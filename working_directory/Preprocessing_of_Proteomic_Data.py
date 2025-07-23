import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import numpy as np
import re
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import math
from collections import defaultdict

# === CONFIGURATION ===
#r"C:/Users/sophiep.WISMAIN/Desktop/
CONFIG = {
    'file_path': r"Preprocessing_of_Proteomic_Data/Data/report.pg_matrix.tsv",
    'output_path': r"Preprocessing_of_Proteomic_Data/output",
    'validity_threshold': 0.7,
    'imputation_params': {'shift': 1.5, 'scale': 0.5},
    'exclude_cols': ['Protein.Ids', 'Protein.Names', 'Genes', 'First.Protein.Description']
}

os.makedirs(CONFIG['output_path'], exist_ok=True)
file_path = CONFIG['file_path']

# === Saving the figures for the summary report ===
all_figures = [] 

# === Step 1: Read raw data, set index, clean column names ===

def load_and_clean_data(file_path):
    """Load data and clean column names."""
    try:
        df = pd.read_csv(file_path, sep='\t', index_col=0, low_memory=False)
        print("File loaded successfully.")
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        raise
    except pd.errors.EmptyDataError:
        print("The file is empty.")
        raise
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        raise

    # Clean column names using regex
    df.columns = df.columns.str.replace(
        r'F:\\WOSP00101-|DIA_|DIA-|DIANN\\|WOSP00101_|CAD20211212chenc_|N20211212chenc_|60min_|\.d$',
        '', regex=True, flags=re.IGNORECASE
    )
    return df

df = load_and_clean_data(CONFIG['file_path'])


# Save cleaned raw data
df.to_csv(os.path.join(CONFIG['output_path'], "step_1_raw_data_cleaned_columns.csv"))

print(f"step_1_raw_data_cleaned_columns.csv has been saved to:\n{CONFIG['output_path']}")

# === Step 2: Overview of the data ===
total_items = df.size
non_missing = df.count().sum()
missing = df.isna().sum().sum()
print(f"Matrix dimensions: {df.shape[0]} genes × {df.shape[1]} samples")
print(f"Total items: {total_items}")
print(f"Non-missing values: {non_missing}")
print(f"Missing values: {missing}")

# === Step 2.1 Protein Count Per Sample ===
# Exclude metadata columns
sample_protein_counts = df.drop(columns=CONFIG['exclude_cols'], errors='ignore').notna().sum()
sample_protein_counts = sample_protein_counts.sort_values(ascending=False)
sample_protein_counts.name = "Protein_Count"

# Plot as bar chart
fig = px.bar(
    sample_protein_counts,
    x=sample_protein_counts.index,
    y="Protein_Count",
    title="Number of Protein IDs Detected per Sample",
    labels={"x": "Sample", "Protein_Count": "Number of Proteins Detected"},
    height=600
)

fig.update_layout(
    xaxis_tickangle=45,
    xaxis_title="Sample",
    yaxis_title="Protein Count",
    bargap=0.3
)

# Save and include in summary HTML
all_figures.append(("Protein Count per Sample", fig))

# === Step 2.2: Heatmap to visualize missing data presence ===
df_presence = df.drop(columns=CONFIG['exclude_cols'], errors='ignore')

presence_matrix = (~df_presence.isna()).astype(int)

colors = ['lightgrey', 'steelblue']
color_labels = {0: "Missing", 1: "Present"}

fig = px.imshow(
    presence_matrix,
    color_continuous_scale=colors,
    x=presence_matrix.columns,
    y=presence_matrix.index,
    aspect="auto"
)

fig.update_coloraxes(showscale=False)

for value, label in color_labels.items():
    fig.add_scatter(
        x=[None],
        y=[None],
        mode='markers',
        marker=dict(size=10, color=colors[value]),
        legendgroup=label,
        showlegend=True,
        name=label
    )

fig.update_layout(
    title="Data Presence: Blue = Present, Grey = Missing",
    xaxis_title="Samples",
    yaxis_title="Genes",
    height=1000
)

all_figures.append(("Missing Data Heatmap", fig))

# === Step 3: Log-transform data ===
df_numeric = df.drop(columns=CONFIG['exclude_cols'], errors='ignore').apply(pd.to_numeric, errors='coerce')
log_df = np.log1p(df_numeric)
log_df.to_csv(os.path.join(CONFIG['output_path'], "step_3_log_transformed.csv"))

print(f"step_3_log_transformed.csv has been saved to:\n{CONFIG['output_path']}")

# === Step 3.1: Histograms of log-transformed expression by sample ===
plot_cols = [col for col in log_df.columns if col not in CONFIG['exclude_cols']]

fig = go.Figure()

for i, sample in enumerate(plot_cols):
    fig.add_trace(go.Histogram(
        x=log_df[sample].dropna(),
        nbinsx=50,
        name=sample,
        visible=(i == 0),
        marker_color='teal',
        opacity=0.75
    ))

dropdown_buttons = [
    dict(label=sample,
         method="update",
         args=[{"visible": [i == j for j in range(len(plot_cols))]},
               {"title": f"Log-Transformed Expression - {sample}"}])
    for i, sample in enumerate(plot_cols)
]

fig.update_layout(
    updatemenus=[dict(
        buttons=dropdown_buttons,
        direction="down",
        x=1.05,
        y=1,
        showactive=True
    )],
    title=f"Log-Transformed Expression - {plot_cols[0]}",
    xaxis_title="Log(Expression)",
    yaxis_title="Frequency",
    showlegend=False,
    bargap=0.1
)

all_figures.append(("log_transformed_histogram_dropdown", fig))

# === Step 4: Add annotations from KEGG and GO ===
kegg_annot_path = r"C:/Users/sophiep.WISMAIN/Desktop/Perseus 2.0.11/bin/conf/annotations/mainAnnot.homo_sapiens.txt"
kegg_df = pd.read_csv(kegg_annot_path, sep='\t', low_memory=False)

log_df_with_annot = log_df.copy()

# Normalize UniProt IDs (uppercase & strip spaces)
log_df_with_annot['UniProt_ID'] = log_df_with_annot.index.str.upper().str.strip()
kegg_df['UniProt'] = kegg_df['UniProt'].str.upper().str.strip()

def build_annotation_map(df, id_col, annot_col):
    mapping = defaultdict(list)
    for _, row in df.iterrows():
        uni_ids = re.split(r'[;,]', row[id_col])
        for uni_id in uni_ids:
            uni_id = uni_id.strip()
            if uni_id:
                annot_val = row[annot_col]
                if pd.notna(annot_val):  # Only add if annotation is not NaN
                    mapping[uni_id].append(str(annot_val))
    # Join multiple annotations with '; '
    return {k: '; '.join(v) for k, v in mapping.items()}


# Build mappings for KEGG and GO annotations
kegg_map = build_annotation_map(kegg_df, 'UniProt', 'KEGG name')
gocc_map = build_annotation_map(kegg_df, 'UniProt', 'GOCC name')
gomf_map = build_annotation_map(kegg_df, 'UniProt', 'GOMF name')

# Map annotations to your data's UniProt IDs
log_df_with_annot['KEGG_Pathway'] = log_df_with_annot['UniProt_ID'].map(kegg_map)
log_df_with_annot['GO_CC'] = log_df_with_annot['UniProt_ID'].map(gocc_map)
log_df_with_annot['GO_MF'] = log_df_with_annot['UniProt_ID'].map(gomf_map)

log_df_with_annot.to_csv(os.path.join(CONFIG['output_path'], "step_4_with_annotations.csv"))

print(f"step_4_with_annotations.csv has been saved to:\n{CONFIG['output_path']}")

# === Step 4.1: Top 20 annotations for each type and plot ===
def top_20_counts(df, col_name, type_name):
    counts = df[col_name].value_counts().reset_index()
    counts.columns = ['Annotation', 'Gene_Count']
    counts['Type'] = type_name
    return counts.head(20)

kegg_counts = top_20_counts(log_df_with_annot, 'KEGG_Pathway', 'KEGG Pathway')
gocc_counts = top_20_counts(log_df_with_annot, 'GO_CC', 'GO Cellular Component')
gomf_counts = top_20_counts(log_df_with_annot, 'GO_MF', 'GO Molecular Function')

combined_counts = pd.concat([kegg_counts, gocc_counts, gomf_counts], ignore_index=True)

fig = px.bar(
    combined_counts,
    x='Annotation',
    y='Gene_Count',
    color='Type',
    facet_col='Type',
    facet_col_wrap=1,
    title='Top 20 Annotations: KEGG, GO CC, GO MF',
    labels={'Annotation': 'Annotation Term', 'Gene_Count': 'Number of Genes'},
    height=1000
)

fig.update_layout(
    xaxis=dict(showticklabels=False),
    showlegend=False,
    title='Top 20 Annotations: KEGG, GO CC, GO MF'
)

all_figures.append(("top20_annotations_barplot", fig))

# === Step 5: Filter genes with at least validity_threshold non-NA values ===
exclude_cols_filter = CONFIG['exclude_cols'] + ['KEGG_Pathway', 'GO_CC', 'GO_MF']

numeric_cols = [col for col in log_df_with_annot.select_dtypes(include='number').columns if col not in exclude_cols_filter]

min_valid = math.ceil(len(numeric_cols) * CONFIG['validity_threshold'])

non_na_counts = log_df_with_annot[numeric_cols].notna().sum(axis=1)

filtered_df = log_df_with_annot[non_na_counts >= min_valid]

print(f"Filtered from {log_df_with_annot.shape[0]} to {filtered_df.shape[0]} genes (rows) based on {int(CONFIG['validity_threshold']*100)}% valid sample values.")
print("Top missingness by sample:")
print(log_df_with_annot[numeric_cols].isna().sum().sort_values(ascending=False).head(60))

filtered_df.to_csv(os.path.join(CONFIG['output_path'], "step_5_filtered_70percent.csv"))

print(f"step_5_filtered_70percent.csv has been saved to:\n{CONFIG['output_path']}")

# === Step 6: Impute missing data by normal distribution ===
def impute_missing_values(df, shift=1.5, scale=0.5):
    """Impute missing values in a DataFrame with normally distributed values 
    based on a left-shifted distribution per column."""
    df_imputed = df.copy()
    
    for col in df.columns:
        missing_mask = df[col].isna()
        if missing_mask.sum() > 0:
            observed_values = df[col][~missing_mask]
            mean_val = observed_values.mean()
            std_val = observed_values.std()

            if std_val == 0 or np.isnan(std_val):
                imputed_values = np.full(missing_mask.sum(), mean_val)
            else:
                imputed_values = np.random.normal(
                    loc=mean_val - shift * std_val,
                    scale=scale * std_val,
                    size=missing_mask.sum()
                )

            df_imputed.loc[missing_mask, col] = imputed_values

    return df_imputed

log_df_imputed = impute_missing_values(log_df, **CONFIG['imputation_params'])
log_df_imputed.to_csv(os.path.join(CONFIG['output_path'], "step_6_imputed.csv"))

print(f"step_6_imputed.csv has been saved to:\n{CONFIG['output_path']}")

# === Step 6.1: Plot histograms of observed vs imputed values ===
plot_cols_hist = [col for col in log_df.columns if col not in CONFIG['exclude_cols']]

fig = go.Figure()

for i, sample in enumerate(plot_cols_hist):
    real_values = log_df[sample].dropna()
    imputed_values = log_df_imputed[sample][log_df[sample].isna()]

    fig.add_trace(go.Histogram(
        x=real_values,
        nbinsx=50,
        name='Observed',
        visible=(i == 0),
        marker_color='teal',
        opacity=0.7
    ))

    fig.add_trace(go.Histogram(
        x=imputed_values,
        nbinsx=50,
        name='Imputed',
        visible=(i == 0),
        marker_color='red',
        opacity=0.7
    ))

dropdown_buttons = []
for i, sample in enumerate(plot_cols_hist):
    visible = [False] * len(plot_cols_hist) * 2
    visible[2*i] = True
    visible[2*i + 1] = True

    dropdown_buttons.append(dict(
        label=sample,
        method="update",
        args=[{"visible": visible},
              {"title": f"Histogram of Observed and Imputed Values - {sample}"}]
    ))

fig.update_layout(
    updatemenus=[dict(buttons=dropdown_buttons, direction="down", x=1.05, y=1.0)],
    title=f"Histogram of Observed and Imputed Values - {plot_cols_hist[0]}",
    xaxis_title="Log(Expression)",
    yaxis_title="Frequency",
    bargap=0.15,
    barmode='overlay',
    showlegend=True
)

all_figures.append(("observed_vs_imputed_histogram", fig))

# === Step 7: Unsupervised clustering with PCA ===
exclude_cols_pca = CONFIG['exclude_cols'] + ['KEGG_Pathway', 'GO_CC', 'GO_MF']

sample_data = log_df_imputed.drop(columns=[col for col in exclude_cols_pca if col in log_df_imputed.columns])

data_T = sample_data.T

scaler = StandardScaler()
scaled_data = scaler.fit_transform(data_T)

scaled_df = pd.DataFrame(scaled_data, index=data_T.index, columns=data_T.columns)

pca = PCA(n_components=2)
pca_result = pca.fit_transform(scaled_df)

pca_df = pd.DataFrame(pca_result, columns=["PC1", "PC2"], index=scaled_df.index)
pca_df.reset_index(inplace=True)
pca_df.rename(columns={"index": "Sample"}, inplace=True)

fig = px.scatter(
    pca_df, x="PC1", y="PC2", text="Sample",
    title="PCA of Samples (Unsupervised Clustering)",
    width=800, height=600
)

fig.update_traces(textposition='top center')
fig.update_layout(
    xaxis_title="Principal Component 1",
    yaxis_title="Principal Component 2",
    showlegend=False
)

all_figures.append(("pca_unsupervised_clustering", fig))

# === Step 8: Generate a summary report with visualizations ===

report_html_parts = []

# Use Plotly's full_html=True on the *first* figure to include JS
for i, (title, fig) in enumerate(all_figures):
    fig_html = pio.to_html(
        fig,
        include_plotlyjs='cdn' if i == 0 else False,  # Only include JS once
        full_html=False
    )
    section_html = f"<h2>{title}</h2>\n{fig_html}<hr style='margin:40px 0;'>"
    report_html_parts.append(section_html)

# Full HTML document
full_report_html = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Proteomics Data Summary Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; }}
        h1 {{ text-align: center; }}
        h2 {{ margin-top: 50px; }}
    </style>
</head>
<body>
    <h1>Summary Report: Proteomics Data</h1>
    {"".join(report_html_parts)}
</body>
</html>
"""

report_path = os.path.join(CONFIG['output_path'], "summary_report_all_plots.html")
with open(report_path, "w", encoding="utf-8") as f:
    f.write(full_report_html)

print(f"✅ Generated plots have been saved in the summary HTML report to:\n{report_path}")