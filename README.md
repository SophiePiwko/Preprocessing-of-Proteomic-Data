# 📊 Preprocessing of Proteomic Data 
After running samples through the mass spectrometer, we align the resulting data against a known database to identify the proteins present in our samples. 
This alignment is typically performed using software called DIA-NN (Data Independent Acquisition Neural Network), 
which generates a matrix listing the peptides identified and quantified in each sample, along with their corresponding intensities.

We usually import this matrix into another program called Perseus for data preprocessing. 
However, Perseus has several significant drawbacks, which motivated me to learn Python and enroll in this course. 
During my master's, I analyzed 37 samples using Perseus, but the software frequently crashed, causing me to lose my work—an incredibly frustrating experience.

While Perseus has a user-friendly interface, it comes at the cost of limited flexibility. 
Users can only perform tasks that are already built into the software, with little room for customization or advanced options.

Since I will be handling many more samples during my PhD, I decided not to rely on Perseus for analysis. 
Instead, I plan to use Python to process my data, which will provide greater control, flexibility, and reliability. 

So in this project, written in Python 3, I would like to streamline and automate the preprocessing steps to apply on my previous samples as well as my future samples. 

# 💡 What does this project do? 
The preprocessing steps include:
1. Cleaning the column names
2. Bar plot with Protein ID for each sample
3. Heatmap to look at the overall missingness of the raw data
4. Log transformation and histograms to see if the data is distributed normally
5. Adding annotations and plotting the top 20 annotations
6. Filtering genes based on a threshold of valid values
7. Imputation based on normal distribution
8. Histogram to visualize the overall distribution and the imputed data
9. PCA plot to see the unsupervised clustering of the data

# 📥 Input
The input is the pg_matrix.tsv file received from the DIA-NN software.

# 👉 Requirements 👈
Download the working_directory folder. This folder includes the raw dataset, a requirements.txt file with the needed packages, a zipped file with the annotations and of course the python code. 

# 📝 Instructions
1. Dowload each file in the working_directory file and save them to a new folder called "Preprocessing_of_Proteomic_Data"
   List of files to be downloaded:
   - Preprocessing_of_Proteomic_Data.py
   - mainAnnot.homo_sapiens.zip
   - report.pg_matrix.tsv
   - requirements.txt
2. Unzip the mainAnnot.homo_sapiens.zip folder. Note: make sure that the unzipped file is a .txt file and not a folder with the .txt file. 
3. Open the terminal and navigate to the created folder using cd path of the folder
4. run pip install -r requirements.txt in the terminal
5. run python Preprocessing_of_Proteomic_Data.py in the terminal

# 📤 Output

1_raw_data_cleaned_columns.csv
This file is created after the column names have been cleaned. It contains the raw data prior to any further processing. Two plots in the .html report are generated from this file:
Protein Count per Sample: An interactive bar plot where hovering over bars reveals detailed counts.
Missing Data Heatmap: An interactive heatmap allowing users to zoom in on individual genes.
These visualizations provide an initial quality assessment of the raw dataset.

2_log_transformed.csv
At this step, the data is log-transformed to approximate a normal distribution—an assumption for downstream analysis. The following plot is generated:
Log-Transformed Histogram (Dropdown): A set of histograms (one per sample), viewable via a dropdown menu. These help confirm the normality of the transformed data distribution.

3_with_annotations.csv
This file contains the dataset enriched with KEGG and GO gene annotations. The corresponding plot is:
Top 20 Annotations Bar Plot: An interactive bar plot where hovering reveals detailed annotation names and the number of associated genes.

4_filtered_70percent.csv
The dataset is filtered to retain only rows with at least 70% valid values—standard practice in many analyses. This threshold is adjustable in the code.

5_imputed.csv
Missing values are imputed using a normal distribution. The following plot is generated:
Observed vs. Imputed Histogram: A histogram showing both observed and imputed values, offering a clear view of the overall data distribution. Imputation parameters (width and shift) can be modified in the code.

Finally, a PCA plot is generated to visualize unsupervised clustering of the dataset, helping to identify sample groupings and potential outliers.

# 📖 Data set used
The data set is publicly available and was downloaded from the Proteomics Identifications Database (Pride). 

Wang, H., Lim, K. P., Kong, W., Gao, H., Wong, B. J. H., Phua, S. X., Guo, T., & Goh, W. W. B. (2023). MultiPro: DDA-PASEF and diaPASEF acquired cell line proteomic datasets with deliberate batch effects. Scientific data, 10(1), 858. https://doi.org/10.1038/s41597-023-02779-8
Link to the PRIDE website: https://www.ebi.ac.uk/pride/archive/projects/PXD041421

# 🔓 Conclusion
This project offers an accessible and simple, yet customizable workflow for preprocessing proteomics data, including visualization of the data throughout the processing steps. 
