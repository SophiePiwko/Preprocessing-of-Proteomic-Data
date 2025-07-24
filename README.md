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
Several .csv files are produced throughout the pipeline, with each file representing the output of a distinct processing stage. The files are numbered, corresponding to the order of processing. In addition, an .html file is created including all the plots generated throughout the pipeline. 
"1_raw_data_cleaned_columns.csv" is generated after the column names and hava been cleaned. This file includes the raw data before any processing. The "Protein Count per Sample" plot and the "Missing Data Heatmap" plot in the html file are generated from this data, giving an overall overview of the raw data. This is helpfull in assessing the quality of the generated data set. Hovering over the bars in the barplot gives more detailed information. The heatmap is interactive, allowing the zoom in on the different genes. Next the data is log transformed and "2_log_transformed.csv" is generated. This file includes the log transformed values of the data set. The "log_transformed_histogram_dropdown" plots in the html file are generated form this .csv file. This allows the user to see that the data is distributed normally which is the assumption for the next steps. Notice that there is a drop down menu to look at the histogram for each sample seperatly. Next the "3_with_annotations.csv" file is generated which includes the KEGG and GO annotations associated with the genes in the data set. In addition, the "top20_annotations_barplot" is created in the html file. Hovering over the plot gives more detailed information about the annotation and how many genes are associated with this annotation. Next the data is filtered based on 70% valid values, creating the "4_filtered_70percent.csv" file. Filtering the data based on 70% valid value is standard and usually the default setting. This value however can be easily adjusted in the code if need be. Next there is data imputation based on normal distribution, generating the "5_imputed.csv" file as well as the "observed_vs_imputed_histogram" plot in the html file. Also here the width and shift of the imputed data can be easily adjusted in the code. The histogram shows the observed as well as the imputed data, allowing to see the overall distribution of the imputed data set at one glance. Lastly the "pca_unsupervised_clustering" is generated showing the unsupervised clustering of the data. 

# 📖 Data set used
The data set is publicly available and was downloaded from the Proteomics Identifications Database (Pride). 

Wang, H., Lim, K. P., Kong, W., Gao, H., Wong, B. J. H., Phua, S. X., Guo, T., & Goh, W. W. B. (2023). MultiPro: DDA-PASEF and diaPASEF acquired cell line proteomic datasets with deliberate batch effects. Scientific data, 10(1), 858. https://doi.org/10.1038/s41597-023-02779-8
Link to the PRIDE website: https://www.ebi.ac.uk/pride/archive/projects/PXD041421

# 🔓 Conclusion
This project offers an accessible and simple, yet customizable workflow for preprocessing proteomics data, including visualization of the data throughout the processing steps. 
