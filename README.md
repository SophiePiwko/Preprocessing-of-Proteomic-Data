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
1. Heatmap to look at the overall distribution and missingness of the data
2. Log transformation
3. Histogram to see if the data is distributed normally
4. Filtering the rows based on a threshold of valid values
5. Imputation and histogram to make sure that the distribution isn't scewed
6. PCA plot to see the unsupervised clustering of the data

# ➡️ Input and Output ⬅️
The input will be the .tsv matrix received from DIA-NN and the output will be a matrix after the mentioned preprocessing steps 
as well as heatmaps, histograms and PCA plots generated during the preprocessing. 

# 🔓 Conclusion
This project offers an accessible and simple, yet customizable workflow for preprocessing proteomics data, including visualization of the data throughout the processing steps. 
