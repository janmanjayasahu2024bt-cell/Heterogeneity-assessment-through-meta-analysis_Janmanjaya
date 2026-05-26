# Heterogeneity-assessment-through-meta-analysis_Janmanjaya
Aim of the project: To achieve higher machine learning-based classification precision for psychiatric disorders using ensemble learning.

Data type: Microarray Gene Expression Datasets and RNA-Seq data

Packages used: Refer to scripts

About the scripts: The gene expression data was processed using the scripts shared. The numbers before the title of the script indicate the sequence in which the analysis was performed. Below are the details of the analysis performed using the R scripts.

Pre-processing: The script includes the creation of a meta-file from the gene expression datasets followed by a random selection of training and testing datasets. The script also includes the processing of an independent dataset and RNA-Seq samples.
Quantile normalization: Training data was independently normalized. While test data, independent dataset, and RNA-Seq samples were normalized with respect to the train data. This was achieved using quantile targets from the respective iterations of train data.
Batch correction: Different batches in train data were batch corrected using comBat. The test, independent dataset, and RNA-Seq samples were batch-corrected with batch-corrected train data as a reference.
Feature selection: Differential gene expression analysis as feature selection method. Feature genes identified from each train data were used to build machine learning models.
Model development and test data prediction: Models were built using different numbers of feature genes and varying kernels. Models were tested using test datasets. Models and test data predictions were saved for further analysis.
Prediction for independent dataset and RNA-Seq data: Quantile normalized, and batch-corrected samples from the independent dataset and RNA-Seq were predicted using microarray-based machine learning models.
Ensemble learning: Test data prediction of SVM and PAM models were ensemble using the Boolean operator “AND”. Ensemble models were evaluated based on precision.
Contact: Vipul : vipulwagh31@gmail.com
