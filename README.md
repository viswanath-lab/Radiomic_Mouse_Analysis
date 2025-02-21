# Radiomic_Mouse_Analysis

File Structure

main: This folder contains the scripts to read mouse dicom data and annotations, extract radiomic features, evaluate features within a subset of features, sort by p-values, train a model and evaluate its performance on other sets

-- mouse2dicom: Script to read dicom directories and convert them to the mha format used in feature extraction

-- feature_extraction_2d_KFMC_scripts: Script that extracts radiomic features from annotations. Requires modification for each dataset and specific annotation ID number. Each annotated slice is treated as seperate.

-- feat_select: Needs to be modified for each dataset being run but stacks the mouse data for each annotated slice and mouse. then for the first dataset which is used for training also performs feature analysis and ranks each feature family by discriminability between groups.

--mouse_classFN: Builds the random forest model used for the rest of the analysis. Trains the model on a subset of training patients and validates on the testing subset. 

--external_validation: Loads the model built in mouse_classFN and applies it to other datasets. 


Feature_Extraction: Contains the code for each of the featue sets extracted in this study

--extract2DFeatureInfo: main module that handles feature extraction subfunction 2DFeatureStats can be changed to adjust which statistics are computed for the features. 


image functions: This folder contains the code needed to convert an image to mha as well as handles generating some additional plots and figures used in the study

# Citation

This code was utilized as part of analysis done for the following paper:

Maneesh Dave, Atul Dev, Rodrigo A Somoza, Nan Zhao, Satish E Viswanath, Pooja Rani Mina, Prathyush Chirra, Verena Carola Obmann, Ganapati H Mahabeleshwar, Paola Menghini, Blythe Durbin-Johnson, Jan Nolta, Christopher Soto, Abdullah Osme, Lam T Khuat, William J Murphy, Arnold I Caplan, Fabio Cominelli, "MSCs mediate long-term efficacy in a Crohn's disease model by sustained anti-inflammatory macrophage programming via efferocytosis", NPJ Regen Med. 2024 Jan 20;9(1):6.  doi: 10.1038/s41536-024-00347-1.

