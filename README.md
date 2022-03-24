# Phylo-GCN
**A Graph Convolutional Network with Phylogenic Tree Embedded for Colorectal Cancer Classification**  
  
Accurate prediction of diseases from metagenomic samples is important for development of non-invasive disease diagnosis. We introduce Phylo-GCN, a graph classification GCN model with phylogenetic tree distance included for colorectal cancer (CRC) prediction. Three public CRC datasets, total 443 samples are collected for training. In each sample, we do sequence analysis and take a format of 2D matrix considering the abundance of taxa and the relative phylogeny tree distance between taxa, we constructed these information into a graph. Then we train a deep learning model to predict whether the sample represents cancer or not, and compare our model with some phylogeny tree distance embedded models and traditional machine learning models.    
## Clone Repository
Users can download the sources of this project by executing the commands below:
``` 
git clone https://github.com/ehjhihlo/Phylo-GCN.git  
cd Phylo-GCN  
```
## data  
The input dataset used for training. For Phylo-GCN, laplacian matrices and label is being used.  
  
**matrix_1137** folder contains the laplacian matrices with all 1137 taxa  
  
**matrix** folder contains the laplacian matrices with all 496 taxa  
  
**label** folder contains the sample and its corresponding label  
  
## train  
The code for training the model, including Phylo-GCN demo and other models.    
## trees
The newick file of phylogeny tree, including all 1137 taxa in this project.  
## utils  
The files, codes for calculating phylogeny distance and laplacian matrix.  
  
**utils_data:** Contains all data used for calculating phylogeny distance and laplacian matrix.  
abundance_GMPR_norm.csv: GMPR normalized abundance table (1137 taxa).  
abundance_GMPR_norm_496.csv: GMPR normalized abundance table (496 taxa).  
taxa_list.txt: The list for all taxa in this project.  
taxa_list_496.txt: The list for 496 taxa in this project.  
species_compare folder: The calculated distance for each 2 taxa.  
taxons folder: The taxa and its abindance in each sample.  
  
**GMPR-normalize.Rmd:**  For GMPR normalization on original abundance table.
  
**get_distance.py:** For phylogeny distance calculation.
  .
**laplacian_matrix_construction.py:**  For laplacian matrix construction.
  
**phylogentic_tree_species.Rmd:**  For phylogeny tree newick file construction.
