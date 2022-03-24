# Phylo-GCN
**A Graph Convolutional Network with Phylogenic Tree Embedded for Colorectal Cancer Classification**  
Accurate prediction of diseases from metagenomic samples is important for development of non-invasive disease diagnosis. We introduce Phylo-GCN, a graph classification GCN model with phylogenetic tree distance included for colorectal cancer (CRC) prediction. Three public CRC datasets, total 443 samples are collected for training. In each sample, we do sequence analysis and take a format of 2D matrix considering the abundance of taxa and the relative phylogeny tree distance between taxa, we constructed these information into a graph. Then we train a deep learning model to predict whether the sample represents cancer or not, and compare our model with some phylogeny tree distance embedded models and traditional machine learning models.    
## Clone Repository
``` 
git clone https://github.com/ehjhihlo/Phylo-GCN/.git  
cd Phylo-GCN  
```
## data  
Including the input dataset used for training.  
## train  
The code for training the model.  
## trees
The newick file of phylogeny tree, including all 1137 taxa in this project.  
## utils  
The files, codes for calculating phylogeny distance and laplacian matrix.  
