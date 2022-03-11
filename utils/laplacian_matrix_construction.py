import csv
import numpy as np
import pandas as pd

taxa = []

sample_list = []
with open('../sample_list.txt') as f:
        for line in f.readlines():
                sample_list.append(line[:-1])

with open('./utils_data/abundance_GMPR_norm.txt','r') as f:
        for line in f.readlines():
                taxa.append(line[:-1])

for index in range(len(sample_list)):
        list_ = [[],[],[],[],[]] # The list includes each two taxa and their distance [[index],[taxa1],[index],[taxa2],[distance]]

        # Construct Adjacency matrix, Diagonal matrix, and the final input Laplacian matrix
        adjacency_matrix = np.zeros((len(taxa),len(taxa)), dtype=np.float64)
        diagonal_matrix = np.zeros((len(taxa),len(taxa)), dtype=np.float64)
        input_matrix = np.zeros((len(taxa),len(taxa)), dtype=np.float64)
        
        with open('./utils_data/species_compare/sample_'+str(index+1)+'_species_new.csv') as f:
                rows = csv.reader(f)
                for i in rows:
                        list_[0].append(int(i[0]))
                        list_[1].append(i[1])
                        list_[2].append(int(i[2]))
                        list_[3].append(i[3])
                        list_[4].append(float(i[4]))

        list2_ = [[],[]] # Append the index and its normalized value [[index],[value]]
        with open('./utils_data/taxons/sample_'+str(index+1)+'_new.csv') as f:
                rows = csv.reader(f)
                for i in rows:
                        list2_[0].append(int(i[0]))
                        list2_[1].append(float(i[2]))

        # Append the distance to list to get the adjacency matrix
        for i in range(len(list_[0])):
                adjacency_matrix[list_[0][i]][list_[2][i]] = list_[4][i]
                adjacency_matrix[list_[2][i]][list_[0][i]] = list_[4][i]               

        # Append the vertex value to get the diagonal matrix
        for j in range(len(list2_[0])):
                diagonal_matrix[list2_[0][j]][list2_[0][j]] = list2_[1][j]

        # Construct the Laplacian matrix by D - A
        input_matrix = diagonal_matrix - adjacency_matrix
    
        # Write the final input matrix 
        with open('../data/Phylo-GCN/matrix/laplacian_matrix_'+str(index+1)+'_species_new.csv', 'w') as f:
                for i in range(len(input_matrix)):
                        for j in range(len(input_matrix)):
                                f.write(str(input_matrix[i][j]))
                                f.write(',')
#                         f.write('\n')