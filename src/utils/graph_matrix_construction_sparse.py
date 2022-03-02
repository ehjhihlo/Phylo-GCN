import csv
import numpy as np
import pandas as pd

# Read taxa list
taxa = []
with open('./CRC_taxa_group1_species_newest.txt','r') as f:
        for line in f.readlines():
                taxa.append(line[:-1])



for index in range(443):

        list_ = [[],[],[],[],[]] # The list includes each two taxa and their distance [[index],[taxa1],[index],[taxa2],[distance]]

        input_matrix = np.zeros((len(taxa),len(taxa)), dtype=np.float64)
        with open('./species_compare/without_log/sample_'+str(index+1)+'_species_new.csv') as f:
                rows = csv.reader(f)
                for i in rows:
                        list_[0].append(int(i[0]))
                        list_[1].append(i[1])
                        list_[2].append(int(i[2]))
                        list_[3].append(i[3])
                        list_[4].append(float(i[4]))

        list2_ = [[],[]] # Append the index and its normalized value [[index],[value]]
        with open('./taxons/without_log/sample_'+str(index+1)+'_new.csv') as f:
                rows = csv.reader(f)
                for i in rows:
                        list2_[0].append(int(i[0]))
                        list2_[1].append(float(i[2]))

        # Append the distance to list
        for i in range(len(list_[0])):
                input_matrix[list_[0][i]][list_[2][i]] = -list_[4][i]
                input_matrix[list_[2][i]][list_[0][i]] = -list_[4][i]

        # Append the vertex value
        for j in range(len(list2_[0])):
                input_matrix[list2_[0][j]][list2_[0][j]] = list2_[1][j]

        # Write the final input matrix 
        with open('./gcn_laplacian_matrix/without_log/laplacian_matrix_'+str(index+1)+'_species_new.csv', 'w') as f:
                for i in range(len(input_matrix)):
                        for j in range(len(input_matrix)):
                                f.write(str(input_matrix[i][j]))
                                f.write(',')
#                         f.write('\n')
