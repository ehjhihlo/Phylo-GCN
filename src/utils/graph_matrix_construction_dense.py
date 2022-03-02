import csv
import numpy as np
import pandas as pd

# Calculating how many species in each sample, and get the abundance value of each species in every sample
count = 0
count_list = []
vertex_value = []
for index in range(443):
#         with open('./taxons/sample_'+str(index+1)+'.csv','r') as f:
#                 for line in f.readlines():
#                         count+=1
#                 count = count-1 #remove first row
#                 count_list.append(count)
#                 count = 0 #計數器歸零
        vertex_value.append([])
        data = pd.read_csv('./taxons/sample_'+str(index+1)+'.csv')
        count_list.append(len(data))
        for i in range(len(data)):
                vertex_value[index].append(data['Value'][i])
# print(count_list)
# print(vertex_value)

for index in range(443):
        distance = []
        with open('./species_compare/sample_'+str(index+1)+'_species.csv','r') as f:
                for line in f.readlines():
                        if(line[-3]==','): # The distance is 5, which contains only one digit
                                distance.append(int(line[-2]))
                        else: # The distance is 10 or larger, containing two digits
                                distance.append(int(line[-3]+line[-2]))

                # Define the adjacency matrix. In each sample, the dimension is (# of species * # of species)
                adjacency = np.zeros((count_list[index],count_list[index]), dtype=np.float64)
                # Fill in the distance into the adjacency matrix 
                counter = 0
                for i in range(count_list[index]): # Fill the value in the upright side of the matrix
                        for j in range(i+1, count_list[index]):
                                adjacency[i][j] = distance[counter]
                                counter+=1
                for i in range(count_list[index]): # Fill the value in the downleft side of the matrix
                        for j in range(count_list[index]):
                                if i > j:
                                        adjacency[i][j] = adjacency[j][i]    

                # Define the adjacency matrix. In each sample, each 2 species connect to each other, so the value is (# of species - 1)
                diagonal = np.zeros((count_list[index],count_list[index]), dtype=np.float64)
                for i in range(count_list[index]):
                        for j in range(count_list[index]):
                                if i == j:
                                        diagonal[i][j] = vertex_value[index][i]

                # Laplacian Eigenvectors
                I = np.eye(count_list[index])
                laplacian = diagonal-adjacency
        #         print(I - ((diagonal**-0.5)*(diagonal-adjacency)*(diagonal**-0.5)))

        # Save the laplacian matrix
        with open('./gcn_laplacian_matrix/laplacian_matrix_'+str(index+1)+'_species.csv', 'w') as f:
                for i in range(len(laplacian)):
                        for j in range(len(laplacian[0])):
                                f.write(str(laplacian[i][j]))
                                f.write(',')
                        f.write('\n')

