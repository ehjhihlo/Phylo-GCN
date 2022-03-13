# PartI: Get the taxa of each sample and it's normalized abundanced

import csv
import numpy as np
import pandas as pd

path = '../sample_list.txt'
data = pd.read_csv("./utils_data/abundance_GMPR_norm.csv")
data.head()

def readSample(path):
        sample_list = []
        with open(path) as f:
                for line in f.readlines():
                        sample_list.append(line[:-1])
        return sample_list

# In each sample, get the taxa and normalized value
def getTaxaValue(): # Version 1
        sample_list = readSample(path)
        for i in range(len(sample_list)):
                # In every sample, check out which taxon has performance on it
                globals()['sample_%s' % (i+1)] = data.loc[data[sample_list[i]] != 0]
                
                select = ['Kingdom','Phylum','Class','Order','Family','Genus','Species',sample_list[i]]
                globals()['sample_%s' % (i+1)] = globals()['sample_%s' % (i+1)].loc[:,select]
                globals()['sample_%s' % (i+1)] = globals()['sample_%s' % (i+1)].rename(columns={sample_list[i]: 'Value'})             
                globals()['sample_%s' % (i+1)].to_csv('./taxons/new/sample_'+str(i+1)+'.csv', index=False)

def getTaxaValue_2(): # Version 2
        sample_list = readSample(path)
        for i in range(len(sample_list)):
                # In every sample, check out which taxon has performance on it
                globals()['sample_%s' % (i+1)] = data.loc[data[sample_list[i]] != 0]
                select = ['Kingdom','Phylum','Class','Order','Family','Genus','Species',sample_list[i]]
                globals()['sample_%s' % (i+1)] = globals()['sample_%s' % (i+1)].loc[:,select]
                globals()['sample_%s' % (i+1)] = globals()['sample_%s' % (i+1)].rename(columns={sample_list[i]: 'Value'})
                
                with open('./utils_data/taxons/sample_'+str(i+1)+'_new.csv', 'w') as f:
                        k_ = globals()['sample_%s' % (i+1)]['Kingdom'].tolist()
                        p_ = globals()['sample_%s' % (i+1)]['Phylum'].tolist()
                        c_ = globals()['sample_%s' % (i+1)]['Class'].tolist()
                        o_ = globals()['sample_%s' % (i+1)]['Order'].tolist()
                        f_ = globals()['sample_%s' % (i+1)]['Family'].tolist()
                        g_ = globals()['sample_%s' % (i+1)]['Genus'].tolist()
                        s_ = globals()['sample_%s' % (i+1)]['Species'].tolist()
                        value_ = globals()['sample_%s' % (i+1)]['Value'].tolist()
                        index_ = globals()['sample_%s' % (i+1)].index.values.tolist()

                        for j in range(len(k_)):
                                f.write(str(index_[j]))
                                f.write(',')
                                f.write(str('k__'+k_[j]+'|p__'+p_[j]+'|c__'+c_[j]+'|o__'+o_[j]+'|f__'+f_[j]+'|g__'+g_[j]+'|s__'+s_[j]))
                                f.write(',')
                                f.write(str(value_[j]))
                                f.write(',')
                                f.write('\n')


def writeSamplelist():
        # sample_list.csv represent the relationship between sample number(sample_1, sample2,...) and sample run number(SRRXXXXXX)
        sample_list = readSample(path)
        with open('../sample_list.csv','w') as f:
                f.write('number, sample name')
                f.write('\n')
                for i in range(len(sample_list)):
                        f.write('sample_%s' % (i+1) + ',')
                        f.write(sample_list[i]+',')
                        f.write('\n')


# PartII: Calculate the phylogeny distance

# Version 1
def getTaxaTreeDistance():
        sample_list = readSample(path)
        taxa_list = []
        tree_distance_list = []
        for index in range(len(sample_list)):
                csvreader = csv.reader(open('./utils_data/taxons/sample_'+str(index+1)+'.csv'))
                matrix = []
                for i in csvreader:
                        matrix.append(i)
                del matrix[0]
                for i in range(len(matrix)):
                        matrix[i] = matrix[i][:-1]
                taxa_list.append(matrix)
                species_list = []

                for i in range(len(matrix)):
                        for j in range(i+1,len(matrix)):
                                list_ = matrix[i] + matrix[j]
                                if len(list_) != 14:
                                        print("Error!")
                                        break
                                species_list.append(list_)

                distance_list = [None for i in range(len(species_list))]

                for i in range(len(species_list)):
                        count = 0 # Calculate the number of 'None' in each taxon
                        # Different Kingdom
                        if species_list[i][0] != species_list[i][7]:
                                for j in range(len(species_list[i])):
                                        if species_list[i][j] == 'None':
                                                count+=1
                                distance_list[i] = (14-count)*5
                                count = 0

                        # Same Kingdom and different Phylum
                        elif species_list[i][1] != species_list[i][8] and species_list[i][1] != 'None' and species_list[i][8] != 'None':
                                for j in range(len(species_list[i])):
                                        if species_list[i][j] == 'None':
                                                count+=1
                                distance_list[i] = (12-count)*5
                                count = 0
                        elif species_list[i][1] != species_list[i][8] and species_list[i][1] == 'None' and species_list[i][8] != 'None':
                                for j in range(7,14):
                                        if species_list[i][j] == 'None':
                                                count+=1
                                distance_list[i] = 30-(5*count)
                        elif species_list[i][1] != species_list[i][8] and species_list[i][1] != 'None' and species_list[i][8] == 'None':
                                for j in range(0,7):
                                        if species_list[i][j] == 'None':
                                                count+=1
                                distance_list[i] = 30-(5*count)

                        # Same Phylum and different Class
                        elif species_list[i][2] != species_list[i][9] and species_list[i][2] != 'None' and species_list[i][9] != 'None':
                                for j in range(len(species_list[i])):
                                        if species_list[i][j] == 'None':
                                                count+=1
                                distance_list[i] = (10-count)*5
                                count = 0
                        elif species_list[i][2] != species_list[i][9] and species_list[i][2] == 'None' and species_list[i][9] != 'None':
                                for j in range(7,14):
                                        if species_list[i][j] == 'None':
                                                count+=1
                                distance_list[i] = 25-5*(count)
                        elif species_list[i][2] != species_list[i][9] and species_list[i][2] != 'None' and species_list[i][9] == 'None':
                                for j in range(0,7):
                                        if species_list[i][j] == 'None':
                                                count+=1
                                distance_list[i] = 25-5*(count)

                        # Same Class and different Order
                        elif species_list[i][3] != species_list[i][10] and species_list[i][3] != 'None' and species_list[i][10] != 'None':
                                for j in range(len(species_list[i])):
                                        if species_list[i][j] == 'None':
                                                count+=1
                                distance_list[i] = (8-count)*5
                                count = 0
                        elif species_list[i][3] != species_list[i][10] and species_list[i][3] == 'None' and species_list[i][10] != 'None':
                                for j in range(7,14):
                                        if species_list[i][j] == 'None':
                                                count+=1
                                distance_list[i] = 20-5*(count)
                        elif species_list[i][3] != species_list[i][10] and species_list[i][3] != 'None' and species_list[i][10] == 'None':
                                for j in range(0,7):
                                        if species_list[i][j] == 'None':
                                                count+=1
                                distance_list[i] = 20-5*(count)

                        # Same Order and different Family
                        elif species_list[i][4] != species_list[i][11] and species_list[i][4] != 'None' and species_list[i][11] != 'None':
                                for j in range(len(species_list[i])):
                                        if species_list[i][j] == 'None':
                                                count+=1
                                distance_list[i] = (6-count)*5
                                count = 0
                        elif species_list[i][4] != species_list[i][11] and species_list[i][4] == 'None' and species_list[i][11] != 'None':
                                for j in range(7,14):
                                        if species_list[i][j] == 'None':
                                                count+=1
                                distance_list[i] = 15-5*(count-1)
                        elif species_list[i][4] != species_list[i][11] and species_list[i][4] != 'None' and species_list[i][11] == 'None':
                                for j in range(0,7):
                                        if species_list[i][j] == 'None':
                                                count+=1
                                distance_list[i] = 15-5*(count-1)

                        # Same Family and different Genus
                        elif species_list[i][5] != species_list[i][12] and species_list[i][5] != 'None' and species_list[i][12] != 'None':
                                for j in range(len(species_list[i])):
                                        if species_list[i][j] == 'None':
                                                count+=1
                                distance_list[i] = (4-count)*5
                                count = 0
                        elif species_list[i][5] != species_list[i][12] and species_list[i][5] == 'None' and species_list[i][12] != 'None':
                                for j in range(7,14):
                                        if species_list[i][j] == 'None':
                                                count+=1
                                distance_list[i] = 10-5*(count-1)
                        elif species_list[i][5] != species_list[i][12] and species_list[i][5] != 'None' and species_list[i][12] == 'None':
                                for j in range(0,7):
                                        if species_list[i][j] == 'None':
                                                count+=1
                                distance_list[i] = 10-5*(count-1)
                        # Same Genus and different Species
                        elif species_list[i][6] != species_list[i][13] and species_list[i][6] != 'None' and species_list[i][13] != 'None':
                                for j in range(len(species_list[i])):
                                        if species_list[i][j] == 'None':
                                                count+=1
                                distance_list[i] = (2-count)*5
                                count = 0
                        elif species_list[i][6] != species_list[i][13] and species_list[i][6] == 'None' and species_list[i][13] != 'None':
                                distance_list[i] = 5
                        elif species_list[i][6] != species_list[i][13] and species_list[i][6] != 'None' and species_list[i][13] == 'None':
                                distance_list[i] = 5

                tree_distance_list.append(distance_list)
        return taxa_list, tree_distance_list

# Get the total tree distance of each sample
def getSumTreeDistance(taxa_list):
        sample_list = readSample(path)
        sum_distance_list = []
        count_none = 0
        for i in range(len(sample_list)):
                sum_distance = 0
                for j in range(len(taxa_list[i])):
                        for k in range(len(taxa_list[i][j])):
                                if taxa_list[i][j][k] == 'None':
                                    count_none += 1
                        sum_distance += (7-count_none)*5
                        count_none = 0
                sum_distance_list.append(sum_distance)
        return sum_distance_list

def normalization(tree_distance_list, sum_distance_list):
        for i in range(len(tree_distance_list)):
                for j in range(len(tree_distance_list[i])):
                        tree_distance_list[i][j] = tree_distance_list[i][j]/sum_distance_list[i]
        
        return tree_distance_list


def writeTreeDistance():
        sample_list = readSample(path)
        new_distance_list = normalization(tree_distance_list, sum_distance_list)
        for index in range(len(sample_list)):
                csvreader = csv.reader(open('../utils_data/taxons/sample_'+str(index+1)+'.csv'))
                matrix = []
                for i in csvreader:
                        matrix.append(i)
                del matrix[0]
                for i in range(len(matrix)):
                        matrix[i] = matrix[i][:-1]
                species_list = []

                for i in range(len(matrix)):
                        for j in range(i+1,len(matrix)):
                                list_ = matrix[i] + matrix[j]
                                if len(list_) != 14:
                                        print("Error!")
                                        break
                                species_list.append(list_)

                with open('./utils_data/species_compare/sample_'+str(index+1)+'_species.csv','w') as f:
                        for i in range(len(new_distance_list[index])):
                                for j in range(len(species_list[i])):
                                        f.write(species_list[i][j])
                                        f.write(',')
                                f.write(str(new_distance_list[index][i]))
                                f.write('\n')



# Version2: In each row, let the whole taxa name be written into one block
def writeTreeDistance_2():
        sample_list = readSample(path)
        for index in range(len(sample_list)):
                csvreader = csv.reader(open('./utils_data/taxons/sample_'+str(index+1)+'_new.csv'))
                matrix = []
                for i in csvreader:
                        matrix.append(i)
                for i in range(len(matrix)):
                        matrix[i] = matrix[i][:-2]

                data2 = []
                for i in range(len(matrix)):
                        for j in range(i+1,len(matrix)):
                                list_ = matrix[i] + matrix[j]
                                data2.append(list_)

                with open('./utils_data/species_compare/sample_'+str(index+1)+'_species_new.csv', 'w') as f:
                        for i in range(len(data2)):
                                for j in range(len(data2[i])):
                                        f.write(data2[i][j])
                                        f.write(',')
                                f.write('\n')

        for index in range(len(sample_list)):
                species_list2 = []
                distance_list2 = []
                with open('./utils_data/species_compare/sample_'+str(index+1)+'_species.csv') as f:
                        rows = csv.reader(f)
                        for i in rows:
                                distance_list2.append(i[-1])

                # Read the newest taxa name(format = k__XXX)        
                with open('./utils_data/species_compare/sample_'+str(index+1)+'_species_new.csv', 'r') as f:
                        rows2 = csv.reader(f)
                        for i in rows2:
                                species_list2.append(i[:-1])

                with open('./utils_data/species_compare/sample_'+str(index+1)+'_species_new.csv', 'w') as f:                
                        for j in range(len(distance_list2)):
                                for k in range(len(species_list2[j])):
                                        f.write(species_list2[j][k])
                                        f.write(',')
                                f.write(str(distance_list2[j]))                        
                                f.write('\n')

# Main
readSample(path)
getTaxaValue()
getTaxaValue_2()
taxa_list, tree_distance_list = getTaxaTreeDistance()
sum_distance_list = getSumTreeDistance(taxa_list)
writeTreeDistance()
writeTreeDistance_2()