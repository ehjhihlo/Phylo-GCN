# -*- coding: utf-8 -*-
"""CRC_SVM

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/15b4feLQOj2nnV-sD-kl315oCeuM3TX1V
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import seaborn as sns
from sklearn import datasets
from sklearn import metrics
from sklearn.preprocessing import LabelBinarizer
from sklearn.model_selection import train_test_split
from sklearn import svm
from sklearn.metrics import confusion_matrix, roc_auc_score, f1_score, precision_score, recall_score
os.environ['TF_CPP_MIN_LOG_LEVEL']= '2'

# !gdown --id '1S9iwczSf6KL5jMSmU20SXKCSD3BUx4o_' --output level-6.csv #GMPR_genus
!gdown --id '1q0yp1iM66BKvqee46bOuSZYwl_SJCTp0' --output level-6.csv #GMPR_species

train = pd.read_csv("level-6.csv")
train.head()
train.info()

labelencoder = LabelEncoder()
train["Diagnosis"] = labelencoder.fit_transform(train["Diagnosis"])
# test["Diagnosis"] = labelencoder.fit_transform(test["Diagnosis"])
# for i in range(len(train)):
#     if train["Diagnosis"][i] == 'Cancer':
#         train["Diagnosis"][i] = str(1)
#     else:
#         train["Diagnosis"][i] = str(0)
train

not_select = ["index", "Diagnosis"]
train_select = train.drop(not_select,axis=1)
df_final_select = train_select

"""#SVM"""

#Use SVM to predict Cancer
x = df_final_select
y = train["Diagnosis"]
# y = np.array(y,dtype=int)
X_train,X_test,y_train,y_test = train_test_split(x,y,test_size=0.2,random_state=0)

clf = svm.SVC()
clf.fit(X_train,y_train)
y_predict = clf.predict(X_test)
score_clf = clf.score(X_test,y_test)
score_clf_train = clf.score(X_train,y_train)
print("train_accuracy = ",score_clf_train*100," %")
print("val_accuracy = ",score_clf*100," %")

mat = confusion_matrix(y_test, y_predict)
sns.heatmap(mat.T, square=True, annot=True, fmt='d', cbar=False)
plt.xlabel('true label')
plt.ylabel('predicted label')
score_recall = recall_score(y_test, y_predict, average=None)
f1score = f1_score(y_test, y_predict, average="macro")
precisionscore = precision_score(y_test, y_predict, average=None)
auc_roc = roc_auc_score(y_test, y_predict)
print("precision = ",precisionscore)
print("recall = ",score_recall)
print("auc_roc = ",auc_roc)
print("f1_score = ",f1score)

with open('SVM_result.csv','w') as f:
    f.write('Precision_Normal,Precision_Cancer,Recall_Normal,Recall_Cancer,Auc_Score,F1_Score,')
    f.write('\n')
    f.write(str(precisionscore[0])+','+str(precisionscore[1])+','+str(score_recall[0])+','+str(score_recall[1])+','+str(auc_roc)+','+str(f1score))
