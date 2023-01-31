# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 14:14:59 2023
The following code is used in Figure 6

@author: Dr. Kamaludin Dingle & Fatme Ghaddar
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from scipy import stats
from sklearn import metrics
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_auc_score

# This code is for analysing RNA structures, and computing classification accuracy based on motifs


def Object2Float(Y):
    Y_float = np.zeros((Y.shape[0],Y.shape[1]))
    for i in range(Y.shape[0]):
        for j in range(Y.shape[1]):
            Y_float[i,j] = np.float(Y[i,j])
    return Y_float


def BootstrapROCAUC(y_true,y_pred):
    print("Original ROC area: {:0.3f}".format(roc_auc_score(y_true, y_pred)))

    n_bootstraps = 1000
   
    bootstrapped_scores = []

    for i in range(n_bootstraps):
        # bootstrap by sampling with replacement on the prediction indices
        indices = np.random.randint(0, len(y_pred) - 1, len(y_pred))
        if len(np.unique(y_true[indices])) < 2:
            # We need at least one positive and one negative sample for ROC AUC
            # to be defined: reject the sample
            continue

        score = roc_auc_score(y_true[indices], y_pred[indices])
        bootstrapped_scores.append(score)

    print('95% confid interval is {:0.2f}'.format(np.percentile(bootstrapped_scores,2.5)),' to {:0.2f}'.format(np.percentile(bootstrapped_scores,97.5)))
    return bootstrapped_scores

####################################################
####################################################
PATH_FIGS = ''

# Choose which length

L = 1000
GC=True
if L==55:
    excel_file = pd.ExcelFile('L55_RFAMnatS246_RandomS1K.xlsx')
    print ('\n+++ Using L=55 ')
elif L==85:
    excel_file = pd.ExcelFile('L85_RFAMnatS497_RandomS1K.xlsx')
    print ('\n+++ Using L=85 ')
elif L==1000 and GC==False:
    excel_file = pd.ExcelFile('L1000_RNACentral_1K_Random_NoJunkNoDuplicates-Copy.xlsx')
    print ('\n+++ Using L=1000 ') #Copy of excel is made bc it wasnt allowing me to run on the old one bc it was open in python
    
elif L==1000 and GC==True:
    excel_file = pd.ExcelFile('L1000_RNACentral_Randomized_NoJunkNoDuplicates_May2022_Updated.xlsx')
    print ('\n+++ Using L=1000 GC Content') #Copy of excel is made bc it wasnt allowing me to run on the old one bc it was open in python

elif L==100:
    excel_file = pd.ExcelFile('L100_RNACentral_30K_Random_NoJunkNoDuplicates.xlsx')
    print ('\n+++ Using L=100 ')
elif L==400:
    
    excel_file = pd.ExcelFile('L400_RNACentral_Random_NoJunkNoDuplicates_31August2022_copyy.xlsx')
    print ('\n+++ Using L=400 ')
if True:
    # Data file
    print ('\n+++ Analysing Normal vs Low BMD +++')
    group1 = 'Natural'
    group2  = 'Random Samples'
    dropcols =[]

if GC:
    # Data file
    print ('\n+++ Analysing Normal vs Low BMD +++')
    group1 = 'Natural'
    group2  = 'Scrambled Natural Samples'
    dropcols =[]

# table
table_raw_all = excel_file.parse()

# drop some columns
table = table_raw_all.drop(dropcols,axis=1)# drop some columns


# count the number of missing values in each cytokine

X = table.values[:,1:]

# make it a float, not an object
X = Object2Float(X)

# scale the data to mean zero and var 1
scaler = preprocessing.StandardScaler().fit(X)
X_scaled = scaler.transform(X)


# the cytokines used in the study are:
cyt_list = [str(c) for c in table.columns[1:]]
print ('\nMotifs used =',cyt_list)
print ('\nNumber of samples=',X_scaled.shape[0])
print ('\nNumber of variables=',X_scaled.shape[1])


# Separate the data out
x1 = X_scaled[table['RNA_Type']=='Natural']
x2 = X_scaled[table['RNA_Type']=='Random ']



#grouplabels=1*(table['Grp (N,L)']=='Normal BMD')
grouplabels = 1*(table['RNA_Type']=='Natural ')

grouplabels = grouplabels.values# dummy variables '1' for L (or OSR) and '0' for N (or OSN).
# In other words, "1" is the more extreme group, and "0" is the less extreme group




#############################################
############  PLSDA Projection   #############
#############################################
# Now do PLS-DA to see in what ways the groups differ, and how strong the diff is
# "Partial least squares discriminant analysis (PLS-DA) is a variant used when the Y is CATEGORICAL." wikipedia
# the weights describe the contribution of each variable to each LV. from https://www.mfitzp.com/article/partial-least-squares-discriminant-analysis-plsda/
from sklearn.cross_decomposition import PLSRegression
plsr = PLSRegression(n_components=2)
plsr.fit(X_scaled,grouplabels)


# Plot PLS-DA 2D plot
plt.figure()
plt.scatter(plsr.x_scores_[grouplabels==0,0],plsr.x_scores_[grouplabels==0,1],c='blue',marker='+',s=80,label=group1,alpha=0.5)# group1, dummy var=0
plt.scatter(plsr.x_scores_[grouplabels==1,0],plsr.x_scores_[grouplabels==1,1],c='orange',marker='o',s=40,label=group2,alpha=0.1)# group2, dummy var=1
plt.xlabel('LV1',fontsize=15)
plt.ylabel('LV2',fontsize=15)
plt.title('L='+str(L),fontsize=15)
plt.legend()
plt.rcParams.update({'font.size': 10})
plt.show()

SAVEIT = True 
if SAVEIT:
    plt.savefig(PATH_FIGS+group1+'_'+group2+'_PLSDA_projection_L+'+str(L)+'.pdf',dpi=300)


#############################################################
###### Variable importances from PLS regression/PLS-DA ######
#############################################################
print ('\nVariable importances (PLSR) ')
VIP_plsr = []
for cyt in cyt_list:
    s = plsr.coef_[cyt_list.index(cyt)]# plsr.coef_ is the coef of linear model Y = X coef_ + Err
    s = np.round(s[0],2)
    VIP_plsr.append(s)
    print (cyt,s)
VIP_plsr = np.array(VIP_plsr)# variable importances, as judged by PLS regression

#  Plot the var imp for PLS regression/PLSDA
#plt.xticks(rotation=90)
plt.figure()
cl_order_plsr = [cyt_list[j] for j in np.argsort(VIP_plsr) if cyt_list[j]!='Length']

VIP_plsr=np.delete(VIP_plsr,0)

plt.bar(x=range(len(cl_order_plsr)),tick_label=cl_order_plsr,height=VIP_plsr[np.argsort(VIP_plsr)],color='red')
plt.ylabel('Variable importances',fontsize=10)
if GC==False:
    plt.title(f" {L} Natural vs. Ramdom Samples")
else:
    plt.title(f" {L} Natural vs. Scrambled Natural")
plt.xticks(fontsize=10, rotation=20, ha='right')
plt.rcParams.update({'font.size': 10})
#plt.title('PLS-DA')
plt.grid('on')
plt.show()

if SAVEIT:
    plt.savefig(PATH_FIGS+group1+'_'+group2+'_PLSDA_VIP_L+'+str(L)+'.pdf',dpi=300)






########################################################################
# Classification performance/ability (ROC AUC) using leave one out
#########################################################################
loo = LeaveOneOut()
group_preds = []
for train, test in loo.split(X_scaled):
    plsr.fit(X_scaled[train],grouplabels[train]) #knearest neighbors  neigh.fit(X, y)
    group_preds.append(plsr.predict(X_scaled[test])[0][0])

# Compare prediction accuracy
rocauc = metrics.roc_auc_score(grouplabels,group_preds)
print ('\nROC AUC for classification ability',np.round(rocauc,3))
print('\n')



########################################################################
# Bootstrap confidence interval - Classification performance/ability (ROC AUC) using leave one out
#########################################################################
BootstrapROCAUC(grouplabels,np.array(group_preds))

excel_file.close()