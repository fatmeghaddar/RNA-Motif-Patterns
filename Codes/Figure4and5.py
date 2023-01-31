# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 14:26:24 2023

@author: Fatme Ghaddar
"""


import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd



randRNA=[]
NNrna=[]
expRNA=[]
NNrnaL=[]


# Program to find most frequent  
# element in a list 
  
def most_frequent(List): 
    counter = 0
    num = List[0] 
      
    for i in List: 
        curr_frequency = List.count(i) 
        if(curr_frequency> counter): 
            counter = curr_frequency 
            num = i 
  
    return num 
  
Ylabel=['Number of Bulges','Number of Loops','Number of Junctions','Number of Helices','Number of Bonds']

# TO EXTRACT MOTIFS



for line in open('Motifs_SS_RNACentral_50to3000_Filtered_Sizes_RandomSamples_20each_NoDuplicates_NoJunks.txt','r').readlines():
  l=line.split()
  
  NNrna.append([float(l[1]),float(l[2]),float(l[3]),float(l[4]),float(l[5]),float(l[6])])
  NNrnaL.append(float(l[1])) #to find the most frequent natural length

for line in open('Motifs_SS_Random_L50upto3000_20each_RNA_Generator_20Feb2022.txt','r').readlines():
    l=line.split()
    randRNA.append([float(l[1]),float(l[2]),float(l[3]),float(l[4]),float(l[5]),float(l[6])])

np.set_printoptions(precision=8,suppress=True)
randRNA=np.array(randRNA)
NNrna=np.array(NNrna)

#Previous lines to get data and store in array for random and natural
L=np.arange(1,3000) #fixed length size
for i in [1,2,3,4,5]: #for moving around motifs
    plt.scatter(NNrna[:,0],NNrna[:,i],c='#ff7f00',s=3,alpha=0.5,label='Natural RNA')
    plt.scatter(randRNA[:,0], randRNA[:,i],s=3, alpha=0.5,c='#0868ac', label='Sampled RNA')

    fit_i_n=np.polyfit(NNrna[:,0],NNrna[:,i],1)
    plt.plot(L,fit_i_n[0]*L+fit_i_n[1],c='#ff7f00',label='Natural RNA Fit') # work in fit and to find std following

    fit_i_r=np.polyfit(randRNA[:,0],randRNA[:,i],1)
    plt.plot(L,fit_i_r[0]*L+fit_i_r[1],c='#0868ac',label='Sampled RNA Fit') # work in fit and to find std following
    print(i)
    plt.rcParams.update({'font.size': 14})
    plt.ylabel(Ylabel[i-1],fontsize=15)
    plt.xlabel('Length',fontsize=15)
    plt.legend() 
    #plt.legend(bbox_to_anchor=(1.1, 1.05)) #place diagrams in correct place
    plt.tight_layout()
    plt.show()



KG=0.074*L-0.377
KP=0.177*L-0.443
mKP=0.1717*L
mean_number_bonds =  0.24*L - 0.7

plt.scatter(L,KG,c='#e41a1c',s=0.05,label='Sampled $K_G$')
plt.scatter(L,KP,c='#a65628',s=0.05,label='Sampled $K_P$')
plt.scatter(L,mKP,c='#999999',s=0.05,label='Analytic $K_P$')
plt.scatter(NNrna[:,0],NNrna[:,4],c='#ff7f00',s=3,alpha=0.5,label='Natural RNA')
plt.scatter(randRNA[:,0], randRNA[:,4],s=3, alpha=0.5,c='#0868ac', label='Sampled RNA')
fit_i_n=np.polyfit(NNrna[:,0],NNrna[:,4],1)
plt.plot(L,fit_i_n[0]*L+fit_i_n[1],c='#ff7f00',label='Natural RNA Fit') # work in fit and to find std following

fit_i_r=np.polyfit(randRNA[:,0],randRNA[:,4],1)
plt.plot(L,fit_i_r[0]*L+fit_i_r[1],c='#0868ac',label='Sampled RNA Fit') # work in fit and to find std following

plt.xlabel('Length')
plt.ylabel('Number of Helices')
plt.legend()

plt.show()

plt.scatter(NNrna[:,0],NNrna[:,5],c='#ff7f00',s=3,alpha=0.5,label='Natural RNA')
plt.scatter(randRNA[:,0], randRNA[:,5],s=3, alpha=0.5,c='#0868ac', label='Sampled RNA')
plt.ylabel('Number of Bonds')
plt.xlabel('Length')
fit_l_bonds=np.polyfit(randRNA[:,0],randRNA[:,5],1)
plt.plot(L,fit_l_bonds[0]*L+fit_l_bonds[1],c='#0868ac')
bonds=randRNA[:,5]
val1=[]
for i in range (len(bonds)):
    v1=bonds[i]-(fit_l_bonds[0]*randRNA[5,0]+fit_l_bonds[1])
    val1.append(v1)
    meanB=np.mean(val1)
    stdB=np.std(val1)
fit_i_n=np.polyfit(NNrna[:,0],NNrna[:,5],1)
plt.plot(L,fit_i_n[0]*L+fit_i_n[1],c='#ff7f00',label='Natural RNA Fit') # work in fit and to find std following
       
plt.scatter(L,mean_number_bonds,c='lightpink',s=0.05,label='Mean Number of Bonds')
plt.legend()
plt.show()

    
