# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 10:07:34 2023
The following code is for figure 7 in the paper, that finds the fits for  natural and randomly sampled RNA, with GC adjusted content RNA.
@author: Fatme Ghaddar
"""


import matplotlib.pyplot as plt
import numpy as np

randRNA=[]
NNrna=[]
expRNA=[]
scrambledRNA=[]
Ylabel=['Number of Bulges','Number of Loops','Number of Junctions','Number of Helices','Number of Bonds']

for line in open('Motifs_SS_RNACentral_50to3000_Filtered_Sizes_RandomSamples_20each_NoDuplicates_NoJunks.txt','r').readlines():
  l=line.split()
  NNrna.append([float(l[1]),float(l[2]),float(l[3]),float(l[4]),float(l[5]),float(l[6])])

for line in open('Motifs_SS_Random_L50upto3000_20each_RNA_Generator_20Feb2022.txt','r').readlines():
  l=line.split()
  randRNA.append([float(l[1]),float(l[2]),float(l[3]),float(l[4]),float(l[5]),float(l[6])])

for line in open('dp-11402-no-dimer_result_with_motifs.txt','r').readlines():
  l=line.split()
  expRNA.append([float(l[1]),float(l[2]),float(l[3]),float(l[4]),float(l[5]),float(l[6])])
for line in open('Motifs_Randomized_RNACentral_50to1500_Filtered_Sizes_RandomSamples_20each_NoDuplicates_NoJunks.txt','r').readlines():
  l=line.split()
  scrambledRNA.append([float(l[1]),float(l[2]),float(l[3]),float(l[4]),float(l[5]),float(l[6])])
randRNA=np.array(randRNA)
NNrna=np.array(NNrna)
scrambledRNA=np.array(scrambledRNA)
expRNA=np.array(expRNA)

for i in [1,2,3,4,5]: #for moving around lengths,loop,junctions,helices and bond
    L=np.arange(1,3000) #fixed length size
  #  plt.xlim([0,1500])
    plt.ylabel(Ylabel[i-1],fontsize=15)
    plt.xlabel('Length',fontsize=15)
    fit_i_N=np.polyfit(NNrna[:,0],NNrna[:,i],1)
    plt.plot(L,fit_i_N[0]*L+fit_i_N[1],c='#ff7f00',label='Natural RNA fit') # work in fit and to find std following
    fit_i_R=np.polyfit(randRNA[:,0],randRNA[:,i],1)
    plt.plot(L,fit_i_R[0]*L+fit_i_R[1],c='#0868ac',label='Sampled RNA fit')
    fit_i_lE=np.polyfit(expRNA[:,0],expRNA[:,i],1)
 #   plt.plot(L,fit_i_lE[0]*L+fit_i_lE[1],c='b',label='Experimental RNA')
    fit_i_lScram=np.polyfit(scrambledRNA[:,0],scrambledRNA[:,i],1)
    plt.plot(L,fit_i_lScram[0]*L+fit_i_lScram[1],c='y',label='GC Adjusted RNA')
    plt.rcParams.update({'font.size': 14})
    
    plt.xlim([0,3000])
 #   plt.legend(bbox_to_anchor=(1.9, 1.05)) #place diagrams in correct place
    plt.legend()   
    plt.tight_layout()
    plt.show()


