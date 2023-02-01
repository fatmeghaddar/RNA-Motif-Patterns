
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 11:24:00 2020
The following code finds the confidence interval for Natural,random and experimental RNA samples by finding the intercept and slope values.
@author: Fatme Ghaddar
"""


import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import numpy.random as npr


randRNA=[]
NNrna=[]
expRNA=[]

randRNA_t=[]
NNrna_t=[]
expRNA_t=[]

randRNA_l=[]
NNrna_l=[]
expRNA_l=[]

Ylabel=['Number of Bulges','Number of Loops','Number of Junctions','Number of Helices','Number of Bonds']

# TO EXTRACT MOTIFS
import scipy.stats

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h



for line in open('Motifs_SS_RNACentral_50to3000_Filtered_Sizes_RandomSamples_20each_NoDuplicates_NoJunks.txt','r').readlines():
  l=line.split()
  NNrna.append([float(l[1]),float(l[2]),float(l[3]),float(l[4]),float(l[5]),float(l[6])])

for line in open('Motifs_SS_Random_L50upto3000_20each_RNA_Generator_20Feb2022.txt','r').readlines():
  l=line.split()
  randRNA.append([float(l[1]),float(l[2]),float(l[3]),float(l[4]),float(l[5]),float(l[6])])

for line in open('dp-11402-no-dimer_result_with_motifs.txt','r').readlines():
  l=line.split()
  expRNA.append([float(l[1]),float(l[2]),float(l[3]),float(l[4]),float(l[5]),float(l[6])])


#  MOTIFS NP ARRAY
randRNA=np.array(randRNA)
NNrna=np.array(NNrna)
expRNA=np.array(expRNA)


natFits=[]

n = len(NNrna)


for j in [1,2,3,4,5]:
    for i in range(1500):

        inds = np.random.randint(n,size=n)

        natFits.append(np.polyfit(NNrna[inds,0],  NNrna[inds,j] ,1 ))



        # Get the corresponding values of 2.5th and 97.5th percentiles
    conf_interval = np.percentile([row[0] for row in natFits],[5.0,95])

        # Print the interval
    print("The confidence interval of natural slope",Ylabel[j-1], ':',[round(x,4) for x in conf_interval])
          # Plot the PDF for bootstrap replicates as histogram
 #   plt.hist([row[1] for row in natFits],bins=10,normed=True)

 # Showing the related percentiles
  #  plt.axvline(x=np.percentile([row[1] for row in natFits],[5.0]),label='5th percentile',c='y')
   # plt.axvline(x=np.percentile([row[1] for row in natFits],[95]), label='95th percentile',c='r')

    #plt.xlabel("Intercept NatRNA")
    #plt.ylabel("Frequency")
    #plt.title("Histogram")
    #plt.legend()
    #plt.show()



        # Get the corresponding values of 2.5th and 97.5th percentiles
    conf_interval = np.percentile([row[1] for row in natFits],[5.0,95])

        # Print the interval
    print("The confidence interval of natural intercept",Ylabel[j-1], ':', [round(x,4) for x in conf_interval])


randFits=[]

n = len(randRNA)


for j in [1,2,3,4,5]:
    for i in range(1500):

        inds = np.random.randint(n,size=n)

        randFits.append(np.polyfit( randRNA[inds,0],  randRNA[inds,j] ,1 ))

      
  #  plt.hist([row[0] for row in randFits],bins=10,normed=True)

 # Showing the related percentiles
   # plt.axvline(x=np.percentile([row[0] for row in randFits],[5.0]),label='5th percentile',c='y')
    #plt.axvline(x=np.percentile([row[0] for row in randFits],[95]), label='95th percentile',c='r')

    #plt.xlabel("Slope RandRNA")
  #  plt.ylabel("Frequency")
   # plt.title("Histogram")
    #plt.legend()
    #plt.show()

        # Get the corresponding values of 2.5th and 97.5th percentiles
    conf_interval = np.percentile([row[0] for row in randFits],[5.0,95])

        # Print the interval
    print("The confidence interval of random slope",Ylabel[j-1], ':',[round(x,4) for x in conf_interval])
    
    #plt.hist([row[1] for row in randFits],bins=10,normed=True)

 # Showing the related percentiles
    #plt.axvline(x=np.percentile([row[1] for row in randFits],[5.0]),label='5th percentile',c='y')
    #plt.axvline(x=np.percentile([row[1] for row in randFits],[95]), label='95th percentile',c='r')

    #plt.xlabel("Intercept RandRNA")
    #plt.ylabel("Frequency")
    #plt.title("Histogram")
    #plt.legend()
    #plt.show()

        # Get the corresponding values of 2.5th and 97.5th percentiles
    conf_interval = np.percentile([row[1] for row in randFits],[5.0,95])

        # Print the interval
    print("The confidence interval of random intercept",Ylabel[j-1], ':',[round(x,4) for x in conf_interval])


expFits=[]

n = len(expRNA)


for j in [1,2,3,4,5]:
    for i in range(1500):

        inds = np.random.randint(n,size=n)

        expFits.append(np.polyfit( expRNA[inds,0],  expRNA[inds,j] ,1 ))


        # Get the corresponding values of 2.5th and 97.5th percentiles
    conf_interval = np.percentile([row[0] for row in expFits],[5.0,95])

        # Print the interval
    print("The confidence interval of exp slope",Ylabel[j-1], ':',[round(x,4) for x in conf_interval])

        # Get the corresponding values of 2.5th and 97.5th percentiles
    conf_interval = np.percentile([row[1] for row in expFits],[5.0,95])

        # Print the interval
    print("The confidence interval of exp intercept",Ylabel[j-1], ':',[round(x,4) for x in conf_interval])



