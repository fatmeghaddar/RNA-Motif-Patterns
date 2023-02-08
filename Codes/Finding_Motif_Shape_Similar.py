# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 10:36:28 2023
The following code uses Shannon's entropy to check if structures with exactly the same motifs show similar shapes
https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.entropy.html
@author: Fatme Ghaddar
"""

# import pandas library
import pandas as pd
import numpy as np
from collections import Counter
import csv
from scipy.stats import entropy
df = pd.read_csv('L100_Random_Motifs_Shapes.csv')



duplicate = df[df.duplicated(['Number_of_Bulges', 'Number_of_Loops','Number_of_Junctions','Number_of_Helices','Number_of_Bonds'])]


duplist= df.groupby(['Number_of_Bulges', 'Number_of_Loops','Number_of_Junctions','Number_of_Helices','Number_of_Bonds'])['Abstract shape Lev5'].apply(np.array).reset_index(name='Abstract shape Lev5')


dup1=[]
dup2=[]
fulldup=[]
allH=[]
all2H=[]
base = 2  # work in units of bits
for i in range (len(duplist['Abstract shape Lev5'])):

  #  if (len(duplist['Abstract shape Lev5'][i])>1):
    print(Counter(duplist['Abstract shape Lev5'][i]))
    
    dup1.append(Counter(duplist['Abstract shape Lev5'][i]))
    for j in Counter(duplist['Abstract shape Lev5'][i]).values():
        dup2.append(j/duplist['Abstract shape Lev5'].str.len()[i]) #getting the frequency of each shape over the total shapes of the structures having exactly the same motifs
        print(j, duplist['Abstract shape Lev5'].str.len()[i], dup2)

    fulldup.append(dup2)
    pk = np.array(dup2)  # each row 
    H = entropy(pk, base=base) #applying Shannon's entropy which gives H == -np.sum(pk * np.log(pk)) / np.log(base)
    
    
    allH.append(H)
    all2H.append(2**H)
    print(dup2)
    print(H)
    print(2**H)
    dup2=[]

    

#Adding all new values to the dataframe as columns

duplist['Shape counts']=dup1

duplist['Total Shapes'] = duplist['Abstract shape Lev5'].str.len()

duplist['Shannons Entropy H']=allH


duplist['2^H']=all2H


# sorting my dataframe

duplistsorted=duplist.sort_values(by='shape list length', ascending=False)

#print(duplistsorted)
duplistsorted.to_csv('L100_Random_Motifs_Shapes_ShannonsEntropy.csv', encoding='utf-8', index=False)


#after this the excel  sheet has been beautified manually









