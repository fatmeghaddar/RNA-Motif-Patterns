# -*- coding: utf-8 -*-
"""
Created on Sun Aug 29 11:30:18 2021

@author: Fatme Ghaddar 
Random generator with probabilities of presence specified.
"""

import numpy as np
import os
import random
pathToOutputFile ='Random_L55_S500K_70GC_30AT_Percent1.txt' 


f=open(pathToOutputFile,'w')

for i in range(500000):
    
    f.write(''.join(np.random.choice(['A','T','C','G'],p=[0.15,0.15,0.35,0.35],size=55)))
    f.write('\n')

f.close()