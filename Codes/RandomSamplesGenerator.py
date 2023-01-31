# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 14:35:37 2023
The following code is to get random samples based on one length
@author: Fatme Ghaddar
"""


import numpy as np
import os
import random
pathToOutputFile ='Random_L60_S100k_RNA_Generator.txt' 


f=open(pathToOutputFile,'w')

for i in range(100000):
    
    f.write(''.join(np.random.choice(['A','G','T','C'],size=60)))
    #f.write(''.join(np.random.choice(['A','G','T','C'],size=random.randint(100,3500)))) #TO GET A RANDOM LENGTH
    f.write('\n')

f.close()

