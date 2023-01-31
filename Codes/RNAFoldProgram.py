# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 14:32:07 2023
The following code is to fold the sequences of Natural and Random data into their corresponding Secondary Structures using Vienna Package.
In order to do so, Vienna Package is needed to be installed within the Python packages. Information on the Vienna Package can be  found in the Methods section of the paper,
@author: Dr. Kamaludin Dingle & Fatme Ghaddar
"""

import subprocess
import time
from subprocess import PIPE
import os
import re
startTime = int(time.time())

projectDirectory = os.path.dirname(os.path.abspath(__file__))
inputFilename='Random_L60_S100k_RNA_Generator'
pathToRNAfold = 'C:\Program Files (x86)\ViennaRNA Package\RNAfold.exe'

pathToInputFile = f'C:/Users/fatoo/Desktop/RNA RESEARCH/ViennaRNAProject/{inputFilename}.txt'

pathToOutputFile = f'{projectDirectory}\\SS_{inputFilename}.txt' 
outputStructures = open(pathToOutputFile, 'w')
readInput = open(pathToInputFile, 'r')

inputAsParam = ''.join(readInput.readlines())

runFold = subprocess.Popen(pathToRNAfold, stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True)
(fullOutput, originalInput) = runFold.communicate(inputAsParam)

output = fullOutput.split()#to add the data into list called output

for i in output:
    if (len(i)>10 and i[0]=='(' ): 
        outputStructures.write(i+'\n')
    elif (len(i)>10 and i[0]=='.' ):
        outputStructures.write(i+'\n')

outputStructures.close()

print('Execution time: ', int(time.time()) - startTime, 'sec')

