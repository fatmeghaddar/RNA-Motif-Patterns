# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 14:31:03 2023
The following code is to extract the sequences from the FASTA Files downloaded from RNACentral website.
@author: Fatme Ghaddar
"""

infile = open('C:/Users/fatoo/Desktop/RNA RESEARCH/Complexity Paper/Subsamples_RNACentral_100_unified_NoDuplicates_NoJunk_uniqueRNA_NoSimilarities80.txt')
newopen2=open('Info_Subsamples_RNACentral_100_unified_NoDuplicates_NoJunk_uniqueRNA_NoSimilarities80.txt','w')
newopen=open('Seqs_Subsamples_RNACentral_100_unified_NoDuplicates_NoJunk_uniqueRNA_NoSimilarities80.txt','w')
for line in infile :
    #replace whats between the quotations with whatever is needed to be removed
    if '>'  not in line : 
        newopen.write(line)
    if '>' in line : 
        newopen2.write(line)
       
newopen.close()
newopen2.close()