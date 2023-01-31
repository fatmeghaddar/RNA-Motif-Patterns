# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 14:38:43 2023

Given an input file of RNA secondary structures (SS), this code returns the 
SS in addition to the number of:
bulges, loops, junctions, helices and bonds


Secstruc package is needed. 

@author: Dr. Kamaludin Dingle and Fatme Ghaddar
"""


import sys
import Secstruc

infilename='Seqs_RNACentral_85_unified_NoDuplicates_NoJunk.txt' 

sys.argv[0]
print ('File to analyse is ',infilename)    

outfilename= 'Motifs_'+infilename
Length=[]
Structures=[]
Bulges=[]
Loops=[]
Junctions=[]
Helices=[]
Bonds=[]
c=0
# Now start to count motifs
print ('Counting motifs...')
for struc in open(infilename, 'r'):
    
    
    SS=struc.strip()
    Structures.append(SS)
    print (SS)
    Length.append(len(SS))
    
    Bulges.append(len(Secstruc.Secstruc(SS).find_bulges()))
    print (len(Secstruc.Secstruc(SS).find_bulges()))
    
    Loops.append(len(Secstruc.Secstruc(SS).find_loops()))
    print (len(Secstruc.Secstruc(SS).find_loops()))
    
    Junctions.append(len(Secstruc.Secstruc(SS).find_junctions()))
    print (len(Secstruc.Secstruc(SS).find_junctions()))
    
    Helices.append(len(Secstruc.Secstruc(SS).find_helices()))
    print (len(Secstruc.Secstruc(SS).find_helices()))
    
    Bonds.append(SS.count('('))
    c=c+1
    print (c)
    if (c%2)==0:
        print (100.0*c/1000,'% completed')

# Now print to file
outputfile=open(outfilename,'a')
#outputfile.write(' \t SS  Length   Bulges     Loops     Junctions     Helices     Bonds \n') 
for indx in range(len(Structures)):
 
    outputfile.write(Structures[indx]+'\t'+str(Length[indx])+'\t'+str(Bulges[indx])+'\t'+str(Loops[indx])+'\t'+str(Junctions[indx])+'\t'+str(Helices[indx])+'\t'+str(Bonds[indx])+'\n')    

    

print ('SS, number of bulges, loops, junctions, helices and bonds in \n',outfilename)   
print ('>>> Done <<<') 
outputfile.close()    

    
