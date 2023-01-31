
"""
Created on Tue Jan 31 14:14:59 2023
The following code is used in Figures 2 and 3.

@author: Dr.Kamaludin Dingle & Fatme Ghaddar
"""



import numpy as np
import matplotlib.pyplot as plt
import collections
import scipy.stats as stats
import os
import pandas as pd



projectDirectory = os.path.dirname(os.path.abspath(__file__))

FIGPATH=projectDirectory

params = [(100,5),
     
          (200,5),
     
          (300,5),
          (400,5)]



for L, Level in params:
    print ('\nL = ',L)
    print ('Level=',Level)

    # Print estimate of number of shapes, from Nebel nebel2009quantitative (Table 1)
    if Level==3:
        est= 1.85* 1.46**L * L**(-1.5)
        print ('\ns_3_'+str(L)+'=log10 ',np.log10(est))
    else:
        est= 2.44* 1.32**L * L**(-1.5)
        print ('\ns_5_'+str(L)+'=log10',np.log10(est))


    # Show rank plot with natural data
    ShowRankPlot = 1
    
    # RNACentral (database natural data)
    FILEPATH_prefix = projectDirectory

    #Sampled_Shapes
    FILEPATH_rand = projectDirectory #to be changed according to the user's path.
    
    if L==100:
        INFILE = FILEPATH_rand + '/Abstract_Random_L100_Lev5_S30k__RNA_Generator.txt'
        Nat_INFILE = FILEPATH_prefix+'/Abstract_L100_Lev5_May.txt'
        
    elif L==200:
        Nat_INFILE = FILEPATH_prefix+'/Abstract_L200_Lev5_March2022.txt'
        INFILE = FILEPATH_rand + '/Abstract_Random_L200_Lev5_S30k__RNA_Generator.txt'
    elif L==300:
        Nat_INFILE = FILEPATH_prefix+'/Abstract_L300_Lev5_March2022.txt'
        #INFILE = FILEPATH_rand + '/Abstract_Random_L'+str(L)+'_Lev'+str(Level)+'_S100k__RNA_Generator.txt'
        INFILE = FILEPATH_rand + '/Abstract_Random_L300_Lev5_S30k__RNA_Generator.txt'
    elif L==400:
        Nat_INFILE = FILEPATH_prefix+'/Abstract_L400_Lev5_March2022.txt'
        #INFILE = FILEPATH_rand + '/Abstract_Random_L'+str(L)+'_Lev'+str(Level)+'_S100k__RNA_Generator.txt'
        INFILE = FILEPATH_rand + '/Abstract_Random_L400_Lev5_S30k__RNA_Generator.txt'
    else:
        print ('Error - no data file provided')


    # Abstract SS in text format
    AbSS = np.genfromtxt(INFILE,dtype=str)# sampled
    Nat_AbSS = np.genfromtxt(Nat_INFILE,dtype=str) # natural data

    c = collections.Counter(AbSS)
    Nat_c = collections.Counter(Nat_AbSS)

    # Outputs and their frequencies
    print ('\nNumber of random samples',len(AbSS))
    print ('\nNumber of natural data',len(Nat_AbSS))
    print ('\nNumber of unique random shapes found = ',len(list(c)))
    print ('\nNumber of unique natural shapes found = ',len(list(Nat_c)))

    # Freqencies
    Freqs = list(c.values())# freqs of each random sampled abstract shape
    Nat_Freqs = list(Nat_c.values()) # freqs of each natural abstract shape
    SortedFreqs = np.sort(Freqs)
    P_rand = 1.0*np.array(Freqs) / sum(Freqs)
    Nat_SortedFreqs = np.sort(Nat_Freqs)
    rand_shape_rank = np.argsort(Freqs)[::-1]
    #P_rand_ranked = P_rand[np.argsort(Freqs)[::-1]]

    # Show rank plot and bias in shape frequencies
    if ShowRankPlot ==1:
        plt.figure()
        therank = np.arange(len(Freqs))+1 # add one to make the rank start at 1, like Nat_rank
        plt.semilogy(therank,1.0*SortedFreqs[::-1]/sum(SortedFreqs[::-1]),'.',ms = 7,mew=1.5,label='Sampled')
        # Find the rank of the natural shapes
        Nat_rank=[len(np.where(np.asarray(list(c.values()))>= c[shp])[0]) for shp in set(Nat_AbSS)] # shp = an abstract shape
        Nat_rank_prob=[1.0*c[shp]/sum(c.values()) for shp in set(Nat_AbSS)]
        plt.scatter(Nat_rank,Nat_rank_prob,facecolors='none',edgecolors='orange',linewidth=3,s=130,label='Natural')
        
        # font sizes
        plt.ylabel(r'$f_{p} ^G$',fontsize=22)
        plt.xlabel('Rank',fontsize=22)
        plt.title("L = "+str(L)+'   Level = '+str(Level),fontsize=18)
        plt.tick_params(axis="x", labelsize=16)
        plt.tick_params(axis="y", labelsize=16)
        plt.legend(fontsize=16)
        # axis limits
        y_min_log10 = np.floor(np.log10(min(1.0*SortedFreqs[::-1]/sum(SortedFreqs[::-1]))))
        plt.ylim(bottom = 10**y_min_log10)
     
        plt.tight_layout()
        plt.show()
        plt.savefig(FIGPATH+'/Rank_Nat_Samp_L'+str(L)+'_Lev'+str(Level)+'.pdf',dpi=300)


    Shapes = list(c) # list of *unique* randomly sampled shapes
    Nat_Shapes = list(Nat_c) # list of *unique* natural shapes from data

    # Now compare frequencies of shapes
    CombinedShapes = []
    Count = []
    Nat_Count = []
    Nat_not_found_by_sampling =[] ######################################## !
    for shape in set(Shapes+Nat_Shapes):# ie collection of all shapes found from sampling AND from nature
        CombinedShapes.append(shape)
        Count.append(c[shape])
        Nat_Count.append(Nat_c[shape])
        if c[shape] == 0:
            Nat_not_found_by_sampling.append(shape)

    print ('\nFraction of unique natural shapes found by random sampling =',len(list(Nat_c))-len(Nat_not_found_by_sampling),'/',len(list(Nat_c)))


    # Probability of shapes
    P = 1.0*np.asarray(Count)/sum(Count)
    Nat_P = 1.0*np.asarray(Nat_Count)/sum(Nat_Count)

    # plot P(x) for natural vs P(x) sampled
    plt.figure()
    plt.scatter(Nat_P[Nat_P>0],P[Nat_P>0],facecolors='none',linewidth=3,s=100,edgecolors='orange',label="L="+str(L)+'   '+'Level='+str(Level))
    # x=y line/prediction
    min_xy = min([min(P[Nat_P>0]),min(Nat_P[Nat_P>0])])
    

    max_xy = max([max(P[Nat_P>0]),max(Nat_P[Nat_P>0])])
    plt.loglog([min_xy,max_xy],[min_xy,max_xy],'-g',label='x = y',alpha=0.6)
    # font sizes
    plt.tick_params(axis="x", labelsize=16)
    plt.tick_params(axis="y", labelsize=16)
    plt.ylabel(r'$f_{p}^G$',fontsize=25)# natural
    plt.xlabel(r'$f_p$',fontsize=25)
    plt.title("L = "+str(L)+'   Level = '+str(Level),fontsize=18)
    #plt.legend(fontsize=16) # no need, because we have title
    plt.tight_layout()
    plt.show()


    # this is the fourth shape from the end of the list

    plt.savefig(FIGPATH+'/PNat_vs_PSamp_L'+str(L)+'_Lev'+str(Level)+'.pdf',dpi=300)

    # Correlation of natural and random probs
    P_nozerosPNatP = np.asarray([P[j] for j in range(len(P)) if P[j]*Nat_P[j]>0]) # ie only if a shape appears in rand and natural at least once
    Nat_P_nozerosPNatP = np.asarray([Nat_P[j] for j in range(len(Nat_P)) if P[j]*Nat_P[j]>0])
    print ('\nCorrelation coeff logP vs logNat_P EXCLUDING zero freqs',
           round(stats.pearsonr(
                                np.log10(P_nozerosPNatP),np.log10(Nat_P_nozerosPNatP)
                                )[0],2), ' p-value ',stats.pearsonr(
                                np.log10(P_nozerosPNatP),np.log10(Nat_P_nozerosPNatP))[1])
