#!/usr/bin/env python
# coding: utf-8

###################################################################################
# Software description: Script for simulating genetic data based on different genetic models
# Author: Anubhav Kaphle
# Date created : 28 AUG 2021
# Tested on Python 3.7.7, Numpy 1.18.5 and Pandas 1.0.4

"""
MIT License

Copyright (c) [2021] [Anubhav Kaphle]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
###################################################################################

import numpy as np
import pandas as pd
import sys
import argparse
import random
import time

#use this function for checking
def checkOR(Ds, Hn, Dh, Hs, OR_data):
    OR = Ds*Hn / (Dh*Hs)
    if np.isclose(OR, OR_data, rtol=0.01):
        print("All good")
    else:
        sys.exit(1)

#Hs is count of risk allele in non-disease state
#Hn is count of non-risk allele in non-disease state
#Ds is count of risk allele in disease state
#Dn is count of non-risk allele in disease state

def CreateGeno(Odd_ratio, Ncases, Ncontrol, maf, model):
    #Compute allele counts for cases and controls
    Hs  = np.rint(2*Ncontrol*maf) 
    #Hs  = np.ceil(2*Ncontrol*maf)
    Hn  = 2*Ncontrol - Hs
    Ds  = np.rint((2*Odd_ratio*Ncases*Hs) / (Hn+Odd_ratio*Hs)) 
    #Ds  = np.ceil((2*Odd_ratio*Ncases*Hs) / (Hn+Odd_ratio*Hs))
    Dn  = 2*Ncases - Ds
    
    case_risk_allele    = np.repeat(1,Ds)
    case_nonrisk_allele = np.repeat(0,Dn)
    cases_allele_pool   = np.concatenate((case_risk_allele, case_nonrisk_allele))
    random.shuffle(cases_allele_pool )
    index               = np.arange(cases_allele_pool.shape[0])
    random.shuffle(index)
    random.shuffle(index)
    index_half          = np.random.choice(index,int(cases_allele_pool.shape[0]/2) , replace=False)
    mask                = np.ones(cases_allele_pool.shape,dtype=bool)
    mask[index_half]    = False
    
    control_risk_allele     = np.repeat(1,Hs)
    control_nonrisk_allele  = np.repeat(0,Hn)
    control_allele_pool     = np.concatenate((control_risk_allele, control_nonrisk_allele))
    random.shuffle(control_allele_pool)
    
    index_2                 = np.arange(control_allele_pool.shape[0])
    random.shuffle(index_2)
    random.shuffle(index_2)
    index_half_2            = np.random.choice(index_2,int(control_allele_pool.shape[0]/2) , replace=False)
    mask_2                  = np.ones(control_allele_pool.shape,dtype=bool)
    mask_2[index_half_2]    = False
    
    sum_case     = cases_allele_pool[index_half] + cases_allele_pool[mask]
    sum_control  = control_allele_pool[index_half_2] + control_allele_pool[mask_2]
     
    #additive model assumes linear increase effect with linear increase in allele count, numerical {0,1,2} coding
    if model.lower() == "additive":     #coding risk allele as 1 and non-risk as 0
        return(np.concatenate((sum_case, sum_control)))

    elif model.lower() == "dominant":  #dominant allele is when even one presence will effect trait, binary coding here
        genotype_case    =  [1 if i in [1,2] else 0 for i in sum_case]
        genotype_control =  [1 if i in [1,2] else 0 for i in sum_control]
        return(np.concatenate((genotype_case,genotype_control)))
       
    elif model.lower() == "recessive": # recessive allele is effective when both are present, binary coding here as well
        genotype_case    =  [1 if i == 2 else 0 for i in sum_case]
        genotype_control =  [1 if i == 2 else 0 for i in sum_control]
        return(np.concatenate((genotype_case,genotype_control)))
    else:
        print("No other models possible. Please check this column in your original file")
        sys.exit(1)

#frequency of effect allele is maf. No. of cases and control in the GWAS for all SNPs is assumed to be constant. 
#This can be changed
def CreateData(pgs_file, outpath=None, Ncases=50000, Ncontrol=50000, model="additive"):
    data      = pd.read_csv(pgs_file,sep='\t')
    snp_id    = list(data.rsID)
    maf       = np.array(data.MAF)
    odd_ratio = np.array(data.OR)
    model     = np.array(data.Model)
    
    #flip allele if OR is <1 and change corresponding MAF by substracting from 1
    for i, j in enumerate(odd_ratio): 
        if j < 1:
            odd_ratio[i]=1/j
            maf[i]=1-maf[i]
        
    #diseased state is labelled 1 and non-diseased state as 0
    indiv_labels = np.concatenate((np.repeat(1,Ncases),np.repeat(0,Ncontrol)))
    indiv_ids    = np.array(["ID_{0}".format(i) for i in np.arange(Ncases+Ncontrol)])
    genotype_data = np.hstack((indiv_ids.reshape(indiv_ids.shape[0],1),indiv_labels.reshape(indiv_labels.shape[0],1)))
    for i in range(maf.shape[0]):
        genotype = CreateGeno(odd_ratio[i], Ncases, Ncontrol, maf[i], model[i])
        print("completed generating data for SNP {0} using genetic model = {1}".format(snp_id[i],model[i]))
        genotype_data=np.hstack((genotype_data,genotype.reshape(genotype.shape[0],1)))
        
        
    Genotype_dataframe = pd.DataFrame(genotype_data,columns=["ID","class"]+snp_id)
    
    outfile = outpath + "/" + "Simulated_genotype_data.txt"
    
    #If file size is too big, use compression in the arguement compression='gzip' or maybe use python HD5 filetype
    Genotype_dataframe.to_csv(outfile
         , sep='\t'
         , header=True
         , index=False
         , chunksize=100000
         , encoding='utf-8')
    
    print("Program ends here")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Software for creating synthetic genetic data based on PRS files and with different genetic models')
    parser.add_argument('--PRS_file', action='store', type=str,
                    help='PRS file', required=True)
    parser.add_argument('--outpath', action='store', type=str,
                    help='Please provide output file path without the last slash', required=True)

    parser.add_argument('--Ncases', action='store', type=int,
                    help='please provide no of cases in the GWAS. Assumed constant for all SNPs ')

    parser.add_argument('--Ncontrol', action='store', type=int,
                    help='please provide no of controls in the GWAS. Assumed constant for all SNPs ')

    args = parser.parse_args()
    
    pgs_file   = args.PRS_file
    outpath    = args.outpath
    Ncases     = args.Ncases
    Ncontrol   = args.Ncontrol
    time_start = time.time()
    CreateData(pgs_file, outpath, Ncases, Ncontrol)
    time_end   = time.time()
    print("\n\nTotal time taken to create the dataset was {0} minutes".format((time_end - time_start)/60))
