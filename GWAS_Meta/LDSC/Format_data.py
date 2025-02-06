#!/usr/bin/env python
# coding: utf-8

# In[1]:

#### formatting data for ldsc - remove MarkerNames col and NAs in SNP col #####


import os
import pandas as pd

directory = "/nobackup/sbcs/lyul1/IBD/metal_update/withUKB/gc_0.01/GC_Final/processed"  # Replace with your desired directory
os.chdir(directory)


def format_gwasmeta(input, output):

        df = pd.read_csv(input, sep="\t")
        
        df = df.drop(columns=['MarkerName'], errors='ignore')  
        
        df = df.dropna(subset=['SNP'])
        
        df.to_csv(output, sep="\t", index=False)


files = [
    ("EASEUR_CD_META1.TBL.SNP.txt", "EASEUR_CD_META1.TBL.SNP2.txt"),
    ("EASEUR_UC_META1.TBL.SNP.txt", "EASEUR_UC_META1.TBL.SNP2.txt"),
    ("EASEUR_IBD_META1.TBL.SNP.txt", "EASEUR_IBD_META1.TBL.SNP2.txt"),
    ("EUR_CD_META1.TBL.SNP.txt", "EUR_CD_META1.TBL.SNP2.txt"),
    ("EUR_UC_META1.TBL.SNP.txt", "EUR_UC_META1.TBL.SNP2.txt"),
    ("EUR_IBD_META1.TBL.SNP.txt", "EUR_IBD_META1.TBL.SNP2.txt")
]

for input, output in files:
    format_gwasmeta(input, output)


# In[2]:

directory = "/nobackup/sbcs/lyul1/IBD/metal_update/withoutUKB/gc_0.01/GC_Final/processed"  # Replace with your desired directory
os.chdir(directory)

def format_gwasmeta(input, output):

        df = pd.read_csv(input, sep="\t")
        
        df = df.drop(columns=['MarkerName'], errors='ignore')  
        
        df = df.dropna(subset=['SNP'])
        
        df.to_csv(output, sep="\t", index=False)


files = [
    ("EASEUR_CD_META1.TBL.SNP.txt", "EASEUR_CD_META1.TBL.SNP2.txt"),
    ("EASEUR_UC_META1.TBL.SNP.txt", "EASEUR_UC_META1.TBL.SNP2.txt"),
    ("EASEUR_IBD_META1.TBL.SNP.txt", "EASEUR_IBD_META1.TBL.SNP2.txt"),
    ("EUR_CD_META1.TBL.SNP.txt", "EUR_CD_META1.TBL.SNP2.txt"),
    ("EUR_UC_META1.TBL.SNP.txt", "EUR_UC_META1.TBL.SNP2.txt"),
    ("EUR_IBD_META1.TBL.SNP.txt", "EUR_IBD_META1.TBL.SNP2.txt")
]

for input, output in files:
    format_gwasmeta(input, output)
