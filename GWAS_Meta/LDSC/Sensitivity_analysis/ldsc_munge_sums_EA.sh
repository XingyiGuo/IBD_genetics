#!/bin/bash

#SBATCH --job-name=ldsc.eas
#SBATCH --error=%x-%j.error
#SBATCH --out=%x-%j.out
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-0:0:0

source activate ldsc

# For overall IBD:
/home/lyul1/ldsc/munge_sumstats.py \
--sumstats /nobackup/sbcs/lyul1/IBD/metal_update/withoutUKB/gc_0.01/GC_Final/processed/EASEUR_CD_META1.TBL.SNP2.txt \
--snp SNP \
--N 25457 \
--chunksize 500000 \
--a1 Allele2 \
--a2 Allele1 \
--out eas_cd \
--merge-alleles /home/lyul1/LDSCref/w_hm3.snplist

/home/lyul1/ldsc/munge_sumstats.py \
--sumstats /nobackup/sbcs/lyul1/IBD/metal_update/withoutUKB/gc_0.01/GC_Final/processed/EASEUR_UC_META1.TBL.SNP2.txt \
--snp SNP \
--N 31677 \
--chunksize 500000 \
--a1 Allele2 \
--a2 Allele1 \
--out eas_uc \
--merge-alleles /home/lyul1/LDSCref/w_hm3.snplist



/home/lyul1/ldsc/munge_sumstats.py \
--sumstats /nobackup/sbcs/lyul1/IBD/metal_update/withoutUKB/gc_0.01/GC_Final/processed/EASEUR_IBD_META1.TBL.SNP2.txt \
--snp SNP \
--N 58232 \
--chunksize 500000 \
--a1 Allele2 \
--a2 Allele1 \
--out eas_ibd \
--merge-alleles /home/lyul1/LDSCref/w_hm3.snplist

