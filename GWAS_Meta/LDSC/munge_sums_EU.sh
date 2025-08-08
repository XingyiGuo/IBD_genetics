#!/bin/bash

#SBATCH --job-name=ldsc.eur
#SBATCH --error=%x-%j.error
#SBATCH --out=%x-%j.out
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-0:0:0

## source activate ldsc # before running ldsc
## module load StdEnv/2020
## module load python/2.7.18
## module load bedtools
## source /home/lyul1/ldsc_env_py27/bin/activate


# For overall IBD:
/home/lyul1/ldsc/munge_sumstats.py \
--sumstats /nobackup/sbcs/keep/data/GuoLab/backup/lyul1/IBD/metal_update/withUKB/gc_0.01/GC_Final/processed/EUR_CD_META1.TBL.SNP2.txt \
--snp SNP \
--N 20093 \
--chunksize 500000 \
--a1 Allele2 \
--a2 Allele1 \
--out /nobackup/sbcs/keep/data/GuoLab/backup/lyul1/IBD/single_cell/sLDSCs/sumstats/eur_cd \
--merge-alleles /home/lyul1/LDSCref/w_hm3.snplist

/home/lyul1/ldsc/munge_sumstats.py \
--sumstats /nobackup/sbcs/keep/data/GuoLab/backup/lyul1/IBD/metal_update/withUKB/gc_0.01/GC_Final/processed/EUR_UC_META1.TBL.SNP2.txt \
--snp SNP \
--N 28485 \
--chunksize 500000 \
--a1 Allele2 \
--a2 Allele1 \
--out /nobackup/sbcs/keep/data/GuoLab/backup/lyul1/IBD/single_cell/sLDSCs/sumstats/eur_uc \
--merge-alleles /home/lyul1/LDSCref/w_hm3.snplist



/home/lyul1/ldsc/munge_sumstats.py \
--sumstats /nobackup/sbcs/keep/data/GuoLab/backup/lyul1/IBD/metal_update/withUKB/gc_0.01/GC_Final/processed/EUR_IBD_META1.TBL.SNP2.txt \
--snp SNP \
--N 49022 \
--chunksize 500000 \
--a1 Allele2 \
--a2 Allele1 \
--out /nobackup/sbcs/keep/data/GuoLab/backup/lyul1/IBD/single_cell/sLDSCs/sumstats/eur_ibd \
--merge-alleles /home/lyul1/LDSCref/w_hm3.snplist
