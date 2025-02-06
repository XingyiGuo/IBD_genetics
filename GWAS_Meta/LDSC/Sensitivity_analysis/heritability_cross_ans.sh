#!/bin/bash

#SBATCH --job-name=GC_corss
#SBATCH --error=%x-%j.error
#SBATCH --out=%x-%j.out
#SBATCH --mem=5G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-0:0:0


/home/lyul1/ldsc/ldsc.py \
--rg eur_cd.sumstats.gz,eas_cd.sumstats.gz \
--ref-ld-chr /home/lyul1/LDSCref/eur_w_ld_chr/ \
--w-ld-chr /home/lyul1/LDSCref/eur_w_ld_chr/ \
--out GC_eas_eur_cd


/home/lyul1/ldsc/ldsc.py \
--rg eur_ibd.sumstats.gz,eas_ibd.sumstats.gz \
--ref-ld-chr /home/lyul1/LDSCref/eur_w_ld_chr/ \
--w-ld-chr /home/lyul1/LDSCref/eur_w_ld_chr/ \
--out GC_eas_eur_ibd


# GC of cd uc
/home/lyul1/ldsc/ldsc.py \
--rg eur_uc.sumstats.gz,eas_uc.sumstats.gz \
--ref-ld-chr /home/lyul1/LDSCref/eur_w_ld_chr/ \
--w-ld-chr /home/lyul1/LDSCref/eur_w_ld_chr/ \
--out GC_eas_eur_uc