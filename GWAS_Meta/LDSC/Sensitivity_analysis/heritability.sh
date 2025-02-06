#!/bin/bash

#SBATCH --job-name=herti
#SBATCH --error=%x-%j.error
#SBATCH --out=%x-%j.out
#SBATCH --mem=5G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-0:0:0


# heritability of ibd EUR
/home/lyul1/ldsc/ldsc.py \
--h2 eur_ibd.sumstats.gz \
--ref-ld-chr /home/lyul1/LDSCref/eur_w_ld_chr/ \
--w-ld-chr /home/lyul1/LDSCref/eur_w_ld_chr/ \
--out h2_eur_ibd

/home/lyul1/ldsc/ldsc.py \
--h2 eur_cd.sumstats.gz \
--ref-ld-chr /home/lyul1/LDSCref/eur_w_ld_chr/ \
--w-ld-chr /home/lyul1/LDSCref/eur_w_ld_chr/ \
--out h2_eur_cd

/home/lyul1/ldsc/ldsc.py \
--h2 eur_uc.sumstats.gz \
--ref-ld-chr /home/lyul1/LDSCref/eur_w_ld_chr/ \
--w-ld-chr /home/lyul1/LDSCref/eur_w_ld_chr/ \
--out h2_eur_uc

# GC of cd uc
/home/lyul1/ldsc/ldsc.py \
--rg eur_cd.sumstats.gz,eur_uc.sumstats.gz \
--ref-ld-chr /home/lyul1/LDSCref/eur_w_ld_chr/ \
--w-ld-chr /home/lyul1/LDSCref/eur_w_ld_chr/ \
--out GC_eur_cd_uc


# heritability of ibd  EAS
/home/lyul1/ldsc/ldsc.py \
--h2 eas_ibd.sumstats.gz \
--ref-ld-chr /home/lyul1/LDSCref/eur_w_ld_chr/ \
--w-ld-chr /home/lyul1/LDSCref/eur_w_ld_chr/ \
--out h2_eas_ibd

/home/lyul1/ldsc/ldsc.py \
--h2 eas_cd.sumstats.gz \
--ref-ld-chr /home/lyul1/LDSCref/eur_w_ld_chr/ \
--w-ld-chr /home/lyul1/LDSCref/eur_w_ld_chr/ \
--out h2_eas_cd


/home/lyul1/ldsc/ldsc.py \
--h2 eas_uc.sumstats.gz \
--ref-ld-chr /home/lyul1/LDSCref/eur_w_ld_chr/ \
--w-ld-chr /home/lyul1/LDSCref/eur_w_ld_chr/ \
--out h2_eas_uc


# GC of cd uc
/home/lyul1/ldsc/ldsc.py \
--rg eas_cd.sumstats.gz,eas_uc.sumstats.gz \
--ref-ld-chr /home/lyul1/LDSCref/eur_w_ld_chr/ \
--w-ld-chr /home/lyul1/LDSCref/eur_w_ld_chr/ \
--out GC_eas_cd_uc
