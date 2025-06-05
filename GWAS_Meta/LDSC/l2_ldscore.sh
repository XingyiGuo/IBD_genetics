#!/bin/bash

#SBATCH --job-name=preldsc
#SBATCH --error=%x-%j.error
#SBATCH --out=%x-%j.out
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-0:0:0




# Set paths
ANODIR=/nobackup/sbcs/keep/data/GuoLab/backup/lyul1/IBD/single_cell/sLDSCs/annotations
PLINK_PREFIX=/nobackup/sbcs/keep/data/GuoLab/backup/lyul1/LDSCref/1000G_EUR_Phase3_plink/1000G.EUR.QC
#OUTPUT_DIR=/nobackup/sbcs/keep/data/GuoLab/backup/lyul1/IBD/single_cell/sLDSCs/pre_l2_ldscore
HM3_DIR=/home/lyul1/LDSCref/hapmap3_snps

mkdir -p "$OUTPUT_DIR"

for celltype in ABS Stem GOB ; do
    for chr in {1..22}; do
        anofile="${ANODIR}/${celltype}.${chr}.annot.gz"

        echo "Processing $celltype on chr$chr"

        /home/lyul1/ldsc/ldsc.py \
            --l2 \
            --bfile "${PLINK_PREFIX}.${chr}" \
            --ld-wind-cm 1 \
            --annot "$anofile" \
            --thin-annot \
            --out "${ANODIR}/${celltype}.${chr}" \
            --print-snps "${HM3_DIR}/hm.${chr}.snp"
    done
done
