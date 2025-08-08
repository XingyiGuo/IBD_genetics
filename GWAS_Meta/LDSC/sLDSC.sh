#!/bin/bash

#SBATCH --job-name=sldsc
#SBATCH --error=%x-%j.error
#SBATCH --out=%x-%j.out
#SBATCH --mem=5G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-0:0:0


ANODIR=/nobackup/sbcs/keep/data/GuoLab/backup/lyul1/IBD/single_cell/sLDSCs/annotations
W_LD_CHR=/nobackup/sbcs/keep/data/GuoLab/backup/lyul1/LDSCref/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
FRQ_DIR=/nobackup/sbcs/keep/data/GuoLab/backup/lyul1/LDSCref/1000G_Phase3_frq/1000G.EUR.QC.
SUMSTS_DIR=/nobackup/sbcs/keep/data/GuoLab/backup/lyul1/IBD/single_cell/sLDSCs/sumstats
OUTPUT_DIR=/nobackup/sbcs/keep/data/GuoLab/backup/lyul1/IBD/single_cell/sLDSCs/results

mkdir -p "$OUTPUT_DIR"

# Loop over cell types and traits ABS Stem GOB 
for celltype in ABS Stem GOB  ; do
  for trait in ibd cd uc ; do

    /home/lyul1/ldsc/ldsc.py \
      --h2 ${SUMSTS_DIR}/eur_${trait}.sumstats.gz \
      --ref-ld-chr ${ANODIR}/${celltype}. \
      --w-ld-chr ${W_LD_CHR} \
      --overlap-annot \
      --print-coefficients \
      --frqfile-chr ${FRQ_DIR} \
      --out ${OUTPUT_DIR}/eur_${trait}_${celltype}

  done
done
