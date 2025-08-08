#!/bin/bash

#SBATCH --job-name=ann
#SBATCH --error=%x-%j.error
#SBATCH --out=%x-%j.out
#SBATCH --mem=2G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-0:0:0


GENSETDIR=/nobackup/sbcs/lyul1/IBD/single_cell/topgenes_update_ENSG
ENSGDIR=/home/lyul1/LDSCref
PLINK_PREFIX=/home/lyul1/LDSCref/1000G_EUR_Phase3_plink
OUTPUT_DIR=/nobackup/sbcs/lyul1/IBD/single_cell/sLDSCs/annotations

mkdir -p "$OUTPUT_DIR"

for cell_base in ABS Stem GOB; do
    geneset_file="${GENSETDIR}/${cell_base}_top_genes_ENSG.txt"

    for chr in {1..22}; do
        output_file="${OUTPUT_DIR}/${cell_base}.${chr}.annot.gz"

        echo "Creating annotation for ${cell_base}, chr${chr}"

        /home/lyul1/ldsc/make_annot.py \
        	--gene-set-file "$geneset_file" \
		    --gene-coord-file "${ENSGDIR}/ENSG_coord.txt" \
		    --windowsize 100000 \
            --bimfile "${PLINK_PREFIX}/1000G.EUR.QC.${chr}.bim" \
            --annot-file "$output_file"

    done
done
