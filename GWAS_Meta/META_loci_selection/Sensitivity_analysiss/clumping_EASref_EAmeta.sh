#!/bin/bash

#SBATCH --job-name=clpEAmeta
#SBATCH --error=%x-%j.error
#SBATCH --out=%x-%j.out
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-2:0:0



# set dir
REF_DIR="/nobackup/sbcs/chenz27/For_others/GuoLab/Reference/1KG_hg38/EAS/processed"
GWAS_DIR="/nobackup/sbcs/lyul1/IBD/metal_update/withoutUKB/gc_0.01/GC_Final/clumping/input"
OUT_DIR="/nobackup/sbcs/lyul1/IBD/metal_update/withoutUKB/gc_0.01/GC_Final/clumping/output_EAS.EURmeta/EASref"
EX_DIR="/nobackup/sbcs/lyul1/Reference/1KG_hg38/allsamples/EAS/processed/duplicateid"

mkdir -p "$EX_DIR"
mkdir -p "$OUT_DIR"

# GWAS result files
GWAS_FILES=( "EASEUR_CD_500k.txt" "EASEUR_UC_500k.txt" "EASEUR_IBD_500k.txt")

# loop through each GWAS file
for gwas_file in "${GWAS_FILES[@]}"; do
    # Extract the base name of the GWAS file (without extension)
    base_name=$(basename "$gwas_file" .txt)

    # loop through each chr
    for chr in {1..22}; do
            # generate a list of variants with duplicate/missing IDs (".")  missing_id_exclude_chr1.txt
       #  awk '$2 == "."' "$REF_DIR/chr${chr}.EUR.bim" > "$EX_DIR/missing_id_exclude_chr${chr}.txt"   
        awk '{if(a[$2]++) print $2}' "$REF_DIR/chr${chr}.EAS.bim" > "$EX_DIR/duplicate_exclude_chr${chr}.txt"



        # LD clumping for the current GWAS file and chr   chr1.EUR.bim (excluded "." id) 
                          
        /home/lyul1/plink --bfile "$REF_DIR/chr${chr}.EAS" \
                  --clump "$GWAS_DIR/$gwas_file" \
                  --clump-p1 5e-8 \
                  --clump-r2 0.6 \
                  --clump-kb 500 \
                  --keep-allele-order \
                  --exclude "$EX_DIR/duplicate_exclude_chr${chr}.txt" \
                  --out "$OUT_DIR/${base_name}_chr${chr}_clumped"

    done
done

# Combine results for each GWAS file across all chromosomes
for gwas_file in "${GWAS_FILES[@]}"; do
    base_name=$(basename "$gwas_file" .txt)
    cat "$OUT_DIR/${base_name}_chr"*"_clumped.clumped" > "$OUT_DIR/${base_name}_all_clumped.clumped"
done