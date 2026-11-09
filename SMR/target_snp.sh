#!/bin/bash
#SBATCH --account=l2_bioinfo1
#SBATCH --job-name=target_snp 
#SBATCH --error=%x-%j.error
#SBATCH --out=%x-%j.out
#SBATCH --mem=60G
##SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1 
#SBATCH --time=2-0:0:0 
#SBATCH --array=0-16%3 

set -euo pipefail


SCRIPT="/data/l2_bioinfo1/lyul1/IBD/SMR/scr_was_gene/01_target_snp.R"

GENE_LIST="/data/l2_bioinfo1/lyul1/IBD/SMR/ST22_genes_update.csv"
GTF="/data/l2_bioinfo1/lyul1/Reference/gencode.v26.annotation.gtf"
PANEL_TSV="/data/l2_bioinfo1/lyul1/IBD/SMR/xqtl_panels.txt"
GWAS_DIR="/data/l2_bioinfo1/lyul1/IBD/metal_update/withUKB/gc_0.01/GC_Final/processed/forsmr_afreq"

OUT_ROOT="/data/l2_bioinfo1/lyul1/IBD/SMR/target_snps"
WINDOW_BP="1000000"

SIF="/labs/accre_public/singularity/RStudio/4.5.0/R-RStudio-x86-64-v2.sif"

DISEASES=(
  "IBD"
  "CD"
  "UC"
)

if [[ ! -f "$SCRIPT" ]]; then
  echo "ERROR: missing R script: $SCRIPT" >&2
  exit 1
fi

if [[ ! -f "$PANEL_TSV" ]]; then
  echo "ERROR: missing panel TSV: $PANEL_TSV" >&2
  exit 1
fi

#panel_n=$(( $(wc -l < "$PANEL_TSV") - 1 ))
panel_n="$(awk 'NR > 1 && NF > 0 {n++} END {print n + 0}' "$PANEL_TSV")"
disease_n="${#DISEASES[@]}"
total=$(( panel_n * disease_n ))

task_id="${SLURM_ARRAY_TASK_ID:-0}"

if (( task_id >= total )); then
  echo "Task ${task_id} outside total jobs ${total}; exiting."
  exit 0
fi

disease_i=$(( task_id / panel_n ))
panel_i=$(( task_id % panel_n + 1 ))

disease="${DISEASES[$disease_i]}"

ancestry="$(
  awk -F'\t' -v row="$panel_i" '
    NR == row + 1 {print $1; exit}
  ' "$PANEL_TSV" | tr -d '\r\n'
)"

panel_label="$(
  awk -F'\t' -v row="$panel_i" '
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        if ($i == "panel_label") panel_col = i
      }
      next
    }
    NR == row + 1 {
      print $panel_col
      exit
    }
  ' "$PANEL_TSV" | tr -d '\r\n'
)"

out_dir="${OUT_ROOT}/${disease}/${ancestry}"

mkdir -p "$out_dir"

echo "Task ID: ${task_id}"
echo "Disease: ${disease}"
echo "Panel index: ${panel_i} / ${panel_n}"
echo "Panel label: ${panel_label}"
echo "Ancestry: ${ancestry}"
echo "Out dir: ${out_dir}"

singularity exec \
  --bind /data/l2_bioinfo1:/data/l2_bioinfo1 \
  --bind /home/lyul1:/mnt_home \
  --bind /cvmfs:/cvmfs \
  "$SIF" \
  Rscript "$SCRIPT" \
    --gene-list "$GENE_LIST" \
    --gtf "$GTF" \
    --panel-tsv "$PANEL_TSV" \
    --panel-index "$panel_i" \
    --disease "$disease" \
    --gwas-dir "$GWAS_DIR" \
    --out-dir "$out_dir" \
    --window-bp "$WINDOW_BP"

echo "Done."