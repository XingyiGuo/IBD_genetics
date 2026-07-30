#!/bin/bash

#SBATCH --mem=60G 
#SBATCH --job-name=03_smr_besd_eur 
#SBATCH --error=%x-%j.error
#SBATCH --out=%x-%j.out
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-0:0:0



set -euo pipefail

/data/l2_bioinfo1/lyul1/tools/smr-1.4.1-linux-x86_64/smr \
  --eqtl-flist /data/l2_bioinfo1/lyul1/IBD/SMR/xQTLs_harmonized/EUR/ABS/UVA_ABS.cis_qtl_pairs.allchr.list \
  --make-besd \
  --out /data/l2_bioinfo1/lyul1/IBD/SMR/xQTLs_harmonized/EUR/ABS/UVA_ABS.cis_qtl_pairs.allchr
  
/data/l2_bioinfo1/lyul1/tools/smr-1.4.1-linux-x86_64/smr \
  --eqtl-flist /data/l2_bioinfo1/lyul1/IBD/SMR/xQTLs_harmonized/EUR/GOB/UVA_GOB.cis_qtl_pairs.allchr.list \
  --make-besd \
  --out /data/l2_bioinfo1/lyul1/IBD/SMR/xQTLs_harmonized/EUR/GOB/UVA_GOB.cis_qtl_pairs.allchr

/data/l2_bioinfo1/lyul1/tools/smr-1.4.1-linux-x86_64/smr \
  --eqtl-flist /data/l2_bioinfo1/lyul1/IBD/SMR/xQTLs_harmonized/EUR/APA/UVA_APA.cis_qtl_pairs.allchr.list \
  --make-besd \
  --out /data/l2_bioinfo1/lyul1/IBD/SMR/xQTLs_harmonized/EUR/APA/UVA_APA.cis_qtl_pairs.allchr

/data/l2_bioinfo1/lyul1/tools/smr-1.4.1-linux-x86_64/smr \
  --eqtl-flist /data/l2_bioinfo1/lyul1/IBD/SMR/xQTLs_harmonized/EUR/Expr/UVA_Expr.cis_qtl_pairs.allchr.list \
  --make-besd \
  --out /data/l2_bioinfo1/lyul1/IBD/SMR/xQTLs_harmonized/EUR/Expr/UVA_Expr.cis_qtl_pairs.allchr

   
/data/l2_bioinfo1/lyul1/tools/smr-1.4.1-linux-x86_64/smr \
  --eqtl-flist /data/l2_bioinfo1/lyul1/IBD/SMR/xQTLs_harmonized/EUR/STM/UVA_STM.cis_qtl_pairs.allchr.list \
  --make-besd \
  --out /data/l2_bioinfo1/lyul1/IBD/SMR/xQTLs_harmonized/EUR/STM/UVA_STM.cis_qtl_pairs.allchr
 
/data/l2_bioinfo1/lyul1/tools/smr-1.4.1-linux-x86_64/smr \
  --eqtl-flist /data/l2_bioinfo1/lyul1/IBD/SMR/xQTLs_harmonized/EUR/AlterSplice/UVA_AlterSplice.cis_qtl_pairs.allchr.list \
  --make-besd \
  --out /data/l2_bioinfo1/lyul1/IBD/SMR/xQTLs_harmonized/EUR/AlterSplice/UVA_AlterSplice.cis_qtl_pairs.allchr

echo "Finished EUR BESD conversion."