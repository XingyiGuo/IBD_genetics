# IBD_genetics

## Overview

![My Image](./Figures/xx.png)

**Step1:** Identify novel GWAS loci for IBD/CD/UC

**Step2:** Identify novel genes associated with IBD/CD/UC through TWAS/AS-WAS/APA-WAS

## Methods
### 1. Identify novel GWAS loci


### 2. Identify novel genes associated with IBD/CD/UC through TWAS/AS-WAS/APA-WAS
We utilized the weight matrix and the summary statistics from IBD/CD/UC GWAS datasets from Europen ancestry and East Asian ancestry. We evaluated the association between gene expression (or AS, APA) and CRC risk under the sTF-TWAS framework.
For more information, please refer to our previous publication [PMID: 36402776](https://pubmed.ncbi.nlm.nih.gov/36402776/) and github repository [sTF-TWAS](https://github.com/XingyiGuo/TF-TWAS).

#### 2.1
Following [GTEx RNA-seq analyses pipelines](https://github.com/broadinstitute/gtex-pipeline), Events were quantile normalized and inverse normal transformed. The potential confounding factors (such as top five principal components (PCs), gender, potential batch effects, PEER factors) for events have been removed using [PEER](https://github.com/PMBio/peer).

Code: run_expression_peer.R


#### 2.2

Code: sTFTWAS_models.R; buildDB.R

#### 2.3 Identification of CRC risk events through associations
We utilized the weight matrix and the summary statistics from CRC GWAS datasets consisting of 186,072 individuals of European ancestry and 72,272 individuals of East Asian ancestry, we evaluated the association between gene expression (or AS, APA) and CRC risk under the sTF-TWAS framework.

- Executive code: \
conda activate [spredixcan](https://github.com/hakyimlab/MetaXcan/blob/master/software/SPrediXcan.py) \
./path/to/src/MetaXcan/software/SPrediXcan.py --model_db_path model.db --covariance model_cov.txt.gz --gwas_folder ./GWAS_SS/ --gwas_file_pattern ".*gz" --snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column BETA  --pvalue_column P --output_file  model.TWAS --verbosity 1


#### 2.4 Fine-mapping TWAS genes at genomic risk regions

- Executive code:


### Data Availability


## Contact
Lesley Lyu: linshuoshuo.lyu@vumc.org \
Qing Li: qing.li@vumc.org \
Xingyi Guo: xingyi.guo@vumc.org \
