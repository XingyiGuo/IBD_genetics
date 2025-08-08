#!/bin/bash

#download SPrediXcan.py from https://github.com/hakyimlab/MetaXcan/tree/master

/src/MetaXcan/software/SPrediXcan.py     --model_db_path ../data/CRC_TF_GTEx.db     --covariance ../data/CRC_TF_GTEx_cov.txt.gz     --gwas_folder ../data/EUR_IBD_perchr/     --gwas_file_pattern ".*".gz     --snp_column SNP     --effect_allele_column A1     --non_effect_allele_column A2     --beta_column BETA      --pvalue_column P     --output_file  CRC_TF_GTEx.EUR_IBD.TWAS     --verbosity 1