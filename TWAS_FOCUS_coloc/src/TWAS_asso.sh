#!/bin/bash

#download SPrediXcan.py from https://github.com/hakyimlab/MetaXcan/tree/master

/src/MetaXcan/software/SPrediXcan.py     --model_db_path ../Example_Data/CRC_TF_GTEx.db     --covariance ../Example_Data/CRC_TF_GTEx_cov.txt.gz     --gwas_folder ../Example_Data/EUR_IBD_perchr/     --gwas_file_pattern ".*".gz     --snp_column SNP     --effect_allele_column A1     --non_effect_allele_column A2     --beta_column BETA      --pvalue_column P     --output_file  CRC_TF_GTEx.EUR_IBD.TWAS     --verbosity 1