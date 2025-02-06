library(data.table)
library(dplyr)
library(tidyr)


###### calculate lambda for each individual gwas ss ########

cal_lambda<- function(trait) {
  
  # read in gwas ss and keep only chr:22
  data<- read.table(trait, header=T, sep="\t", stringsAsFactors = F)
  data$P.value<- as.numeric(data$P.value)

  # calculate z score
  #data$z<- data$Effect/data$StdErr
  
  #data$z<- data$beta/data$sebeta
  
  # (1) Convert your output to chi-squared values
  # For z-scores, just square them
  #chisq <- data$z^2
  
  # For chi-squared values, keep as is
  #chisq <- data$chisq
  
  # For p-values, calculate chi-squared statistic
  
  chisq <- qchisq(1-data$P.value,1)

  
  # (2) Calculate lambda gc (λgc)
  median(chisq)/qchisq(0.5,1)
  
}

cal_lambda2<- function(trait) {
  
  # read in gwas ss and keep only chr:22
  data<- read.table(trait, header=T, sep="\t", stringsAsFactors = F)
  
  # calculate z score
  #data$z<- data$Effect/data$StdErr
  
  #data$z<- data$beta/data$sebeta
  
  # (1) Convert your output to chi-squared values
  # For z-scores, just square them
  #chisq <- data$z^2
  
  # For chi-squared values, keep as is
  #chisq <- data$chisq
  
  # For p-values, calculate chi-squared statistic
  # chisq <- qchisq(1-data$P.value,1)
  chisq <- qchisq(1-data$pval,1)
  
  # (2) Calculate lambda gc (λgc)
  median(chisq)/qchisq(0.5,1)
  
}

eascd<- cal_lambda("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2023NG/processed/ibd_EAS_SiKJ_meta_CD_markername_rmdup.txt")
easuc<- cal_lambda("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2023NG/processed/ibd_EAS_SiKJ_meta_UC_markername_rmdup.txt")
easibd<- cal_lambda("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2023NG/processed/ibd_EAS_SiKJ_meta_IBD_markername_rmdup.txt")

mvp_cd<- cal_lambda2("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2023MVP/processed/Phe_555_1.EUR.CD_markername_rmdup.txt")
mvp_uc<- cal_lambda2("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2023MVP/processed/Phe_555_2.EUR.UC_markername_rmdup.txt")
mvp_ibd<- cal_lambda2("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2023MVP/processed/Phe_555.EUR.IBD_markername_rmdup.txt")

finn_cd<- cal_lambda2("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_FinnGen/R12/processed/finngen_R12_K11_CD_STRICT2_markername_rmdup_chr1.22.txt")
finn_uc<- cal_lambda2("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_FinnGen/R12/processed/finngen_R12_K11_UC_STRICT2_markername_rmdup_chr1.22.txt")
finn_ibd<- cal_lambda2("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_FinnGen/R12/processed/finngen_R12_K11_IBD_STRICT_markername_rmdup_chr1.22.txt")

ukb_cd<- cal_lambda("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_Pan_UKB/processed/555.1_CD_build38_markername_rmdup_chr1.22.txt")
ukb_uc<- cal_lambda("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_Pan_UKB/processed/555.2_UC_build38_markername_rmdup_chr1.22.txt")
ukb_ibd<- cal_lambda("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_Pan_UKB/processed/555_IBD_build38_markername_rmdup_chr1.22.txt")

ng17_cd<- cal_lambda("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2017NG/processed/cd_build38_markername_rmdup_maf_chr1.22.txt")
ng17_uc<- cal_lambda("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2017NG/processed/uc_build38_markername_rmdup_maf_chr1.22.txt")
ng17_ibd<- cal_lambda("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2017NG/processed/ibd_build38_markername_rmdup_maf_chr1.22.txt")


lambda<- c(mvp_cd, mvp_uc, mvp_ibd, eascd, easuc, easibd, finn_cd, finn_uc, finn_ibd, ukb_cd, ukb_uc, ukb_ibd, ng17_cd, ng17_uc, ng17_ibd) # cd, uc, ibd,  "cd", "uc", "ibd", 
colname<- c("mvp_cd", "mvp_uc", "mvp_ibd", "eascd", "easuc", "easibd", "finn_cd", "finn_uc", "finn_ibd", "ukb_cd", "ukb_uc", "ukb_ibd", "ng17_cd", "ng17_uc", "ng17_ibd")
rb<- rbind(colname,  lambda) %>% as.data.frame()
rb

write.csv(rb, paste0("/nobackup/sbcs/lyul1/IBD/metal_update/withUKB/gc_0.01/select/lambda_individual.csv"), col.names = F, row.names = F)

