library(data.table)
library(dplyr)
library(tidyr)


###### calculate lambda for gwas meta ss ########

filt<- "/nobackup/sbcs/lyul1/IBD/metal_update/withoutUKB/gc_0.01/GC_Final/final_gc_"

cal_lambda<- function(trait) {
  
  # read in gwas ss
  data<- read.table(trait, header=T, sep="\t", stringsAsFactors = F)
  data$Effect<- as.numeric(data$Effect)
  data$StdErr<- as.numeric(data$StdErr)
  
  # calculate z score
  data$z<- data$Effect/data$StdErr
  
  # (1) Convert your output to chi-squared values
  # For z-scores, just square them
  chisq <- data$z^2
  
  # (2) Calculate lambda gc (λgc)
  median(chisq)/qchisq(0.5,1)
  
}



cd<- cal_lambda(paste0(filt, "EUR_CD_META1.TBL"))
uc<- cal_lambda(paste0(filt, "EUR_UC_META1.TBL"))
ibd<- cal_lambda(paste0(filt, "EUR_IBD_META1.TBL"))
eascd<- cal_lambda(paste0(filt, "EASEUR_CD_META1.TBL"))
easuc<- cal_lambda(paste0(filt, "EASEUR_UC_META1.TBL"))
easibd<- cal_lambda(paste0(filt, "EASEUR_IBD_META1.TBL"))

lambda<- c(cd, uc, ibd, eascd, easuc, easibd) # cd, uc, ibd,  "cd", "uc", "ibd", 
colname<- c("cd", "uc", "ibd", "eascd", "easuc", "easibd")
rb<- rbind(colname,  lambda) %>% as.data.frame()
rb

write.csv(rb, paste0("/nobackup/sbcs/lyul1/IBD/metal_update/withoutUKB/gc_0.01/select/lambda_meta_", "final_gc_0.01", ".csv"), col.names = F, row.names = F)
