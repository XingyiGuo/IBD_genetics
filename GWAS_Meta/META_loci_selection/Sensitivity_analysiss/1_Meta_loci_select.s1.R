library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)
library(readxl)

####### selecting significant snps (p<5e-08) from GWAS meta #########

setwd("/nobackup/sbcs/lyul1/IBD/metal_update/withoutUKB/gc_0.01")

####### read in previous gwas significant snps #######
kn_eascd_sig<- fread("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/EASEUR_cd_5e-8_nodup.txt")  %>% subset(CHR %in% 1:22)
colnames(kn_eascd_sig)[9]<- "MarkerName"

kn_easuc_sig<- fread("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/EASEUR_uc_5e-8_nodup.txt")  %>% subset(CHR %in% 1:22)
colnames(kn_easuc_sig)[9]<- "MarkerName"

kn_easibd_sig<- fread("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/EASEUR_ibd_5e-8_nodup.txt")  %>% subset(CHR %in% 1:22)
colnames(kn_easibd_sig)[9]<- "MarkerName"
kn_easibd_all_sig<- rbind(kn_eascd_sig, kn_easuc_sig, kn_easibd_sig)  %>% subset(CHR %in% 1:22)

####### Meta Results #######
## EUR ## Allele 1 = EA; Allele 2 = NEA
# processed 
filt<- "/nobackup/sbcs/lyul1/IBD/metal_update/withoutUKB/gc_0.01/GC_Final/processed/"
# cd
cd<-  read.table(paste0(filt, "EUR_CD_META1.TBL.txt"), header=T, sep="\t", stringsAsFactors = F) 
# uc
uc<-  read.table(paste0(filt, "EUR_UC_META1.TBL.txt"), header=T, sep="\t", stringsAsFactors = F)
# ibd
ibd<-  read.table(paste0(filt, "EUR_IBD_META1.TBL.txt"), header=T, sep="\t", stringsAsFactors = F)

## EASEUR ##
# cd
eascd<- read.table(paste0(filt, "EASEUR_CD_META1.TBL.txt"), header = T, sep = "\t", stringsAsFactors = F)
# uc
easuc<- read.table(paste0(filt, "EASEUR_UC_META1.TBL.txt"), header = T, sep = "\t", stringsAsFactors = F)
# ibd
easibd<- read.table(paste0(filt, "EASEUR_IBD_META1.TBL.txt"), header = T, sep = "\t", stringsAsFactors = F)

select_sig<- function(meta, known_snp, trait) {

  # select 1: significant snps (5e-08) in gwas meta
  
  meta_sig1<- filter(meta, meta$P.value<5e-08) 
  
  meta_sig1<- meta_sig1 %>% mutate(
    MarkerName_dup=MarkerName) %>% 
    separate(MarkerName_dup, into = c("CHR", "BP"), sep = ":") %>% 
    arrange(CHR, BP)
  
  # select 2: significant snps (5e-08) in gwas meta but not in previous gawss
  
  meta_sig2 <- anti_join(meta_sig1, known_snp, by = "MarkerName")
  
  return(list(sig1 = meta_sig1, sig2 = meta_sig2))
}

sig_cd<- select_sig(cd, kn_eascd_sig, "EUR_CD_META")
sig_uc<- select_sig(uc, kn_easuc_sig, "EUR_UC_META")
sig_ibd<- select_sig(ibd, kn_easibd_all_sig, "EUR_IBD_META")

sig_eascd<- select_sig(eascd, kn_eascd_sig, "EASEUR_CD_META")
sig_easuc<- select_sig(easuc, kn_easuc_sig, "EASEUR_UC_META")
sig_easibd<- select_sig(easibd, kn_easibd_all_sig, "EASEUR_IBD_META")

# write results to an excel
wr_tables <- function(lists, traits) {
  wb <- createWorkbook()
  
  for (name in names(lists)) {
    addWorksheet(wb, sheetName = name)
    writeData(wb, sheet = name, x = lists[[name]])
  }
  
  saveWorkbook(wb, paste0("/nobackup/sbcs/lyul1/IBD/metal_update/withoutUKB/gc_0.01/select/final_gc_0.01_", traits, "_5e-08.xlsx"), overwrite = TRUE)
}
wr_tables(sig_cd, "EUR_CD")
wr_tables(sig_uc, "EUR_UC")
wr_tables(sig_ibd, "EUR_IBD")

wr_tables(sig_eascd, "EASEUR_CD")
wr_tables(sig_easuc, "EASEUR_UC")
wr_tables(sig_easibd, "EASEUR_IBD")

