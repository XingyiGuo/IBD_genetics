library(data.table)
library(dplyr)
library(tidyr)

##### Selecting previously significant (p<5e-08) snps from 4 EUR GWASs, 1EAS GWAS, and 1EASEUR meta GWAS.

setwd("/nobackup/sbcs/lyul1/GWAS_SS/IBD")

############### EUR IBD ####################
# mvp
mvp<- "/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2023MVP/processed/"
mvp_cd<- fread(paste0(mvp,"Phe_555_1.EUR.CD_markername.txt"), fill = T) 
colnames(mvp_cd)[1]<- "SNP"
colnames(mvp_cd)[2]<- "CHR"
colnames(mvp_cd)[3]<- "POS"
mvp_uc<- fread(paste0(mvp,"Phe_555_2.EUR.UC_markername.txt"), fill = T) 
colnames(mvp_uc)[1]<- "SNP"
colnames(mvp_uc)[2]<- "CHR"
colnames(mvp_uc)[3]<- "POS"
mvp_ibd<- fread(paste0(mvp,"Phe_555.EUR.IBD_markername.txt"), fill = T) 
colnames(mvp_ibd)[1]<- "SNP"
colnames(mvp_ibd)[2]<- "CHR"
colnames(mvp_ibd)[3]<- "POS"

select_mvp<- function(gwas, trait) {
  gwas_filtered<- filter(gwas, gwas$pval < 5e-08)
  
  write.table(gwas_filtered, paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/mvp_", trait, "_5e-8.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  return(gwas_filtered)
}

mvp_cd<- select_mvp(mvp_cd , "Phe_555_1.EUR.CD") 
mvp_uc<- select_mvp(mvp_uc, "Phe_555_2.EUR.UC")
mvp_ibd<- select_mvp(mvp_ibd, "Phe_555.EUR.IBD")

# finn
finn<- "/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_FinnGen/R12/processed/"
finn_cd<- fread(paste0(finn, "finngen_R12_K11_CD_STRICT2_markername.txt"), fill = T)
colnames(finn_cd)[5]<- "SNP"
colnames(finn_cd)[1]<- "CHR"
colnames(finn_cd)[2]<- "POS"
finn_uc<- fread(paste0(finn, "finngen_R12_K11_UC_STRICT2_markername.txt"), fill = T) 
colnames(finn_uc)[5]<- "SNP"
colnames(finn_uc)[1]<- "CHR"
colnames(finn_uc)[2]<- "POS"
finn_ibd<- fread(paste0(finn, "finngen_R12_K11_IBD_STRICT_markername.txt"), fill = T)
colnames(finn_ibd)[5]<- "SNP"
colnames(finn_ibd)[1]<- "CHR"
colnames(finn_ibd)[2]<- "POS"

select_finn<- function(gwas, trait) {
  gwas_filtered<- filter(gwas, gwas$pval < 5e-08)
  
  write.table(gwas_filtered, paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/", trait, "_5e-8.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  return(gwas_filtered)
}

finn_cd<- select_finn(finn_cd , "finngen_R12_K11_CD_STRICT2")
finn_uc<- select_finn(finn_uc, "finngen_R12_K11_UC_STRICT2")
finn_ibd<- select_finn(finn_ibd, "finngen_R12_K11_IBD_STRICT")


# ng17 
ng17<- "/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2017NG/" 
ng17_cd<- fread(paste0(ng17, "cd_build37_40266_20161107.txt"), fill = T)
ng17_uc<- fread(paste0(ng17, "uc_build37_45975_20161107.txt"), fill = T)
ng17_ibd<- fread(paste0(ng17, "ibd_build37_59957_20161107.txt"), fill = T)

select_ng<- function(gwas, trait) {
  gwas_filtered<- filter(gwas, gwas$P.value < 5e-08)
  
  write.table(gwas_filtered, paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/ng_",  trait, "_5e-8.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  return(gwas_filtered)
}

ng17_cd<- select_ng(ng17_cd , "cd_build37")
ng17_uc<- select_ng(ng17_uc, "uc_build37")
ng17_ibd<- select_ng(ng17_ibd, "ibd_build37")


##### map to get significant snps with grch38 chr and pos
cd38<- fread("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2017NG/processed/cd_build38_markername.txt", fill = T)
uc38<- fread("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2017NG/processed/uc_build38_markername.txt", fill = T)
ibd38<- fread("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2017NG/processed/ibd_build38_markername.txt", fill = T)

merge_grch<- function(gwas, hg, trait) {
  gwas<- merge(gwas, hg[,c(1,18,19,25)], by="MarkerName")
  
  write.table(gwas, paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/ng_",  trait, "_5e-8.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  return(gwas)
}

ng17_cd<- merge_grch(ng17_cd, cd38, "cd_build38")
ng17_uc<- merge_grch(ng17_uc, uc38, "uc_build38")
ng17_ibd<- merge_grch(ng17_ibd, ibd38, "ibd_build38")

# By checking build37 and build38 significant snps, one snp was not successfully converted to build38 during liftover, so it was manually added to the file
ng17_ibd1<- fread("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/ng_ibd_build37_5e-8.txt", fill = T) 
diff_rows <- anti_join(ng17_ibd1, ng17_ibd[,1], by = "MarkerName") # manually add 9:139242225_T_C in ibd_build38_5e-8.txt
chrpos <- data.frame(CHR = 9, POS = 136347773, markernames = "9:136347773:C:T")
add<- cbind(diff_rows, chrpos) 
ng17_ibd<- rbind(ng17_ibd,add)
write.table(ng17_ibd, paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/ng_",  "ibd_build38", "_5e-8_2.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# ng17_cd<- fread("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/ng_cd_build38_5e-8.txt", fill = T) %>% select(, -16)
# ng17_uc<- fread("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/ng_uc_build38_5e-8.txt", fill = T) %>% select(, -16) 
# ng17_ibd<- fread("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/ng_ibd_build38_5e-8_2.txt", fill = T) %>% select(, -16) 

# map rsids
cd_rsid<- fread("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2017NG/processed/cd_build37_N_SNP.txt", fill = T)
uc_rsid<- fread("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2017NG/processed/uc_build37_N_SNP.txt", fill = T)
ibd_rsid<- fread("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2017NG/processed/ibd_build37_N_SNP.txt", fill = T)

map_rsid_ng17<- function(gwas, rsid, trait) {
  
  gwas<- merge(gwas, rsid[,c(5,24)], by = "MarkerName", all.x = T)
  colnames(gwas)[19]<- "SNP"
  
  write.table(gwas, paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/ng_",  trait, "_5e-8.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  return(gwas)
}

ng17_cd<- map_rsid_ng17(ng17_cd, cd_rsid, "cd_build38")
ng17_uc<- map_rsid_ng17(ng17_uc, uc_rsid, "uc_build38")
ng17_ibd<- map_rsid_ng17(ng17_ibd, ibd_rsid, "ibd_build38")

select_ukb<- function(trait) {
  
  gwas<- fread(paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_Pan_UKB/processed/", trait, "_markername.txt"), fill = T, header = T) 
  
  gwas$P.value<- as.numeric(gwas$P.value)
  gwas_filtered<- filter(gwas, gwas$P.value < 5e-08)
  gwas_filtered$SNP<- NA
  
  write.table(gwas_filtered, paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/ukb_",  trait, "_5e-8.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  return(gwas_filtered)
}

ukb_cd<- select_ukb("555.1_CD_build38")
ukb_uc<- select_ukb("555.2_UC_build38")
ukb_ibd<- select_ukb("555_IBD_build38")

####### combine significant snps from mvp, finn, and ng17 together ########

header<- c("SNP","CHR", "POS", "REF", "EA", "beta", "se", "pval", "markernames")

ukb_cd<- ukb_cd %>% 
  select(, c("SNP","CHR", "BP", "ref", "alt", "beta_EUR", "se_EUR", "P.value", "markernames")) %>%  setNames(header) %>% 
  mutate(dataset = "UKB")

ukb_uc<- ukb_uc %>% 
  select(, c("SNP","CHR", "BP", "ref", "alt", "beta_EUR", "se_EUR", "P.value", "markernames")) %>%  setNames(header) %>% 
  mutate(dataset = "UKB")

ukb_ibd<- ukb_ibd %>% 
  select(, c("SNP","CHR", "BP", "ref", "alt", "beta_EUR", "se_EUR", "P.value", "markernames")) %>%  setNames(header) %>% 
  mutate(dataset = "UKB")


mvp_cd<-  mvp_cd %>% 
  select(, c("SNP","CHR", "POS", "ref", "ea", "beta", "se", "pval", "markernames")) %>%  setNames(header) %>% 
  mutate(dataset = "mvp")
mvp_uc<-  mvp_uc %>% 
  select(, c("SNP","CHR", "POS", "ref", "ea", "beta", "se", "pval", "markernames")) %>%  setNames(header) %>% 
  mutate(dataset = "mvp")
mvp_ibd<- mvp_ibd %>% 
  select(, c("SNP","CHR", "POS", "ref", "ea", "beta", "se", "pval", "markernames")) %>%  setNames(header) %>% 
  mutate(dataset = "mvp")

finn_cd<- finn_cd %>% 
  select(, c("SNP","CHR", "POS", "ref",  "alt","beta", "sebeta", "pval", "markernames")) %>%  setNames(header) %>% 
  mutate(dataset = "finn")
finn_uc<-  finn_uc %>% 
  select(, c("SNP","CHR", "POS", "ref",  "alt","beta", "sebeta", "pval", "markernames")) %>%  setNames(header) %>% 
  mutate(dataset = "finn")
finn_ibd<-  finn_ibd %>% 
  select(, c("SNP","CHR", "POS", "ref",  "alt","beta", "sebeta", "pval", "markernames")) %>%  setNames(header) %>% 
  mutate(dataset = "finn")

ng17_cd<-  ng17_cd %>% 
  select(, c("SNP","CHR", "POS", "Allele1", "Allele2", "Effect", "StdErr", "P.value", "markernames"))  %>%  setNames(header) %>% 
  mutate(dataset = "10.1038/ng.3359.")
ng17_uc<-  ng17_uc %>% 
  select(, c("SNP","CHR", "POS", "Allele1", "Allele2", "Effect", "StdErr", "P.value", "markernames")) %>%  setNames(header) %>% 
  mutate(dataset = "10.1038/ng.3359.")
ng17_ibd<-  ng17_ibd %>% 
  select(, c("SNP","CHR", "POS", "Allele1", "Allele2", "Effect", "StdErr", "P.value", "markernames")) %>%  setNames(header) %>% 
  mutate(dataset = "10.1038/ng.3359.")

cd<- rbind(mvp_cd, finn_cd, ng17_cd, ukb_cd)
write.table(cd, "/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/cd_5e-8.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

uc<- rbind(mvp_uc, finn_uc, ng17_uc, ukb_uc)
write.table(uc, "/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/uc_5e-8.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

ibd<- rbind(mvp_ibd, finn_ibd, ng17_ibd, ukb_ibd)
write.table(ibd, "/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/ibd_5e-8.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


############### EAS EUR Meta IBD ####################
# eas
eas<- "/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2023NG/processed/"
eas_cd<- fread(paste0(eas,"ibd_EAS_SiKJ_meta_CD_markername.txt"), fill = T) 
eas_uc<- fread(paste0(eas,"ibd_EAS_SiKJ_meta_UC_markername.txt"), fill = T) 
eas_ibd<- fread(paste0(eas,"ibd_EAS_SiKJ_meta_IBD_markername.txt"), fill = T) 
# eas eur meta
eas_eur_cd<- fread(paste0(eas,"ibd_EAS_EUR_SiKJEF_meta_CD_markername.txt"), fill = T) 
eas_eur_uc<- fread(paste0(eas,"ibd_EAS_EUR_SiKJEF_meta_UC_markername.txt"), fill = T) 
eas_eur_ibd<- fread(paste0(eas,"ibd_EAS_EUR_SiKJEF_meta_IBD_markername.txt"), fill = T) 


select_eas<- function(gwas, trait, rsid) {
  
  colnames(gwas)[10]<- "P.value"
  
  gwas_filtered<- filter(gwas, gwas$P.value < 5e-08)
  
  gwas_filtered$SNP<- NA
  
  write.table(gwas_filtered, paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/EAS_", trait, "_5e-8.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  return(gwas_filtered)
}

eas_cd<- select_eas(eas_cd , "CD")
eas_uc<- select_eas(eas_uc, "UC")
eas_ibd<- select_eas(eas_ibd, "IBD")

eas_eur_cd<- select_eas(eas_eur_cd, "EUR.meta_CD")
eas_eur_uc<- select_eas(eas_eur_uc, "EUR.meta_UC")
eas_eur_ibd<- select_eas(eas_eur_ibd, "EUR.meta_IBD")


# read in significant snps from eas/easeur gwas
# eas_cd<- read.table(paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/","EAS_CD_5e-8.txt"), header = T, sep = "\t", stringsAsFactors = F)
# eas_uc<- read.table(paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/", "EAS_UC_5e-8.txt"), header = T, sep = "\t", stringsAsFactors = F)
# eas_ibd<- read.table(paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/", "EAS_IBD_5e-8.txt"), header = T, sep = "\t", stringsAsFactors = F)

# eas_eur_cd<- read.table(paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/","EAS_EUR.meta_CD_5e-8.txt"), header = T, sep = "\t", stringsAsFactors = F)
# eas_eur_uc<- read.table(paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/", "EAS_EUR.meta_UC_5e-8.txt"), header = T, sep = "\t", stringsAsFactors = F)
# eas_eur_ibd<- read.table(paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/", "EAS_EUR.meta_IBD_5e-8.txt"), header = T, sep = "\t", stringsAsFactors = F)

# cd<- read.table(paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/","cd_5e-8.txt"), header = T, sep = "\t", stringsAsFactors = F)
# uc<- read.table(paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/","uc_5e-8.txt"), header = T, sep = "\t", stringsAsFactors = F)
# ibd<- read.table(paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/","ibd_5e-8.txt"), header = T, sep = "\t", stringsAsFactors = F)

# cd
eas_cd<- eas_cd %>% 
  select(, c("SNP","CHR", "BP", "Allele1", "Allele2", "Effect", "StdErr", "P.value", "markernames")) %>%  setNames(header) %>% 
  mutate(dataset = "10.1038/s41588-023-01384-0.eas")
eas_eur_cd<- eas_eur_cd %>% 
  select(, c("SNP","CHR", "BP", "Allele1", "Allele2", "Effect", "StdErr", "P.value", "markernames")) %>%  setNames(header) %>% 
  mutate(dataset = "10.1038/s41588-023-01384-0.eas.eur")

cd1<- rbind(cd, eas_cd)
write.table(cd1, "/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/EAS+EUR_cd_5e-8.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

cd<- rbind(cd, eas_cd, eas_eur_cd)
write.table(cd, "/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/EASEUR_cd_5e-8.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# uc
eas_uc<- eas_uc %>% 
  select(, c("SNP","CHR", "BP", "Allele1", "Allele2", "Effect", "StdErr", "P.value", "markernames")) %>%  setNames(header) %>% 
  mutate(dataset = "10.1038/s41588-023-01384-0.eas")
eas_eur_uc<- eas_eur_uc %>% 
  select(, c("SNP","CHR", "BP", "Allele1", "Allele2", "Effect", "StdErr", "P.value", "markernames")) %>%  setNames(header) %>% 
  mutate(dataset = "10.1038/s41588-023-01384-0.eas.eur")

uc1<- rbind(uc, eas_cd)
write.table(uc1, "/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/EAS+EUR_uc_5e-8.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

uc<- rbind(uc, eas_uc, eas_eur_uc)
write.table(uc, "/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/EASEUR_uc_5e-8.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# ibd
eas_ibd<- eas_ibd %>% 
  select(, c("SNP","CHR", "BP", "Allele1", "Allele2", "Effect", "StdErr", "P.value", "markernames")) %>%  setNames(header) %>% 
  mutate(dataset = "10.1038/s41588-023-01384-0.eas")
eas_eur_ibd<- eas_eur_ibd %>% 
  select(, c("SNP","CHR", "BP", "Allele1", "Allele2", "Effect", "StdErr", "P.value", "markernames")) %>%  setNames(header) %>% 
  mutate(dataset = "10.1038/s41588-023-01384-0.eas.eur")

ibd1<- rbind(ibd, eas_ibd)
write.table(ibd1, "/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/EAS+EUR_ibd_5e-8.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

ibd<- rbind(ibd, eas_ibd, eas_eur_ibd)
write.table(ibd, "/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/EASEUR_ibd_5e-8.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# remove duplicate SNPs
read_sig<- function(trait) {
  
  data<- fread(paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/", trait, "_5e-8.txt"), fill = T)
  data<- as.data.frame(data)
  data <- data %>%
    mutate(
      markernames = as.character(markernames), # Convert MarkerName to character
      pval = as.numeric(pval)               # Convert pval to numeric
    ) %>%
    arrange(markernames, pval) %>%           # Sort by MarkerName and pval
    group_by(markernames) %>%
    slice_min(order_by = pval, n = 1) %>%   # Keep the row with the smallest pval
    ungroup()
  
  write.table(data, paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/", trait, "_5e-8_nodup.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  return(data)
}

# EUR
kn_cd_sig <- read_sig("cd")
kn_uc_sig <- read_sig("uc")
kn_ibd_sig <- read_sig("ibd")
# CD and UC aassociated SNPs are also considered as IBD-associated
kn_ibd_all_sig<- rbind(kn_cd_sig, kn_uc_sig, kn_ibd_sig) %>% # Combine data frames
  mutate(
    markernames = as.character(markernames),  # Ensure markernames is character
    pval = as.numeric(pval)                  # Ensure pval is numeric
  ) %>%
  arrange(markernames, pval) %>%             # Sort by markernames and pval
  group_by(markernames) %>%
  slice_min(order_by = pval, n = 1) %>%      # Select row with the smallest pval
  ungroup()
write.table(kn_ibd_all_sig, paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/", "kn_ibd_all_sig", "_5e-8_nodup.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# EAS
kn_eascd_sig <- read_sig("EASEUR_cd")
kn_easuc_sig <- read_sig("EASEUR_uc")
kn_easibd_sig <- read_sig("EASEUR_ibd")

kn_easibd_all_sig<- rbind(kn_eascd_sig, kn_easuc_sig, kn_easibd_sig) %>%
  mutate(
    markernames = as.character(markernames),  # Ensure markernames is character
    pval = as.numeric(pval)                  # Ensure pval is numeric
  ) %>%
  arrange(markernames, pval) %>%             # Sort by markernames and pval
  group_by(markernames) %>%
  slice_min(order_by = pval, n = 1) %>%      # Select row with the smallest pval
  ungroup()
write.table(kn_easibd_all_sig, paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/", "kn_easibd_all_sig", "_5e-8_nodup.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
