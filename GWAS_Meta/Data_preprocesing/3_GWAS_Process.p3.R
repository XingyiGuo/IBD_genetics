library(data.table)
library(dplyr)
library(tidyr)


###### map rsid, allele freq in EAS ########

### EAS meta
setwd<- "/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2023NG/processed"

# loop through chr 1 to 22
ref_process_eas <- function(chr_id) {
  
  # read bim files
  path_bim <- sprintf("/nobackup/sbcs/chenz27/For_others/GuoLab/Reference/1KG_hg38/EAS/processed/chr%d.EAS.bim", chr_id)
  
  ref_bim <- fread(path_bim)
  
  header1<- c("CHR",	"SNP",	"GENO",	"BP",	"ALT",	"REF")
  colnames(ref_bim)<- header1
  ref_bim<-  ref_bim %>% select(, -c("GENO"))
  ref_bim$CHR<- as.integer(ref_bim$CHR)
  ref_bim$BP<- as.integer(ref_bim$BP)  
  
  ref_bim <- ref_bim %>%
    rowwise() %>%
    mutate(Allele1_sorted = pmin(ALT, REF),
           Allele2_sorted = pmax(ALT, REF)) %>%
    ungroup()

    
  return(ref_bim)
  
}


map_ref_to_gwas_eas <- function(gwas, chr_id) {
  
  # read processed ref data for the given chromosome
  ref <- ref_process_eas(chr_id)
  
  # filter GWAS data for the current chromosome
  gwas$CHR<- as.integer(gwas$CHR)
  gwas$BP<- as.integer(gwas$BP)  
  gwas_chr <- gwas %>%
    filter(CHR == chr_id)
  
  # merge bim.frq with gwas meta results
  merged_data <- merge(gwas_chr, ref, by = c("CHR", "BP", "Allele1_sorted", "Allele2_sorted"), all.x = T)
  
  return(merged_data)
  
}

res_process_eas <- function(merged, trait) {

  gwas <- fread(merged, header = TRUE, sep = "\t", fill = T)
  
  # map by chr
  mapped <- lapply(1:22, function(chr_id) map_ref_to_gwas_eas(gwas, chr_id)) %>%
    bind_rows() %>% 
    arrange(CHR, BP)
  
  write.table(mapped, trait, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

  return(mapped)
}
  

res_process_eas("ibd_EAS_SiKJ_meta_CD_markername.txt", "ibd_EAS_SiKJ_meta_CD_markername_SNP.txt")
res_process_eas("ibd_EAS_SiKJ_meta_UC_markername.txt", "ibd_EAS_SiKJ_meta_UC_markername_SNP.txt")
res_process_eas("ibd_EAS_SiKJ_meta_IBD_markername.txt", "ibd_EAS_SiKJ_meta_IBD_markername_SNP.txt")
res_process_eas("ibd_EAS_EUR_SiKJEF_meta_CD_markername.txt", "ibd_EAS_EUR_SiKJEF_meta_CD_markername_SNP.txt")
res_process_eas("ibd_EAS_EUR_SiKJEF_meta_UC_markername.txt", "ibd_EAS_EUR_SiKJEF_meta_UC_markername_SNP.txt")
res_process_eas("ibd_EAS_EUR_SiKJEF_meta_IBD_markername.txt", "ibd_EAS_EUR_SiKJEF_meta_IBD_markername_SNP.txt")

# remove the X chrs
# EAS
                   
# east asian IBD | Allele1: reference allele Allele2: effect allele
# eas<- "/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2023NG/processed/"

eas_ibd <- fread(paste0(eas,"ibd_EAS_SiKJ_meta_IBD_markername_rmdup.txt"), fill = T)   %>% subset(CHR %in% 1:22)
write.table(eas_ibd, paste0(eas, "ibd_EAS_SiKJ_meta_IBD_markername_rmdup_chr1.22.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

eas_cd<- fread(paste0(eas, "ibd_EAS_SiKJ_meta_CD_markername_rmdup.txt"), fill = T)   %>% subset(CHR %in% 1:22)
write.table(eas_cd, paste0(eas, "ibd_EAS_SiKJ_meta_CD_markername_rmdup_chr1.22.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

eas_uc <- fread(paste0(eas,"ibd_EAS_SiKJ_meta_UC_markername_rmdup.txt"), fill = T)   %>% subset(CHR %in% 1:22)
write.table(eas_uc, paste0(eas, "ibd_EAS_SiKJ_meta_UC_markername_rmdup_chr1.22.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
