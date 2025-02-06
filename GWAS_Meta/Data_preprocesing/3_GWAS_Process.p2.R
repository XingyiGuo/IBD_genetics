
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

### Allele frequency/rsid is missing in a gwas meta-analysis, so 1KG was used to map MAF/rsid #########
###### map rsid, allele freq for gwas ########

#### IBDGC
setwd("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2017NG/processed")

# loop through chr 1 to 22
ref_process <- function(chr_id) {
  # read bim files
  path_bim <- sprintf("/nobackup/sbcs/chenz27/For_others/GuoLab/Reference/1KG_hg38/EUR/processed/chr%d.EUR.bim", chr_id)
  
  ref_bim <- fread(path_bim)
  
  header1<- c("CHR",	"SNP",	"GENO",	"BP",	"ALT",	"REF")
  colnames(ref_bim)<- header1
  ref_bim<-  ref_bim %>% select(, -c("GENO"))
  
  ref_bim <- ref_bim %>%
    rowwise() %>%
    mutate(Allele1_sorted = pmin(ALT, REF),
           Allele2_sorted = pmax(ALT, REF)) %>%
    ungroup()
  
  # frq file to get maf
  path_frq <- sprintf("/nobackup/sbcs/chenz27/For_others/GuoLab/Reference/1KG_hg38/EUR/processed/chr%d.EUR.frq", chr_id)
  
  ref_frq <- fread(path_frq)
  
  header2<- c("CHR","SNP","ALT","REF","MAF.EUR","NCHROBS.EUR")
  colnames(ref_frq)<- header2
  
  ref_frq<- ref_frq %>% select(, c("MAF.EUR","NCHROBS.EUR"))
  
  # merged_bim_frq<- merge(ref_bim, ref_frq, by = c("SNP","ALT","REF"), all.x = T, all.y = T, allow.cartesian = T)
  cb_bim_frq<- cbind(ref_bim, ref_frq)
  
  return(cb_bim_frq)
  
}

  
map_ref_to_gwas <- function(gwas, chr_id) {
  
  # read processed ref data for the given chr
  ref <- ref_process(chr_id)
  
  gwas$CHR<- as.integer(gwas$CHR)
  gwas$BP<- as.integer(gwas$BP)  
  
  gwas_chr <- gwas %>%
    filter(CHR == chr_id)
  
  gwas_chr<- gwas_chr %>% mutate(MarkerName_dup = MarkerName)  %>% 
    separate(MarkerName_dup, into = c("CHR", "BP", "Allele1_sorted", "Allele2_sorted"), sep = ":")
  
  gwas_chr <- gwas_chr %>%
    mutate(Allele1 = toupper(Allele1),
           Allele2 = toupper(Allele2))
  
  gwas_chr <- gwas_chr %>%
    rowwise() %>%
    mutate(Allele1_sorted = pmin(Allele1, Allele2),
           Allele2_sorted = pmax(Allele1, Allele2)) %>%
    ungroup()
  
  colnames(gwas_chr)[19]<- "BP"
  
  # merge bim.frq with gwas meta results
  merged_data <- merge(gwas_chr, ref, by = c("CHR", "BP", "Allele1_sorted", "Allele2_sorted"), all.x = T)

  # maf
  merged_data$maf <- ifelse(merged_data$MAF.EUR < 0.5, merged_data$MAF.EUR, 1 - merged_data$MAF.EUR)
  merged_data$CHR<- as.integer(merged_data$CHR)
  merged_data$BP<- as.integer(merged_data$BP)
  
  
  return(merged_data)
  
}

res_process <- function(merged, trait) {

  gwas <- read.table(merged, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # map by chr
  mapped <- lapply(1:22, function(chr_id) map_ref_to_gwas(gwas, chr_id)) %>%
    bind_rows() %>% 
    arrange(CHR, BP)
  
  write.table(mapped, trait, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  return(mapped)
}


res_process("cd_build38_markername_rmdup.txt", "cd_build38_markername_rmdup_maf.txt")
res_process("uc_build38_markername_rmdup.txt", "uc_build38_markername_rmdup_maf.txt")
res_process("ibd_build38_markername_rmdup.txt", "ibd_build38_markername_rmdup_maf.txt")

# remove the X chrs
# IBDGC
ng17<- "/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2017NG/processed/"
ng17_ibd<- fread(paste0(ng17, "ibd_build38_markername_rmdup_maf.txt"), fill = T)   %>% subset(CHR %in% 1:22)
write.table(ng17_ibd, paste0(ng17, "ibd_build38_markername_rmdup_maf_chr1.22.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

ng17_cd<- fread(paste0(ng17, "cd_build38_markername_rmdup_maf.txt"), fill = T)   %>% subset(CHR %in% 1:22)
write.table(ng17_cd, paste0(ng17, "cd_build38_markername_rmdup_maf_chr1.22.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

ng17_uc<- fread(paste0(ng17, "uc_build38_markername_rmdup_maf.txt"), fill = T)  %>% subset(CHR %in% 1:22)
write.table(ng17_uc, paste0(ng17, "uc_build38_markername_rmdup_maf_chr1.22.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

