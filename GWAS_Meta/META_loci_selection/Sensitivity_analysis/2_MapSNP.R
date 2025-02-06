library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)
library(readxl)


######## map rsids, alt, ref alleles, and maf to gwas meta significant snps ######
setwd("/nobackup/sbcs/lyul1/IBD/metal_update/withoutUKB/gc_0.01")

####### read meta significant results #######
# eur
sig_cd <- read_excel("./select/final_gc_0.01_EUR_CD_5e-08.xlsx", sheet = 1) %>% as.data.frame() # final_gc_0.01_EUR_CD_5e-08.xlsx
sig_uc <- read_excel("./select/final_gc_0.01_EUR_UC_5e-08.xlsx", sheet = 1) %>% as.data.frame()
sig_ibd <- read_excel("./select/final_gc_0.01_EUR_IBD_5e-08.xlsx", sheet = 1) %>% as.data.frame()
# eas eur
sig_eascd <- read_excel("./select/final_gc_0.01_EASEUR_CD_5e-08.xlsx", sheet = 1) %>% as.data.frame()
sig_easuc <- read_excel("./select/final_gc_0.01_EASEUR_UC_5e-08.xlsx", sheet = 1) %>% as.data.frame()
sig_easibd <- read_excel("./select/final_gc_0.01_EASEUR_IBD_5e-08.xlsx", sheet = 1) %>% as.data.frame()


############ map rsid, ref, alt alleles, freq in EUR ref data ##############

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
  
  # filter GWAS data for the current chr
  gwas <- gwas %>%
    mutate(Allele1 = toupper(Allele1),
           Allele2 = toupper(Allele2))
  
  gwas <- gwas %>%
    rowwise() %>%
    mutate(Allele1_sorted = pmin(Allele1, Allele2),
           Allele2_sorted = pmax(Allele1, Allele2)) %>%
    ungroup()

  # colnames(gwas)[13]<- "BP"
  
  gwas_chr <- gwas %>%
    filter(CHR == chr_id)
  
  # merge bim.frq with gwas meta results
  merged_data <- merge(gwas_chr, ref, by = c("CHR", "BP", "Allele1_sorted", "Allele2_sorted"), all.x = T)
  
  merged_data<- merged_data %>% select(, c("CHR", "BP", "ALT", "REF", "SNP", "MarkerName", "Allele1", 
                                           "Allele2", "Effect", "StdErr", "P.value",  "Direction",
                                           "HetISq",	"HetChiSq",	"HetDf",	"HetPVal", "MAF.EUR", "NCHROBS.EUR"))
  
  return(merged_data)
  
}

res_process <- function(merged, trait) {

  gwas<- merged
  
  # map by chr
  mapped <- lapply(1:22, function(chr_id) map_ref_to_gwas(gwas, chr_id)) %>%
    bind_rows() %>% 
    arrange(CHR, BP)
  
  write.csv(mapped, trait, row.names = F)
  
  return(mapped)

}


mapped_cd<- res_process(sig_cd, "./select/final_gc_0.01_EUR_CD_5e-08.csv")
mapped_uc<- res_process(sig_uc, "./select/final_gc_0.01_EUR_UC_5e-08.csv")
mapped_ibd<- res_process(sig_ibd, "./select/final_gc_0.01_EUR_IBD_5e-08.csv")

mapped_eascd <- res_process(sig_eascd, "./select/final_gc_0.01_EASEUR_CD_5e-08_eurfrq.csv")
mapped_easuc <- res_process(sig_easuc,  "./select/final_gc_0.01_EASEUR_UC_5e-08_eurfrq.csv")
mapped_easibd <- res_process(sig_easibd, "./select/final_gc_0.01_EASEUR_IBD_5e-08_eurfrq.csv")



############ map rsid, ref, alt alleles, freq in EAS ref data ##############

# loop through chr 1 to 22
ref_process <- function(chr_id) {
  # read bim files
  path_bim <- sprintf("/nobackup/sbcs/chenz27/For_others/GuoLab/Reference/1KG_hg38/EAS/processed/chr%d.EAS.bim", chr_id)
  
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
  path_frq <- sprintf("/nobackup/sbcs/chenz27/For_others/GuoLab/Reference/1KG_hg38/EAS/processed/chr%d.EAS.frq", chr_id)
  
  ref_frq <- fread(path_frq)
  
  header2<- c("CHR","SNP","ALT","REF","MAF.EAS","NCHROBS.EAS")
  colnames(ref_frq)<- header2
  
  ref_frq<- ref_frq %>% select(, c("MAF.EAS","NCHROBS.EAS"))

  cb_bim_frq<- cbind(ref_bim, ref_frq)
  
  return(cb_bim_frq)
  
}


map_ref_to_gwas <- function(gwas, chr_id) {
  
  # read processed ref data for the given chr
  ref <- ref_process(chr_id)
  
  # filter GWAS data for the current chr
  gwas <- gwas %>%
    mutate(Allele1 = toupper(Allele1),
           Allele2 = toupper(Allele2))
  
  gwas <- gwas %>%
    rowwise() %>%
    mutate(Allele1_sorted = pmin(Allele1, Allele2),
           Allele2_sorted = pmax(Allele1, Allele2)) %>%
    ungroup()
  
  colnames(gwas)[13]<- "BP"
  
  gwas_chr <- gwas %>%
    filter(CHR == chr_id)
  
  # merge bim.frq with gwas meta results
  merged_data <- merge(gwas_chr, ref, by = c("CHR", "BP", "Allele1_sorted", "Allele2_sorted"), all.x = T)
  
  merged_data<- merged_data %>% select(, c("CHR", "BP", "ALT", "REF", "SNP", "MarkerName", "Allele1", 
                                           "Allele2", "Effect", "StdErr", "P.value", "Direction", 
                                           "HetISq",	"HetChiSq",	"HetDf",	"HetPVal","MAF.EAS", "NCHROBS.EAS"))
  
  return(merged_data)
  
}

res_process <- function(merged, trait) {
  
  gwas<- merged
  
  # map by chr
  mapped <- lapply(1:22, function(chr_id) map_ref_to_gwas(gwas, chr_id)) %>%
    bind_rows() %>% 
    arrange(CHR, BP)
  
  write.csv(mapped, trait, row.names = F)
  
  return(mapped)
  
}

mapped_eascd <- res_process(sig_eascd, "./select/final_gc_0.01_EASEUR_CD_5e-08_easfrq.csv")
mapped_easuc <- res_process(sig_easuc,  "./select/final_gc_0.01_EASEUR_UC_5e-08_easfrq.csv")
mapped_easibd <- res_process(sig_easibd, "./select/final_gc_0.01_EASEUR_IBD_5e-08_easfrq.csv")
