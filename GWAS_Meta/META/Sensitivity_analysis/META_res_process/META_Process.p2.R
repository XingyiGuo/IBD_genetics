library(data.table)
library(dplyr)
library(tidyr)



###### map rsid, allele freq for GWAS meta results ########

setwd("/nobackup/sbcs/lyul1/IBD/metal_update/withoutUKB/gc_0.01/GC_Final/processed")

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
  
  gwas_chr <- gwas %>%
    mutate(MarkerName_dup = MarkerName) %>%
    separate(MarkerName_dup, into = c("CHR", "BP", "Allele1_sorted", "Allele2_sorted"), sep = ":") %>%
    filter(CHR == chr_id)
  
  # merge bim.frq with gwas meta results
  merged_data <- merge(gwas_chr, ref, by = c("CHR", "BP", "Allele1_sorted", "Allele2_sorted"), all.x = T)
  
  merged_data<- merged_data %>% select(, c("CHR", "BP", "ALT", "REF", "SNP", "MarkerName", "Allele1", 
                                           "Allele2", "Effect", "StdErr", "P.value", "Direction", 
                                           "HetISq",	"HetChiSq",	"HetDf",	"HetPVal", "MAF.EUR", "NCHROBS.EUR")) 

  
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

res_process("EASEUR_CD_META1.TBL.txt", "EASEUR_CD_META1.TBL.SNP.txt")
res_process("EASEUR_UC_META1.TBL.txt", "EASEUR_UC_META1.TBL.SNP.txt")
res_process("EASEUR_IBD_META1.TBL.txt", "EASEUR_IBD_META1.TBL.SNP.txt")

res_process("EUR_CD_META1.TBL.txt", "EUR_CD_META1.TBL.SNP.txt")
res_process("EUR_UC_META1.TBL.txt", "EUR_UC_META1.TBL.SNP.txt")
res_process("EUR_IBD_META1.TBL.txt", "EUR_IBD_META1.TBL.SNP.txt")

