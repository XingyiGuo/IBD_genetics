library(data.table)
library(dplyr)
library(tidyr)

########## Liftover grch37->38 (has been done using online tools ############
### used liftover online tool to get the successful conversion and failed output file
### map the grch 38 with original gwas ss #######

# select the success conversion
setwd("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2017NG")
ng17<- "/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2017NG/processed/"

select_success<- function(uc, df, success, trait) {
  uc <- uc %>%
    mutate(MarkerName_dupl = MarkerName) %>%
    separate(MarkerName_dupl, into = c("chr_pos", "alleles"), sep = "_")
  
  colnames(df)[1]<- "chr_pos"
  df <- df %>%
    separate(chr_pos, into = c("chr_pos", "null"), sep = "-") %>%
    select(chr_pos) %>%
    mutate(chr_pos = gsub("^chr", "", chr_pos))
  
  matching_index<- which(uc$chr_pos %in% df$chr_pos)
  uc_filtered <- uc[-matching_index, ] %>% as.data.frame()

  colnames(success)[1]<- "hg38chr_pos"
  combined<- cbind(uc_filtered, success) 
  
  combined <- combined %>%
    mutate(
      CHR = as.integer(sub("chr(.*):.*", "\\1", hg38chr_pos)),
      POS = as.integer(sub(".*:(.*)-.*", "\\1", hg38chr_pos))
    )
  
  combined <- combined %>%
    select(-hg38chr_pos)
  
  write.table(combined, paste0(ng17, trait,"_build38.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
}

ibd_list=c("cd_build37_40266_20161107.txt","ibd_build37_59957_20161107.txt", "uc_build37_45975_20161107.txt")
ng17_cd <-  fread(ibd_list[1])
failed_cd <- read.table("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2017NG/processed/liftover/cd_failed.txt", header = FALSE, sep = "\t")
success_cd<- read.table("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2017NG/processed/liftover/cd_success_hglft_genome_3910e_571cb0.bed", header = F, sep = "\t", stringsAsFactors = F) %>% as.data.frame() 

ng17_uc <-  fread(ibd_list[3])
failed_uc <- read.table("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2017NG/processed/liftover/uc_failed.txt", header = FALSE, sep = "\t")
success_uc<- read.table("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2017NG/processed/liftover/uc_success_hglft_genome_3cc68_574470.bed", header = F, sep = "\t", stringsAsFactors = F) %>% as.data.frame() 

ng17_ibd <- fread(ibd_list[2])
failed_ibd <- read.table("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2017NG/processed/liftover/ibd_failed.txt", header = FALSE, sep = "\t")
success_ibd<- read.table("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2017NG/processed/liftover/ibd_success_hglft_genome_2e90a_592350.bed", header = F, sep = "\t", stringsAsFactors = F) %>% as.data.frame() 

select_success(ng17_cd, failed_cd, success_cd, "cd")
select_success(ng17_uc, failed_uc, success_uc, "uc")
select_success(ng17_ibd, failed_ibd, success_ibd, "ibd")

###### Pan-UKB
setwd("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_Pan_UKB/processed/liftover/")

# select the success conversion
select_success<- function(trait, part) {
  
  df<- fread(paste0("chrpos_", trait, "_", part, ".bed"), header = F) # success converted: grch 37 chrpos38_555.1_CD_p1.bed
  
  colnames(df)[1]<- "hg38chr_pos"
  df <- df %>%
    separate(hg38chr_pos, into = c("hg38chr_pos", "null"), sep = "-") %>%
    select(hg38chr_pos) %>%
    mutate(hg38chr_pos = gsub("^chr", "", hg38chr_pos))
  
  original<- fread(paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_Pan_UKB/processed/", trait, "_", part, ".txt"), header = T)
  original$chr_pos <- paste0("chr", original$chr, ":", original$pos, "-", original$pos)
  
  failed<- fread(paste0("chrpos_", trait, "_", part, "_failed38.txt"), header = F)
  colnames(failed)[1]<- "chr_pos"
  
  matching_index<- which(original$chr_pos %in% failed$chr_pos) # grch38
  success_conv_grch37 <- original[-matching_index, ] %>% as.data.frame()
  
  colnames(success_conv_grch37)[11]<- "hg37chr_pos"
  combined<- cbind(df, success_conv_grch37) 
  
  combined <- combined %>%
    select(-hg37chr_pos)
  
  combined <- combined %>% separate(hg38chr_pos, into = c("CHR", "BP"), sep = ":")
  
  combined$chr<- as.integer(combined$chr)
  combined$BP<- as.integer(combined$BP)
  
  write.table(combined, paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_Pan_UKB/processed/", trait, "_", part, "_build38.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  return(combined)
}

cd_381<- select_success("555.1_CD", "p1")
uc_381<- select_success("555.2_UC", "p1")
ibd_381<- select_success("555_IBD", "p1")

cd_382<- select_success("555.1_CD", "p2")
uc_382<- select_success("555.2_UC", "p2")
ibd_382<- select_success("555_IBD", "p2")


#### combine two parts of gwas into a single file
read_38<- function(trait, part) {
  
  data<- fread(paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_Pan_UKB/processed/", trait, "_", part, "_build38.txt"), header = T)
  
  return(data)
}

cd_381<- read_38("555.1_CD", "p1")
uc_381<- read_38("555.2_UC", "p1")
ibd_381<- read_38("555_IBD", "p1")

cd_382<- read_38("555.1_CD", "p2")
uc_382<- read_38("555.2_UC", "p2")
ibd_382<- read_38("555_IBD", "p2")


cd_38<- rbind(cd_381, cd_382)
uc_38<- rbind(uc_381, uc_382)
ibd_38<- rbind(ibd_381, ibd_382)

cd_38<- filter(cd_38, cd_38$CHR != "CHR")
uc_38<- filter(uc_38, uc_38$CHR != "CHR")
ibd_38<- filter(ibd_38, ibd_38$CHR != "CHR")

write.table(cd_38, paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_Pan_UKB/processed/", "555.1_CD", "_build38.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(uc_38, paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_Pan_UKB/processed/", "555.2_UC", "_build38.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(ibd_38, paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_Pan_UKB/processed/", "555_IBD", "_build38.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



