library(data.table)
library(dplyr)
library(tidyr)

########## format chr:pos for liftover ucsc online conversion ############
# IBDGC GWAS SS
setwd("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2017NG")
ng17<- "/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2017NG/processed/"

format_pos<- function(trait) {

  data<- fread(paste0(trait, "_20161107.txt"), header = T)
  
  data <- data %>%
   separate(MarkerName, into = c("chr_pos", "alleles"), sep = "_") %>%
     separate(chr_pos, into = c("chr", "pos"), sep = ":")
  data<- data[,-3]
  data$chrpos <- paste0("chr", cd$chr, ":",
                       format(cd$pos, big.mark = ",", scientific = FALSE),
                       "-",
                       format(cd$pos, big.mark = ",", scientific = FALSE))

     write.table(data, paste0(ng17, "chrpos_",trait, ".txt"), sep = "\t", row.names = FALSE, col.names = F, quote = FALSE)

   }

format_pos("cd_build37_40266")
format_pos("uc_build37_45975")
format_pos("ibd_build37_59957")



# Pan-UKB GWAS SS
setwd("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_Pan_UKB")
ukb<- "/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_Pan_UKB/processed/liftover/"

format_pos<- function(trait) {
  
  data<- fread(paste0("phecode-", trait,"-both_sexes.tsv"), header = T)
  
  data<- data %>% filter( chr %in% 1:22)
  
  data<- data[, c("chr", "pos")]
  data$chrpos <- paste0("chr", data$chr, ":", data$pos, "-", data$pos)
  
  data<- data[,-(1:2)]
  
  write.table(data, paste0(ukb, "chrpos_",trait, ".txt"), sep = "\t", row.names = FALSE, col.names = F, quote = FALSE)
  
  rws<- nrow(data)
  #r1<- nrow(data)/2
  #r2<- 2* nrow(data)/3 
  
  data1<- data[1:14029018, ]
  data2<- data[14029019:28058035, ]
  
  write.table(data1, paste0(ukb, "chrpos_",trait, "_p1.txt"), sep = "\t", row.names = FALSE, col.names = F, quote = FALSE)
  write.table(data2, paste0(ukb, "chrpos_",trait, "_p2.txt"), sep = "\t", row.names = FALSE, col.names = F, quote = FALSE)  
  
}


format_pos("555_IBD")
format_pos("555.1_CD")
format_pos("555.2_UC")
