library(data.table)
library(dplyr)
library(tidyr)


###### merge back heterogenity test info for gwas meta results ##########

setwd("/nobackup/sbcs/lyul1/IBD/metal_update/withoutUKB/gc_0.01/GC_Final/processed")

filt1<- "/nobackup/sbcs/lyul1/IBD/metal_update/withoutUKB/gc_0.01/maf0.01_gc_" 
# cd
cd<-  read.table(paste0(filt1, "EUR_CD_META1.TBL"), header=T, sep="\t", stringsAsFactors = F) %>% select(., c(1,7:11))
# uc
uc<-  read.table(paste0(filt1, "EUR_UC_META1.TBL"), header=T, sep="\t", stringsAsFactors = F) %>% select(., c(1,7:11))
# ibd
ibd<-  read.table(paste0(filt1, "EUR_IBD_META1.TBL"), header=T, sep="\t", stringsAsFactors = F) %>% select(., c(1,7:11))

## EASEUR ##
# cd
eascd<- read.table(paste0(filt1, "EASEUR_CD_META1.TBL"), header = T, sep = "\t", stringsAsFactors = F) %>% select(., c(1, 7:11))
# uc
easuc<- read.table(paste0(filt1, "EASEUR_UC_META1.TBL"), header = T, sep = "\t", stringsAsFactors = F) %>% select(., c(1, 7:11))
# ibd
easibd<- read.table(paste0(filt1, "EASEUR_IBD_META1.TBL"), header = T, sep = "\t", stringsAsFactors = F) %>% select(., c(1, 7:11))


filt<- "/nobackup/sbcs/lyul1/IBD/metal_update/withoutUKB/gc_0.01/GC_Final/final_gc_"
# cd
cd_final<-  read.table(paste0(filt, "EUR_CD_META1.TBL"), header=T, sep="\t", stringsAsFactors = F) %>% select(., -(7:11))
# uc
uc_final<-  read.table(paste0(filt, "EUR_UC_META1.TBL"), header=T, sep="\t", stringsAsFactors = F) %>% select(., -(7:11))
# ibd
ibd_final<-  read.table(paste0(filt, "EUR_IBD_META1.TBL"), header=T, sep="\t", stringsAsFactors = F) %>% select(., -(7:11))

## EASEUR ##
# cd
eascd_final<- read.table(paste0(filt, "EASEUR_CD_META1.TBL"), header = T, sep = "\t", stringsAsFactors = F) %>% select(., -(7:11))
# uc
easuc_final<- read.table(paste0(filt, "EASEUR_UC_META1.TBL"), header = T, sep = "\t", stringsAsFactors = F) %>% select(., -(7:11))
# ibd
easibd_final<- read.table(paste0(filt, "EASEUR_IBD_META1.TBL"), header = T, sep = "\t", stringsAsFactors = F) %>% select(., -(7:11))

merge_gws<- function( trait) {
  
  org<-  read.table(paste0(filt1, trait), header=T, sep="\t", stringsAsFactors = F) %>% select(., c(1,7:11))
  
  final<-  read.table(paste0(filt, trait), header=T, sep="\t", stringsAsFactors = F) %>% select(., -(7:11))
  
  merge<-  merge(final, org, by = "MarkerName")
  
  write.table(merge, paste0(trait, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

}

merge_gws("EASEUR_CD_META1.TBL")
merge_gws("EASEUR_UC_META1.TBL")
merge_gws("EASEUR_IBD_META1.TBL")
merge_gws("EUR_CD_META1.TBL")
merge_gws("EUR_UC_META1.TBL")
merge_gws("EUR_IBD_META1.TBL")
