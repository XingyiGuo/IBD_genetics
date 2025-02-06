library(data.table)
library(dplyr)
library(tidyr)

########## format original gwas ss (calculate beta se, maf / create markernames) for meta analysis #################

# EUR IBD
# mvp
mvp<- "/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2023MVP/"
mvp_cd<- fread(paste0(mvp,"MVP_R4.1000G_AGR.Phe_555_1.EUR.GIA.dbGaP.txt"), fill = T) # Phe_555_1	binary	AfAdE	PheCodes	Digestive System	Regional enteritis
mvp_uc<- fread(paste0(mvp,"MVP_R4.1000G_AGR.Phe_555_2.EUR.GIA.dbGaP.txt"), fill = T) # Phe_555_2	binary	AfAdE	PheCodes	Digestive System	Ulcerative colitis
mvp_ibd<- fread(paste0(mvp,"MVP_R4.1000G_AGR.Phe_555.EUR.GIA.dbGaP.txt"), fill = T)  # Phe_555	binary	AfAdE	PheCodes	Digestive System	Inflammatory bowel disease and other gasteroenteritis and colitis


mvp_process <- function(data, trait) {
  # ci
  data <- data %>% 
    separate(ci, into = c("lci", "uci"), sep = ",")
  # beta se
  data$beta <- log(data$or)
  # data$se <- abs(log(data$or) / qnorm(data$pval / 2)) # se=abs(log(or)/qnorm(p/2)) 
  data$se <- (log(as.numeric(data$uci)) - log(as.numeric(data$lci))) / (2 * 1.96)
  # nea
  data$nea <- ifelse(data$alt == data$ea, data$ref, data$alt)
  # maf
  data$maf <- ifelse(data$af < 0.5, data$af, 1 - data$af)
  
  # markernames (chr:pos:a1_sorted:a2_sorted)  
  data <- data %>%
    rowwise() %>%
    mutate(Allele1_sorted = pmin(ref, alt),
           Allele2_sorted = pmax(ref, alt)) %>%
    ungroup()
  
  data$pos <- format(as.integer(as.numeric(data$pos)), scientific = FALSE)
  data$pos<- as.integer(data$pos)
  data$chrom<- as.integer(data$chrom)
  
  data<- data %>% mutate(
    markernames = paste0(chrom, ":", pos, ":", Allele1_sorted, ":", Allele2_sorted)
  )
  
  write.table(data, paste0(mvp, "processed/", trait, "_markername.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  # in case duplicate snps, remove the one with higher p.val
  data_rmdup <- data %>%
    arrange(markernames, pval) %>%
    group_by(markernames) %>%
    slice(1) %>%
    ungroup()
  
  write.table(data_rmdup, paste0(mvp, "processed/", trait, "_markername_rmdup.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

mvp_process(mvp_cd, "Phe_555_1.EUR.CD")
mvp_process(mvp_uc, "Phe_555_2.EUR.UC")
mvp_process(mvp_ibd, "Phe_555.EUR.IBD")

rm(mvp_cd)
rm(mvp_uc)
rm(mvp_ibd)

# finn
finn<- "/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_FinnGen/R12/"
header<- c("chrom","pos", "ref","alt","rsids", "nearest_genes", "pval","mlogp", "beta","sebeta","af_alt","af_alt_cases","af_alt_controls")
finn_cd<- read.table(paste0(finn, "finngen_R12_K11_CD_STRICT2"),header=F, sep="\t", stringsAsFactors = F)  %>%  setNames(header)
finn_uc<- read.table(paste0(finn, "finngen_R12_K11_UC_STRICT2"),header=F, sep="\t", stringsAsFactors = F) %>%  setNames(header)
finn_ibd<- read.table(paste0(finn, "finngen_R12_K11_IBD_STRICT"),header=F, sep="\t", stringsAsFactors = F)  %>%  setNames(header)
# calculate N
# cd
finn_cd <- finn_cd %>%
  mutate(
    num_cases = 2489,
    num_controls = 497622,
    num_samples = num_cases + num_controls
  )
# uc
finn_uc <- finn_uc %>%
  mutate(
    num_cases = 7220,
    num_controls = 492160,
    num_samples = num_cases + num_controls
  )
# ibd
finn_ibd <- finn_ibd %>%
  mutate(
    num_cases = 10960,
    num_controls = 489388,
    num_samples = num_cases + num_controls
  )

# maf + markername
finn_process <- function(data, trait) {
  
  data$rsids[data$rsids == ""]<-NA # if "", metal will ignore this cell and shift cells left, so replace it with NA
  
  # maf
  data$maf_alt <- ifelse(data$af_alt < 0.5, data$af_alt, 
                         1 - data$af_alt)
  
  data$maf_alt_controls <- ifelse(data$af_alt_controls < 0.5, data$af_alt_controls, 
                                  1 - data$af_alt_controls)
  
  # markername (chr:pos:a1_sorted:a2_sorted)
  data <- data %>%
    rowwise() %>%
    mutate(Allele1_sorted = pmin(ref, alt),
           Allele2_sorted = pmax(ref, alt)) %>%
    ungroup()
  
  data$pos <- format(as.integer(as.numeric(data$pos)), scientific = FALSE)
  data$pos<- as.integer(data$pos)
  data$chrom<- as.integer(data$chrom)
  
  data<- data %>% mutate(
    markernames = paste0(chrom, ":", pos, ":", Allele1_sorted, ":", Allele2_sorted)
  )
  
  write.table(data, paste0(finn,"processed/", trait, "_markername.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  data_rmdup <- data %>%
    arrange(markernames, pval) %>%
    group_by(markernames) %>%
    slice(1) %>%
    ungroup()
  
  write.table(data_rmdup, paste0(finn, "processed/", trait, "_markername_rmdup.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

  data_rmdup<- data_rmdup %>% filter( CHR %in% 1:22)
  
  write.table(data_rmdup, paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_FinnGen/R12/processed/", trait, "_markername_rmdup_chr1.22.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

}

finn_process(finn_cd, "finngen_R12_K11_CD_STRICT2")
finn_process(finn_uc, "finngen_R12_K11_UC_STRICT2")
finn_process(finn_ibd, "finngen_R12_K11_IBD_STRICT")

rm(finn_ibd)
rm(finn_uc)
rm(finn_cd)

# ng17
# read in grch 38 converted ng17 gwas ss ### this processing code was provided in liftover.R file
ng17<- "/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2017NG/processed/" 
ng17_cd<- fread(paste0(ng17, "cd_build38.txt"), fill = T)
ng17_uc<- fread(paste0(ng17, "uc_build38.txt"), fill = T)
ng17_ibd<- fread(paste0(ng17, "ibd_build38.txt"), fill = T)

ng17_cd <- ng17_cd %>%
  mutate(
    num_cases = 12194,
    num_controls = 28072,
    num_samples = num_cases + num_controls
  )

# uc
ng17_uc <- ng17_uc %>%
  mutate(
    num_cases = 12366,
    num_controls = 33609,
    num_samples = num_cases + num_controls
  )

# ibd
ng17_ibd <- ng17_ibd %>%
  mutate(
    num_cases = 25042	,
    num_controls = 34915,
    num_samples = num_cases + num_controls
  )

ng17_process <- function(data, trait) {
  
  data <- data %>%
    mutate(Allele1 = toupper(Allele1),
           Allele2 = toupper(Allele2))
  
  data <- data %>%
    rowwise() %>%
    mutate(Allele1_sorted = pmin(Allele1, Allele2),
           Allele2_sorted = pmax(Allele1, Allele2)) %>%
    ungroup()
  
  data$POS <- format(as.integer(as.numeric(data$POS)), scientific = FALSE)
  data$POS<- as.integer(data$POS)
  data$CHR<- as.integer(data$CHR)
  
  # markername (chr:pos:a1_sorted:a2_sorted)
  data<- data %>% mutate(
    markernames = paste0(CHR, ":", POS, ":", Allele1_sorted, ":", Allele2_sorted)
  )
  
  write.table(data, paste0(ng17, trait, "_markername.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  data_rmdup <- data %>%
    arrange(markernames, P.value) %>%
    group_by(markernames) %>%
    slice(1) %>%
    ungroup()
  
  write.table(data_rmdup, paste0(ng17, trait, "_markername_rmdup.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

}

ng17_process(ng17_cd, "cd_build38")
ng17_process(ng17_uc, "uc_build38")
ng17_process(ng17_ibd, "ibd_build38")

rm(ng17_cd)
rm(ng17_uc)
rm(ng17_ibd)

# Pan-UKB
cd_38<- fread(paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_Pan_UKB/processed/", "555.1_CD", "_build38.txt"), header = T)
uc_38<- fread(paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_Pan_UKB/processed/", "555.2_UC", "_build38.txt"), header = T)
ibd_38<- fread(paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_Pan_UKB/processed/", "555_IBD", "_build38.txt"), header = T)

cd_38 <- cd_38 %>%
  mutate(
    n_cases_EUR = 2008,
    n_controls_EUR = 330270,
    n_samples_EUR = n_cases_EUR + n_controls_EUR
  )

# uc
uc_38 <- uc_38 %>%
  mutate(
    n_cases_EUR= 3670,
    n_controls_EUR= 330270,
    n_samples_EUR = n_cases_EUR + n_controls_EUR
  )

# ibd
ibd_38 <- ibd_38 %>%
  mutate(
    n_cases_EUR = 5183, 
    n_controls_EUR = 330270,
    n_samples_EUR = n_cases_EUR + n_controls_EUR
  )

ukb_process <- function(data, trait) {
  
  data$P.value<- 10^(-(data$neglog10_pval_EUR))
  
  # maf
  data$maf_alt_cases <- ifelse(data$af_cases_EUR < 0.5, data$af_cases_EUR, 
                               1 - data$af_cases_EUR)
  
  data$maf_alt_controls <- ifelse(data$af_controls_EUR < 0.5, data$af_controls_EUR, 
                                  1 - data$af_controls_EUR)
  
  data$maf_alt_est<- (data$maf_alt_cases + data$maf_alt_controls) / 2
  
  data <- data %>%
    rowwise() %>%
    mutate(Allele1_sorted = pmin(ref, alt),
           Allele2_sorted = pmax(ref, alt)) %>%
    ungroup()
  
  data$BP <- format(as.integer(as.numeric(data$BP)), scientific = FALSE)
  data$BP<- as.integer(data$BP)
  data$CHR<- as.integer(data$CHR)
  
  #markername (chr:pos:a1_sorted:a2_sorted)
  data<- data %>% mutate(
    markernames = paste0(CHR, ":", BP, ":", Allele1_sorted, ":", Allele2_sorted)
  )
  
  write.table(data, paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_Pan_UKB/processed/", trait, "_markername.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  data_rmdup <- data %>%
    arrange(markernames, P.value) %>%
    group_by(markernames) %>%
    slice(1) %>%
    ungroup()
  
  write.table(data_rmdup, paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_Pan_UKB/processed/", trait, "_markername_rmdup.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  data_rmdup<- data_rmdup %>% filter( CHR %in% 1:22)
  
  write.table(data_rmdup, paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_Pan_UKB/processed/", trait, "_markername_rmdup_chr1.22.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
}

ukb_process(cd_38, "555.1_CD_build38")
ukb_process(uc_38, "555.2_UC_build38")
ukb_process(ibd_38, "555_IBD_build38")


# east asian IBD | Allele1: reference allele Allele2: effect allele
eas<- "/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2023NG/liu-2022-east-asian-gwas/summary-stats/"
# EAS only GWAS
eas_cd<- fread(paste0(eas, "ibd_EAS_SiKJ_meta_CD.TBL.txt"), fill = T) 
eas_uc <- fread(paste0(eas,"ibd_EAS_SiKJ_meta_UC.TBL.txt"), fill = T) 
eas_ibd <- fread(paste0(eas,"ibd_EAS_SiKJ_meta_IBD.TBL.txt"), fill = T) 
# EAS EUR only GWAS
eas_eur_cd <- fread(paste0(eas,"ibd_EAS_EUR_SiKJEF_meta_CD.TBL.txt.gz"), fill = T) %>% select(, c(1:17))
eas_eur_uc <- fread(paste0(eas,"ibd_EAS_EUR_SiKJEF_meta_UC.TBL.txt.gz"), fill = T) %>% select(, c(1:17))
eas_eur_ibd <- fread(paste0(eas,"ibd_EAS_EUR_SiKJEF_meta_IBD.TBL.txt.gz"), fill = T) %>% select(, c(1:17))

# cd
eas_cd <- eas_cd %>%
  mutate(
    num_cases = 7372,
    num_controls = 15456,
    num_samples = num_cases + num_controls
  )
# cd eas eur
eas_eur_cd <- eas_eur_cd %>%
  mutate(
    num_cases = 20873,
    num_controls = 346719,
    num_samples = num_cases + num_controls
  )

# uc
eas_uc <- eas_uc %>%
  mutate(
    num_cases = 6862 ,
    num_controls = 15456,
    num_samples = num_cases + num_controls
  )
# eas_eur_uc
eas_eur_uc <- eas_eur_uc %>%
  mutate(
    num_cases = 23252 ,
    num_controls = 352256,
    num_samples = num_cases + num_controls
  )

# ibd
eas_ibd <- eas_ibd %>%
  mutate(
    num_cases = 14393 ,
    num_controls = 15456,
    num_samples = num_cases + num_controls
  )
# ibd eas eur
eas_eur_ibd <- eas_eur_ibd %>%
  mutate(
    num_cases = 45106 ,
    num_controls = 353562,
    num_samples = num_cases + num_controls
  )
# 20873 CD and 346719 Control ibd_EAS_EUR_SiKJEF_meta_CD.TBL.txt.gz # 23252 UC and 352256 Control ibd_EAS_EUR_SiKJEF_meta_UC.TBL.txt.gz # 45106 IBD and 353562 Control ibd_EAS_EUR_SiKJEF_meta_IBD.TBL.txt.gz

eas_process <- function(data, trait) {
  colnames(data)[10]<- "P.value"
  # maf
  data$maf <- ifelse(data$Freq1 < 0.5, data$Freq1, 
                     1 - data$Freq1)
  
  # markername (chr:pos:a1_sorted:a2_sorted)
  data <- data %>%
    mutate(Allele1 = toupper(Allele1),
           Allele2 = toupper(Allele2))
  
  data <- data %>%
    rowwise() %>%
    mutate(Allele1_sorted = pmin(Allele1, Allele2),
           Allele2_sorted = pmax(Allele1, Allele2)) %>%
    ungroup()
  
  data$BP <- format(as.integer(as.numeric(data$BP)), scientific = FALSE)
  data$BP<- as.integer(data$BP)
  data$CHR <- gsub("X", "23", data$CHR)
  data$CHR<- as.integer(data$CHR)
  
  data<- data %>% mutate(
    markernames = paste0(CHR, ":", BP, ":", Allele1_sorted, ":", Allele2_sorted)
  )
  
  write.table(data, paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2023NG/processed/", trait, "_markername.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  data_rmdup <- data %>%
    arrange(markernames, P.value) %>%
    group_by(markernames) %>%
    slice(1) %>%
    ungroup()
  
  write.table(data_rmdup, paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2023NG/processed/", trait, "_markername_rmdup.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  data_rmdup<- data_rmdup %>% subset( CHR %in% 1:22)
  
  write.table(data_rmdup, paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2023NG/processed/", trait, "_markername_rmdup_chr1.22.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

}


eas_process(eas_cd, "ibd_EAS_SiKJ_meta_CD")
eas_process(eas_uc, "ibd_EAS_SiKJ_meta_UC")
eas_process(eas_ibd, "ibd_EAS_SiKJ_meta_IBD")
eas_process(eas_eur_cd, "ibd_EAS_EUR_SiKJEF_meta_CD")
eas_process(eas_eur_uc, "ibd_EAS_EUR_SiKJEF_meta_UC")
eas_process(eas_eur_ibd, "ibd_EAS_EUR_SiKJEF_meta_IBD")


