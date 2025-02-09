library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)
library(readxl)

##### applying additional filters (eg., distance, MAF) and selecting loci; format data for clumping; selecting lead variants from clumped results #####

setwd("/nobackup/sbcs/lyul1/IBD/metal_update/withoutUKB/gc_0.01")
filt<- "./GC_Final/processed/"

####### read in previous gwas significant snps #######  ***** please note here, the reference previous snp list for EUREAS meta are based on 3EUR gwas, 1EAS only gwas, and 1EASEUR gwas meta.
kn_eascd_sig<- fread("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/EASEUR_cd_5e-8_nodup.txt") %>% subset(CHR %in% 1:22)
colnames(kn_eascd_sig)[9]<- "MarkerName"
kn_easuc_sig<- fread("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/EASEUR_uc_5e-8_nodup.txt") %>% subset(CHR %in% 1:22)
colnames(kn_easuc_sig)[9]<- "MarkerName"

kn_easibd_sig<- fread("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/EASEUR_ibd_5e-8_nodup.txt") %>% subset(CHR %in% 1:22)
kn_easibd_all_sig<- fread("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_sign/kn_easibd_all_sig_5e-8_nodup.txt") %>% subset(CHR %in% 1:22)
colnames(kn_easibd_all_sig)[9]<- "MarkerName"
####### read processed meta results #######
## EUR ## Allele 1 = EA; Allele 2 = NEA
# cd
# cd
cd<-  read.csv("./select/final_gc_0.01_EUR_CD_5e-08.csv",header=T) 
# uc
uc<-  read.csv("./select/final_gc_0.01_EUR_UC_5e-08.csv",header=T)
# ibd
ibd<-  read.csv("./select/final_gc_0.01_EUR_IBD_5e-08.csv",header=T) # bc original data has chr 23, while ref data dont have chr, so the number of snp in csv file is smaller.
# ibd1 <- ibd %>% distinct(MarkerName, .keep_all = TRUE)

## EASEUR ## merge eas & eur minor allele frequency
# cd
eufrq_cd<- read.csv("./select/final_gc_0.01_EASEUR_CD_5e-08_eurfrq.csv", header = T)
eafrq_cd<- read.csv("./select/final_gc_0.01_EASEUR_CD_5e-08_easfrq.csv", header = T)
eascd<- merge(eufrq_cd, eafrq_cd[, c(6,17,18)], by = "MarkerName")
#write.csv(eascd, "./select/final_gc_0.01_EASEUR_CD_5e-08_eureasfrq.csv", row.names = F)

# uc
eufrq_uc<- read.csv("./select/final_gc_0.01_EASEUR_UC_5e-08_eurfrq.csv", header = T)
eafrq_uc<- read.csv("./select/final_gc_0.01_EASEUR_UC_5e-08_easfrq.csv", header = T)
easuc<- merge(eufrq_uc, eafrq_uc[, c(6,17,18)], by = "MarkerName")
#write.csv(easuc, "./select/final_gc_0.01_EASEUR_UC_5e-08_eureasfrq.csv", row.names = F)

# ibd
eufrq_ibd<- read.csv("./select/final_gc_0.01_EASEUR_IBD_5e-08_eurfrq.csv", header = T)
eafrq_ibd<- read.csv("./select/final_gc_0.01_EASEUR_IBD_5e-08_easfrq.csv", header = T)
easibd<- merge(eufrq_ibd, eafrq_ibd[, c(6,17,18)], by = "MarkerName")
#write.csv(easibd, "./select/final_gc_0.01_EASEUR_IBD_5e-08_eureasfrq.csv", row.names = F)


select_sig<- function(meta, known_snp) {
  
  # select 1: significant snps (5e-08) in gwas meta (have been done - res in cd/uc/ibd)
  
  # select 2: significant snps (5e-08) in gwas meta but not in previous gawss
  meta_sig1<- meta
  
  meta_sig2 <- anti_join(meta_sig1, known_snp, by = "MarkerName")
  
  return(list(sig1 = meta_sig1, sig2 = meta_sig2))
}

sig_cd<- select_sig(cd, kn_eascd_sig)
sig_uc<- select_sig(uc, kn_easuc_sig)
sig_ibd<- select_sig(ibd, kn_easibd_all_sig)

sig_eascd<- select_sig(eascd, kn_eascd_sig)
sig_easuc<- select_sig(easuc, kn_easuc_sig)
sig_easibd<- select_sig(easibd, kn_easibd_all_sig)


# select 3: check whether meta.sign.snps are within 500kb of any SNP in gwas sig.snps
window <- 500000
select_sig_500kb <- function(meta_sig, known_snp) {
  
  # check whether meta.sign.snps are within 500kb of any SNP in gwas sig.snps
  meta_sig$within_500kb <- apply(meta_sig, 1, function(row) {
    
    any(abs(as.integer(known_snp$POS) - as.integer(row["BP"])) <= 500000 & as.integer(known_snp$CHR) == as.integer(row["CHR"]))
  })
  
  
  # check whether meta.sign.snps are within 1000kb of any SNP in gwas sig.snps
  meta_sig$within_1mb <- apply(meta_sig, 1, function(row) {
    
    any(abs(as.integer(known_snp$POS) - as.integer(row["BP"])) <= 1000000 & as.integer(known_snp$CHR) == as.integer(row["CHR"]))
  })
  
  # obtained common and rare variants based on MAF(EUR) and MAF.EAS
  if ("MAF.EUR" %in% colnames(meta_sig)) {
    meta_sig$common <- ifelse(meta_sig$MAF.EUR > 0.01, "common", "rare")
  }
  
  if ("MAF.EAS" %in% colnames(meta_sig)) {
    meta_sig$common <- ifelse(meta_sig$MAF.EUR > 0.01, "common", "rare")
    meta_sig$common.EAS <- ifelse(meta_sig$MAF.EAS > 0.01, "common", "rare")
    meta_sig$common.EUREAS <- ifelse(meta_sig$MAF.EAS > 0.01 & meta_sig$MAF.EUR > 0.01, "common", "rare")
  }
  
  return(meta_sig)
}

sig_cd[["sig2"]] <- select_sig_500kb(sig_cd[["sig2"]], kn_eascd_sig)

sig_uc[["sig2"]] <- select_sig_500kb(sig_uc[["sig2"]], kn_easuc_sig)

sig_ibd[["sig2"]] <- select_sig_500kb(sig_ibd[["sig2"]], kn_easibd_all_sig)

sig_eascd[["sig2"]] <- select_sig_500kb(sig_eascd[["sig2"]], kn_eascd_sig)

sig_easuc[["sig2"]] <- select_sig_500kb(sig_easuc[["sig2"]], kn_easuc_sig)

sig_easibd[["sig2"]] <- select_sig_500kb(sig_easibd[["sig2"]], kn_easibd_all_sig)

# select 4: snps that are not within 500kb  & select 5: common snp (maf>0.01)
select_novel<- function(sig)  {
  
  meta_sig_novel_500kb<-  filter(sig[["sig2"]], sig[["sig2"]]$within_500kb == FALSE)
  meta_sig_novel_500kb_common<-  filter(meta_sig_novel_500kb, meta_sig_novel_500kb$common == "common")
  meta_sig_novel_500kb_rare<-  filter(meta_sig_novel_500kb, meta_sig_novel_500kb$common == "rare")
  
  if  ("MAF.EAS" %in% colnames(meta_sig_novel_500kb)) {
    meta_sig_novel_500kb_common.EAS<-  filter(meta_sig_novel_500kb, meta_sig_novel_500kb$common.EAS == "common")
    meta_sig_novel_500kb_rare.EAS<-  filter(meta_sig_novel_500kb, meta_sig_novel_500kb$common.EAS == "rare")
    
    meta_sig_novel_500kb_common.EUREAS<-  filter(meta_sig_novel_500kb, meta_sig_novel_500kb$common.EUREAS == "common")
    meta_sig_novel_500kb_rare.EUREAS<-  filter(meta_sig_novel_500kb, meta_sig_novel_500kb$common.EUREAS == "rare")
    
  }
  
  meta_sig_novel_1mb<-  filter(sig[["sig2"]], sig[["sig2"]]$within_1mb == FALSE)
  meta_sig_novel_1mb_common<-  filter(meta_sig_novel_1mb, meta_sig_novel_1mb$common == "common")
  meta_sig_novel_1mb_rare<-  filter(meta_sig_novel_1mb, meta_sig_novel_1mb$common == "rare")
  
  if  ("MAF.EAS" %in% colnames(meta_sig_novel_1mb)) {
    meta_sig_novel_1mb_common.EAS<-  filter(meta_sig_novel_1mb, meta_sig_novel_1mb$common.EAS == "common")
    meta_sig_novel_1mb_rare.EAS<-  filter(meta_sig_novel_1mb, meta_sig_novel_1mb$common.EAS == "rare")
    
    meta_sig_novel_1mb_common.EUREAS<-  filter(meta_sig_novel_1mb, meta_sig_novel_1mb$common.EUREAS == "common")
    meta_sig_novel_1mb_rare.EUREAS<-  filter(meta_sig_novel_1mb, meta_sig_novel_1mb$common.EUREAS == "rare")
    
  }
  
  meta_sig_list <- append(sig, list(sig_novel_500kb= meta_sig_novel_500kb, sig_novel_500kb_common= meta_sig_novel_500kb_common,
                                    sig_novel_500kb_rare= meta_sig_novel_500kb_rare,
                                    sig_novel_1mb= meta_sig_novel_1mb, sig_novel_1mb_common= meta_sig_novel_1mb_common,
                                    sig_novel_1mb_rare= meta_sig_novel_1mb_rare))
  # append list if MAF.EAS exists
  if ("MAF.EAS" %in% colnames(meta_sig_novel_500kb)) {
    meta_sig_list <- append(meta_sig_list, list(sig_novel_500kb_common.EAS = meta_sig_novel_500kb_common.EAS,
                                                sig_novel_500kb_rare.EAS = meta_sig_novel_500kb_rare.EAS,
                                                sig_novel_500kb_common.EUREAS = meta_sig_novel_500kb_common.EUREAS,
                                                sig_novel_500kb_rare.EUREAS = meta_sig_novel_500kb_rare.EUREAS))
  }
  
  if ("MAF.EAS" %in% colnames(meta_sig_novel_1mb)) {
    meta_sig_list <- append(meta_sig_list, list(sig_novel_1mb_common.EAS = meta_sig_novel_1mb_common.EAS,
                                                sig_novel_1mb_rare.EAS = meta_sig_novel_1mb_rare.EAS,
                                                sig_novel_1mb_common.EUREAS = meta_sig_novel_1mb_common.EUREAS,
                                                sig_novel_1mb_rare.EUREAS = meta_sig_novel_1mb_rare.EUREAS))
  }
  
  return(meta_sig_list)
  
}

sig_cd <- select_novel(sig_cd)

sig_uc <- select_novel(sig_uc)

sig_ibd <- select_novel(sig_ibd)

sig_eascd <- select_novel(sig_eascd)

sig_easuc <- select_novel(sig_easuc)

sig_easibd <- select_novel(sig_easibd)


# write results to an excel
wr_tables <- function(lists, traits) {
  wb <- createWorkbook()
  
  for (name in names(lists)) {
    addWorksheet(wb, sheetName = name)
    writeData(wb, sheet = name, x = lists[[name]])
  }
  
  saveWorkbook(wb, paste0("./select/",traits, "_5e-08_SNP.xlsx"), overwrite = TRUE)
}
wr_tables(sig_cd, "EUR_CD")
wr_tables(sig_uc, "EUR_UC")
wr_tables(sig_ibd, "EUR_IBD")

wr_tables(sig_eascd, "EASEUR_CD")
wr_tables(sig_easuc, "EASEUR_UC")
wr_tables(sig_easibd, "EASEUR_IBD")

#### format for clumping ####
dir<- "/nobackup/sbcs/lyul1/IBD/metal_update/withoutUKB/gc_0.01/GC_Final/clumping/input/"
dir.create(dir, recursive = TRUE, showWarnings = FALSE)
format_dt<- function(sig, trait) {
  col_select<- sig %>% select(, c("CHR","BP","Allele1","REF","SNP","P.value" ))
  colnames(col_select)[3]<- "A1"
  colnames(col_select)[4]<- "A2"
  colnames(col_select)[6]<- "P"
  
  write.table(col_select, paste0(dir, trait), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
}

format_dt(sig_cd[["sig_novel_500kb"]], "EUR_CD_500k.txt")
format_dt(sig_uc[["sig_novel_500kb"]], "EUR_UC_500k.txt")
format_dt(sig_ibd[["sig_novel_500kb"]], "EUR_IBD_500k.txt")

format_dt(sig_eascd[["sig_novel_500kb"]], "EASEUR_CD_500k.txt")
format_dt(sig_easuc[["sig_novel_500kb"]], "EASEUR_UC_500k.txt")
format_dt(sig_easibd[["sig_novel_500kb"]], "EASEUR_IBD_500k.txt")

stop("Execution stopped")


# performed LD clumping in shell scr

####### Read processed meta results ########
read_ex_to_list<- function(excel_file)  {
  sheets <- excel_sheets(excel_file)
  
  list <- lapply(sheets, function(sheet) {
    read_excel(excel_file, sheet = sheet)
  })
  
  names(list) <- sheets
  
  return(list)
}

sig_cd <- read_ex_to_list("./select/EUR_CD_5e-08_SNP.xlsx")
sig_uc <- read_ex_to_list("./select/EUR_UC_5e-08_SNP.xlsx")
sig_ibd <- read_ex_to_list("./select/EUR_IBD_5e-08_SNP.xlsx")

sig_eascd<- read_ex_to_list("./select/EASEUR_CD_5e-08_SNP.xlsx")
sig_easuc<- read_ex_to_list("./select/EASEUR_UC_5e-08_SNP.xlsx")
sig_easibd<- read_ex_to_list("./select/EASEUR_IBD_5e-08_SNP.xlsx")


####### Clumped results #########
eur<- "/nobackup/sbcs/lyul1/IBD/metal_update/withoutUKB/gc_0.01/GC_Final/clumping/output_EUR.meta/"
eas.eur<- "/nobackup/sbcs/lyul1/IBD/metal_update/withoutUKB/gc_0.01/GC_Final/clumping/output_EAS.EURmeta/EURref/"
eas.eas<- "/nobackup/sbcs/lyul1/IBD/metal_update/withoutUKB/gc_0.01/GC_Final/clumping/output_EAS.EURmeta/EASref/"

read_data<- function(dir, trait) {
  data<- fread(paste0(dir, trait),fill = T)
  data<- filter(data, data$CHR!="CHR")
  data<- filter(data, data$CHR!="")
  
  data$CHR<- as.integer(data$CHR)
  data$BP<- as.integer(data$BP)
  data$P<- as.numeric(data$P)
  
  data<- data %>% arrange(CHR, BP)
  
  return(data)
}
cd.eur<- read_data(eur, "EUR_CD_500k_all_clumped.clumped")
uc.eur<- read_data(eur, "EUR_UC_500k_all_clumped.clumped")
ibd.eur<- read_data(eur, "EUR_IBD_500k_all_clumped.clumped")
# eas eur meta / eur ref data
cd_euref.eas<- read_data(eas.eur, "EASEUR_CD_500k_all_clumped.clumped")
uc_euref.eas<- read_data(eas.eur, "EASEUR_UC_500k_all_clumped.clumped")
ibd_euref.eas<- read_data(eas.eur, "EASEUR_IBD_500k_all_clumped.clumped")
# eas eur meta / eas ref data
cd_earef.eas<- read_data(eas.eas,"EASEUR_CD_500k_all_clumped.clumped")
uc_earef.eas<- read_data(eas.eas,"EASEUR_UC_500k_all_clumped.clumped")
ibd_earef.eas<- read_data(eas.eas,"EASEUR_IBD_500k_all_clumped.clumped")


# check lead variants check within 500kb
# merge back with effect size, p.val, etc.
window<- 500000
check_lead_var <- function(clumped, meta_res, window) {
  
  #clumped<- filter(clumped, clumped$common == "common") #EUR ref
  #clumped<- filter(clumped, clumped$common.EAS == "common") #EAS ref
  clumped<- clumped %>% arrange(CHR, BP)
  
  merged<- merge(clumped, meta_res, by = c("CHR", "BP", "SNP"), all.x = T) %>% arrange(CHR, BP)
  lead.var<- merged %>% select(., c("CHR",	"BP",	"SNP",	"ALT",	"REF",	"MarkerName",	"Allele1",	"Allele2",	
                                    "Effect",	"StdErr",	"P.value",	"Direction",	"HetISq",	"HetChiSq",	"HetDf",	"HetPVal",	"MAF.EUR"))
  
  filter_snps <- function(df, window) {
    df <- df %>% arrange(P.value)  # Sort by p-value (lowest first)
    keep <- rep(TRUE, nrow(df))    # Track which SNPs to keep
    
    for (i in seq_len(nrow(df))) {
      if (keep[i]) {  # If SNP is still marked as keep
        # Find all SNPs within the window
        within_window <- which(abs(df$BP - df$BP[i]) <= window)
        # Remove SNPs in the window except the one with the lowest p-value
        keep[within_window] <- FALSE  
        keep[i] <- TRUE  # Keep the current SNP
      }
    }
    
    return(df[keep, ])
  }
  
  # Apply filtering per chromosome
  lead.var.clean <- lead.var %>%
    group_by(CHR) %>%
    group_split() %>%
    lapply(filter_snps, window = window) %>%
    bind_rows()
  
  return(list(LeadVariants.clean = lead.var.clean, LeadVariants = lead.var, Raw = merged))
  
}

#
lead_var.cd <- check_lead_var(cd.eur, cd, window)
lead_var.uc <- check_lead_var(uc.eur, uc, window)
lead_var.ibd <- check_lead_var(ibd.eur, ibd, window)

lead_var.cd_euref.eas <- check_lead_var(cd_euref.eas, eascd, window)
lead_var.uc_euref.eas <- check_lead_var(uc_euref.eas, easuc, window)
lead_var.ibd_euref.eas <- check_lead_var(ibd_euref.eas, easibd, window)

lead_var.cd_earef.eas <- check_lead_var(cd_earef.eas, eascd, window)
lead_var.uc_earef.eas <- check_lead_var(uc_earef.eas, easuc, window)
lead_var.ibd_earef.eas <- check_lead_var(ibd_earef.eas, easibd, window)

# write results to an excel
wr_tables.lead <- function(lists, traits) {
  wb <- createWorkbook()
  
  for (name in names(lists)) {
    addWorksheet(wb, sheetName = name)
    writeData(wb, sheet = name, x = lists[[name]])
  }
  
  saveWorkbook(wb, paste0(traits, "_5e-08_SNP_Novel_clumped_500checked.xlsx"), overwrite = TRUE)
}
wr_tables.lead(lead_var.cd, "./select/EUR_CD")
wr_tables.lead(lead_var.uc, "./select/EUR_UC")
wr_tables.lead(lead_var.ibd, "./select/EUR_IBD")

wr_tables.lead(lead_var.cd_euref.eas, "./select/EASEUR_CD.euref")
wr_tables.lead(lead_var.uc_euref.eas, "./select/EASEUR_UC.euref")
wr_tables.lead(lead_var.ibd_euref.eas, "./select/EASEUR_IBD.euref")

wr_tables.lead(lead_var.cd_earef.eas, "./select/EASEUR_CD.earef")
wr_tables.lead(lead_var.uc_earef.eas, "./select/EASEUR_UC.earef")
wr_tables.lead(lead_var.ibd_earef.eas, "./select/EASEUR_IBD.earef")

### merge clumped results from eas panel and eur panel to select lead variants (loci) #####
merge_twopanel<- function(eas, eur, trait) {
  
  eas1 <- eas[["LeadVariants.clean"]] %>%
    mutate(Phenotype.loci = trait) %>%
    mutate(Source.loci = "EUREAS meta") %>%
    mutate(Ref.panel = "EAS Ref")
  
  eur1 <- eur[["LeadVariants.clean"]] %>%
    mutate(Phenotype.loci = trait) %>%
    mutate(Source.loci = "EUREAS meta") %>%
    mutate(Ref.panel = "EUR Ref")
  
  panel<- rbind(eas1, eur1) # %>% arrange(CHR, BP)
  
  panel$Reference <- ifelse(duplicated(panel$SNP) | duplicated(panel$SNP, fromLast = TRUE), "EUR/EAS Ref", panel$Ref.panel)
  
  panel_clean<- panel %>% distinct(MarkerName, .keep_all = TRUE)
  
  return(list(panel = panel, panel_clean = panel_clean))
  
}

panel.cd<- merge_twopanel(lead_var.cd_euref.eas, lead_var.cd_earef.eas, "CD")
panel.uc<- merge_twopanel(lead_var.uc_euref.eas, lead_var.uc_earef.eas, "UC")
panel.ibd<- merge_twopanel(lead_var.ibd_euref.eas, lead_var.ibd_earef.eas, "IBD")

# list <- list(eur.cd = ovl.cd, eur.uc =ovl.uc, eur.ibd =ovl.ibd, eur.cd.eas =ovl.cd1, eur.uc.eas= ovl.uc1, eur.ibd.eas= ovl.ibd1)
# write results to an excel
wr_tables.lead <- function(lists, traits) {
  wb <- createWorkbook()
  
  for (name in names(lists)) {
    addWorksheet(wb, sheetName = name)
    writeData(wb, sheet = name, x = lists[[name]])
  }
  
  saveWorkbook(wb, paste0("./select/EASEURref.panel.compare.", traits, ".xlsx"), overwrite = TRUE)
}

wr_tables.lead(panel.cd, "CD")
wr_tables.lead(panel.uc, "UC")
wr_tables.lead(panel.ibd, "IBD")

# check 500kb distance again for loci identified from eur and eas reference panel
window<- 500000
check_lead_var.panel <- function(panel, window) {
  
  filter_snps <- function(df, window) {
    df <- df %>% arrange(P.value)  # Sort by p-value (lowest first)
    keep <- rep(TRUE, nrow(df))    # Track which SNPs to keep
    
    for (i in seq_len(nrow(df))) {
      if (keep[i]) {  # If SNP is still marked as keep
        # Find all SNPs within the window
        within_window <- which(abs(df$BP - df$BP[i]) <= window)
        # Remove SNPs in the window except the one with the lowest p-value
        keep[within_window] <- FALSE  
        keep[i] <- TRUE  # Keep the current SNP
      }
    }
    
    return(df[keep, ])
  }
  
  lead.var.clean <- panel %>%
    group_by(CHR) %>%
    group_split() %>%
    lapply(filter_snps, window = window) %>%
    bind_rows()
  
}

panel.cd.check <- check_lead_var.panel(panel.cd[["panel_clean"]], window)
panel.uc.check <- check_lead_var.panel(panel.uc[["panel_clean"]], window)
panel.ibd.check <- check_lead_var.panel(panel.ibd[["panel_clean"]], window)

# write results to an excel
wr_tables.lead <- function(lists, traits) {
  wb <- createWorkbook()
  
  for (name in names(lists)) {
    addWorksheet(wb, sheetName = name)
    writeData(wb, sheet = name, x = lists[[name]])
  }
  
  saveWorkbook(wb, paste0("./select/EASEURref.panel.compare.", traits, "_500checked.xlsx"), overwrite = TRUE)
}

wr_tables.lead(panel.cd.check, "EASEUR_CD")
wr_tables.lead(panel.uc.check, "EASEUR_UC")
wr_tables.lead(panel.ibd.check, "EASEUR_IBD")


## check the overlap between EUR GWAS meta and EAS GWAS meta (combined EUR/EAS panel)
check_overlap<- function(EUR_loci, EAS_loci, trait) {
  
  add_EAS_loci <- anti_join(EAS_loci, EUR_loci, by = "MarkerName")
  EUR_loci$Phenotype.loci<- trait 
  EUR_loci$Source.loci<- "EUR meta"
  EUR_loci$Ref.panel<- "EUR Ref"
  EUR_loci$Reference<- "EUR"
  combined<- rbind(EUR_loci, add_EAS_loci)
  
  filter_snps <- function(df, window) {
    df <- df %>% arrange(P.value)  # sort by p-value (lowest first)
    keep <- rep(TRUE, nrow(df))     
    
    for (i in seq_len(nrow(df))) {
      if (keep[i]) {   
        # get all snps within the window
        within_window <- which(abs(df$BP - df$BP[i]) <= window)
        # remove selected snps (only those with Source.loci == "EUREAS meta")
        to_remove <- within_window[df$Source.loci[within_window] == "EUREAS meta"]
        # remove the selected SNPs, but keep the current one
        keep[to_remove] <- FALSE
        keep[i] <- TRUE  # Keep the current SNP
      }
    }
    
    return(df[keep, ])
  }
  
  combined.500k <- combined %>%
    group_by(CHR) %>%
    group_split() %>%
    lapply(filter_snps, window = window) %>%
    bind_rows()
  
  add_EAS_loci.500k <- filter(combined.500k, combined.500k$Source.loci == "EUREAS meta")
  
  return(list(all.EUREAS.500k = combined.500k, additional.500k = add_EAS_loci.500k, EUR.only = EUR_loci, 
              all.EUREAS = combined, additional = add_EAS_loci, EAS.only = EAS_loci))
  
}

overlap_cd<- check_overlap(lead_var.cd[["LeadVariants.clean"]], panel.cd.check, "CD")
overlap_uc<- check_overlap(lead_var.uc[["LeadVariants.clean"]], panel.uc.check, "UC")
overlap_ibd<- check_overlap(lead_var.ibd[["LeadVariants.clean"]], panel.ibd.check, "IBD")

# write results to an excel
wr_tables.lead <- function(lists, traits) {
  wb <- createWorkbook()
  
  for (name in names(lists)) {
    addWorksheet(wb, sheetName = name)
    writeData(wb, sheet = name, x = lists[[name]])
  }
  
  saveWorkbook(wb, paste0(traits, "_5e-08_SNP_All_loci_panel.check.addEAS.final.xlsx"), overwrite = TRUE)
}

wr_tables.lead(overlap_cd, "./select/EASEUR_CD")
wr_tables.lead(overlap_uc, "./select/EASEUR_UC")
wr_tables.lead(overlap_ibd, "./select/EASEUR_IBD")


# Novel EUR
combined_EUR<- rbind(overlap_cd[["EUR.only"]], overlap_uc[["EUR.only"]], overlap_ibd[["EUR.only"]])
# Additional EAS
add_EAS<- rbind(overlap_cd[["additional.500k"]], overlap_uc[["additional.500k"]], overlap_ibd[["additional.500k"]])
# Combined the EAS EUR
combined_EASEUR<- rbind(overlap_cd[["all.EUREAS.500k"]], overlap_uc[["all.EUREAS.500k"]], overlap_ibd[["all.EUREAS.500k"]])
T12<- list(novel_EUR = combined_EUR, additional = add_EAS, combined_all = combined_EASEUR)

wb <- createWorkbook()
for (name in names(T12)) {
  addWorksheet(wb, sheetName = name)
  writeData(wb, sheet = name, x = T12[[name]])
}

saveWorkbook(wb, "./select/Table1_2.xlsx", overwrite = TRUE)


##### cross check with gwas loci ######
read_ex_to_list<- function(excel_file)  {
  sheets <- excel_sheets(excel_file)
  
  list <- lapply(sheets, function(sheet) {
    read_excel(excel_file, sheet = sheet)
  })
  
  names(list) <- sheets
  
  return(list)
}

overlap_cd <- read_ex_to_list("./select/EASEUR_CD_5e-08_SNP_All_loci_panel.check.addEAS.final.xlsx")
overlap_uc <- read_ex_to_list("./select/EASEUR_UC_5e-08_SNP_All_loci_panel.check.addEAS.final.xlsx")
overlap_ibd <- read_ex_to_list("./select/EASEUR_IBD_5e-08_SNP_All_loci_panel.check.addEAS.final.xlsx")

excel_file<- "/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2023NG/41588_2023_1384_MOESM4_ESM (1).xlsx"
sheets <- excel_sheets(excel_file)
pre_loci <- lapply(sheets, function(sheet) {
  read_excel(excel_file, sheet = sheet)
})
names(pre_loci) <- sheets
easmeta<- pre_loci[["ST2"]][,c(1:5,9,11,12)]
colnames(easmeta)[1]<- "CHR"
colnames(easmeta)[2]<- "BP"

easeur<- pre_loci[["ST8"]]
colnames(easeur)<- easeur[1,]
easeur<- easeur[-1, c(2:6,12,14,16,17)]
colnames(easeur)[1]<- "CHR"
colnames(easeur)[2]<- "BP"

easeurcb<- rbind(easeur[, -8], easmeta)

check.cd<- merge(easeur, overlap_cd[["all.EUREAS"]], by = c("CHR", "BP"))
check.uc<- merge(easeur, overlap_uc[["all.EUREAS"]], by = c("CHR", "BP"))
check.ibd<- merge(easeur, overlap_ibd[["all.EUREAS"]], by = c("CHR", "BP"))

combined_check<- rbind(check.cd, check.uc, check.ibd)
write.csv(combined_check, "./select/check.with.latest.EASEURmeta.csv", row.names = F)
