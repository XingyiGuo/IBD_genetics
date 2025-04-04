suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(tidyr))
suppressMessages(library(purrr))
suppressMessages(library(readxl))
suppressMessages(library(stringr))


####################################
####Part 1: define some dictionaries###
####################################

sample_size_dict <- list(
  Breast_ss=247173,Colon_ss=125478,Lung_ss=85716,Ovary_ss=63347,Pancreas_ss=14998,Prostate_ss=140306,
  deCODE_pqtl=35559,UKBPPP_pqtl=34557,
  wb_mqtl=670, OneK1K_sceqtl=100, 
  Breast_sqtl=396, Breast_eqtl=396, Colon_sqtl=368, Colon_eqtl=368, Lung_sqtl=515, Lung_eqtl=515,
  Ovary_sqtl=167, Ovary_eqtl=167, Pancreas_sqtl=305, Pancreas_eqtl=305, Prostate_sqtl=221, Prostate_eqtl=221,
  CD_EUR_ss=867409, UC_EUR_ss=869870, IBD_EUR_ss=874262, 
  CD_EAS_ss=22828, UC_EAS_ss=22318, IBD_EAS_ss=29849, 
  UVA_Expr=423,UVA_APA=423,UVA_AlterSplice=423,UVA_ABS=423, UVA_STM=423,
  ACCC_Expr=364,ACCC_APA=364,ACCC_AlterSplice=364
)

GWAS_ss_dict <- list(
"Breast"="breast_dbSNPs_impute_summary_statistic_polyfun.txt",
"Colon"="colorectal_dbSNPs_impute_summary_statistic_polyfun.txt",
"Lung"="lung_dbSNPs_impute_summary_statistic_polyfun.txt",
"Ovary"="ovarian_dbSNPs_impute_summary_statistic_polyfun.txt",
"Pancreas"="pancreatic_dbSNPs_impute_summary_statistic_polyfun.txt",
"Prostate"="prostate_dbSNPs_impute_summary_statistic_polyfun.txt"
)

disase_tissue_dict <- list(
"Breast"="Breast_Mammary_Tissue",
"Colon"="Colon_Transverse",
"Lung"="Lung",
"Ovary"="Ovary",
"Pancreas"="Pancreas",
"Prostate"="Prostate"
)

####################################
####Part 2: define some functions###
####################################

get_cisgene_annotation <- function(gene_annot_file, chrom, pos, flank)
{
	pos_flank_start=as.numeric(pos)-as.numeric(flank)
	pos_flank_end=as.numeric(pos)+as.numeric(flank)
    gene_df <- read.table(gene_annot_file,header=F,stringsAsFactors =F,sep="\t",fill=T)
	gene_df$V4<-as.numeric(gene_df$V4)
	gene_df$V5<-as.numeric(gene_df$V5)
	gene_df$V3 <- as.character(gene_df$V3)
	gene_df$V7 <- as.character(gene_df$V7)
    gene_df_tmp <- filter(gene_df,(V3=="gene")&(V1%in%(paste0('chr',chrom))))
	gene_df1_1 <- filter(gene_df_tmp, 
						V7 == "+" & 
						((V4 >= pos_flank_start & V4 <= pos_flank_end) | 
						(V5 >= pos_flank_start & V5 <= pos_flank_end)))
	gene_df1_2 <- filter(gene_df_tmp, 
						V7 == "-" & 
						((V5 >= pos_flank_start & V5 <= pos_flank_end) | 
						(V4 >= pos_flank_start & V4 <= pos_flank_end)))
	gene_df1 <- rbind(gene_df1_1, gene_df1_2)
	if(nrow(gene_df1)>0){    
		geneid <- str_extract(gene_df1[,9], "ENSG\\d+.\\d+")
		genename <- gsub("gene_name (\\S+);","\\1",str_extract(gene_df1[,9], "gene_name (\\S+);"), perl=T)
		gene_used <- as.data.frame(cbind(geneid,genename,gene_df1[,c(1,4,5,3,7)]))
		colnames(gene_used) <- c("geneid","genename","chr","start","end","anno","direction")
		return(gene_used)
	}
}

get_gwas_variants_for_leadvariant <- function(GWAS_ss_root, disease, chr, pos, flank_length){
	disease_gwas=data.frame(fread(paste0(GWAS_ss_root,GWAS_ss_dict[disease])))
	disease_gwas$SNP <- as.character(disease_gwas$SNP)
	disease_gwas$BETA <- as.numeric(disease_gwas$BETA)
	disease_gwas$SE <- as.numeric(disease_gwas$SE)
	disease_gwas$CHR <- as.numeric(disease_gwas$CHR)
	colnames(disease_gwas)[3]<-"POS"
	disease_gwas$POS <- as.numeric(disease_gwas$POS)
	disease_gwas$A1FREQ <- as.numeric(disease_gwas$A1FREQ)
	disease_gwas$maf <- ifelse(disease_gwas$A1FREQ < 0.5, disease_gwas$A1FREQ, 1 - disease_gwas$A1FREQ)
	disease_gwas<-filter(disease_gwas, maf>0 & maf<1 & PVALUE>0 & PVALUE<1)
	flank_start=as.numeric(pos)-flank_length
	flank_end=as.numeric(pos)+flank_length
	gwas_variants_in_flank <- filter(disease_gwas, CHR == as.numeric(chr) & POS >= flank_start & POS <= flank_end)
	gwas_variants_in_flank <- gwas_variants_in_flank[!duplicated(gwas_variants_in_flank$SNP), ]
	rownames(gwas_variants_in_flank)=gwas_variants_in_flank$SNP
	return(as.data.frame(gwas_variants_in_flank))
}

get_gwas_variants_for_leadvariant2 <- function(GWAS_ss_root, xqtl_prefix, disease, chr, pos, flank_length){
	if(grepl("UVA", xqtl_prefix)){
		disease_gwas_raw=data.frame(fread(paste0(GWAS_ss_root,"EUR_",disease,"_perchr/chr", chr, ".gwas.ss.gz")))
	}else if(grepl("ACCC", xqtl_prefix)){
		disease_gwas_raw=data.frame(fread(paste0(GWAS_ss_root,"EAS_",disease,"_perchr/chr", chr, ".gwas.ss.gz")))
	}
	disease_gwas<-disease_gwas_raw[disease_gwas_raw$SNP!=".",]
	disease_gwas$SNP <- as.character(disease_gwas$SNP)
	disease_gwas$BETA <- as.numeric(disease_gwas$BETA)
	disease_gwas$SE <- as.numeric(disease_gwas$SE)
	disease_gwas$CHR <- as.numeric(disease_gwas$CHR)
	disease_gwas$BP <- as.numeric(disease_gwas$BP)
	disease_gwas$FRQ <- as.numeric(disease_gwas$FRQ)
	disease_gwas$maf <- ifelse(disease_gwas$FRQ < 0.5, disease_gwas$FRQ, 1 - disease_gwas$FRQ)
	disease_gwas<-filter(disease_gwas, maf>0 & maf<1 & P>0 & P<1)
	flank_start=as.numeric(pos)-flank_length
	flank_end=as.numeric(pos)+flank_length
	gwas_variants_in_flank <- filter(disease_gwas, CHR == as.numeric(chr) & BP >= flank_start & BP <= flank_end)
	gwas_variants_in_flank <- gwas_variants_in_flank[!duplicated(gwas_variants_in_flank$SNP), ]
	rownames(gwas_variants_in_flank)=gwas_variants_in_flank$SNP
	return(as.data.frame(gwas_variants_in_flank))
}

get_proteins_in_same_chrom_with_leadvariant <- function(pqtls_root, pqtl_db, chr_){
	if(pqtl_db=="UKBPPP"){
		protein_loopup_table=as.data.frame(fread(paste0(pqtls_root,"olink_protein_map_3k_v1.tsv")))
		poteins_in_same_chrom_df=filter(protein_loopup_table, chr==chr_)
		poteins_in_same_chrom=as.vector(poteins_in_same_chrom_df$OlinkID)
	}else if(pqtl_db=="deCODE"){
		protein_loopup_table=as.data.frame(fread(paste0(pqtls_root,"deCODE-SMP-lookup.csv")))
		poteins_in_same_chrom_df=filter(protein_loopup_table, chr==chr_)
		poteins_in_same_chrom=as.vector(poteins_in_same_chrom_df$file)
	}
	return(as.data.frame(poteins_in_same_chrom))
}

matching_xqtls_and_gwas<-function(gwas_var, xqtls_var){
	merged_df <- merge(gwas_var, xqtls_var, by = "SNP", suffixes = c("_g", "_x"))
	if(nrow(merged_df)>0){
		for (i in 1:nrow(merged_df)) {
		#Make sure xQTL A1 match to GWAS A1. If not, switch A1 and A2 in xQTL and correct beta
		if (toupper(merged_df$A2_x[i]) == toupper(merged_df$A1_g[i])) {
		  merged_df$BETA_x[i] <- -merged_df$BETA_x[i]
		  tmp <- toupper(merged_df$A1_x[i])
		  merged_df$A1_x[i] <- toupper(merged_df$A2_x[i])
		  merged_df$A2_x[i] <- tmp
		}
	}
	}
	return(merged_df)
}