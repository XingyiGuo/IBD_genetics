suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(tidyr))
suppressMessages(library("coloc"))
suppressMessages(library(purrr))
suppressMessages(library(readxl))
suppressMessages(library(arrow))
args <- commandArgs(trailingOnly = TRUE)

################################
####Step 1: Define variables####
################################

flank_length <- 1000000
xqtls <- c("UVA_Expr","UVA_APA", "UVA_AlterSplice","ACCC_Expr","ACCC_AlterSplice","UVA_ABS","UVA_STM")  
#ACCC_AlterSplice  ACCC_APA  ACCC_Expr ACCC_ABS ACCC_GOB ACCC_STM 
#UVA_AlterSplice  UVA_APA  UVA_Expr UVA_ABS UVA_GOB UVA_STM
xqtl_prefix <- xqtls[as.numeric(args[1])]
GWAS_ss_root="/Example_Data/"
xqtls_root=paste0("/Example_Data/",xqtl_prefix, "/rsformat/")

tag="knownGWASLoci"  #"metal_dec" for novel GWAS loci; "knownGWASLoci" for known GWAS loci
out_root="/IBD_GWAS_TWAS/eQTLs_coloc/colocres/knowGWASLoci/"

################################
####Step 2: Load functions####
################################
source("./coloc_gwas_xqtl_support.R")

################################
####Step 3: coloc ##############
################################
leadv_lmres_df<- data.frame(read_excel(paste0("/Example_Data/lmres_",tag,"/CRC_UVA_ACCC_ExprAlterSAPA_ABSGOBSTM_lmres.xlsx"), sheet = xqtl_prefix, col_names = TRUE))
#Select rows with BFsig
leadv_lmres_df_sig <- leadv_lmres_df[!is.na(leadv_lmres_df$BF) & (leadv_lmres_df$BF<0.05), ]
if(dim(leadv_lmres_df_sig)[1]==0){
cat("No sig ",xqtl_prefix," are found")
q()}

#Define a list to store resutls
stash <- matrix(integer(0), nrow = 0, ncol = 13) %>% as.data.frame()
names(stash) <- c(
   "disease", "lead_variant", "rs", "CHR","POS", "gene_id", "gene_name","nsnps",
   "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf") 

for (row_index in 1:dim(leadv_lmres_df_sig)[1]) {
    cat(row_index,"/",dim(leadv_lmres_df_sig)[1])
	##Load GWAS lead variants and the flank gwas variants 
	lead_variant = leadv_lmres_df_sig[row_index,"leadvariant"]
	rs=leadv_lmres_df_sig[row_index,"rs"]
	disease = leadv_lmres_df_sig[row_index,"Trait"]
	gene_id = leadv_lmres_df_sig[row_index,"gene_id"]
	gene_name = leadv_lmres_df_sig[row_index,"gene_name"]
	chr = gsub("chr", "", strsplit(lead_variant, ":")[[1]][1])
	pos = strsplit(lead_variant,":")[[1]][2]
	
	##Load GWAS ss, SNP(rs) as rownames
	gwas_variants_flank_leadvariant = get_gwas_variants_for_leadvariant2(GWAS_ss_root, xqtl_prefix, disease, chr, pos, flank_length)
	
	##Load xqtls
	xqtls=as.data.frame(read_parquet(paste0(xqtls_root, xqtl_prefix, ".cis_qtl_pairs.chr", chr,".parquet")))
	if(ncol(xqtls)==12){xqtls<-xqtls[, -12]}
	##Find xqtls overlapped with gwas variants, SNP(rs) as rownames
	#phenotype_id   variant_id  start_distance        af  ma_samples  ma_count  pval_nominal    slope  slope_se
	
	qtls.ss.geneid.raw <- filter(xqtls, str_detect(phenotype_id, gene_id))
	if(grepl("Expr", xqtl_prefix)){
		colnames(qtls.ss.geneid.raw) <- c("phenotype_id", "variant_id","start_distance", "af","ma_samples","ma_count","P","BETA","SE", "SNP")
	}else{
		colnames(qtls.ss.geneid.raw) <- c("phenotype_id", "variant_id","start_distance","end_distance","af","ma_samples","ma_count","P","BETA","SE", "SNP")
	}
	
	qtls.ss.geneid <- qtls.ss.geneid.raw[!duplicated(qtls.ss.geneid.raw$SNP) & !is.na(qtls.ss.geneid.raw$SNP), ]
	qtls.ss.geneid_comp <- qtls.ss.geneid[complete.cases(qtls.ss.geneid), ]
	qtls.ss.geneid <- filter(qtls.ss.geneid, P>0 & P<1)
	rownames(qtls.ss.geneid) <- qtls.ss.geneid$SNP
	qtls.ss.geneid <- qtls.ss.geneid %>% separate(variant_id, into = c("CHR", "BP", "A2", "A1", "B"), sep = "_")

	###Match xqtls and GWAS betas
	match_xqtl_gwas = matching_xqtls_and_gwas(gwas_variants_flank_leadvariant, qtls.ss.geneid)
		
	###Conduct coloc when over 50? variants overlapped between GWAS and pQTLs
	coloc_res=c("NA","NA","NA","NA","NA","NA")
	if(nrow(match_xqtl_gwas)>=50){
			ds1 = list(beta=match_xqtl_gwas$BETA_g, varbeta=(match_xqtl_gwas$SE_g)^2, 
					   pvalues=match_xqtl_gwas$P_g, snp=match_xqtl_gwas$SNP, 
					   type = "cc", N = as.numeric(sample_size_dict[[paste0(disease,"_ss")]]))
			ds2 = list(beta=match_xqtl_gwas$BETA_x, varbeta=(match_xqtl_gwas$SE_x)^2, 
					   pvalues=match_xqtl_gwas$P_x, snp=match_xqtl_gwas$SNP,
					   type = "quant", N = as.numeric(sample_size_dict[xqtl_prefix]))
			result <- coloc.abf(ds1, ds2, MAF=match_xqtl_gwas$maf)
			##Store coloc results if exits
			coloc_res=unlist(result$summary)
			if(coloc_res[6]>=0.5){fwrite(result$result, paste0(out_root,xqtl_prefix,"_",disease,"_",lead_variant,"_",gene_id,".coloc.csv"))}
	}else{
		cat("Skipping ",gene_id,", less than 50 overlap btw gwas region and xqtl.\n")
	}
	stash[nrow(stash) + 1, ] <- c(disease, lead_variant, rs, chr, pos, gene_id, gene_name, coloc_res)
}
fwrite(stash, paste0(out_root,xqtl_prefix,"_",as.character(flank_length/1000),"K.csv"))	
