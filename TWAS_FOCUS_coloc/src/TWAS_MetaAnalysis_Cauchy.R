#conda activate tftwas
library(ACAT)
library(data.table)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)

disList=c("cd","ibd","uc")
ModelTypeList=c("BULK_models","SC_models")
ModelType=ModelTypeList[as.numeric(args[1])]

if(ModelType=="BULK_models"){
	root="/IBD_MBD/IBD_GWAS_TWAS/TWAS/BULK_models/TFTWAS_res_metal_dec/"
	ModelSubTypeList=c("Expr","APA","AlterSplice")
}else if(ModelType=="SC_models"){
	root="/IBD_MBD/IBD_GWAS_TWAS/TWAS/SC_models/TFTWAS_res_metal_dec/"
	ModelSubTypeList=c("ABS","GOB","STM")
}

ModelSubTypeIndex=as.numeric(args[2]) #1[Expr], 2[APA], 3[AlterSplice]
ModelSubType=ModelSubTypeList[ModelSubTypeIndex]

bulk_EAS_models_list=c("CRC_TF_ACCC","CRC_TF_ACCC_APA","CRC_Asian_Splicing_hg38")
bulk_EUR_models_list=c("CRC_UVA","CRC_UVA_APA","CRC_UVA_Splice","CRC_TF_GTEx","CRC_TF_GTEx_APA","CRC_TF_GTEx_Splice")

sc_EAS_models_list=c("CRC_ASIAN_ABS","CRC_ASIAN_GOB","CRC_ASIAN_STM")
sc_EUR_models_list=c("CRC_EUR_ABS","CRC_EUR_GOB","CRC_EUR_STM","GTEX_colontrans_ABS","GTEX_colontrans_GOB","GTEX_colontrans_STM")

GWAS_SS_EAS=c("EAS_CD", "EAS_IBD", "EAS_UC")
GWAS_SS_EUR=c("EUR_CD", "EUR_IBD","EUR_UC")

for(i in 1:3){ #cd, ibd, uc
	if(ModelType=="BULK_models"){
		x1_filename=paste0(root, bulk_EUR_models_list[0+ModelSubTypeIndex], ".", GWAS_SS_EUR[i], ".TWAS") #uva
		x2_filename=paste0(root, bulk_EUR_models_list[3+ModelSubTypeIndex], ".", GWAS_SS_EUR[i], ".TWAS") #gtex
		x3_filename=paste0(root, bulk_EAS_models_list[0+ModelSubTypeIndex], ".", GWAS_SS_EAS[i], ".TWAS") #asian
	}
	else if(ModelType=="SC_models"){
		x1_filename=paste0(root, sc_EUR_models_list[0+ModelSubTypeIndex], ".", GWAS_SS_EUR[i], ".TWAS") #uva
		x2_filename=paste0(root, sc_EUR_models_list[3+ModelSubTypeIndex], ".", GWAS_SS_EUR[i], ".TWAS") #gtex
		x3_filename=paste0(root, sc_EAS_models_list[0+ModelSubTypeIndex], ".", GWAS_SS_EAS[i], ".TWAS") #asian
	}
	x1_df=as.data.frame(fread(x1_filename))
	x2_df=as.data.frame(fread(x2_filename))
	x3_df=as.data.frame(fread(x3_filename))

	#keep genes with pred_perf_r2>0.01
	x1_df = filter(x1_df, pred_perf_r2>0.01)
	x2_df = filter(x2_df, pred_perf_r2>0.01)
	x3_df = filter(x3_df, pred_perf_r2>0.01)

	#add suffix
	colnames(x1_df)=paste0(colnames(x1_df), "_UVA")
	colnames(x2_df)=paste0(colnames(x2_df), "_GTEx") 
	colnames(x3_df)=paste0(colnames(x3_df), "_ASIAN")

	#rename gene names
	if((ModelType=="BULK_models") && (ModelSubTypeIndex==3)){ #CRC_Asian_Splicing_hg38
	
		##Rename gene chr5:40716498:40728345:clu_39972_-:ENSG00000113638.13 to chr5:40716498:40728345:ENSG00000113638
		x1_df$gene_hg38_UVA <- gsub(":(clu_.*:)?(ENSG[0-9]+)\\..*$", ":\\2", x1_df$gene_UVA)
		x2_df$gene_hg38_GTEx <- gsub(":(clu_.*:)?(ENSG[0-9]+)\\..*$", ":\\2", x2_df$gene_GTEx)
		
		names(x1_df)[names(x1_df) == "gene_hg38_UVA"] <- "gene_key"
		names(x2_df)[names(x2_df) == "gene_hg38_GTEx"] <- "gene_key"
		names(x3_df)[names(x3_df) == "gene_hg38_ASIAN"] <- "gene_key"
		##merge two EUR dataframe p values toether
		x1_x2_df = merge(x1_df, x2_df, by="gene_key", all = TRUE)
		#Cannot have NAs in the p-values!, 
		x1_x2_df <- x1_x2_df %>%
		rowwise() %>%
		mutate(META_UVA_GTEx_P = ifelse(is.na(pvalue_UVA), pvalue_GTEx,
						  ifelse(is.na(pvalue_GTEx), pvalue_UVA,
								 ACAT(c(pvalue_UVA, pvalue_GTEx)))))		
		# merge EURs with Asian
		x1_x2_x3_df = merge(x1_x2_df, x3_df, by="gene_key", all = TRUE)
		x1_x2_x3_df <- x1_x2_x3_df %>%
		rowwise() %>%
		mutate(META_UVA_GTEx_ASIAN_P = ifelse(is.na(META_UVA_GTEx_P), pvalue_ASIAN,
									   ifelse(is.na(pvalue_ASIAN), META_UVA_GTEx_P,
									   ACAT(c(pvalue_ASIAN, META_UVA_GTEx_P)))))
	}else if((ModelType=="BULK_models") && (ModelSubTypeIndex==2)){
		##Rename gene ENST00000106952.3 to be ENST00000106952
		x1_df$gene_hg38_UVA <- gsub("\\..*", "", x1_df$gene_UVA)
		x2_df$gene_hg38_GTEx <- gsub("\\..*", "", x2_df$gene_GTEx)
		x3_df$gene_hg38_ASIAN <- gsub("\\..*", "", x3_df$gene_ASIAN)
		
		print(head(x1_df))
		print(head(x2_df))
		print(head(x3_df))
		
		names(x1_df)[names(x1_df) == "gene_hg38_UVA"] <- "gene_key"
		names(x2_df)[names(x2_df) == "gene_hg38_GTEx"] <- "gene_key"
		names(x3_df)[names(x3_df) == "gene_hg38_ASIAN"] <- "gene_key"
		##merge two EUR dataframe p values toether
		x1_x2_df = merge(x1_df, x2_df, by="gene_key", all = TRUE)
		#Cannot have NAs in the p-values!, 
		x1_x2_df <- x1_x2_df %>%
		rowwise() %>%
		mutate(META_UVA_GTEx_P = ifelse(is.na(pvalue_UVA), pvalue_GTEx,
						  ifelse(is.na(pvalue_GTEx), pvalue_UVA,
								 ACAT(c(pvalue_UVA, pvalue_GTEx)))))		
		# merge EURs with Asian
		x1_x2_x3_df = merge(x1_x2_df, x3_df, by="gene_key", all = TRUE)
		x1_x2_x3_df <- x1_x2_x3_df %>%
		rowwise() %>%
		mutate(META_UVA_GTEx_ASIAN_P = ifelse(is.na(META_UVA_GTEx_P), pvalue_ASIAN,
									   ifelse(is.na(pvalue_ASIAN), META_UVA_GTEx_P,
									   ACAT(c(pvalue_ASIAN, META_UVA_GTEx_P)))))
	}else{ #merge by gene symbol
		names(x1_df)[names(x1_df) == "gene_name_UVA"] <- "gene_name_key"
		names(x2_df)[names(x2_df) == "gene_name_GTEx"] <- "gene_name_key"
		names(x3_df)[names(x3_df) == "gene_name_ASIAN"] <- "gene_name_key"
		##merge two EUR dataframe p values toether
		x1_x2_df = merge(x1_df, x2_df, by="gene_name_key", all = TRUE)
		x1_x2_df <- x1_x2_df %>%
		rowwise() %>%
		mutate(META_UVA_GTEx_P = ifelse(is.na(pvalue_UVA), pvalue_GTEx,
						  ifelse(is.na(pvalue_GTEx), pvalue_UVA,
								 ACAT(c(pvalue_UVA, pvalue_GTEx)))))		
		# merge EURs with Asian
		x1_x2_x3_df = merge(x1_x2_df, x3_df, by="gene_name_key", all = TRUE)
		x1_x2_x3_df <- x1_x2_x3_df %>%
		rowwise() %>%
		mutate(META_UVA_GTEx_ASIAN_P = ifelse(is.na(META_UVA_GTEx_P), pvalue_ASIAN,
									   ifelse(is.na(pvalue_ASIAN), META_UVA_GTEx_P,
									   ACAT(c(pvalue_ASIAN, META_UVA_GTEx_P)))))
	}

	fwrite(x1_x2_x3_df, paste0("/IBD_MBD/IBD_GWAS_TWAS/TWAS/Cauchy_meta/CRC_UVA_GTEX_ASIAN_ONEGWAS_",disList[i],"_",ModelType,"_",ModelSubType,".csv"))
}

