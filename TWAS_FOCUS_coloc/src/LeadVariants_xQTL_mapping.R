# conda activate tftwas_coloc
library("methods")
library("glmnet")
library("stringr")
library("data.table")
library("dplyr")
library("readxl")

setwd("/IBD_MBD/IBD_GWAS_TWAS/eQTLs_coloc/")
"%&%" <- function(a,b) paste(a,b, sep = "")

###### I. paramter ###### 
args <- commandArgs(trailingOnly=T)
eths=c("UVA","ACCC")
profile_types=c("Expr","AlterSplice", "APA","ABS","GOB","STM")

cur_eth=eths[as.numeric(args[1])] #args[1]
cur_profile=profile_types[as.numeric(args[2])] #args[2]
chrom <- as.numeric(args[3]) #args[3]

prefix <- "CRC_" %&% cur_eth %&% "_" %&% cur_profile

###### II. input files ######
#Lead variants
leadV_folder="/IBD_MBD/IBD_GWAS_TWAS/eQTLs_coloc/"
if(cur_eth=="UVA"){
	leadV_df <- as.data.frame(read_excel(paste0(leadV_folder, "gwas_loci.xlsx"), sheet = "EUR",skip = 1))
}else if(cur_eth=="ACCC"){
	leadV_df <- as.data.frame(read_excel(paste0(leadV_folder, "gwas_loci.xlsx"), sheet = "EUREAS",skip = 1))
}

# genotype file
if (cur_eth == "UVA"){
  genotype_prefix <- paste0("/UVA/genotype/chr", chrom,".gt.gt")
  # vcf head
  vcfhead <- read.table("/UVA/genotype/vcf.head",header = F, stringsAsFactors = F)
  vcfhead <- gsub("(VM\\d+)_\\S+","\\1",vcfhead,perl=T)
  # gene annotation
  gene_annot_file <- "/ref/gencode/gencode.v26.annotation.gtf"
  if(cur_profile=="Expr"){
    profile_prefix = "/UVA/expression/463samples_423white_QN.expression.PEER_inverse.txt" 
  }else if(cur_profile=="AlterSplice"){
    profile_prefix = "/UVA/expression/UVA.leafcutter.bed.quantile_norm.txt_PEER_inverse.txt_01032023"
  }else if(cur_profile=="APA"){
    profile_prefix <-"/UVA/expression/UVA_CRC_APA.txt_QN_05042023" # quantile normalization
    # peer  factors for APA
    peer_file <- "/UVA/expression/UVA_CRC_APA.txt_QN_05042023_PEER_factors.txt"
  }else if(cur_profile=="ABS"){
	profile_prefix = "/IBD_MBD/IBD_GWAS_TWAS/eQTLs_coloc/UVA_ABS.txt" 
  }else if(cur_profile=="GOB"){
	profile_prefix = "/IBD_MBD/IBD_GWAS_TWAS/eQTLs_coloc/UVA_GOB.txt" 
  }else if(cur_profile=="STM"){
	profile_prefix = "/IBD_MBD/IBD_GWAS_TWAS/eQTLs_coloc/UVA_STM.txt" 
  }
}else if (cur_eth == "ACCC"){
  genotype_prefix <- paste0("/genotype4/merge/chr", chrom,".gt.gt")
  # vcf head
  vcfhead <- read.table("/genotype4/merge/vcf.head",header = F, stringsAsFactors = F)
  vcfhead <- gsub("WG\\d+-\\S+-(CNUHHCRC_\\d+)@\\d+","\\1",vcfhead, perl=TRUE)
  vcfhead <- gsub('CNUHH-CRC-(\\d+)_QC-3158-Cai_\\S+','CNUHHCRC_\\1',vcfhead, perl=T)
  vcfhead <- gsub('FMMU-(CC\\d+N)_QC-3158-Cai_\\S+','FMMU_\\1',vcfhead, perl =T)
  vcfhead <- gsub('(HCES2-\\d+)_QC-3158-Cai_\\S+','\\1', vcfhead, perl=T)
  # gene annotation
  gene_annot_file <- "/ref/gencode/gencode.v19.annotation.gtf"
  if(cur_profile=="Expr"){
    profile_prefix = "/accc/expression/381samples_used_364sample_QN.expression_PEER_inverse.txt"
  }else if(cur_profile=="AlterSplice"){
    profile_prefix = "/accc/expression/ACCC_364samples.leafcutter_quantile_norm_PEER_inverse.txt"
  }else if(cur_profile=="APA"){
    profile_prefix <-"/accc/expression/Asian_CRC_APA.txt_QN_05102023" # quantile normalization
    # peer  factors for APA
    peer_file <- "/accc/expression/Asian_CRC_APA.txt_QN_05102023_PEER_factors.txt"
  }else if(cur_profile=="ABS"){
	profile_prefix = "/IBD_MBD/IBD_GWAS_TWAS/eQTLs_coloc/ACCC_ABS.txt" 
  }else if(cur_profile=="GOB"){
	profile_prefix = "/IBD_MBD/IBD_GWAS_TWAS/eQTLs_coloc/ACCC_GOB.txt" 
  }else if(cur_profile=="STM"){
	profile_prefix = "/IBD_MBD/IBD_GWAS_TWAS/eQTLs_coloc/ACCC_STM.txt" 
  }
}

###### III. set function  ######
get_gene_annotation <- function(gene_annot_file, chrom, bp, flank){
    gene_df <- read.table(gene_annot_file,header=F,stringsAsFactors =F,sep="\t",fill=T)
    gene_df1 <- filter(gene_df,V3 %in% "gene")
    geneid <- str_extract(gene_df1[,9], "ENSG\\d+.\\d+")
    genename <- gsub("gene_name (\\S+);","\\1",str_extract(gene_df1[,9], "gene_name (\\S+);"), perl=T)
    gene_used <- as.data.frame(cbind(geneid,genename,gene_df1[,c(1,4,5,3)]))
    colnames(gene_used) <- c("geneid","genename","chr","start","end","anno")
    gtf_used <- filter(gene_used, gene_used[,3] %in% ('chr' %&% chrom))
	variant_flank_start=max(as.numeric(bp)-flank, 0 )
	variant_flank_end=as.numeric(bp)+flank
	gtf_used2 <- filter(gtf_used, 
					   ((gtf_used[,4] > variant_flank_start) & (gtf_used[,4] < variant_flank_end)) | 
					   ((gtf_used[,5] > variant_flank_start) & (gtf_used[,5] < variant_flank_end)))
    gtf_used2
}

get_transcript_annotation <- function(gene_annot_file, chrom, bp, flank){
    gene_df <- read.table(gene_annot_file,header=F,stringsAsFactors =F,sep="\t",fill=T)
    gene_df1 <- filter(gene_df,V3 %in% "transcript")
    geneid <- str_extract(gene_df1[,9], "ENST\\d+.\\d+")
    genename <- gsub("gene_name (\\S+);","\\1",str_extract(gene_df1[,9], "gene_name (\\S+);"), perl=T)
    gene_used <- as.data.frame(cbind(geneid, genename,gene_df1[,c(1,4,5,3)]))
    colnames(gene_used) <- c("geneid","genename","chr","start","end","anno")
    gtf_used <- filter(gene_used, gene_used[,3] %in% ('chr' %&% chrom))
	variant_flank_start=max(as.numeric(bp)-flank, 0 )
	variant_flank_end=as.numeric(bp)+flank
	gtf_used2 <- filter(gtf_used, 
					   ((gtf_used[,4] > variant_flank_start) & (gtf_used[,4] < variant_flank_end)) | 
					   ((gtf_used[,5] > variant_flank_start) & (gtf_used[,5] < variant_flank_end)))
    gtf_used2
}

get_gene_expression <- function(profile_prefix, gene_annot) {
	data_df <-  as.data.frame(fread(profile_prefix, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
	# Set the first row as column names
	colnames(data_df) <- as.character(unlist(data_df[1,])) #ENSG
	# Remove the first row from the data
	data_df <- data_df[-1,]
	# Set the first column as row names
	rownames(data_df) <- as.character(unlist(data_df[,1])) #sampleID
	# Remove the first row from the data
	data_df <- data_df[,-1]
	data_df <- data_df[order(row.names(data_df)), ]
	# Return gene dataframe with selected geneid
	data_df_selected_geneid <- data_df %>% select(one_of(intersect(gene_annot$geneid, colnames(data_df))))
	data_df_selected_geneid
}

get_splice_expression <- function(profile_prefix, gene_annot) {
	data_df <-  as.data.frame(fread(profile_prefix, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
	# Set the first row as column names
	colnames(data_df) <- as.character(unlist(data_df[1,])) #ENSG
	# Remove the first row from the data
	data_df <- data_df[-1,]
	# Set the first column as row names
	rownames(data_df) <- as.character(unlist(data_df[,1])) #sampleID
	# Remove the first row from the data
	data_df <- data_df[,-1]
	data_df <- data_df[order(row.names(data_df)), ]
	##gene_annot contains genes in the +/- 1M regions of lead variants
	mat = matrix(ncol = 0, nrow = nrow(data_df)) 
	selected_splices=data.frame(mat) 
	for(row_index in 1:nrow(gene_annot)){
		#Locate the spliceID which contains target ENSG 
		selected_elements <- grep(gene_annot[row_index, "geneid"], colnames(data_df), value = TRUE)
		#Select columns with target spliceID
		selected_splice_for_one_geneid=data.frame(data_df[, colnames(data_df) %in% selected_elements])
		cat(gene_annot[row_index, "geneid"], ":", dim(selected_splice_for_one_geneid)[2], ", ")
		if(dim(selected_splice_for_one_geneid)[2]==1){
			colnames(selected_splice_for_one_geneid)=selected_elements
		}
		selected_splices=cbind(selected_splices, selected_splice_for_one_geneid)
	} 
	rownames(selected_splices)=rownames(data_df)
	selected_splices
}

get_apa_expression <- function(profile_prefix, trans_annot) {
	data_df_trans <-  as.data.frame(fread(profile_prefix, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
	data_df=t(data_df_trans)
	# Set the first row as column names
	colnames(data_df) <- as.character(unlist(data_df[1,])) #ENSG
	# Remove the first row from the data
	data_df <- data_df[-1,]
	# Set the first column as row names
	rownames(data_df) <- as.character(unlist(data_df[,1])) #sampleID
	# Remove the first row from the data
	data_df <- data_df[,-1]
	data_df <- data_df[order(row.names(data_df)), ]
	row_names=rownames(data_df)
	# Convert the content back to numeric
	data_df <- apply(data_df, 2, as.numeric)
	# Convert back to a data frame if needed (optional)
	data_df <- as.data.frame(data_df)
	rownames(data_df)=row_names
	##trans_annot contains genes in the +/- 1M regions of lead variants
	mat = matrix(ncol = 0, nrow = nrow(data_df)) 
	selected_apa=data.frame(mat) 
	for(row_index in 1:nrow(trans_annot)){
		#Locate the APA which contains target ENST 
		selected_elements <- grep(trans_annot[row_index, "geneid"], colnames(data_df), value = TRUE)
		#Select columns with target APA
		selected_apa_for_one_geneid=data.frame(data_df[, colnames(data_df) %in% selected_elements])
		cat(trans_annot[row_index, "geneid"], ":", dim(selected_apa_for_one_geneid)[2], ", ")
		if(dim(selected_apa_for_one_geneid)[2]==1){
			colnames(selected_apa_for_one_geneid)=selected_elements
		}
		selected_apa=cbind(selected_apa, selected_apa_for_one_geneid)
	}
	rownames(selected_apa)=rownames(data_df)
	selected_apa
}

get_maf_filtered_genotype <- function(genotype_file, vcfhead, maf, samples) {
  gt_df<- fread(genotype_file,sep="\t",header=F)
  colnames(gt_df) <- vcfhead
  gt_df <- as.data.frame(gt_df)
  row.names(gt_df) <- gt_df$ID
  gt_df <- gt_df[,-c(1:5)]
  gt_df1 <- as.data.frame(t(gt_df))
  gt_df1 <- gt_df1[row.names(gt_df1) %in% samples,]
  effect_allele_freqs <- colMeans(gt_df1) / 2
  gt_df1 <- gt_df1[,which((effect_allele_freqs >= maf) & (effect_allele_freqs <= 1-maf))]
  gt_df1 <- gt_df1[order(row.names(gt_df1)), ]
  gt_df1
}

format_output <- function(model_coefficients){
	# Extracting values for the Intercept
	int_beta <- model_coefficients["(Intercept)", "Estimate"]
	int_se <- model_coefficients["(Intercept)", "Std. Error"]
	int_t <- model_coefficients["(Intercept)", "t value"]
	int_p <- model_coefficients["(Intercept)", "Pr(>|t|)"]
	# Extracting values for the slope (x)
	slope_beta <- model_coefficients["one_chrom_gt_df_variant", "Estimate"]
	slope_se <- model_coefficients["one_chrom_gt_df_variant", "Std. Error"]
	slope_t <- model_coefficients["one_chrom_gt_df_variant", "t value"]
	slope_p <- model_coefficients["one_chrom_gt_df_variant", "Pr(>|t|)"]
	#return res
	return(c(slope_beta, slope_se, slope_t, slope_p, int_beta, int_se, int_t, int_p))
}

###### IV. run analysis ######
# parameter for analysis 
maf=0.01
flank=1000000
# output file
model_summary_file <- './lmres_metal_dec/' %&% prefix %&% '_chr' %&% chrom %&% '_model_summaries.txt'
model_summary_cols <- c('Trait', 'leadvariant', 'rs', 'gene_id', 'gene_name', 'slope_beta', 'slope_se', 'slope_t', 'slope_p', 'int_beta', 'int_se', 'int_t', 'int_p')
write(model_summary_cols, file = model_summary_file, ncol = 13, sep = '\t')

if(cur_profile=="APA"){
	data_df_trans <-  as.data.frame(fread(profile_prefix, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
	data_df=t(data_df_trans)
	# Set the first row as column names
	colnames(data_df) <- as.character(unlist(data_df[1,])) #ENSG
	# Remove the first row from the data
	data_df <- data_df[-1,]
	# Set the first column as row names
	exp_samples <- as.character(unlist(data_df[,1])) #sampleID
	# Peer factors for APA
    peer <- read.table(peer_file,header=T,stringsAsFactors=F,row.names=1,sep="\t")
    peer <- peer[order(row.names(peer)), ]
    peer <- as.data.frame(peer)
}else{
	data_df <-  as.data.frame(fread(profile_prefix, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
	# Set the first row as column names
	colnames(data_df) <- as.character(unlist(data_df[1,])) #ENSG
	# Remove the first row from the data
	data_df <- data_df[-1,]
	# Set the first column as row names
	exp_samples <- as.character(unlist(data_df[,1])) #sampleID
}
one_chrom_gt_df <- get_maf_filtered_genotype(genotype_prefix, vcfhead, maf, exp_samples) #genotype order by sample names
gt_exp_samples <- intersect(exp_samples,row.names(one_chrom_gt_df))

# inite analysis
one_chrom_leadV_df=filter(leadV_df, leadV_df[,"CHR"]==chrom)
for(row_index in 1:nrow(one_chrom_leadV_df)){
	if(dim(one_chrom_leadV_df)[1]>0){
		bp=one_chrom_leadV_df[row_index, "BP"]
		alt=one_chrom_leadV_df[row_index, "ALT"]
		ref=one_chrom_leadV_df[row_index, "REF"]
		SNP=one_chrom_leadV_df[row_index, "SNP"]
		Trait=one_chrom_leadV_df[row_index, "Phenotype"]
		cat(chrom, " ", bp, " ", SNP)
		# Load annot file, chr, pos, ENGX, gene_name
		if(cur_profile=="APA"){
			one_chrom_gene_annot <- get_transcript_annotation(gene_annot_file, chrom, bp, flank) #Find genes located in the +/- flank regions of the lead variant
		}else{
			one_chrom_gene_annot <- get_gene_annotation(gene_annot_file, chrom, bp, flank)
		}
		cat("\nNumber of annot (ENSG/ENST) is ", dim(one_chrom_gene_annot)[1], "\n")
		# Load profile values
		if(cur_profile=="APA"){
			one_chrom_profile_df <- get_apa_expression(profile_prefix, one_chrom_gene_annot) #apa profile order by sample names
		}else if(cur_profile=="AlterSplice"){
			one_chrom_profile_df <- get_splice_expression(profile_prefix, one_chrom_gene_annot) #splice order by sample names
		}else{ #Expr, ABS, GOB, STM
			one_chrom_profile_df <- get_gene_expression(profile_prefix, one_chrom_gene_annot) #gene expression order by sample names
		}
		# Count n_tests as 
		tests <- colnames(one_chrom_profile_df)
		n_tests <- length(one_chrom_profile_df)
		cat("\nNumber of test is", n_tests)
		if(dim(one_chrom_profile_df)[2]==1){
			one_chrom_profile_df <- one_chrom_profile_df[row.names(one_chrom_profile_df) %in% gt_exp_samples,colnames(one_chrom_profile_df), drop = FALSE]
			one_chrom_profile_df <- one_chrom_profile_df[match(gt_exp_samples, row.names(one_chrom_profile_df)), colnames(one_chrom_profile_df), drop = FALSE]
		}else{
			one_chrom_profile_df <- one_chrom_profile_df[row.names(one_chrom_profile_df) %in% gt_exp_samples,] #cause pb when only 1 col present
			one_chrom_profile_df <- one_chrom_profile_df[match(gt_exp_samples, row.names(one_chrom_profile_df)), ]
		}
		# Check if the variants present in the one_chrom_gt_df file
		variants_list=colnames(one_chrom_gt_df)
		if(cur_eth=="UVA"){
			index1 <- which(grepl(paste0("chr",chrom,":",bp,":",ref,":",alt), variants_list))
			index2 <- which(grepl(paste0("chr",chrom,":",bp,":",alt,":",ref), variants_list))
		}else if(cur_eth=="ACCC"){
			index1 <- which(grepl(paste0(chrom,":",bp,";",chrom,":",bp,":",ref,":",alt), variants_list))
			index2 <- which(grepl(paste0(chrom,":",bp,";",chrom,":",bp,":",alt,":",ref), variants_list))
		}
		if(length(index1)==1){
			index=index1
		}else if(length(index2)==1){
			index=index2
		}else{
			index=0
		}
		lead_variant <- variants_list[index]
		if(length(lead_variant)!=0){ #lead variants exists in gt
			for(gene_index in 1:n_tests){
				one_chrom_gt_df_variant=one_chrom_gt_df[1:nrow(one_chrom_gt_df), lead_variant]
				cat(gene_index, "/", n_tests, "\n")
				testID <- tests[gene_index]
				#testID
				#ENSG00000227232.5 for Expr, 
				#chr1:14829:14970:clu_50317:ENSG00000227232.5 for AlterSplice
				#ENST00000377921.7:RSU1_chr10:16590611-16593496 for APA,
				if(cur_profile=="APA"){
					ENSX = str_extract(testID, "ENST\\d+.\\d+")
				}else{
					ENSX = str_extract(testID, "ENSG\\d+.\\d+")
				}
				gene_name <- as.character(one_chrom_gene_annot$genename[one_chrom_gene_annot$geneid == ENSX])
				model_summary <- c(Trait, lead_variant, SNP, testID, gene_name, 0,0,0,0,0,0,0,0)
				profile_df_one_gene=one_chrom_profile_df[1:nrow(one_chrom_profile_df), testID]
				#Run lm analyses
				if(length(profile_df_one_gene)>0){ 
					if(cur_profile=="APA"){
						peer_used <- peer
						na_id <- which(is.na(profile_df_one_gene))
						if (length(na_id) > 0){
							profile_df_one_gene <- profile_df_one_gene[-na_id]
							peer_used <- peer_used[-na_id,]
							one_chrom_gt_df_variant <- one_chrom_gt_df_variant[-na_id]
						}
						res <- residuals(lm(profile_df_one_gene ~ ., peer_used))
						adj_expression <- as.numeric(qnorm( rank(res,ties.method="r") / (length(res)+1))) # inverse normalization
						model <- lm(adj_expression ~ one_chrom_gt_df_variant)
					}else{
						model <- lm(profile_df_one_gene ~ one_chrom_gt_df_variant)
					}
				} #Else, model_summary is default
				model_summary <- c(Trait, lead_variant, SNP, testID, gene_name, format_output(summary(model)$coefficients))
				write(model_summary, file = model_summary_file, append = TRUE, ncol = 13, sep = '\t')
			}
		}else{
			cat("\nNo variants from gt match chr",chrom,":",bp,":",ref,":",alt)
		}
	}else{
		cat("\nNo lead variants from ", chrom)
	}
}
