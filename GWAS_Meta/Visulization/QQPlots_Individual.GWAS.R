library(qqman)
library(data.table)
library(dplyr)
library(tidyr)


setwd<- "/nobackup/sbcs/lyul1/IBD/metal_update/withUKB/gc_0.01/GC_Final/Plots/QQplot/ind"

# Finn
traits<- c("CD_STRICT2", "UC_STRICT2", "IBD_STRICT")
traits

plots<- function (trait) {
   inputfilename <- paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_FinnGen/R12/processed/finngen_R12_K11_", trait, "_markername_rmdup_chr1.22.txt") # finngen_R12_K11_CD_STRICT2_markername_rmdup_chr1.22.txt
    outfilename <- paste0("Fin_", trait, ".jpg")
    fig_title<- paste0("QQ Plot for GWAS of Fin_", trait)


    gwas.ss<- read.table(inputfilename, header=T, sep="\t", stringsAsFactors = F)

    if(is.na(outfilename)){
    print('Error: No output file name has been assigned!')
    }else {
    tiff(outfilename, width = 12, height = 8, units = "in", res = 300) #use .pdf as the extension in the filename
  
    qq(gwas.ss$pval , 
         main = fig_title,
        col = "blue4")
  
  
  dev.off()

}

}

plots(traits[1])
plots(traits[2])
plots(traits[3])


# MVP
traits<- c("Phe_555_1.EUR.CD", "Phe_555_2.EUR.UC", "Phe_555.EUR.IBD")

traits

plots<- function (trait) {
   inputfilename <- paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2023MVP/processed/", trait, "_markername_rmdup.txt") # Phe_555_2.EUR.UC_markername_rmdup.txt
    outfilename <- paste0("MVP_", trait, ".jpg")
    fig_title<- paste0("QQ Plot for GWAS of MVP_", trait)


    gwas.ss<- read.table(inputfilename, header=T, sep="\t", stringsAsFactors = F)

    if(is.na(outfilename)){
    print('Error: No output file name has been assigned!')
    }else {
    tiff(outfilename, width = 12, height = 8, units = "in", res = 300) #use .pdf as the extension in the filename
  
    qq(gwas.ss$pval , 
         main = fig_title,
        col = "blue4")
  
  
  dev.off()

}

}

plots(traits[1])
plots(traits[2])
plots(traits[3])

# IBDGC

traits<- c("cd", "uc", "ibd")

traits

plots<- function (trait) {
   inputfilename <- paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2017NG/processed/", trait, "_build38_markername_rmdup_maf_chr1.22.txt") # cd_build38_markername_rmdup_maf_chr1.22.txt
    outfilename <- paste0("NG17_", trait, ".jpg")
    fig_title<- paste0("QQ Plot for GWAS of NG17_", trait)


    gwas.ss<- read.table(inputfilename, header=T, sep="\t", stringsAsFactors = F)

    if(is.na(outfilename)){
    print('Error: No output file name has been assigned!')
    }else {
    tiff(outfilename, width = 12, height = 8, units = "in", res = 300) #use .pdf as the extension in the filename
  
    qq(gwas.ss$P.value , 
         main = fig_title,
        col = "blue4")
  
  
  dev.off()

}

}

plots(traits[1])
plots(traits[2])
plots(traits[3])


# Pan-UKB
traits<- c("555.1_CD", "555.2_UC", "555_IBD")

traits

plots<- function (trait) {
   inputfilename <- paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_Pan_UKB/processed/", trait, "_build38_markername_rmdup_chr1.22.txt") # 555_IBD_build38_markername_rmdup_chr1.22.txt 
    outfilename <- paste0("UKB_", trait, ".jpg")
    fig_title<- paste0("QQ Plot for GWAS of UKB_", trait)


    gwas.ss<- read.table(inputfilename, header=T, sep="\t", stringsAsFactors = F)

    if(is.na(outfilename)){
    print('Error: No output file name has been assigned!')
    }else {
    tiff(outfilename, width = 12, height = 8, units = "in", res = 300) #use .pdf as the extension in the filename
  
    qq(gwas.ss$P.value , 
         main = fig_title,
        col = "blue4")
  
  
  dev.off()

}

}

plots(traits[1])
plots(traits[2])
plots(traits[3])


# EAS NG23
traits<- c("CD", "UC", "IBD")

traits

plots<- function (trait) {
   inputfilename <- paste0("/nobackup/sbcs/lyul1/GWAS_SS/IBD/IBD_2023NG/processed/ibd_EAS_SiKJ_meta_", trait, "_markername_rmdup_chr1.22.txt") # ibd_EAS_SiKJ_meta_UC_markername_rmdup_chr1.22.txt
    outfilename <- paste0("NG23_", trait, ".jpg")
    fig_title<- paste0("QQ Plot for GWAS of NG23_", trait)


    gwas.ss<- read.table(inputfilename, header=T, sep="\t", stringsAsFactors = F)

    if(is.na(outfilename)){
    print('Error: No output file name has been assigned!')
    }else {
    tiff(outfilename, width = 12, height = 8, units = "in", res = 300) #use .pdf as the extension in the filename
  
    qq(gwas.ss$P.value , 
         main = fig_title,
        col = "blue4")
  
  
  dev.off()

}

}

plots(traits[1])
plots(traits[2])
plots(traits[3])


