library(qqman)
library(data.table)
library(dplyr)
library(tidyr)


setwd<- "/nobackup/sbcs/lyul1/IBD/metal_update/withUKB/gc_0.01/GC_Final/Plots/QQplot"

traits<- c("final_gc_EUR_CD", "final_gc_EUR_UC", "final_gc_EUR_IBD", "final_gc_EASEUR_CD", "final_gc_EASEUR_UC", "final_gc_EASEUR_IBD")

traits

plots<- function (trait) {
    inputfilename <- paste0("/nobackup/sbcs/lyul1/IBD/metal_update/withUKB/gc_0.01/GC_Final/", trait, "_META1.TBL")
    outfilename <- paste0("final_gc_", trait, ".jpg")
    fig_title<- paste0("QQ Plot for GWAS of ", trait)


    gwas.ss<- read.table(inputfilename, header=T, sep="\t", stringsAsFactors = F)

    if(is.na(outfilename)){
    print('Error: No output file name has been assigned!')
    }else {
    tiff(outfilename, width = 12, height = 8, units = "in", res = 300) #use .pdf as the extension in the filename
  
    qq(gwas.ss$P.value, 
         main = fig_title,
        col = "blue4")
  
  
  dev.off()

}

}

plots(traits[1])
plots(traits[2])
plots(traits[3])
plots(traits[4])
plots(traits[5])
plots(traits[6])

 
