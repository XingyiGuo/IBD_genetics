
library(data.table)
library(dplyr)
library(tidyr)


outdir <- "/data/l2_bioinfo1/sbcs/lyul1/IBD/metal_update/withUKB/gc_0.01/GC_Final/Plots/QQplot/ind"
indir <- "/data/l2_bioinfo1/lyul1/GWAS_SS/IBD/IBD_2023NG/processed"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd(outdir)

traits <- c("CD", "UC", "IBD")

plots <- function(trait) {
  inputfilename <- file.path(
    indir,
    paste0("ibd_EAS_SiKJ_meta_", trait, "_markername_rmdup_chr1.22.txt")
  )

  outfilename <- file.path(outdir, paste0("NG23_", trait, ".tiff"))
  fig_title <- paste0("QQ Plot for GWAS of NG23_", trait)

  gwas.ss <- read.table(
    inputfilename,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )

  pvals <- gwas.ss$P.value
  pvals <- pvals[!is.na(pvals) & pvals > 0 & pvals <= 1]

  observed <- -log10(sort(pvals))
  expected <- -log10(ppoints(length(pvals)))

  tiff(
    filename = outfilename,
    width = 12,
    height = 8,
    units = "in",
    res = 300,
    compression = "lzw",
    type = "cairo"
  )

  plot(
    expected,
    observed,
    pch = 20,
    cex = 0.7,
    col = "#0072B2",
    main = fig_title,
    xlab = expression(Expected~~-log[10](P)),
    ylab = expression(Observed~~-log[10](P))
  )

  abline(0, 1, lty = 2, col = "black", lwd = 1)

  dev.off()

  message("Saved: ", outfilename)
}

for (trait in traits) {
  plots(trait)
}
