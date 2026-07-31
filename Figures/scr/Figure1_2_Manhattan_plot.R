library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)



setwd("/data/l2_bioinfo1/lyul1/IBD/metal_update/withUKB/gc_0.01/GC_Final/Plots/Manhattan")

alltraits <- read.csv(
  "/data/l2_bioinfo1/lyul1/IBD/metal_update/withUKB/gc_0.01/select/lead.var.for.plots_update.csv",
  header = TRUE
)

inputfile <- "_META1.TBL.SNP.txt"

mk.mhtplt <- function(trait, y_max = 30) {
  
  inputfilename <- paste0(trait, inputfile)
  
  gwas.ss <- fread(
    paste0(
      "/data/l2_bioinfo1/lyul1/IBD/metal_update/withUKB/gc_0.01/GC_Final/processed/",
      inputfilename
    ),
    fill = TRUE
  )
  
  gwas.ss$BP <- as.numeric(gwas.ss$BP)
  colnames(gwas.ss)[11] <- "P"
  gwas.ss$P <- as.numeric(gwas.ss$P)
  
  gwas.ss <- gwas.ss %>%
    filter(CHR %in% 1:22, !is.na(BP), !is.na(P), P > 0)
  
  lead.var <- filter(alltraits, Trait == trait)
  snpsOfInterest <- lead.var$SNP
  
  # Color-blind-safe alternating chromosome colors.
  chrom_colors <- c("#0072B2", "#999999")
  col_scheme <- rep(chrom_colors, length.out = 22)
  
  outfilename <- paste0(trait, "_Manhattan_truncated.tiff")
  
  don <- gwas.ss %>%
    group_by(CHR) %>%
    summarise(chr_len = max(BP), .groups = "drop") %>%
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    left_join(gwas.ss, ., by = "CHR") %>%
    arrange(CHR, BP) %>%
    mutate(
      BPcum = BP + tot,
      logp = -log10(P),
      logp_plot = pmin(logp, y_max),
      is_truncated = logp > y_max,
      is_highlight = SNP %in% snpsOfInterest
    )
  
  axisdf <- don %>%
    group_by(CHR) %>%
    summarise(center = (max(BPcum) + min(BPcum)) / 2, .groups = "drop")
  
  p <- ggplot(don, aes(x = BPcum, y = logp_plot)) +
    
    geom_point(
      data = subset(don, !is_truncated),
      aes(color = as.factor(CHR)),
      alpha = 0.75,
      size = 1.1
    ) +
    
    geom_point(
      data = subset(don, is_truncated),
      aes(color = as.factor(CHR)),
      shape = 17,
      alpha = 0.9,
      size = 1.6
    ) +
    
    scale_color_manual(values = col_scheme) +
    
    geom_hline(
      yintercept = -log10(5e-8),
      linetype = "dashed",
      color = "black", #"#D55E00"
      linewidth = 0.4
    ) +
    
    geom_point(
      data = subset(don, is_highlight),
      color = "#E69F00",
      size = 2.0
    ) +
    
    scale_x_continuous(
      label = axisdf$CHR,
      breaks = axisdf$center,
      expand = c(0.01, 0.01)
    ) +
    
    scale_y_continuous(
      limits = c(0, y_max),
      breaks = seq(0, y_max, by = 5),
      expand = c(0, 0)
    ) +
    
    labs(
      x = "Chromosome",
      y = expression(-log[10](P)),
      title = paste0("Manhattan Plot for GWAS Meta-analysis of ", trait),
      caption = paste0(
        "Y-axis truncated at -log10(P) = ", y_max,
        "; variants exceeding this value are shown as triangles at the upper boundary."
      )
    ) +
    
    theme_bw(base_family = "Arial") +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.minor.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.caption = element_text(hjust = 0, size = 9),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  ggsave(
    filename = outfilename,
    plot = p,
    width = 12,
    height = 8,
    units = "in",
    dpi = 500,
    device = "tiff",
    compression = "lzw",
      type = "cairo"
  )
}


mk.mhtplt("EUR_CD")
mk.mhtplt("EUR_UC")
mk.mhtplt("EUR_IBD")
mk.mhtplt("EASEUR_CD")
mk.mhtplt("EASEUR_UC")
#mk.mhtplt("EASEUR_IBD")
