library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)



setwd<- "/nobackup/sbcs/lyul1/IBD/metal_update/withUKB/gc_0.01/GC_Final/Plots/Manhattan"

# lead novel variants
alltraits<- read.csv("/nobackup/sbcs/lyul1/IBD/metal_update/withUKB/gc_0.01/select/lead.var.for.plots.csv", header=T)
alltraits <- filter(alltraits, !(alltraits$SNP == "rs11168249" & alltraits$Trait == "EUR_UC"))
alltraits <- filter(alltraits, !(alltraits$SNP == "rs12523195" & alltraits$Trait == "EUR_UC"))
alltraits <- filter(alltraits, !(alltraits$SNP == "rs7253302" & alltraits$Trait == "EUR_IBD"))
alltraits <- filter(alltraits, !(alltraits$SNP == "rs11168249" & alltraits$Trait == "EUR_IBD"))
#alltraits <- alltraits[alltraits$NCHROBS.EAS == "" | is.na(alltraits$NCHROBS.EAS), ]

inputfile<- "_META1.TBL.SNP.txt"

mk.mhtplt<- function (trait) {

inputfilename<- paste0(trait, inputfile)
gwas.ss<- fread(paste0("/nobackup/sbcs/lyul1/IBD/metal_update/withUKB/gc_0.01/GC_Final/processed/", inputfilename), fill = T)   #%>% subset(CHR %in% 1:22)
gwas.ss$BP<- as.numeric(gwas.ss$BP)
colnames(gwas.ss)[11]<- "P"
gwas.ss$P<- as.numeric(gwas.ss$P)

lead.var<- filter(alltraits, alltraits$Trait == trait)
snpsOfInterest<- lead.var$SNP

chrom_colors<- c("grey", "skyblue")

chrom_colors<- c("#789DBC", "#FFE3E3")

chrom_colors<- c("#16423C" , "#6A9C89")

col_scheme <- rep(chrom_colors, length.out = 22)

outfilename<- paste0(trait, "Manhattan1.tiff")

if(is.na(outfilename)){
  print('Error: No output file name has been assigned!')
} else {
  
  # Prepare the dataset
  don <- gwas.ss %>%
    
    # Compute chromosome size
    group_by(CHR) %>%
    summarise(chr_len=max(BP)) %>%
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(gwas.ss, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate(BPcum = BP + tot) %>%
    
    # Add highlight and annotation information
    mutate(is_highlight = ifelse(SNP %in% snpsOfInterest, "yes", "no"))

  # Prepare X axis
  axisdf <- don %>%
    group_by(CHR) %>%
    summarize(center = (max(BPcum) + min(BPcum)) / 2)

  # Make the plot
  p <- ggplot(don, aes(x = BPcum, y = -log10(P))) +
    
    # Show all points
    geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.2) +
    scale_color_manual(values = col_scheme) +
    
    # Custom X axis:
    scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) +
    scale_y_continuous(expand = c(0, 0), breaks = c(seq(0, 150, by = 10))) + 
    # Expand y-axis limit to 150
    # scale_y_continuous(expand = c(0, 0), limits = c(0, 150)) +

    # scale_y_continuous(expand = c(0, 0)) +  # Remove space between plot area and x-axis
    # Remove vertical grid lines (major and minor) and keep only the horizontal line for y-intercept
        theme(
    panel.grid.major.y = element_blank(),   # Remove major horizontal grid lines
    panel.grid.minor.y = element_blank(),   # Remove minor horizontal grid lines
    panel.grid.major.x = element_blank(),   # Remove major vertical grid lines
    panel.grid.minor.x = element_blank()    # Remove minor vertical grid lines
        ) +
        
    # Add horizontal line at 5e-08
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red") + # -log10(5e-8)


    # Add highlighted points
    geom_point(data = subset(don, is_highlight == "yes"), color = "orange", size = 2) +
    
    # Add label using ggrepel to avoid overlapping
    #geom_label_repel(data = subset(don, is_highlight == "yes"), aes(label = SNP), size = 3, family = "Arial") +
    
    # Customize the theme:
    theme_bw() +
    theme(
      text = element_text(family = "Arial"),  
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    
    # Change x-axis label
    xlab("Chromosome")  +
    
    # Add figure title
    ggtitle(paste0("Manhattan Plot for GWAS Meta-analysis of ", trait))
    
  
  # Save the plot as a TIFF file
  ggsave(filename = outfilename, plot = p, width = 12, height = 8, units = "in", dpi = 500, device = "tiff")
}


}


mk.mhtplt("EUR_CD")
mk.mhtplt("EUR_UC")
mk.mhtplt("EUR_IBD")
mk.mhtplt("EASEUR_CD")
mk.mhtplt("EASEUR_UC")
mk.mhtplt("EASEUR_IBD")
