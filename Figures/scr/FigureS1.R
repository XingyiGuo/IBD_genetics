library(data.table)
library(dplyr)
library(ggplot2)

indir <- paste0(
  "/data/l2_bioinfo1/lyul1/IBD/metal_update/withUKB/",
  "gc_0.01/GC_Final/processed"
)

filtered_indir <- file.path(
  indir,
  "rm_loci"
)

outdir <- paste0(
  "/data/l2_bioinfo1/lyul1/IBD/metal_update/withUKB/",
  "gc_0.01/GC_Final/Plots/QQplot/rm_loci"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

traits <- c(
  "EUR_CD",
  "EUR_UC",
  "EUR_IBD",
  "EASEUR_CD",
  "EASEUR_UC",
  "EASEUR_IBD"
)
qq_data <- function(pvals, label, trait) {
  pvals <- suppressWarnings(
    as.numeric(pvals)
  )

  pvals <- pvals[
    !is.na(pvals) &
      is.finite(pvals) &
      pvals > 0 &
      pvals <= 1
  ]

  n <- length(pvals)

  if (n == 0L) {
    return(
      data.frame(
        Expected = numeric(0),
        Observed = numeric(0),
        Groups = character(0),
        Trait = character(0),
        stringsAsFactors = FALSE
      )
    )
  }

  expected <- -log10(ppoints(n))
  observed <- -log10(sort(pvals))

  data.frame(
    Expected = expected,
    Observed = observed,
    Groups = label,
    Trait = trait,
    stringsAsFactors = FALSE
  )
}

for (trait in traits) {
  message("")
  message("Processing ", trait)

  full_gwas_file <- file.path(
    indir,
    paste0(
      trait,
      "_META1.TBL.SNP.txt"
    )
  )

  filtered_gwas_file <- file.path(
    filtered_indir,
    paste0(
      trait,
      "_META1.TBL.SNP_rm_pre_loci.300000.txt"
    )
  )

  if (!file.exists(full_gwas_file)) {
    warning(
      "Full GWAS file was not found: ",
      full_gwas_file
    )
    next
  }

  if (!file.exists(filtered_gwas_file)) {
    warning(
      "Filtered GWAS file was not found: ",
      filtered_gwas_file
    )
    next
  }


  gwas_full <- fread(
    full_gwas_file,
    fill = TRUE,
    showProgress = FALSE
  )

  gwas_filtered <- fread(
    filtered_gwas_file,
    fill = TRUE,
    showProgress = FALSE
  )

  # Standardize the P-value column name
  if (!"P.value" %in% names(gwas_full) &&
      "P" %in% names(gwas_full)) {
    setnames(
      gwas_full,
      old = "P",
      new = "P.value"
    )
  }

  if (!"P.value" %in% names(gwas_filtered) &&
      "P" %in% names(gwas_filtered)) {
    setnames(
      gwas_filtered,
      old = "P",
      new = "P.value"
    )
  }

  if (!"P.value" %in% names(gwas_full)) {
    warning(
      "No P.value or P column was found in: ",
      full_gwas_file
    )
    next
  }

  if (!"P.value" %in% names(gwas_filtered)) {
    warning(
      "No P.value or P column was found in: ",
      filtered_gwas_file
    )
    next
  }

  # Generate QQ data
  qq_full <- qq_data(
    pvals = gwas_full$P.value,
    label = "Full GWAS meta-analysis",
    trait = trait
  )

  qq_filtered <- qq_data(
    pvals = gwas_filtered$P.value,
    label = "After removing known loci",
    trait = trait
  )

  qq_combined <- bind_rows(
    qq_full,
    qq_filtered
  )

  if (nrow(qq_combined) == 0L) {
    warning(
      "No valid P values were found for ",
      trait
    )
    next
  }

  # Set the order of groups in the legend
  qq_combined$Groups <- factor(
    qq_combined$Groups,
    levels = c(
      "Full GWAS meta-analysis",
      "After removing known loci"
    )
  )


group_colors <- c(
  "Full GWAS meta-analysis" = "#E69F00",
  "After removing known loci" = "#0072B2"
)

  qq_plot <- ggplot(
    qq_combined,
    aes(
      x = Expected,
      y = Observed,
      color = Groups
    )
  ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "black",
      linewidth = 0.5
    ) +
    geom_point(
      alpha = 0.7,
      size = 0.8
    ) +
    scale_color_manual(
      values = group_colors,
      breaks = names(group_colors),
      drop = FALSE
    ) +
    labs(
      title = paste0(
        "QQ Plot for GWAS of ",
        trait
      ),
      subtitle = paste0(
        "Before and after removing variants within 300 kb ",
        "of known IBD loci"
      ),
      x = "Expected -log10(P)",
      y = "Observed -log10(P)",
      color = NULL
    ) +
    theme_minimal(
      base_size = 12
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 14,
        face = "bold"
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 10
      ),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "right",
      legend.text = element_text(size = 10),
      plot.background = element_rect(
        fill = "white",
        color = NA
      ),
      panel.background = element_rect(
        fill = "white",
        color = NA
      )
    )


  if (capabilities("cairo")) {
    output_file <- file.path(
      outdir,
      paste0(
        "QQ_Plot_",
        trait,
        ".300000.png"
      )
    )

    ggsave(
      filename = output_file,
      plot = qq_plot,
      device = grDevices::png,
      type = "cairo-png",
      width = 7,
      height = 6,
      units = "in",
      dpi = 300,
      bg = "white"
    )
  } else {
    warning(
      "Cairo is unavailable. Saving ",
      trait,
      " as PDF instead."
    )

    output_file <- file.path(
      outdir,
      paste0(
        "QQ_Plot_",
        trait,
        ".300000.pdf"
      )
    )

    ggsave(
      filename = output_file,
      plot = qq_plot,
      device = grDevices::pdf,
      width = 7,
      height = 6,
      units = "in",
      bg = "white"
    )
  }
