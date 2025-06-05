### select top 10% genes acorss each cell type for stratified- LDSC #####

#select top genes
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(purrr)
#install.packages("BiocManager")
#BiocManager::install("biomaRt")
library(biomaRt)

setwd("/nobackup/sbcs/lyul1/IBD/single_cell/ibd")
ibd_matrix<- read.csv("/nobackup/sbcs/lyul1/IBD/single_cell/ibd_expression_matrix.csv")
ibd_matrix <- as.data.frame(ibd_matrix)
colnames(ibd_matrix) <- make.names(colnames(ibd_matrix))
# Extract gene names and sample data
long_expr <- pivot_longer(ibd_matrix,
                          cols = -Gene,
                          names_to = "sample_celltype_condition",
                          values_to = "expression")
split_info <- strsplit(long_expr$sample_celltype_condition, "_")
long_expr$sample_id <- sapply(split_info, `[`, 1)
long_expr$cell_type <- sapply(split_info, `[`, 2)
long_expr$condition <- sapply(split_info, `[`, 3)

# Replace "Control" as 0, "IBD" as 1
long_expr$condition <- ifelse(long_expr$condition == "Control", 0, 1)

long_expr <- long_expr %>%
  mutate(cell_group = case_when(
    cell_type %in% c("Early.Colonocyte", "Early.Enterocyte", 
                     "Intermediate.Colonocyte", "Intermediate.Enterocyte", "Mature.Colonocyte",
                     "Mature.Enterocyte") ~ "ABS",
    cell_type %in% c("Mature.Goblet", "Early.Goblet", 
                     "Goblet.Proliferating") ~ "GOB",
    TRUE ~ cell_type  # keep original cell type if not in ABS or GOB
  ))

write.csv(long_expr, "/nobackup/sbcs/lyul1/IBD/single_cell/ibd_expression_matrix_long.csv", row.names = F)


# Step 1 Compute t-statistics per gene per cell type
t_stats <- long_expr %>%
  group_by(Gene, cell_group) %>%  # Include cell_group here
  summarise(
    t_stat = tryCatch(
      t.test(expression ~ condition)$statistic,
      error = function(e) NA
    ),
    .groups = "drop"
  )

# Step 2 Rank genes by t-stat and select top 10% per cell type
top_genes <- t_stats %>%
  group_by(cell_group) %>%
  arrange(desc(t_stat), .by_group = TRUE) %>%
  mutate(rank = row_number(), total = n(), pct_rank = rank / total) %>%
  filter(pct_rank <= 0.10) %>%
  ungroup()
write.table(top_genes,
            file = file.path("/nobackup/sbcs/lyul1/IBD/single_cell/topgenes_update/", paste0("top_genes_all",".txt")),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

unique_groups <- unique(top_genes$cell_group)

for (group in unique_groups) {
  gene_list <- top_genes %>%
    dplyr::filter(cell_group == group) %>%
    dplyr::pull(Gene) %>%
    unique()
  
  group_clean <- gsub("[^A-Za-z0-9]", "_", group)
  

  write.table(gene_list,
              file = file.path("/nobackup/sbcs/lyul1/IBD/single_cell/topgenes_update/", paste0( group_clean, "_top_genes",".txt")),
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE)
}
# save gene names from top_genes
 write.table(unique(top_genes$gene),
             file = "/nobackup/sbcs/lyul1/IBD/single_cell/topgenes_update/top_genes_set.txt",
             row.names = FALSE,
             col.names = FALSE,
             quote = FALSE)



ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

#
cell_types <- c("ABS", "Stem", "GOB")
geneset_dir <- "/nobackup/sbcs/lyul1/IBD/single_cell/topgenes_update"
output_dir <- "/nobackup/sbcs/lyul1/IBD/single_cell/topgenes_update_ENSG"

dir.create(output_dir, showWarnings = FALSE)

# Loop through each cell type
for (cell in cell_types) {
  input_file <- file.path(geneset_dir, paste0(cell, "_top_genes.txt"))
  output_file <- file.path(output_dir, paste0(cell, "_top_genes_ENSG.txt"))
  
  gene_symbols <- read_lines(input_file)
  
  # Map to Ensembl IDs
  mapping <- getBM(
    attributes = c("hgnc_symbol", "ensembl_gene_id"),
    filters = "hgnc_symbol",
    values = gene_symbols,
    mart = ensembl
  )
  
  # save Ensembl IDs only
  write_lines(unique(mapping$ensembl_gene_id), output_file)
  
  cat(paste("Mapped", nrow(mapping), "genes for", cell, "\n"))
}
