library(readxl)
library(openxlsx)
library(enrichR)

##------------------------------------------
## Enrichr setup
##------------------------------------------
listEnrichrSites()
setEnrichrSite("Enrichr")

#dbs <- listEnrichrDbs()
dbs <- c(
    "Reactome_Pathways_2024",
    "KEGG_2021_Human",
    "MSigDB_Hallmark_2020")

##------------------------------------------
## Input / Output
##------------------------------------------
input_file <- "/SupplementaryTables.xlsx"
output_file <- "/enriched_genes.xlsx"
sheet_name <- "ST22"

##------------------------------------------
## Read gene list
##------------------------------------------
gene_df <- read_excel(
    input_file,
    sheet = sheet_name,
    skip = 1
)

gene_list <- unique(trimws(na.omit(gene_df$`Gene Name`)))
cat("Number of unique genes:", length(gene_list), "\n")

##------------------------------------------
## Function
##------------------------------------------
run_enrichr <- function(gene_list,
                        databases,
                        padj_cutoff = 0.05) {
    enrichment_results <- enrichr(gene_list, databases)
    if (is.null(enrichment_results)) {
        stop("Enrichr query failed.")
    }
    filtered_results <- list()
    for (db in names(enrichment_results)) {
        result_df <- enrichment_results[[db]]
        if (!"Adjusted.P.value" %in% colnames(result_df))
            next
        result_df <- result_df[
            result_df$Adjusted.P.value < padj_cutoff,
        ]
        if (nrow(result_df) == 0)
            next
        result_df$Source <- db
        filtered_results[[db]] <- result_df
    }
    if (length(filtered_results) == 0)
        return(NULL)
    combined_df <- do.call(rbind, filtered_results)

    desired_columns <- c(
        "Term",
        "Overlap",
        "P.value",
        "Adjusted.P.value",
        "Old.P.value",
        "Old.Adjusted.P.value",
        "Odds.Ratio",
        "Combined.Score",
        "Genes",
        "Source"
    )
    desired_columns <- intersect(desired_columns, colnames(combined_df))
    combined_df <- combined_df[, desired_columns]
    combined_df <- combined_df[
        order(
            combined_df$Source,
            combined_df$Adjusted.P.value
        ),
    ]
    rownames(combined_df) <- NULL
    return(combined_df)
}

##------------------------------------------
## Run enrichment
##------------------------------------------
enrichment_df <- run_enrichr(
    gene_list = gene_list,
    databases = dbs
)

##------------------------------------------
## Save results
##------------------------------------------
if (is.null(enrichment_df)) {
    message("No significant enrichment found.")
} else {
    write.xlsx(
        enrichment_df,
        file = output_file,
        overwrite = TRUE,
        rowNames = FALSE
    )
    cat(
        "Saved",
        nrow(enrichment_df),
        "significant enrichment terms to:\n",
        output_file,
        "\n"
    )
}