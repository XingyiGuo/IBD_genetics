#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (!length(hit)) return(default)
  if (hit == length(args)) stop("Missing value for ", flag)
  args[[hit + 1]]
}

need_arg <- function(flag) {
  value <- get_arg(flag)
  if (is.null(value) || value == "") stop("Required argument missing: ", flag)
  value
}

clean_text <- function(x) {
  trimws(gsub("[\r\n]", "", as.character(x)))
}

clean_chr <- function(x) {
  sub("^chr", "", clean_text(x))
}

clean_id <- function(x) {
  x <- clean_text(x)
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  fifelse(is.na(x) | x == "", "locus", x)
}

strip_ensg_version <- function(x) {
  sub("\\..*$", "", clean_text(x))
}

extract_ensg <- function(x) {
  x <- clean_text(x)
  out <- rep(NA_character_, length(x))

  for (i in seq_along(x)) {
    m <- regexpr("ENSG[0-9]+(\\.[0-9]+)?", x[i], perl = TRUE)
    if (m[1] > 0) {
      out[i] <- regmatches(x[i], m)
    }
  }

  strip_ensg_version(out)
}

escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x, perl = TRUE)
}

has_gene_symbol <- function(x, gene_name) {
  x <- clean_text(x)
  gene_name <- clean_text(gene_name)

  if (length(gene_name) != 1 || is.na(gene_name) || gene_name == "") {
    return(rep(FALSE, length(x)))
  }

  pattern <- paste0("(^|[^A-Za-z0-9])", escape_regex(gene_name), "([^A-Za-z0-9]|$)")
  grepl(pattern, x, perl = TRUE)
}

normalize_qtl_type <- function(x) {
  x <- clean_text(x)
  fifelse(x %in% c("AS", "AlterSplice"), "AlterSplice", x)
}

infer_qtl_type_from_panel <- function(panel_label, qtl_type = NA_character_) {
  qtl_type <- normalize_qtl_type(qtl_type)
  panel_label <- clean_text(panel_label)

  out <- qtl_type
  missing <- is.na(out) | out == ""

  out[missing & grepl("AlterSplice|_AS$|_AS_", panel_label, ignore.case = TRUE)] <- "AlterSplice"
  out[missing & grepl("Expr|eQTL", panel_label, ignore.case = TRUE)] <- "Expr"
  out[missing & grepl("APA", panel_label, ignore.case = TRUE)] <- "APA"
  out[missing & grepl("ABS", panel_label, ignore.case = TRUE)] <- "ABS"
  out[missing & grepl("GOB", panel_label, ignore.case = TRUE)] <- "GOB"
  out[missing & grepl("STM", panel_label, ignore.case = TRUE)] <- "STM"

  out
}

pick_col <- function(dt, candidates, label) {
  hit <- candidates[candidates %in% names(dt)]
  if (!length(hit)) {
    stop("Missing ", label, " column. Found: ", paste(names(dt), collapse = ", "))
  }
  hit[1]
}

gwas_summary_for <- function(gwas_dir, ancestry, disease) {
  ancestry <- clean_text(ancestry)
  disease <- clean_text(disease)

  if (ancestry == "EAS") {
    return(file.path(gwas_dir, paste0("EASEUR_", disease, ".alt.harmo.eas.ma")))
  }
  if (ancestry == "EUR") {
    return(file.path(gwas_dir, paste0("EUR_", disease, ".alt.harmo.ma")))
  }

  stop("Unknown ancestry: ", ancestry)
}

extract_gtf_attr <- function(attribute, key) {
  pattern <- paste0(key, " \"([^\"]+)\"")
  hit <- regexec(pattern, attribute)
  value <- regmatches(attribute, hit)
  vapply(value, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
}

read_gtf_genes <- function(gtf_file) {
  gtf <- fread(gtf_file, sep = "\t", header = FALSE, skip = "#", showProgress = FALSE)
  if (ncol(gtf) < 9) stop("GTF has fewer than 9 columns: ", gtf_file)

  setnames(gtf, c(
    "seqname", "source", "type", "start", "end",
    "score", "strand", "frame", "attribute"
  ))

  gtf <- gtf[type == "gene"]
  gtf[, gene_id := extract_gtf_attr(attribute, "gene_id")]
  gtf[, gene_name := extract_gtf_attr(attribute, "gene_name")]
  gtf
}

make_regions <- function(gene_file, gtf_file, disease, window_bp, gene_id_filter = "") {
  st22 <- fread(gene_file, encoding = "UTF-8")
  names(st22) <- clean_text(gsub("\ufeff", "", names(st22)))

  required <- c("Disease", "Gene ID", "Gene Name", "Model", "Chr")
  missing <- setdiff(required, names(st22))
  if (length(missing)) stop("Missing ST22 columns: ", paste(missing, collapse = ", "))

  st22 <- st22[clean_text(Disease) == clean_text(disease)]
  if (gene_id_filter != "") {
    st22 <- st22[strip_ensg_version(`Gene ID`) == strip_ensg_version(gene_id_filter)]
  }

  if (!nrow(st22)) stop("No ST22 genes found for disease: ", disease)

  st22[, gene_id_clean := strip_ensg_version(`Gene ID`)]
  st22[, Chr := clean_chr(Chr)]

  gtf <- read_gtf_genes(gtf_file)
  gtf[, gene_id_clean := strip_ensg_version(gene_id)]
  gtf[, Chr := clean_chr(seqname)]

  gtf_keep <- gtf[, .(
    gene_id_clean,
    Chr,
    gtf_gene_id = clean_text(gene_id),
    gtf_gene_name = clean_text(gene_name),
    gtf_start = as.integer(start),
    gtf_end = as.integer(end),
    strand = clean_text(strand)
  )]

  annotated <- merge(
    st22,
    gtf_keep,
    by = c("gene_id_clean", "Chr"),
    all.x = TRUE
  )

  if (!"Start" %in% names(annotated) || !"End" %in% names(annotated)) {
    stop("ST22 table must contain Start and End columns for GTF-missing fallback.")
  }

  annotated[, st22_start := as.integer(Start)]
  annotated[, st22_end := as.integer(End)]
  annotated[, has_gtf := !is.na(strand) & !is.na(gtf_start) & !is.na(gtf_end)]

  if (!any(annotated$has_gtf)) {
    warning("No ST22 genes matched GTF for disease ", disease, "; using ST22 Start/End fallback.")
  } else if (any(!annotated$has_gtf)) {
    warning(sum(!annotated$has_gtf), " ST22 genes did not match GTF; using ST22 Start/End fallback for those genes.")
  }

  annotated <- annotated[
    has_gtf |
      (!is.na(st22_start) & !is.na(st22_end))
  ]
  if (!nrow(annotated)) {
    stop("No genes have either a GTF match or usable ST22 Start/End coordinates for disease: ", disease)
  }

  annotated[, center_bp := fifelse(
    has_gtf,
    fifelse(strand == "+", gtf_start, gtf_end),
    as.integer(floor((pmin(st22_start, st22_end) + pmax(st22_start, st22_end)) / 2))
  )]
  annotated[, start_bp := pmax(1L, as.integer(center_bp - window_bp))]
  annotated[, end_bp := as.integer(center_bp + window_bp)]
  annotated[has_gtf == FALSE, `:=`(
    start_bp = pmax(1L, as.integer(pmin(st22_start, st22_end) - window_bp)),
    end_bp = as.integer(pmax(st22_start, st22_end) + window_bp)
  )]
  annotated[, anchor := fifelse(has_gtf, "TSS", "ST22_gene_body_fallback")]
  annotated[, locus_id := clean_id(gene_id_clean)]
  annotated[, region_i := .I]

  unique(
    annotated[, .(
      region_i,
      locus_id,
      Disease,
      Gene_Name = `Gene Name`,
      Gene_ID = `Gene ID`,
      gene = gene_id_clean,
      chr = Chr,
      strand = fifelse(has_gtf, strand, NA_character_),
      center_bp,
      start_bp,
      end_bp,
      anchor
    )],
    by = "locus_id"
  )
}

read_panel_row <- function(panel_tsv, panel_index = NULL, panel_label = NULL) {
  panels <- fread(panel_tsv, showProgress = FALSE)
  names(panels) <- clean_text(gsub("\ufeff", "", names(panels)))

  required <- c("ancestry", "panel_label", "xqtl_list")
  missing <- setdiff(required, names(panels))
  if (length(missing)) stop("Missing panel TSV columns: ", paste(missing, collapse = ", "))

  panels[, ancestry := clean_text(ancestry)]
  panels[, panel_label := clean_text(panel_label)]
  if ("qtl_type" %in% names(panels)) {
    panels[, qtl_type := infer_qtl_type_from_panel(panel_label, qtl_type)]
  } else {
    panels[, qtl_type := infer_qtl_type_from_panel(panel_label)]
  }
  panels[, xqtl_list := clean_text(xqtl_list)]

  if (!is.null(panel_label) && panel_label != "") {
    wanted_panel_label <- clean_text(panel_label)
    panel <- panels[panel_label == wanted_panel_label]
  } else {
    idx <- as.integer(panel_index)
    if (is.na(idx) || idx < 1 || idx > nrow(panels)) {
      stop("--panel-index must be between 1 and ", nrow(panels))
    }
    panel <- panels[idx]
  }

  if (nrow(panel) != 1) {
    stop(
      "Panel selection did not return exactly one row. panel_label=",
      panel_label,
      "; matched rows=",
      nrow(panel)
    )
  }

  panel
}

read_gwas_snps <- function(gwas_file) {
  if (!file.exists(gwas_file)) stop("GWAS file not found: ", gwas_file)

  gwas <- fread(gwas_file, showProgress = FALSE)
  names(gwas) <- clean_text(names(gwas))

  snp_col <- pick_col(gwas, c("SNP", "snp", "rsid", "RSID", "MarkerName"), "GWAS SNP")
  unique(clean_text(gwas[[snp_col]]))
}

prepare_xqtl_list <- function(xqtl_list_file) {
  if (!file.exists(xqtl_list_file)) stop("xQTL list not found: ", xqtl_list_file)

  xqtl <- fread(xqtl_list_file, showProgress = FALSE)
  names(xqtl) <- clean_text(names(xqtl))

  path_col <- pick_col(xqtl, c("PathOfEsd", "path", "esd", "ESD"), "PathOfEsd")

  if (!"ProbeID" %in% names(xqtl)) xqtl[, ProbeID := NA_character_]
  if (!"Gene" %in% names(xqtl)) xqtl[, Gene := NA_character_]

  list_dir <- dirname(xqtl_list_file)

  xqtl[, esd_file := clean_text(get(path_col))]
  xqtl[!file.exists(esd_file), esd_file := file.path(list_dir, esd_file)]
  xqtl <- xqtl[file.exists(esd_file)]

  xqtl[, probe_clean := strip_ensg_version(ProbeID)]
  xqtl[, gene_clean := strip_ensg_version(Gene)]
  xqtl[, probe_ensg := extract_ensg(ProbeID)]
  xqtl[, gene_ensg := extract_ensg(Gene)]
  xqtl[, path_ensg := extract_ensg(basename(esd_file))]

  xqtl
}

same_gene_xqtl_rows <- function(xqtl, gene_id_clean, gene_name) {
  gene_id_clean <- clean_text(gene_id_clean)
  gene_name <- clean_text(gene_name)

  xqtl[
    probe_clean == gene_id_clean |
      gene_clean == gene_id_clean |
      probe_ensg == gene_id_clean |
      gene_ensg == gene_id_clean |
      path_ensg == gene_id_clean |
      has_gene_symbol(ProbeID, gene_name) |
      has_gene_symbol(Gene, gene_name) |
      has_gene_symbol(basename(esd_file), gene_name)
  ]
}

pick_header_col <- function(cols, candidates, label) {
  hit <- candidates[candidates %in% cols]
  if (!length(hit)) {
    stop("Missing ", label, " column. Found: ", paste(cols, collapse = ", "))
  }
  hit[1]
}

read_best_snp_from_esd <- function(esd_file, region, gwas_snps) {
  cols <- tryCatch(names(fread(esd_file, nrows = 0, showProgress = FALSE)), error = function(e) NULL)
  if (is.null(cols)) return(NULL)
  cols <- clean_text(cols)

  chr_col <- pick_header_col(cols, c("Chr", "chr", "CHR"), "ESD chr")
  bp_col <- pick_header_col(cols, c("Bp", "bp", "BP", "pos", "POS"), "ESD bp")
  snp_col <- pick_header_col(cols, c("SNP", "snp", "rsid", "RSID"), "ESD SNP")
  p_col <- pick_header_col(cols, c("p", "P", "pvalue", "Pvalue", "pval"), "ESD p")

  esd <- tryCatch(
    fread(esd_file, select = unique(c(chr_col, bp_col, snp_col, p_col)), showProgress = FALSE),
    error = function(e) NULL
  )
  if (is.null(esd) || !nrow(esd)) return(NULL)
  names(esd) <- clean_text(names(esd))

  esd <- esd[, .(
    target_chr = clean_chr(get(chr_col)),
    target_bp = as.integer(get(bp_col)),
    target_snp = clean_text(get(snp_col)),
    target_p = as.numeric(get(p_col))
  )]

  esd <- esd[
    target_chr == clean_chr(region$chr) &
      target_bp >= as.integer(region$start_bp) &
      target_bp <= as.integer(region$end_bp) &
      target_snp %in% gwas_snps &
      !is.na(target_p) &
      target_p > 0 &
      target_p <= 1
  ]

  if (!nrow(esd)) return(NULL)

  esd[order(target_p)][1]
}

select_targets_same_gene_gwas <- function(regions, panel, gwas_snps) {
  xqtl_list_file <- clean_text(panel$xqtl_list[1])

  if (!file.exists(xqtl_list_file)) {
    out <- copy(regions)
    out[, `:=`(
      target_snp = NA_character_,
      target_p = NA_real_,
      target_chr = NA_character_,
      target_bp = NA_integer_,
      source_esd = NA_character_,
      n_same_gene_esd = 0L,
      error = "missing_xqtl_list"
    )]
    return(out)
  }

  xqtl <- prepare_xqtl_list(xqtl_list_file)
  res <- vector("list", nrow(regions))

  for (i in seq_len(nrow(regions))) {
    region <- regions[i]

    candidates <- same_gene_xqtl_rows(
      xqtl = xqtl,
      gene_id_clean = region$gene,
      gene_name = region$Gene_Name
    )

    if (!nrow(candidates)) {
      res[[i]] <- cbind(
        copy(region),
        data.table(
          target_snp = NA_character_,
          target_p = NA_real_,
          target_chr = NA_character_,
          target_bp = NA_integer_,
          source_esd = NA_character_,
          n_same_gene_esd = 0L,
          error = "no_same_gene_esd"
        )
      )
      next
    }

    hits <- lapply(unique(candidates$esd_file), function(f) {
      h <- read_best_snp_from_esd(f, region, gwas_snps)
      if (is.null(h)) return(NULL)
      h[, source_esd := f]
      h
    })

    hits <- rbindlist(hits, fill = TRUE)

    if (!nrow(hits)) {
      res[[i]] <- cbind(
        copy(region),
        data.table(
          target_snp = NA_character_,
          target_p = NA_real_,
          target_chr = NA_character_,
          target_bp = NA_integer_,
          source_esd = NA_character_,
          n_same_gene_esd = length(unique(candidates$esd_file)),
          error = "no_snp_in_region_and_gwas"
        )
      )
      next
    }

    best <- hits[order(target_p)][1]

    res[[i]] <- cbind(
      copy(region),
      data.table(
        target_snp = best$target_snp,
        target_p = best$target_p,
        target_chr = best$target_chr,
        target_bp = best$target_bp,
        source_esd = best$source_esd,
        n_same_gene_esd = length(unique(candidates$esd_file)),
        error = ""
      )
    )

    if (i %% 100 == 0) {
      message("Processed regions: ", i, " / ", nrow(regions))
      gc()
    }
  }

  rbindlist(res, fill = TRUE)
}

gene_file <- need_arg("--gene-list")
gtf_file <- need_arg("--gtf")
panel_tsv <- need_arg("--panel-tsv")
disease <- need_arg("--disease")
out_dir <- need_arg("--out-dir")
gwas_dir <- get_arg("--gwas-dir", "")
gwas_summary <- get_arg("--gwas-summary", "")
window_bp <- as.integer(get_arg("--window-bp", "1000000"))
panel_index <- get_arg("--panel-index", NULL)
panel_label <- get_arg("--panel-label", NULL)
gene_id_filter <- get_arg("--gene-id", "")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

panel <- read_panel_row(panel_tsv, panel_index = panel_index, panel_label = panel_label)

ancestry <- clean_text(panel$ancestry[1])
panel_label_value <- clean_text(panel$panel_label[1])

if (gwas_summary == "") {
  if (gwas_dir == "") stop("Provide either --gwas-summary or --gwas-dir")
  gwas_summary <- gwas_summary_for(gwas_dir, ancestry, disease)
}

message("Disease: ", disease)
message("Ancestry: ", ancestry)
message("Panel: ", panel_label_value)
message("GWAS: ", gwas_summary)
if (gene_id_filter != "") message("Gene filter: ", gene_id_filter)

regions <- make_regions(
  gene_file = gene_file,
  gtf_file = gtf_file,
  disease = disease,
  window_bp = window_bp,
  gene_id_filter = gene_id_filter
)

regions_file <- file.path(out_dir, paste0(disease, ".", ancestry, ".", clean_id(panel_label_value), ".regions.tsv"))
fwrite(regions, regions_file, sep = "\t")
message("Wrote regions: ", regions_file)

gwas_snps <- read_gwas_snps(gwas_summary)
message("GWAS SNPs loaded: ", length(gwas_snps))

targets <- select_targets_same_gene_gwas(
  regions = regions,
  panel = panel,
  gwas_snps = gwas_snps
)

targets[, `:=`(
  ancestry = ancestry,
  panel_label = panel_label_value,
  source = if ("source" %in% names(panel)) clean_text(panel$source[1]) else NA_character_,
  qtl_type = if ("qtl_type" %in% names(panel)) clean_text(panel$qtl_type[1]) else NA_character_,
  xqtl_list = clean_text(panel$xqtl_list[1]),
  gwas_summary = gwas_summary
)]

out_file <- file.path(
  out_dir,
  paste0(disease, ".", ancestry, ".", clean_id(panel_label_value), ".target_snps.tsv")
)

fwrite(targets, out_file, sep = "\t")

message("Wrote target SNPs: ", out_file)
message("Rows: ", nrow(targets))
message("Rows with target SNP: ", sum(!is.na(targets$target_snp) & targets$target_snp != ""))
message("Errors:")
print(targets[, .N, by = error][order(-N)])