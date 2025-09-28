## ================================
## Top-20 Up/Down tables (HGNC), sex-corrected
## Works for multiple comparisons.
## Outputs only:
##   top20_up_<label>.csv
##   top20_down_<label>.csv
## ================================

setwd("/Users/kasiadur/Desktop/analysis - part1")

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

options(expressions = 5e5)

wd <- "/Users/kasiadur/Desktop/analysis - part1"
meta_file   <- "gene - metadata.txt"
counts_file <- "gene_counts.txt"
setwd(wd)

meta   <- read.delim(meta_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
counts <- read.delim(counts_file, header = TRUE, row.names = 1, check.names = FALSE)

## --- helpers (OFFLINE) ---

id_to_symbol <- function(ids) {
  ids <- as.character(ids)
  sym <- rep(NA_character_, length(ids))
  is_entrez <- grepl("^[0-9]+$", ids)
  is_ensg   <- grepl("^ENSG", ids)
  if (any(is_entrez)) {
    m <- mapIds(org.Hs.eg.db, keys = unique(ids[is_entrez]),
                keytype = "ENTREZID", column = "SYMBOL", multiVals = "first")
    sym[is_entrez] <- unname(m[ids[is_entrez]])
  }
  if (any(is_ensg)) {
    m <- mapIds(org.Hs.eg.db, keys = unique(ids[is_ensg]),
                keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
    sym[is_ensg] <- ifelse(is.na(sym[is_ensg]), unname(m[ids[is_ensg]]), sym[is_ensg])
  }
  left <- which(is.na(sym) | sym == "")
  if (length(left) > 0) sym[left] <- paste0(ids[left], " (unmapped)")
  sym
}

is_protein_coding <- function(ids) {
  ids <- as.character(ids)
  res <- rep(FALSE, length(ids))
  is_entrez <- grepl("^[0-9]+$", ids)
  is_ensg   <- grepl("^ENSG", ids)
  if (any(is_entrez)) {
    gt <- mapIds(org.Hs.eg.db, keys = unique(ids[is_entrez]),
                 keytype = "ENTREZID", column = "GENETYPE", multiVals = "first")
    gt <- unname(gt); names(gt) <- unique(ids[is_entrez])
    res[is_entrez] <- gt[ids[is_entrez]] == "protein-coding"
  }
  if (any(is_ensg)) {
    gt <- mapIds(org.Hs.eg.db, keys = unique(ids[is_ensg]),
                 keytype = "ENSEMBL", column = "GENETYPE", multiVals = "first")
    gt <- unname(gt); names(gt) <- unique(ids[is_ensg])
    res[is_ensg] <- gt[ids[is_ensg]] == "protein-coding"
  }
  res[is.na(res)] <- FALSE
  res
}

fmt_table <- function(df) {
  df %>%
    mutate(
      log2FC = sprintf("%.3f", as.numeric(log2FC)),
      pvalue = ifelse(is.na(pvalue), NA, formatC(pvalue, format = "e", digits = 2)),
      padj   = ifelse(is.na(padj),   NA, formatC(padj,   format = "e", digits = 2))
    )
}

make_groups <- function(x) {
  ifelse(
    grepl("^DS", x),
    gsub("(_[anp]\\d+)_\\d+_S\\d+", "\\1", x),
    gsub("(_[anp]\\d+)_\\d+$", "\\1", x)
  )
}

empty_tbl <- function() {
  fmt_table(tibble::tibble(GeneSymbol=character(), ID=character(),
                           log2FC=numeric(), pvalue=numeric(), padj=numeric()))
}

safe_write_two <- function(up_tbl, down_tbl, label) {
  write.csv(up_tbl,   paste0("top20_up_",   label, ".csv"),   row.names = FALSE)
  write.csv(down_tbl, paste0("top20_down_", label, ".csv"),   row.names = FALSE)
}

## ---------- main ----------
run_top20 <- function(celltype, patients, label, disease_by_patient = NULL) {
  mm <- subset(meta, Celltype == celltype & Patient %in% patients)
  if (nrow(mm) == 0) {
    warning("No samples for ", label, " — writing empty tables.")
    safe_write_two(empty_tbl(), empty_tbl(), label)
    return(invisible(list(up = empty_tbl(), down = empty_tbl())))
  }
  mm$Group <- make_groups(mm$Row.names)
  
  if (!all(mm$Row.names %in% colnames(counts))) {
    warning("Missing columns for ", label, " — writing empty tables.")
    safe_write_two(empty_tbl(), empty_tbl(), label)
    return(invisible(list(up = empty_tbl(), down = empty_tbl())))
  }
  
  sub_counts <- counts[, mm$Row.names, drop = FALSE]
  collapsed_counts <- sapply(unique(mm$Group), function(g) {
    reps <- mm$Row.names[mm$Group == g]
    rowSums(sub_counts[, reps, drop = FALSE])
  })
  collapsed_counts <- as.data.frame(collapsed_counts)
  colnames(collapsed_counts) <- unique(mm$Group)
  
  keep <- rowSums(collapsed_counts >= 10) >= 2
  collapsed_counts <- collapsed_counts[keep, , drop = FALSE]
  
  ## protein-coding only (offline)
  if (nrow(collapsed_counts) == 0) {
    warning("No genes after count filter for ", label, " — writing empty tables.")
    safe_write_two(empty_tbl(), empty_tbl(), label)
    return(invisible(list(up = empty_tbl(), down = empty_tbl())))
  }
  
  pc_mask <- is_protein_coding(rownames(collapsed_counts))
  collapsed_counts <- collapsed_counts[pc_mask, , drop = FALSE]
  if (nrow(collapsed_counts) == 0) {
    warning("No protein-coding genes left for ", label, " — writing empty tables.")
    safe_write_two(empty_tbl(), empty_tbl(), label)
    return(invisible(list(up = empty_tbl(), down = empty_tbl())))
  }
  
  m2 <- mm[!duplicated(mm$Group), c("Group","Disease","Patient")]
  rownames(m2) <- m2$Group
  if (!is.null(disease_by_patient)) {
    if (!all(m2$Patient %in% names(disease_by_patient))) {
      warning("Missing disease mapping for ", label, " — writing empty tables.")
      safe_write_two(empty_tbl(), empty_tbl(), label)
      return(invisible(list(up = empty_tbl(), down = empty_tbl())))
    }
    m2$Disease <- unname(disease_by_patient[m2$Patient])
  }
  
  if (all(m2$Disease %in% c("control","carrier"))) {
    m2$Disease <- factor(m2$Disease, levels = c("control","carrier"))
  } else {
    m2$Disease <- factor(m2$Disease, levels = sort(unique(m2$Disease)))
  }
  
  sex_map <- c(d8="male", d9="female", d36="male", d37="male", d10="female",
               c11="male", c35="female", c1="unknown")
  m2$Sex <- factor(ifelse(m2$Patient %in% names(sex_map), sex_map[m2$Patient], "unknown"))
  
  if (length(unique(m2$Disease)) < 2) {
    warning("Need two groups for DE (", label, ") — writing empty tables.")
    safe_write_two(empty_tbl(), empty_tbl(), label)
    return(invisible(list(up = empty_tbl(), down = empty_tbl())))
  }
  
  dds <- DESeqDataSetFromMatrix(collapsed_counts, m2[, c("Disease","Sex")], design = ~ Sex + Disease)
  dds <- DESeq(dds)
  res <- results(dds)
  
  df <- as.data.frame(res)
  df$ID <- rownames(df)
  df$GeneSymbol <- id_to_symbol(df$ID)
  
  ## thresholds per your layout
  padj_cut <- 0.05
  lab <- tolower(label)
  is_astro      <- (celltype == "a")
  is_15_vs_ctrl <- grepl("1\\.5[_ ]?vs[_ ]?ctrl", lab)
  # Rule:
  # - Astrocytes: ±1 (any comparison)
  # - Any 3.0_vs_ctrl or 1.5_vs_3.0: ±1
  # - Neurons/NPCs 1.5_vs_ctrl: ±0.56
  lfc_cut <- if (is_astro) {
    1
  } else if (grepl("3\\.0[_ ]?vs[_ ]?ctrl", lab) || grepl("1\\.5[_ ]?vs[_ ]?3\\.0", lab)) {
    1
  } else if (is_15_vs_ctrl) {
    0.56
  } else {
    1
  }
  
  up_tbl <- df %>%
    filter(!is.na(padj), padj < padj_cut, log2FoldChange >  lfc_cut) %>%
    arrange(desc(log2FoldChange)) %>%
    slice_head(n = 20) %>%
    transmute(GeneSymbol, ID, log2FC = log2FoldChange, pvalue, padj) %>%
    fmt_table()
  
  down_tbl <- df %>%
    filter(!is.na(padj), padj < padj_cut, log2FoldChange < -lfc_cut) %>%
    arrange(log2FoldChange) %>%
    slice_head(n = 20) %>%
    transmute(GeneSymbol, ID, log2FC = log2FoldChange, pvalue, padj) %>%
    fmt_table()
  
  # Always write, even if 0 rows
  safe_write_two(up_tbl, down_tbl, label)
  message("Wrote: top20_up_", label, " / top20_down_", label,
          "  (padj<", padj_cut, ", |LFC|>", lfc_cut, "; up=", nrow(up_tbl), ", down=", nrow(down_tbl), ")")
  invisible(list(up = up_tbl, down = down_tbl))
}

## ---- RUN ----
run_top20("a", c("d8","d9","c11"),                    "astro_1.5_vs_ctrl")
run_top20("n", c("d8","d9","c11"),                    "neuron_1.5_vs_ctrl")
run_top20("p", c("d8","d9","c11"),                    "npc_1.5_vs_ctrl")
run_top20("p", c("d8","d9","c11","c35","c1"),         "1.5_vs_ctrl")
run_top20("p", c("d36","d37","d10","c11","c35","c1"), "3.0_vs_ctrl")
run_top20("p", c("d8","d9","d36","d37","d10"),        "1.5_vs_3.0",
          disease_by_patient = c(d8="1.5Mb", d9="1.5Mb", d36="3.0Mb", d37="3.0Mb", d10="3.0Mb"))
