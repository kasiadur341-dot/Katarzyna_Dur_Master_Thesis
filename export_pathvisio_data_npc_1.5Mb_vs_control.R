# ==========================================================
# DESeq2 -> PathVisio export (Entrez-only)
# ==========================================================

setwd("/Users/kasiadur/Desktop/analysis - part1")

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ashr)
  library(tibble)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(biomaRt)
})

# Convert gene IDs to trimmed, integer-like Entrez strings; blank -> NA
normalize_entrez <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA
  x[!is.na(x)] <- as.character(as.integer(x[!is.na(x)]))
  x
}

# Optional: Map Entrez -> (Ensembl, HGNC) for coverage stats only (not exported)
map_entrez_to_symbols <- function(entrez_ids) {
  entrez_ids <- unique(as.character(entrez_ids))
  m <- suppressWarnings(AnnotationDbi::select(
    org.Hs.eg.db,
    keys = entrez_ids,
    keytype = "ENTREZID",
    columns = c("ENTREZID", "ENSEMBL", "SYMBOL")
  )) %>%
    dplyr::distinct(ENTREZID, .keep_all = TRUE) %>%
    dplyr::transmute(
      entrezgene_id   = as.character(ENTREZID),
      ensembl_gene_id = as.character(ENSEMBL),
      hgnc_symbol     = as.character(SYMBOL)
    )
  m
}

# -------------------- Core analysis for one cohort --------------------
run_deseq_for_cohort <- function(patients, outfile) {
  
  # 1) Load input tables
  meta   <- read.delim("gene - metadata.txt", header = TRUE, stringsAsFactors = FALSE)
  counts <- read.delim("gene_counts.txt", row.names = 1, check.names = FALSE)
  
  # 2) Keep only NPC samples for selected patients
  npc_meta <- subset(meta, Celltype == "p" & Patient %in% patients)
  
  # 3) Derive replicate Group ID from sample names
  npc_meta$Group <- ifelse(
    grepl("^DS", npc_meta$Row.names),
    gsub("(_[anp]\\d+)_\\d+_S\\d+", "\\1", npc_meta$Row.names),
    gsub("(_p\\d+)_\\d+$", "\\1", npc_meta$Row.names)
  )
  
  # 4) Collapse replicate columns by Group (sum counts)
  npc_counts <- counts[, npc_meta$Row.names, drop = FALSE]
  collapsed_counts <- sapply(unique(npc_meta$Group), function(g) {
    reps <- npc_meta$Row.names[npc_meta$Group == g]
    rowSums(npc_counts[, reps, drop = FALSE])
  })
  counts_collapsed <- as.data.frame(collapsed_counts)
  
  # 5) Gene-level filtering: keep genes with >=10 counts in >=2 samples
  keep_genes <- rowSums(counts_collapsed >= 10) >= 2
  counts_filtered <- counts_collapsed[keep_genes, , drop = FALSE]
  
  # 6) Collapse metadata to one row per Group & set factors
  meta_collapsed <- npc_meta[!duplicated(npc_meta$Group), ]
  rownames(meta_collapsed) <- meta_collapsed$Group
  meta_collapsed <- meta_collapsed[, c("Disease", "Celltype", "Patient")]
  
  # Make sure control is reference level
  meta_collapsed$Disease <- factor(meta_collapsed$Disease, levels = c("control", "carrier"))
  
  # Add Sex factor from Patient (hard-coded map)
  meta_collapsed$Sex <- factor(dplyr::case_when(
    meta_collapsed$Patient %in% c("d8", "d36", "d37", "c11") ~ "male",
    meta_collapsed$Patient %in% c("d9", "d10", "c35")       ~ "female",
    TRUE                                                    ~ "unknown"
  ))
  
  # 7) Build DESeqDataSet, run DESeq2, and shrink LFCs for Disease_carrier_vs_control
  dds <- DESeqDataSetFromMatrix(
    countData = counts_filtered,
    colData   = meta_collapsed,
    design    = ~ Sex + Disease
  )
  dds <- DESeq(dds)
  
  res <- lfcShrink(
    dds,
    coef = "Disease_carrier_vs_control",
    res  = results(dds),
    type = "ashr"
  )
  
  # 8) Add Entrez IDs from rownames, normalized as strings
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene_id_raw") %>%
    dplyr::mutate(entrezgene_id = normalize_entrez(gene_id_raw)) %>%
    tibble::as_tibble()
  
  # 9) Optional: print coverage stats against org.Hs.eg.db (not exported)
  cov <- map_entrez_to_symbols(res_df$entrezgene_id)
  cat(
    "Rows:", nrow(res_df),
    " Entrez present:", sum(!is.na(res_df$entrezgene_id)),
    " Ensembl available:", sum(!is.na(cov$ensembl_gene_id)), "\n"
  )
  
  # 10) Export PathVisio-ready table (ENTREZ ONLY)
  out_tbl <- res_df %>%
    dplyr::transmute(
      entrezgene_id = normalize_entrez(entrezgene_id),
      log2FoldChange, pvalue, padj
    )
  
  write.table(
    out_tbl,
    file      = outfile,
    sep       = ",",
    row.names = FALSE,
    col.names = TRUE,
    quote     = FALSE
  )
  
  # Return full result (with entrezgene_id) for later merging
  res_df
}

# -------------------- Run two cohorts --------------------

deseq_1p5Mb <- run_deseq_for_cohort(
  patients = c("d8", "d9", "c11", "c35", "c1"),
  outfile  = "PathVisio_Full_NPC_1.5Mb_vs_Control.csv"
)

deseq_3p0Mb <- run_deseq_for_cohort(
  patients = c("d36", "d37", "d10", "c11", "c35", "c1"),
  outfile  = "PathVisio_Full_NPC_3.0Mb_vs_Control.csv"
)

# -------------------- Merge cohorts (Entrez-only) --------------------

deseq_1p5Mb <- deseq_1p5Mb %>% dplyr::mutate(source = "1.5Mb")
deseq_3p0Mb <- deseq_3p0Mb %>% dplyr::mutate(source = "3.0Mb")

merged <- dplyr::full_join(
  deseq_1p5Mb %>%
    dplyr::select(
      