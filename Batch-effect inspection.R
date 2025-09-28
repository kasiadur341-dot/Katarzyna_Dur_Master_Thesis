## =====================================================
## Batch-effect inspection for 6 comparisons
## =====================================================

# Clean env
rm(list = ls())

# ---- Settings ----
main_dir <- "/Users/kasiadur/Desktop/analysis - part1"
out_dir  <- file.path(main_dir, "batch_inspection")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- Libraries ----
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
})
if (!requireNamespace("vegan", quietly = TRUE)) install.packages("vegan")
suppressPackageStartupMessages(library(vegan))

# Optional: clear stuck graphics device (RStudio quirk)
while (!is.null(dev.list())) dev.off()

# ---- Load data ----
meta   <- read.delim(file.path(main_dir, "gene - metadata.txt"),
                     header = TRUE, stringsAsFactors = FALSE)
counts <- read.delim(file.path(main_dir, "gene_counts.txt"),
                     row.names = 1, check.names = FALSE)

# ---- Helper: build & run batch check for one subset ----
run_batch_check <- function(celltype, patients, comp_name,
                            disease_relabel = NULL) {
  message("\n=== Running batch check for: ", comp_name, " ===")
  
  # 1) Subset metadata (by celltype + patients)
  this_meta <- subset(meta, Celltype == celltype & Patient %in% patients)
  
  # If re-labelling disease (e.g., "1.5Mb" vs "3.0Mb"), apply here
  if (!is.null(disease_relabel)) {
    # disease_relabel is a function that maps Patient -> new disease label
    this_meta$Disease <- disease_relabel(this_meta$Patient)
  }
  
  # 2) Collapse technical replicates exactly like your pipeline
  this_meta$Group <- ifelse(
    grepl("^DS", this_meta$Row.names),
    gsub("(_[anp]\\d+)_\\d+_S\\d+", "\\1", this_meta$Row.names),
    gsub("(_p\\d+)_\\d+$", "\\1", this_meta$Row.names)
  )
  # Counts subset
  stopifnot(all(this_meta$Row.names %in% colnames(counts)))
  this_counts <- counts[, this_meta$Row.names, drop = FALSE]
  
  # Collapse
  collapsed_counts <- sapply(unique(this_meta$Group), function(g) {
    reps <- this_meta$Row.names[this_meta$Group == g]
    rowSums(this_counts[, reps, drop = FALSE])
  })
  collapsed_counts <- as.data.frame(collapsed_counts)
  
  # 3) Filter low counts (>=10 in at least 2 samples)
  keep_genes <- rowSums(collapsed_counts >= 10) >= 2
  counts_filtered <- collapsed_counts[keep_genes, , drop = FALSE]
  
  # 4) Collapsed metadata
  this_meta_collapsed <- this_meta[!duplicated(this_meta$Group), ]
  rownames(this_meta_collapsed) <- this_meta_collapsed$Group
  this_meta_collapsed <- this_meta_collapsed[, c("Disease", "Celltype", "Patient")]
  
  # Order Disease levels (put the 1st unique as reference)
  if (is.null(disease_relabel)) {
    this_meta_collapsed$Disease <- factor(this_meta_collapsed$Disease,
                                          levels = unique(c("control","carrier",
                                                            this_meta_collapsed$Disease)))
  } else {
    # Keep the order encountered
    this_meta_collapsed$Disease <- factor(this_meta_collapsed$Disease,
                                          levels = unique(this_meta_collapsed$Disease))
  }
  
  # 5) Add Sex (manual mapping you provided)
  this_meta_collapsed$Sex <- dplyr::case_when(
    this_meta_collapsed$Patient == "d8"  ~ "male",
    this_meta_collapsed$Patient == "d9"  ~ "female",
    this_meta_collapsed$Patient == "d36" ~ "male",
    this_meta_collapsed$Patient == "d37" ~ "male",
    this_meta_collapsed$Patient == "d10" ~ "female",
    this_meta_collapsed$Patient == "c11" ~ "male",
    this_meta_collapsed$Patient == "c35" ~ "female",
    this_meta_collapsed$Patient == "c1"  ~ "unknown",
    TRUE ~ NA_character_
  )
  this_meta_collapsed$Sex <- factor(this_meta_collapsed$Sex)
  
  # 6) Dynamic design (drop constant factors so DESeq2 doesn't error)
  design_vars <- character(0)
  if (nlevels(this_meta_collapsed$Sex)     > 1) design_vars <- c(design_vars, "Sex")
  if (nlevels(this_meta_collapsed$Disease) > 1) design_vars <- c(design_vars, "Disease")
  
  design_formula <- if (length(design_vars) == 0) {
    ~ 1
  } else {
    as.formula(paste("~", paste(design_vars, collapse = " + ")))
  }
  
  # 7) Build DESeq2 object (no full DE fit needed for PCA)
  dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                                colData = this_meta_collapsed,
                                design = design_formula)
  dds <- estimateSizeFactors(dds)
  
  # 8) Add DiffRound from sample IDs (e.g., DS..._a1 -> "a1")
  DiffRound <- sub(".*_([anp]\\d+)$", "\\1", colnames(dds))
  colData(dds)$DiffRound <- factor(DiffRound)
  
  # Sanity table
  tab_DR <- with(as.data.frame(colData(dds)), table(DiffRound, Disease))
  print(tab_DR)
  
  # 9) PCA (VST), colored by Disease + DiffRound
  vsd <- vst(dds, blind = FALSE)
  p_pca <- plotPCA(vsd, intgroup = c("Disease","DiffRound")) +
    ggtitle(paste("PCA -", comp_name, "(Disease + DiffRound)"))
  ggsave(file.path(out_dir, paste0(comp_name, "_PCA.pdf")),
         p_pca, width = 7, height = 5, useDingbats = FALSE)
  
  # 10) PERMANOVA (DiffRound + Disease + Sex as available)
  dist_mat <- dist(t(assay(vsd)))
  
  # Build PERMANOVA formula mirroring available columns
  perm_terms <- c()
  if ("DiffRound" %in% colnames(colData(dds))) perm_terms <- c(perm_terms, "DiffRound")
  if ("Disease"   %in% colnames(colData(dds)) && nlevels(colData(dds)$Disease) > 1) perm_terms <- c(perm_terms, "Disease")
  if ("Sex"       %in% colnames(colData(dds)) && nlevels(colData(dds)$Sex)     > 1) perm_terms <- c(perm_terms, "Sex")
  
  perm_formula <- if (length(perm_terms) == 0) as.formula("dist_mat ~ 1")
  else as.formula(paste("dist_mat ~", paste(perm_terms, collapse = " + ")))
  
  perm <- adonis2(perm_formula,
                  data = as.data.frame(colData(dds)),
                  permutations = 999, by = "margin")
  capture.output(perm, file = file.path(out_dir, paste0(comp_name, "_PERMANOVA.txt")))
  print(perm)
  
  # Console p-values
  get_p <- function(term) if (term %in% rownames(perm)) signif(perm$`Pr(>F)`[term], 3) else NA
  cat(sprintf("%s — PERMANOVA p-values | DiffRound: %s | Disease: %s | Sex: %s\n",
              comp_name, get_p("DiffRound"), get_p("Disease"), get_p("Sex")))
}

## ---- Run all 6 comparisons ----

# 1) Astrocytes: 1.5Mb carriers vs controls
run_batch_check("a", c("d8","d9","c11"), "Astrocytes")

# 2) Neurons: 1.5Mb carriers vs controls
run_batch_check("n", c("d8","d9","c11"), "Neurons")

# 3) NPCs: 1.5Mb carriers vs controls
run_batch_check("p", c("d8","d9","c11"), "NPCs")

# 4) NPCs: all 1.5Mb carriers vs all controls
run_batch_check("p", c("d8","d9","c11","c35","c1"), "NPCs_1.5_vs_Control")

# 5) NPCs: 3.0Mb carriers vs controls
run_batch_check("p", c("d36","d37","d10","c11","c35","c1"), "NPCs_3.0_vs_Control")

# 6) NPCs: 1.5Mb vs 3.0Mb (no controls) — relabel Disease by Patient
relabel_1v3 <- function(patient_vec) ifelse(patient_vec %in% c("d8","d9"), "1.5Mb", "3.0Mb")
run_batch_check("p", c("d8","d9","d36","d37","d10"), "NPCs_1.5_vs_3.0",
                disease_relabel = relabel_1v3)

cat("\nBatch-effect inspection complete.\nOutputs saved to:\n", out_dir, "\n")
