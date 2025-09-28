# Clean environment
rm(list = ls())

# Set working directory
setwd("/Users/kasiadur/Desktop/analysis - part1")

# Load libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(biomaRt)
  library(tidyr)
  library(RCy3)
  library(ashr)
  library(ReactomePA)
})

# Load metadata and counts
meta <- read.delim("gene - metadata.txt", header = TRUE, stringsAsFactors = FALSE)
raw_counts <- read.delim("gene_counts.txt", header = TRUE, check.names = FALSE)
ensembl_ids <- sub("\\..*", "", raw_counts$Row.names)
stopifnot(!anyDuplicated(ensembl_ids))
rownames(raw_counts) <- ensembl_ids
counts <- raw_counts[, -1]

# Subset metadata for neurons from 1.5Mb carriers and control
neurons_meta <- subset(meta, Celltype == "n" & Patient %in% c("d8", "d9", "c11"))
neurons_meta$Disease <- ifelse(grepl("^c", neurons_meta$Patient), "control", "carrier")

# Create group IDs
neurons_meta$Group <- ifelse(
  grepl("^DS", neurons_meta$Row.names),
  gsub("(_[anp]\\d+)_\\d+_S\\d+", "\\1", neurons_meta$Row.names),
  gsub("(_p\\d+)_\\d+$", "\\1", neurons_meta$Row.names)
)

# Collapse replicates
neurons_counts <- counts[, neurons_meta$Row.names]
groups <- unique(neurons_meta$Group)
collapsed_counts <- sapply(groups, function(g) rowSums(neurons_counts[, neurons_meta$Row.names[neurons_meta$Group == g], drop = FALSE]))
neurons_counts_collapsed <- as.data.frame(collapsed_counts)

# Filter genes
keep_genes <- rowSums(neurons_counts_collapsed >= 10) >= 2
neurons_counts_filtered <- neurons_counts_collapsed[keep_genes, ]

# Collapse metadata
neurons_meta_collapsed <- neurons_meta[!duplicated(neurons_meta$Group), ]
rownames(neurons_meta_collapsed) <- neurons_meta_collapsed$Group
neurons_meta_collapsed <- neurons_meta_collapsed[, c("Disease", "Celltype", "Patient")]
neurons_meta_collapsed$Disease <- factor(neurons_meta_collapsed$Disease, levels = c("control", "carrier"))

# Add Sex
neurons_meta_collapsed$Sex <- dplyr::case_when(
  neurons_meta_collapsed$Patient == "d8"  ~ "male",
  neurons_meta_collapsed$Patient == "d9"  ~ "female",
  neurons_meta_collapsed$Patient == "c11" ~ "male"
)
neurons_meta_collapsed$Sex <- factor(neurons_meta_collapsed$Sex)

# DESeq2
dds <- DESeqDataSetFromMatrix(neurons_counts_filtered, neurons_meta_collapsed, design = ~ Sex + Disease)
dds <- DESeq(dds)
res <- lfcShrink(dds, coef = "Disease_carrier_vs_control", type = "ashr")
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$gene_clean <- sub("\\..*", "", res_df$gene)

# Map symbols + Entrez
mart <- useEnsembl("ensembl", "hsapiens_gene_ensembl", mirror = "useast")
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
  filters = "ensembl_gene_id",
  values = res_df$gene_clean,
  mart = mart
)
res_df <- merge(res_df, gene_map, by.x = "gene_clean", by.y = "ensembl_gene_id", all.x = TRUE)

# Significant genes
sig_genes_entrez <- na.omit(res_df$entrezgene_id[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1])

# DEG table for Cytoscape coloring
degs_data <- res_df[, c("entrezgene_id", "log2FoldChange")]
colnames(degs_data) <- c("ID", "log2FC")
degs_data <- na.omit(degs_data)

# Run enrichments
edo_go <- enrichGO(gene = sig_genes_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.05)
edo_kegg <- enrichKEGG(gene = sig_genes_entrez, organism = "hsa", pvalueCutoff = 0.05)
edo_reactome <- enrichPathway(gene = sig_genes_entrez, organism = "human", pvalueCutoff = 0.05)
edo_wikipathways <- enrichWP(gene = sig_genes_entrez, organism = "Homo sapiens", pvalueCutoff = 0.05)

# ✅ Add readable label for Cytoscape node display
edo_go@result$Label <- edo_go@result$Description

# Standardize enrichment output
standardize_enrich <- function(df, source_name) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df <- as.data.frame(df)
  df$Source <- source_name
  required_cols <- c("ID", "Description", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
  missing_cols <- setdiff(required_cols, colnames(df))
  df[missing_cols] <- NA
  df <- df[, c(required_cols, "Source")]
  return(df)
}

# Merge all results
all_enrich <- do.call(rbind, list(
  standardize_enrich(edo_go, "GO"),
  standardize_enrich(edo_kegg, "KEGG"),
  standardize_enrich(edo_reactome, "Reactome"),
  standardize_enrich(edo_wikipathways, "WikiPathways")
))

# Clean GO terms from descriptions
all_enrich$Description <- gsub("^GO:\\d+\\s*", "", all_enrich$Description)
write.csv(all_enrich, "Neurons_all_enrichment_combined.csv", row.names = FALSE)

# Cytoscape visualization
if (nrow(edo_go) > 0) {
  source("cyemapplot_ora.R")
  edo_sim <- pairwise_termsim(edo_go)
  cy.emapplot(
    edo_sim,
    analysis.name = "emapplot-ora-withdata",
    degs_data = degs_data,
    show_category = 30,
    min_edge = 0.3
  )
}
