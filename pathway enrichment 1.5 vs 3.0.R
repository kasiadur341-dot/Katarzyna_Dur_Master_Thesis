# -------- Pathway Enrichment ONLY: NPCs (1.5Mb vs 3.0Mb) --------
# Clean environment
rm(list = ls())

# Set working directory
setwd("/Users/kasiadur/Desktop/analysis - part1")

# Load required libraries
suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(gprofiler2)
  library(biomaRt)
  library(stringr)
  if (!requireNamespace("ashr", quietly = TRUE)) install.packages("ashr")
  library(ashr)
})

# ----- Load metadata and gene counts
meta <- read.delim("gene - metadata.txt", header = TRUE, stringsAsFactors = FALSE)
raw_counts <- read.delim("gene_counts.txt", header = TRUE, check.names = FALSE)

# Ensure unique Ensembl IDs (strip version)
ensembl_ids <- sub("\\..*", "", raw_counts$Row.names)
if (anyDuplicated(ensembl_ids) != 0) stop("Ensembl IDs are not unique.")
rownames(raw_counts) <- ensembl_ids
counts <- raw_counts[, -1, drop = FALSE]

# ----- Select NPCs: 1.5Mb vs 3.0Mb carriers (no controls)
# Patients: 1.5Mb = d8, d9 ; 3.0Mb = d36, d37, d10
npc_meta <- subset(meta, Celltype == "p" & Patient %in% c("d8","d9","d36","d37","d10"))
npc_meta$Disease <- ifelse(npc_meta$Patient %in% c("d8","d9"), "1.5Mb", "3.0Mb")

# Create group IDs to collapse technical replicates
npc_meta$Group <- ifelse(
  grepl("^DS", npc_meta$Row.names),
  gsub("(_[anp]\\d+)_\\d+_S\\d+", "\\1", npc_meta$Row.names),
  gsub("(_p\\d+)_\\d+$", "\\1", npc_meta$Row.names)
)

# Subset counts and collapse replicates
stopifnot(all(npc_meta$Row.names %in% colnames(counts)))
npc_counts <- counts[, npc_meta$Row.names, drop = FALSE]
replicate_groups <- unique(npc_meta$Group)

collapsed_counts <- sapply(replicate_groups, function(g) {
  reps <- npc_meta$Row.names[npc_meta$Group == g]
  rowSums(npc_counts[, reps, drop = FALSE])
})
npc_counts_collapsed <- as.data.frame(collapsed_counts)

# Filter low-expression genes (>=10 counts in at least 2 samples)
keep_genes <- rowSums(npc_counts_collapsed >= 10) >= 2
npc_counts_filtered <- npc_counts_collapsed[keep_genes, , drop = FALSE]

# Collapsed metadata
npc_meta_collapsed <- npc_meta[!duplicated(npc_meta$Group), ]
rownames(npc_meta_collapsed) <- npc_meta_collapsed$Group
npc_meta_collapsed <- npc_meta_collapsed[, c("Disease", "Celltype", "Patient")]
npc_meta_collapsed$Disease <- factor(npc_meta_collapsed$Disease, levels = c("1.5Mb","3.0Mb"))

# Add Sex per patient
npc_meta_collapsed$Sex <- dplyr::case_when(
  npc_meta_collapsed$Patient == "d8"  ~ "male",
  npc_meta_collapsed$Patient == "d9"  ~ "female",
  npc_meta_collapsed$Patient == "d36" ~ "male",
  npc_meta_collapsed$Patient == "d37" ~ "male",
  npc_meta_collapsed$Patient == "d10" ~ "female",
  TRUE ~ "unknown"
)
npc_meta_collapsed$Sex <- factor(npc_meta_collapsed$Sex)

# ----- Differential expression (just to get gene list)
dds <- DESeqDataSetFromMatrix(
  countData = npc_counts_filtered,
  colData   = npc_meta_collapsed,
  design    = ~ Sex + Disease
)
dds <- DESeq(dds)

# Shrink logFC for 3.0Mb vs 1.5Mb
res <- lfcShrink(dds, contrast = c("Disease","3.0Mb","1.5Mb"), type = "ashr")
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$gene_clean <- sub("\\..*", "", res_df$gene)
res_df <- res_df[order(res_df$padj), ]

# ----- Restrict to protein-coding genes (biomaRt)
if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt", update = FALSE)
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")

biotype_info <- getBM(
  attributes = c("ensembl_gene_id", "gene_biotype"),
  filters    = "ensembl_gene_id",
  values     = res_df$gene_clean,
  mart       = mart
)
res_df <- merge(res_df, biotype_info, by.x = "gene_clean", by.y = "ensembl_gene_id", all.x = TRUE)
res_df <- subset(res_df, gene_biotype == "protein_coding")

# Significant gene set for enrichment
sig_up   <- subset(res_df, padj < 0.05 & log2FoldChange >  1)
sig_down <- subset(res_df, padj < 0.05 & log2FoldChange < -1)

sig_genes       <- unique(c(sig_up$gene, sig_down$gene))
sig_genes_clean <- sub("\\..*", "", sig_genes)

if (length(sig_genes_clean) < 5) {
  stop("Not enough significant genes to run pathway enrichment.")
}

# ----- g:Profiler pathway enrichment (Reactome, KEGG, WikiPathways only)
comp_label    <- "NPCs (1.5Mb vs 3.0Mb)"
top_n_terms   <- 20
fdr_threshold <- 0.05

gp_res <- gprofiler2::gost(
  query             = sig_genes_clean,
  organism          = "hsapiens",
  sources           = c("REAC","KEGG","WP"),
  user_threshold    = fdr_threshold,
  correction_method = "fdr",
  significant       = TRUE
)

if (is.null(gp_res) || is.null(gp_res$result) || nrow(gp_res$result) == 0) {
  stop("No significant pathway enrichment terms returned.")
}

gp_df <- as.data.frame(gp_res$result)
gp_df <- subset(gp_df, source %in% c("REAC","KEGG","WP"))
gp_df$p_value <- as.numeric(gp_df$p_value)
gp_df$source  <- factor(gp_df$source,
                        levels = c("REAC","KEGG","WP"),
                        labels = c("Reactome","KEGG","WikiPathways"))

# Flatten list-columns for CSV
list_cols <- vapply(gp_df, is.list, logical(1))
if (any(list_cols)) {
  gp_df[list_cols] <- lapply(gp_df[list_cols], function(col)
    vapply(col, function(x) paste(x, collapse = ", "), character(1))
  )
}
gp_df$p_value <- as.numeric(gp_df$p_value)

# Save full enrichment table
write.csv(gp_df, "npc_1.5Mb_vs_3.0Mb_pathway_enrichment.csv", row.names = FALSE)

# ----- Pathway dot plot (Top N by FDR p-value)
plot_df <- head(gp_df[order(gp_df$p_value, gp_df$term_name), ], top_n_terms)
plot_df$term_name_wrapped <- stringr::str_wrap(plot_df$term_name, width = 50)
plot_df$term_name_wrapped <- factor(plot_df$term_name_wrapped,
                                    levels = rev(plot_df$term_name_wrapped))

enrich_theme <- theme_minimal(base_size = 14) +
  theme(
    axis.text.y  = element_text(size = 12, face = "bold"),
    axis.text.x  = element_text(size = 12),
    axis.title   = element_text(size = 14, face = "bold"),
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 11),
    plot.margin  = margin(15, 15, 15, 15)
  )

pathway_plot <- ggplot(
  plot_df,
  aes(x = term_name_wrapped, y = -log10(p_value),
      size = intersection_size, color = source)
) +
  geom_point() +
  coord_flip() +
  enrich_theme +
  labs(
    title = paste0("Pathway Enrichment – ", comp_label),
    x = "Pathway",
    y = expression(-log[10](FDR~p~value)),
    size = "Gene Count",
    color = "Source"
  )

print(pathway_plot)
ggsave("Pathway_Enrichment_Top20_NPC_1.5Mb_vs_3.0Mb.png",
       pathway_plot, width = 11.7, height = 8.3, dpi = 300)
ggsave("Pathway_Enrichment_Top20_NPC_1.5Mb_vs_3.0Mb.pdf",
       pathway_plot, width = 11.7, height = 8.3)
# ---------------------------------------------------------------------------
