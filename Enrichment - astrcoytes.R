# Clean environment
rm(list = ls())

# Set working directory
setwd("/Users/kasiadur/Desktop/analysis - part1")

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(EnhancedVolcano)
  library(dplyr)
  library(factoextra)
  library(ggrepel)
  library(gprofiler2)
  library(biomaRt)
  library(cowplot)
  library(stringr)
})

# Load metadata
meta <- read.delim("gene - metadata.txt", header = TRUE, stringsAsFactors = FALSE)

# Load and process gene counts
raw_counts <- read.delim("gene_counts.txt", header = TRUE, check.names = FALSE)
ensembl_ids <- sub("\\..*", "", raw_counts$Row.names)
if (anyDuplicated(ensembl_ids) != 0) stop("Ensembl IDs are not unique.")
rownames(raw_counts) <- ensembl_ids
counts <- raw_counts[, -1]

# Subset metadata for astrocytes from 1.5Mb carriers and control
astro_meta <- subset(meta, Celltype == "a" & Patient %in% c("d8", "d9", "c11"))
astro_meta$Disease <- ifelse(grepl("^c", astro_meta$Patient), "control", "carrier")

# Create group IDs for collapsing
astro_meta$Group <- ifelse(
  grepl("^DS", astro_meta$Row.names),
  gsub("(_[anp]\\d+)_\\d+_S\\d+", "\\1", astro_meta$Row.names),
  gsub("(_p\\d+)_\\d+$", "\\1", astro_meta$Row.names)
)

# Filter and collapse technical replicates
astro_counts <- counts[, astro_meta$Row.names]
replicate_groups <- unique(astro_meta$Group)
collapsed_counts <- matrix(nrow = nrow(astro_counts), ncol = length(replicate_groups))
rownames(collapsed_counts) <- rownames(astro_counts)
colnames(collapsed_counts) <- replicate_groups

for (i in seq_along(replicate_groups)) {
  group <- replicate_groups[i]
  reps <- astro_meta$Row.names[astro_meta$Group == group]
  collapsed_counts[, i] <- rowSums(astro_counts[, reps, drop = FALSE])
}

astro_counts_collapsed <- as.data.frame(collapsed_counts)

# Filter out low-expression genes
keep_genes <- rowSums(astro_counts_collapsed >= 10) >= 2
astro_counts_filtered <- astro_counts_collapsed[keep_genes, ]

# Metadata for collapsed samples
astro_meta_collapsed <- astro_meta[!duplicated(astro_meta$Group), ]
rownames(astro_meta_collapsed) <- astro_meta_collapsed$Group
astro_meta_collapsed <- astro_meta_collapsed[, c("Disease", "Celltype", "Patient")]
astro_meta_collapsed$Disease <- factor(astro_meta_collapsed$Disease, levels = c("control", "carrier"))

# Add sex information manually
astro_meta_collapsed$Sex <- dplyr::case_when(
  astro_meta_collapsed$Patient == "d8" ~ "male",
  astro_meta_collapsed$Patient == "d9" ~ "female",
  astro_meta_collapsed$Patient == "c11" ~ "male" 
)
astro_meta_collapsed$Sex <- factor(astro_meta_collapsed$Sex)


# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = astro_counts_filtered,
                              colData = astro_meta_collapsed,
                              design = ~ Sex + Disease)
dds <- DESeq(dds)

# Shrink log fold changes
if (!requireNamespace("ashr", quietly = TRUE)) install.packages("ashr")
library(ashr)
res <- results(dds)
res <- lfcShrink(dds, coef = "Disease_carrier_vs_control", res = res, type = "ashr")

# Convert result to data frame
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$gene_clean <- sub("\\..*", "", res_df$gene)
res_df <- res_df[order(res_df$padj), ]


# Protein coding genes only - specification 
if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt", update = FALSE)
library(biomaRt)
# Get gene biotype information from Ensembl
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")

biotype_info <- getBM(attributes = c("ensembl_gene_id", "gene_biotype"),
                      filters = "ensembl_gene_id",
                      values = res_df$gene_clean,
                      mart = mart)
# Merge with DE results
res_df <- merge(res_df, biotype_info, by.x = "gene_clean", by.y = "ensembl_gene_id", all.x = TRUE)
# Filter to protein-coding genes only
res_df <- subset(res_df, gene_biotype == "protein_coding")


# Significant gene subsets
sig_p_astro <- subset(res_df, pvalue < 0.05 & !is.na(pvalue))
sig_padj_astro <- subset(res_df, padj < 0.05 & !is.na(padj))
sig_up_astro <- subset(res_df, padj < 0.05 & log2FoldChange > 1)
sig_down_astro <- subset(res_df, padj < 0.05 & log2FoldChange < -1)

# Map Ensembl to HGNC
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt", update = FALSE)
library(biomaRt)

mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = res_df$gene_clean,
                  mart = mart)

res_df <- merge(res_df, gene_map, by.x = "gene_clean", by.y = "ensembl_gene_id", all.x = TRUE)

sig_up_astro$symbol <- res_df$hgnc_symbol[match(sig_up_astro$gene_clean, res_df$gene_clean)]
sig_down_astro$symbol <- res_df$hgnc_symbol[match(sig_down_astro$gene_clean, res_df$gene_clean)]

# ---- Barplots: Top 10 up & downregulated genes ----
top_up <- head(sig_up_astro[order(-sig_up_astro$log2FoldChange), ], 10)
top_up$label <- ifelse(
  is.na(top_up$symbol) | top_up$symbol == "",
  paste0("novel_", substr(top_up$gene_clean, 10, 15)),  # e.g. novel_29087
  top_up$symbol
)
top_up$label <- factor(top_up$label, levels = rev(top_up$label))

top_down <- head(sig_down_astro[order(sig_down_astro$log2FoldChange), ], 10)
top_down$label <- ifelse(
  is.na(top_down$symbol) | top_down$symbol == "",
  paste0("novel_", substr(top_down$gene_clean, 10, 15)),
  top_down$symbol
)
top_down$label <- factor(top_down$label, levels = rev(top_down$label))

pub_theme <- theme_minimal(base_size = 14) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 13, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

p_up <- print(ggplot(top_up, aes(x = label, y = log2FoldChange)) +
  geom_bar(stat = "identity", fill = "#d73027", width = 0.6) +
  coord_flip() +
  pub_theme +
  labs(title = "Upregulated", y = "Log2 Fold Change"))

p_down <- print(ggplot(top_down, aes(x = label, y = log2FoldChange)) +
  geom_bar(stat = "identity", fill = "#4575b4", width = 0.6) +
  coord_flip() +
  pub_theme +
  labs(title = "Downregulated", y = "Log2 Fold Change"))

final_plot <- plot_grid(p_up, p_down, labels = NULL, ncol = 2, align = "hv")
titled_plot <- ggdraw() +
  draw_label("Top 10 Up- and Down-regulated Genes – Astrocytes (1.5Mb vs c11)",
             fontface = 'bold', size = 16, x = 0.5, y = 0.98, hjust = 0.5) +
  draw_plot(final_plot, y = 0, height = 0.95)

print(final_plot)
ggsave("Top10_Astro_DEG_barplots_1.5Mb_vs_c11.pdf", titled_plot, width = 11.7, height = 8.3)

# ---- Pathway enrichment ----
sig_genes <- unique(c(sig_up_astro$gene, sig_down_astro$gene))
sig_genes_clean <- sub("\\..*", "", sig_genes)

gp_res <- gost(
  query = sig_genes_clean,
  organism = "hsapiens",
  sources = c("GO:BP", "REAC", "KEGG", "WP"),
  user_threshold = 0.05,
  correction_method = "fdr",
  significant = TRUE
)

if (!is.null(gp_res) && !is.null(gp_res$result) && nrow(gp_res$result) > 0) {
  gp_df <- as.data.frame(gp_res$result)
  gp_df[] <- lapply(names(gp_df), function(colname) {
    col <- gp_df[[colname]]
    if (is.list(col)) {
      sapply(col, paste, collapse = ", ")
    } else if (colname == "p_value") {
      signif(col, 3)  
    } else if (is.numeric(col)) {
      round(col, 2)   
    } else {
      col
    }
  })
  colnames(gp_df) <- names(gp_res$result)
  
  write.csv(gp_df, "astro_1.5Mb_vs_c11_gprofiler_enrichment.csv", row.names = FALSE)
  
  enrich_theme <- theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11),
      plot.margin = margin(15, 15, 15, 15)
    )
  
  go_bp_only <- subset(gp_df, source == "GO:BP")
  go_bp_only$term_name_wrapped <- str_wrap(go_bp_only$term_name, width = 40)
  top_go <- head(go_bp_only[order(go_bp_only$p_value), ], 10)
  top_go$term_name_wrapped <- factor(top_go$term_name_wrapped, levels = rev(top_go$term_name_wrapped))
  
  go_plot <- print(ggplot(top_go, aes(x = term_name_wrapped, y = -log10(p_value))) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
    coord_flip() +
    enrich_theme +
    labs(title = "Top 10 GO:BP Terms – Astrocytes (1.5Mb vs c11)",
         x = "GO Biological Process", y = "-log10(p-value)"))
  
  ggsave("GO_BP_Top10_Astro_1.5Mb_vs_c11.pdf", go_plot, width = 11.7, height = 8.3)
  
  top_pathways <- head(gp_df[order(gp_df$p_value), ], 20)
  top_pathways$term_name_wrapped <- str_wrap(top_pathways$term_name, width = 40)
  top_pathways$term_name_wrapped <- factor(top_pathways$term_name_wrapped, levels = rev(top_pathways$term_name_wrapped))
  
  pathway_plot <- print(ggplot(top_pathways, aes(x = term_name_wrapped, y = -log10(p_value),
                                           size = intersection_size, color = source)) +
    geom_point() +
    coord_flip() +
    enrich_theme +
    labs(title = "Top 20 Enriched Pathways – Astrocytes (1.5Mb vs c11)",
         x = "Pathway / Term", y = "-log10(p-value)", size = "Gene Count", color = "Source"))
  
  ggsave("Enriched_Pathways_Top20_Astro_1.5Mb_vs_c11.pdf", pathway_plot, width = 11.7, height = 8.3)
} else {
  cat("No enrichment results returned by g:Profiler.\n")
}
