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
  library(ggfortify)
})

# Load metadata and counts
meta <- read.delim("gene - metadata.txt", header = TRUE, stringsAsFactors = FALSE)
counts <- read.delim("gene_counts.txt", row.names = 1, check.names = FALSE)

# Filter metadata for NPCs (celltype 'p') for 3.0Mb carriers vs controls
npc_meta <- subset(meta, Celltype == "p" & Patient %in% c("d36", "d37", "d10", "c11", "c35", "c1"))

# Create group IDs for collapsing
npc_meta$Group <- ifelse(
  grepl("^DS", npc_meta$Row.names),
  gsub("(_[anp]\\d+)_\\d+_S\\d+", "\\1", npc_meta$Row.names),
  gsub("(_p\\d+)_\\d+$", "\\1", npc_meta$Row.names)
)

# Filter count matrix
npc_counts <- counts[, npc_meta$Row.names]

# Collapse technical replicates
replicate_groups <- unique(npc_meta$Group)
collapsed_counts <- matrix(nrow = nrow(npc_counts), ncol = length(replicate_groups))
rownames(collapsed_counts) <- rownames(npc_counts)
colnames(collapsed_counts) <- replicate_groups

for (i in seq_along(replicate_groups)) {
  group <- replicate_groups[i]
  reps <- npc_meta$Row.names[npc_meta$Group == group]
  collapsed_counts[, i] <- rowSums(npc_counts[, reps, drop = FALSE])
}

npc_counts_collapsed <- as.data.frame(collapsed_counts)

# Filter lowly expressed genes
keep_genes <- rowSums(npc_counts_collapsed >= 10) >= 2
npc_counts_filtered <- npc_counts_collapsed[keep_genes, ]

# Create metadata
npc_meta_collapsed <- npc_meta[!duplicated(npc_meta$Group), ]
rownames(npc_meta_collapsed) <- npc_meta_collapsed$Group
npc_meta_collapsed <- npc_meta_collapsed[, c("Disease", "Celltype", "Patient")]
npc_meta_collapsed$Disease <- factor(npc_meta_collapsed$Disease, levels = c("control", "carrier"))

# Add Sex manually based on Patient ID
npc_meta_collapsed$Sex <- dplyr::case_when(
  npc_meta_collapsed$Patient == "d36" ~ "male",
  npc_meta_collapsed$Patient == "d37" ~ "male",
  npc_meta_collapsed$Patient == "d10" ~ "female",
  npc_meta_collapsed$Patient == "c11" ~ "male",
  npc_meta_collapsed$Patient == "c35" ~ "female",
  npc_meta_collapsed$Patient == "c1"  ~ "unknown"
)
npc_meta_collapsed$Sex <- factor(npc_meta_collapsed$Sex)


# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = npc_counts_filtered,
                              colData = npc_meta_collapsed,
                              design = ~ Sex + Disease)
dds <- DESeq(dds)

# Shrink LFC
if (!requireNamespace("ashr", quietly = TRUE)) install.packages("ashr")
library(ashr)
res <- results(dds)
res <- lfcShrink(dds, coef = "Disease_carrier_vs_control", res = res, type = "ashr")

# Convert to data frame
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]

# Install and load biomaRt if needed
if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")
library(biomaRt)

# Connect to Ensembl BioMart
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")

# Prepare Entrez IDs for mapping
entrez_ids <- as.character(rownames(res_df))

# Split Entrez IDs into chunks of 1000
chunks <- split(entrez_ids, ceiling(seq_along(entrez_ids) / 1000))

# Query biomaRt in chunks
id_map_list <- lapply(chunks, function(ids) {
  tryCatch({
    getBM(
      attributes = c("entrezgene_id", "ensembl_gene_id", "hgnc_symbol"),
      filters = "entrezgene_id",
      values = ids,
      mart = mart
    )
  }, error = function(e) {
    message("Chunk failed: ", conditionMessage(e))
    return(NULL)
  })
})

# Combine results
id_map <- do.call(rbind, id_map_list)

# Merge with res_df
res_df$entrezgene_id <- rownames(res_df)
res_df <- merge(res_df, id_map, by = "entrezgene_id", all.x = TRUE)

# Save full DGE table with Ensembl & HGNC (optional)
write.csv(res_df[, c("ensembl_gene_id", "log2FoldChange", "pvalue", "padj", "hgnc_symbol")],
          "npc_dge_3.0Mb_with_ensembl.csv", row.names = FALSE)

# Clean version without HGNC
dge_clean_3 <- res_df[, c("ensembl_gene_id", "log2FoldChange", "pvalue", "padj")]
write.csv(dge_clean_3, "npc_dge_3.0Mb_clean_no_HGNC.csv", row.names = FALSE)


# Load both DGE tables
dge_1.5 <- read.csv("npc_dge_clean_no_HGNC.csv")
dge_3.0 <- read.csv("npc_dge_3.0Mb_clean_no_HGNC.csv")

# Rename columns to identify comparison
colnames(dge_1.5)[2:4] <- c("log2FC_1.5Mb", "pval_1.5Mb", "padj_1.5Mb")
colnames(dge_3.0)[2:4] <- c("log2FC_3.0Mb", "pval_3.0Mb", "padj_3.0Mb")

# 1. Ensure Ensembl IDs are character type
dge_1.5$ensembl_gene_id <- as.character(dge_1.5$ensembl_gene_id)
dge_3.0$ensembl_gene_id <- as.character(dge_3.0$ensembl_gene_id)

# 2. Remove NA Ensembl IDs
dge_1.5 <- dge_1.5[!is.na(dge_1.5$ensembl_gene_id), ]
dge_3.0 <- dge_3.0[!is.na(dge_3.0$ensembl_gene_id), ]

# 3. Optionally remove duplicated Ensembl IDs (if present)
dge_1.5 <- dge_1.5[!duplicated(dge_1.5$ensembl_gene_id), ]
dge_3.0 <- dge_3.0[!duplicated(dge_3.0$ensembl_gene_id), ]

# 4. Merge safely
merged_dge <- merge(dge_1.5, dge_3.0, by = "ensembl_gene_id", all = TRUE)


# Save merged file
write.csv(merged_dge, "npc_dge_merged_1.5_vs_3.0_vs_control.csv", row.names = FALSE, quote = FALSE)

# Significance filters
sig_p_npc_3 <- subset(res_df, pvalue < 0.05 & !is.na(pvalue))
sig_padj_npc_3 <- subset(res_df, padj < 0.05 & !is.na(padj))
sig_up_npc_3 <- subset(res_df, padj < 0.05 & log2FoldChange > 1)
sig_down_npc_3 <- subset(res_df, padj < 0.05 & log2FoldChange < -1)

# Save tables
write.csv(sig_p_npc_3, "npc_3.0_vs_control_sig_pval.csv", row.names = FALSE)
write.csv(sig_padj_npc_3, "npc_3.0_vs_control_sig_adj_pval.csv", row.names = FALSE)
write.csv(sig_up_npc_3, "npc_3.0_vs_control_upregulated.csv", row.names = FALSE)
write.csv(sig_down_npc_3, "npc_3.0_vs_control_downregulated.csv", row.names = FALSE)

# Summary
cat("Gene count summary for NPCs (3.0 vs control):\n")
cat("1. Raw p-value < 0.05: ", nrow(sig_p_npc_3), "genes\n")
cat("2. Adjusted p-value < 0.05: ", nrow(sig_padj_npc_3), "genes\n")
cat("3. Upregulated (log2FC > 1 & adj.p < 0.05): ", nrow(sig_up_npc_3), "genes\n")
cat("4. Downregulated (log2FC < -1 & adj.p < 0.05): ", nrow(sig_down_npc_3), "genes\n")

# Create summary row and add to global summary list
summary_row <- data.frame(
  Comparison = "3.0_vs_Control",
  Significant_padj_0.05 = nrow(sig_padj_npc),
  Upregulated_log2FC_gt1 = nrow(sig_up_npc),
  Downregulated_log2FC_lt_neg1 = nrow(sig_down_npc)
)

if (!exists("summary_list")) summary_list <- list()
summary_list[["3.0Mb_vs_Control"]] <- summary_row


# Volcano Plot
plot(res_df$log2FoldChange, -log10(res_df$pvalue),
     pch = 20, col = "gray",
     main = "Volcano Plot − NPCs (3.0Mb vs Control)",
     xlab = "Log2 Fold Change",
     ylab = "−log10(p−value)",
     xlim = c(-6, 6))
abline(h = -log10(0.05), col = "darkgray", lty = 2)
with(subset(res_df, padj < 0.05 & log2FoldChange < -1),
     points(log2FoldChange, -log10(pvalue), col = "blue", pch = 20))
with(subset(res_df, padj < 0.05 & log2FoldChange > 1),
     points(log2FoldChange, -log10(pvalue), col = "red", pch = 20))

# Boxplots
boxplot(log(npc_counts_filtered + 1),
        main = "Raw log−counts (Collapsed) − 3.0 Vs. Control",
        ylab = "log1p(Counts)",
        las = 2, col = "gray")

norm_counts <- counts(dds, normalized = TRUE)
boxplot(log(norm_counts + 1),
        main = "Normalized log−counts − 3.0 Vs. Control",
        ylab = "log1p(Normalized Counts)",
        las = 2, col = "gray")

# ---- HEATMAP ----
library(pheatmap); library(RColorBrewer)

# --- choose genes for the heatmap (top by adjusted p-value; fallback = most variable) ---
top_n <- 100
res_df_nona <- res_df[!is.na(res_df$padj), ]
res_df_nona <- res_df_nona[order(res_df_nona$padj), ]
top_genes <- head(res_df_nona$gene, top_n)

if (length(top_genes) == 0) {
  if (!requireNamespace("matrixStats", quietly = TRUE)) install.packages("matrixStats")
  library(matrixStats)
  v_all <- vst(dds, blind = TRUE)
  vars <- rowVars(assay(v_all))
  top_genes <- names(sort(vars, decreasing = TRUE))[seq_len(min(top_n, length(vars)))]
}

v <- vst(dds, blind = TRUE)
vmat <- assay(v)[top_genes, , drop = FALSE]

cap <- 10
vmat_cap <- pmin(vmat, cap)

pal <- colorRampPalette(c("white", "#c6dbef", "#6baed6", "#2171b5", "#08306b"))(100)
bk  <- seq(0, cap, length.out = 101)

annotation_col <- data.frame(Group = npc_meta_collapsed$Disease)
rownames(annotation_col) <- rownames(npc_meta_collapsed)
annotation_col <- annotation_col[colnames(vmat_cap), , drop = FALSE]

pheatmap(
  vmat_cap,
  color = pal, breaks = bk,
  cluster_rows = TRUE, cluster_cols = TRUE,
  annotation_col = annotation_col,
  show_colnames = TRUE, show_rownames = FALSE,
  fontsize_col = 10,
  legend_breaks = c(0,2,4,6,8,10),
  legend_labels = c("0","2","4","6","8","10"),
  main = "Heatmap of Top DE Genes — 3.0Mb vs Control (VST expression, log2-like scale, capped at 10)",
  border_color = NA,
  use_raster = TRUE,
  filename = "3.0Mb vs Control_heatmap.pdf",
  width = 8, height = 10
)

# PCA
library(ggfortify)

pca <- prcomp(t(log1p(norm_counts)))

group_info <- npc_meta_collapsed$Disease[match(rownames(t(norm_counts)), rownames(npc_meta_collapsed))]

autoplot(pca, data = data.frame(group = group_info), colour = 'group', label = TRUE, label.size = 3) +
  scale_color_manual(values = c("control" = "#4575b4", "carrier" = "#d73027")) +
  ggtitle("PCA − 3.0Mb vs Control") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


# Elbow Plot
data_scaled <- t(scale(log1p(norm_counts)))
elbow_plot <- fviz_nbclust(data_scaled, kmeans, method = "wss", k.max = 6) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Elbow Plot – NPCs (3.0Mb vs control)",
    x = "Number of Clusters (k)",
    y = "Within-Cluster Sum of Squares"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 12)
  )
print(elbow_plot)

# Dendrogram 
library(ggplot2)
library(ggdendro)

# data → clustering
norm_counts <- counts(dds, normalized = TRUE)
hc <- hclust(dist(t(log1p(norm_counts))), method = "ward.D2")
dd <- dendro_data(hc)

y_top <- 240  # fixed scale for all comparisons

p_dend <- ggplot() +
  geom_segment(data = dd$segments,
               aes(x = x, y = y, xend = xend, yend = yend)) +
  labs(title = "Sample Clustering Dendrogram – 3.0 vs Control",
       x = NULL, y = "Height") +
  scale_y_continuous(
    limits = c(0, y_top),
    breaks = seq(0, y_top, by = 20),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    breaks = dd$labels$x,
    labels = dd$labels$label,
    expand = c(0.01, 0.02)
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid   = element_blank(),
    panel.border = element_blank(),
    axis.line.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
    plot.title   = element_text(hjust = 0.5, face = "bold")
  )

ggsave("3.0 vs Control_dendrogram.pdf", plot = p_dend, width = 8, height = 6, units = "in")

# ===========================================
# CYTOSCAPE FILE
# ===========================================
# Prepare Cytoscape input for 3.0 vs Control
cyto_3.0 <- res_df[, c("hgnc_symbol", "log2FoldChange", "pvalue", "padj")]
colnames(cyto_3.0) <- c("GeneSymbol", "log2FC", "pvalue", "padj")

# Remove rows without GeneSymbol
cyto_3.0 <- cyto_3.0[!is.na(cyto_3.0$GeneSymbol) & cyto_3.0$GeneSymbol != "", ]

# Save file
write.csv(cyto_3.0, "cytoscape_input_3.0_vs_control.csv", row.names = FALSE)
cat("Cytoscape input file created: cytoscape_input_3.0_vs_control.csv\n")
