# Clean environment
rm(list = ls())

# Set working directory
setwd("/Users/kasiadur/Desktop/analysis - part1")

# Load required libraries
suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(EnhancedVolcano)
  library(dplyr)
  library(factoextra)
  library(ggrepel)
  library(gprofiler2)
  library(gridExtra)
  library(biomaRt)
  library(cowplot)
  library(stringr)
})

# Load metadata and gene counts
meta <- read.delim("gene - metadata.txt", header = TRUE, stringsAsFactors = FALSE)
raw_counts <- read.delim("gene_counts.txt", header = TRUE, check.names = FALSE)
ensembl_ids <- sub("\\..*", "", raw_counts$Row.names)
if (anyDuplicated(ensembl_ids) != 0) stop("Ensembl IDs are not unique.")
rownames(raw_counts) <- ensembl_ids
counts <- raw_counts[, -1, drop = FALSE]

# Filter metadata for NPCs: 3.0Mb carriers vs controls
npc_meta <- subset(meta, Celltype == "p" & Patient %in% c("d36", "d37", "d10", "c11", "c35", "c1"))
npc_meta$Disease <- ifelse(grepl("^c", npc_meta$Patient), "control", "carrier")

# Create group IDs
npc_meta$Group <- ifelse(
  grepl("^DS", npc_meta$Row.names),
  gsub("(_[anp]\\d+)_\\d+_S\\d+", "\\1", npc_meta$Row.names),
  gsub("(_p\\d+)_\\d+$", "\\1", npc_meta$Row.names)
)

# Subset and collapse technical replicates
stopifnot(all(npc_meta$Row.names %in% colnames(counts)))
npc_counts <- counts[, npc_meta$Row.names, drop = FALSE]
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

# Filter low-expression genes
keep_genes <- rowSums(npc_counts_collapsed >= 10) >= 2
npc_counts_filtered <- npc_counts_collapsed[keep_genes, , drop = FALSE]

# Metadata for collapsed samples
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

# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = npc_counts_filtered,
                              colData   = npc_meta_collapsed,
                              design    = ~ Sex + Disease)
dds <- DESeq(dds)

# Shrink log fold changes
if (!requireNamespace("ashr", quietly = TRUE)) install.packages("ashr")
library(ashr)
res <- results(dds)
res <- lfcShrink(dds, coef = "Disease_carrier_vs_control", res = res, type = "ashr")

# Results to data frame
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$gene_clean <- sub("\\..*", "", res_df$gene)  # keep original ID without version
res_df <- res_df[order(res_df$padj), ]

# Get biotype info and filter for protein-coding genes only
if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt", update = FALSE)
library(biomaRt)
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")

biotype_info <- getBM(
  attributes = c("ensembl_gene_id", "gene_biotype"),
  filters    = "ensembl_gene_id",
  values     = res_df$gene_clean,
  mart       = mart
)

res_df <- merge(res_df, biotype_info, by.x = "gene_clean", by.y = "ensembl_gene_id", all.x = TRUE)
res_df <- subset(res_df, gene_biotype == "protein_coding")

# Significant subsets (thresholds per your previous conventions)
sig_up_npc   <- subset(res_df, padj < 0.05 & log2FoldChange >  1)
sig_down_npc <- subset(res_df, padj < 0.05 & log2FoldChange < -1)

# ---- Pathway enrichment (PATHWAY-ONLY: Reactome, KEGG, WikiPathways) ----
comp_label    <- "NPCs (3.0Mb vs Control)"
top_n_terms   <- 20
fdr_threshold <- 0.05

sig_genes <- unique(c(sig_up_npc$gene, sig_down_npc$gene))
sig_genes_clean <- sub("\\..*", "", sig_genes)

if (length(sig_genes_clean) < 5) {
  stop("Not enough significant genes to run pathway enrichment.")
}

gp_res <- gprofiler2::gost(
  query             = sig_genes_clean,
  organism          = "hsapiens",
  sources           = c("REAC", "KEGG", "WP"),  # pathway-only
  user_threshold    = fdr_threshold,
  correction_method = "fdr",
  significant       = TRUE
)

if (is.null(gp_res) || is.null(gp_res$result) || nrow(gp_res$result) == 0) {
  stop("No significant pathway enrichment terms returned.")
}

gp_df <- as.data.frame(gp_res$result)
# Ensure only pathway sources (safety)
gp_df <- subset(gp_df, source %in% c("REAC","KEGG","WP"))
gp_df$p_value <- as.numeric(gp_df$p_value)
gp_df$source  <- factor(gp_df$source,
                        levels = c("REAC","KEGG","WP"),
                        labels  = c("Reactome","KEGG","WikiPathways"))

# --- Flatten list-columns so write.csv works ---
list_cols <- vapply(gp_df, is.list, logical(1))
if (any(list_cols)) {
  gp_df[list_cols] <- lapply(gp_df[list_cols], function(col)
    vapply(col, function(x) paste(x, collapse = ", "), character(1))
  )
}
# ensure p_value is numeric (already set above, but safe)
gp_df$p_value <- as.numeric(gp_df$p_value)

# Save full pathway enrichment table
write.csv(gp_df, "npc_3.0Mb_vs_Control_pathway_enrichment.csv", row.names = FALSE)

# ---- Pathway dot plot ----
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
ggsave("Pathway_Enrichment_Top20_NPC_3.0Mb_vs_Control.png",
       pathway_plot, width = 11.7, height = 8.3, dpi = 300)
ggsave("Pathway_Enrichment_Top20_NPC_3.0Mb_vs_Control.pdf",
       pathway_plot, width = 11.7, height = 8.3)

# ---- (Optional) Top-10 DEG barplots (kept from your script; remove if not needed) ----
top_up <- head(sig_up_npc[order(-sig_up_npc$log2FoldChange), ], 10)
top_up$label <- ifelse(is.na(top_up$hgnc_symbol) | top_up$hgnc_symbol == "",
                       paste0("novel_", substr(top_up$gene_clean, 10, 15)),
                       top_up$hgnc_symbol)
top_up$label <- factor(top_up$label, levels = rev(top_up$label))

top_down <- head(sig_down_npc[order(sig_down_npc$log2FoldChange), ], 10)
top_down$label <- ifelse(is.na(top_down$hgnc_symbol) | top_down$hgnc_symbol == "",
                         paste0("novel_", substr(top_down$gene_clean, 10, 15)),
                         top_down$hgnc_symbol)
top_down$label <- factor(top_down$label, levels = rev(top_down$label))

pub_theme <- theme_minimal(base_size = 14) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 12, face = "bold"),
    axis.text.x  = element_text(size = 12),
    plot.title   = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 13, face = "bold"),
    plot.margin  = margin(10, 10, 10, 10)
  )

p_up <- ggplot(top_up, aes(x = label, y = log2FoldChange)) +
  geom_bar(stat = "identity", fill = "#d73027", width = 0.6) +
  coord_flip() + pub_theme + labs(title = "Upregulated", y = "Log2 Fold Change")

p_down <- ggplot(top_down, aes(x = label, y = log2FoldChange)) +
  geom_bar(stat = "identity", fill = "#4575b4", width = 0.6) +
  coord_flip() + pub_theme + labs(title = "Downregulated", y = "Log2 Fold Change")

final_plot <- cowplot::plot_grid(p_up, p_down, labels = NULL, ncol = 2, align = "hv")
titled_plot <- cowplot::ggdraw() +
  cowplot::draw_label("Top 10 Up- and Down-regulated Genes – NPCs (3.0Mb vs Control)",
                      fontface = 'bold', size = 16, x = 0.5, y = 0.98, hjust = 0.5) +
  cowplot::draw_plot(final_plot, y = 0, height = 0.95)

ggsave("Top10_NPC_DEG_barplots_3.0Mb_vs_Control_A4.pdf", titled_plot, width = 11.7, height = 8.3)
