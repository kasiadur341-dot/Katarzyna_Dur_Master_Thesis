# ===========================================
# SETUP
# ===========================================
setwd("/Users/kasiadur/Desktop/analysis - part1")  # change if needed

suppressPackageStartupMessages({
  library(pheatmap)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(grid)
  library(biomaRt)
})

# ===========================================
# STEP 1: Full 22q11.2 gene list from Ensembl
# ===========================================
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# GRCh38 coordinates of 22q11.2 region (approx 18–22 Mb)
region_chr   <- "22"
region_start <- 18924718L
region_end   <- 21500000L

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

genes_22q <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype",
                 "start_position", "end_position"),
  filters    = "chromosomal_region",
  values     = sprintf("%s:%d:%d", region_chr, region_start, region_end),
  mart       = mart
)

# keep genes with HGNC symbols and prep for matching
genes_22q <- genes_22q[genes_22q$hgnc_symbol != "", ]
genes_22q$hgnc_symbol_upper <- toupper(genes_22q$hgnc_symbol)
pathway_genes <- unique(genes_22q$hgnc_symbol_upper)

# uppercase for matching
genes_22q$hgnc_symbol_upper <- toupper(genes_22q$hgnc_symbol)

# final vector for downstream filtering
pathway_genes <- unique(genes_22q$hgnc_symbol_upper)

cat("Number of genes in 22q11.2 list:", length(pathway_genes), "\n")

# ===========================================
# STEP 2: Load DE results
# (expects columns: ensembl_gene_id, log2FC_1.5Mb, log2FC_3.0Mb, padj_1.5Mb, padj_3.0Mb)
# ===========================================
de <- read.csv("merged_for_pathvisio.csv", stringsAsFactors = FALSE)

# ===========================================
# STEP 3: Map Ensembl IDs to HGNC
# ===========================================
mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = de$ensembl_gene_id, mart = mart
)
de <- merge(de, mapping, by = "ensembl_gene_id", all.x = TRUE)
de$hgnc_symbol_upper <- toupper(de$hgnc_symbol)

# ===========================================
# STEP 4: Filter for pathway genes & save
# ===========================================
de_pathway <- de %>% dplyr::filter(hgnc_symbol_upper %in% pathway_genes)
write.csv(de_pathway, "expression_for_22q11_pathway.csv", row.names = FALSE)

# ===========================================
# STEP 5: Summarize significance (alpha = 0.05)
# ===========================================
alpha <- 0.05
sig_15  <- subset(de_pathway, padj_1.5Mb < alpha)
sig_30  <- subset(de_pathway, padj_3.0Mb < alpha)
up_30   <- subset(sig_30,  `log2FC_3.0Mb` > 1)
down_30 <- subset(sig_30,  `log2FC_3.0Mb` < -1)

write.csv(up_30,   "upregulated_pathway_genes_3Mb.csv",   row.names = FALSE)
write.csv(down_30, "downregulated_pathway_genes_3Mb.csv", row.names = FALSE)

# ===========================================
# STEP 6: Build matrix for Heatmap (auto-detect all log2FC_* columns)
# ===========================================
log2_cols <- grep("^log2FC", names(de_pathway), value = TRUE)

fc_matrix <- de_pathway |>
  dplyr::select(hgnc_symbol_upper, dplyr::all_of(log2_cols)) |>
  tibble::column_to_rownames("hgnc_symbol_upper") |>
  as.matrix()

# Heatmap with fixed scale (-1 to 1)
pheatmap(
  fc_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "none",
  color = colorRampPalette(c("blue", "white", "red"))(201),
  breaks = seq(-1, 1, length.out = 202),
  main = "Expression of 22q11.2 Genes (log2FC)",
  fontsize_row = 10,
  fontsize_col = 12,
  legend = TRUE,
  angle_col = 45,
  filename = "heatmap_22q11_fixedscale.pdf",
  width = 8, height = 10
)

# ===========================================
# STEP 7: Bar plot of significant Up/Down
# ===========================================
counts <- list(
  `1.5Mb` = data.frame(
    Condition = "1.5Mb",
    Up   = sum(sig_15$`log2FC_1.5Mb`  > 1, na.rm = TRUE),
    Down = sum(sig_15$`log2FC_1.5Mb`  < -1, na.rm = TRUE)
  ),
  `3.0Mb` = data.frame(
    Condition = "3.0Mb",
    Up   = sum(sig_30$`log2FC_3.0Mb`  > 1, na.rm = TRUE),
    Down = sum(sig_30$`log2FC_3.0Mb`  < -1, na.rm = TRUE)
  )
) |> dplyr::bind_rows()

counts_long <- counts |>
  tidyr::pivot_longer(cols = c("Up","Down"), names_to = "Direction", values_to = "Count") |>
  dplyr::mutate(Direction = factor(Direction, levels = c("Up","Down")))

y_max <- max(counts_long$Count, na.rm = TRUE)
y_limit <- if (is.finite(y_max) && y_max > 0) y_max * 1.15 else 1

p <- ggplot(counts_long, aes(x = Condition, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
  geom_text(aes(label = Count),
            position = position_dodge(width = 0.7),
            vjust = -0.6, size = 6) +
  labs(
    title = "Up/Down-Regulated Genes in 22q11.2 Pathway",
    y = "Number of Genes", x = "Condition"
  ) +
  theme_minimal(base_size = 16) +
  scale_fill_manual(name = "Direction", values = c("Up" = "red", "Down" = "blue")) +
  scale_y_continuous(limits = c(0, y_limit), expand = expansion(mult = c(0.02, 0.08))) +
  coord_cartesian(clip = "off") +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(t = 10, b = 10)),
    plot.margin = margin(t = 24, r = 20, b = 20, l = 24)
  )

ggsave("barplot_22q11_fixed.pdf", p, width = 8, height = 6)

cat("Saved:\n - heatmap_22q11_fixedscale.pdf\n - barplot_22q11_fixed.pdf\n")
