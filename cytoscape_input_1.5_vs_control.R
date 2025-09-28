# ---- CSV for Cytoscape ----
suppressPackageStartupMessages({
  library(DESeq2); library(dplyr); library(ashr); library(biomaRt)
})

# paths
setwd("/Users/kasiadur/Desktop/analysis - part1")
meta   <- read.delim("gene - metadata.txt", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
counts <- read.delim("gene_counts.txt", header = TRUE, row.names = 1, check.names = FALSE)

# subset NPCs and patients
npc_meta <- subset(meta, Celltype == "p" & Patient %in% c("d8","d9","c11","c35","c1"))

# collapse replicates
npc_meta$Group <- ifelse(
  grepl("^DS", npc_meta$Row.names),
  gsub("(_[anp]\\d+)_\\d+_S\\d+", "\\1", npc_meta$Row.names),
  gsub("(_p\\d+)_\\d+$", "\\1", npc_meta$Row.names)
)
npc_counts <- counts[, npc_meta$Row.names, drop = FALSE]
collapsed_counts <- sapply(unique(npc_meta$Group), function(g){
  reps <- npc_meta$Row.names[npc_meta$Group == g]
  rowSums(npc_counts[, reps, drop = FALSE])
})
collapsed_counts <- as.data.frame(collapsed_counts)
colnames(collapsed_counts) <- unique(npc_meta$Group)

# filter low counts
keep <- rowSums(collapsed_counts >= 10) >= 2
collapsed_counts <- collapsed_counts[keep, , drop = FALSE]

# collapsed metadata
npc_meta_collapsed <- npc_meta[!duplicated(npc_meta$Group), ]
rownames(npc_meta_collapsed) <- npc_meta_collapsed$Group
npc_meta_collapsed <- npc_meta_collapsed[, c("Disease","Celltype","Patient")]
npc_meta_collapsed$Disease <- factor(npc_meta_collapsed$Disease, levels = c("control","carrier"))
npc_meta_collapsed$Sex <- dplyr::case_when(
  npc_meta_collapsed$Patient == "d8"  ~ "male",
  npc_meta_collapsed$Patient == "d9"  ~ "female",
  npc_meta_collapsed$Patient == "c11" ~ "male",
  npc_meta_collapsed$Patient == "c35" ~ "female",
  TRUE ~ "unknown"
)
npc_meta_collapsed$Sex <- factor(npc_meta_collapsed$Sex)

# DESeq2
dds <- DESeqDataSetFromMatrix(countData = collapsed_counts,
                              colData = npc_meta_collapsed,
                              design = ~ Sex + Disease)
dds <- DESeq(dds)
res <- results(dds)
res <- lfcShrink(dds, coef = "Disease_carrier_vs_control", res = res, type = "ashr")
res_df <- as.data.frame(res)
res_df$input_id <- rownames(res_df)

suppressPackageStartupMessages({
  library(dplyr)
  library(biomaRt)
})

# res_df = your DESeq2 results (rownames = Entrez IDs)
res_df$Entrez <- as.character(rownames(res_df))   # force character

# connect to Ensembl
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# map Entrez → Ensembl
id_map <- getBM(
  attributes = c("entrezgene_id", "ensembl_gene_id"),
  filters    = "entrezgene_id",
  values     = res_df$Entrez,
  mart       = mart
)

id_map$entrezgene_id <- as.character(id_map$entrezgene_id)  # also character

# merge mapping
res_annot <- res_df %>%
  left_join(id_map, by = c("Entrez" = "entrezgene_id"))

# build Cytoscape table
cyto_tbl <- res_annot %>%
  transmute(
    Ensembl = ensembl_gene_id,
    log2FC  = log2FoldChange,
    pvalue,
    padj
  ) %>%
  filter(!is.na(Ensembl)) %>%
  distinct(Ensembl, .keep_all = TRUE)

write.csv(cyto_tbl, "cytoscape_input_1.5_vs_control.csv", row.names = FALSE)
