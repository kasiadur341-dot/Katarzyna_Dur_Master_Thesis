# ================================
# Mislabel check: NPC 1.5Mb vs 3.0Mb
# ================================

# --- Setup ---
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(matrixStats)
  library(limma)
  library(caret)
  library(glmnet)
  library(MVN)
  library(ggrepel)
  library(ggfortify)
})

# ---------- I/O ----------
outdir <- "mislabel_check"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

meta_path   <- "gene - metadata.txt"
counts_path <- "gene_counts.txt"
interval_path <- "22q11_interval_genes.tsv"   # optional

# ---------- Load ----------
meta   <- read.delim(meta_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
counts <- read.delim(counts_path, row.names = 1, check.names = FALSE)

# ---------- Subset to NPCs & carriers ----------
# Keep NPCs ("p") for patients used in your 1.5 vs 3.0 analysis
npc_meta <- meta %>%
  filter(Celltype == "p", Patient %in% c("d8","d9","d36","d37","d10")) %>%
  mutate(Disease = ifelse(Patient %in% c("d8","d9"), "1.5Mb", "3.0Mb"))

# align counts with metadata
stopifnot(all(npc_meta$Row.names %in% colnames(counts)))
npc_counts <- counts[, npc_meta$Row.names]
npc_meta   <- npc_meta %>% arrange(match(Row.names, colnames(npc_counts)))
stopifnot(identical(npc_meta$Row.names, colnames(npc_counts)))

# ---------- Normalize (DESeq2 VST) ----------
dds <- DESeqDataSetFromMatrix(countData = npc_counts,
                              colData   = npc_meta[, c("Row.names","Patient","Disease")],
                              design    = ~ 1)
dds <- estimateSizeFactors(dds)
v   <- vst(dds, blind = TRUE)
X   <- assay(v)                    # genes x samples
Xl2 <- log1p(counts(dds, normalized = TRUE))  # for PCA alt.

# ---------- 1) Correlations & clustering ----------
cors <- cor(X, method = "spearman")
pheatmap(cors,
         annotation_col = npc_meta[,c("Patient","Disease")] %>% `rownames<-`(npc_meta$Row.names),
         main = "Sample–sample Spearman correlations (VST)",
         filename = file.path(outdir, "01_correlation_heatmap.pdf"),
         width = 6, height = 5)

hc <- hclust(as.dist(1 - cors), method = "average")
pdf(file.path(outdir, "01_dendrogram.pdf"), width = 6, height = 4)
plot(hc, main = "Hierarchical clustering (1 - Spearman)", xlab = "", sub = "")
dev.off()

# ---------- 2) Distances to centroids ----------
# ---- helper: group centroids (rows = groups, cols = genes) ----
centroid_by <- function(M, groups) {
  groups <- as.factor(groups)
  grp_levels <- levels(groups)
  cent <- lapply(grp_levels, function(g) {
    idx <- which(groups == g)
    rowMeans(M[, idx, drop = FALSE])
  })
  cent <- do.call(rbind, cent)
  rownames(cent) <- grp_levels
  cent
}

# Build disease centroids
cent_disease <- centroid_by(X, npc_meta$Disease)   # rows: 1.5Mb, 3.0Mb ; cols: genes

# Distances of each sample to each disease centroid (Euclidean in VST space)
dist_to_grp <- t(apply(X, 2, function(s) {
  # subtract sample from each centroid (row) gene-wise, then Euclidean norm
  diffs <- sweep(cent_disease, 2, s, FUN = "-")   # same dims as cent_disease
  sqrt(rowSums(diffs^2))
}))

# Name the two columns with the *row* names of the centroid matrix (the groups)
colnames(dist_to_grp) <- paste0("dist_", rownames(cent_disease))
dist_to_grp <- as.data.frame(dist_to_grp)

# Add labels and derived fields
dist_to_grp$sample    <- colnames(X)
dist_to_grp$true      <- npc_meta$Disease
# figure out which centroid is closer
grp_names <- rownames(cent_disease)
dist_to_grp$closer_to <- ifelse(dist_to_grp[[paste0("dist_", grp_names[1])]] <
                                  dist_to_grp[[paste0("dist_", grp_names[2])]],
                                grp_names[1], grp_names[2])
# distance margin between the two groups
d1 <- dist_to_grp[[paste0("dist_", grp_names[1])]]
d2 <- dist_to_grp[[paste0("dist_", grp_names[2])]]
dist_to_grp$margin <- abs(d1 - d2)

# (optional) sanity check
print(head(dist_to_grp, 3))

# ---------- 3) Dosage signature over 22q11.2 interval (optional) ----------
if (file.exists(interval_path)) {
  interval_genes <- read.delim(interval_path, header = TRUE, stringsAsFactors = FALSE)
  gcol <- intersect(colnames(interval_genes), c("gene_id","gene","symbol","Gene","ENSEMBL","ENSEMBL_ID"))
  if (length(gcol) == 0) stop("Interval file found but no gene ID column recognized.")
  id_col <- gcol[1]
  gset <- intersect(rownames(X), unique(interval_genes[[id_col]]))
  message("Interval genes detected: ", length(gset), " found in matrix.")
} else {
  gset <- character(0)
  message("Interval file not found; will fall back to top variable genes for classifier.")
}

dosage_tbl <- NULL
if (length(gset) >= 10) {
  dosage_score <- colMeans(X[gset, , drop = FALSE])
  dosage_tbl <- data.frame(sample = names(dosage_score),
                           dosage_score = as.numeric(dosage_score)) %>%
    left_join(npc_meta[,c("Row.names","Patient","Disease")], by = c("sample" = "Row.names"))
  write.csv(dosage_tbl, file.path(outdir, "03_interval_dosage_scores.csv"), row.names = FALSE)
  
  p_dos <- ggplot(dosage_tbl, aes(x = Disease, y = dosage_score, color = Disease)) +
    geom_jitter(width = 0.1, size = 2) +
    geom_boxplot(outlier.shape = NA, alpha = 0.2) +
    theme_minimal() + labs(title = "22q11.2 interval dosage score (VST mean)")
  ggsave(file.path(outdir, "03_interval_dosage_scores.pdf"), p_dos, width = 5, height = 4)
}

# ---------- 4) LOOCV classifier for Disease ----------

set.seed(7)

# Make levels safe for caret
Y <- factor(npc_meta$Disease)
levels(Y) <- c("Mb15","Mb30")   # instead of 1.5Mb / 3.0Mb

# Features: top variable genes
vars <- matrixStats::rowVars(X)
top  <- names(sort(vars, decreasing = TRUE))[seq_len(min(500, nrow(X)))]
Xsub <- t(X[top, ])

# Train with LOOCV
ctrl <- trainControl(method = "LOOCV", classProbs = TRUE, summaryFunction = twoClassSummary)
fit  <- train(x = Xsub, y = Y, method = "glmnet",
              trControl = ctrl, metric = "ROC", family = "binomial")

# Predictions
pred <- predict(fit, Xsub)
prob <- predict(fit, Xsub, type = "prob")

# Build output table, map back to original labels
map_back <- c(Mb15 = "1.5Mb", Mb30 = "3.0Mb")
pred_tbl <- data.frame(
  sample   = rownames(Xsub),
  true     = npc_meta$Disease,
  pred     = map_back[as.character(pred)],
  prob_1.5 = prob[,"Mb15"],
  prob_3.0 = prob[,"Mb30"]
)

write.csv(pred_tbl, file.path(outdir, "04_classifier_predictions.csv"), row.names = FALSE)

# ---------- 5) Within-donor outlier check (PC-space Mahalanobis with ridge) ----------
# Build PCA once on all samples (VST matrix X: genes x samples)
pca_all <- prcomp(t(X), center = TRUE, scale. = TRUE)
K <- min(10, ncol(pca_all$x))   # number of PCs to use (tune if you like)
PC <- pca_all$x[, 1:K, drop = FALSE]   # samples x K
rownames(PC) <- colnames(X)            # sample names

ridge_mahal <- function(M, lam = 1e-3) {
  # M: n x K (rows = samples from one donor)
  mu  <- colMeans(M)
  S   <- cov(M)
  K   <- ncol(M)
  Srg <- S + diag(lam, K)       # ridge regularization
  as.numeric(mahalanobis(M, center = mu, cov = Srg))
}

outlier_rows <- list()
for (don in unique(npc_meta$Patient)) {
  idx <- which(npc_meta$Patient == don)
  if (length(idx) < 2) next                # need at least 2 reps
  S <- PC[npc_meta$Row.names[idx], , drop = FALSE]  # donor's PC scores
  m <- ridge_mahal(S, lam = 1e-3)
  z <- as.numeric(scale(m))
  outlier_rows[[don]] <- data.frame(
    sample = rownames(S),
    donor  = don,
    mahal  = m,
    zscore = z
  )
}
outlier_tbl <- dplyr::bind_rows(outlier_rows)
readr::write_csv(outlier_tbl, file.path(outdir, "05_within_donor_mahalanobis_PC.csv"))

# flag suspicious ones
sus <- subset(outlier_tbl, abs(zscore) > 3)
if (nrow(sus) == 0) {
  cat("\nWithin-donor outliers (PC-space mahalanobis): none with |z| > 3.\n")
} else {
  cat("\nWithin-donor outliers (PC-space mahalanobis, |z| > 3):\n")
  print(sus)
}

# ---------- Compact on-screen summary ----------
cat("\n=== MISLABEL CHECK SUMMARY ===\n")
# 1) Which samples are closer to the 'other' Disease centroid?
flip <- dist_to_grp %>%
  mutate(is_flip = closer_to != true) %>%
  arrange(desc(margin))
print(flip[, c("sample","true","closer_to","margin")] )

# 2) Classifier disagreements
disagree <- pred_tbl %>% filter(as.character(true) != as.character(pred))
if (nrow(disagree) == 0) {
  cat("\nClassifier: all samples predicted as their labeled group (LOOCV).\n")
} else {
  cat("\nClassifier disagreements (possible mislabels):\n")
  print(disagree)
}

# 3) Within-donor outliers (z > 3)
if (!is.null(outlier_tbl)) {
  sus <- outlier_tbl %>% filter(abs(zscore) > 3)
  if (nrow(sus) == 0) {
    cat("\nWithin-donor outliers: none with |z| > 3.\n")
  } else {
    cat("\nWithin-donor outliers (|z| > 3):\n")
    print(sus)
  }
}

cat(paste0("\nResults written to: ", normalizePath(outdir), "\n"))
