# ICA 3: exploratory analysis and quality-control visualisations

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
})

processed_dir <- file.path("data", "processed")
fig_dir <- "figures"
dir.create(fig_dir, showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)

expr_file <- file.path(processed_dir, "expression_matrix_fpkm.rds")
meta_file <- file.path(processed_dir, "sample_metadata.rds")
if (!file.exists(expr_file) || !file.exists(meta_file)) {
  stop("Run code/02_prepare_expression.R first.")
}

expr <- readRDS(expr_file)
meta <- readRDS(meta_file)
meta <- meta[match(colnames(expr), meta$sample), ]
stopifnot(all(meta$sample == colnames(expr)))

log_expr <- log2(expr + 1)

# Long-format expression distribution.
long_expr <- as.data.frame(log_expr) %>%
  rownames_to_column("gene_id") %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "log2_fpkm") %>%
  left_join(meta, by = "sample")

p_box <- ggplot(long_expr, aes(x = sample, y = log2_fpkm,
                               fill = interaction(genotype, condition))) +
  geom_boxplot(outlier.size = 0.15) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "Expression distribution across samples",
       x = "Sample", y = "log2(FPKM + 1)", fill = "Group")
ggsave(file.path(fig_dir, "expression_boxplot.png"), p_box,
       width = 10, height = 6, dpi = 300)

p_density <- ggplot(long_expr, aes(x = log2_fpkm, colour = sample)) +
  geom_density(linewidth = 0.5, na.rm = TRUE) +
  theme_bw() +
  labs(title = "Expression density by sample",
       x = "log2(FPKM + 1)", y = "Density")
ggsave(file.path(fig_dir, "expression_density.png"), p_density,
       width = 10, height = 6, dpi = 300)

# Sample correlation.
cor_mat <- cor(log_expr, use = "pairwise.complete.obs", method = "pearson")
png(file.path(fig_dir, "sample_correlation.png"),
    width = 1800, height = 1600, res = 200)
pheatmap(cor_mat, annotation_col = data.frame(
  Genotype = meta$genotype,
  Condition = meta$condition,
  row.names = meta$sample
), main = "Sample-to-sample Pearson correlation")
dev.off()

# PCA on genes with highest variance.
gene_var <- apply(log_expr, 1, var, na.rm = TRUE)
n_top <- min(500, sum(is.finite(gene_var) & gene_var > 0))
if (n_top < 3) stop("Too few variable genes for PCA")
top_genes <- names(sort(gene_var[is.finite(gene_var) & gene_var > 0], decreasing = TRUE))[seq_len(n_top)]
pca <- prcomp(t(log_expr[top_genes, , drop = FALSE]), scale. = TRUE)
pca_df <- as.data.frame(pca$x[, 1:2]) %>%
  rownames_to_column("sample") %>%
  left_join(meta, by = "sample")
variance_pct <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = genotype,
                            shape = condition, label = sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.8, size = 2.8, show.legend = FALSE) +
  theme_bw() +
  labs(title = "PCA of rice samples",
       x = paste0("PC1 (", round(variance_pct[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_pct[2], 1), "%)"))
ggsave(file.path(fig_dir, "pca_samples.png"), p_pca,
       width = 9, height = 6, dpi = 300)

write.csv(pca_df, file.path("output", "pca_coordinates.csv"),
          row.names = FALSE)
message("Exploratory plots saved in figures/. Review them before interpreting.")
