# ICA 3: limma differential expression on log2(FPKM + 1)
# Primary contrasts: Drought vs Control separately within IR64 and N22.
# FPKM is not raw counts; do not use DESeq2 or voom on this matrix.

suppressPackageStartupMessages({
  library(limma)
  library(tidyverse)
})

processed_dir <- file.path("data", "processed")
out_dir <- "output"
fig_dir <- "figures"
dir.create(out_dir, showWarnings = FALSE)
dir.create(fig_dir, showWarnings = FALSE)

expr <- readRDS(file.path(processed_dir, "expression_matrix_fpkm.rds"))
meta <- readRDS(file.path(processed_dir, "sample_metadata.rds"))
meta <- meta[match(colnames(expr), meta$sample), ]
stopifnot(all(meta$sample == colnames(expr)))

log_expr <- log2(expr + 1)

run_contrast <- function(genotype_name) {
  keep <- meta$genotype == genotype_name
  m <- droplevels(factor(meta$condition[keep],
                         levels = c("Control", "Drought")))
  x <- log_expr[, keep, drop = FALSE]
  expressed <- rowSums(expr[, keep, drop = FALSE] > 1) >= 2
  x <- x[expressed, , drop = FALSE]
  if (nrow(x) < 2) stop("Insufficient expressed genes in ", genotype_name)

  if (sum(m == "Control") < 2 || sum(m == "Drought") < 2) {
    stop("At least two samples per condition are required for ", genotype_name)
  }

  design <- model.matrix(~ 0 + m)
  colnames(design) <- levels(m)
  fit <- lmFit(x, design)
  contrast <- makeContrasts(Drought_vs_Control = Drought - Control,
                            levels = design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

  tab <- topTable(fit2, coef = "Drought_vs_Control",
                  number = Inf, adjust.method = "BH", sort.by = "P") %>%
    rownames_to_column("gene_id") %>%
    mutate(significance = case_when(
      adj.P.Val < 0.05 & logFC > 0 ~ "Upregulated",
      adj.P.Val < 0.05 & logFC < 0 ~ "Downregulated",
      TRUE ~ "Not significant"
    ))

  write.csv(tab,
            file.path(out_dir, paste0("limma_", genotype_name, ".csv")),
            row.names = FALSE)

  # Volcano plot.
  p <- ggplot(tab, aes(x = logFC, y = -log10(pmax(adj.P.Val, 1e-300)),
                       colour = significance)) +
    geom_point(alpha = 0.65, size = 1) +
    scale_colour_manual(values = c(
      Upregulated = "#D55E00",
      Downregulated = "#0072B2",
      `Not significant` = "grey70"
    )) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_bw() +
    labs(title = paste(genotype_name, "drought vs control"),
         x = "log2 fold change", y = "-log10 adjusted p-value",
         colour = "Result")
  ggsave(file.path(fig_dir, paste0("volcano_", genotype_name, ".png")),
         p, width = 8, height = 6, dpi = 300)

  tab
}

res_ir64 <- run_contrast("IR64")
res_n22 <- run_contrast("N22")

# Combine top significant genes for a descriptive heatmap.
sig_ir64 <- res_ir64 %>% filter(adj.P.Val < 0.05) %>% slice_head(n = 25)
sig_n22 <- res_n22 %>% filter(adj.P.Val < 0.05) %>% slice_head(n = 25)
selected <- unique(c(sig_ir64$gene_id, sig_n22$gene_id))

if (length(selected) >= 2) {
  selected <- selected[seq_len(min(50, length(selected)))]
  mat <- log_expr[selected, , drop = FALSE]
  png(file.path(fig_dir, "heatmap_top_genes.png"),
      width = 1800, height = 1600, res = 200)
  pheatmap::pheatmap(
    mat,
    scale = "row",
    annotation_col = data.frame(
      Genotype = meta$genotype,
      Condition = meta$condition,
      row.names = meta$sample
    ),
    main = "Selected differentially expressed genes"
  )
  dev.off()
} else {
  message("Fewer than two significant genes found; heatmap not generated.")
}

summary_table <- tibble(
  genotype = c("IR64", "N22"),
  genes_tested = c(nrow(res_ir64), nrow(res_n22)),
  significant_FDR_0_05 = c(sum(res_ir64$adj.P.Val < 0.05, na.rm = TRUE),
                           sum(res_n22$adj.P.Val < 0.05, na.rm = TRUE)),
  upregulated = c(sum(res_ir64$significance == "Upregulated"),
                  sum(res_n22$significance == "Upregulated")),
  downregulated = c(sum(res_ir64$significance == "Downregulated"),
                    sum(res_n22$significance == "Downregulated"))
)
write.csv(summary_table, file.path(out_dir, "differential_expression_summary.csv"),
          row.names = FALSE)
print(summary_table)
message("Differential expression results saved in output/.")
