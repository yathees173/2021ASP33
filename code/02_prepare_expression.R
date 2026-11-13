#combine the 12 supplied original GEO FPKM files.

suppressPackageStartupMessages(library(tidyverse))

raw_dir <- file.path("data", "raw")
processed_dir <- file.path("data", "processed")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(raw_dir, pattern = "FPKM\\.txt\\.gz$", full.names = TRUE,
                    recursive = TRUE)
if (length(files) != 12) {
  stop("Expected 12 original FPKM files; found ", length(files),
       ". Check data/raw/ and do not mix unrelated files.")
}

parse_metadata <- function(path) {
  nm <- basename(path)
  m <- stringr::str_match(
    nm, "^(GSM[0-9]+)_(IR64|N22)_(CT|DS)_REP([0-9]+)_FPKM\\.txt\\.gz$"
  )
  if (any(is.na(m[1, ]))) stop("Unexpected filename: ", nm)
  tibble(
    sample = paste(m[2], m[3], m[4], paste0("REP", m[5]), sep = "_"),
    gsm = m[2],
    genotype = m[3],
    condition = recode(m[4], CT = "Control", DS = "Drought"),
    replicate = as.integer(m[5]),
    file = path
  )
}
metadata <- purrr::map_dfr(files, parse_metadata) %>%
  arrange(genotype, condition, replicate)

# Confirm a balanced 2 genotype x 2 condition x 3 replicate design.
group_counts <- metadata %>% count(genotype, condition)
print(group_counts)
if (nrow(group_counts) != 4 || any(group_counts$n != 3)) {
  stop("Expected 3 samples in each of the four genotype/condition groups.")
}
if (anyDuplicated(metadata$gsm)) stop("Duplicate GEO sample accessions detected.")

read_one <- function(path) {
  d <- read.delim(gzfile(path), header = TRUE, sep = "\t",
                  skip = 0, blank.lines.skip = TRUE,
                  check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("gene_id", "gene_short_name")
  if (!all(required %in% names(d))) {
    stop("Expected gene_id and gene_short_name columns in ", basename(path))
  }
  expr_cols <- setdiff(names(d), required)
  if (length(expr_cols) != 1) {
    stop("Expected exactly one FPKM column in ", basename(path))
  }
  tibble(
    gene_id = as.character(d$gene_id),
    gene_short_name = as.character(d$gene_short_name),
    fpkm = suppressWarnings(as.numeric(d[[expr_cols]]))
  ) %>%
    filter(!is.na(gene_id), gene_id != "", !is.na(fpkm)) %>%
    distinct(gene_id, .keep_all = TRUE)
}

sample_data <- purrr::map(metadata$file, read_one)
names(sample_data) <- metadata$sample

# Preserve both the transcript-assembly identifier and gene annotation.
annotation <- bind_rows(sample_data) %>%
  select(gene_id, gene_short_name) %>%
  distinct(gene_id, .keep_all = TRUE)

# Join one sample at a time by gene identifier (not by row order).
expr_wide <- tibble(gene_id = sort(unique(unlist(
  lapply(sample_data, function(d) d$gene_id)))))
for (sample_name in names(sample_data)) {
  one <- sample_data[[sample_name]] %>% select(gene_id, fpkm)
  names(one)[2] <- sample_name
  expr_wide <- left_join(expr_wide, one, by = "gene_id")
}

expr_mat <- expr_wide %>%
  column_to_rownames("gene_id") %>%
  as.matrix()
storage.mode(expr_mat) <- "numeric"
expr_mat[is.na(expr_mat)] <- 0

# Reorder columns to match metadata exactly.
expr_mat <- expr_mat[, metadata$sample, drop = FALSE]

write.csv(metadata, file.path(processed_dir, "sample_metadata.csv"),
          row.names = FALSE)
write.csv(annotation, file.path(processed_dir, "gene_annotation.csv"),
          row.names = FALSE)
write.csv(data.frame(gene_id = rownames(expr_mat), expr_mat,
                     check.names = FALSE),
          file.path(processed_dir, "expression_matrix_fpkm.csv"),
          row.names = FALSE)
saveRDS(expr_mat, file.path(processed_dir, "expression_matrix_fpkm.rds"))
saveRDS(metadata, file.path(processed_dir, "sample_metadata.rds"))

summary <- tibble(
  metric = c("original_files", "samples", "genes_union",
             "genes_with_any_zero_after_union", "negative_values"),
  value = c(length(files), ncol(expr_mat), nrow(expr_mat),
            sum(rowSums(expr_mat == 0) > 0), sum(expr_mat < 0))
)
write.csv(summary, file.path(processed_dir, "preprocessing_summary.csv"),
          row.names = FALSE)
print(summary)
message("Prepared matrix: ", nrow(expr_mat), " genes x ",
        ncol(expr_mat), " samples.")

