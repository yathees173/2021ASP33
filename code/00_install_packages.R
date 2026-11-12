# Install packages
# Run this script once in RStudio.

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c("tidyverse", "ggplot2", "pheatmap", "RColorBrewer",
                   "knitr", "rmarkdown", "statmod")
missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_cran)) install.packages(missing_cran)

bioc_packages <- c("GEOquery", "limma")
missing_bioc <- bioc_packages[
  !vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

message("Package installation check completed.")
