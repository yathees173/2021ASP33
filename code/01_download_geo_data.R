# Optional: download supplementary data from NCBI GEO if needed.
# The supplied original dataset files are already included in data/raw/.
# Run this only if you need to restore/re-download the source files.

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery")
}
GEOquery::getGEOSuppFiles(
  GEO = "GSE183894",
  baseDir = "data/geo_download",
  makeDirectory = TRUE,
  fetch_files = TRUE
)

