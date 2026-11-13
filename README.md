# Differential Gene Expression Analysis of Crop Drought Stress Response

**Coursework:** Bioinformatics  
**Author:** Yathees Thiraviyarajah  
**Student ID:** 2021/ASP/33  
**Organism:** *Oryza sativa* (rice)  
**GEO accession:** [GSE183894](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183894)  
**Analysis environment:** R 4.3.2; Windows 11  
**Main statistical package:** `limma`

---

## Abstract

Drought stress can alter crop growth and gene expression. This project analyses public rice (*Oryza sativa*) data from GEO series GSE183894, comparing drought-simulation and control samples from IR64 and N22. Twelve FPKM files were combined into a gene-by-sample matrix, followed by quality assessment, PCA, correlation analysis and differential expression modelling. Since the data are processed FPKM values rather than raw counts, the analysis uses log2(FPKM + 1) and `limma`, not DESeq2 or `voom`. The union matrix contains 39,718 genes. Separate comparisons tested 23,019 genes in IR64 and 23,316 in N22. At adjusted p-value < 0.05, 6,157 genes were significant in IR64 and 7,520 in N22. These results are candidate expression changes, not proof that any gene causes drought tolerance.

## 1. Introduction

Water limitation affects plant growth and can alter genes involved in signalling, metabolism and stress responses. Rice is an important crop for studying these processes. The selected GEO study describes IR64 as drought-sensitive and N22 as drought-tolerant. Comparing drought-simulation and control samples within each cultivar allows their expression responses to be examined.

**Research question:** Which genes show evidence of differential expression between drought-treated and control rice seedlings in IR64 and N22?

**Aim:** To analyse public rice gene-expression data and identify candidate genes whose expression differs between drought-simulation and control conditions in the two cultivars.

**Objectives:**

1. Obtain and document the public GSE183894 expression dataset.
2. Inspect the supplied sample files and construct a consistent gene-expression matrix.
3. Check sample distributions, correlations and major expression patterns using exploratory visualisations.
4. Apply `limma` to log-transformed FPKM values for separate drought-versus-control comparisons in IR64 and N22.
5. Summarise and visualise the differential-expression results and discuss their limitations.

## 2. Dataset and experimental design

The data are from NCBI GEO series [GSE183894](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183894), profiling IR64 and N22 rice cultivars. GEO sample descriptions indicate that 21-day-old seedlings were exposed for six hours to reverse-osmosis water (control) or 20% PEG 6000 (drought-simulation treatment). PEG imposes osmotic stress but does not reproduce every feature of field drought.

The supplied project includes 12 compressed tab-delimited files (`.txt.gz`). The experimental design is balanced, with three biological replicates for each cultivar-condition combination:

| Cultivar | Control | Drought simulation | Total |
|---|---:|---:|---:|
| IR64 | 3 | 3 | 6 |
| N22 | 3 | 3 | 6 |
| **Total** | **6** | **6** | **12** |

The files are stored in `data/raw/`. Sample identifiers and group labels are parsed from filenames and recorded in `data/processed/sample_metadata.csv`. The preparation summary records 12 original files and 12 retained samples. The combined matrix contains 39,718 unique gene identifiers in the union across the sample files. The preprocessing summary reports 13,988 genes with at least one zero after the union and no negative expression values.

### Data type and analytical implications

The supplied values are FPKM, normalised for transcript length and sequencing depth, not raw read counts. The workflow therefore uses `log2(FPKM + 1)` and `limma`; DESeq2 and `voom` are not applied. FPKM has limitations for between-sample comparison, so these results analyse the supplied processed matrix and do not replace a raw-count analysis.

## 3. Methods

### 3.1 Software and packages

The project was run using R 4.3.2 on Windows 11. The main packages include `tidyverse` for data import and manipulation, `ggplot2` for visualisation, `pheatmap` for heatmaps, and Bioconductor's `limma` package for differential expression. `GEOquery` is included for optional retrieval of supplementary GEO files, although the original 12 files are already supplied in `data/raw/` and do not need to be downloaded again.

### 3.2 Preprocessing

The preprocessing script reads the compressed `.txt.gz` files directly. It checks that exactly 12 FPKM files are available and parses each filename to determine the GEO sample accession, cultivar, condition and replicate. Gene expression values are combined by gene identifier to form a matrix with genes as rows and samples as columns. Gene short names are stored separately in `data/processed/gene_annotation.csv`. The resulting matrix and metadata are saved in both CSV and RDS formats so that they can be inspected and reused by later scripts.

The processing summary indicates 12 input files and 12 samples. The union matrix contains 39,718 genes. A zero in the union matrix indicates that the corresponding gene has a zero expression value in at least one sample; it does not mean that the gene is absent from every sample. No negative expression values were reported. The project preserves the original compressed files and writes derived data to `data/processed/`.

### 3.3 Exploratory analysis and quality assessment

Expression values were transformed using `log2(FPKM + 1)`. Boxplots, density plots, a Pearson correlation heatmap and PCA based on up to 500 variable genes were generated. In the supplied PCA, PC1 explains 74.1% and PC2 20.7% of variance. Cultivars separate mainly along PC1, while treatment patterns are visible along PC2. Replicates generally cluster by cultivar and condition in the correlation heatmap. These are descriptive checks, not grounds for removing samples without independent evidence.