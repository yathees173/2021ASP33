# Run this file from the project root after installing dependencies.
# The 12 original GSE183894 FPKM files are already in data/raw/.

source("code/02_prepare_expression.R")
source("code/03_exploratory_analysis.R")
source("code/04_differential_expression_limma.R")
source("code/05_session_info.R")
source("code/06_render_report.R")
message("Pipeline finished. Inspect output/ and figures/ before submission.")

