# Run from project root in RStudio or command line.
if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Install rmarkdown first with code/00_install_packages.R")
}

# Ensure pandoc is available when running outside of RStudio
if (!rmarkdown::pandoc_available()) {
  potential_paths <- c(
    "K:/Installed Application/RStudio/resources/app/bin/quarto/bin/tools",
    "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools",
    "C:/Program Files/RStudio/bin/pandoc"
  )
  for (p in potential_paths) {
    if (dir.exists(p)) {
      Sys.setenv(RSTUDIO_PANDOC = p)
      break
    }
  }
}

dir.create("output", showWarnings = FALSE)
rmarkdown::render(
  input = "code/Project_Report.Rmd",
  output_format = "html_document",
  output_file = "Project_Report.html",
  output_dir = normalizePath("output"),
  knit_root_dir = normalizePath("code"),
  quiet = FALSE
)
message("Rendered HTML: output/Project_Report.html")
