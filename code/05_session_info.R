# Save software and package versions for reproducibility.
dir.create("output", showWarnings = FALSE)
writeLines(capture.output(sessionInfo()),
           con = file.path("output", "sessionInfo.txt"))
message("R session information saved to output/sessionInfo.txt")
