# Runs all stages in order, from the raw NHANES files to the results, figures, and Table 1.
# Usage, from the repository root: Rscript run_all.R
# Optional: STUDY4_BUILD sets the directory for intermediate files (default: _build/).
# Runtime is about 5 to 10 minutes; the jackknife, bootstrap, and imputation steps take most of it.

t0 <- Sys.time()
rscript <- file.path(R.home("bin"), "Rscript")
run <- function(cmd, args) {
  message("\n=== ", paste(c(basename(cmd), args), collapse = " "))
  st <- system2(cmd, args)
  if (!identical(as.integer(st), 0L)) stop("Failed: ", paste(args, collapse = " "))
}
for (s in c(
  "01_data_derivation.R", "02_cohort.R", "03_primary_analysis.R", "04_sensitivity_analysis.R",
  "05_standardization.R", "05b_standardized_curves.R", "05c_function_and_prediction.R",
  "05d_guideline_categories.R", "06_figures.R", "07_tables.R", "08_reproducibility_check.R"
))
  run(rscript, c("--vanilla", file.path("code", s)))
writeLines(capture.output(sessionInfo()), file.path("output", "sessionInfo.txt"))
message(sprintf(
  "\nAll stages completed in %.1f minutes",
  as.numeric(difftime(Sys.time(), t0, units = "mins"))
))
