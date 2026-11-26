# Stage 8. Checks: (1) collects all analysis outputs in output/csv/analysis_outputs/ and compares each
# file with the reported results in results/analysis_outputs/; (2) checks that the Table 1 values equal
# the analysis output; (3) runs code/independent_check.R. Reports are written to output/qc/.

source("code/00_setup.R")
suppressPackageStartupMessages(library(dplyr))
dir.create(file.path(P$csv, "analysis_outputs"), showWarnings = FALSE)
outs <- c(
  list.files(P$out_orig, "\\.csv$", full.names = TRUE),
  list.files(P$out_add, "\\.csv$", full.names = TRUE)
)
stopifnot(!anyDuplicated(basename(outs)))
invisible(file.copy(outs, file.path(P$csv, "analysis_outputs", basename(outs)), overwrite = TRUE))
message(length(outs), " analysis outputs from stages 01-05 copied to output/csv/analysis_outputs")
# Compare every output with the reported results in results/analysis_outputs/.
ref_dir <- file.path(PKG, "results", "analysis_outputs")
if (dir.exists(ref_dir)) {
  ref_files <- list.files(ref_dir, "\\.csv$")
  sha <- function(f) if (file.exists(f)) digest::digest(
    file = f,
    algo = "sha256"
  ) else NA_character_
  cmp <- data.frame(
    file = ref_files,
    reported_sha256 = vapply(file.path(ref_dir, ref_files), sha, ""),
    reproduced_sha256 = vapply(file.path(P$csv, "analysis_outputs", ref_files), sha, ""),
    row.names = NULL
  )
  participant_file <- "guideline_analysis_variables.csv"
  cmp <- rbind(cmp, data.frame(file = participant_file,
    reported_sha256 = "80f14c61e29ec1bcb4f0835d9eb3603a843f83b35bc4ae73720a33cfa42758fb",
    reproduced_sha256 = unname(sha(file.path(P$csv, "analysis_outputs", participant_file)))))
  cmp$identical <- !is.na(cmp$reproduced_sha256) & cmp$reported_sha256 == cmp$reproduced_sha256
  write.csv(cmp, file.path(P$qc, "comparison_with_reported_results.csv"), row.names = FALSE)
  message(
    sum(cmp$identical),
    " of ",
    nrow(cmp),
    " analysis outputs identical to results/analysis_outputs"
  )
  if (!all(cmp$identical)) warning("Some outputs differ from results/analysis_outputs; see output/qc/comparison_with_reported_results.csv")
}
# Table 1 QC: recomputed values must equal the original analysis output for shared rows.
t1 <- read.csv(file.path(P$csv, "table1_characteristics_long.csv"), stringsAsFactors = FALSE)
orig <- res("baseline_characteristics")
map <- c(
  age = "age", grip = "grip", height = "height", bmi = "bmi", sbp = "sbp", tc = "tc", hdl = "hdl",
  egfr = "egfr", risk = "risk", diabetes = "diabetes", smoking = "smoking", bp_tx = "bp_tx",
  statin = "statin", grip_low = "grip_low", female = "female"
)
qc <- bind_rows(lapply(
  names(map),
  function(nm) bind_rows(lapply(c("primary", "mobility"), function(k) {
    oc <- if (k == "primary") "Core age 50-79" else "Mobility age 60-79"
    o <- orig[orig$cohort == oc & orig$variable == map[[nm]], ]
    n <- t1[t1$cohort == k & t1$variable == nm, ]
    scale <- if (o$kind == "proportion_SE") 100 else 1
    data.frame(
      cohort = k, variable = nm, original = o$estimate * scale, recomputed = n$estimate,
      abs_diff = abs(o$estimate * scale - n$estimate)
    )
  }))
))
write.csv(qc, file.path(P$qc, "Table1_recomputation_check.csv"), row.names = FALSE)
stopifnot(max(qc$abs_diff) < 1e-8)

message("Table 1 recomputation check: max abs difference ", signif(max(qc$abs_diff), 3))
# Independent re-derivation.
out <- file.path(P$qc, "Independent_Check")
log <- file.path(P$logs, "independent_check.R.log")
status <- system2(file.path(R.home("bin"), "Rscript"), c(
  "--vanilla", shQuote(file.path(PKG, "code", "independent_check.R")),
  shQuote(BUILD), shQuote(out)
), stdout = log, stderr = log)
if (!identical(as.integer(status), 0L)) stop("Independent check failed; see ", log)
print(read.csv(file.path(out, "independent_check_results.csv")))
