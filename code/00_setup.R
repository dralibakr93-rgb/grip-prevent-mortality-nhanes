# Shared setup for all stages: paths, required packages, and staging of the input files.
# Every stage is run from the repository root, for example: Rscript code/03_primary_analysis.R
# The scripts in code/analysis/ use paths relative to a build directory (default: _build/).
# stage_build() copies the raw files from data/raw/ into that directory, and run_analysis()
# runs one of those scripts there in a clean R session.

PKG <- normalizePath(Sys.getenv("STUDY4_PKG", getwd()), mustWork = TRUE)
if (!dir.exists(file.path(PKG, "code")) || !dir.exists(file.path(PKG, "data")))
  stop("Run from the repository root (the folder that contains code/ and data/), or set STUDY4_PKG.")
BUILD <- Sys.getenv("STUDY4_BUILD", file.path(PKG, "_build"))
dir.create(BUILD, showWarnings = FALSE, recursive = TRUE)
BUILD <- normalizePath(BUILD)

P <- list(
  raw = file.path(PKG, "data", "raw"),
  coef = file.path(PKG, "data", "reference"),
  analytic = file.path(PKG, "output", "analytic_dataset"),
  analysis = file.path(PKG, "code", "analysis"),
  rlib = file.path(PKG, "environment", "R_library"),
  out_orig = file.path(BUILD, "study", "05_Outputs"),
  out_add = file.path(BUILD, "study", "additional", "Results"),
  derived = file.path(BUILD, "study", "03_Derived"),
  derived_add = file.path(BUILD, "study", "additional", "Derived"),
  csv = file.path(PKG, "output", "csv"),
  fig_main = file.path(PKG, "output", "figures", "main"),
  fig_supp = file.path(PKG, "output", "figures", "supplementary"),
  fig_src = file.path(PKG, "output", "figures", "ggplot_objects"),
  tab_main = file.path(PKG, "output", "tables"),
  qc = file.path(PKG, "output", "qc"),
  logs = file.path(PKG, "output", "qc", "logs")
)
for (d in P[c("csv", "fig_main", "fig_supp", "fig_src", "tab_main", "qc", "logs", "analytic")])
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

required_pkgs <- c(
  "haven", "dplyr", "survey", "survival", "Matrix", "splines", "Hmisc", "mice", "mitools",
  "ggplot2", "patchwork", "scales", "svglite", "ragg", "officer", "flextable", "jsonlite", "digest"
)
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing R packages: ", paste(missing, collapse = ", "))

# Stage public inputs into the build directory.
stage_build <- function() {
  S4 <- file.path(BUILD, "study")
  for (d in c(
    "01_Sources", "02_Raw/2011", "02_Raw/2013", "02_Raw/Mortality", "03_Derived", "04_Code",
    "05_Outputs", "08_Audit"
  ))
    dir.create(file.path(S4, d), recursive = TRUE, showWarnings = FALSE)
  cp <- function(from, to) stopifnot(file.copy(from, to, overwrite = TRUE, copy.date = TRUE))
  for (f in list.files(file.path(P$raw, "2011_2012"), "\\.xpt$", full.names = TRUE))
    if (!grepl("^ALB_CR", basename(f))) cp(f, file.path(S4, "02_Raw/2011", basename(f)))
  for (f in list.files(file.path(P$raw, "2013_2014"), "\\.xpt$", full.names = TRUE))
    if (!grepl("^ALB_CR", basename(f))) cp(f, file.path(S4, "02_Raw/2013", basename(f)))
  for (f in list.files(file.path(P$raw, "Mortality_2019_LMF"), "\\.dat$", full.names = TRUE))
    cp(f, file.path(S4, "02_Raw/Mortality", basename(f)))
  cp(file.path(P$raw, "Reference", "RXQ_DRUG.xpt"), file.path(S4, "02_Raw", "RXQ_DRUG.xpt"))
  cp(
    file.path(P$coef, "preventr_base_10yr_coefficients.csv"),
    file.path(S4, "01_Sources", "preventr_base_10yr_coefficients.csv")
  )
  cp(file.path(P$coef, "calculator_check_41003.csv"), file.path(S4, "01_Sources", "calculator_check_41003.csv"))
  for (f in list.files(P$analysis, "\\.R$", full.names = TRUE)) cp(
    f,
    file.path(S4, "04_Code", basename(f))
  )
  # Urine albumin/creatinine inputs used by 04_additional_sensitivity.R.
  S3 <- file.path(BUILD, "study", "02_Raw", "urine")
  dir.create(file.path(S3, "2011"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(S3, "2013"), recursive = TRUE, showWarnings = FALSE)
  cp(file.path(P$raw, "2011_2012", "ALB_CR_G.xpt"), file.path(S3, "2011", "ALB_CR_G.xpt"))
  cp(file.path(P$raw, "2013_2014", "ALB_CR_H.xpt"), file.path(S3, "2013", "ALB_CR_H.xpt"))
  # Bundled preventr 0.12.0 used by 04_additional_sensitivity.R.
  RL <- file.path(BUILD, "R_library")
  dir.create(RL, recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(P$rlib, "preventr"), RL, recursive = TRUE, overwrite = TRUE, copy.date = TRUE)
  invisible(TRUE)
}

# Run one analysis script in a clean R session with the build directory as working directory.
run_analysis <- function(script) {
  log <- file.path(P$logs, paste0(script, ".log"))
  rscript <- file.path(R.home("bin"), "Rscript")
  t0 <- Sys.time()
  status <- withr_wd(
    BUILD,
    system2(rscript, c("--vanilla", shQuote(file.path("study", "04_Code", script))),
      stdout = log, stderr = log
    )
  )
  if (!identical(as.integer(status), 0L)) stop(script, " failed; see ", log)
  message(sprintf(
    "%s finished in %.1f s",
    script,
    as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ))
  invisible(log)
}
withr_wd <- function(dir, expr) {
  old <- setwd(dir)
  on.exit(setwd(old))
  force(expr)
}

# Locate an analysis output CSV written by stages 01-05.
res_path <- function(name) {
  f <- c(file.path(P$out_orig, paste0(name, ".csv")), file.path(P$out_add, paste0(name, ".csv")))
  f <- f[file.exists(f)]
  if (!length(f)) stop("Result file not found: ", name, " (run stages 01-05 first)")
  f[1]
}
res <- function(name) read.csv(res_path(name), stringsAsFactors = FALSE, check.names = FALSE)

# Prior-year overnight hospitalization (HUQ071)
# and ADL/IADL limitation (PFQ061 items; any of codes 2-4 = limitation, all items code 1 = none, else unclassified).
add_function_hospital_vars <- function(d) {
  huq <- rbind(
    as.data.frame(haven::read_xpt(file.path(P$raw, "2011_2012", "HUQ_G.xpt")))[c("SEQN", "HUQ071")],
    as.data.frame(haven::read_xpt(file.path(P$raw, "2013_2014", "HUQ_H.xpt")))[c("SEQN", "HUQ071")]
  )
  stopifnot(!anyDuplicated(huq$SEQN))
  d$HUQ071 <- as.integer(huq$HUQ071[match(d$SEQN, huq$SEQN)])
  d$prior_hosp <- ifelse(d$HUQ071 %in% 1, 1, ifelse(d$HUQ071 %in% 2, 0, NA))
  code <- function(v) ifelse(v %in% 2:4, 1, ifelse(v %in% 1, 0, NA))
  anylim <- function(vars) {
    m <- sapply(vars, function(v) code(d[[v]]))
    ifelse(rowSums(m == 1, na.rm = TRUE) > 0, 1, ifelse(rowSums(!is.na(m)) == length(vars), 0, NA))
  }
  d$adl <- anylim(c("PFQ061H", "PFQ061J", "PFQ061K", "PFQ061L"))
  d$iadl <- anylim(c("PFQ061A", "PFQ061F", "PFQ061G"))
  d
}

set.seed(20261003)
