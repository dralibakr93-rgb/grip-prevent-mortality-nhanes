# Stage 2. Survey-weighted characteristics of the primary and older cohorts (Table 1) and export of the
# analytic dataset (output/analytic_dataset/). No outcome model is fitted here.
# Survey design: pooled examination weight (WTMEC2YR/2), masked strata and PSUs, and the full design
# defined before the cohort is selected as a domain.

source("code/00_setup.R")
suppressPackageStartupMessages({
  library(survey)
  library(dplyr)
})
options(survey.lonely.psu = "fail")

d <- readRDS(file.path(P$derived, "study4_master.rds"))
d$race_eth <- factor(
  dplyr::case_when(
    d$RIDRETH3 %in% c(1, 2) ~ "Hispanic",
    d$RIDRETH3 %in% 3 ~ "Non-Hispanic White",
    d$RIDRETH3 %in% 4 ~ "Non-Hispanic Black",
    d$RIDRETH3 %in% 6 ~ "Non-Hispanic Asian",
    d$RIDRETH3 %in% 7 ~ "Other or multiracial"
  ),
  levels = c(
    "Non-Hispanic White",
    "Non-Hispanic Black",
    "Hispanic",
    "Non-Hispanic Asian",
    "Other or multiracial"
  )
)
d$education <- factor(
  dplyr::case_when(
    d$DMDEDUC2 %in% 1:2 ~ "Less than high school",
    d$DMDEDUC2 %in% 3 ~ "High school or equivalent",
    d$DMDEDUC2 %in% 4:5 ~ "More than high school"
  ),
  levels = c("Less than high school", "High school or equivalent", "More than high school")
)
d$pir <- d$INDFMPIR
d$female <- as.numeric(d$sex == 2)
d <- add_function_hospital_vars(d) # prior-year hospitalization; ADL/IADL items
# Base PREVENT 10-year ASCVD risk, same inputs and coefficient file as total CVD.
sf <- readRDS(file.path(P$derived, "scoring_functions.rds"))
d$risk_ascvd <- local({
  z <- sf$terms(d)
  out <- rep(NA_real_, nrow(d))
  for (sx in 1:2) {
    i <- which(d$sex == sx & d$complete_score & d$range_ok)
    out[i] <- plogis(z[i, , drop = FALSE] %*% sf$coefs[[paste0(c("male", "female")[sx], "_ascvd")]])
  }
  out
})

S <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~w,
  nest = TRUE,
  data = d[d$design_ok, ]
)
cohorts <- list(primary = subset(S, core), mobility = subset(S, mobility_cohort))

rows <- list()
add <- function(...) rows[[length(rows) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)
for (k in names(cohorts)) {
  des <- cohorts[[k]]
  v <- des$variables
  add(
    cohort = k, section = "Counts", variable = "participants", level = "", type = "count",
    n = nrow(v), estimate = nrow(v), spread = NA, n_missing = 0
  )
  add(
    cohort = k, section = "Counts", variable = "deaths", level = "", type = "count",
    n = sum(v$death), estimate = sum(v$death), spread = NA, n_missing = 0
  )
  add(
    cohort = k, section = "Counts", variable = "median_followup_years", level = "", type = "median",
    n = nrow(v), estimate = median(v$time), spread = NA, n_missing = 0
  )
  cont <- c(
    age = "age", grip = "grip", height = "height", bmi = "bmi", sbp = "sbp", tc = "tc", hdl = "hdl",
    egfr = "egfr", risk = "risk", risk_ascvd = "risk_ascvd", pir = "pir"
  )
  for (nm in names(cont)) {
    x <- cont[[nm]]
    ok <- !is.na(v[[x]])
    dd <- subset(des, !is.na(get(x)))
    m <- coef(svymean(as.formula(paste0("~", x)), dd))
    s <- sqrt(coef(svyvar(as.formula(paste0("~", x)), dd)))
    add(
      cohort = k, section = "Continuous", variable = nm, level = "", type = "mean_sd",
      n = sum(ok), estimate = unname(m), spread = unname(s), n_missing = sum(!ok)
    )
  }
  for (sx in 1:2) {
    dd <- subset(des, sex == sx)
    add(
      cohort = k,
      section = "Continuous",
      variable = paste0("grip_", c("men", "women")[sx]),
      level = "",
      type = "mean_sd", n = sum(v$sex == sx), estimate = unname(coef(svymean(~grip, dd))),
      spread = unname(sqrt(coef(svyvar(~grip, dd)))), n_missing = 0
    )
  }
  bin <- c(
    "female",
    "diabetes",
    "smoking",
    "bp_tx",
    "statin",
    "grip_low",
    "prior_hosp",
    if (k == "mobility") "mobility"
  )
  for (x in bin) {
    ok <- !is.na(v[[x]])
    dd <- subset(des, !is.na(get(x)))
    add(
      cohort = k, section = "Binary", variable = x, level = "1", type = "n_pct",
      n = sum(v[[x]][ok] == 1),
      estimate = 100 * unname(coef(svymean(as.formula(paste0("~", x)), dd))),
      spread = NA, n_missing = sum(!ok)
    )
  }
  for (x in c("race_eth", "education")) {
    ok <- !is.na(v[[x]])
    dd <- subset(des, !is.na(get(x)))
    p <- coef(svymean(as.formula(paste0("~", x)), dd))
    for (lv in levels(v[[x]])) {
      add(
        cohort = k, section = "Categorical", variable = x, level = lv, type = "n_pct",
        n = sum(v[[x]] == lv, na.rm = TRUE), estimate = 100 * unname(p[paste0(x, lv)]), spread = NA,
        n_missing = sum(!ok)
      )
    }
  }
}
t1 <- bind_rows(rows)
write.csv(t1, file.path(P$csv, "table1_characteristics_long.csv"), row.names = FALSE)

# Analytic dataset: full positive-weight examination design sample with domain flags, so that
# survey variance for any domain can be reproduced (domains must not be subset before design).
keep <- c(
  "SEQN",
  "cycle",
  "SDMVSTRA",
  "SDMVPSU",
  "WTMEC2YR",
  "w",
  "design_ok",
  "clinical",
  "pregrip",
  "core",
  "mobility_cohort",
  "age",
  "sex",
  "race_eth",
  "education",
  "pir",
  "height",
  "bmi",
  "sbp",
  "tc",
  "hdl",
  "scr",
  "egfr",
  "diabetes",
  "diabetes_a1c",
  "LBXGH",
  "smoking",
  "bp_tx",
  "statin",
  "risk",
  "logitrisk",
  "risk_ascvd",
  "grip",
  "grip_any",
  "nhand1",
  "nhand2",
  "both_hands",
  "grip_status",
  "grip_low",
  "lower_z",
  "lower5",
  "hand_pain",
  "PFQ054",
  "PFQ061B",
  "PFQ061C",
  "mobility",
  "PFQ061A",
  "PFQ061F",
  "PFQ061G",
  "PFQ061H",
  "PFQ061J",
  "PFQ061K",
  "PFQ061L",
  "adl",
  "iadl",
  "HUQ071",
  "prior_hosp",
  "MCQ220",
  "linked",
  "months",
  "time",
  "death",
  "cause", "death_cvd"
)
ad <- d[d$design_ok, keep]
ad$race_eth <- as.character(ad$race_eth)
ad$education <- as.character(ad$education)
write.csv(ad,
  file.path(P$analytic, "study4_analytic_design_sample.csv"),
  row.names = FALSE,
  na = ""
)
saveRDS(ad, file.path(P$analytic, "study4_analytic_design_sample.rds"))
message(
  "Analytic dataset rows: ",
  nrow(ad),
  "; primary cohort: ",
  sum(ad$core),
  "; mobility cohort: ",
  sum(ad$mobility_cohort)
)
