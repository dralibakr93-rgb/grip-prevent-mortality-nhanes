# Reads the raw NHANES files, derives the analysis variables, computes PREVENT base risk, applies the
# eligibility criteria, and writes the master dataset, cohort flow, and grip-strength scaling.
# Run through code/01_data_derivation.R (working directory: the build directory).

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(survey)
})
options(survey.lonely.psu = "fail")
B <- "study"
out <- function(
  x,
  n
) write.csv(x,
  file.path(B, "05_Outputs", paste0(n, ".csv")),
  row.names = FALSE,
  na = ""
)
getv <- function(d, n) if (n %in% names(d)) as.numeric(d[[n]]) else rep(NA_real_, nrow(d))
bpmean <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) > 1) x <- x[-1]
  if (length(x)) mean(x) else NA_real_
}
egfr <- function(scr, age, sex) {
  k <- ifelse(sex == 2, .7, .9)
  a <- ifelse(sex == 2, -.241, -.302)
  142 * pmin(scr / k, 1)^a * pmax(scr / k, 1)^(-1.2) * .9938^age * ifelse(sex == 2, 1.012, 1)
}
mobility <- function(
  a,
  b,
  c
) case_when(
  a %in% 1 | b %in% 2:4 | c %in% 2:4 ~ 1,
  a %in% 2 & b %in% 1 & c %in% 1 ~ 0,
  TRUE ~ NA_real_
)
stopifnot(
  bpmean(c(120, 130, 140)) == 135,
  bpmean(c(NA, 130)) == 130,
  is.na(bpmean(c(0, NA))),
  abs(egfr(.9, 55, 1) - 100.8632) < .001
)
stopifnot(identical(
  mobility(c(1, 2, 2, NA), c(NA, 1, 5, NA), c(NA, 1, 1, NA)),
  c(1, 0, NA_real_, NA_real_)
))
dictionary <- as.data.frame(read_xpt(file.path(B, "02_Raw/RXQ_DRUG.xpt")))
cc <- grep("^RXD[DI]CI[1-6][ABC]$", names(dictionary), value = TRUE)
dictionary$statin <- rowSums(dictionary[cc] == 173, na.rm = TRUE) > 0
headers <- list()
medlogs <- list()
integrity <- list()
build <- function(year, suf) {
  read <- function(comp) {
    d <- as.data.frame(read_xpt(file.path(B, "02_Raw", year, paste0(comp, "_", suf, ".xpt"))))
    headers[[paste(year, comp)]] <<- data.frame(
      year = year,
      component = comp,
      variable = names(d),
      label = vapply(d, function(x) {
        a <- attr(x, "label")
        if (is.null(a)) "" else a
      }, character(1))
    )
    if (comp != "RXQ_RX") stopifnot(!anyDuplicated(d$SEQN))
    integrity[[paste(year, comp)]] <<- data.frame(
      year = year,
      component = comp,
      rows = nrow(d),
      unique_ID = length(unique(d$SEQN))
    )
    d
  }
  d <- read("DEMO")
  d$year <- year
  for (comp in c(
    "BMX",
    "BPX",
    "BPQ",
    "SMQ",
    "DIQ",
    "MCQ",
    "TCHOL",
    "HDL",
    "BIOPRO",
    "GHB",
    "MGX",
    "PFQ"
  )) {
    z <- read(comp)
    keep <- setdiff(names(z), setdiff(names(d), "SEQN"))
    d <- left_join(d, z[, keep, drop = FALSE], by = "SEQN")
  }
  rx <- read("RXQ_RX")
  drug <- toupper(trimws(rx$RXDDRUG))
  j <- match(trimws(rx$RXDDRGID), trimws(dictionary$RXDDRGID))
  hit <- grepl(
    "ATORVASTATIN|FLUVASTATIN|LOVASTATIN|PITAVASTATIN|PRAVASTATIN|ROSUVASTATIN|SIMVASTATIN|CERIVASTATIN",
    drug
  ) | dictionary$statin[j] %in% TRUE
  unknown <- is.na(drug) | drug == "" | grepl(
    "^(55555|77777|99999)$|UNKNOWN|UNSPECIFIED|UNIDENTIFIED",
    drug
  )
  r <- data.frame(SEQN = rx$SEQN, use = rx$RXDUSE, drug = drug, hit = hit, unknown = unknown)
  medlogs[[as.character(year)]] <<- r %>%
    count(drug, hit, unknown, name = "records") %>%
    mutate(year = year)
  r <- r %>%
    group_by(SEQN) %>%
    summarise(
      statin = case_when(any(hit) ~ 1, all(use %in% 2) ~ 0, all(use %in% 1) & !any(unknown) ~ 0, TRUE ~ NA_real_),
      .groups = "drop"
    )
  d <- left_join(d, r, by = "SEQN")
  d$age <- d$RIDAGEYR
  d$sex <- d$RIAGENDR
  d$bmi <- d$BMXBMI
  d$height <- d$BMXHT
  d$sbp <- apply(d[paste0("BPXSY", 1:4)], 1, bpmean)
  d$tc <- d$LBXTC
  d$hdl <- d$LBDHDD
  d$scr <- d$LBXSCR
  d$egfr <- egfr(d$scr, d$age, d$sex)
  d$diabetes <- case_when(
    d$DIQ010 %in% 1 | d$DIQ050 %in% 1 | getv(d, "DIQ070") %in% 1 ~ 1,
    d$DIQ010 %in% c(2, 3) ~ 0,
    TRUE ~ NA_real_
  )
  d$diabetes_a1c <- case_when(
    d$diabetes %in% 1 | d$LBXGH >= 6.5 ~ 1,
    d$diabetes %in% 0 & !is.na(d$LBXGH) & d$LBXGH < 6.5 ~ 0,
    TRUE ~ NA_real_
  )
  d$smoking <- case_when(
    d$SMQ020 == 1 & d$SMQ040 %in% c(1, 2) ~ 1,
    d$SMQ020 == 2 | (d$SMQ020 == 1 & d$SMQ040 == 3) ~ 0,
    TRUE ~ NA_real_
  )
  d$bp_tx <- case_when(
    d$BPQ050A == 1 ~ 1,
    d$BPQ050A == 2 ~ 0,
    is.na(d$BPQ050A) & (d$BPQ020 == 2 | (d$BPQ020 == 1 & d$BPQ040A == 2)) ~ 0,
    TRUE ~ NA_real_
  )
  cv <- d[paste0("MCQ160", LETTERS[2:6])]
  d$cvd_no <- rowSums(cv == 2, na.rm = TRUE) == 5
  d$cvd_yes <- rowSums(cv == 1, na.rm = TRUE) > 0
  d$mobility <- mobility(d$PFQ054, d$PFQ061B, d$PFQ061C)
  trial <- unlist(lapply(1:2, function(h) paste0("MGXH", h, "T", 1:3)))
  m <- as.matrix(d[trial])
  eff <- as.matrix(d[paste0(trial, "E")])
  m[is.na(eff) | eff != 1 | !is.finite(m) | m <= 0] <- NA
  d$nhand1 <- rowSums(!is.na(m[, 1:3]))
  d$nhand2 <- rowSums(!is.na(m[, 4:6]))
  d$both_hands <- d$nhand1 > 0 & d$nhand2 > 0
  d$grip_any <- apply(m, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
  d$grip <- ifelse(d$both_hands, d$grip_any, NA_real_)
  d$grip_low <- ifelse(!is.na(d$grip), as.integer(d$grip < ifelse(d$sex == 1, 27, 16)), NA_real_)
  d$hand_pain <- d$MGQ070 %in% 1 | d$MGQ100 %in% 1
  d$grip_status <- case_when(
    d$both_hands ~ "bilateral valid maximal effort",
    !is.na(d$grip_any) ~ "one hand valid maximal effort",
    d$MGDEXSTS %in% 3 ~ "not tested mixed reasons",
    d$MGDEXSTS %in% c(1, 2) ~ "tested but no qualifying effort",
    TRUE ~ "component unavailable"
  )
  d$w <- d$WTMEC2YR / 2
  d$design_ok <- !is.na(d$w) & d$w > 0 & !is.na(d$SDMVSTRA) & !is.na(d$SDMVPSU)
  input <- c(
    "age",
    "sex",
    "tc",
    "hdl",
    "sbp",
    "egfr",
    "bmi",
    "diabetes",
    "smoking",
    "bp_tx",
    "statin"
  )
  d$complete_score <- complete.cases(d[input])
  d$range_ok <- with(
    d,
    tc >= 130 & tc <= 320 & hdl >= 20 & hdl <= 100 & sbp >= 90 & sbp <= 200 & egfr >= 15 & egfr <= 140 & bmi >= 18.5 & bmi <= 39.9
  )
  d$range_ok[is.na(d$range_ok)] <- FALSE
  lines <- readLines(file.path(
    B,
    "02_Raw/Mortality",
    paste0("NHANES_", year, "_", year + 1, "_MORT_2019_PUBLIC.dat")
  ))
  numfield <- function(a, b) {
    x <- trimws(substr(lines, a, b))
    x[x %in% c("", ".")] <- NA_character_
    stopifnot(all(is.na(x) | grepl("^[0-9]+$", x)))
    as.numeric(x)
  }
  mm <- data.frame(
    SEQN = numfield(1, 6),
    eligible = numfield(15, 15),
    death = numfield(16, 16),
    cause = numfield(17, 19),
    months = numfield(46, 48)
  )
  stopifnot(!anyDuplicated(mm$SEQN))
  d <- left_join(d, mm, by = "SEQN")
  d$linked <- d$eligible %in% 1 & d$death %in% c(0, 1) & !is.na(d$months) & d$months >= 0
  d$time <- ifelse(d$months == 0, .5, d$months) / 12
  d$death_cvd <- ifelse(d$linked, as.integer(d$death == 1 & d$cause %in% c(1, 5)), NA_real_)
  d$death_other <- ifelse(d$linked,
    as.integer(d$death == 1 & !is.na(d$cause) & !(d$cause %in% c(1, 5))),
    NA_real_
  )
  d$clinical <- d$design_ok & d$age >= 50 & d$age <= 79 & d$cvd_no & !(d$RIDEXPRG %in% 1)
  d$pregrip <- d$clinical & d$complete_score & d$range_ok & d$linked & !is.na(d$height)
  d$core <- d$pregrip & d$both_hands
  d$mobility_cohort <- d$core & d$age >= 60 & !is.na(d$mobility)
  d$sex_factor <- factor(d$sex, levels = c(1, 2), labels = c("Men", "Women"))
  d$cycle <- factor(year)
  stopifnot(!anyDuplicated(d$SEQN))
  d
}
d <- bind_rows(build(2011, "G"), build(2013, "H"))
stopifnot(!anyDuplicated(d$SEQN))
coefs <- read.csv(file.path(B, "01_Sources/preventr_base_10yr_coefficients.csv"))
terms <- function(d) {
  a <- (d$age - 55) / 10
  n <- (d$tc - d$hdl) * .02586 - 3.5
  h <- (d$hdl * .02586 - 1.3) / .3
  s1 <- (pmin(d$sbp, 110) - 110) / 20
  s2 <- (pmax(d$sbp, 110) - 130) / 20
  b1 <- (pmin(d$bmi, 30) - 25) / 5
  b2 <- (pmax(d$bmi, 30) - 30) / 5
  e1 <- (pmin(d$egfr, 60) - 60) / -15
  e2 <- (pmax(d$egfr, 60) - 90) / -15
  cbind(
    a,
    n,
    h,
    s1,
    s2,
    d$diabetes,
    d$smoking,
    b1,
    b2,
    e1,
    e2,
    d$bp_tx,
    d$statin,
    d$bp_tx * s2,
    d$statin * n,
    a * n,
    a * h,
    a * s2,
    a * d$diabetes,
    a * d$smoking,
    a * b2,
    a * e1,
    1
  )
}
score <- function(d) {
  z <- terms(d)
  out <- rep(NA_real_, nrow(d))
  for (sx in 1:2) {
    i <- which(d$sex == sx & d$complete_score & d$range_ok)
    out[i] <- plogis(z[i, ,
      drop = FALSE
    ] %*% coefs[[if (sx == 1) "male_total_cvd" else "female_total_cvd"]])
  }
  out
}
d$risk <- score(d)
d$logitrisk <- qlogis(d$risk)
# Compare against archived outputs of the independent calculator for the original random sample.
set.seed(41003)
ii <- sample(which(d$core), min(100, sum(d$core)))
reference <- read.csv(file.path(B, "01_Sources/calculator_check_41003.csv"))
stopifnot(nrow(reference) == length(ii))
dif <- abs(d$risk[ii] - reference$total_cvd)
stopifnot(max(dif) < 1e-10, all(d$risk[d$core] > 0 & d$risk[d$core] < 1), !anyNA(d$risk[d$core]))
out(
  data.frame(test = "independent source arithmetic", profiles = length(ii), maximum_difference = max(dif)),
  "score_validation"
)
S <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~w,
  nest = TRUE,
  data = d[d$design_ok, ]
)
C <- subset(S, core)
scaling <- bind_rows(lapply(1:2, function(sx) {
  z <- subset(C, sex == sx)
  data.frame(
    sex = sx,
    mean_kg = as.numeric(svymean(~grip, z)),
    sd_kg = sqrt(as.numeric(svyvar(~grip, z)))
  )
}))
d$lower_z <- (scaling$mean_kg[match(d$sex, scaling$sex)] - d$grip) / scaling$sd_kg[match(
  d$sex,
  scaling$sex
)]
d$lower5 <- -d$grip / 5
d$joint <- factor(ifelse(d$mobility_cohort, paste0(ifelse(d$mobility == 0, "No limitation", "Limitation"), " / ", ifelse(d$grip_low == 0, "preserved strength", "low strength")), NA),
  levels = c(
    "No limitation / preserved strength",
    "No limitation / low strength",
    "Limitation / preserved strength",
    "Limitation / low strength"
  )
)
stages <- list(
  all_respondents = rep(TRUE, nrow(d)),
  positive_MEC_design = d$design_ok,
  age_50_79 = d$design_ok & d$age >= 50 & d$age <= 79,
  no_CVD_not_pregnant = d$clinical,
  complete_score = d$clinical & d$complete_score,
  score_range = d$clinical & d$complete_score & d$range_ok,
  linked_height = d$pregrip,
  core = d$core,
  mobility_age60_79 = d$mobility_cohort
)
flow <- bind_rows(lapply(names(stages), function(n) {
  q <- stages[[n]]
  data.frame(
    stage = n,
    N = sum(q),
    deaths = sum(d$death[q], na.rm = TRUE),
    cvd_deaths = sum(d$death_cvd[q], na.rm = TRUE)
  )
}))
out(flow, "cohort_flow")
out(scaling, "strength_scaling")
out(
  data.frame(cohort = c("core", "mobility"), N = c(sum(d$core), sum(d$mobility_cohort)), deaths = c(sum(d$death[d$core]), sum(d$death[d$mobility_cohort])), at_risk_5 = c(sum(d$time[d$core] >= 5), sum(d$time[d$mobility_cohort] >= 5)), cvd_deaths = c(sum(d$death_cvd[d$core]), sum(d$death_cvd[d$mobility_cohort])), median_followup = c(median(d$time[d$core]), median(d$time[d$mobility_cohort]))),
  "cohort_summary"
)
out(
  d %>% filter(mobility_cohort) %>% group_by(joint) %>% summarise(N = n(), deaths = sum(death), .groups = "drop"),
  "joint_group_counts"
)
out(
  d %>% filter(pregrip) %>% group_by(grip_status) %>% summarise(N = n(), deaths = sum(death), mean_age = mean(age), .groups = "drop"),
  "grip_selection"
)
out(
  data.frame(variable = c("age", "height", "bmi", "sbp", "tc", "hdl", "scr", "diabetes", "smoking", "bp_tx", "statin", "grip", "mobility"), missing = sapply(c("age", "height", "bmi", "sbp", "tc", "hdl", "scr", "diabetes", "smoking", "bp_tx", "statin", "grip", "mobility"), function(v) sum(is.na(d[[v]][d$clinical]))), denominator = sum(d$clinical)),
  "missingness"
)
out(
  data.frame(full_design_df = degf(S), core_design_df = degf(C), core_strata = length(unique(C$variables$SDMVSTRA)), core_psu = nrow(unique(C$variables[c("SDMVSTRA", "SDMVPSU")])), zero_months = sum(d$core & d$months %in% 0), unknown_cause = sum(d$core & d$death %in% 1 & is.na(d$cause))),
  "survey_audit"
)
out(bind_rows(headers), "source_variable_dictionary")
out(bind_rows(integrity), "source_key_checks")
out(bind_rows(medlogs), "medication_audit")
saveRDS(d, file.path(B, "03_Derived/study4_master.rds"))
saveRDS(
  list(score = score, terms = terms, coefs = coefs),
  file.path(B, "03_Derived/scoring_functions.rds")
)
capture.output(sessionInfo(), file = file.path(B, "08_Audit/sessionInfo_build.txt"))
print(flow)
print(scaling)
print(read.csv(file.path(B, "05_Outputs/joint_group_counts.csv")))
