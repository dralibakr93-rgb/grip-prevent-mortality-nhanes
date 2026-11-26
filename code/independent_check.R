# Independent re-derivation of the primary cohort, the PREVENT scores, and the primary hazard ratio.
# It re-reads the raw files with separate code, parses the mortality files separately, and computes
# PREVENT with preventr::estimate_risk() instead of the transcribed coefficients.

suppressPackageStartupMessages({
  library(haven)
  library(survey)
  library(survival)
  library(splines)
})
options(survey.lonely.psu = "fail")
args <- commandArgs(trailingOnly = TRUE)
REPRO <- args[1] # path to .../05_Reproducibility
OUT <- args[2] # output directory
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
RAW <- file.path(REPRO, "study", "02_Raw")
.libPaths(c(file.path(REPRO, "R_library"), .libPaths()))
suppressPackageStartupMessages(library(preventr))

rd <- function(cyc, comp) {
  sfx <- if (cyc == 2011) "G" else "H"
  as.data.frame(read_xpt(file.path(RAW, cyc, sprintf("%s_%s.xpt", comp, sfx))))
}
drugdict <- as.data.frame(read_xpt(file.path(RAW, "RXQ_DRUG.xpt")))
classcols <- grep("^RXD[DI]CI[1-6][ABC]$", names(drugdict), value = TRUE)
statin_ids <- drugdict$RXDDRGID[apply(drugdict[classcols], 1, function(r) any(r %in% 173))]
statin_words <- c(
  "ATORVASTATIN", "FLUVASTATIN", "LOVASTATIN", "PITAVASTATIN", "PRAVASTATIN",
  "ROSUVASTATIN", "SIMVASTATIN", "CERIVASTATIN"
)

one_cycle <- function(cyc) {
  keep <- function(df, vars) df[, intersect(c("SEQN", vars), names(df)), drop = FALSE]
  demo <- keep(
    rd(cyc, "DEMO"),
    c("RIDAGEYR", "RIAGENDR", "RIDEXPRG", "WTMEC2YR", "SDMVPSU", "SDMVSTRA")
  )
  parts <- list(
    keep(rd(cyc, "BMX"), c("BMXBMI", "BMXHT")),
    keep(rd(cyc, "BPX"), paste0("BPXSY", 1:4)),
    keep(rd(cyc, "BPQ"), c("BPQ020", "BPQ040A", "BPQ050A")),
    keep(rd(cyc, "SMQ"), c("SMQ020", "SMQ040")),
    keep(rd(cyc, "DIQ"), c("DIQ010", "DIQ050", "DIQ070")),
    keep(rd(cyc, "MCQ"), c(paste0("MCQ160", c("B", "C", "D", "E", "F")), "MCQ220")),
    keep(rd(cyc, "TCHOL"), "LBXTC"), keep(rd(cyc, "HDL"), "LBDHDD"),
    keep(rd(cyc, "BIOPRO"), "LBXSCR"), keep(rd(cyc, "GHB"), "LBXGH"),
    keep(rd(cyc, "MGX"), c(
      paste0("MGXH", rep(1:2, each = 3), "T", 1:3),
      paste0("MGXH", rep(1:2, each = 3), "T", 1:3, "E")
    )),
    keep(rd(cyc, "PFQ"), c("PFQ054", "PFQ061B", "PFQ061C"))
  )
  x <- demo
  for (p in parts) x <- merge(x, p, by = "SEQN", all.x = TRUE, sort = FALSE)
  # Statins: affirmative by ingredient name or drug-class 173; known non-use otherwise.
  rx <- rd(cyc, "RXQ_RX")
  nm <- toupper(trimws(rx$RXDDRUG))
  is_statin <- grepl(
    paste(statin_words, collapse = "|"),
    nm
  ) | trimws(rx$RXDDRGID) %in% trimws(statin_ids)
  unknown <- nm == "" | is.na(nm) | grepl(
    "^(55555|77777|99999)$",
    nm
  ) | grepl(
    "UNKNOWN|UNSPECIFIED|UNIDENTIFIED",
    nm
  )
  agg <- do.call(
    rbind,
    lapply(split(data.frame(is_statin, unknown, use = rx$RXDUSE), rx$SEQN), function(g)
      data.frame(statin = if (any(g$is_statin)) 1 else if (all(g$use %in% 2)) 0 else if (all(g$use %in% 1) && !any(g$unknown)) 0 else NA))
  )
  agg$SEQN <- as.numeric(rownames(agg))
  x <- merge(x, agg, by = "SEQN", all.x = TRUE, sort = FALSE)
  # Mortality: independent read_fwf-style parser using documented 2019 LMF positions.
  lmf <- read.fwf(file.path(RAW, "Mortality", sprintf("NHANES_%d_%d_MORT_2019_PUBLIC.dat", cyc, cyc + 1)),
    widths = c(6, 8, 1, 1, 3, 1, 1, 21, 3, 3), colClasses = "character",
    col.names = c(
      "SEQN",
      "f1",
      "ELIGSTAT",
      "MORTSTAT",
      "UCOD",
      "DIAB",
      "HYP",
      "f2",
      "PERMTH_INT",
      "PERMTH_EXM"
    )
  )
  lmf <- data.frame(
    SEQN = as.numeric(lmf$SEQN), ELIGSTAT = suppressWarnings(as.numeric(lmf$ELIGSTAT)),
    MORTSTAT = suppressWarnings(as.numeric(lmf$MORTSTAT)),
    UCOD = suppressWarnings(as.numeric(lmf$UCOD)),
    PERMTH_EXM = suppressWarnings(as.numeric(lmf$PERMTH_EXM))
  )
  x <- merge(x, lmf, by = "SEQN", all.x = TRUE, sort = FALSE)
  x$cycle <- cyc
  x
}
d <- rbind(one_cycle(2011), one_cycle(2013))
stopifnot(!anyDuplicated(d$SEQN))

# Derived variables, coded from the manuscript description.
sbp_mat <- as.matrix(d[paste0("BPXSY", 1:4)])
d$sbp <- apply(sbp_mat, 1, function(v) {
  v <- v[!is.na(v) & v > 0]
  if (length(v) >= 2) mean(v[-1]) else if (length(v) == 1) v else NA
})
fem <- d$RIAGENDR == 2
k <- ifelse(fem, 0.7, 0.9)
a <- ifelse(fem, -0.241, -0.302)
d$egfr <- 142 * pmin(
  d$LBXSCR / k,
  1
)^a * pmax(
  d$LBXSCR / k,
  1
)^-1.2 * 0.9938^d$RIDAGEYR * ifelse(fem,
  1.012,
  1
)
d$dm <- ifelse(d$DIQ010 %in% 1 | d$DIQ050 %in% 1 | d$DIQ070 %in% 1,
  1,
  ifelse(d$DIQ010 %in% c(2, 3), 0, NA)
)
d$smk <- ifelse(d$SMQ020 %in% 1 & d$SMQ040 %in% c(1, 2),
  1,
  ifelse(d$SMQ020 %in% 2 | (d$SMQ020 %in% 1 & d$SMQ040 %in% 3), 0, NA)
)
d$bptx <- ifelse(d$BPQ050A %in% 1, 1, ifelse(d$BPQ050A %in% 2, 0,
  ifelse(is.na(d$BPQ050A) & (d$BPQ020 %in% 2 | (d$BPQ020 %in% 1 & d$BPQ040A %in% 2)), 0, NA)
))
cvd <- as.matrix(d[paste0("MCQ160", c("B", "C", "D", "E", "F"))])
d$no_cvd <- rowSums(cvd == 2, na.rm = TRUE) == 5
trials <- paste0("MGXH", rep(1:2, each = 3), "T", 1:3)
g <- as.matrix(d[trials])
e <- as.matrix(d[paste0(trials, "E")])
g[!(e %in% 1) | is.na(g) | g <= 0] <- NA
h1 <- rowSums(!is.na(g[, 1:3])) > 0
h2 <- rowSums(!is.na(g[, 4:6])) > 0
d$grip <- ifelse(h1 & h2, apply(g, 1, function(v) max(v, na.rm = TRUE)), NA)
d$w <- d$WTMEC2YR / 2
d$design_ok <- !is.na(d$w) & d$w > 0 & !is.na(d$SDMVPSU) & !is.na(d$SDMVSTRA)

# Sequential eligibility (manuscript Methods / Figure 1).
inp <- c(
  "RIDAGEYR",
  "RIAGENDR",
  "LBXTC",
  "LBDHDD",
  "sbp",
  "egfr",
  "BMXBMI",
  "dm",
  "smk",
  "bptx",
  "statin"
)
st <- list()
st$all <- rep(TRUE, nrow(d))
st$design <- d$design_ok
st$age <- st$design & d$RIDAGEYR >= 50 & d$RIDAGEYR <= 79
st$clinical <- st$age & d$no_cvd & !(d$RIDEXPRG %in% 1)
st$complete <- st$clinical & complete.cases(d[inp])
inr <- with(d, LBXTC >= 130 & LBXTC <= 320 & LBDHDD >= 20 & LBDHDD <= 100 & sbp >= 90 & sbp <= 200 &
  egfr >= 15 & egfr <= 140 & BMXBMI >= 18.5 & BMXBMI <= 39.9)
st$range <- st$complete & inr %in% TRUE
st$linked <- st$range & d$ELIGSTAT %in% 1 & d$MORTSTAT %in% 0:1 & !is.na(d$PERMTH_EXM) & !is.na(d$BMXHT)
st$core <- st$linked & !is.na(d$grip)
flow <- data.frame(
  stage = names(st),
  N = sapply(st, sum),
  deaths = sapply(st, function(s) sum(d$MORTSTAT[s] %in% 1))
)
write.csv(flow, file.path(OUT, "independent_flow.csv"), row.names = FALSE)

# PREVENT via the package's own estimator (base model, 10-year total CVD).
core <- d[st$core, ]
pr <- function(i) {
  r <- core[i, ]
  if (r$sbp > 180) return(NA_real_) # preventr limits SBP to 180; the AHA equations allow 200
  out <- estimate_risk(
    age = r$RIDAGEYR, sex = if (r$RIAGENDR == 1) "male" else "female", sbp = r$sbp,
    bp_tx = r$bptx == 1, total_c = r$LBXTC, hdl_c = r$LBDHDD, statin = r$statin == 1,
    dm = r$dm == 1, smoking = r$smk == 1, egfr = r$egfr, bmi = r$BMXBMI,
    model = "base", time = "10yr", quiet = TRUE
  )
  as.numeric(out$total_cvd)
}
core$prevent_pkg <- vapply(seq_len(nrow(core)), pr, numeric(1))
write.csv(data.frame(SEQN = core$SEQN, prevent_pkg = core$prevent_pkg),
  file.path(OUT, "independent_prevent.csv"),
  row.names = FALSE
)

# Compare with the derived master file (participant-level agreement).
m <- readRDS(file.path(REPRO, "study", "03_Derived", "study4_master.rds"))
mc <- m[m$core, ]
same_ids <- setequal(mc$SEQN, core$SEQN)
j <- match(core$SEQN, mc$SEQN)
risk_diff <- abs(core$prevent_pkg - mc$risk[j])
cmp <- data.frame(
  check = c(
    "Primary cohort SEQN sets identical", "Primary cohort N (independent)", "Deaths (independent)",
    "Max |grip difference| (kg)", "Max |SBP difference| (mm Hg)", "Max |eGFR difference|",
    "Participants with SBP >180 (preventr cannot score)",
    "Max |PREVENT difference|, package vs study (SBP<=180)",
    "Package rounding: PREVENT returned to 3 decimals"
  ),
  value = c(
    same_ids, nrow(core), sum(core$MORTSTAT == 1),
    max(abs(core$grip - mc$grip[j])),
    max(abs(core$sbp - mc$sbp[j])),
    max(abs(core$egfr - mc$egfr[j])),
    sum(core$sbp > 180), max(risk_diff, na.rm = TRUE), "see note"
  )
)
# Survey-weighted descriptive checks.
des <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~w,
  nest = TRUE,
  data = d[d$design_ok, ]
)
dd <- des$variables
dd$core <- dd$SEQN %in% core$SEQN
des$variables <- dd
C <- subset(des, core)
cmp <- rbind(cmp, data.frame(
  check = c("Weighted mean age", "Weighted % women", "Core design df (PSU - strata)", "Sum of pooled weights, full design (millions)"),
  value = c(
    round(coef(svymean(~RIDAGEYR, C)), 2),
    round(100 * coef(svymean(~ I(RIAGENDR == 2), C))[2], 2),
    degf(C),
    round(sum(d$w[d$design_ok]) / 1e6, 2)
  )
))
# Independent primary model with the same prespecified structure.
C$variables$logitrisk <- qlogis(mc$risk[match(C$variables$SEQN, mc$SEQN)])
sexw <- sapply(1:2, function(s) {
  z <- subset(C, RIAGENDR == s)
  c(coef(svymean(~grip, z)), sqrt(coef(svyvar(~grip, z))))
})
C$variables$lower_z <- (sexw[1, C$variables$RIAGENDR] - C$variables$grip) / sexw[
  2,
  C$variables$RIAGENDR
]
rk <- ns(C$variables$logitrisk, df = 3)
ak <- ns(C$variables$RIDAGEYR, df = 3)
C$variables$time <- C$variables$PERMTH_EXM / 12
C$variables$death <- C$variables$MORTSTAT
f <- Surv(time, death) ~ lower_z + rk + ak + factor(RIAGENDR) + factor(cycle) + BMXHT + BMXBMI
C$variables$rk <- rk
C$variables$ak <- ak
fit <- svycoxph(f, design = C)
b <- coef(fit)["lower_z"]
se <- sqrt(vcov(fit)["lower_z", "lower_z"])
p <- length(coef(fit))
ci <- function(df) exp(b + c(-1, 1) * qt(0.975, df) * se)
cmp <- rbind(cmp, data.frame(
  check = c(
    "Sex-specific weighted SD (men, women)",
    "Primary HR (independent)",
    "95% CI, df = design df - p + 1 (22)",
    "95% CI, df = design df (32), alternative convention", "95% CI, normal approximation"
  ),
  value = c(
    paste(round(sexw[2, ], 3), collapse = ", "), round(exp(b), 4),
    paste(round(ci(degf(C) - p + 1), 3), collapse = "-"),
    paste(round(ci(degf(C)), 3), collapse = "-"),
    paste(round(exp(b + c(-1, 1) * 1.959964 * se), 3), collapse = "-")
  )
))
write.csv(cmp, file.path(OUT, "independent_check_results.csv"), row.names = FALSE)
print(flow)
print(cmp, right = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_independent_check.txt"))
