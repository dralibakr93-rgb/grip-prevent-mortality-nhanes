# Stage 5d. Analyses within the guideline PREVENT categories:
#   A. low grip strength and functional limitation by PREVENT-ASCVD category and below or above the
#      treatment thresholds, with 5-year mortality (Figures 4, S7, and S8);
#   B. the creatinine-based kidney term: kidney-neutral PREVENT and an eGFR spline;
#   C. PREVENT-ASCVD and PREVENT heart failure as alternative adjustment scores;
#   D. selection at each exclusion step (Figure S1).

source("code/00_setup.R")
suppressPackageStartupMessages({
  library(survey)
  library(survival)
  library(splines)
})
options(survey.lonely.psu = "fail")
OUT <- file.path(P$csv, "analysis_outputs")
dir.create(OUT, showWarnings = FALSE)
wcsv <- function(x, n) write.csv(x, file.path(OUT, paste0(n, ".csv")), row.names = FALSE)

A <- readRDS(file.path(P$derived, "analysis_models.rds"))
d <- A$data
sf <- readRDS(file.path(P$derived, "scoring_functions.rds"))
d <- add_function_hospital_vars(d)

# A1. PREVENT outcomes (and B1 kidney-neutral)
score_col <- function(x, col) {
  z <- sf$terms(x)
  out <- rep(NA_real_, nrow(x))
  for (sx in 1:2) {
    i <- which(x$sex == sx & x$complete_score & x$range_ok)
    out[i] <- plogis(z[i, , drop = FALSE] %*% sf$coefs[[paste0(c("male", "female")[sx], "_", col)]])
  }
  out
}
d$r_tot <- score_col(d, "total_cvd")
d$r_asc <- score_col(d, "ascvd")
d$r_hf <- score_col(d, "heart_failure")
stopifnot(max(abs(d$r_tot - d$risk), na.rm = TRUE) < 1e-12)
dn <- d
dn$egfr <- 90
d$r_neut <- score_col(dn, "total_cvd")
# The reference file follows the draw order of the original calculator check.
set.seed(20261118)
ii <- sample(which(d$core), 100)
reference <- read.csv(file.path(P$coef, "calculator_check_20261118.csv"))
stopifnot(nrow(reference) == length(ii))
chk <- rbind(abs(d$r_asc[ii] - reference$ascvd), abs(d$r_hf[ii] - reference$heart_failure))
stopifnot(max(chk) < 1e-10)
wcsv(data.frame(
  check = c("PREVENT-ASCVD vs extracted calculator arithmetic", "PREVENT-HF vs extracted calculator arithmetic"),
  profiles = length(ii), maximum_difference = apply(chk, 1, max)
), "ascvd_hf_score_validation")

# definitions
quant <- res("risk_strength_contrasts")
d$low_grip <- as.integer(d$grip <= quant$q25[match(d$sex, quant$sex)])
d$asc5 <- as.integer(d$r_asc >= 0.05)
d$tot75 <- as.integer(d$r_tot >= 0.075)
d$asc_cat <- as.character(cut(d$r_asc,
  c(0, .03, .05, .10, 1),
  right = FALSE,
  labels = c("<3%", "3% to <5%", "5% to <10%", ">=10%")
))
d$older <- d$core & d$age >= 60
d$anyfx <- ifelse(d$mobility %in% 1 | d$adl %in% 1 | d$iadl %in% 1, 1,
  ifelse(d$mobility %in% 0 & d$adl %in% 0 & d$iadl %in% 0, 0, NA)
)
stopifnot(
  sum(d$core) == 2650,
  sum(d$core & d$asc5 %in% 1) == 1318,
  sum(d$core & d$low_grip %in% 1) == 836
)
# Derived variables used by the guideline-category and kidney-term analyses, for the analysis workbook.
gv <- d[
  d$design_ok,
  c("SEQN", "r_tot", "r_asc", "r_hf", "r_neut", "asc_cat", "asc5", "tot75", "low_grip", "anyfx")
]
names(gv) <- c(
  "SEQN",
  "prevent_total_cvd",
  "prevent_ascvd",
  "prevent_hf",
  "prevent_total_cvd_kidney_neutral",
  "ascvd_category",
  "ascvd_ge5", "total_cvd_ge7_5", "low_grip", "any_functional_limitation"
)
wcsv(gv, "guideline_analysis_variables")

S <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~w,
  nest = TRUE,
  data = d[d$design_ok, ]
)
R <- as.svrepdesign(S, type = "JKn", mse = TRUE)
base <- "r1+r2+r3+a1+a2+a3+sex_factor+cycle+height+bmi"

# A4. Burden within guideline categories
measures <- list(
  low_grip = "core",
  mobility = "older",
  adl = "older",
  iadl = "older",
  anyfx = "older"
)
cats <- list(
  asc_cat = c("<3%", "3% to <5%", "5% to <10%", ">=10%"),
  asc5 = c(0, 1),
  tot75 = c(0, 1)
)
burden <- list()
Sv <- S$variables
for (m in names(measures)) {
  for (cv in names(cats)) for (lv in cats[[cv]]) {
    k <- Sv[[measures[[m]]]] & !is.na(Sv[[m]]) & Sv[[cv]] %in% lv
    D <- subset(S, k)
    p <- svyciprop(as.formula(paste0("~I(", m, "==1)")), D, method = "logit")
    tt <- svytotal(as.formula(paste0("~I(", m, "==1)")), D)
    df <- degf(D)
    tot <- coef(tt)[[2]]
    se <- SE(tt)[[2]]
    h <- qt(.975, df) * se / tot # log-scale CI (always positive)
    burden[[length(burden) + 1]] <- data.frame(
      measure = m, cohort = measures[[m]], category_var = cv, category = as.character(lv),
      N = sum(k),
      cases = sum(Sv[[m]][k] == 1),
      prevalence = coef(p)[[1]],
      prev_lower = confint(p)[1],
      prev_upper = confint(p)[2],
      pop_millions = tot / 1e6,
      pop_lower = tot * exp(-h) / 1e6,
      pop_upper = tot * exp(h) / 1e6,
      df = df
    )
  }
  # Share of adults with the measure who fall below each treatment threshold.
  for (cv in c("asc5", "tot75")) {
    k <- Sv[[measures[[m]]]] & Sv[[m]] %in% 1
    p <- svyciprop(as.formula(paste0("~I(", cv, "==0)")), subset(S, k), method = "logit")
    burden[[length(burden) + 1]] <- data.frame(
      measure = m,
      cohort = measures[[m]],
      category_var = paste0("share_below_", cv),
      category = "below threshold",
      N = sum(k),
      cases = sum(Sv[[cv]][k] == 0),
      prevalence = coef(p)[[1]],
      prev_lower = confint(p)[1],
      prev_upper = confint(p)[2],
      pop_millions = NA, pop_lower = NA, pop_upper = NA, df = degf(subset(S, k))
    )
  }
}
wcsv(do.call(rbind, burden), "guideline_burden")

# A5-A7. Discordance groups
km5 <- function(time, death, w) { # weighted Kaplan-Meier cumulative mortality at 5 years
  k <- w > 0
  f <- survfit(Surv(time[k], death[k]) ~ 1, weights = w[k])
  1 - summary(f, times = 5, extend = TRUE)$surv
}
Rv <- R$variables
rw <- weights(R, type = "analysis")
fw <- weights(R, type = "sampling")
discord <- function(label, dom, fx, thr) {
  k <- Rv[[dom]] & !is.na(Rv[[fx]]) & !is.na(Rv[[thr]])
  grp <- list(
    below_preserved = c(0, 0),
    below_low = c(0, 1),
    above_preserved = c(1, 0),
    above_low = c(1, 1)
  )
  idx <- lapply(grp, function(g) k & Rv[[thr]] == g[1] & Rv[[fx]] == g[2])
  est <- function(w) {
    r <- vapply(idx, function(i) km5(Rv$time[i], Rv$death[i], w[i]), 0)
    c(r,
      rd_below = r[["below_low"]] - r[["below_preserved"]],
      rd_above = r[["above_low"]] - r[["above_preserved"]],
      additive_interaction = (r[["above_low"]] - r[["above_preserved"]]) - (r[["below_low"]] - r[["below_preserved"]]),
      discordance = r[["below_low"]] - r[["above_preserved"]]
    )
  }
  pt0 <- est(fw)
  reps <- t(vapply(seq_len(ncol(rw)), function(j) est(rw[, j]), numeric(length(pt0))))
  df <- degf(subset(S, S$variables[[dom]] & !is.na(S$variables[[fx]]) & !is.na(S$variables[[thr]])))
  v <- R$scale * colSums(sweep(sweep(reps, 2, pt0, "-")^2, 1, R$rscales, "*"))
  q <- qt(.975, df)
  cnt <- t(vapply(idx, function(i) c(
    N = sum(i), deaths = sum(Rv$death[i]), cvd_deaths = sum(Rv$death_cvd[i], na.rm = TRUE),
    wt_mean_age = weighted.mean(Rv$age[i], fw[i]), pop_millions = sum(fw[i]) / 1e6
  ), numeric(5)))
  out <- data.frame(
    analysis = label, cohort = dom, function_measure = fx, threshold = thr, estimand = names(pt0),
    estimate = pt0,
    SE = sqrt(v),
    lower = pt0 - q * sqrt(v),
    upper = pt0 + q * sqrt(v),
    df = df,
    row.names = NULL
  )
  out <- merge(out, data.frame(estimand = rownames(cnt), cnt), all.x = TRUE, sort = FALSE)
  out$other_deaths <- out$deaths - out$cvd_deaths
  # Mortality is reported only for groups with at least 10 deaths.
  small <- !is.na(out$deaths) & out$deaths < 10
  out$reportable <- !small
  out[order(match(out$estimand, names(pt0))), ]
}
dg <- rbind(
  discord("grip, ASCVD 5%", "core", "low_grip", "asc5"),
  discord("grip, total CVD 7.5%", "core", "low_grip", "tot75"),
  discord("reported function, ASCVD 5%", "older", "anyfx", "asc5"),
  discord("reported function, total CVD 7.5%", "older", "anyfx", "tot75")
)
# A contrast is reportable only if every group it uses is reportable.
uses <- list(
  rd_below = c("below_low", "below_preserved"), rd_above = c("above_low", "above_preserved"),
  additive_interaction = c("below_preserved", "below_low", "above_preserved", "above_low"),
  discordance = c("below_low", "above_preserved")
)
for (a in unique(dg$analysis)) for (e in names(uses)) {
  ok <- all(dg$reportable[dg$analysis == a & dg$estimand %in% uses[[e]]])
  dg$reportable[dg$analysis == a & dg$estimand == e] <- ok
}
wcsv(dg, "guideline_discordance")

# Supporting adjusted Cox models: low grip x threshold, primary covariates.
ci_terms <- function(m, des, terms) {
  df <- degf(des) - length(coef(m)) + 1
  V <- vcov(m)
  b <- coef(m)
  q <- qt(.975, df)
  do.call(rbind, lapply(names(terms), function(nm) {
    L <- terms[[nm]]
    e <- sum(L * b[names(L)])
    se <- sqrt(drop(t(L) %*% V[names(L), names(L)] %*% L))
    data.frame(
      term = nm,
      HR = exp(e),
      lower = exp(e - q * se),
      upper = exp(e + q * se),
      P = 2 * pt(-abs(e / se), df),
      df = df
    )
  }))
}
adj <- list()
for (thr in c("asc5", "tot75")) {
  D <- subset(S, core)
  m <- svycoxph(as.formula(paste0("Surv(time,death) ~ low_grip*", thr, "+", base)), design = D)
  it <- paste0("low_grip:", thr)
  r <- ci_terms(
    m,
    D,
    list(
      low_grip_below = c(low_grip = 1),
      low_grip_above = setNames(c(1, 1), c("low_grip", it)),
      interaction = setNames(1, it)
    )
  )
  adj[[thr]] <- cbind(threshold = thr, N = m$n, deaths = m$nevent, r)
}
wcsv(do.call(rbind, adj), "guideline_discordance_adjusted_hr")

# B, C. Alternative adjustment bases
addb <- function(x, var, nm) {
  b <- ns(qlogis(x[[var]][x$core]), df = 3)
  p <- predict(b, qlogis(x[[var]]))
  for (j in 1:3) x[[paste0(nm, j)]] <- p[, j]
  x
}
d <- addb(d, "r_neut", "kn")
d <- addb(d, "r_asc", "as")
d <- addb(d, "r_hf", "hf")
eb <- ns(d$egfr[d$core], df = 3)
ep <- predict(eb, d$egfr)
for (j in 1:3) d[[paste0("eg", j)]] <- ep[, j]
S <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~w,
  nest = TRUE,
  data = d[d$design_ok, ]
)
rest <- "a1+a2+a3+sex_factor+cycle+height+bmi"
bases <- c(
  "Kidney-neutral PREVENT (eGFR set to 90)" = paste("kn1+kn2+kn3+", rest),
  "Kidney-neutral PREVENT plus eGFR spline" = paste("kn1+kn2+kn3+eg1+eg2+eg3+", rest),
  "PREVENT-ASCVD in place of total CVD" = paste("as1+as2+as3+", rest),
  "PREVENT-HF in place of total CVD" = paste("hf1+hf2+hf3+", rest)
)
alt <- list()
for (nm in names(bases)) {
  for (spec in list(
    list(c = "core", rhs = "lower_z", t = "lower_z", lab = "Primary"),
    list(
      c = "mobility_cohort",
      rhs = "lower_z+mobility",
      t = c("lower_z", "mobility"),
      lab = "Mobility D"
    )
  )) {
    D <- subset(S, get(spec$c))
    m <- svycoxph(as.formula(paste("Surv(time,death) ~", spec$rhs, "+", bases[[nm]])), design = D)
    df <- degf(D) - length(coef(m)) + 1
    se <- sqrt(diag(vcov(m)))
    q <- qt(.975, df)
    for (t in spec$t) alt[[length(alt) + 1]] <- data.frame(
      adjustment = nm, model = spec$lab, term = t, N = m$n, deaths = m$nevent,
      HR = exp(coef(m)[[t]]),
      lower = exp(coef(m)[[t]] - q * se[[t]]),
      upper = exp(coef(m)[[t]] + q * se[[t]]),
      P = 2 * pt(-abs(coef(m)[[t]] / se[[t]]), df), df = df
    )
  }
}
wcsv(do.call(rbind, alt), "alternative_adjustment_sensitivity")

# B4. Creatinine, eGFR, and the kidney term by sex-specific grip quartile (weighted quartiles, primary cohort).
Sv <- S$variables
C <- subset(S, core)
qs <- lapply(
  1:2,
  function(sx) as.numeric(coef(svyquantile(~grip, subset(C, sex == sx), c(.25, .5, .75))))
)
qcut <- function(g, sx) findInterval(g, qs[[sx]], left.open = TRUE) + 1
Sv$grip_q <- NA_integer_
for (sx in 1:2) {
  i <- Sv$core & Sv$sex == sx
  Sv$grip_q[i] <- qcut(Sv$grip[i], sx)
}
S$variables$grip_q <- Sv$grip_q
S$variables$kterm_pp <- 100 * (Sv$r_tot - Sv$r_neut)
C <- subset(S, core)
b4 <- do.call(rbind, lapply(1:4, function(g) {
  D <- subset(C, grip_q == g)
  data.frame(
    grip_quartile = g,
    N = nrow(D$variables),
    creatinine_mg_dl = coef(svymean(~scr, D))[[1]],
    egfr = coef(svymean(~egfr, D))[[1]],
    egfr_below60_pct = 100 * coef(svymean(~ I(egfr < 60), D))[[2]],
    kidney_term_pp = coef(svymean(~kterm_pp, D))[[1]],
    prevent_total_pct = 100 * coef(svymean(~r_tot, D))[[1]],
    prevent_neutral_pct = 100 * coef(svymean(~r_neut, D))[[1]],
    age = coef(svymean(~age, D))[[1]]
  )
}))
b4$grip_quartile_label <- c("Q1 (weakest)", "Q2", "Q3", "Q4 (strongest)")
wcsv(b4, "kidney_term_by_grip_quartile")

# D. Selection at each exclusion gate
dd <- d[d$design_ok, ]
tgt <- dd$clinical
gates <- list(
  "Incomplete PREVENT inputs" = tgt & !dd$complete_score,
  "PREVENT inputs outside equation ranges" = tgt & dd$complete_score & !dd$range_ok,
  "No eligible mortality linkage or no height" = tgt & dd$complete_score & dd$range_ok & !dd$pregrip,
  "Grip not tested" = dd$pregrip & dd$grip_status == "not tested mixed reasons",
  "Valid maximal-effort grip from 1 hand only" = dd$pregrip & dd$grip_status == "one hand valid maximal effort",
  "No qualifying maximal-effort trial" = dd$pregrip & dd$grip_status == "tested but no qualifying effort",
  "Included (primary cohort)" = dd$core
)
stopifnot(sum(sapply(gates, sum)) == sum(tgt))
gates <- c(
  gates[-length(gates)],
  list("All excluded candidates" = tgt & !dd$core),
  gates[length(gates)]
)
sel <- do.call(rbind, lapply(names(gates), function(g) {
  k <- gates[[g]]
  L <- k & dd$linked
  wm <- function(x) if (sum(k & !is.na(x))) weighted.mean(
    x[k & !is.na(x)],
    dd$w[k & !is.na(x)]
  ) else NA
  data.frame(
    group = g,
    N = sum(k),
    linked = sum(L),
    deaths = sum(dd$death[L] == 1),
    crude_mortality_pct = 100 * mean(dd$death[L] == 1),
    wt_mean_age = wm(dd$age),
    wt_female_pct = 100 * wm(as.numeric(dd$sex == 2)),
    wt_mean_bmi = wm(dd$bmi),
    wt_diabetes_pct = 100 * wm(dd$diabetes), median_followup_years = median(dd$time[L])
  )
}))
wcsv(sel, "selection_by_gate")
writeLines(capture.output(sessionInfo()), file.path(P$logs, "sessionInfo_guideline_categories.txt"))
message("Guideline-category analyses written to ", OUT)
