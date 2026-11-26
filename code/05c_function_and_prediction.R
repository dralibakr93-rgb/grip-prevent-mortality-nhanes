# Stage 5c. Reported function, prior hospitalization, and mortality prediction:
#   A. ADL and IADL limitation at ages 60-79 years: survey logistic models and standardized
#      prevalence (Figure 3C);
#   B. overnight hospitalization in the prior year: prevalence and sensitivity models;
#   C. incremental mortality prediction: weighted C statistic, 5-year AUC, NRI, IDI, and calibration
#      slope, with jackknife CIs and survey-bootstrap optimism correction (Figures 2B and S6).

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
d <- add_function_hospital_vars(d)
d$older <- d$core & d$age >= 60
stopifnot(
  sum(d$older & !is.na(d$adl)) == 1523,
  sum(d$older & d$adl %in% 1) == 217
) # expected analytic counts

design <- function(x) svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~w,
  nest = TRUE,
  data = x[x$design_ok, ]
)
S <- design(d)
base <- "r1+r2+r3+a1+a2+a3+sex_factor+cycle+height+bmi"
sc <- res("strength_scaling")
quant <- res("risk_strength_contrasts")
z_at <- function(
  x,
  q
) (sc$mean_kg[match(x$sex, sc$sex)] - quant[[q]][match(x$sex, quant$sex)]) / sc$sd_kg[match(
  x$sex,
  sc$sex
)]
R <- as.svrepdesign(S, type = "JKn", mse = TRUE)
jk_ci <- function(Rs, point, reps, df) {
  v <- Rs$scale * colSums(sweep(sweep(reps, 2, point, "-")^2, 1, Rs$rscales, "*"))
  data.frame(
    estimate = point,
    SE = sqrt(v),
    lower = point - qt(0.975, df) * sqrt(v),
    upper = point + qt(0.975, df) * sqrt(v),
    df = df
  )
}

# A. ADL / IADL (cross-sectional)
fd <- list()
fs <- list()
for (o in c("adl", "iadl")) {
  D <- subset(S, older & !is.na(get(o)))
  f <- as.formula(paste(o, "~ lower_z +", base))
  m <- svyglm(f, design = D, family = quasibinomial())
  df <- degf(D) - length(coef(m)) + 1
  b <- coef(m)[["lower_z"]]
  se <- SE(m)[["lower_z"]]
  q <- qt(0.975, df)
  pv <- svymean(as.formula(paste0("~", o)), D)
  fd[[o]] <- data.frame(
    outcome = o, N = nrow(D$variables), cases = sum(D$variables[[o]]),
    weighted_prevalence = coef(pv)[[1]],
    OR = exp(b),
    lower = exp(b - q * se),
    upper = exp(b + q * se),
    P = 2 * pt(-abs(b / se), df), beta = b, SE = se, df = df
  )
  # Standardized prevalence at the 25th vs 75th sex-specific percentiles, JKn refits.
  Rs <- subset(R, older & !is.na(get(o)))
  dat <- Rs$variables
  fw <- weights(Rs, type = "sampling")
  rw <- weights(Rs, type = "analysis")
  calc <- function(w) {
    k <- w > 0
    z <- dat[k, ]
    z$ww <- w[k] / mean(w[k]) # rescaled as in svyglm (glm diverges otherwise)
    g <- glm(f, data = z, weights = ww, family = quasibinomial())
    stopifnot(g$converged)
    pr <- vapply(c("q25", "q75"), function(qn) {
      zz <- z
      zz$lower_z <- z_at(zz, qn)
      weighted.mean(predict(g, newdata = zz, type = "response"), zz$ww)
    }, 0)
    c(pr, pr[1] - pr[2], coef(g)[["lower_z"]])
  }
  pt0 <- calc(fw)
  stopifnot(abs(pt0[4] - b) < 1e-6)
  reps <- t(vapply(seq_len(ncol(rw)), function(j) calc(rw[, j]), numeric(4)))
  ci <- jk_ci(Rs, pt0[1:3], reps[, 1:3, drop = FALSE], df)
  fs[[o]] <- cbind(
    outcome = o,
    estimand = c("q25_prevalence", "q75_prevalence", "difference"),
    ci,
    replicates = nrow(reps)
  )
}
wcsv(do.call(rbind, fd), "functional_disability_estimates")
wcsv(do.call(rbind, fs), "functional_disability_standardized")

# B. Prior hospitalization
hp <- list()
for (k in c("core", "mobility_cohort")) {
  D <- subset(S, get(k) & !is.na(prior_hosp))
  v <- D$variables
  hp[[k]] <- data.frame(
    cohort = k, N_cohort = sum(d[[k]]), known = nrow(v), unknown = sum(d[[k]]) - nrow(v),
    hospitalized = sum(v$prior_hosp), weighted_prevalence = coef(svymean(~prior_hosp, D))[[1]],
    deaths_hospitalized = sum(v$death[v$prior_hosp == 1]),
    deaths_not = sum(v$death[v$prior_hosp == 0])
  )
}
wcsv(do.call(rbind, hp), "prior_hospitalization_prevalence")
fitcox <- function(label, D, rhs, terms) {
  m <- svycoxph(as.formula(paste("Surv(time,death)~", rhs)), design = D)
  df <- degf(D) - length(coef(m)) + 1
  se <- sqrt(diag(vcov(m)))
  q <- qt(0.975, df)
  do.call(rbind, lapply(terms, function(t) data.frame(
    model = label, term = t, N = m$n, deaths = m$nevent, HR = exp(coef(m)[[t]]),
    lower = exp(coef(m)[[t]] - q * se[[t]]),
    upper = exp(coef(m)[[t]] + q * se[[t]]),
    P = 2 * pt(-abs(coef(m)[[t]] / se[[t]]), df),
    df = df
  )))
}
hs <- rbind(
  fitcox("Primary: excluding prior-year hospitalization", subset(S, core & prior_hosp %in% 0), paste("lower_z+", base), "lower_z"),
  fitcox("Mobility D: excluding prior-year hospitalization", subset(S, mobility_cohort & prior_hosp %in% 0), paste("lower_z+mobility+", base), c("lower_z", "mobility")),
  fitcox(
    "Primary: adjusted for prior-year hospitalization",
    subset(S, core & !is.na(prior_hosp)),
    paste("lower_z+prior_hosp+", base),
    c("lower_z", "prior_hosp")
  )
)
wcsv(hs, "prior_hospitalization_sensitivity")

# C. Incremental mortality prediction
wauc <- function(score, case, ctrl, w) { # weighted AUC (Mann-Whitney) for cases vs controls
  sc1 <- score[case]
  w1 <- w[case]
  sc0 <- score[ctrl]
  w0 <- w[ctrl]
  gt <- outer(sc1, sc0, ">") + 0.5 * outer(sc1, sc0, "==")
  sum(w1 * (gt %*% w0)) / (sum(w1) * sum(w0))
}
perf <- function(fits, dat, w) { # metrics of fitted models evaluated on data with weights w
  k <- w > 0
  z <- dat[k, ]
  ww <- w[k]
  case <- z$death == 1 & z$time <= 5
  ctrl <- !case & z$time >= 5 # alive at 5 y; only those censored before 5 y excluded
  lps <- lapply(fits, function(f) predict(f, newdata = z, type = "lp", reference = "sample"))
  r5 <- mapply(function(f, lp) {
    bh <- basehaz(f, centered = TRUE)
    H <- tail(bh$hazard[bh$time <= 5], 1)
    -expm1(-H * exp(lp))
  }, fits, lps, SIMPLIFY = FALSE)
  out <- c()
  for (nm in names(fits)) {
    cc <- concordance(Surv(z$time, z$death) ~ lps[[nm]], weights = ww, reverse = TRUE)$concordance
    out[paste0("C:", nm)] <- cc
    out[paste0("AUC5:", nm)] <- wauc(lps[[nm]], case, ctrl, ww)
  }
  for (nm in setdiff(names(fits), "M0")) {
    out[paste0("dC:", nm)] <- out[paste0("C:", nm)] - out["C:M0"]
    out[paste0("dAUC5:", nm)] <- out[paste0("AUC5:", nm)] - out["AUC5:M0"]
    up <- r5[[nm]] > r5[["M0"]]
    dn <- r5[[nm]] < r5[["M0"]]
    wm <- function(x, s) sum(ww[s] * x[s]) / sum(ww[s])
    out[paste0("NRI:", nm)] <- (wm(up, case) - wm(dn, case)) + (wm(dn, ctrl) - wm(up, ctrl))
    out[paste0("IDI:", nm)] <- (wm(r5[[nm]] - r5[["M0"]], case)) - (wm(r5[[nm]] - r5[["M0"]], ctrl))
  }
  out
}
slopes <- function(fits, dat, w_eval) { # calibration slope of each model's lp on evaluation weights
  k <- w_eval > 0
  z <- dat[k, ]
  ww <- w_eval[k] / mean(w_eval[k])
  vapply(fits, function(f) {
    z$lp <- predict(f, newdata = z, type = "lp", reference = "sample")
    coef(coxph(Surv(time, death) ~ lp, data = z, weights = ww))[["lp"]]
  }, 0)
}
fitset <- function(specs, dat, w) {
  k <- w > 0
  z <- dat[k, ]
  z$ww <- w[k] / mean(w[k])
  lapply(
    specs,
    function(rhs) coxph(as.formula(paste("Surv(time,death) ~", rhs)),
      data = z,
      weights = ww,
      ties = "efron",
      x = FALSE
    )
  )
}

run_prediction <- function(label, domain, specs) {
  Rs <- subset(R, get(domain))
  dat <- Rs$variables
  fw <- weights(Rs, type = "sampling")
  rw <- weights(Rs, type = "analysis")
  fits <- fitset(specs, dat, fw)
  app <- perf(fits, dat, fw)
  reps <- t(vapply(
    seq_len(ncol(rw)),
    function(j) perf(fitset(specs, dat, rw[, j]), dat, rw[, j]),
    numeric(length(app))
  ))
  ci <- jk_ci(Rs, app, reps, 22) # t(22) reference distribution
  # Survey bootstrap optimism (Rao-Wu rescaling bootstrap on the full design, then domain).
  set.seed(20261008)
  B <- subset(as.svrepdesign(S, type = "subbootstrap", replicates = 200), get(domain))
  stopifnot(identical(B$variables$SEQN, dat$SEQN))
  bw <- weights(B, type = "analysis")
  opt <- matrix(NA_real_, ncol(bw), length(app))
  sl <- matrix(NA_real_, ncol(bw), length(specs))
  for (j in seq_len(ncol(bw))) {
    fb <- fitset(specs, dat, bw[, j])
    opt[j, ] <- perf(fb, dat, bw[, j]) - perf(fb, dat, fw)
    sl[j, ] <- slopes(fb, dat, fw)
  }
  res_tab <- data.frame(
    cohort = label, metric = sub(":.*", "", names(app)), model = sub(".*:", "", names(app)),
    apparent = app, ci[, c("SE", "lower", "upper", "df")], optimism = colMeans(opt),
    corrected = app - colMeans(opt),
    bootstrap_replicates = nrow(opt),
    jk_replicates = nrow(reps),
    row.names = NULL
  )
  cal <- data.frame(
    cohort = label,
    metric = "calibration_slope",
    model = names(specs),
    apparent = 1,
    SE = NA,
    lower = NA,
    upper = NA,
    df = NA,
    optimism = 1 - colMeans(sl),
    corrected = colMeans(sl),
    bootstrap_replicates = nrow(sl),
    jk_replicates = NA
  )
  rbind(res_tab, cal)
}
specs_core <- list(M0 = base, M1 = paste("lower_z+", base), PREVENT_only = "r1+r2+r3")
pc <- run_prediction("primary", "core", specs_core)
specs_mob <- list(
  M0 = base,
  Mobility = paste("mobility+", base),
  Grip = paste("lower_z+", base),
  Both = paste("lower_z+mobility+", base)
)
pm <- run_prediction("mobility", "mobility_cohort", specs_mob)
wcsv(rbind(pc, pm), "prediction_metrics")
# Check: weighted Cox coefficient equals the primary survey-Cox estimate.
fchk <- fitset(
  list(M1 = specs_core$M1),
  subset(R, core)$variables,
  weights(subset(R, core), type = "sampling")
)
stopifnot(abs(coef(fchk$M1)[["lower_z"]] - coef(A$models[["Primary: continuous strength"]])[["lower_z"]]) < 1e-6)
writeLines(capture.output(sessionInfo()), file.path(P$logs, "sessionInfo_function_prediction.txt"))
message("Function and prediction analyses written to ", OUT)
