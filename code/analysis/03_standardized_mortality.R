# Risk-stratified and cohort-specific spline models, selection summaries, and standardized 5-year
# mortality at the 25th vs 75th grip-strength percentile with full jackknife refits.
# Run through code/04_sensitivity_analysis.R (working directory: the build directory).

suppressPackageStartupMessages({
  library(survey)
  library(survival)
  library(splines)
  library(dplyr)
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
A <- readRDS(file.path(B, "03_Derived/analysis_models.rds"))
d <- A$data
S <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~w,
  nest = TRUE,
  data = d[d$design_ok, ]
)
C <- subset(S, core)
M <- subset(S, mobility_cohort)
base <- "r1+r2+r3+a1+a2+a3+sex_factor+cycle+height+bmi"
newmodels <- list()
results <- list()
diags <- list()
run <- function(label, dd, f) {
  m <- svycoxph(as.formula(f), dd, x = TRUE, y = TRUE, model = TRUE)
  df <- degf(dd) - length(coef(m)) + 1
  se <- sqrt(diag(vcov(m)))
  q <- qt(.975, df)
  t <- intersect(c("lower_z", "mobility"), names(coef(m)))
  results[[label]] <<- data.frame(
    model = label,
    term = t,
    beta = coef(m)[t],
    SE = se[t],
    HR = exp(coef(m)[t]),
    lower = exp(coef(m)[t] - q * se[t]),
    upper = exp(coef(m)[t] + q * se[t]),
    P = 2 * pt(-abs(coef(m)[t] / se[t]), df),
    df = df,
    N = nrow(m$model),
    deaths = m$nevent,
    row.names = NULL
  )
  diags[[label]] <<- data.frame(model = label, condition_number = kappa(scale(m$x), exact = TRUE))
  newmodels[[label]] <<- m
  invisible(m)
}
for (cc in c("core", "mobility_cohort")) {
  i <- d[[cc]]
  bk <- unique(quantile(d$risk[i], c(0, .25, .5, .75, 1)))
  d$risk_band <- cut(d$risk, breaks = c(-Inf, bk[2:4], Inf), labels = FALSE, include.lowest = TRUE)
  SS <- svydesign(
    ids = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~w,
    nest = TRUE,
    data = d[d$design_ok, ]
  )
  D <- subset(SS, get(cc))
  run(
    paste("Post-diagnostic: risk-stratified", cc),
    D,
    paste(
      "Surv(time,death)~lower_z+",
      if (cc == "mobility_cohort") "mobility+" else "",
      base,
      "+strata(risk_band)"
    )
  )
}
dm <- M$variables
rr <- ns(dm$logitrisk, df = 3)
aa <- ns(dm$age, df = 3)
for (j in 1:3) {
  M$variables[[paste0("r", j)]] <- rr[, j]
  M$variables[[paste0("a", j)]] <- aa[, j]
}
run(
  "Post-diagnostic: mobility-specific spline knots",
  M,
  paste("Surv(time,death)~lower_z+mobility+", base)
)
out(bind_rows(results), "post_diagnostic_estimates")
out(bind_rows(diags), "post_diagnostic_conditioning")
saveRDS(newmodels, file.path(B, "03_Derived/post_diagnostic_models.rds"))
# Selection assessment: missingness and availability, not an unweighted comparison P test.
d$selection_group <- ifelse(d$core,
  "Included core",
  ifelse(d$clinical, "Excluded clinical candidates", NA)
)
sel <- list()
for (g in c("Included core", "Excluded clinical candidates")) {
  xx <- d[d$design_ok & d$selection_group %in% g, ]
  des <- subset(S, SEQN %in% xx$SEQN)
  for (v in c("age", "bmi", "sbp", "egfr", "diabetes", "smoking")) {
    z <- svymean(as.formula(paste0("~", v)), des, na.rm = TRUE)
    sel[[paste(g, v)]] <- data.frame(
      group = g,
      variable = v,
      N = nrow(xx),
      available = sum(!is.na(xx[[v]])),
      estimate = coef(z)[1],
      SE = SE(z)[1]
    )
  }
}
out(bind_rows(sel), "selection_comparison")
# Distinguish questionnaire routing from missing responses in the age-eligible mobility cohort.
out(
  data.frame(group = c("Core younger than 60: outside mobility target", "Core 60-79 classified", "Core 60-79 unclassified"), N = c(sum(d$core & d$age < 60), sum(d$mobility_cohort), sum(d$core & d$age >= 60 & is.na(d$mobility)))),
  "mobility_routing"
)
# Weighted contrast values are fixed across replicates, as documented.
c <- C$variables
quant <- bind_rows(lapply(1:2, function(sx) {
  v <- svyquantile(~grip, subset(C, sex == sx), c(.25, .75), ci = FALSE)
  data.frame(sex = sx, q25 = coef(v)[1], q75 = coef(v)[2])
}))
out(quant, "risk_strength_contrasts")
c$cell <- interaction(c$sex,
  cut(c$age, c(49, 59, 69, 79)),
  cut(c$bmi, c(18.49, 25, 30, 39.91), right = FALSE),
  drop = TRUE
)
supports <- bind_rows(lapply(split(c, c$cell), function(z) {
  q <- quant[quant$sex == z$sex[1], ]
  lo <- quantile(z$grip, .05)
  hi <- quantile(z$grip, .95)
  data.frame(
    cell = as.character(z$cell[1]),
    N = nrow(z),
    q05 = lo,
    q95 = hi,
    support = nrow(z) >= 30 & q$q25 >= lo & q$q75 <= hi
  )
}))
c$supported <- supports$support[match(c$cell, supports$cell)]
stopifnot(sum(c$supported) > 0, mean(c$time >= 5) > .90, min(c$time[c$death == 0]) >= 4.8)
out(supports, "risk_support_cells")
out(
  data.frame(core_N = nrow(c), supported_N = sum(c$supported), supported_deaths = sum(c$death[c$supported]), supported_weight_fraction = sum(c$w[c$supported]) / sum(c$w), at_risk_5 = sum(c$time >= 5)),
  "risk_support_summary"
)
C$variables <- c
R <- as.svrepdesign(S, type = "JKn", mse = TRUE)
R <- subset(R, core)
rw <- weights(R, type = "analysis")
fw <- weights(R, type = "sampling")
f <- as.formula(paste("Surv(time,death)~lower_z+", base))
sc <- read.csv(file.path(B, "05_Outputs/strength_scaling.csv"))
predrisk <- function(w) {
  idx <- which(w > 0)
  z <- c[idx, ]
  z$ww <- w[idx]
  fit <- coxph(f, data = z, weights = ww, ties = "efron", x = TRUE, model = TRUE)
  stopifnot(all(is.finite(coef(fit))))
  bh <- basehaz(fit, centered = TRUE)
  h <- tail(bh$hazard[bh$time <= 5], 1)
  avg <- function(col) {
    zz <- z[z$supported, ]
    zz$grip <- quant[[col]][match(zz$sex, quant$sex)]
    zz$lower_z <- (sc$mean_kg[match(zz$sex, sc$sex)] - zz$grip) / sc$sd_kg[match(zz$sex, sc$sex)]
    lp <- predict(fit, newdata = zz, type = "lp", reference = "sample")
    r <- -expm1(-h * exp(lp))
    weighted.mean(r, zz$ww)
  }
  p25 <- avg("q25")
  p75 <- avg("q75")
  c(q25_risk = p25, q75_risk = p75, difference = p25 - p75, strength_beta = coef(fit)["lower_z"])
}
p <- predrisk(fw)
stopifnot(abs(p[4] - coef(A$models[["Primary: continuous strength"]])["lower_z"]) < 1e-7)
rep <- matrix(NA_real_, ncol(rw), length(p))
logs <- list()
for (k in seq_len(ncol(rw))) {
  w <- character()
  err <- ""
  v <- tryCatch(withCallingHandlers(predrisk(rw[, k]), warning = function(e) {
    w <<- c(w, conditionMessage(e))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    err <<- conditionMessage(e)
    rep(NA_real_, length(p))
  })
  rep[k, ] <- v
  logs[[k]] <- data.frame(
    replicate = k,
    success = all(is.finite(v)),
    warnings = paste(w, collapse = " | "),
    error = err
  )
}
out(bind_rows(logs), "jackknife_diagnostics")
stopifnot(all(is.finite(rep)))
vv <- R$scale * colSums(sweep((sweep(rep, 2, p, "-"))^2, 1, R$rscales, "*"))
df <- degf(C) - length(coef(A$models[["Primary: continuous strength"]])) + 1
rr <- data.frame(
  estimand = names(p)[1:3],
  estimate = p[1:3],
  SE = sqrt(vv[1:3]),
  lower = p[1:3] - qt(.975, df) * sqrt(vv[1:3]),
  upper = p[1:3] + qt(.975, df) * sqrt(vv[1:3]),
  df = df,
  replicates = nrow(rep)
)
out(rr, "standardized_5yr_risks")
out(data.frame(replicate = seq_len(nrow(rep)), rep), "jackknife_estimates")
out(
  data.frame(scale = R$scale, replicates = ncol(rw), min_rscale = min(R$rscales), max_rscale = max(R$rscales), mse = R$mse, inference_df = df),
  "jackknife_scaling"
)
# Direct survival-prediction cross-check for the full-sample baseline-hazard calculation.
cx <- coxph(f, data = c, weights = w, x = TRUE, model = TRUE, ties = "efron")
nn <- c[1, , drop = FALSE]
ss <- summary(survfit(cx, newdata = nn), times = 5)$surv
hh <- tail(basehaz(cx, centered = TRUE)$hazard[basehaz(cx, centered = TRUE)$time <= 5], 1)
pp <- exp(-hh * exp(predict(cx, newdata = nn, type = "lp", reference = "sample")))
stopifnot(abs(ss - pp) < 1e-10)
out(
  data.frame(test = "baseline hazard versus survfit at five years", absolute_difference = abs(ss - pp)),
  "survival_prediction_validation"
)
capture.output(sessionInfo(), file = file.path(B, "08_Audit/sessionInfo_validation.txt"))
print(bind_rows(results))
print(rr)
print(read.csv(file.path(B, "05_Outputs/risk_support_summary.csv")))
