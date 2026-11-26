# Additional sensitivity analyses: BMI range relaxed, grip-strength spline with nonlinearity tests,
# flexible age, 3 trials per hand, full PREVENT (UACR and HbA1c), education, income and cancer history,
# income imputation, and proportional-hazards tests by model term.
# Run through code/04_sensitivity_analysis.R (working directory: the build directory).

suppressPackageStartupMessages({
  library(survey)
  library(survival)
  library(splines)
  library(dplyr)
  library(haven)
  library(Hmisc)
  library(mice)
  library(mitools)
})
options(survey.lonely.psu = "fail")
set.seed(410031)
B <- "study"
V <- file.path(B, "additional")
for (p in c("Results", "Audit", "Derived", "Sources")) dir.create(file.path(V, p),
  recursive = TRUE,
  showWarnings = FALSE
)
out <- function(
  x,
  n
) write.csv(x,
  file.path(V, "Results", paste0(n, ".csv")),
  row.names = FALSE,
  na = ""
)
A <- readRDS(file.path(B, "03_Derived/analysis_models.rds"))
d <- A$data
sc <- read.csv(file.path(B, "05_Outputs/strength_scaling.csv"))
sf <- readRDS(file.path(B, "03_Derived/scoring_functions.rds"))
coefs <- sf$coefs
terms <- sf$terms
score_any <- function(x) {
  z <- terms(x)
  v <- rep(NA_real_, nrow(x))
  for (s in 1:2) {
    i <- which(x$sex == s)
    v[i] <- plogis(z[i, ,
      drop = FALSE
    ] %*% coefs[[if (s == 1) "male_total_cvd" else "female_total_cvd"]])
  }
  v
}
design <- function(x) svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~w,
  nest = TRUE,
  data = x[x$design_ok, ]
)
base <- "r1+r2+r3+a1+a2+a3+sex_factor+cycle+height+bmi"
models <- list()
est <- list()
diag <- list()
fit <- function(label, des, rhs) {
  ww <- character()
  m <- withCallingHandlers(svycoxph(as.formula(paste("Surv(time,death)~", rhs)), des, x = TRUE, y = TRUE, model = TRUE),
    warning = function(w) {
      ww <<- c(ww, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  df <- degf(des) - length(coef(m)) + 1
  stopifnot(df > 5, all(is.finite(coef(m))), all(is.finite(vcov(m))))
  se <- sqrt(diag(vcov(m)))
  tt <- intersect(names(coef(m)), c("lower_z", "mobility", "g1", "g2"))
  q <- qt(.975, df)
  est[[label]] <<- data.frame(
    model = label,
    term = tt,
    HR = exp(coef(m)[tt]),
    lower = exp(coef(m)[tt] - q * se[tt]),
    upper = exp(coef(m)[tt] + q * se[tt]),
    P = 2 * pt(-abs(coef(m)[tt] / se[tt]), df),
    beta = coef(m)[tt],
    SE = se[tt],
    df = df,
    N = nrow(m$model),
    deaths = m$nevent,
    row.names = NULL
  )
  diag[[label]] <<- data.frame(
    model = label,
    N = nrow(m$model),
    deaths = m$nevent,
    parameters = length(coef(m)),
    design_df = degf(des),
    inference_df = df,
    warnings = paste(ww, collapse = " | "),
    min_cov_eigen = min(eigen(vcov(m), symmetric = TRUE, only.values = TRUE)$values)
  )
  models[[label]] <<- m
  m
}
S <- design(d)
C <- subset(S, core)
# Audit inputs individually and sequentially; do not infer missing eligibility.
lims <- list(
  tc = c(130, 320),
  hdl = c(20, 100),
  sbp = c(90, 200),
  egfr = c(15, 140),
  bmi = c(18.5, 39.9)
)
out(bind_rows(lapply(names(lims), function(v) {
  x <- d[[v]][d$clinical]
  l <- lims[[v]]
  data.frame(
    variable = v,
    N = length(x),
    missing = sum(is.na(x)),
    below = sum(x < l[1], na.rm = TRUE),
    above = sum(x > l[2], na.rm = TRUE)
  )
})), "range_exclusions")
nonbmi <- with(
  d,
  tc >= 130 & tc <= 320 & hdl >= 20 & hdl <= 100 & sbp >= 90 & sbp <= 200 & egfr >= 15 & egfr <= 140
)
d$bmi_relaxed <- d$clinical & d$complete_score & nonbmi & d$linked & !is.na(d$height) & d$both_hands
d$bmi_relaxed[is.na(d$bmi_relaxed)] <- FALSE
stopifnot(all(d$bmi_relaxed[d$core]), sum(d$bmi_relaxed) == 2879)
e <- d
e$risk <- score_any(e)
e$logitrisk <- qlogis(e$risk)
rb <- predict(A$risk_basis, e$logitrisk)
for (j in 1:3) e[[paste0("r", j)]] <- rb[, j]
fit("BMI range relaxed", subset(design(e), bmi_relaxed), paste("lower_z+", base))
out(
  data.frame(group = c("Original core", "BMI-only additions", "BMI-relaxed sensitivity", "Clinical candidates with missing score inputs and otherwise valid grip/linkage/height"), N = c(sum(d$core), sum(d$bmi_relaxed & !d$core), sum(d$bmi_relaxed), sum(d$clinical & !d$complete_score & d$both_hands & d$linked & !is.na(d$height))), deaths = c(sum(d$death[d$core]), sum(d$death[d$bmi_relaxed & !d$core]), sum(d$death[d$bmi_relaxed]), sum(d$death[d$clinical & !d$complete_score & d$both_hands & d$linked & !is.na(d$height)]))),
  "selection_counts"
)
# Three total RCS knots; two columns including linear term.
kn <- as.numeric(coef(svyquantile(~lower_z, C, c(.1, .5, .9), ci = FALSE)))
g <- rcspline.eval(d$lower_z, knots = kn, inclx = TRUE)
d$g1 <- g[, 1]
d$g2 <- g[, 2]
S <- design(d)
C <- subset(S, core)
nl <- fit("Grip restricted cubic spline", C, paste("g1+g2+", base))
dfnl <- degf(C) - length(coef(nl)) + 1
tests <- lapply(
  list(overall = ~ g1 + g2, nonlinear = ~g2),
  function(f) regTermTest(nl, f, method = "Wald", df = dfnl)
)
out(bind_rows(lapply(names(tests), function(n) {
  z <- tests[[n]]
  data.frame(test = n, F = as.numeric(z$Ftest), num_df = z$df, den_df = z$ddf, P = as.numeric(z$p))
})), "nonlinearity_tests")
out(data.frame(probability = c(.1, .5, .9), lower_z = kn), "grip_spline_knots")
grid <- seq(quantile(d$lower_z[d$core], .01), quantile(d$lower_z[d$core], .99), length.out = 161)
G <- sweep(
  rcspline.eval(grid, knots = kn, inclx = TRUE),
  2,
  as.numeric(rcspline.eval(0, knots = kn, inclx = TRUE)),
  "-"
)
ii <- match(c("g1", "g2"), names(coef(nl)))
lp <- as.vector(G %*% coef(nl)[ii])
se <- sqrt(rowSums((G %*% vcov(nl)[ii, ii]) * G))
out(
  data.frame(lower_z = grid, strength_z = -grid, HR = exp(lp), lower = exp(lp - qt(.975, dfnl) * se), upper = exp(lp + qt(.975, dfnl) * se)),
  "grip_spline_curve"
)
# Individual clinical factors with the same flexible age adjustment as the primary model.
fit(
  "Clinical factors with matched flexible age",
  C,
  "lower_z+a1+a2+a3+sex_factor+cycle+height+bmi+sbp+tc+hdl+egfr+diabetes+smoking+bp_tx+statin"
)
fit(
  "Three maximal-effort trials per hand",
  subset(S, core & nhand1 == 3 & nhand2 == 3),
  paste("lower_z+", base)
)
# Full PREVENT with UACR and HbA1c, missing SDI, on a common sample.
.libPaths(c("R_library", .libPaths()))
ns <- asNamespace("preventr")
fullcoef <- get("full_10yr", ns)
write.csv(fullcoef, file.path(V, "Sources/preventr_full_10yr_coefficients.csv"), row.names = FALSE)
capture.output(get("prep_terms", ns), file = file.path(V, "Sources/preventr_prep_terms.txt"))
for (yy in c(2011, 2013)) {
  suffix <- if (yy == 2011) "G" else "H"
  p <- file.path("study/02_Raw/urine", yy, paste0("ALB_CR_", suffix, ".xpt"))
  file.copy(p, file.path(V, "Sources", basename(p)), overwrite = TRUE)
  u <- as.data.frame(read_xpt(p))
  if (!"uacr" %in% names(d)) d$uacr <- NA_real_
  jj <- match(d$SEQN, u$SEQN)
  ok <- !is.na(jj)
  d$uacr[ok] <- u$URXUMA[jj[ok]] / u$URXUCR[jj[ok]] * 100
  if ("URDACT" %in% names(u)) stopifnot(max(abs(u$URXUMA / u$URXUCR * 100 - u$URDACT),
    na.rm = TRUE
  ) < .02)
}
d$enhanced_ok <- d$core & !is.na(d$uacr) & d$uacr >= .1 & d$uacr <= 25000 & !is.na(d$LBXGH) & d$LBXGH >= 4.5 & d$LBXGH <= 15
d$risk_full <- NA_real_
i <- which(d$enhanced_ok)
x <- d[i, ]
Z <- cbind(
  terms(x)[, 1:22],
  0,
  0,
  1,
  log(x$uacr),
  0,
  ifelse(x$diabetes == 1, x$LBXGH - 5.3, 0),
  ifelse(x$diabetes == 0, x$LBXGH - 5.3, 0),
  0,
  1
)
for (s in 1:2) {
  j <- which(x$sex == s)
  d$risk_full[i[j]] <- plogis(Z[j, ,
    drop = FALSE
  ] %*% fullcoef[[if (s == 1) "male_total_cvd" else "female_total_cvd"]])
}
set.seed(410031)
ch <- sample(i, 100)
err <- sapply(ch, function(j) {
  a <- d[j, ]
  pd <- data.frame(
    age = a$age,
    total_c = a$tc,
    hdl_c = a$hdl,
    sbp = a$sbp,
    dm = a$diabetes,
    smoking = a$smoking,
    bmi = a$bmi,
    egfr = a$egfr,
    bp_tx = a$bp_tx,
    statin = a$statin,
    uacr = a$uacr,
    hba1c = a$LBXGH
  )
  pt <- get("prep_terms", ns)(pd, model = "full", units = "mg/dL")
  pt <- pt[names(pt) != "age_squared"]
  p <- plogis(sum(pt * fullcoef[[if (a$sex == 1) "male_total_cvd" else "female_total_cvd"]]))
  abs(p - a$risk_full)
})
stopifnot(max(err) < 1e-12)
out(
  data.frame(test = "Vectorized full-score terms versus preventr term preparation", N = 100, maximum_difference = max(err), package_version = as.character(packageVersion("preventr"))),
  "enhanced_score_validation"
)
S <- design(d)
E <- subset(S, enhanced_ok)
fit("Base PREVENT common biomarker sample", E, paste("lower_z+", base))
ek <- ns(qlogis(d$risk_full[d$enhanced_ok]), df = 3)
er <- predict(ek, qlogis(d$risk_full))
for (j in 1:3) d[[paste0("e", j)]] <- er[, j]
fit(
  "Full PREVENT with UACR and HbA1c",
  subset(design(d), enhanced_ok),
  paste("lower_z+", gsub("r([123])", "e\\1", base))
)
out(
  data.frame(N = sum(d$enhanced_ok), deaths = sum(d$death[d$enhanced_ok]), weighted_base = weighted.mean(d$risk[i], d$w[i]), weighted_full = weighted.mean(d$risk_full[i], d$w[i]), missing_SDI = TRUE),
  "enhanced_score_summary"
)
# Limited additional confounding. Education uses 3 levels to conserve design df.
d$education <- factor(ifelse(d$DMDEDUC2 %in% 1:2, "Below high school", ifelse(d$DMDEDUC2 %in% 3, "High school", ifelse(d$DMDEDUC2 %in% 4:5, "Above high school", NA))),
  levels = c("Above high school", "High school", "Below high school")
)
d$pir <- d$INDFMPIR
d$cancer <- ifelse(d$MCQ220 %in% 1:2, as.numeric(d$MCQ220 == 1), NA_real_)
d$expanded_ok <- d$core & complete.cases(d[c("education", "pir", "cancer")])
S <- design(d)
X <- subset(S, expanded_ok)
fit("Base adjustment in expanded complete cases", X, paste("lower_z+", base))
expanded <- paste("lower_z+", base, "+education+pir+cancer")
fit("Additional education income and cancer", X, expanded)
out(
  data.frame(variable = c("education", "pir", "cancer"), missing = sapply(d[c("education", "pir", "cancer")], function(x) sum(is.na(x[d$core]))), N = sum(d$core)),
  "expanded_missingness"
)
# Imputation sensitivity analysis; it does not expand score or grip eligibility.
# Impute income via PMM; exclude the single unknown cancer response rather than infer cancer.
midom <- d$core & !is.na(d$education) & !is.na(d$cancer)
miidx <- which(midom)
x <- d[miidx, ]
nafit <- survfit(Surv(time, death) ~ 1, data = x)
x$nelson_aalen <- stepfun(nafit$time, c(0, nafit$cumhaz))(x$time)
mi <- x[c(
  "pir",
  "age",
  "lower_z",
  "height",
  "bmi",
  "sbp",
  "tc",
  "hdl",
  "egfr",
  "diabetes",
  "smoking",
  "bp_tx",
  "statin",
  "education",
  "cancer",
  "death",
  "nelson_aalen",
  "r1",
  "r2",
  "r3",
  "a1",
  "a2",
  "a3"
)]
mi$sex <- factor(x$sex)
mi$log_weight <- log(x$w)
mi$design_cell <- interaction(x$SDMVSTRA, x$SDMVPSU, drop = TRUE)
method <- rep("", ncol(mi))
names(method) <- names(mi)
method["pir"] <- "pmm"
pm <- matrix(0, ncol(mi), ncol(mi), dimnames = list(names(mi), names(mi)))
pm["pir", setdiff(names(mi), "pir")] <- 1
# Linear age is already spanned by the fixed age spline; avoid redundant predictors.
pm["pir", "age"] <- 0
imp <- mice(mi,
  m = 40,
  maxit = 10,
  method = method,
  predictorMatrix = pm,
  seed = 410031,
  printFlag = FALSE,
  remove.collinear = FALSE
)
out(
  if (is.null(imp$loggedEvents)) data.frame(status = "No logged events") else imp$loggedEvents,
  "imputation_logged_events"
)
imods <- list()
idata <- list()
for (k in 1:40) {
  z <- d
  v <- complete(imp, k)
  z$pir[miidx] <- v$pir
  z$mi_domain <- midom
  dd <- subset(design(z), mi_domain)
  m <- svycoxph(as.formula(paste("Surv(time,death)~", expanded)), dd, x = TRUE, model = TRUE)
  imods[[k]] <- m
  idata[[k]] <- data.frame(
    imputation = k,
    imputed_mean = mean(v$pir[is.na(mi$pir)]),
    imputed_sd = sd(v$pir[is.na(mi$pir)]),
    beta = coef(m)["lower_z"],
    SE = sqrt(vcov(m)["lower_z", "lower_z"])
  )
}
dfcomp <- degf(subset(S, core)) - length(coef(imods[[1]])) + 1
pooled <- MIcombine(lapply(imods, coef), lapply(imods, vcov), df.complete = dfcomp)
j <- match("lower_z", names(pooled$coefficients))
beta <- pooled$coefficients[j]
se <- sqrt(pooled$variance[j, j])
df <- pooled$df[j]
out(
  data.frame(model = "Income multiply imputed expanded adjustment", N = length(miidx), deaths = sum(d$death[miidx]), imputations = 40, HR = exp(beta), lower = exp(beta - qt(.975, df) * se), upper = exp(beta + qt(.975, df) * se), P = 2 * pt(-abs(beta / se), df), beta = beta, SE = se, df = df, missing_information = pooled$missinfo[j]),
  "imputation_result"
)
out(bind_rows(idata), "imputation_diagnostics")
out(
  data.frame(iteration = rep(1:10, 40), chain = rep(1:40, each = 10), mean = as.vector(imp$chainMean["pir", , ]), variance = as.vector(imp$chainVar["pir", , ])),
  "imputation_traces"
)
saveRDS(imp, file.path(V, "Derived/income_imputations.rds"))
# Joint-term residual diagnostics use identical original bases and coefficient checks.
ph <- list()
phobj <- list()
for (coh in c("core", "mobility_cohort")) {
  des <- subset(S, get(coh))
  cv <- function(x) paste0("c(", paste(sprintf("%.17g", x), collapse = ","), ")")
  f <- as.formula(paste(
    "Surv(time,death)~lower_z+",
    if (coh == "mobility_cohort") "mobility+" else "",
    sprintf(
      "ns(logitrisk,knots=%s,Boundary.knots=%s)+ns(age,knots=%s,Boundary.knots=%s)+sex_factor+cycle+height+bmi",
      cv(attr(A$risk_basis, "knots")),
      cv(attr(A$risk_basis, "Boundary.knots")),
      cv(attr(A$age_basis, "knots")),
      cv(attr(A$age_basis, "Boundary.knots"))
    )
  ))
  m <- svycoxph(f, des, x = TRUE, y = TRUE, model = TRUE)
  ref <- A$models[[if (coh == "core") "Primary: continuous strength" else "Mobility D: strength and limitation"]]
  stopifnot(abs(coef(m)["lower_z"] - coef(ref)["lower_z"]) < 1e-8)
  q <- cox.zph(m, transform = "km", terms = TRUE)
  tt <- rownames(q$table)
  tt[grepl("^ns\\(logitrisk", tt)] <- "PREVENT risk spline block"
  tt[grepl("^ns\\(age", tt)] <- "Age spline block"
  ph[[coh]] <- data.frame(cohort = coh, term = tt, q$table, row.names = NULL)
  phobj[[coh]] <- q
}
out(bind_rows(ph), "PH_joint_terms")
saveRDS(phobj, file.path(V, "Derived/ph_joint.rds"))
out(bind_rows(est), "additional_sensitivity_estimates")
out(bind_rows(diag), "additional_sensitivity_diagnostics")
saveRDS(
  list(data = d, models = models, grip_knots = kn, enhanced_basis = ek, MI_models = imods),
  file.path(V, "Derived/additional_models.rds")
)
capture.output(sessionInfo(), file = file.path(V, "Audit/sessionInfo_additional.txt"))
print(bind_rows(est))
print(read.csv(file.path(V, "Results/imputation_result.csv")))
print(read.csv(file.path(V, "Results/nonlinearity_tests.csv")))
