# Survey-weighted Cox models (primary model, older-cohort models A-D, sensitivity models), with
# estimates, model diagnostics, proportional-hazards tests, and weighted descriptive statistics.
# Run through code/03_primary_analysis.R (working directory: the build directory).

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
d <- readRDS(file.path(B, "03_Derived/study4_master.rds"))
sc <- read.csv(file.path(B, "05_Outputs/strength_scaling.csv"))
# All basis knots use covariates only, without outcomes, and are fixed thereafter.
rk <- ns(d$logitrisk[d$core], df = 3)
ak <- ns(d$age[d$core], df = 3)
addbasis <- function(x) {
  r <- predict(rk, x$logitrisk)
  a <- predict(ak, x$age)
  for (j in 1:3) {
    x[[paste0("r", j)]] <- r[, j]
    x[[paste0("a", j)]] <- a[, j]
  }
  x
}
d <- addbasis(d)
out(
  data.frame(variable = rep(c("logit_PREVENT", "age"), each = 4), role = rep(c("knot1", "knot2", "boundary_low", "boundary_high"), 2), value = c(attr(rk, "knots"), attr(rk, "Boundary.knots"), attr(ak, "knots"), attr(ak, "Boundary.knots"))),
  "spline_specification"
)
design <- function(x) svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~w,
  nest = TRUE,
  data = x[x$design_ok, ]
)
S <- design(d)
C <- subset(S, core)
M <- subset(S, mobility_cohort)
base <- "r1+r2+r3+a1+a2+a3+sex_factor+cycle+height+bmi"
models <- list()
estimates <- list()
diagnostics <- list()
warn <- list()
fit <- function(name, des, rhs, event = "Surv(time,death)", targets = c("lower_z", "mobility")) {
  f <- as.formula(paste(event, "~", rhs))
  ww <- character()
  m <- withCallingHandlers(svycoxph(f, design = des, x = TRUE, y = TRUE, model = TRUE),
    warning = function(w) {
      ww <<- c(ww, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  p <- length(coef(m))
  dd <- degf(des) - p + 1
  stopifnot(dd > 5, all(is.finite(coef(m))), all(is.finite(vcov(m))))
  se <- sqrt(diag(vcov(m)))
  q <- qt(.975, dd)
  pv <- 2 * pt(-abs(coef(m) / se), dd)
  estimates[[name]] <<- data.frame(
    model = name,
    term = names(coef(m)),
    beta = unname(coef(m)),
    SE = se,
    HR = exp(coef(m)),
    lower = exp(coef(m) - q * se),
    upper = exp(coef(m) + q * se),
    P = pv,
    df = dd,
    N = nrow(m$model),
    deaths = m$nevent,
    row.names = NULL
  )
  xx <- scale(m$x)
  keep <- apply(xx, 2, function(x) all(is.finite(x)))
  diagnostics[[name]] <<- data.frame(
    model = name,
    N = nrow(m$model),
    events = m$nevent,
    parameters = p,
    design_df = degf(des),
    inference_df = dd,
    condition_number = kappa(xx[, keep, drop = FALSE], exact = TRUE),
    min_cov_eigen = min(eigen(vcov(m), symmetric = TRUE, only.values = TRUE)$values),
    warnings = paste(ww, collapse = " | ")
  )
  models[[name]] <<- m
  warn[[name]] <<- ww
  invisible(m)
}
fit("Primary: continuous strength", C, paste("lower_z+", base))
fit("Mobility A: risk and covariates", M, base)
fit("Mobility B: add reported limitation", M, paste("mobility+", base))
fit("Mobility C: add measured strength", M, paste("lower_z+", base))
fit("Mobility D: strength and limitation", M, paste("lower_z+mobility+", base))
fit("Sensitivity: per 5 kg lower", C, paste("lower5+", base), targets = "lower5")
fit("Sensitivity: attained age", C, paste("lower_z+", base), "Surv(age,age+time,death)")
# Landmark eligibility uses observed survival and complete follow-up beyond two years.
d$landmark <- d$core & d$time > 2
d$time_lm <- d$time - 2
SL <- design(d)
L <- subset(SL, landmark)
fit("Sensitivity: two-year landmark", L, paste("lower_z+", base), "Surv(time_lm,death)")
one <- d
one$lower_z <- (sc$mean_kg[match(one$sex, sc$sex)] - one$grip_any) / sc$sd_kg[match(
  one$sex,
  sc$sex
)]
one$onehand_cohort <- one$pregrip & !is.na(one$grip_any)
O <- subset(design(one), onehand_cohort)
fit("Sensitivity: any valid hand", O, paste("lower_z+", base))
fit("Sensitivity: no reported hand pain", subset(S, core & !hand_pain), paste("lower_z+", base))
aa <- d
sf <- readRDS(file.path(B, "03_Derived/scoring_functions.rds"))
coefs <- sf$coefs
terms <- sf$terms
score <- sf$score
aa$diabetes <- aa$diabetes_a1c
aa$complete_score <- aa$complete_score & !is.na(aa$diabetes)
aa$risk <- score(aa)
aa$logitrisk <- qlogis(aa$risk)
aa <- addbasis(aa)
aa$a1c_cohort <- aa$core & !is.na(aa$risk)
fit(
  "Sensitivity: HbA1c-inclusive diabetes",
  subset(design(aa), a1c_cohort),
  paste("lower_z+", base)
)
# Prespecified component replacement; linear components avoid excessive degrees of freedom.
fit(
  "Sensitivity: individual clinical factors",
  C,
  "lower_z+age+sex_factor+cycle+height+bmi+sbp+tc+hdl+egfr+diabetes+smoking+bp_tx+statin"
)
out(bind_rows(estimates), "model_estimates_all")
out(bind_rows(diagnostics), "model_diagnostics")
focal <- bind_rows(estimates) %>% filter(term %in% c("lower_z", "lower5", "mobility"))
out(focal, "focal_estimates")
# Conventional residual diagnostics: approximate only, no design-based P claim.
ph <- lapply(c("Primary: continuous strength", "Mobility D: strength and limitation"), function(n) {
  z <- cox.zph(models[[n]], transform = "km")
  saveRDS(
    z,
    file.path(B, "03_Derived", if (grepl("Primary", n)) "ph_primary.rds" else "ph_mobility.rds")
  )
  data.frame(model = n, term = rownames(z$table), z$table, row.names = NULL)
})
out(bind_rows(ph), "PH_exploratory")
# Weighted descriptive characteristics; estimates and SD, not standard errors of the mean.
tab <- list()
for (coh in c("Core age 50-79", "Mobility age 60-79")) {
  des <- if (grepl("Core", coh)) C else M
  for (v in c("age", "grip", "height", "bmi", "sbp", "tc", "hdl", "egfr", "risk")) {
    z <- svymean(as.formula(paste0("~", v)), des)
    s <- sqrt(as.numeric(svyvar(as.formula(paste0("~", v)), des)))
    tab[[paste(coh, v)]] <- data.frame(
      cohort = coh,
      variable = v,
      kind = "mean_SD",
      estimate = as.numeric(z),
      spread = s,
      N = nrow(des$variables)
    )
  }
  for (v in c("diabetes", "smoking", "bp_tx", "statin", "grip_low")) {
    z <- svymean(as.formula(paste0("~", v)), des)
    tab[[paste(coh, v)]] <- data.frame(
      cohort = coh,
      variable = v,
      kind = "proportion_SE",
      estimate = as.numeric(z),
      spread = SE(z)[1],
      N = nrow(des$variables)
    )
  }
  z <- svymean(~ I(sex == 2), des)
  tab[[paste(coh, "female")]] <- data.frame(
    cohort = coh,
    variable = "female",
    kind = "proportion_SE",
    estimate = coef(z)[2],
    spread = SE(z)[2],
    N = nrow(des$variables)
  )
}
out(bind_rows(tab), "baseline_characteristics")
out(
  data.frame(cohort = "Mobility age 60-79", prevalence = coef(svymean(~mobility, M))[1], SE = SE(svymean(~mobility, M))[1], N = nrow(M$variables), with_limitation = sum(M$variables$mobility)),
  "mobility_prevalence"
)
# Continuous-exposure central distribution supports, not effect estimates.
support <- bind_rows(lapply(
  split(d[d$core, ], interaction(d$sex_factor[d$core], cut(d$age[d$core], c(49, 59, 69, 79)), drop = TRUE)),
  function(x) data.frame(
    sex = as.character(x$sex_factor[1]),
    age_min = min(x$age),
    age_max = max(x$age),
    N = nrow(x),
    deaths = sum(x$death),
    grip_p10 = quantile(x$grip, .1),
    grip_p90 = quantile(x$grip, .9),
    risk_p10 = quantile(x$risk, .1),
    risk_p90 = quantile(x$risk, .9)
  )
))
out(support, "covariate_support")
saveRDS(
  list(models = models, risk_basis = rk, age_basis = ak, data = d),
  file.path(B, "03_Derived/analysis_models.rds")
)
capture.output(sessionInfo(), file = file.path(B, "08_Audit/sessionInfo_analysis.txt"))
print(focal)
print(bind_rows(diagnostics))
print(bind_rows(ph))
