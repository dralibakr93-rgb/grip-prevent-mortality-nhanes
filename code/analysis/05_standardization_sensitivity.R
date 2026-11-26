# Standardized 5-year mortality under alternative target populations and models, each with
# jackknife refits. Run through code/05_standardization.R (working directory: the build directory).

suppressPackageStartupMessages({
  library(survey)
  library(survival)
  library(dplyr)
  library(Hmisc)
})
options(survey.lonely.psu = "fail")
B <- "study"
V <- file.path(B, "additional")
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
ad <- readRDS(file.path(V, "Derived/additional_models.rds"))
S <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~w,
  nest = TRUE,
  data = d[d$design_ok, ]
)
C <- subset(S, core)
c <- C$variables
sc <- read.csv(file.path(B, "05_Outputs/strength_scaling.csv"))
quant <- read.csv(file.path(B, "05_Outputs/risk_strength_contrasts.csv"))
sup <- read.csv(file.path(B, "05_Outputs/risk_support_cells.csv"))
c$cell <- interaction(c$sex,
  cut(c$age, c(49, 59, 69, 79)),
  cut(c$bmi, c(18.49, 25, 30, 39.91), right = FALSE),
  drop = TRUE
)
c$original_support <- sup$support[match(c$cell, sup$cell)]
c$age_cell <- interaction(c$sex, cut(c$age, c(49, 59, 69, 79)), drop = TRUE)
ss <- bind_rows(lapply(split(c, c$age_cell), function(z) {
  q <- quant[quant$sex == z$sex[1], ]
  lo <- quantile(z$grip, .05)
  hi <- quantile(z$grip, .95)
  data.frame(
    cell = as.character(z$age_cell[1]),
    N = nrow(z),
    q05 = lo,
    q95 = hi,
    support = nrow(z) >= 30 & q$q25 >= lo & q$q75 <= hi
  )
}))
c$sex_age_support <- ss$support[match(c$age_cell, ss$cell)]
c$full_target <- TRUE
bk <- quantile(c$risk, c(.25, .5, .75))
c$risk_band <- cut(c$risk, c(-Inf, bk, Inf), labels = FALSE)
g <- rcspline.eval(c$lower_z, knots = ad$grip_knots, inclx = TRUE)
c$g1 <- g[, 1]
c$g2 <- g[, 2]
out(ss, "alternative_support_cells")
targets <- c("original_support", "sex_age_support", "full_target")
out(bind_rows(lapply(
  c("original_support", "outside_original_support", "sex_age_support", "full_target"),
  function(t) {
    i <- if (t == "outside_original_support") !c$original_support else c[[t]]
    z <- c[i, ]
    data.frame(
      target = t,
      N = nrow(z),
      deaths = sum(z$death),
      deaths_5 = sum(z$death == 1 & z$time <= 5),
      weighted_fraction = sum(z$w) / sum(c$w),
      age_unweighted = mean(z$age),
      age_weighted = weighted.mean(z$age, z$w),
      PREVENT_unweighted = mean(z$risk),
      PREVENT_weighted = weighted.mean(z$risk, z$w),
      grip_weighted = weighted.mean(z$grip, z$w)
    )
  }
)), "support_target_comparison")
C$variables <- c
R <- subset(as.svrepdesign(S, type = "JKn", mse = TRUE), core)
rw <- weights(R, type = "analysis")
fw <- weights(R, type = "sampling")
base <- "r1+r2+r3+a1+a2+a3+sex_factor+cycle+height+bmi"
specs <- list(
  original = as.formula(paste("Surv(time,death)~lower_z+", base)),
  risk_stratified = as.formula(paste("Surv(time,death)~lower_z+", base, "+strata(risk_band)")),
  nonlinear_grip = as.formula(paste("Surv(time,death)~g1+g2+", base))
)
calc <- function(w, f, kind) {
  z <- c[w > 0, ]
  z$ww <- w[w > 0]
  m <- coxph(f, z, weights = ww, ties = "efron", x = TRUE, model = TRUE)
  stopifnot(all(is.finite(coef(m))))
  bh <- basehaz(m, centered = TRUE)
  if (kind == "risk_stratified") {
    bhh <- sapply(split(bh, bh$strata), function(b) {
      a <- b$hazard[b$time <= 5]
      if (length(a)) tail(a, 1) else 0
    })
    lev <- as.character(sort(unique(z$risk_band)))
    stopifnot(length(lev) == length(bhh))
    names(bhh) <- lev
  } else hh <- tail(bh$hazard[bh$time <= 5], 1)
  val <- c()
  for (t in targets) {
    v <- sapply(c("q25", "q75"), function(q) {
      a <- z[z[[t]], ]
      a$grip <- quant[[q]][match(a$sex, quant$sex)]
      a$lower_z <- (sc$mean_kg[match(a$sex, sc$sex)] - a$grip) / sc$sd_kg[match(a$sex, sc$sex)]
      gg <- rcspline.eval(a$lower_z, knots = ad$grip_knots, inclx = TRUE)
      a$g1 <- gg[, 1]
      a$g2 <- gg[, 2]
      lp <- predict(m, newdata = a, type = "lp", reference = "sample")
      h <- if (kind == "risk_stratified") unname(bhh[as.character(a$risk_band)]) else hh
      weighted.mean(-expm1(-h * exp(lp)), a$ww)
    })
    val <- c(val, setNames(c(v, v[1] - v[2]), paste(t, c("q25", "q75", "difference"), sep = ":")))
  }
  list(values = val, fit = m, bh = bh)
}
ans <- list()
logs <- list()
checks <- list()
reps <- list()
for (k in names(specs)) {
  f <- specs[[k]]
  a <- calc(fw, f, k)
  point <- a$values
  # Independent survfit check for each baseline stratum or one observation.
  ix <- if (k == "risk_stratified") sapply(
    sort(unique(c$risk_band)),
    function(b) which(c$risk_band == b)[1]
  ) else 1
  for (j in ix) {
    new <- c[j, , drop = FALSE]
    bh <- a$bh
    if (k == "risk_stratified") {
      levelsbh <- unique(as.character(bh$strata))
      z <- bh[as.character(bh$strata) == levelsbh[match(new$risk_band, sort(unique(c$risk_band)))], ]
    } else z <- bh
    hv <- tail(z$hazard[z$time <= 5], 1)
    manual <- exp(-hv * exp(predict(a$fit, newdata = new, type = "lp", reference = "sample")))
    direct <- as.numeric(summary(survfit(a$fit, newdata = new), times = 5)$surv)
    stopifnot(abs(manual - direct) < 1e-9)
    checks[[paste(k, j)]] <- data.frame(
      model = k,
      stratum = if (k == "risk_stratified") new$risk_band else 0,
      difference = abs(manual - direct)
    )
  }
  mat <- matrix(NA_real_, ncol(rw), length(point), dimnames = list(NULL, names(point)))
  for (j in seq_len(ncol(rw))) {
    ww <- character()
    err <- ""
    val <- tryCatch(withCallingHandlers(calc(rw[, j], f, k)$values, warning = function(w) {
      ww <<- c(ww, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) {
      err <<- conditionMessage(e)
      rep(NA_real_, length(point))
    })
    mat[j, ] <- val
    logs[[paste(k, j)]] <- data.frame(
      model = k,
      replicate = j,
      finite = all(is.finite(val)),
      warnings = paste(ww, collapse = " | "),
      error = err
    )
  }
  stopifnot(all(is.finite(mat)))
  vv <- R$scale * colSums(sweep(sweep(mat, 2, point, "-")^2, 1, R$rscales, "*"))
  df <- degf(C) - length(coef(a$fit)) + 1
  spl <- strsplit(names(point), ":", fixed = TRUE)
  ans[[k]] <- data.frame(
    model = k,
    target = vapply(spl, "[", "", 1),
    estimand = vapply(spl, "[", "", 2),
    estimate = point,
    SE = sqrt(vv),
    lower = point - qt(.975, df) * sqrt(vv),
    upper = point + qt(.975, df) * sqrt(vv),
    df = df,
    replicates = nrow(mat),
    row.names = NULL
  )
  reps[[k]] <- mat
  cat("Completed absolute-risk specification:", k, "\n")
}
out(bind_rows(ans), "absolute_risk_sensitivities")
out(bind_rows(logs), "absolute_risk_replication_log")
out(bind_rows(checks), "stratified_prediction_validation")
ref <- read.csv(file.path(B, "05_Outputs/standardized_5yr_risks.csv"))
stopifnot(max(abs(ans$original$estimate[ans$original$target == "original_support"] - ref$estimate)) < 1e-10)
saveRDS(reps, file.path(V, "Derived/risk_replicates.rds"))
capture.output(sessionInfo(), file = file.path(V, "Audit/sessionInfo_absolute_risk.txt"))
print(bind_rows(ans))
