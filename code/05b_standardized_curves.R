# Stage 5b. Standardized cumulative mortality over follow-up at the 25th and 75th percentiles of grip
# strength (Figure 3A and 3B), and survey-weighted Kaplan-Meier mortality by sex-specific third of grip
# strength (Figure S2). Uses the primary model, contrast values, target population, and jackknife
# procedure of code/analysis/03_standardized_mortality.R; the 5-year values are checked against
# standardized_5yr_risks.csv.

source("code/00_setup.R")
suppressPackageStartupMessages({
  library(survey)
  library(survival)
  library(splines)
})
options(survey.lonely.psu = "fail")

A <- readRDS(file.path(P$derived, "analysis_models.rds"))
d <- A$data
S <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~w,
  nest = TRUE,
  data = d[d$design_ok, ]
)
C <- subset(S, core)
c <- C$variables
sc <- res("strength_scaling")
quant <- res("risk_strength_contrasts")
sup <- res("risk_support_cells")
c$cell <- interaction(c$sex,
  cut(c$age, c(49, 59, 69, 79)),
  cut(c$bmi, c(18.49, 25, 30, 39.91), right = FALSE),
  drop = TRUE
)
c$supported <- sup$support[match(as.character(c$cell), sup$cell)]
stopifnot(sum(c$supported) == res("risk_support_summary")$supported_N)
R <- subset(as.svrepdesign(S, type = "JKn", mse = TRUE), core)
rw <- weights(R, type = "analysis")
fw <- weights(R, type = "sampling")
base <- "r1+r2+r3+a1+a2+a3+sex_factor+cycle+height+bmi"
f <- as.formula(paste("Surv(time,death)~lower_z+", base))
grid <- sort(unique(round(
  c(seq(0, 8, by = 1 / 12), 0:8),
  10
))) # monthly grid; whole years included exactly

curves <- function(w) {
  idx <- which(w > 0)
  z <- c[idx, ]
  z$ww <- w[idx]
  fit <- coxph(f, data = z, weights = ww, ties = "efron", x = TRUE, model = TRUE)
  bh <- basehaz(fit, centered = TRUE)
  H <- vapply(grid, function(t) {
    h <- bh$hazard[bh$time <= t]
    if (length(h)) tail(h, 1) else 0
  }, numeric(1))
  zz <- z[z$supported, ]
  one <- function(col) {
    zz$grip <- quant[[col]][match(zz$sex, quant$sex)]
    zz$lower_z <- (sc$mean_kg[match(zz$sex, sc$sex)] - zz$grip) / sc$sd_kg[match(zz$sex, sc$sex)]
    lp <- predict(fit, newdata = zz, type = "lp", reference = "sample")
    risk <- -expm1(-outer(exp(lp), H))
    colSums(risk * zz$ww) / sum(zz$ww)
  }
  lo <- one("q25")
  hi <- one("q75")
  c(lo, hi, lo - hi)
}
point <- curves(fw)
reps <- t(vapply(seq_len(ncol(rw)), function(k) curves(rw[, k]), numeric(length(point))))
stopifnot(all(is.finite(reps)))
vv <- R$scale * colSums(sweep(sweep(reps, 2, point, "-")^2, 1, R$rscales, "*"))
df <- degf(C) - length(coef(A$models[["Primary: continuous strength"]])) + 1
q <- qt(0.975, df)
G <- length(grid)
out <- data.frame(
  time = rep(grid, 3),
  series = rep(c("lower_strength_q25", "higher_strength_q75", "difference"), each = G),
  estimate = point,
  SE = sqrt(vv),
  lower = point - q * sqrt(vv),
  upper = point + q * sqrt(vv),
  df = df,
  replicates = nrow(reps)
)
# Exact agreement with the prespecified 5-year estimates.
ref <- res("standardized_5yr_risks")
chk <- c(
  out$estimate[out$time == 5 & out$series == "lower_strength_q25"] - ref$estimate[ref$estimand == "q25_risk"],
  out$estimate[out$time == 5 & out$series == "higher_strength_q75"] - ref$estimate[ref$estimand == "q75_risk"],
  out$SE[out$time == 5 & out$series == "difference"] - ref$SE[ref$estimand == "difference"]
)
stopifnot(max(abs(chk)) < 1e-10)
out$at_risk <- vapply(out$time, function(t) sum(c$time >= t), numeric(1))
dir.create(file.path(P$csv, "analysis_outputs"), showWarnings = FALSE)
write.csv(out,
  file.path(P$csv, "analysis_outputs", "standardized_mortality_curves.csv"),
  row.names = FALSE
)
message(
  "Standardized curves: 5-year values reproduce the prespecified estimates (max |diff| ",
  signif(max(abs(chk)), 2), ")"
)

# Descriptive survey-weighted Kaplan-Meier cumulative mortality by sex-specific weighted grip tertile.
tt <- do.call(rbind, lapply(1:2, function(sx) {
  qv <- as.numeric(coef(svyquantile(~grip, subset(C, sex == sx), c(1 / 3, 2 / 3), ci = FALSE)))
  data.frame(sex = sx, t1 = qv[1], t2 = qv[2])
}))
C$variables$tertile <- with(
  C$variables,
  ifelse(grip <= tt$t1[match(sex, tt$sex)], 1, ifelse(grip <= tt$t2[match(sex, tt$sex)], 2, 3))
)
C$variables$tertile <- factor(C$variables$tertile,
  levels = c(1, 2, 3),
  labels = c("Lowest third", "Middle third", "Highest third")
)
km <- svykm(Surv(time, death) ~ tertile, design = C, se = TRUE)
kmd <- do.call(rbind, lapply(names(km), function(k) {
  s <- km[[k]]
  data.frame(tertile = k, time = s$time, surv = s$surv, varlog = s$varlog)
}))
kmd$cum_mortality <- 1 - kmd$surv
se_s <- sqrt(kmd$varlog) * kmd$surv
kmd$lower <- pmax(0, 1 - pmin(1, kmd$surv + 1.96 * se_s))
kmd$upper <- pmin(1, 1 - pmax(0, kmd$surv - 1.96 * se_s))
risk_tab <- do.call(rbind, lapply(levels(C$variables$tertile), function(k) {
  v <- C$variables[C$variables$tertile == k, ]
  data.frame(
    tertile = k, time = 0:8, at_risk = vapply(0:8, function(t) sum(v$time >= t), numeric(1)),
    deaths = sum(v$death), n = nrow(v)
  )
}))
write.csv(kmd, file.path(P$csv, "analysis_outputs", "km_by_grip_tertile.csv"), row.names = FALSE)
write.csv(risk_tab,
  file.path(P$csv, "analysis_outputs", "km_numbers_at_risk.csv"),
  row.names = FALSE
)
write.csv(tt, file.path(P$csv, "analysis_outputs", "grip_tertile_cutpoints.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(P$logs, "sessionInfo_standardized_curves.txt"))
