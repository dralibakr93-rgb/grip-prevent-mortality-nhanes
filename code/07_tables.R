# Stage 7. Table 1 (CSV and Word) and the formatted values behind the supplementary figures
# (output/csv/supplementary_values/). In footnotes, ^a^ marks a superscript letter.

invisible(Sys.setlocale("LC_ALL", "en_US.UTF-8"))
source("code/00_setup.R")
suppressPackageStartupMessages({
  library(flextable)
  library(officer)
})
# Footnote formatting: ^x^ marks a superscript footnote letter; figure and table citations are set in bold.
CITE_RE <- "(?<![A-Za-z])((?:Figures?|Tables?) S?[0-9]+[A-E]?(?:\u2013S?[0-9]+[A-E]?)?)(?![0-9A-Za-z])"
rich_fpar <- function(x, props) {
  parts <- regmatches(x, gregexpr("\\^[^^]+\\^|[^^]+", x))[[1]]
  runs <- list()
  add <- function(txt, p) if (nzchar(txt)) runs[[length(runs) + 1]] <<- officer::ftext(txt, p)
  for (p in parts) {
    if (grepl("^\\^.+\\^$", p)) {
      add(gsub("\\^", "", p), update(props, vertical.align = "superscript"))
      next
    }
    m <- gregexpr(CITE_RE, p, perl = TRUE)[[1]]
    if (m[1] == -1) {
      add(p, props)
      next
    }
    pos <- 1
    for (k in seq_along(m)) {
      s <- m[k]
      e <- s + attr(m, "match.length")[k] - 1
      if (s > pos) add(substr(p, pos, s - 1), props)
      add(substr(p, s, e), update(props, bold = TRUE))
      pos <- e + 1
    }
    if (pos <= nchar(p)) add(substr(p, pos, nchar(p)), props)
  }
  do.call(officer::fpar, runs)
}
strip_sup <- function(x) gsub("\\^([^^]+)\\^", "\\1", x) # remove ^x^ markup for CSV values

NBSP <- "\u00a0"
EN <- "\u2013"
fmt_n <- function(x) vapply(
  round(x),
  function(v) if (is.na(v)) "" else formatC(v, big.mark = ",", format = "d"),
  character(1)
)
fx <- function(x, d = 2) formatC(x, format = "f", digits = d)
ci <- function(h, l, u, d = 2) sprintf("%s (%s%s%s)", fx(h, d), fx(l, d), EN, fx(u, d))
pfmt <- function(p) ifelse(p < 0.001, "<0.001", ifelse(p < 0.01, fx(p, 3), fx(p, 2)))
FONT <- "Times New Roman"

style_ft <- function(ft, header_rows = 1) {
  ft <- font(ft, fontname = FONT, part = "all")
  ft <- fontsize(ft, size = 9.5, part = "all")
  ft <- bold(ft, part = "header")
  ft <- bg(ft, bg = "#EDEDED", part = "header")
  ft <- border_remove(ft)
  ft <- hline_top(ft, border = fp_border(color = "black", width = 1), part = "header")
  ft <- hline_bottom(ft, border = fp_border(color = "black", width = 0.75), part = "header")
  ft <- hline_bottom(ft, border = fp_border(color = "black", width = 1), part = "body")
  ft <- valign(ft, valign = "top", part = "all")
  ft <- padding(ft,
    padding.top = 1.5,
    padding.bottom = 1.5,
    padding.left = 3,
    padding.right = 3,
    part = "all"
  )
  ft <- align(ft, j = 1, align = "left", part = "all")
  if (ncol_keys(ft) > 1) ft <- align(ft, j = 2:ncol_keys(ft), align = "center", part = "all")
  set_table_properties(ft, layout = "autofit", width = 1)
}
section_rows <- function(ft, rows) { # bold, shaded section-header rows spanning the table
  if (!length(rows)) return(ft)
  ft <- merge_h_range(ft, i = rows, j1 = 1, j2 = ncol_keys(ft))
  ft <- bold(ft, i = rows, part = "body")
  ft <- bg(ft, i = rows, bg = "#F7F7F7", part = "body")
  align(ft, i = rows, j = 1, align = "left", part = "body")
}
indent_rows <- function(
  ft,
  rows
) if (length(rows)) padding(ft,
  i = rows,
  j = 1,
  padding.left = 12,
  part = "body"
) else ft
TABLES <- list()
# Only Table 1 is a formatted table. The other summaries are written as CSV files.
register <- function(id, title, df, footnote, ft, dir = NULL) {
  if (id != "Table_1") {
    dir.create(file.path(P$csv, "supplementary_values"), showWarnings = FALSE)
    write.csv(df, file.path(P$csv, "supplementary_values", paste0(id, ".csv")), row.names = FALSE, fileEncoding = "UTF-8")
    return(invisible(NULL))
  }
  write.csv(df, file.path(P$csv, paste0(id, ".csv")), row.names = FALSE, fileEncoding = "UTF-8")
  TABLES[[id]] <<- list(id = id, title = title, ft = ft, footnote = footnote, data = df)
  doc <- read_docx()
  doc <- body_add_fpar(doc, fpar(ftext(
    paste0(gsub("_", " ", id), ". ", title),
    fp_text(font.family = FONT, font.size = 10, bold = TRUE)
  )))
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_fpar(doc, rich_fpar(footnote, fp_text(font.family = FONT, font.size = 8.5)))
  print(doc, target = file.path(dir, paste0(id, ".docx")))
}

# TABLE 1
t1 <- read.csv(file.path(P$csv, "table1_characteristics_long.csv"), stringsAsFactors = FALSE)
g <- function(k, var, lv = "") t1[t1$cohort == k & t1$variable == var & t1$level == lv, ]
cell <- function(k, var, lv = "", d = 1, scale = 1) {
  r <- if (lv == "" && var %in% c(
    "female",
    "diabetes",
    "smoking",
    "bp_tx",
    "statin",
    "grip_low",
    "mobility",
    "prior_hosp"
  )) g(
    k,
    var,
    "1"
  ) else g(
    k,
    var,
    lv
  )
  if (!nrow(r)) return("Not assessed")
  switch(r$type,
    count = fmt_n(r$estimate),
    median = fx(r$estimate, 2),
    mean_sd = sprintf("%s (%s)", fx(r$estimate * scale, d), fx(r$spread * scale, d)),
    n_pct = sprintf("%s (%s)", fmt_n(r$n), fx(r$estimate, 1))
  )
}
spec <- list(
  list("Participants, n", "participants"), list("Deaths, n", "deaths"),
  list("Median follow-up, y", "median_followup_years"),
  list("Demographic characteristics", NA),
  list("Age, y", "age"), list("Women", "female"),
  list("Race and ethnicity", NA, sub = TRUE),
  list("Non-Hispanic White", "race_eth", "Non-Hispanic White", indent = TRUE),
  list("Non-Hispanic Black", "race_eth", "Non-Hispanic Black", indent = TRUE),
  list("Hispanic", "race_eth", "Hispanic", indent = TRUE),
  list("Non-Hispanic Asian", "race_eth", "Non-Hispanic Asian", indent = TRUE),
  list("Other or multiracial", "race_eth", "Other or multiracial", indent = TRUE),
  list("Education", NA, sub = TRUE),
  list("Less than high school", "education", "Less than high school", indent = TRUE),
  list("High school or equivalent", "education", "High school or equivalent", indent = TRUE),
  list("More than high school", "education", "More than high school", indent = TRUE),
  list("Income-to-poverty ratio^a^", "pir"),
  list("Physical function", NA),
  list("Grip strength in men, kg", "grip_men"), list("Grip strength in women, kg", "grip_women"),
  list("Reported mobility limitation^b^", "mobility"),
  list("Cardiovascular risk factors", NA),
  list("PREVENT 10-y total CVD risk, %", "risk", scale = 100),
  list("PREVENT 10-y ASCVD risk, %", "risk_ascvd", scale = 100),
  list("Systolic blood pressure, mm Hg", "sbp"), list("Total cholesterol, mg/dL", "tc"),
  list("HDL cholesterol, mg/dL", "hdl"), list("eGFR, mL/min per 1.73 m\u00b2", "egfr"),
  list("Body-mass index, kg/m\u00b2", "bmi"), list("Height, cm", "height"),
  list("Diagnosed diabetes", "diabetes"), list("Current smoking", "smoking"),
  list("Antihypertensive treatment", "bp_tx"), list("Statin use", "statin"),
  list("Recent health care use", NA),
  list("Overnight hospital stay in the prior year^c^", "prior_hosp")
)
rows <- lapply(spec, function(s) {
  if (is.na(s[[2]])) return(data.frame(
    Characteristic = s[[1]],
    Primary = "",
    Mobility = "",
    kind = if (isTRUE(s$sub)) "sub" else "section"
  ))
  lv <- if (length(s) >= 3 && is.character(s[[3]])) s[[3]] else ""
  sc <- if (!is.null(s$scale)) s$scale else 1
  data.frame(
    Characteristic = s[[1]], Primary = cell("primary", s[[2]], lv, scale = sc),
    Mobility = cell("mobility", s[[2]], lv, scale = sc),
    kind = if (isTRUE(s$indent)) "indent" else "row"
  )
})
T1 <- do.call(rbind, rows)
sup1 <- grep("\\^[a-z]\\^$", T1$Characteristic)
sup1_let <- sub(".*\\^([a-z])\\^$", "\\1", T1$Characteristic[sup1])
T1$Characteristic <- sub(
  "\\^[a-z]\\^$",
  "",
  T1$Characteristic
) # footnote letters added as superscripts below
pir_m <- g("primary", "pir")$n_missing
pir_mm <- g("mobility", "pir")$n_missing
n_p <- fmt_n(g("primary", "participants")$estimate)
n_m <- fmt_n(g("mobility", "participants")$estimate)
df1 <- T1[, 1:3]
names(df1) <- c(
  "Characteristic", sprintf("Primary cohort, age 50%s79 y (n=%s)", EN, n_p),
  sprintf("Older cohort, age 60%s79 y (n=%s)", EN, n_m)
)
ft1 <- flextable(df1) |> style_ft()
ft1 <- section_rows(ft1, which(T1$kind == "section"))
ft1 <- bold(ft1, i = which(T1$kind == "sub"), j = 1, part = "body")
ft1 <- indent_rows(ft1, which(T1$kind %in% c("indent")))
for (k in seq_along(sup1)) ft1 <- compose(ft1,
  i = sup1[k],
  j = 1,
  part = "body",
  value = as_paragraph(T1$Characteristic[sup1[k]], as_sup(sup1_let[k]))
)
fn1 <- paste0(
  "Values are unweighted number (survey-weighted percentage) or survey-weighted mean (SD); participant and death counts ",
  "and median follow-up are unweighted. The older cohort is a subset of the primary cohort. ASCVD indicates atherosclerotic cardiovascular disease; CVD, cardiovascular ",
  "disease; eGFR, estimated glomerular filtration rate; HDL, high-density lipoprotein; and PREVENT, Predicting Risk of ",
  "Cardiovascular Disease Events. ^a^Income-to-poverty ratio was missing for ",
  pir_m,
  " participants in the primary cohort and ",
  pir_mm,
  " in the older cohort. ^b^Defined only in participants aged 60",
  EN,
  "79 years. ^c^Self-reported overnight hospital stay in the ",
  "12 months before the interview (before grip testing)",
  local({
    m1 <- g("primary", "prior_hosp", "1")$n_missing
    m2 <- g("mobility", "prior_hosp", "1")$n_missing
    if (m1 == 0 && m2 == 0) "; status was known for all participants." else
      paste0(
        "; status was unknown for ",
        m1,
        " participants in the primary cohort and ",
        m2,
        " in the older cohort."
      )
  })
)
register("Table_1", "Characteristics of the Primary and Older Cohorts", df1, fn1, ft1, P$tab_main)

# SUPPLEMENTARY VALUES
# input availability and ranges among clinical candidates
re <- res("range_exclusions")
mi <- res("missingness")
lims <- c(
  tc = "130\u2013320 mg/dL",
  hdl = "20\u2013100 mg/dL",
  sbp = "90\u2013200 mm Hg",
  egfr = "15\u2013140 mL/min per 1.73 m\u00b2",
  bmi = "18.5\u201339.9 kg/m\u00b2"
)
lab <- c(
  tc = "Total cholesterol",
  hdl = "HDL cholesterol",
  sbp = "Systolic blood pressure",
  egfr = "eGFR",
  bmi = "Body-mass index",
  age = "Age",
  height = "Height",
  scr = "Serum creatinine",
  diabetes = "Diagnosed diabetes",
  smoking = "Current smoking",
  bp_tx = "Antihypertensive treatment",
  statin = "Statin use",
  grip = "Bilateral maximal-effort grip strength",
  mobility = "Mobility status"
)
den <- unique(mi$denominator)
s1a <- data.frame(
  Input = lab[re$variable], Range = lims[re$variable], Missing = fmt_n(re$missing),
  Below = fmt_n(re$below), Above = fmt_n(re$above)
)
oth <- mi[!mi$variable %in% re$variable, ]
s1b <- data.frame(
  Input = lab[oth$variable],
  Range = ifelse(oth$variable == "age", "50\u201379 y (eligibility)", "No range restriction"),
  Missing = fmt_n(oth$missing), Below = "", Above = ""
)
s1 <- rbind(
  data.frame(Input = "Inputs with range restrictions", Range = "", Missing = "", Below = "", Above = ""),
  s1a,
  data.frame(Input = "Other inputs and exposures", Range = "", Missing = "", Below = "", Above = ""),
  s1b
)
names(s1) <- c("Variable", "Study range", "Missing, n", "Below range, n", "Above range, n")
ft <- flextable(s1) |> style_ft()
ft <- section_rows(ft, c(1, nrow(s1a) + 2))
ft <- indent_rows(ft, setdiff(seq_len(nrow(s1)), c(1, nrow(s1a) + 2)))
register(
  "selection_input_availability",
  "Input Availability and Range Exclusions Among Clinical Candidates",
  s1,
  paste0(
    "Counts among the ",
    fmt_n(den),
    " clinically eligible candidates (age 50",
    EN,
    "79 years, no self-reported cardiovascular disease, ",
    "not known to be pregnant). Counts overlap and should not be summed. Out-of-range values were excluded, not truncated. ",
    "Grip includes participants without a qualifying maximal-effort trial in each hand; mobility includes questionnaire routing ",
    "below age 60 years, not only item nonresponse. eGFR indicates estimated glomerular filtration rate; and HDL, high-density lipoprotein."
  ),
  ft
)

# fixed spline bases
ss <- res("spline_specification")
s2 <- data.frame(
  Variable = ifelse(ss$variable == "logit_PREVENT", "Logit PREVENT risk", "Age, y"),
  Feature = c(knot1 = "Interior knot 1", knot2 = "Interior knot 2", boundary_low = "Lower boundary knot", boundary_high = "Upper boundary knot")[ss$role],
  Value = ifelse(ss$variable == "logit_PREVENT", fx(ss$value, 3), fx(ss$value, 0)),
  Risk = ifelse(ss$variable == "logit_PREVENT", paste0(fx(100 * plogis(ss$value), 1), "%"), "")
)
names(s2) <- c("Variable", "Basis feature", "Value", "Equivalent PREVENT risk")
register(
  "spline_specification", "Fixed Natural Spline Specification for PREVENT Risk and Age", s2,
  paste0(
    "Natural cubic splines with 3 df; knots were derived from covariate distributions in the primary cohort before any ",
    "exposure-outcome model was fitted and were reused, unchanged, in the mobility models, sensitivity analyses, and jackknife ",
    "replicates. A sensitivity analysis re-derived knots within the older cohort."
  ),
  flextable(s2) |> style_ft()
)

# selection at each exclusion gate, with crude mortality of each excluded group
sg <- read.csv(file.path(P$csv, "analysis_outputs", "selection_by_gate.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
stopifnot(sg$N[sg$group == "Included (primary cohort)"] == 2650)
s3 <- data.frame(
  Group = sg$group, N = fmt_n(sg$N), Linked = fmt_n(sg$linked),
  Deaths = ifelse(sg$linked > 0, sprintf("%s (%s)", fmt_n(sg$deaths), fx(sg$crude_mortality_pct, 1)), "NA"),
  Age = fx(sg$wt_mean_age, 1),
  Women = fx(sg$wt_female_pct, 1),
  BMI = fx(sg$wt_mean_bmi, 1),
  DM = fx(sg$wt_diabetes_pct, 1)
)
names(s3) <- c(
  "Group, in order of exclusion",
  "n",
  "Linked, n",
  "Deaths, n (% of linked)",
  "Age, y",
  "Women, %",
  "Body-mass index, kg/m\u00b2",
  "Diagnosed diabetes, %"
)
ft <- flextable(s3) |> style_ft()
ft <- fontsize(ft, size = 8.5, part = "all")
ft <- bold(ft,
  i = which(sg$group %in% c("All excluded candidates", "Included (primary cohort)")),
  part = "body"
)
ft <- hline(ft,
  i = which(sg$group == "No qualifying maximal-effort trial"),
  border = fp_border(color = "black", width = 0.5),
  part = "body"
)
register(
  "selection_by_exclusion_step",
  "Selection at Each Exclusion Gate Among Clinically Eligible Candidates",
  s3,
  paste0(
    "Clinically eligible candidates were examined participants aged 50",
    EN,
    "79 years with a positive examination weight, without self-reported cardiovascular disease, and not known to be pregnant (n=",
    fmt_n(sum(sg$N[!sg$group %in% c("All excluded candidates")])),
    "). Groups are mutually exclusive and are listed in the order in which the eligibility criteria were applied. ",
    "Deaths are counted among participants with eligible mortality linkage; percentages are crude (unweighted) proportions over each group's own follow-up and are not adjusted for age. ",
    "Age, sex, body-mass index, and diabetes are survey-weighted among participants with available values. \"Grip not tested\" includes several reasons and does not mean physical inability alone. NA, not applicable."
  ),
  ft
)

# grip disposition
gs <- res("grip_selection")
lab4 <- c(
  "bilateral valid maximal effort" = "Valid maximal-effort trial in each hand (primary cohort)",
  "one hand valid maximal effort" = "Valid maximal-effort trial in 1 hand only",
  "not tested mixed reasons" = "Not tested (any reason)",
  "tested but no qualifying effort" = "Tested, no maximal-effort trial"
)
s4 <- data.frame(
  Status = lab4[gs$grip_status],
  N = fmt_n(gs$N),
  Deaths = fmt_n(gs$deaths),
  Age = fx(gs$mean_age, 1)
)
names(s4) <- c("Grip-testing status", "n", "Deaths, n", "Mean age, y (unweighted)")
register(
  "grip_disposition",
  "Grip-Testing Disposition Among Participants Eligible Before the Grip Requirement",
  s4,
  paste0(
    "Population: the ",
    fmt_n(sum(gs$N)),
    " participants who met all eligibility criteria other than grip testing. ",
    "\"Not tested\" includes several reasons and should not be interpreted as physical inability alone."
  ),
  flextable(s4) |> style_ft()
)

# adjusted associations
fe <- res("focal_estimates")
pick <- function(m, t) {
  r <- fe[fe$model == m & fe$term == t, ]
  stopifnot(nrow(r) == 1)
  r
}
s5r <- list(
  c("Primary cohort", "Lower grip strength, per SD", "Primary: continuous strength", "lower_z"),
  c("Primary cohort", "Lower grip strength, per 5 kg", "Sensitivity: per 5 kg lower", "lower5"),
  c("Older cohort, model B", "Mobility limitation", "Mobility B: add reported limitation", "mobility"),
  c("Older cohort, model C", "Lower grip strength, per SD", "Mobility C: add measured strength", "lower_z"),
  c("Older cohort, model D", "Lower grip strength, per SD", "Mobility D: strength and limitation", "lower_z"),
  c(
    "Older cohort, model D",
    "Mobility limitation",
    "Mobility D: strength and limitation",
    "mobility"
  )
)
s5 <- do.call(rbind, lapply(s5r, function(x) {
  r <- pick(x[3], x[4])
  data.frame(
    Model = x[1], Exposure = x[2], N = fmt_n(r$N), Deaths = fmt_n(r$deaths),
    Beta = sprintf("%s (%s)", fx(r$beta, 3), fx(r$SE, 3)),
    HR = ci(r$HR, r$lower, r$upper),
    P = pfmt(r$P),
    df = r$df
  )
}))
names(s5) <- c("Model", "Exposure", "n", "Deaths", "\u03b2 (SE)", "HR (95% CI)", "P value", "df")
register(
  "adjusted_associations",
  "Adjusted Associations of Grip Strength and Mobility Limitation With All-Cause Mortality",
  s5,
  "All models adjust for natural splines of logit PREVENT risk and age, sex, survey cycle, height, and body-mass index. Grip contrasts are per 1 sex-specific SD lower strength unless labeled per 5 kg; mobility limitation is present vs absent. Models B to D use the same participants aged 60\u201379 years; model A (covariates only) has no functional-exposure estimate. SEs are design-based; CIs and P values use t distributions with the df shown (design df minus coefficients plus 1). HR indicates hazard ratio.",
  flextable(s5) |> style_ft()
)

# grip spline knots and tests
kn <- res("grip_spline_knots")
nl <- res("nonlinearity_tests")
sc <- res("strength_scaling")
s6 <- rbind(
  data.frame(
    Item = paste0("Knot at weighted ", 100 * kn$probability, "th percentile"),
    Value = sprintf(
      "%s SD (men %s kg; women %s kg)", fx(-kn$lower_z, 2),
      fx(sc$mean_kg[sc$sex == 1] - kn$lower_z * sc$sd_kg[sc$sex == 1], 1),
      fx(sc$mean_kg[sc$sex == 2] - kn$lower_z * sc$sd_kg[sc$sex == 2], 1)
    )
  ),
  data.frame(
    Item = c("Overall association, Wald F test", "Nonlinear component, Wald F test"),
    Value = sprintf(
      "F=%s (%d, %d df); P%s", fx(nl$F, 2), nl$num_df, nl$den_df,
      ifelse(nl$P < 0.001, "<0.001", paste0("=", pfmt(nl$P)))
    )
  )
)
names(s6) <- c("Specification or test", "Value")
register(
  "grip_spline_tests",
  "Restricted Cubic Spline for Grip Strength: Knots and Design-Based Tests",
  s6,
  "Three knots at the survey-weighted 10th, 50th, and 90th percentiles of standardized strength, expressed as SD from the sex-specific mean (positive values indicate stronger grip) with sex-specific kilogram equivalents. The model retained the original covariate splines. A nonsignificant nonlinear term does not establish exact linearity.",
  flextable(s6) |> style_ft()
)

# joint descriptive groups
jg <- res("joint_group_counts")
s7 <- data.frame(Group = gsub(" / ", ", ", jg$joint), N = fmt_n(jg$N), Deaths = fmt_n(jg$deaths))
names(s7) <- c("Mobility status and descriptive strength group", "n", "Deaths, n")
register(
  "joint_groups", "Descriptive Joint Mobility and Strength Groups in the Older Cohort", s7,
  "Unweighted counts among participants aged 60\u201379 years with classifiable mobility. Low strength: <27 kg in men or <16 kg in women; preserved strength means above these descriptive cutoffs. Because the low-strength groups were sparse, no adjusted categorical mortality comparison was fitted.",
  flextable(s7) |> style_ft()
)

# standardized five-year mortality, primary and alternative specifications
ar <- res("absolute_risk_sensitivities")
st <- res("support_target_comparison")
mlab <- c(
  original = "Primary model (linear strength)",
  risk_stratified = "Risk-quartile-stratified baseline hazards",
  nonlinear_grip = "Spline-modeled strength"
)
tlab <- c(
  original_support = sprintf("Supported target (n=%s)", fmt_n(st$N[st$target == "original_support"])),
  full_target = sprintf("Full primary cohort (n=%s)", fmt_n(st$N[st$target == "full_target"]))
)
s8 <- do.call(
  rbind,
  lapply(names(mlab), function(m) do.call(rbind, lapply(names(tlab), function(t) {
    z <- ar[ar$model == m & ar$target == t, ]
    e <- function(k) z[z$estimand == k, ]
    data.frame(
      Model = mlab[[m]], Target = tlab[[t]],
      Low = ci(100 * e("q25")$estimate, 100 * e("q25")$lower, 100 * e("q25")$upper),
      High = ci(100 * e("q75")$estimate, 100 * e("q75")$lower, 100 * e("q75")$upper),
      Diff = ci(
        100 * e("difference")$estimate,
        100 * e("difference")$lower,
        100 * e("difference")$upper
      )
    )
  })))
)
prim <- res("standardized_5yr_risks")
stopifnot(abs(ar$estimate[ar$model == "original" & ar$target == "original_support" & ar$estimand == "difference"] -
  prim$estimate[prim$estimand == "difference"]) < 1e-12)
names(s8) <- c(
  "Model",
  "Standardization target",
  "Lower strength, % (95% CI)",
  "Higher strength, % (95% CI)",
  "Difference, percentage points (95% CI)"
)
qq <- res("risk_strength_contrasts")
register(
  "standardized_mortality_alternatives",
  "Standardized 5-Year All-Cause Mortality at Lower and Higher Grip Strength",
  s8,
  paste0(
    "Lower and higher strength are the sex-specific weighted 25th and 75th percentiles (men ",
    fx(qq$q25[qq$sex == 1], 1),
    " and ",
    fx(qq$q75[qq$sex == 1], 1),
    " kg; women ",
    fx(qq$q25[qq$sex == 2], 1),
    " and ",
    fx(qq$q75[qq$sex == 2], 1),
    " kg). ",
    "Every model was fitted in all ",
    fmt_n(st$N[st$target == "full_target"]),
    " primary-cohort participants; the first row is the ",
    "main estimate and the others are sensitivity analyses. CIs use ",
    unique(ar$replicates),
    " stratified ",
    "jackknife replicates in which coefficients, baseline hazards, and target weights were refitted; contrast values, spline knots, ",
    "risk-quartile boundaries, and target membership were fixed. A sex-age support rule retained all participants and therefore ",
    "coincides with the full-cohort target. Full-cohort predictions extend beyond the original sex-age-BMI support rule. These are ",
    "model-based adjusted predictions, not effects of changing strength."
  ),
  flextable(s8) |> style_ft()
)

# standardization target comparison
tl9 <- c(
  original_support = "Supported standardization target",
  outside_original_support = "Outside the supported target",
  full_target = "Full primary cohort"
)
s9 <- do.call(rbind, lapply(names(tl9), function(t) {
  z <- st[st$target == t, ]
  data.frame(
    Group = tl9[[t]],
    N = fmt_n(z$N),
    Deaths = fmt_n(z$deaths),
    D5 = fmt_n(z$deaths_5),
    W = fx(100 * z$weighted_fraction, 1),
    Age = fx(z$age_weighted, 1),
    Risk = fx(100 * z$PREVENT_weighted, 1),
    Grip = fx(z$grip_weighted, 1)
  )
}))
names(s9) <- c(
  "Group",
  "n",
  "Deaths, n",
  "Deaths within 5 y, n",
  "Weighted share, %",
  "Mean age, y",
  "Mean PREVENT risk, %",
  "Mean grip strength, kg"
)
register(
  "standardization_target",
  "Participants Inside and Outside the Supported Standardization Target",
  s9,
  "Counts are unweighted; means and shares are survey weighted. Support required sex, age-decade, and body-mass-index cells with at least 30 participants in which both assigned strength values lay within the cell's 5th to 95th percentiles. All participants contributed to the fitted mortality model regardless of target membership.",
  flextable(s9) |> style_ft()
)

# planned sensitivity analyses
s10r <- list(
  c("Attained age as time scale", "Sensitivity: attained age"),
  c("Two-year landmark", "Sensitivity: two-year landmark"),
  c("Any valid hand (1-hand testing allowed)", "Sensitivity: any valid hand"),
  c("Excluding recent hand pain or stiffness", "Sensitivity: no reported hand pain"),
  c("Glycated hemoglobin-inclusive diabetes definition", "Sensitivity: HbA1c-inclusive diabetes"),
  c("Individual risk factors replace PREVENT", "Sensitivity: individual clinical factors")
)
s10 <- do.call(rbind, lapply(s10r, function(x) {
  r <- pick(x[2], "lower_z")
  data.frame(
    Analysis = x[1],
    N = fmt_n(r$N),
    Deaths = fmt_n(r$deaths),
    HR = ci(r$HR, r$lower, r$upper),
    P = pfmt(r$P),
    df = r$df
  )
}))
names(s10) <- c("Analysis", "n", "Deaths", "HR (95% CI)", "P value", "df")
register(
  "planned_sensitivity", "Sensitivity Analyses for Grip Strength", s10,
  "HR per 1 sex-specific SD lower grip strength, with the original sex-specific scaling. The attained-age model used age at examination as delayed entry; the landmark analysis restarted follow-up at 2 years among participants still under observation; recent hand pain or stiffness refers to either hand in the previous 7 days; the individual-factor model replaced the PREVENT and age splines with linear age and the PREVENT inputs. The per-5-kg contrast is shown in Table S5. These correlated estimates are not independent replications. HR indicates hazard ratio.",
  flextable(s10) |> style_ft()
)

# proportional-hazards diagnostics
ph <- res("PH_joint_terms")
tlab11 <- c(
  lower_z = "Lower grip strength",
  mobility = "Mobility limitation",
  "PREVENT risk spline block" = "PREVENT risk spline (3 df)",
  "Age spline block" = "Age spline (3 df)",
  sex_factor = "Sex",
  cycle = "Survey cycle",
  height = "Height",
  bmi = "Body-mass index",
  GLOBAL = "Global test"
)
s11 <- rbind(
  data.frame(Term = "Primary cohort", Chi = "", df = "", P = ""),
  with(ph[ph$cohort == "core", ], data.frame(Term = tlab11[term], Chi = fx(chisq, 2), df = as.character(df), P = pfmt(p))),
  data.frame(Term = "Older cohort, model D", Chi = "", df = "", P = ""),
  with(
    ph[ph$cohort == "mobility_cohort", ],
    data.frame(Term = tlab11[term], Chi = fx(chisq, 2), df = as.character(df), P = pfmt(p))
  )
)
names(s11) <- c("Term", "\u03c7\u00b2", "df", "P value")
hdr11 <- which(s11$df == "")
ft <- flextable(s11) |> style_ft()
ft <- section_rows(ft, hdr11)
ft <- indent_rows(ft, setdiff(seq_len(nrow(s11)), hdr11))
register(
  "PH_tests", "Exploratory Proportional-Hazards Diagnostics", s11,
  "Scaled Schoenfeld residual tests (Kaplan-Meier time transform), with each spline tested as a single multi-df term. These conventional diagnostics do not account for the complex survey design and are exploratory. Smoothed residual plots are shown in Figure S4.",
  ft
)

# incremental mortality prediction
pm <- read.csv(file.path(P$csv, "analysis_outputs", "prediction_metrics.csv"),
  stringsAsFactors = FALSE
)
fe5 <- res("focal_estimates")
pm_ <- function(co, me, mo) {
  r <- pm[pm$cohort == co & pm$metric == me & pm$model == mo, ]
  stopifnot(nrow(r) <= 1)
  r
}
cell5 <- function(co, me, mo, d = 3, corrected = TRUE) {
  r <- pm_(co, me, mo)
  if (!nrow(r)) return("NA")
  if (me == "calibration_slope") return(fx(r$corrected, 2))
  fs <- function(x) {
    v <- fx(x, d)
    v <- sub("^-(0\\.0+)$", "\\1", v)
    gsub("-", "\u2212", v)
  } # true minus; no "-0.000"
  out <- sprintf("%s\n(%s to %s)", fs(r$apparent), fs(r$lower), fs(r$upper))
  if (corrected) out <- sprintf("%s\n[%s]", out, fs(r$corrected))
  out
}
mods5 <- list(
  primary = c(
    PREVENT_only = "PREVENT spline only (reference)", M0 = "Base model, no functional measure",
    M1 = "Base model + grip strength"
  ),
  mobility = c(
    M0 = "Base model, no functional measure", Mobility = "+ Mobility limitation",
    Grip = "+ Grip strength", Both = "+ Grip strength and mobility limitation"
  )
)
fp5 <- fe5[fe5$model == "Primary: continuous strength" & fe5$term == "lower_z", ]
fm5 <- fe5[fe5$model == "Mobility D: strength and limitation" & fe5$term == "lower_z", ]
hdr5 <- c(
  primary = sprintf("Primary cohort, age 50%s79 y (n=%s; %s deaths)", EN, fmt_n(fp5$N), fmt_n(fp5$deaths)),
  mobility = sprintf(
    "Older cohort, age 60%s79 y (n=%s; %s deaths)",
    EN,
    fmt_n(fm5$N),
    fmt_n(fm5$deaths)
  )
)
s5p <- do.call(rbind, lapply(names(mods5), function(co) rbind(
  data.frame(Model = hdr5[[co]], C = "", dC = "", AUC = "", dAUC = "", NRI = "", IDI = "", Slope = ""),
  do.call(rbind, lapply(names(mods5[[co]]), function(mo) {
    inc <- !(mo %in% c("M0", "PREVENT_only"))
    data.frame(
      Model = mods5[[co]][[mo]],
      C = cell5(co, "C", mo),
      dC = if (inc) cell5(co, "dC", mo) else "NA",
      AUC = cell5(co, "AUC5", mo), dAUC = if (inc) cell5(co, "dAUC5", mo) else "NA",
      NRI = if (inc) cell5(co, "NRI", mo, 2, FALSE) else "NA",
      IDI = if (inc) cell5(co, "IDI", mo, 3, FALSE) else "NA",
      Slope = cell5(co, "calibration_slope", mo)
    )
  }))
)))
names(s5p) <- c(
  "Model",
  "Harrell C",
  "Change in C vs base",
  "5-y AUC",
  "Change in 5-y AUC vs base",
  "Continuous NRI at 5 y",
  "IDI at 5 y",
  "Calibration slope"
)
hdr5i <- which(s5p$Model %in% hdr5) # section rows (columns already renamed)
ft <- flextable(s5p) |> style_ft()
ft <- fontsize(ft, size = 8.5, part = "all")
ft <- section_rows(ft, hdr5i)
ft <- indent_rows(ft, setdiff(seq_len(nrow(s5p)), hdr5i))
ft <- compose(ft,
  part = "header",
  i = 1,
  j = 8,
  value = as_paragraph("Calibration slope", as_sup("a"))
)
register(
  "incremental_prediction",
  "Incremental All-Cause Mortality Prediction With Grip Strength and Mobility Limitation",
  s5p,
  paste0(
    "Exploratory secondary analyses. The base model included the PREVENT and age splines, sex, survey cycle, height, and body-mass index; the mobility-cohort models were fitted in identical participants. ",
    "Values are examination-weighted apparent estimates with 95% CIs from ",
    unique(na.omit(pm$jk_replicates)),
    " stratified jackknife replicate refits (t distribution, 22 df); values in brackets are optimism-corrected using ",
    unique(pm$bootstrap_replicates),
    " survey (Rao",
    EN,
    "Wu rescaling) bootstrap replicates. The 5-year AUC compares deaths within 5 years with participants known to be alive at 5 years; participants censored before 5 years are excluded from the 5-year metrics. ",
    "Continuous NRI and IDI compare 5-year predicted risks with those of the base model. ^a^Optimism-corrected calibration slope (apparent value 1 by construction). ",
    "These metrics address all-cause mortality, not cardiovascular events. AUC indicates area under the receiver operating characteristic curve; IDI, integrated discrimination improvement; NA, not applicable; and NRI, net reclassification improvement."
  ),
  ft
)

# physical function within guideline PREVENT categories
gb <- read.csv(file.path(P$csv, "analysis_outputs", "guideline_burden.csv"),
  stringsAsFactors = FALSE
)
qq7 <- res("risk_strength_contrasts")
meas6 <- c(
  low_grip = "Low grip strength^a^ (age 50\u201379 y)",
  mobility = "Mobility limitation (age 60\u201379 y)",
  adl = "ADL limitation (age 60\u201379 y)",
  iadl = "IADL limitation (age 60\u201379 y)",
  anyfx = "Any of the 3 limitations (age 60\u201379 y)"
)
gb_ <- function(m, cv, lv) {
  r <- gb[gb$measure == m & gb$category_var == cv & gb$category == lv, ]
  stopifnot(nrow(r) == 1)
  r
}
cell6 <- function(r, pop = TRUE) paste0(
  sprintf("%s/%s\n%s (%s%s%s)", fmt_n(r$cases), fmt_n(r$N), fx(100 * r$prevalence, 1), fx(100 * r$prev_lower, 1), EN, fx(100 * r$prev_upper, 1)),
  if (pop) sprintf(
    "\n%s (%s%s%s)",
    fx(r$pop_millions, 1),
    fx(r$pop_lower, 1),
    EN,
    fx(r$pop_upper, 1)
  ) else ""
)
rows6 <- list(
  list("PREVENT-ASCVD 10-year risk", NULL),
  list("<3%", c("asc_cat", "<3%")),
  list("3% to <5%", c("asc_cat", "3% to <5%")),
  list("5% to <10%", c("asc_cat", "5% to <10%")),
  list("\u226510%", c("asc_cat", ">=10%")),
  list("PREVENT total CVD 10-year risk", NULL),
  list("<7.5%", c("tot75", "0")),
  list("\u22657.5%", c("tot75", "1")),
  list("Share of adults with the measure who are below the treatment threshold", NULL),
  list("PREVENT-ASCVD <5% (2026 dyslipidemia guideline)", c("share_below_asc5", "below threshold"), FALSE),
  list(
    "PREVENT total CVD <7.5% (2025 blood pressure guideline)",
    c("share_below_tot75", "below threshold"),
    FALSE
  )
)
s6 <- do.call(rbind, lapply(rows6, function(z) {
  if (is.null(z[[2]])) return(as.data.frame(as.list(setNames(c(z[[1]], rep("", length(meas6))), c("Category", names(meas6)))),
    stringsAsFactors = FALSE
  ))
  pop <- if (length(z) >= 3) z[[3]] else TRUE
  as.data.frame(as.list(setNames(c(z[[1]], vapply(names(meas6), function(m) cell6(gb_(m, z[[2]][1], z[[2]][2]), pop), "")), c("Category", names(meas6)))),
    stringsAsFactors = FALSE
  )
}))
hdr6 <- which(vapply(rows6, function(z) is.null(z[[2]]), TRUE))
names(s6) <- c("PREVENT category", strip_sup(meas6))
ft <- flextable(s6) |> style_ft()
ft <- fontsize(ft, size = 8, part = "all")
ft <- section_rows(ft, hdr6)
ft <- indent_rows(ft, setdiff(seq_len(nrow(s6)), hdr6))
ft <- compose(ft,
  part = "header",
  i = 1,
  j = 2,
  value = as_paragraph("Low grip strength", as_sup("a"), " (age 50\u201379 y)")
)
register(
  "function_within_PREVENT_categories",
  "Low Grip Strength and Functional Limitation Within Guideline PREVENT Risk Categories",
  s6,
  paste0(
    "Each cell shows participants with the measure/participants in the category (unweighted), the survey-weighted prevalence, % (95% CI), and the weighted number of US adults, millions (95% CI); ",
    "in the last section, the weighted percentage of adults with the measure whose PREVENT estimate is below the treatment threshold. Population numbers are average annual estimates for 2011",
    EN,
    "2014 ",
    "and refer to adults who met the cohort criteria (no self-reported cardiovascular disease, complete PREVENT inputs, and bilateral grip testing), not to all US adults of these ages. ",
    "Reported limitations were assessed only at ages 60",
    EN,
    "79 years. ^a^At or below the sex-specific weighted 25th percentile of the primary cohort (",
    fx(qq7$q25[qq7$sex == 1], 1),
    " kg in men and ",
    fx(qq7$q25[qq7$sex == 2], 1),
    " kg in women). ADL indicates activities of daily living; ASCVD, atherosclerotic cardiovascular disease; CVD, cardiovascular disease; and IADL, instrumental activities of daily living."
  ),
  ft
)

# five-year mortality by function within treatment thresholds
gd <- read.csv(file.path(P$csv, "analysis_outputs", "guideline_discordance.csv"),
  stringsAsFactors = FALSE
)
ga <- read.csv(file.path(P$csv, "analysis_outputs", "guideline_discordance_adjusted_hr.csv"),
  stringsAsFactors = FALSE
)
mfx <- function(x) {
  v <- fx(x, 1)
  v <- sub("^-(0\\.0)$", "\\1", v)
  gsub("-", "\u2212", v)
} # true minus; no "-0.0"
pct_ci <- function(r) sprintf(
  "%s (%s to %s)",
  mfx(100 * r$estimate),
  mfx(100 * r$lower),
  mfx(100 * r$upper)
)
sec7 <- list(
  list("grip, ASCVD 5%", "Grip strength, age 50\u201379 y; PREVENT-ASCVD threshold 5%", "asc5", "low grip", "preserved grip", "PREVENT-ASCVD"),
  list("grip, total CVD 7.5%", "Grip strength, age 50\u201379 y; PREVENT total CVD threshold 7.5%", "tot75", "low grip", "preserved grip", "PREVENT total CVD"),
  list("reported function, ASCVD 5%", "Reported function, age 60\u201379 y; PREVENT-ASCVD threshold 5%", NA, "any limitation", "no limitation", "PREVENT-ASCVD"),
  list(
    "reported function, total CVD 7.5%",
    "Reported function, age 60\u201379 y; PREVENT total CVD threshold 7.5%",
    NA,
    "any limitation",
    "no limitation",
    "PREVENT total CVD"
  )
)
s7 <- do.call(rbind, lapply(sec7, function(z) {
  g <- gd[gd$analysis == z[[1]], ]
  thr <- if (grepl("ASCVD", z[[1]])) "5%" else "7.5%"
  gl <- c(
    below_preserved = sprintf("Below %s, %s", thr, z[[5]]),
    below_low = sprintf("Below %s, %s", thr, z[[4]]),
    above_preserved = sprintf("At or above %s, %s", thr, z[[5]]),
    above_low = sprintf("At or above %s, %s", thr, z[[4]])
  )
  cl <- c(
    rd_below = sprintf("Difference below %s (%s minus %s)", thr, z[[4]], z[[5]]),
    rd_above = sprintf("Difference at or above %s", thr),
    additive_interaction = "Difference between the 2 differences",
    discordance = sprintf("Below %s with %s minus at or above %s with %s", thr, z[[4]], thr, z[[5]])
  )
  rr <- function(e) {
    r <- g[g$estimand == e, ]
    stopifnot(nrow(r) == 1)
    r
  }
  out <- data.frame(
    Group = z[[2]],
    N = "",
    Deaths = "",
    Age = "",
    Pop = "",
    Mort = "",
    stringsAsFactors = FALSE
  )
  for (e in names(gl)) {
    r <- rr(e)
    out <- rbind(out, data.frame(
      Group = gl[[e]],
      N = fmt_n(r$N),
      Deaths = sprintf("%s (%s/%s)", fmt_n(r$deaths), fmt_n(r$cvd_deaths), fmt_n(r$other_deaths)),
      Age = fx(r$wt_mean_age, 1),
      Pop = fx(r$pop_millions, 1),
      Mort = if (r$reportable) pct_ci(r) else "Not estimated (<10 deaths)"
    ))
  }
  for (e in names(cl)) {
    r <- rr(e)
    out <- rbind(out, data.frame(
      Group = cl[[e]], N = "", Deaths = "", Age = "", Pop = "",
      Mort = if (r$reportable) paste(pct_ci(r), "points") else "Not estimated"
    ))
  }
  if (!is.na(z[[3]])) {
    h <- ga[ga$threshold == z[[3]], ]
    hb <- h[h$term == "low_grip_below", ]
    ha <- h[h$term == "low_grip_above", ]
    hi <- h[h$term == "interaction", ]
    out <- rbind(out, data.frame(
      Group = sprintf("Adjusted HR for low grip: below %s; at or above %s; P for product term", thr, thr),
      N = "",
      Deaths = "",
      Age = "",
      Pop = "",
      Mort = sprintf(
        "%s; %s; P=%s",
        ci(hb$HR, hb$lower, hb$upper),
        ci(ha$HR, ha$lower, ha$upper),
        pfmt(hi$P)
      )
    ))
  }
  out
}))
names(s7) <- c(
  "Group",
  "n",
  "Deaths, n (cardiovascular/other)",
  "Weighted mean age, y",
  "US adults, millions",
  "5-year mortality, % (95% CI)"
)
hdr7 <- which(s7$Group %in% vapply(sec7, function(z) z[[2]], ""))
ft <- flextable(s7) |> style_ft()
ft <- fontsize(ft, size = 8, part = "all")
ft <- section_rows(ft, hdr7)
ft <- indent_rows(ft, setdiff(seq_len(nrow(s7)), hdr7))
for (i in grep("^Adjusted HR", s7$Group)) ft <- compose(ft,
  i = i,
  j = 1,
  part = "body",
  value = as_paragraph(s7$Group[i], as_sup("a"))
)
register(
  "mortality_within_treatment_thresholds",
  "Five-Year All-Cause Mortality by Physical Function Within Guideline PREVENT Treatment Thresholds",
  s7,
  paste0(
    "Five-year mortality is the survey-weighted Kaplan",
    EN,
    "Meier cumulative mortality (crude, not adjusted for age), with 95% CIs from ",
    61,
    " stratified jackknife replicates (t distribution, ",
    unique(gd$df),
    " df); differences are in percentage points. Mortality is not estimated for groups with fewer than 10 deaths. Low grip strength was at or below the sex-specific weighted 25th percentile; ",
    "reported function was limitation in mobility, ADL, or IADL. Cardiovascular deaths were those assigned to heart disease or cerebrovascular disease. ",
    "^a^Survey-weighted Cox model in the primary cohort with low grip strength, the threshold indicator, and their product, adjusted for the PREVENT total-CVD and age splines, sex, survey cycle, height, and body-mass index. ",
    "ADL indicates activities of daily living; ASCVD, atherosclerotic cardiovascular disease; CVD, cardiovascular disease; HR, hazard ratio; and IADL, instrumental activities of daily living."
  ),
  ft
)

# post-diagnostic models
pd <- res("post_diagnostic_estimates")
cn <- res("post_diagnostic_conditioning")
md <- res("model_diagnostics")
s12r <- list(
  c("Risk-quartile-stratified baseline hazards, primary cohort", "Post-diagnostic: risk-stratified core", "lower_z", "Lower grip strength"),
  c("Risk-quartile-stratified baseline hazards, older cohort", "Post-diagnostic: risk-stratified mobility_cohort", "lower_z", "Lower grip strength"),
  c("Risk-quartile-stratified baseline hazards, older cohort", "Post-diagnostic: risk-stratified mobility_cohort", "mobility", "Mobility limitation"),
  c("Mobility-specific spline knots", "Post-diagnostic: mobility-specific spline knots", "lower_z", "Lower grip strength"),
  c(
    "Mobility-specific spline knots",
    "Post-diagnostic: mobility-specific spline knots",
    "mobility",
    "Mobility limitation"
  )
)
s12 <- do.call(rbind, lapply(s12r, function(x) {
  r <- pd[pd$model == x[2] & pd$term == x[3], ]
  data.frame(
    Model = x[1],
    Exposure = x[4],
    N = fmt_n(r$N),
    Deaths = fmt_n(r$deaths),
    HR = ci(r$HR, r$lower, r$upper),
    CN = fx(cn$condition_number[cn$model == x[2]], 1)
  )
}))
names(s12) <- c("Model", "Exposure", "n", "Deaths", "HR (95% CI)", "Condition number")
register(
  "post_diagnostic", "Sensitivity Models for Proportional Hazards and Conditioning", s12,
  paste0(
    "Stratified models allowed separate baseline hazards across PREVENT-risk quartiles while retaining continuous risk adjustment. Mobility-specific knots re-derived the PREVENT and age splines within the older cohort to address numerical conditioning; the scaled design-matrix condition number was ",
    fx(md$condition_number[md$model == "Mobility D: strength and limitation"], 1),
    " with the original bases and ",
    fx(md$condition_number[md$model == "Primary: continuous strength"], 1),
    " in the primary model. The original models remain the principal estimates. HR indicates hazard ratio."
  ),
  flextable(s12) |> style_ft()
)

# eligibility audit
sl <- res("selection_counts")
lab13 <- c(
  "Original core" = "Primary cohort", "BMI-only additions" = "Added by relaxing only the BMI range",
  "BMI-relaxed sensitivity" = "BMI-relaxed sensitivity cohort",
  "Clinical candidates with missing score inputs and otherwise valid grip/linkage/height" = "Missing score inputs, otherwise eligible (not analyzed)"
)
s13 <- data.frame(Population = lab13[sl$group], N = fmt_n(sl$N), Deaths = fmt_n(sl$deaths))
names(s13) <- c("Population", "n", "Deaths, n")
register(
  "eligibility_audit", "Eligibility Audit", s13,
  "Relaxing the body-mass-index (BMI) range retained all other eligibility rules, measured BMI as a covariate, and the original spline bases and strength scaling; the base total-CVD PREVENT equation has zero coefficients for BMI terms. Participants with missing score inputs could not be assessed against input ranges and were not imputed into the cohort.",
  flextable(s13) |> style_ft()
)

# sensitivity estimates
rv <- res("additional_sensitivity_estimates")
mres <- res("imputation_result")
s14r <- list(
  c("BMI range restriction relaxed", "BMI range relaxed"),
  c("Individual risk factors with age spline", "Clinical factors with matched flexible age"),
  c("Three maximal-effort trials per hand", "Three maximal-effort trials per hand"),
  c("Base PREVENT, biomarker-eligible sample", "Base PREVENT common biomarker sample"),
  c("Full PREVENT with UACR and HbA1c, same sample", "Full PREVENT with UACR and HbA1c"),
  c("Base adjustment, socioeconomic complete cases", "Base adjustment in expanded complete cases"),
  c(
    "Plus education, income, and cancer history, same sample",
    "Additional education income and cancer"
  )
)
s14 <- do.call(rbind, lapply(s14r, function(x) {
  r <- rv[rv$model == x[2] & rv$term == "lower_z", ]
  data.frame(
    Analysis = x[1],
    N = fmt_n(r$N),
    Deaths = fmt_n(r$deaths),
    HR = ci(r$HR, r$lower, r$upper),
    P = pfmt(r$P),
    df = fx(r$df, 0)
  )
}))
s14 <- rbind(s14, data.frame(
  Analysis = sprintf("Plus education, income, and cancer history; income imputed (m=%d)", mres$imputations),
  N = fmt_n(mres$N),
  Deaths = fmt_n(mres$deaths),
  HR = ci(mres$HR, mres$lower, mres$upper),
  P = pfmt(mres$P),
  df = fx(mres$df, 1)
))
names(s14) <- c("Analysis", "n", "Deaths", "HR (95% CI)", "P value", "df")
register(
  "additional_sensitivity", "Additional Sensitivity Analyses for Grip Strength", s14,
  paste0(
    "HR per 1 original sex-specific SD lower grip strength. Paired base and expanded models use identical participants. Full PREVENT used the bundled preventr 0.12.0 coefficients with the specified missing social-deprivation term and its own spline knots. Income imputation used predictive mean matching; Rubin's rules with a finite complete-data df gave the df shown (fraction of missing information for grip, ",
    fx(mres$missing_information, 2),
    "). Imputation did not expand score or grip eligibility, and 1 participant with unresolved cancer history was excluded. Sensitivity analyses are exploratory. HbA1c indicates glycated hemoglobin; HR, hazard ratio; and UACR, urine albumin-to-creatinine ratio."
  ),
  flextable(s14) |> style_ft()
)

saveRDS(TABLES, file.path(BUILD, "tables_registry.rds"))
message("Tables written: ", paste(names(TABLES), collapse = ", "))
