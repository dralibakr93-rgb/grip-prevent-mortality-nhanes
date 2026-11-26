# Stage 4. Standardized 5-year mortality at the 25th vs 75th sex-specific percentile of grip strength,
# with jackknife (JKn) confidence intervals; risk-stratified and cohort-specific spline models; and
# selection summaries (code/analysis/03_standardized_mortality.R). Then additional sensitivity analyses:
# BMI range relaxed, grip-strength spline, flexible age, 3 trials per hand, full PREVENT, education,
# income and cancer history, and income imputation (code/analysis/04_additional_sensitivity.R).

source("code/00_setup.R")
run_analysis("03_standardized_mortality.R")
run_analysis("04_additional_sensitivity.R")
