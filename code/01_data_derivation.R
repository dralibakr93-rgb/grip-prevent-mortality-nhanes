# Stage 1. Builds the participant-level dataset from the raw NHANES 2011-2014 and linked mortality
# files: derived variables, PREVENT base risk, eligibility flags, cohort flow (Figure 1), and the
# sex-specific mean and SD of grip strength used to scale the exposure.
# Runs code/analysis/01_build_cohort.R.

source("code/00_setup.R")
stage_build()
run_analysis("01_build_cohort.R")
