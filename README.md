# Grip strength, reported function, and all-cause mortality within guideline PREVENT risk categories: NHANES 2011-2014

Analysis code and data for a survey-weighted cohort study of US adults aged 50 to 79 years without
cardiovascular disease. The study relates measured grip strength and reported mobility and daily
function to all-cause mortality through 2019, after adjustment for PREVENT-estimated cardiovascular
risk, and describes low grip strength and functional limitation within the PREVENT-ASCVD categories
of the 2026 ACC/AHA dyslipidemia guideline.

This submission version contains the R analysis scripts, variable dictionaries, aggregate results,
and instructions for downloading the exact public-use input files from NCHS. The scripts build
the participant-level datasets locally. The manuscript’s principal estimates and their output files
are listed below.

## Repository structure

| Path | Contents |
|---|---|
| `run_all.R` | Runs all stages in order. |
| `code/` | Stage scripts `00` to `08`, run in numerical order by `run_all.R`. |
| `code/analysis/` | Core estimation scripts called by stages 1, 3, 4, and 5. |
| `code/independent_check.R` | Separate re-derivation of the cohort, PREVENT scores, and primary estimate. |
| `data/raw/` | Local download destination for NHANES 2011–2014 and the 2019 public-use linked mortality files. |
| `data/raw_data_manifest.txt` | Size, SHA-256 checksum, and NCHS download address of every raw file. |
| `data/download_raw_data.R` | Checks raw files against the manifest and downloads them when run with `--download`. |
| `data/reference/` | PREVENT coefficient tables from the preventr R package (version 0.12.0). |
| `output/analytic_dataset/` | Local output for the derived design-sample data. |
| `data/dictionary/` | Definitions of derived analytic variables and their NHANES sources. |
| `results/` | The results as reported in the manuscript (see `results/README.md`). |
| `environment/` | R package versions, R session information, and the bundled preventr 0.12.0 package. |

`data/README.md` describes, step by step, how the analytic dataset is built from the raw files.

## Reproducing the analysis

Requirements: R 4.6.0 with the packages listed in `environment/R_package_versions.csv` (survey 4.5,
survival 3.8-6, mice 3.19.0, Hmisc 5.2-5, and others). The PREVENT package preventr 0.12.0 is
bundled in `environment/R_library/` and is loaded from there.

From the repository root:

```
Rscript data/download_raw_data.R --download   # download and verify the 35 public-use inputs
Rscript run_all.R                    # 2. run all stages (about 5 to 10 minutes)
```

All outputs are written to `output/` (not part of the repository) and intermediate files to
`_build/`. Stage 8 compares every reproduced analysis output with the corresponding file in
`results/analysis_outputs/` and writes the comparison to
`output/qc/comparison_with_reported_results.csv`. The comparison covers 67 archived CSVs and the recorded SHA-256 of the participant-level guideline file. Other versions of R or of the packages may give small numerical
differences.

Each stage can also be run on its own, in order, for example `Rscript code/03_primary_analysis.R`.

## Analysis stages

| Stage | Script | What it does | Main outputs (`output/csv/analysis_outputs/`) |
|---|---|---|---|
| 0 | `code/00_setup.R` | Paths, package checks, staging of the raw files | - |
| 1 | `code/01_data_derivation.R` | Derived variables, PREVENT base risk, eligibility, cohort flow, grip-strength scaling | `cohort_flow.csv`, `strength_scaling.csv` |
| 2 | `code/02_cohort.R` | Survey-weighted characteristics (Table 1); export of the analytic dataset | `table1_characteristics_long.csv` |
| 3 | `code/03_primary_analysis.R` | Survey-weighted Cox models: primary model, older-cohort models A-D, sensitivity models, Schoenfeld tests | `focal_estimates.csv`, `model_estimates_all.csv`, `PH_exploratory.csv` |
| 4 | `code/04_sensitivity_analysis.R` | Standardized 5-year mortality with jackknife CIs; additional sensitivity analyses; income imputation | `standardized_5yr_risks.csv`, `additional_sensitivity_estimates.csv`, `imputation_result.csv` |
| 5 | `code/05_standardization.R` | Standardized mortality under alternative targets and models | `absolute_risk_sensitivities.csv` |
| 5b | `code/05b_standardized_curves.R` | Standardized mortality over follow-up; Kaplan-Meier curves by grip third | `standardized_mortality_curves.csv`, `km_by_grip_tertile.csv` |
| 5c | `code/05c_function_and_prediction.R` | ADL and IADL limitation; prior hospitalization; incremental prediction | `functional_disability_estimates.csv`, `prediction_metrics.csv` |
| 5d | `code/05d_guideline_categories.R` | Function within guideline PREVENT categories and thresholds; kidney term; alternative scores; selection | `guideline_burden.csv`, `guideline_discordance.csv`, `kidney_term_by_grip_quartile.csv`, `selection_by_gate.csv` |
| 6 | `code/06_figures.R` | Figures 1-4 and S1-S9 | `output/figures/` |
| 7 | `code/07_tables.R` | Table 1 and the formatted values behind the supplementary figures | `output/tables/`, `output/csv/` |
| 8 | `code/08_reproducibility_check.R` | Comparison with `results/`; Table 1 check; independent re-derivation | `output/qc/` |

## Where the main estimates come from

| Estimate in the manuscript | File in `results/analysis_outputs/` | Row |
|---|---|---|
| Cohort flow: 19,931 participants, 2,650 in the primary cohort (192 deaths), 1,512 in the older cohort (148 deaths) | `cohort_flow.csv` | stages `all_respondents`, `core`, `mobility_age60_79` |
| HR per SD lower grip strength, 1.72 (95% CI, 1.32-2.23); per 5 kg, 1.51 (1.25-1.82) | `focal_estimates.csv` | `Primary: continuous strength`; `Sensitivity: per 5 kg lower` |
| Ages 60-79, joint model: grip strength HR 1.45 (1.13-1.87), mobility limitation HR 2.37 (1.52-3.71) | `focal_estimates.csv` | `Mobility D: strength and limitation` |
| Standardized 5-year mortality 3.69% vs 1.90%; difference 1.79 percentage points (0.86-2.72) | `standardized_5yr_risks.csv` | `q25_risk`, `q75_risk`, `difference` |
| Change in C statistic: mobility limitation 0.026 (0.002 to 0.049); grip strength 0.014 (-0.007 to 0.034); both 0.034 (0.013 to 0.055) | `prediction_metrics.csv` | cohort `mobility`, metric `dC` |
| ADL limitation, OR per SD lower grip strength 1.58 (1.28-1.96) | `functional_disability_estimates.csv` | outcome `adl` |
| Low grip strength 16% (PREVENT-ASCVD <3%) to 46% (>=10%); 42.1% of adults with low grip strength below the 5% threshold | `guideline_burden.csv` | measure `low_grip`, `asc_cat` and `share_below_asc5` |
| 5-year mortality by grip strength below and above the statin threshold | `guideline_discordance.csv`, `guideline_discordance_adjusted_hr.csv` | analysis `grip, ASCVD 5%` |
| Kidney term by grip-strength quartile | `kidney_term_by_grip_quartile.csv` | all rows |

## Statistical methods in brief

All estimates use the NHANES survey design: the pooled 4-year examination weight (WTMEC2YR/2),
masked variance strata and PSUs, and the full design with the cohort as a domain. Cox models were
fitted with `survey::svycoxph`; confidence intervals and P values use t distributions with the design
degrees of freedom minus the number of coefficients plus 1. Models adjust for natural cubic splines
(3 df) of logit PREVENT 10-year total cardiovascular risk and of age, sex, survey cycle, height, and
body-mass index. Grip strength is expressed per sex-specific SD lower strength. Standardized
mortality and prediction metrics use 61 stratified jackknife (JKn) replicates; optimism correction
uses 200 survey bootstrap replicates. Random seeds are fixed in the scripts.

## Data use

NHANES and the NCHS linked mortality files are public-use data released by the National Center for
Health Statistics. They contain no direct identifiers. Users must follow the NCHS data use
restrictions, including the prohibition on attempting to identify any participant.

For the calculator cross-checks, `data/reference/calculator_check_*.csv` contains the original
reference risks in draw order for the two seeded samples of 100 profiles. The tests recompute
the risks and compare them with these values. The public release does not redistribute the AHA
web calculator source. The independent cohort and primary estimate check also runs against preventr.

The drug classifications originate from the Multum Lexicon. See the [NCHS drug information codebook](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1988/DataFiles/RXQ_DRUG.htm).

The bundled preventr package and its coefficient tables are distributed under the MIT license.
See [the preventr license](environment/preventr-MIT.txt).

## Citation

If you use this code, please cite the accompanying article (citation to be added on publication).
