# Data

## Sources

All inputs are public-use files released by the National Center for Health Statistics (NCHS).
Download them into `raw/` with `Rscript data/download_raw_data.R --download` from the repository root. The file `raw_data_manifest.txt` lists the size, SHA-256 checksum,
and download address of each file.

| Folder in `raw/` | Files | Content |
|---|---|---|
| `2011_2012/` | DEMO_G, BMX_G, BPX_G, BPQ_G, SMQ_G, DIQ_G, MCQ_G, TCHOL_G, HDL_G, BIOPRO_G, GHB_G, MGX_G, PFQ_G, RXQ_RX_G, ALB_CR_G, HUQ_G | NHANES 2011-2012 examination, laboratory, and questionnaire files (XPT) |
| `2013_2014/` | the same components with suffix _H | NHANES 2013-2014 |
| `Reference/` | RXQ_DRUG | Drug information file; drug class 173 identifies statins |
| `Mortality_2019_LMF/` | NHANES_2011_2012_MORT_2019_PUBLIC.dat, NHANES_2013_2014_MORT_2019_PUBLIC.dat | 2019 public-use linked mortality files, follow-up through December 31, 2019 |

`reference/` holds the PREVENT base and full 10-year coefficient tables and the term-preparation
code of the preventr R package, version 0.12.0. NCHS publishes the [cycle-specific codebooks](https://wwwn.cdc.gov/nchs/nhanes/).

## From the raw files to the analytic dataset

The steps below are carried out by `code/analysis/01_build_cohort.R` (stage 1) and
`code/02_cohort.R` (stage 2). The two cycles are processed separately and then combined.

1. **Merge.** For each cycle, the examination and questionnaire files are merged with the
   demographic file by participant identifier (SEQN). Each file is checked for duplicate SEQN.
2. **Blood pressure.** Systolic pressure is the mean of the available readings (BPXSY1-4),
   excluding the first reading when more than 1 is available.
3. **Kidney function.** eGFR is computed from serum creatinine (LBXSCR) with the 2021 CKD-EPI
   creatinine equation without race.
4. **Risk factors.** Diabetes: reported diagnosis (DIQ010), insulin (DIQ050), or diabetes pills
   (DIQ070). Current smoking: at least 100 lifetime cigarettes (SMQ020) and smoking every day or
   some days (SMQ040). Antihypertensive treatment: currently taking prescribed medicine for high
   blood pressure (BPQ050A). Statin use: any statin by drug name in the prescription file, or drug
   class 173 in the drug information file.
5. **Cardiovascular disease.** Reported congestive heart failure, coronary heart disease, angina,
   heart attack, or stroke (MCQ160B-F). Absence requires a "no" answer to all 5 questions.
6. **Grip strength.** Up to 3 trials per hand (MGXH1T1-3, MGXH2T1-3). A trial is used only when
   its effort code indicates maximal effort. At least 1 such trial in each hand is required, and
   grip strength is the highest value over both hands, in kg.
7. **Mobility limitation** (asked at ages 60 years and older). Present when the participant needs
   special equipment to walk (PFQ054) or reports some difficulty, much difficulty, or inability to
   walk a quarter mile (PFQ061B) or climb 10 steps (PFQ061C). Absent when all 3 answers are negative.
8. **Daily activities.** ADL limitation: any difficulty walking between rooms, getting in and out
   of bed, using a fork or knife or drinking from a cup, or dressing (PFQ061H, J, K, L). IADL
   limitation: any difficulty managing money, doing house chores, or preparing meals
   (PFQ061A, F, G). Absent only when every item is answered "no difficulty".
9. **PREVENT risk.** Base PREVENT 10-year total cardiovascular, ASCVD, and heart failure risks are
   computed by sex from age, total and HDL cholesterol, systolic pressure, eGFR, BMI, diabetes,
   current smoking, antihypertensive treatment, and statin use, with the coefficients in
   `reference/`. The computed scores are checked against archived outputs from an independent transcription
   of the AHA calculator and against `preventr::estimate_risk()`.
10. **Mortality.** Participants with eligible linkage and known vital status are followed from the
    examination to death or December 31, 2019. Follow-up time is months divided by 12; a value of
    0 months is set to 0.5 months. Deaths from heart disease or cerebrovascular disease (underlying
    cause codes 1 and 5) are classed as cardiovascular.
11. **Survey design.** The pooled examination weight is `w = WTMEC2YR / 2`; strata are `SDMVSTRA`
    and PSUs `SDMVPSU`.

## Eligibility

Counts are participants (deaths), from `results/analysis_outputs/cohort_flow.csv`.

| Step | Criterion | Remaining |
|---|---|---|
| 1 | NHANES 2011-2014 participants | 19,931 (1,095) |
| 2 | Positive examination weight and valid design variables (the design sample) | 19,151 (1,024) |
| 3 | Aged 50 to 79 years | 4,609 (544) |
| 4 | No reported cardiovascular disease and not pregnant | 3,853 (355) |
| 5 | All PREVENT inputs available | 3,380 (290) |
| 6 | PREVENT inputs within the study ranges: total cholesterol 130-320 mg/dL, HDL cholesterol 20-100 mg/dL, systolic pressure 90-200 mm Hg, eGFR 15-140 mL/min/1.73 m2, BMI 18.5-39.9 kg/m2 | 2,917 (216) |
| 7 | Eligible mortality linkage and measured height | 2,911 (216) |
| 8 | At least 1 maximal-effort grip trial in each hand: **primary cohort** (`core`) | 2,650 (192) |
| 9 | Aged 60 to 79 years with classifiable mobility status: **older cohort** (`mobility_cohort`) | 1,512 (148) |

## Exposure definitions used in the models

- **Grip strength per SD lower:** `lower_z = (sex-specific mean - grip) / sex-specific SD`, with the
  survey-weighted mean and SD of the primary cohort (`results/analysis_outputs/strength_scaling.csv`).
- **Per 5 kg lower:** `lower5 = -grip / 5`.
- **Low grip strength** (guideline-category analyses): grip strength at or below the sex-specific
  survey-weighted 25th percentile of the primary cohort, 39.9 kg in men and 24.7 kg in women
  (`results/analysis_outputs/risk_strength_contrasts.csv`).
- **Functional limitation:** mobility, ADL, or IADL limitation at ages 60 to 79 years; absent only
  when all 3 are classifiable and absent.
- **Guideline PREVENT categories:** 10-year PREVENT-ASCVD risk <3%, 3% to <5%, 5% to <10%, and
  >=10%. Statin threshold: PREVENT-ASCVD >=5%. Blood-pressure threshold: PREVENT total
  cardiovascular risk >=7.5%.

## Locally generated analytic datasets

Stage 2 writes the design sample (19,151 participants) to
`output/analytic_dataset/study4_analytic_design_sample.csv`. Stage 5d writes
`output/csv/analysis_outputs/guideline_analysis_variables.csv`. Join the files by `SEQN`.
The two CSV dictionaries in `data/dictionary/` describe the fields, units and derivations.

**Survey variance.** Standard errors must be estimated on the full design sample with the cohort
as a domain, for example `subset(svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~w,
nest = TRUE, data = d), core)`. Fitting survey models to a cohort extract alone gives incorrect
standard errors.

## Checking the raw files

```
Rscript data/download_raw_data.R              # check every file against the manifest
Rscript data/download_raw_data.R --download   # download any missing file from NCHS
```

NCHS occasionally reissues public files. A downloaded file whose checksum differs from the
manifest is reported, because the results were produced with the listed versions.
