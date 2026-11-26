# Reported results

These files are the results reported in the manuscript and its supplement, as produced by the code
in this repository.

| Path | Content |
|---|---|
| `analysis_outputs/` | 67 reference outputs from stages 1 to 5d. The participant-level guideline file is verified separately by SHA-256. |
| `supplementary_values/` | Formatted values behind the supplementary figures and the Supplemental Methods. |
| `Table_1.csv`, `table1_characteristics_long.csv` | Values of Table 1, formatted and in long form. |
| `figure1_flow_values.csv` | Counts shown in Figure 1. |

Running `Rscript run_all.R` writes the same files to `output/csv/`. Stage 8 compares each file in
`analysis_outputs/` with its reproduced copy and records the result in
`output/qc/comparison_with_reported_results.csv`.
