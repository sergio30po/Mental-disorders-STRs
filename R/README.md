# R Scripts Folder

This directory contains the R scripts used in the analysis workflow for the *Mental-disorders-STRs* project.

The reviewer-revision workflow is organized sequentially from `01_Environment.R` through `09_Reviewer_master_results.R`.

---

## ⚙️ Script Overview

| Script | Description | Main Purpose |
|---|---|---|
| **01_Environment.R** | Loads packages, imports the mental-disorder and control datasets, preprocesses variables, and creates shared analysis objects. | Initializes the analytical environment and derived datasets used by downstream scripts. |
| **02_Demographic_analysis.R** | Performs descriptive and comparative analyses of demographic and clinical variables across diagnostic groups. | Generates cohort-level descriptive statistics and group comparisons. |
| **03_Genotype_stats.R** | Tests intermediate-allele frequencies in *HTT*, *ATXN1*, and *ATXN2* with explicit case/reference coding. | Produces oriented ORs, predefined testing families, adjusted p-values, and reviewer-facing master tables. |
| **04_CAG_repeat_sizes.R** | Derives within-participant short/long alleles and evaluates continuous CAG repeat distributions. | Performs omnibus/pairwise sensitivity analyses with explicit multiplicity control and generates separate short/long plus 2D short×long figures. |
| **05_Regression_models.R** | Fits the direct *HTT* IA BD-I vs control logistic models and prespecified continuous-CAG multinomial gene-block models. | Evaluates covariate robustness of the *HTT* BD-I signal and global continuous-CAG contributions without stepwise/AIC gene selection. |
| **06_Survival_age_analysis.R** | Performs age-at-onset and observed disease-duration analyses. | Replaces the previous survival interpretation with prespecified onset analyses and exploratory cross-sectional duration analyses. |
| **07_Enrichr-KG.R** | Validates and summarizes the Enrichr-KG network and target-centered first-degree context. | Retains the network only as contextual/hypothesis-generating gene-level information, not repeat-specific functional evidence. |
| **08_HTT_meta_analysis.R** | Harmonizes published *HTT* intermediate-allele counts with the present cohort and performs fixed/random-effects synthesis. | Reports study effects, pooled ORs, heterogeneity metrics, eligibility decisions, and a forest plot. |
| **09_Reviewer_master_results.R** | Reads final outputs from scripts 03–08 and standardizes them into a common reviewer-facing schema. | Produces a master results table, testing-family summary, correction-stable subset, and nominal-only subset without redefining p-value families. |

---

## 🧩 Reviewer-Revision Inference Framework

- **Intermediate-allele comparisons (`03`)**: explicit case/reference OR direction and predefined testing families.
- **Continuous repeat-length comparisons (`04`)**: short/long alleles derived within participant; Holm or BH-FDR according to the prespecified family.
- **Regression sensitivity (`05`)**: direct *HTT* BD-I model adjusted for sex, age, and APOE ε4; continuous gene inference uses global multinomial gene-block likelihood-ratio tests.
- **Clinical outcomes (`06`)**: age-at-onset and observed disease duration are handled as separate analysis families. Disease duration is not treated as a survival endpoint.
- **Network context (`07`)**: contextual only; no repeat-QTL or CAG-specific mechanistic inference.
- **External synthesis (`08`)**: current primary synthesis is all-BD vs control because the published Ferrari et al. data do not provide BD-I-specific IA carrier counts.
- **Master integration (`09`)**: carries forward the already computed effects and adjusted p-values from scripts 03–08; no new multiplicity correction is introduced.

Machine-readable reviewer-revision outputs are written to:

```text
results/reviewer_revision/
```

Generated reviewer-requested figures are written to:

```text
figures/
```

---

## ▶️ Workflow Notes

1. Start from a fresh R session for a full reproducibility run.
2. Run scripts sequentially from `01_Environment.R` through `09_Reviewer_master_results.R`.
3. `01_Environment.R` interactively requests the primary input files and initializes shared objects.
4. Scripts `03`–`08` generate reviewer-facing results in `results/reviewer_revision/`.
5. The revised `06` workflow contains no Cox, Kaplan-Meier, log-rank, or other time-to-event inference for disease duration.
6. The revised `05` workflow contains no stepwise, greedy, or AIC-based gene-selection inference.

---

## 📦 Main Dependencies

Core packages used across the reviewer-revision workflow include:

`tidyverse`, `dplyr`, `readxl`, `writexl`, `ggplot2`, `rstatix`, `nnet`, `lmtest`, `broom`, `splines`, `igraph`, `visNetwork`, `purrr`, `readr`, and `patchwork`.

---

**Author:** Sergio Pérez-Oliveira  
**Last update:** October 2026
