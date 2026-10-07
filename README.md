# Mental-disorders-STRs

This repository contains the analysis code for the scientific study titled:

**Exploring the role of CAG repeats in *HTT*, *ATXN1* and *ATXN2* genes in the genetic architecture of mental disorders: schizophrenia and bipolar disorder.**

### 📋 Authors

- Sergio Pérez-Oliveira<sup>1,2,3</sup>
- Olaya Fernández-Álvarez<sup>4</sup>
- Manuel Menéndez-González<sup>1,5,6</sup>
- Pilar Sierra<sup>7,8,9</sup>
- Belén Arranz<sup>10</sup>
- Gemma Safont<sup>11</sup>
- Pablo García-González<sup>12</sup>
- Maitee Rosende-Roca<sup>12,13</sup>
- Mercè Boada<sup>12,13</sup>
- Agustín Ruiz<sup>12,13</sup>
- Paz García-Portilla<sup>1,14,15,16</sup>
- Victoria Álvarez<sup>1,3,15</sup>

1 Health Research Institute of the Principality of Asturias (ISPA), Oviedo, Spain

2 University of Oviedo, Oviedo, Spain

3 Genetics Laboratory, Central University Hospital of Asturias (HUCA), Oviedo, Spain

4 Asociación Parkinson Asturias (APARKAS), Oviedo, Spain

5 Department of Neurology, Central University Hospital of Asturias (HUCA), Oviedo, Spain

6 Department of Medicine, University of Oviedo, Oviedo, Spain

7 Department of Psychiatry and Psychology, La Fe University and Polytechnic Hospital, Valencia, Spain

8 Department of Medicine, University of Valencia, Valencia, Spain

9 Mental Health Research Group, La Fe Health Research Institute, Valencia, Spain

10 Parc Sanitari Sant Joan de Déu; CIBERSAM, Barcelona, Spain

11 Psychiatry Department, Hospital Universitari Mútua Terrassa, Barcelona; ISIC Medical Center, Barcelona; Universitat de Barcelona; CIBERSAM

12 Ace Alzheimer Center Barcelona, Universitat Internacional de Catalunya, 08028 Barcelona, Spain

13 Networking Research Center on Neurodegenerative Diseases (CIBERNED), Instituto de Salud Carlos III, 28029 Madrid, Spain

14 Department of Psychiatry, University of Oviedo, Oviedo, Spain

15 Health Service of the Principality of Asturias (SESPA), Oviedo, Spain

16 Biomedical Research Networking Center in Mental Health (CIBERSAM), Oviedo, Spain

---

### 🧠 Project Description

The project evaluates associations between **CAG repeat variation** in *HTT*, *ATXN1*, and *ATXN2* and **schizophrenia (SCZ)** and **bipolar disorder (BD)**. The revised analytical workflow includes categorical intermediate-allele comparisons, continuous repeat-length sensitivity analyses, covariate-adjusted regression models, age-at-onset and cross-sectional disease-duration analyses, contextual gene-network visualization, a harmonized external meta-analysis of *HTT* intermediate alleles, and a reviewer-facing master integration of the final inferential results.

The analyses are observational. Gene-network results are used as contextual, hypothesis-generating information and are not interpreted as repeat-specific functional validation.

---

### 📊 Statistical Analysis

- Analyses were executed in **R 4.6.0**.
- Odds ratios are reported with an explicit case/reference direction.
- Multiple-testing families are defined in code before adjustment.
- Holm adjustment is used for prespecified family-wise comparisons where applicable.
- BH-FDR is used for explicitly exploratory analysis families.
- Continuous CAG analyses use within-participant short/long allele ordering.
- The main direct *HTT* BD-I sensitivity model is adjusted for sex, age, and APOE ε4 status.
- Disease duration is treated as a cross-sectional observed-duration variable; no Cox, Kaplan-Meier, or log-rank analysis is used in the revised workflow.
- The *HTT* meta-analysis harmonizes intermediate-allele definitions and reports study-specific effects, pooled estimates, and heterogeneity metrics.

---

### 📁 Repository Structure

```text
Mental-disorders-STRs/
├── data/                           # Input datasets and network files
├── R/                              # Analysis scripts
│   ├── 01_Environment.R
│   ├── 02_Demographic_analysis.R
│   ├── 03_Genotype_stats.R
│   ├── 04_CAG_repeat_sizes.R
│   ├── 05_Regression_models.R
│   ├── 06_Survival_age_analysis.R
│   ├── 07_Enrichr-KG.R
│   ├── 08_HTT_meta_analysis.R
│   └── 09_Reviewer_master_results.R
├── figures/                        # Generated figures
├── results/
│   ├── BD.rds / BD.xlsx
│   ├── SCZ.rds / SCZ.xlsx
│   ├── DT.rds / DT.xlsx
│   └── reviewer_revision/          # Reviewer-facing testing plans and results
├── LICENSE
└── README.md
```

#### Analysis scripts

| Script | Main role |
|---|---|
| `01_Environment.R` | Package setup, data import, preprocessing, and creation of shared analysis objects. |
| `02_Demographic_analysis.R` | Demographic and clinical descriptive analyses. |
| `03_Genotype_stats.R` | Intermediate-allele frequency comparisons with explicit OR orientation and predefined testing families. |
| `04_CAG_repeat_sizes.R` | Continuous CAG sensitivity analyses, short/long allele QC, multiplicity control, and 1D/2D allele-size figures. |
| `05_Regression_models.R` | Direct *HTT* BD-I logistic sensitivity analysis and prespecified continuous-CAG multinomial gene-block tests. |
| `06_Survival_age_analysis.R` | Age-at-onset and cross-sectional observed disease-duration analyses; no survival/time-to-event inference. |
| `07_Enrichr-KG.R` | Contextual, hypothesis-generating gene-network analysis. |
| `08_HTT_meta_analysis.R` | Harmonized *HTT* intermediate-allele meta-analysis using currently available published aggregate data. |
| `09_Reviewer_master_results.R` | Integrates final reviewer-revision outputs into a single master table without redefining or recalculating testing families. |

The `results/reviewer_revision/` directory contains machine-readable testing plans, QC tables, study-specific estimates, adjusted results, and analysis-decision outputs generated by the reviewer-revision scripts.

---

### ▶️ Execution

The scripts are intended to be run sequentially from `01_` through `09_`.

`01_Environment.R` interactively requests the main input datasets and initializes the shared analysis objects. Subsequent scripts use those objects and generate outputs in `results/` and `figures/`.

For a clean reproducibility check, start a fresh R session, remove previously generated outputs if appropriate, and rerun the scripts in numerical order.

---

### 📌 Citation

If you use this code in your work, please cite the corresponding paper (reference will be added upon publication).

---

### 📎 Availability

The full analysis code is openly available at [MENTAL DISORDERS STRs](https://github.com/sergio30po/Mental-disorders-STRs).

---

### 📜 License

This project is licensed under the MIT License – see the [LICENSE](./LICENSE.txt) file for details.
