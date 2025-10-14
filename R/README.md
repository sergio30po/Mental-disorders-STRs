# R Scripts Folder

This directory contains all R scripts used in the analysis workflow for the project *Mental-disorders-STRs*.  
Each script corresponds to a specific analytical step, ensuring modularity and reproducibility of the pipeline.

---

## ⚙️ Script Overview

| Script | Description | Main Purpose |
|---------|--------------|---------------|
| **01_Environment.R** | Loads required R packages, defines global options, and sets the working environment. | Initializes the analytical environment and ensures package consistency. |
| **02_Demographic_analysis.R** | Performs descriptive and comparative analyses of demographic and clinical variables across diagnostic groups. | Generates summary statistics and group comparisons. |
| **03_Genotype_stats.R** | Conducts statistical evaluation of STR genotypes in *HTT*, *ATXN1*, and *ATXN2* genes. | Tests allelic distributions and differences between clinical groups. |
| **04_CAG_repeat_sizes.R** | Compares CAG repeat length distributions and calculates effect sizes. | Includes Kruskal-Wallis tests, Dunn post hoc comparisons, and effect size classification. |
| **05_Regression_models.R** | Fits binomial, multinomial, and linear regression models using genetic, clinical, and lifestyle predictors. | Identifies significant associations and interactions influencing disease risk and clinical features. |
| **06_Survival_age_analysis.R** | Performs survival-type and age-related analyses, including onset and disease duration models. | Evaluates the relationship between genetic variables and clinical progression. |
| **07_Enrichr-KG.R** | Integrates results from gene enrichment and knowledge graph analyses (enrich-KG). | Constructs and visualizes functional interaction networks using edge and node data. |
| **Pipeline_Mental-disorders-STRs.Rmd** | Complete pipeline | Complete R Markdown for all analysis. |

---

## 🧩 Workflow Notes

1. Scripts are designed to be executed sequentially (from `01_` to `07_`).  
2. Input data files are automatically retrieved from the `data/` directory.  
3. Output files, plots, and results are stored in the `results/` folder.  
4. Each script can be run independently, provided that the environment (`01_Environment.R`) has been initialized.  
5. The pipeline supports reproducible execution within RStudio or command-line R sessions.

---

## 📦 Dependencies

Main R packages used across scripts include:  
`tidyverse`, `ggplot2`, `dplyr`, `DescTools`, `rstatix`, `survival`, `nnet`, and `igraph`.

---

**Author:** Sergio Pérez-Oliveira  
**Last update:** October 2025
