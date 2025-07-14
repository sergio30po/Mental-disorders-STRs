# Mental-disorders-STRs

This repository contains the analysis code for the scientific study titled:

**Exploring the role of CAG repeats in *HTT*, *ATXN1* and *ATXN2* genes in the genetic architecture of mental disorders: schizophrenia and bipolar disorder.**

### 📋 Authors

- Sergio Pérez-Oliveira<sup>1,2,3</sup>
- Olaya Fernández-Álvarez<sup>7</sup>
- Manuel Menéndez-González
- Paz García-Portilla<sup>4,5,6</sup>
- Victoria Álvarez<sup>1,3</sup>

<sup>1</sup> Instituto de Investigación Sanitaria del Principado de Asturias (ISPA)  
<sup>2</sup> Programa de Doctorado en Biomedicina y Oncología Molecular, Universidad de Oviedo  
<sup>3</sup> Servicio de Neurología, Hospital Universitario Central de Asturias  
<sup>4</sup> Departamento de Psiquiatría, Universidad de Oviedo  
<sup>5</sup> CIBERSAM  
<sup>6</sup> Área de Salud Mental, Hospital Universitario Central de Asturias  
<sup>7</sup> Unidad de Análisis Bioinformático, Universidad de Oviedo

---

### 🧠 Project Description

The project aims to evaluate the role of intermediate-length **CAG repeats** in *HTT*, *ATXN1*, and *ATXN2* genes in the etiology of **schizophrenia (SCZ)** and **bipolar disorder (BD)**. We perform association analyses, regression models, and effect size calculations using R to investigate potential modifying roles of these STR loci in mental illness phenotypes.

---

### 📊 Statistical Analysis

- All statistical analyses were performed using **R version 4.3.2**.
- Visualizations were generated using **R** and **GraphPad Prism version 10.0**.
- A **p-value < 0.05** was considered statistically significant.
- Effect sizes, regression models (binomial/multinomial), and multiple testing corrections were applied where appropriate.

---

### 📁 Repository Structure

📁 Mental-disorders-STRs/

├── 📁 data/                        # Input datasets

│   ├── [Controls.xlsx](data/Controls.xlsx)              # Genetic data from control individuals

│   ├── [Mental_disorders.xlsx](data/Mental_disorders.xlsx)      # Clinical and genetic data from patients

│   ├── [Variables.xlsx](data/Variables.xlsx)            # Covariates and metadata

│   ├── [edges.tsv](data/edges.tsv)                 # Edgelist for gene/miRNA interaction network

│   └── [nodes.tsv](data/nodes.tsv)                   # Nodelist for gene/miRNA interaction network

├── 📁 R/                           # R scripts for analysis and visualization

│   ├── [01_Environment.R](R/01_Environment.R)           # Package setup and data import

│   ├── [02_Demographic_analysis.R](R/02_Demographic_analysis.R)  # Demographic summary and statistics

│   ├── [03_Genotype_stats.R](R/03_Genotype_stats.R)        # Analysis of genotype distributions

│   ├── [04_CAG_repeat_sizes.R](R/04_CAG_repeat_sizes.R)      # Descriptive analysis of CAG repeat sizes

│   ├── [05_Regression_models.R](R/05_Regression_models.R)     # Binomial and multinomial regression models

│   ├── [06_Survival_age_analysis.R](R/06_Survival_age_analysis.R) # Survival and age-at-onset analyses

│   └── [07_Enrichr-KG.R](R/07_Enrichr-KG.R)            # Functional enrichment using Enrichr and knowledge graphs


├── 📁 results/                    # Output figures and summary tables

│   ├── [BD.xlsx](results/BD.xlsx)               # Biolar disorder dataset

│   ├── [SCH.xlsx](results/SCH.xlsx)            # Shizophrenia dataset

│   └── [DT.xlsx](results/DT.xlsx)             # Complete dataset

├── LICENSE                       # License file for the project

└── README.md                     # Project description and instructions


---

### 📌 Citation

If you use this code in your work, please cite the corresponding paper (reference will be added upon publication).

---

### 📎 Availability

The full analysis code is openly available here at [MENTAL DISORDERS STRs](https://github.com/sergio30po/Mental-disorders-STRs).

---

### 📜 License

This project is licensed under the MIT License – see the [LICENSE](./LICENSE.txt) file for details.
