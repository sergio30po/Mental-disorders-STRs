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

The project aims to evaluate the role of intermediate-length **CAG repeats** in *HTT*, *ATXN1*, and *ATXN2* genes in the etiology of **schizophrenia (SCZ)** and **bipolar disorder (BD)**. We perform association analyses, regression models, and effect size calculations using R to investigate potential modifying roles of these STR loci in mental illness phenotypes.

---

### 📊 Statistical Analysis

- All statistical analyses were performed using **R version 4.4.3**.
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

│   └── [Pipeline_Mental-disorders-STRs.Rmd](R/Pipeline_Mental-disorders-STRs.Rmd)            # Complete pipeline

├── 📁 results/                    # Output figures and summary tables

│   ├── [BD.xlsx](results/BD.xlsx)               # Biolar disorder dataset

│   ├── [SCZ.xlsx](results/SCZ.xlsx)            # Shizophrenia dataset

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
