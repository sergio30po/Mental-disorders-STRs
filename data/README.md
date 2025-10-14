# Data Folder

This directory contains the anonymized clinical and genetic datasets used in the analyses of the project *Mental-disorders-STRs*.

---

## 📊 Contents

| File | Description | Format | Notes |
|------|--------------|--------|-------|
| **Controls.xlsx** | Clinical and genetic data from healthy control individuals. | Excel (.xlsx) | Used as reference group in comparative and regression analyses. |
| **Mental_disorders.xlsx** | Clinical and genetic data from patients diagnosed with schizophrenia and bipolar disorder. | Excel (.xlsx) | Includes anonymized identifiers, diagnostic subtype, and STR genotypes for *HTT*, *ATXN1*, and *ATXN2*. |
| **Variables.xlsx** | Variable definitions, data type, and coding system used across the study. | Excel (.xlsx) | Describes demographic, clinical, and genetic variables. |
| **edges.tsv** | Edge list representing relationships between genes identified in the enrichment analysis (enrich-KG). | TSV (tab-separated) | Used to reconstruct functional interaction networks. |
| **nodes.tsv** | Node table associated with `edges.tsv`. | TSV (tab-separated) | Contains gene identifiers, categories, and enrichment annotations. |

---

## ⚙️ Usage

- These files are directly imported by the R scripts in the `R/` folder.  
- All data are **fully anonymized**, but include essential genetic and diagnostic information for statistical modeling.  
- The datasets allow replication of descriptive, comparative, and enrichment analyses performed in the manuscript.  
- Do **not edit** the data manually; all preprocessing is automated within the pipeline.

---

## 🔒 Data Policy

All data are anonymized and comply with ethical and confidentiality standards for human genetic research.  
The repository does **not contain identifiable information**.  
Clinical and genetic variables are included solely for research reproducibility.

---

## 🧩 Related Scripts

- `01_data_loading_and_filtering.R`: loads and formats datasets.  
- `02_descriptive_analysis.R`: performs descriptive and comparative statistics.  
- `03_DE_analysis.R`: includes analyses of differential expansion or distribution patterns.  
- `05_Regression_models.R`: fits regression models using genetic and clinical predictors.  
- `enrich_KG_analysis.R`: generates and visualizes the enrichment network using `edges.tsv` and `nodes.tsv`.

---

**Author:** Sergio Pérez-Oliveira  
**Last update:** October 2025
