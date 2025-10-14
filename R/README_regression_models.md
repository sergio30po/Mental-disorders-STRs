# Regression Models for STR-Psychiatric Study

This repository contains the complete R script used to perform logistic and linear regression analyses assessing the association between CAG repeat lengths in the *HTT*, *ATXN1*, and *ATXN2* genes and various psychiatric and clinical outcomes.

## Repository structure

- `05_Regression_models.R`: Full script containing all regression models evaluated in the study.

## Summary of analyses included

The script covers:

1. **Diagnostic risk models** (BD vs SCZ):
   - Combined model including SCA1 and SCA2 alleles
   - Individual models for HTT, ATXN1, and ATXN2

2. **Multinomial model**:
   - Comparison between Control, BD, and SCZ using STR and lifestyle variables

3. **Bipolar disorder subtypes**:
   - Logistic regression models for clinical subtype (e.g., CD vs non-CD) stratified by gene

4. **Cognitive deterioration in BD**:
   - Logistic regression model using HTT genotype and covariates

5. **Age at onset models**:
   - Linear regression models for BD and SCZ evaluating STR effects on age at disease onset

## Notes

- Only statistically significant or relevant findings were reported in the main manuscript.
- All other models (including those with non-significant results) are included here for transparency and reproducibility.
- Covariates such as sex, smoking status, and coffee consumption were included when available.

## How to use

To run the analysis:
1. Load your R environment and required packages.
2. Source the environment file when prompted in the script.
3. Execute the script sequentially or run individual model blocks as needed.

---

For questions or feedback, please contact the repository author.

