# Script name: 02_Demographic_analysis.R
# ==============================================================================
# Title: Descriptive analysis of psychiatric cohorts (SCZ, BD, Controls)
# Author: Sergio Pérez Oliveira
# Description: This script performs descriptive comparisons across psychiatric 
#              groups (SCZ, BD, Controls) and genotype frequency analyses for 
#              HTT, ATXN1 and ATXN2 loci across subgroups (e.g., binary subtype, 
#              DCO, controls). Analyses include demographic, clinical, lifestyle 
#              variables and STR genotype distributions.
# Inputs:
#   - Manually selected environment file with custom functions (.R)
#   - Manually selected data file to update DT, SCH_CONTROLES, BP_CONTROLES
#   - Dataframes: DT, MENTAL, BD, SCH, SCH_CONTROLES, BP_CONTROLES
# Outputs:
#   - Descriptive statistics (mean ± SD, proportions)
#   - Non-parametric tests (Wilcoxon, Kruskal-Wallis, Dunn)
#   - Pairwise tests (Wilcoxon, Fisher's exact with Holm correction)
#   - Effect sizes (rank biserial)
#   - Genotype frequency tables and pairwise comparisons
#   - Missing data summary
# ==============================================================================

# Load environment ----
Env_path <- file.choose()
source(Env_path)
rm(Env_path)


# Age at onset ----
mean_sd(MENTAL, "PATHOLOGY", "ONSET_AGE", "Age at onset")
wilcox.test(ONSET_AGE ~ PATHOLOGY, data = MENTAL, exact = TRUE)
rank_biserial(ONSET_AGE ~ PATHOLOGY, data = MENTAL)

# Age at death / last visit ----
mean_sd(DT, "PATHOLOGY", "AGE", "Age at death / last visit")
kruskal.test(AGE ~ PATHOLOGY, data = DT)
dunn <- dunnTest(DT$AGE ~ DT$PATHOLOGY, method = "holm")
print(dunn, dunn.test.results = TRUE)

# Rank biserial effect sizes by pairwise comparison
levels <- unique(DT$PATHOLOGY)
comparisons <- combn(levels, 2, simplify = FALSE)
for (pair in comparisons) {
  cat("\nComparison:", pair[1], "vs", pair[2], "\n")
  sub_data <- subset(DT, PATHOLOGY %in% pair)
  sub_data$PATHOLOGY <- factor(sub_data$PATHOLOGY, levels = pair)
  print(rank_biserial(AGE ~ PATHOLOGY, data = sub_data))
}

# Disease duration ----
mean_sd(MENTAL, "PATHOLOGY", "DURATION", "Disease duration")
wilcox.test(DURATION ~ PATHOLOGY, data = MENTAL, exact = TRUE)
pairwise.wilcox.test(x = MENTAL$DURATION, g = MENTAL$PATHOLOGY, p.adjust.method = "holm")
rank_biserial(DURATION ~ PATHOLOGY, data = MENTAL)


# Sex distribution ----
Table <- TABLE(DT, "PATHOLOGY", "SEX", "Sex distribution")
rowPercents(Table)
my_function(Table)
pairwise_fisher_test(Table, p.adjust.method = "holm", conf.int = TRUE, detailed = TRUE)

# Coffee consumption ----
coffee_tab <- TABLE(MENTAL, "PATHOLOGY", "COFFEE", "Coffee consumption by pathology")
rowPercents(coffee_tab)
my_function(coffee_tab)
pairwise_fisher_test(coffee_tab, p.adjust.method = "holm", conf.int = TRUE, detailed = TRUE)

# Smoking status ----
smoking_tab <- TABLE(DT, "PATHOLOGY", "SMOKER", "Smoking by pathology")
rowPercents(smoking_tab)
my_function(smoking_tab)
pairwise_fisher_test(smoking_tab, p.adjust.method = "holm", conf.int = TRUE, detailed = TRUE)

# Cognitive decline in BD ----
cd_tab <- table(BD$CD)
round(prop.table(cd_tab) * 100, 2)
cd_tab
cd_bin_tab <- table(BD$CD_BINARY)
round(prop.table(cd_bin_tab) * 100, 2)
cd_bin_tab

# Subtype of pathology ----
subtype_tab <- table(BD$PATHOLOGY_TYPE)
round(prop.table(subtype_tab) * 100, 2)
subtype_tab
subtype_bin_tab <- table(BD$PATHOLOGY_TYPE_BINARY)
round(prop.table(subtype_bin_tab) * 100, 2)
subtype_bin_tab

sch_type_tab <- table(SCH$PATHOLOGY_TYPE)
round(prop.table(sch_type_tab) * 100, 2)
sch_type_tab
sch_type_bin_tab <- table(SCH$PATHOLOGY_TYPE_BINARY)
round(prop.table(sch_type_bin_tab) * 100, 2)
sch_type_bin_tab


# Missing data ----
DT %>%
  group_by(PATHOLOGY) %>%
  summarise(across(everything(), ~ sum(is.na(.))))

MENTAL %>%
  group_by(PATHOLOGY) %>%
  summarise(across(everything(), ~ sum(is.na(.))))


# Session info ----
sessionInfo()