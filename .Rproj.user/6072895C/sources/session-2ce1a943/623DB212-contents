# Script name: 06_Survival_age_analysis.R
# ==============================================================================
# Title: Clinical and survival analysis by genotype in psychiatric cohorts.

# Author: Sergio Pérez Oliveira

# Description: This script evaluates the association between STR-based genotype 
#              classifications (HTT_CODE, SCA1_CODE, SCA2_CODE) and clinical 
#              variables (age at onset, disease duration) in schizophrenia (SCH) 
#              and bipolar disorder (BD) patients. Analyses include subgroup 
#              comparisons within BD (e.g., BD-I, BD-CD) and survival models.

# Inputs:
#   - Manually selected environment file with custom functions (.R)
#   - Manually selected data file with clinical and genotype variables
#   - Dataframes: DT, BD, SCH, BD_I, BD_CD, BD_noCD

# Outputs:
#   - Descriptive statistics for age at onset and disease duration (mean ± SD)
#   - Non-parametric tests (Wilcoxon, Kruskal-Wallis, Dunn with Holm correction)
#   - Subgroup comparisons within BD and SCH by genotype category
#   - Survival objects using Surv() from survival package
#   - Cox proportional hazards models (with covariates: SEX, SMOKER, genotype)
#   - Log-rank tests and Kaplan-Meier plots (e.g., for SCA2_CODE in BD)
# ==============================================================================

# Load environment ----
Env_path <- file.choose()
source(Env_path)
rm(Env_path)

# Subsetting by BD-I and CD
BD_I <- subset(BD, BD$PATHOLOGY_TYPE_BINARY == "BD-I")
BD_CD <- subset(BD, BD$CD_BINARY == "CD")
BD_NOCD <- subset(BD, BD$CD_BINARY == "NO")

# Subsetting by CD in SCH
SCH_CD <- subset(SCH, SCH$CD_BINARY == "CD")
SCH_NOCD <- subset(SCH, SCH$CD_BINARY == "NO")

# AGE AT ONSET CORRELATIONS ----

# BD - HTT
mean_sd(BD, "HTT_CODE", "ONSET_AGE", "Age at onset BD")
kruskal.test(ONSET_AGE ~ HTT_CODE, data = BD)
dunn <- dunnTest(BD$ONSET_AGE ~ BD$HTT_CODE, method = "holm")
print(dunn, dunn.test.results = TRUE)

mean_sd(BD_I, "HTT_CODE", "ONSET_AGE", "Age at onset BD-I")
kruskal.test(ONSET_AGE ~ HTT_CODE, data = BD_I)

mean_sd(BD_CD, "HTT_CODE", "ONSET_AGE", "Age at onset BD-CD")
kruskal.test(ONSET_AGE ~ HTT_CODE, data = BD_CD)

mean_sd(BD_NOCD, "HTT_CODE", "ONSET_AGE", "Age at onset BD-NO-CD")
kruskal.test(ONSET_AGE ~ HTT_CODE, data = BD_NOCD)

# SCH - HTT
mean_sd(SCH, "HTT_CODE", "ONSET_AGE", "Age at onset SCH")
wilcox.test(ONSET_AGE ~ HTT_CODE, data = SCH)

mean_sd(SCH_CD, "HTT_CODE", "ONSET_AGE", "Age at onset SCH_CD")
kruskal.test(ONSET_AGE ~ HTT_CODE, data = SCH_CD)

mean_sd(SCH_NOCD, "HTT_CODE", "ONSET_AGE", "Age at onset SCH_NOCD")
kruskal.test(ONSET_AGE ~ HTT_CODE, data = SCH_NOCD)

# BD - SCA1
mean_sd(BD, "SCA1_CODE", "ONSET_AGE", "Age at onset BD")
wilcox.test(ONSET_AGE ~ SCA1_CODE, data = BD)

mean_sd(BD_I, "SCA1_CODE", "ONSET_AGE", "Age at onset BD-I")
kruskal.test(ONSET_AGE ~ SCA1_CODE, data = BD_I)

mean_sd(BD_CD, "SCA1_CODE", "ONSET_AGE", "Age at onset BD-CD")
kruskal.test(ONSET_AGE ~ SCA1_CODE, data = BD_CD)

mean_sd(BD_NOCD, "SCA1_CODE", "ONSET_AGE", "Age at onset BD-NO-CD")
kruskal.test(ONSET_AGE ~ SCA1_CODE, data = BD_NOCD)

# SCH - SCA1
mean_sd(SCH, "SCA1_CODE", "ONSET_AGE", "Age at onset SCH")
wilcox.test(ONSET_AGE ~ SCA1_CODE, data = SCH)

mean_sd(SCH_CD, "SCA1_CODE", "ONSET_AGE", "Age at onset SCH_CD")
kruskal.test(ONSET_AGE ~ SCA1_CODE, data = SCH_CD)

mean_sd(SCH_NOCD, "SCA1_CODE", "ONSET_AGE", "Age at onset SCH_NOCD")
kruskal.test(ONSET_AGE ~ SCA1_CODE, data = SCH_NOCD)

# BD - SCA2
mean_sd(BD, "SCA2_CODE", "ONSET_AGE", "Age at onset BD")
wilcox.test(ONSET_AGE ~ SCA2_CODE, data = BD)

mean_sd(BD_I, "SCA2_CODE", "ONSET_AGE", "Age at onset BD-I")
wilcox.test(ONSET_AGE ~ SCA2_CODE, data = BD_I)

mean_sd(BD_CD, "SCA2_CODE", "ONSET_AGE", "Age at onset BD-CD")
wilcox.test(ONSET_AGE ~ SCA2_CODE, data = BD_CD)

mean_sd(BD_NOCD, "SCA2_CODE", "ONSET_AGE", "Age at onset BD-NO-CD")
kruskal.test(ONSET_AGE ~ SCA2_CODE, data = BD_NOCD)

# SCH - SCA2
mean_sd(SCH, "SCA2_CODE", "ONSET_AGE", "Age at onset SCH")
kruskal.test(ONSET_AGE ~ SCA2_CODE, data = SCH)

mean_sd(SCH_CD, "SCA2_CODE", "ONSET_AGE", "Age at onset SCH_CD")
kruskal.test(ONSET_AGE ~ SCA2_CODE, data = SCH_CD)

mean_sd(SCH_NOCD, "SCA2_CODE", "ONSET_AGE", "Age at onset SCH_NOCD")
kruskal.test(ONSET_AGE ~ SCA2_CODE, data = SCH_NOCD)

# DURATION CORRELATIONS ----

# BD - HTT
mean_sd(BD, "HTT_CODE", "DURATION", "Duration BD")
kruskal.test(DURATION ~ HTT_CODE, data = BD)

mean_sd(BD_I, "HTT_CODE", "DURATION", "Duration BD-I")
kruskal.test(DURATION ~ HTT_CODE, data = BD_I)

mean_sd(BD_CD, "HTT_CODE", "DURATION", "Duration BD-CD")
kruskal.test(DURATION ~ HTT_CODE, data = BD_CD)

mean_sd(BD_NOCD, "HTT_CODE", "DURATION", "Duration BD-NO-CD")
kruskal.test(DURATION ~ HTT_CODE, data = BD_NOCD)

# SCH - HTT
mean_sd(SCH, "HTT_CODE", "DURATION", "Duration SCH")
wilcox.test(DURATION ~ HTT_CODE, data = SCH)

mean_sd(SCH_CD, "HTT_CODE", "DURATION", "Duration SCH-CD")
kruskal.test(DURATION ~ HTT_CODE, data = SCH_CD)

mean_sd(SCH_NOCD, "HTT_CODE", "DURATION", "Duration SCH-NO-CD")
kruskal.test(DURATION ~ HTT_CODE, data = SCH_NOCD)

# BD - SCA1
mean_sd(BD, "SCA1_CODE", "DURATION", "Duration BD")
wilcox.test(DURATION ~ SCA1_CODE, data = BD)

mean_sd(BD_I, "SCA1_CODE", "DURATION", "Duration BD-I")
kruskal.test(DURATION ~ SCA1_CODE, data = BD_I)

mean_sd(BD_CD, "SCA1_CODE", "DURATION", "Duration BD-CD")
kruskal.test(DURATION ~ SCA1_CODE, data = BD_CD)

mean_sd(BD_NOCD, "SCA1_CODE", "DURATION", "Duration BD-NO-CD")
kruskal.test(DURATION ~ SCA1_CODE, data = BD_NOCD)

# SCH - SCA1
mean_sd(SCH, "SCA1_CODE", "DURATION", "Duration SCH")
wilcox.test(DURATION ~ SCA1_CODE, data = SCH)

mean_sd(SCH_CD, "SCA1_CODE", "DURATION", "Duration SCH-CD")
kruskal.test(DURATION ~ SCA1_CODE, data = SCH_CD)

mean_sd(SCH_NOCD, "SCA1_CODE", "DURATION", "Duration SCH-NO-CD")
kruskal.test(DURATION ~ SCA1_CODE, data = SCH_NOCD)

# BD - SCA2
mean_sd(BD, "SCA2_CODE", "DURATION", "Duration BD")
wilcox.test(DURATION ~ SCA2_CODE, data = BD)

mean_sd(BD_I, "SCA2_CODE", "DURATION", "Duration BD-I")
wilcox.test(DURATION ~ SCA2_CODE, data = BD_I)

mean_sd(BD_CD, "SCA2_CODE", "DURATION", "Duration BD-CD")
wilcox.test(DURATION ~ SCA2_CODE, data = BD_CD)

mean_sd(BD_NOCD, "SCA2_CODE", "DURATION", "Duration BD-NoCD")
wilcox.test(DURATION ~ SCA2_CODE, data = BD_NOCD)

# SCH - SCA2
mean_sd(SCH, "SCA2_CODE", "DURATION", "Duration SCH")
kruskal.test(DURATION ~ SCA2_CODE, data = SCH)

mean_sd(SCH_CD, "SCA2_CODE", "DURATION", "Duration SCH-CD")
kruskal.test(DURATION ~ SCA2_CODE, data = SCH_CD)

mean_sd(SCH_NOCD, "SCA2_CODE", "DURATION", "Duration SCH-NO-CD")
kruskal.test(DURATION ~ SCA2_CODE, data = SCH_NOCD)

# SURVIVAL CURVES ----

# BD
surv_object <- Surv(time = BD$DURATION, event = rep(1, length(BD$DURATION)))
cox_model <- coxph(surv_object ~ SEX + SMOKER + SCA2_CODE + HTT_CODE, data = BD)
summary(cox_model)

survdiff(surv_object ~ COFFEE, data = BD, rho = 0)
survdiff(surv_object ~ SEX, data = BD, rho = 0)
survdiff(surv_object ~ SMOKER, data = BD, rho = 0)
survdiff(surv_object ~ SCA2_CODE, data = BD, rho = 0)
cox_model <- coxph(surv_object ~ SCA2_CODE, data = BD)
summary(cox_model)
survdiff(surv_object ~ HTT_CODE, data = BD, rho = 0)
survdiff(surv_object ~ SCA1_CODE, data = BD, rho = 0)

surv_object <- Surv(time = BD_NOCD$DURATION, event = rep(1, length(BD_NOCD$DURATION)))
survdiff(surv_object ~ SCA2_CODE, data = BD_NOCD, rho = 0)
survdiff(surv_object ~ HTT_CODE, data = BD_NOCD, rho = 0)
survdiff(surv_object ~ SCA1_CODE, data = BD_NOCD, rho = 0)

surv_object <- Surv(time = BD_CD$DURATION, event = rep(1, length(BD_CD$DURATION)))
survdiff(surv_object ~ SCA2_CODE, data = BD_CD, rho = 0)
cox_model <- coxph(surv_object ~ SCA2_CODE, data = BD_CD)
summary(cox_model)
survdiff(surv_object ~ HTT_CODE, data = BD_CD, rho = 0)
survdiff(surv_object ~ SCA1_CODE, data = BD_CD, rho = 0)

surv_object <- Surv(time = BD$DURATION, event = rep(1, length(BD$DURATION)))
gen.km <- survfit(surv_object ~ SCA2_CODE, data = BD, type = "kaplan-meier", error = "tsiatis", conf.type = "log-log", conf.int = 0.95)
plot_km <- ggsurvplot(
  fit = gen.km,
  data = BD,
  conf.int = TRUE,
  pval = TRUE,
  risk.table = TRUE,
  title = "ATXN2 Kaplan-Meier Curve in BD",
  xlab = "Time (years)",
  ylab = "Survival probability",
  legend.title = "ATXN2 Genotype",
  legend.labs = c("Normal", "IAs")
)
print(plot_km)

# SCH
surv_object <- Surv(time = SCH$DURATION, event = rep(1, length(SCH$DURATION)))
cox_model <- coxph(surv_object ~ SEX + SMOKER + SCA2_CODE + HTT_CODE, data = SCH)
summary(cox_model)

survdiff(surv_object ~ SEX, data = SCH, rho = 0)
survdiff(surv_object ~ SMOKER, data = SCH, rho = 0)
survdiff(surv_object ~ SCA2_CODE, data = SCH, rho = 0)
survdiff(surv_object ~ HTT_CODE, data = SCH, rho = 0)
survdiff(surv_object ~ SCA1_CODE, data = SCH, rho = 0)

surv_object <- Surv(time = SCH_NOCD$DURATION, event = rep(1, length(SCH_NOCD$DURATION)))
survdiff(surv_object ~ SCA2_CODE, data = SCH_NOCD, rho = 0)
survdiff(surv_object ~ HTT_CODE, data = SCH_NOCD, rho = 0)
survdiff(surv_object ~ SCA1_CODE, data = SCH_NOCD, rho = 0)

surv_object <- Surv(time = SCH_CD$DURATION, event = rep(1, length(SCH_CD$DURATION)))
survdiff(surv_object ~ SCA2_CODE, data = SCH_CD, rho = 0)
survdiff(surv_object ~ HTT_CODE, data = SCH_CD, rho = 0)
survdiff(surv_object ~ SCA1_CODE, data = SCH_CD, rho = 0)


# Session info ----
sessionInfo()