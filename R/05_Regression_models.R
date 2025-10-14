# Script name: 05_Regression_models.R
# ==============================================================================
# Title: Psychiatric Disorder Risk and Age-of-Onset Regression Modeling.

# Author: Sergio Pérez Oliveira

# Description: This script performs regression analyses to investigate:
#              1) Association between STR genotypes (HTT, ATXN1, ATXN2) and:
#                 - Risk of psychiatric disorders (BD, SCZ vs controls)
#                 - Risk of clinical subtypes (e.g., bipolar subtypes)
#                 - Risk of cognitive deterioration (CD)
#              2) Association between STRs and age-of-onset in psychiatric disorders
#              Includes covariate adjustment, model selection, and effect size calculations.

# Inputs:
#   - Environment file with custom functions (manually selected)
#   - Data file with preprocessed clinical/genetic data
#   - Dataframes: DT, MENTAL, BD, SCH

# Outputs:
#   - Logistic regression models for disorder risk
#   - Linear regression models for age-of-onset
#   - Model summaries, stepwise selection results, and effect sizes

# ==============================================================================

# Load environment  ----
Env_path <- file.choose()
source(Env_path)
rm(Env_path)


# MENTAL DISORDER RISK MODEL----

# Multinomial model (3 categories)
DT$PATHOLOGY <- relevel(DT$PATHOLOGY, ref = "CONTROL")

modelo_multi <- multinom(PATHOLOGY ~ ., data = na.omit(DT))
summary(modelo_multi)
step(modelo_multi, direction = "backward", trace = 0)

modelo_multi <- multinom(PATHOLOGY ~ SEX + AGE + SMOKER + ALLELE2_SCA1 + 
                           SCA1_CODE + ALLELE2_SCA2 + SCA2_CODE, 
                         data = na.omit(DT))
summary(modelo_multi)

# Calculate p-values
z <- summary(modelo_multi)$coefficients / summary(modelo_multi)$standard.errors
pvals <- 2 * (1 - pnorm(abs(z)))
pvals
exp(coef(modelo_multi))
confint(modelo_multi)

# Binomial model (2 categories)
MENTAL$PATHOLOGY <- relevel(MENTAL$PATHOLOGY, ref = "BD")
MLM <- select(MENTAL, -ONSET_AGE, -AGE, -DURATION, -PATHOLOGY_TYPE, -CD)

# Full model 
model <- glm(PATHOLOGY ~ ., family = binomial, data = na.omit(MLM))
summary(model)

# Backward selection
step(model, direction = "backward", trace = 0)

# Final model
mod <- glm(PATHOLOGY ~ SEX + SCA1_CODE + ALLELE2_SCA2 + SCA2_CODE, 
           data = MLM, family = binomial)
summary(mod)
lrtest(mod)
confint(mod)
exp(coef(mod))

# HTT RISK MODEL ----
DTLM <- MENTAL %>%
  filter(!is.na(ALLELE1_HTT), !is.na(ALLELE2_HTT), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER))

# Full model with covariates
model <- glm(PATHOLOGY ~ SEX + COFFEE + SMOKER + 
               ALLELE1_HTT + I(ALLELE1_HTT^2) + 
               ALLELE2_HTT + I(ALLELE2_HTT^2) + 
               ALLELE1_HTT:ALLELE2_HTT,
             data = DTLM, family = binomial)
summary(model)
step(model, direction = "backward", trace = 0)

# Final model
mod <- glm(PATHOLOGY ~ SEX, 
           data = MENTAL, family = binomial)
summary(mod)
lrtest(mod)
confint(mod)
exp(coef(mod))

# Model without covariates
model_nocov <- glm(PATHOLOGY ~ ALLELE1_HTT + I(ALLELE1_HTT^2) + 
                     ALLELE2_HTT + I(ALLELE2_HTT^2) + 
                     ALLELE1_HTT:ALLELE2_HTT,
                   data = DTLM, family = binomial)
summary(model_nocov)
summary(step(model_nocov, direction = "backward", trace = 0))

# Cognitive Deterioration (CD) model----
BD$CD_BINARY <- relevel(BD$CD_BINARY, ref = "CD")
BDLM <- BD %>%
  filter(!is.na(ALLELE1_HTT), !is.na(ALLELE2_HTT), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(CD_BINARY))

model_CD <- glm(CD_BINARY ~ SEX + COFFEE + SMOKER + 
                  ALLELE1_HTT + I(ALLELE1_HTT^2) + 
                  ALLELE2_HTT + I(ALLELE2_HTT^2) + 
                  ALLELE1_HTT:ALLELE2_HTT,
                data = BDLM, family = binomial())
summary(model_CD)
summary(step(model_CD, direction = "backward", trace = 0))

SCH$CD_BINARY <- relevel(SCH$CD_BINARY, ref = "CD")
SCHLM <- SCH %>%
  filter(!is.na(ALLELE1_HTT), !is.na(ALLELE2_HTT), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(CD_BINARY))

model_CD <- glm(CD_BINARY ~ SEX + COFFEE + SMOKER + 
                  ALLELE1_HTT + I(ALLELE1_HTT^2) + 
                  ALLELE2_HTT + I(ALLELE2_HTT^2) + 
                  ALLELE1_HTT:ALLELE2_HTT,
                data = SCHLM, family = binomial())
summary(model_CD)
summary(step(model_CD, direction = "backward", trace = 0))

#Subtypes----
# Bipolar Subtype
BDLM <- BD %>%
  filter(!is.na(ALLELE1_HTT), !is.na(ALLELE2_HTT), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(PATHOLOGY_TYPE_BINARY))

model_BD <- glm(PATHOLOGY_TYPE_BINARY ~ SEX + COFFEE + SMOKER + 
                  ALLELE1_HTT + I(ALLELE1_HTT^2) + 
                  ALLELE2_HTT + I(ALLELE2_HTT^2) + 
                  ALLELE1_HTT:ALLELE2_HTT,
                data = BDLM, family = binomial())
summary(model_BD)
summary(step(model_BD, direction = "backward", trace = 0))

# Schizophrenia Subtype
SCHLM <- SCH %>%
  filter(!is.na(ALLELE1_HTT), !is.na(ALLELE2_HTT), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(PATHOLOGY_TYPE_BINARY))
model_SCH <- glm(PATHOLOGY_TYPE_BINARY ~ ALLELE1_HTT + I(ALLELE1_HTT^2) + 
                   ALLELE2_HTT + I(ALLELE2_HTT^2) + 
                   ALLELE1_HTT:ALLELE2_HTT,
                 data = SCHLM, family = binomial())
summary(model_SCH)
summary(step(model_SCH, direction = "backward", trace = 0))


# ATXN1 RISK MODEL ----
DTLM <- MENTAL %>%
  filter(!is.na(ALLELE1_SCA1), !is.na(ALLELE2_SCA1), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER))

# Full model with covariates
model_ATXN1 <- glm(PATHOLOGY ~ SEX + COFFEE + SMOKER + 
                     ALLELE1_SCA1 + I(ALLELE1_SCA1^2) + 
                     ALLELE2_SCA1 + I(ALLELE2_SCA1^2) + 
                     ALLELE1_SCA1:ALLELE2_SCA1,
                   data = DTLM, family = binomial)
summary(model_ATXN1)
summary(step(model_ATXN1, direction = "backward", trace = 0))

# Final model
mod_ATXN1 <- glm(PATHOLOGY ~ SEX + ALLELE1_SCA1 + ALLELE2_SCA1 + 
                   ALLELE1_SCA1:ALLELE2_SCA1, 
                 data = MENTAL, family = binomial)
summary(mod_ATXN1)
lrtest(mod_ATXN1)
confint(mod_ATXN1)
exp(coef(mod_ATXN1))

# Model without covariates
model_ATXN1_nocov <- glm(PATHOLOGY ~ ALLELE1_SCA1 + I(ALLELE1_SCA1^2) + 
                           ALLELE2_SCA1 + I(ALLELE2_SCA1^2) + 
                           ALLELE1_SCA1:ALLELE2_SCA1,
                         data = DTLM, family = binomial)
summary(model_ATXN1_nocov)
summary(step(model_ATXN1_nocov, direction = "backward", trace = 0))

# Cognitive Deterioration (CD) model----
BDLM_ATXN1 <- BD %>%
  filter(!is.na(ALLELE1_SCA1), !is.na(ALLELE2_SCA1), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(CD_BINARY))

model_CD_ATXN1 <- glm(CD_BINARY ~ SEX + COFFEE + SMOKER + 
                        ALLELE1_SCA1 + I(ALLELE1_SCA1^2) + 
                        ALLELE2_SCA1 + I(ALLELE2_SCA1^2) + 
                        ALLELE1_SCA1:ALLELE2_SCA1,
                      data = BDLM_ATXN1, family = binomial())
summary(model_CD_ATXN1)
summary(step(model_CD_ATXN1, direction = "backward", trace = 0))

SCHLM_ATXN1 <- SCH %>%
  filter(!is.na(ALLELE1_SCA1), !is.na(ALLELE2_SCA1), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(CD_BINARY))

model_CD_ATXN1 <- glm(CD_BINARY ~ SEX + COFFEE + SMOKER + 
                        ALLELE1_SCA1 + I(ALLELE1_SCA1^2) + 
                        ALLELE2_SCA1 + I(ALLELE2_SCA1^2) + 
                        ALLELE1_SCA1:ALLELE2_SCA1,
                      data = SCHLM_ATXN1, family = binomial())
summary(model_CD_ATXN1)
summary(step(model_CD_ATXN1, direction = "backward", trace = 0))
#Subtypes----
# Bipolar Subtype
BDLM_ATXN1 <- BD %>%
  filter(!is.na(ALLELE1_SCA1), !is.na(ALLELE2_SCA1), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(PATHOLOGY_TYPE_BINARY))

model_BD_ATXN1 <- glm(PATHOLOGY_TYPE_BINARY ~ SEX + COFFEE + SMOKER + 
                        ALLELE1_SCA1 + I(ALLELE1_SCA1^2) + 
                        ALLELE2_SCA1 + I(ALLELE2_SCA1^2) + 
                        ALLELE1_SCA1:ALLELE2_SCA1,
                      data = BDLM_ATXN1, family = binomial())
summary(model_BD_ATXN1)
summary(step(model_BD_ATXN1, direction = "backward", trace = 0))

# Schizophrenia Subtype
SCHLM_ATXN1 <- SCH %>%
  filter(!is.na(ALLELE1_SCA1), !is.na(ALLELE2_SCA1), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(PATHOLOGY_TYPE_BINARY))

model_SCH_ATXN1 <- glm(PATHOLOGY_TYPE_BINARY ~ SEX + COFFEE + SMOKER + 
                         ALLELE1_SCA1 + I(ALLELE1_SCA1^2) + 
                         ALLELE2_SCA1 + I(ALLELE2_SCA1^2) + 
                         ALLELE1_SCA1:ALLELE2_SCA1,
                       data = SCHLM_ATXN1, family = binomial())
summary(model_SCH_ATXN1)
summary(step(model_SCH_ATXN1, direction = "backward", trace = 0))


# ATXN2 RISK MODEL ----
DTLM_ATXN2 <- MENTAL %>%
  filter(!is.na(ALLELE1_SCA2), !is.na(ALLELE2_SCA2), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER))

# Full model with covariates
model_ATXN2 <- glm(PATHOLOGY ~ SEX + COFFEE + SMOKER + 
                     ALLELE1_SCA2 + I(ALLELE1_SCA2^2) + 
                     ALLELE2_SCA2 + I(ALLELE2_SCA2^2) + 
                     ALLELE1_SCA2:ALLELE2_SCA2,
                   data = DTLM_ATXN2, family = binomial)
summary(model_ATXN2)
summary(step(model_ATXN2, direction = "backward", trace = 0))

# Final model
mod_ATXN2 <- glm(PATHOLOGY ~ SEX + ALLELE1_SCA2 + I(ALLELE1_SCA2^2) + 
                   ALLELE2_SCA2 + I(ALLELE2_SCA2^2),
                 data = MENTAL, family = binomial)
summary(mod_ATXN2)
lrtest(mod_ATXN2)
confint(mod_ATXN2)
exp(coef(mod_ATXN2))
summary(step(mod_ATXN2, direction = "backward", trace = 0))

# Model without covariates
model_ATXN2_nocov <- glm(PATHOLOGY ~ ALLELE1_SCA2 + I(ALLELE1_SCA2^2) + 
                           ALLELE2_SCA2 + I(ALLELE2_SCA2^2) + 
                           ALLELE1_SCA2:ALLELE2_SCA2,
                         data = DTLM_ATXN2, family = binomial)
summary(model_ATXN2_nocov)
summary(step(model_ATXN2_nocov, direction = "backward", trace = 0))

# Cognitive Deterioration (CD) model----
BDLM_ATXN2 <- BD %>%
  filter(!is.na(ALLELE1_SCA2), !is.na(ALLELE2_SCA2), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(CD_BINARY))

model_CD_ATXN2 <- glm(CD_BINARY ~ SEX + COFFEE + SMOKER + 
                        ALLELE1_SCA2 + I(ALLELE1_SCA2^2) + 
                        ALLELE2_SCA2 + I(ALLELE2_SCA2^2) + 
                        ALLELE1_SCA2:ALLELE2_SCA2,
                      data = BDLM_ATXN2, family = binomial())
summary(model_CD_ATXN2)
summary(step(model_CD_ATXN2, direction = "backward", trace = 0))

SCHLM_ATXN2 <- SCH %>%
  filter(!is.na(ALLELE1_SCA2), !is.na(ALLELE2_SCA2), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(CD_BINARY))

model_CD_ATXN2 <- glm(CD_BINARY ~ SEX + COFFEE + SMOKER + 
                        ALLELE1_SCA2 + I(ALLELE1_SCA2^2) + 
                        ALLELE2_SCA2 + I(ALLELE2_SCA2^2) + 
                        ALLELE1_SCA2:ALLELE2_SCA2,
                      data = SCHLM_ATXN2, family = binomial())
summary(model_CD_ATXN2)
summary(step(model_CD_ATXN2, direction = "backward", trace = 0))
#Subtypes----
# Bipolar Subtype
BDLM_ATXN2 <- BD %>%
  filter(!is.na(ALLELE1_SCA2), !is.na(ALLELE2_SCA2), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(PATHOLOGY_TYPE_BINARY))

model_BD_ATXN2 <- glm(PATHOLOGY_TYPE_BINARY ~ SEX + COFFEE + SMOKER + 
                        ALLELE1_SCA2 + I(ALLELE1_SCA2^2) + 
                        ALLELE2_SCA2 + I(ALLELE2_SCA2^2) + 
                        ALLELE1_SCA2:ALLELE2_SCA2,
                      data = BDLM_ATXN2, family = binomial())
summary(model_BD_ATXN2)
summary(step(model_BD_ATXN2, direction = "backward", trace = 0))

# Schizophrenia Subtype
SCHLM_ATXN2 <- SCH %>%
  filter(!is.na(ALLELE1_SCA2), !is.na(ALLELE2_SCA2), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(PATHOLOGY_TYPE_BINARY))

model_SCH_ATXN2 <- glm(PATHOLOGY_TYPE_BINARY ~ SEX + COFFEE + SMOKER + 
                         ALLELE1_SCA2 + I(ALLELE1_SCA2^2) + 
                         ALLELE2_SCA2 + I(ALLELE2_SCA2^2) + 
                         ALLELE1_SCA2:ALLELE2_SCA2,
                       data = SCHLM_ATXN2, family = binomial())
summary(model_SCH_ATXN2)
summary(step(model_SCH_ATXN2, direction = "backward", trace = 0))



# ==============================================================================

# AGE-ONSET MODELING ----


# HTT and Age-of-Onset ----

# Bipolar disorder
model_onset_BD <- lm(ONSET_AGE ~ SEX + COFFEE + SMOKER + 
                       ALLELE1_HTT + I(ALLELE1_HTT^2) + 
                       ALLELE2_HTT + I(ALLELE2_HTT^2) + 
                       ALLELE1_HTT:ALLELE2_HTT, 
                     data = na.omit(BD))
summary(model_onset_BD)
summary(step(model_onset_BD, direction = "backward", trace = 0))

# Schizophrenia
model_onset_SCH <- lm(ONSET_AGE ~ SEX + COFFEE + SMOKER + 
                        ALLELE1_HTT + I(ALLELE1_HTT^2) + 
                        ALLELE2_HTT + I(ALLELE2_HTT^2) + 
                        ALLELE1_HTT:ALLELE2_HTT, 
                      data = na.omit(SCH))
summary(model_onset_SCH)
summary(step(model_onset_SCH, direction = "backward", trace = 0))

# ATXN1 and Age-of-Onset ----

# Bipolar Disorder
model_onset_BD_ATXN1 <- lm(ONSET_AGE ~ SEX + COFFEE + SMOKER + 
                             ALLELE1_SCA1 + I(ALLELE1_SCA1^2) + 
                             ALLELE2_SCA1 + I(ALLELE2_SCA1^2) + 
                             ALLELE1_SCA1:ALLELE2_SCA1,
                           data = na.omit(BD))
summary(model_onset_BD_ATXN1)
summary(step(model_onset_BD_ATXN1, direction = "backward", trace = 0))

# Schizophrenia
model_onset_SCH_ATXN1 <- lm(ONSET_AGE ~ SEX + COFFEE + SMOKER + 
                              ALLELE1_SCA1 + I(ALLELE1_SCA1^2) + 
                              ALLELE2_SCA1 + I(ALLELE2_SCA1^2) + 
                              ALLELE1_SCA1:ALLELE2_SCA1,
                            data = na.omit(SCH))
summary(model_onset_SCH_ATXN1)
summary(step(model_onset_SCH_ATXN1, direction = "backward", trace = 0))


# ATXN2 and Age-of-Onset ----

# Bipolar Disorder
model_onset_BD_ATXN2 <- lm(ONSET_AGE ~ SEX + COFFEE + SMOKER + 
                             ALLELE1_SCA2 + I(ALLELE1_SCA2^2) + 
                             ALLELE2_SCA2 + I(ALLELE2_SCA2^2) + 
                             ALLELE1_SCA2:ALLELE2_SCA2,
                           data = na.omit(BD))
summary(model_onset_BD_ATXN2)
summary(step(model_onset_BD_ATXN2, direction = "backward", trace = 0))

# Schizophrenia
model_onset_SCH_ATXN2 <- lm(ONSET_AGE ~ SEX + COFFEE + SMOKER + 
                              ALLELE1_SCA2 + I(ALLELE1_SCA2^2) + 
                              ALLELE2_SCA2 + I(ALLELE2_SCA2^2) + 
                              ALLELE1_SCA2:ALLELE2_SCA2,
                            data = na.omit(SCH))
summary(model_onset_SCH_ATXN2)
summary(step(model_onset_SCH_ATXN2, direction = "backward", trace = 0))

# Session info ----
sessionInfo()