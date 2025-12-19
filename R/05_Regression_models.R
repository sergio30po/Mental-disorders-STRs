#Script name: 05_Regression_models.R
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
library(broom)

# MENTAL DISORDER RISK MODEL----

# Multinomial model

DT$PATHOLOGY <- relevel(DT$PATHOLOGY, ref = "CONTROL")
DT <- DT %>%
  filter(
    SCA1_CODE != "EXPANDED",
    SCA2_CODE != "EXPANDED",
    HTT_CODE  != "EXPANDED"
  ) %>%
  mutate(
    SCA1_CODE = droplevels(SCA1_CODE),
    SCA2_CODE = droplevels(SCA2_CODE),
    HTT_CODE  = droplevels(HTT_CODE)
  )
vars <- c("PATHOLOGY", "SEX", "AGE", "ALLELE1_SCA1", 
          "ALLELE2_SCA1", "SCA1_CODE","ALLELE1_SCA2","ALLELE2_SCA2", "SCA2_CODE","ALLELE1_HTT","ALLELE2_HTT","HTT_CODE")
DT_clean <- DT %>% filter(if_all(all_of(vars), ~ !is.na(.)))
model_multi <- multinom(PATHOLOGY ~ SEX+AGE+ALLELE1_SCA1+ALLELE2_SCA1+SCA1_CODE+ALLELE1_SCA2+ALLELE2_SCA2+SCA2_CODE+ALLELE1_HTT+ALLELE2_HTT+HTT_CODE, data = DT_clean)
summary(model_multi)
step(model_multi, direction = "backward", trace = 0)
model <- multinom(PATHOLOGY ~SEX + AGE + SCA1_CODE + ALLELE2_SCA2 + 
                    SCA2_CODE + HTT_CODE, 
                         data = DT_clean)
summary(model)
tidy(model, exponentiate = TRUE, conf.int = TRUE)

#Binomial model
MENTAL$PATHOLOGY <- relevel(MENTAL$PATHOLOGY, ref = "BD")
MENTAL <- MENTAL %>%
  filter(
    SCA1_CODE != "EXPANDED",
    SCA2_CODE != "EXPANDED",
    HTT_CODE  != "EXPANDED"
  ) %>%
  mutate(
    SCA1_CODE = droplevels(SCA1_CODE),
    SCA2_CODE = droplevels(SCA2_CODE),
    HTT_CODE  = droplevels(HTT_CODE)
  )
MLM <- select(MENTAL, -ONSET_AGE, -DURATION, -PATHOLOGY_TYPE, -CD,-CD_BINARY)

# Full model 
model <- glm(PATHOLOGY ~ ., family = binomial, data = na.omit(MLM))
summary(model)

#Backward selection
step(model, direction = "backward", trace = 0)

# Final model
mod <- glm(PATHOLOGY ~SEX + AGE + SMOKER + SCA1_CODE + ALLELE2_SCA2 + 
             SCA2_CODE, 
           data = MLM, family = binomial)
summary(mod)
lmtest::lrtest(mod)
confint(mod)
exp(coef(mod))
tidy(mod, exponentiate = TRUE, conf.int = TRUE)

# HTT RISK MODEL ----

DTLM <- MENTAL %>%
  filter(!is.na(ALLELE1_HTT), !is.na(ALLELE2_HTT), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(AGE))

model_HTT_full <- glm(PATHOLOGY ~ SEX + AGE + COFFEE + SMOKER + 
                        ALLELE1_HTT + I(ALLELE1_HTT^2) + 
                        ALLELE2_HTT + I(ALLELE2_HTT^2) + 
                        ALLELE1_HTT:ALLELE2_HTT,
                      data = DTLM, family = binomial)

step(model_HTT_full, direction = "backward", trace = 0)

# Final model
mod <- glm(PATHOLOGY ~ SEX+AGE, 
           data = MENTAL, family = binomial)
summary(mod)
tidy(mod, exponentiate = TRUE, conf.int = TRUE)

# CD models for HTT----

#BD
BD$CD_BINARY <- relevel(BD$CD_BINARY, ref = "CD")
BDLM <- BD %>%
  filter(!is.na(ALLELE1_HTT), !is.na(ALLELE2_HTT), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(CD_BINARY), !is.na(AGE))

model_CD_full <- glm(CD_BINARY ~ SEX + COFFEE + SMOKER + AGE +
                  ALLELE1_HTT + I(ALLELE1_HTT^2) + 
                  ALLELE2_HTT + I(ALLELE2_HTT^2) + 
                  ALLELE1_HTT:ALLELE2_HTT,
                data = BDLM, family = binomial())
summary(model_CD_full)
step(model_CD_full, direction = "backward", trace = 0)

mod <- glm(CD_BINARY ~ AGE + I(ALLELE1_HTT^2), 
           data = BD, family = binomial)
summary(mod)
tidy(mod, exponentiate = TRUE, conf.int = TRUE)

#SCH
SCH$CD_BINARY <- relevel(SCH$CD_BINARY, ref = "CD")
SCHLM <- SCH %>%
  filter(!is.na(ALLELE1_HTT), !is.na(ALLELE2_HTT), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(CD_BINARY), !is.na(AGE))

model_CD_full <- glm(CD_BINARY ~ SEX + COFFEE + SMOKER + AGE +
                  ALLELE1_HTT + I(ALLELE1_HTT^2) + 
                  ALLELE2_HTT + I(ALLELE2_HTT^2) + 
                  ALLELE1_HTT:ALLELE2_HTT,
                data = SCHLM, family = binomial())
summary(model_CD_full)
step(model_CD_full, direction = "backward", trace = 0)

mod <- glm(CD_BINARY ~ ALLELE1_HTT + ALLELE2_HTT, 
           data = SCH, family = binomial)
summary(mod)
tidy(mod, exponentiate = TRUE, conf.int = TRUE)

#Subtypes models for HTT----

#BD
BDLM <- BD %>%
  filter(!is.na(ALLELE1_HTT), !is.na(ALLELE2_HTT), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(PATHOLOGY_TYPE_BINARY), !is.na(AGE))

model_BD_full <- glm(PATHOLOGY_TYPE_BINARY ~ SEX + COFFEE + SMOKER + AGE +
                  ALLELE1_HTT + I(ALLELE1_HTT^2) + 
                  ALLELE2_HTT + I(ALLELE2_HTT^2) + 
                  ALLELE1_HTT:ALLELE2_HTT,
                data = BDLM, family = binomial())
summary(model_BD_full)
step(model_BD_full, direction = "backward", trace = 0)

mod <- glm(PATHOLOGY_TYPE_BINARY ~ SEX + COFFEE + AGE + ALLELE1_HTT + 
             I(ALLELE1_HTT^2) + ALLELE2_HTT + ALLELE1_HTT:ALLELE2_HTT, 
           data = BD, family = binomial)
summary(mod)
tidy(mod, exponentiate = TRUE, conf.int = TRUE)

#SCH
SCHLM <- SCH %>%
  filter(!is.na(ALLELE1_HTT), !is.na(ALLELE2_HTT), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(PATHOLOGY_TYPE_BINARY), !is.na(AGE))

model_SCH_full <- glm(PATHOLOGY_TYPE_BINARY ~ SEX + COFFEE + SMOKER + AGE + 
                  ALLELE1_HTT + I(ALLELE1_HTT^2) + 
                   ALLELE2_HTT + I(ALLELE2_HTT^2) + 
                   ALLELE1_HTT:ALLELE2_HTT,
                 data = SCHLM, family = binomial())
summary(model_SCH_full)
step(model_SCH_full, direction = "backward", trace = 0)
mod <- glm(PATHOLOGY_TYPE_BINARY ~  SEX + AGE + ALLELE1_SCA1, 
           data = SCH, family = binomial)
summary(mod)
tidy(mod, exponentiate = TRUE, conf.int = TRUE)

# ATXN1 RISK MODEL ----

DTLM <- MENTAL %>%
  filter(!is.na(ALLELE1_SCA1), !is.na(ALLELE2_SCA1), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(AGE))

model_ATXN1_full <- glm(PATHOLOGY ~ SEX + COFFEE + SMOKER + AGE + 
                     ALLELE1_SCA1 + I(ALLELE1_SCA1^2) + 
                     ALLELE2_SCA1 + I(ALLELE2_SCA1^2) + 
                     ALLELE1_SCA1:ALLELE2_SCA1,
                   data = DTLM, family = binomial)
summary(model_ATXN1_full)
step(model_ATXN1_full, direction = "backward", trace = 0)

# Final model
mod_ATXN1 <- glm(PATHOLOGY ~ SEX + AGE + ALLELE1_SCA1 + ALLELE2_SCA1 + 
                   ALLELE1_SCA1:ALLELE2_SCA1, 
                 data = MENTAL, family = binomial)
summary(mod_ATXN1)
tidy(mod_ATXN1, exponentiate = TRUE, conf.int = TRUE)

# CD model for ATXN1----

#BD
BDLM_ATXN1 <- BD %>%
  filter(!is.na(ALLELE1_SCA1), !is.na(ALLELE2_SCA1), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(CD_BINARY),!is.na(AGE)) 

model_CD_full <- glm(CD_BINARY ~ SEX + COFFEE + SMOKER + AGE +
                        ALLELE1_SCA1 + I(ALLELE1_SCA1^2) + 
                        ALLELE2_SCA1 + I(ALLELE2_SCA1^2) + 
                        ALLELE1_SCA1:ALLELE2_SCA1,
                      data = BDLM_ATXN1, family = binomial())
summary(model_CD_full)
step(model_CD_full, direction = "backward", trace = 0)
mod <- glm(CD_BINARY ~ AGE + ALLELE1_SCA1 + I(ALLELE1_SCA1^2), 
           data = BD, family = binomial)
summary(mod)
tidy(mod, exponentiate = TRUE, conf.int = TRUE)

#SCH
SCHLM_ATXN1 <- SCH %>%
  filter(!is.na(ALLELE1_SCA1), !is.na(ALLELE2_SCA1), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(CD_BINARY),!is.na(AGE)) 

model_CD_full <- glm(CD_BINARY ~ SEX + COFFEE + SMOKER + AGE + 
                        ALLELE1_SCA1 + I(ALLELE1_SCA1^2) + 
                        ALLELE2_SCA1 + I(ALLELE2_SCA1^2) + 
                        ALLELE1_SCA1:ALLELE2_SCA1,
                      data = SCHLM_ATXN1, family = binomial())
summary(model_CD_full)
step(model_CD_full, direction = "backward", trace = 0)
mod <- glm(CD_BINARY ~ ALLELE2_SCA1 + I(ALLELE2_SCA1^2), 
           data = SCH, family = binomial)
summary(mod)
tidy(mod, exponentiate = TRUE, conf.int = TRUE)

#Subtypes for ATXN1----

#BD
BDLM_ATXN1 <- BD %>%
  filter(!is.na(ALLELE1_SCA1), !is.na(ALLELE2_SCA1), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(PATHOLOGY_TYPE_BINARY), !is.na(AGE))

model_BD_full <- glm(PATHOLOGY_TYPE_BINARY ~ SEX + COFFEE + SMOKER + AGE +
                        ALLELE1_SCA1 + I(ALLELE1_SCA1^2) + 
                        ALLELE2_SCA1 + I(ALLELE2_SCA1^2) + 
                        ALLELE1_SCA1:ALLELE2_SCA1,
                      data = BDLM_ATXN1, family = binomial())
summary(model_BD_full)
step(model_BD_full, direction = "backward", trace = 0)
mod <- glm(PATHOLOGY_TYPE_BINARY ~  SEX + AGE + ALLELE1_SCA1, 
           data = BD, family = binomial)
summary(mod)
tidy(mod, exponentiate = TRUE, conf.int = TRUE)

#SCH
SCHLM_ATXN1 <- SCH %>%
  filter(!is.na(ALLELE1_SCA1), !is.na(ALLELE2_SCA1), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(PATHOLOGY_TYPE_BINARY), !is.na(AGE))

model_SCH_full <- glm(PATHOLOGY_TYPE_BINARY ~ SEX + COFFEE + SMOKER + AGE +
                         ALLELE1_SCA1 + I(ALLELE1_SCA1^2) + 
                         ALLELE2_SCA1 + I(ALLELE2_SCA1^2) + 
                         ALLELE1_SCA1:ALLELE2_SCA1,
                       data = SCHLM_ATXN1, family = binomial())
summary(model_SCH_full)
step(model_SCH_full, direction = "backward", trace = 0)
mod <- glm(PATHOLOGY_TYPE_BINARY ~ AGE + ALLELE1_SCA1 + ALLELE2_SCA1 + 
             I(ALLELE2_SCA1^2) + ALLELE1_SCA1:ALLELE2_SCA1, 
           data = SCH, family = binomial)
summary(mod)
tidy(mod, exponentiate = TRUE, conf.int = TRUE)

# ATXN2 RISK MODEL ----

DTLM_ATXN2 <- MENTAL %>%
  filter(!is.na(ALLELE1_SCA2), !is.na(ALLELE2_SCA2), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(AGE))

# Full model with covariates
model_ATXN2_full <- glm(PATHOLOGY ~ SEX + COFFEE + SMOKER + AGE +
                     ALLELE1_SCA2 + I(ALLELE1_SCA2^2) + 
                     ALLELE2_SCA2 + I(ALLELE2_SCA2^2) + 
                     ALLELE1_SCA2:ALLELE2_SCA2,
                   data = DTLM_ATXN2, family = binomial)
summary(model_ATXN2_full)
step(model_ATXN2_full, direction = "backward", trace = 0)

# Final model
mod_ATXN2 <- glm(PATHOLOGY ~ SEX + AGE + ALLELE1_SCA2 + I(ALLELE1_SCA2^2) + 
                   ALLELE2_SCA2 + I(ALLELE2_SCA2^2),
                 data = MENTAL, family = binomial)
summary(mod_ATXN2)
tidy(mod_ATXN2, exponentiate = TRUE, conf.int = TRUE)

library(boot)

# Define cost function (e.g. misclassification error)
cost_function <- function(actual, predicted_probs) {
  mean(abs(actual - (predicted_probs > 0.5)))
}

# LOOCV
cv_ATXN2 <- cv.glm(data = DTLM_ATXN2, glmfit = mod_ATXN2, cost = cost_function, K = nrow(DTLM_ATXN2))

# Print LOOCV error estimate
cv_ATXN2$delta

install.packages("ggeffects")
library(ggeffects)
library(ggplot2)

# Predicted effect of short allele (ALLELE1_SCA2)
plot_short <- ggpredict(mod_ATXN2, terms = "ALLELE1_SCA2 [all]") %>%
  plot() +
  ggtitle("Effect of Short Allele (ALLELE1_SCA2)") +
  theme_minimal()

# Predicted effect of long allele (ALLELE2_SCA2)
plot_long <- ggpredict(mod_ATXN2, terms = "ALLELE2_SCA2 [all]") %>%
  plot() +
  ggtitle("Effect of Long Allele (ALLELE2_SCA2)") +
  theme_minimal()

# Display both
library(patchwork)
plot_short + plot_long

install.packages("visreg")
library(visreg)

# Center covariates to better isolate interaction
mod_centered <- glm(PATHOLOGY ~ scale(AGE) + SEX +
                      ALLELE1_SCA2 * ALLELE2_SCA2 + 
                      I(ALLELE1_SCA2^2) + I(ALLELE2_SCA2^2),
                    data = DTLM_ATXN2, family = binomial)

# 3D interaction plot
visreg2d(mod_centered, "ALLELE1_SCA2", "ALLELE2_SCA2", 
         scale = "response", plot.type = "persp",
         main = "Predicted Probability by Allele1 and Allele2 (ATXN2)",
         xlab = "Short allele", ylab = "Long allele", zlab = "Probability")

# 2D contour plot (alternative)
visreg2d(mod_centered, "ALLELE1_SCA2", "ALLELE2_SCA2", 
         scale = "response", plot.type = "image",
         main = "Contour plot of allele interaction (ATXN2)")

# CD model for ATXN2----

#BD
BDLM_ATXN2 <- BD %>%
  filter(!is.na(ALLELE1_SCA2), !is.na(ALLELE2_SCA2), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(CD_BINARY),!is.na(AGE))

model_CD_full <- glm(CD_BINARY ~ SEX + COFFEE + SMOKER + AGE +
                        ALLELE1_SCA2 + I(ALLELE1_SCA2^2) + 
                        ALLELE2_SCA2 + I(ALLELE2_SCA2^2) + 
                        ALLELE1_SCA2:ALLELE2_SCA2,
                      data = BDLM_ATXN2, family = binomial())
summary(model_CD_full)
step(model_CD_full, direction = "backward", trace = 0)
mod <- glm(CD_BINARY ~ AGE, 
           data = BD, family = binomial)
summary(mod)
tidy(mod, exponentiate = TRUE, conf.int = TRUE)

#SCH
SCHLM_ATXN2 <- SCH %>%
  filter(!is.na(ALLELE1_SCA2), !is.na(ALLELE2_SCA2), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(CD_BINARY),!is.na(AGE))

model_CD_full <- glm(CD_BINARY ~ SEX + COFFEE + SMOKER + AGE +
                        ALLELE1_SCA2 + I(ALLELE1_SCA2^2) + 
                        ALLELE2_SCA2 + I(ALLELE2_SCA2^2) + 
                        ALLELE1_SCA2:ALLELE2_SCA2,
                      data = SCHLM_ATXN2, family = binomial())
summary(model_CD_full)
step(model_CD_full, direction = "backward", trace = 0)
mod <- glm(CD_BINARY ~ ALLELE1_SCA2 + I(ALLELE1_SCA2^2), 
           data = SCH, family = binomial)
summary(mod)
tidy(mod, exponentiate = TRUE, conf.int = TRUE)

#Subtypes for ATXN2----

#BD
BDLM_ATXN2 <- BD %>%
  filter(!is.na(ALLELE1_SCA2), !is.na(ALLELE2_SCA2), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(PATHOLOGY_TYPE_BINARY),!is.na(AGE))

model_BD_full <- glm(PATHOLOGY_TYPE_BINARY ~ SEX + COFFEE + SMOKER + AGE +
                        ALLELE1_SCA2 + I(ALLELE1_SCA2^2) + 
                        ALLELE2_SCA2 + I(ALLELE2_SCA2^2) + 
                        ALLELE1_SCA2:ALLELE2_SCA2,
                      data = BDLM_ATXN2, family = binomial())
summary(model_BD_full)
step(model_BD_full, direction = "backward", trace = 0)

#Schizophrenia Subtype
SCHLM_ATXN2 <- SCH %>%
  filter(!is.na(ALLELE1_SCA2), !is.na(ALLELE2_SCA2), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(PATHOLOGY_TYPE_BINARY),!is.na(AGE))

model_SCH_full <- glm(PATHOLOGY_TYPE_BINARY ~ SEX + COFFEE + SMOKER + AGE +
                         ALLELE1_SCA2 + I(ALLELE1_SCA2^2) + 
                         ALLELE2_SCA2 + I(ALLELE2_SCA2^2) + 
                         ALLELE1_SCA2:ALLELE2_SCA2,
                       data = SCHLM_ATXN2, family = binomial())
summary(model_SCH_full)
step(model_SCH_full, direction = "backward", trace = 0)
mod <- glm(PATHOLOGY_TYPE_BINARY ~ SEX + AGE + ALLELE1_SCA2 + 
             I(ALLELE1_SCA2^2) + ALLELE2_SCA2 + I(ALLELE2_SCA2^2) + ALLELE1_SCA2:ALLELE2_SCA2, 
           data = SCH, family = binomial)
summary(mod)
tidy(mod, exponentiate = TRUE, conf.int = TRUE)


# ==============================================================================

# AGE-ONSET MODELING ----

# HTT and Age-of-Onset ----

#BD
BDLM <- BD %>%
  filter(!is.na(ALLELE1_HTT), !is.na(ALLELE2_HTT), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(ONSET_AGE))
model_onset_BD <- lm(ONSET_AGE ~ SEX + COFFEE + SMOKER + 
                       ALLELE1_HTT + I(ALLELE1_HTT^2) + 
                       ALLELE2_HTT + I(ALLELE2_HTT^2) + 
                       ALLELE1_HTT:ALLELE2_HTT, 
                     data = na.omit(BDLM))
summary(model_onset_BD)
summary(step(model_onset_BD, direction = "backward", trace = 0))
mod <- lm(ONSET_AGE ~ I(ALLELE1_HTT^2), 
           data = BD)
summary(mod)
tidy(mod, exponentiate = TRUE, conf.int = TRUE)

#SCH
SCHLM <- SCH %>%
  filter(!is.na(ALLELE1_HTT), !is.na(ALLELE2_HTT), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(ONSET_AGE))
model_onset_SCH <- lm(ONSET_AGE ~ SEX + COFFEE + SMOKER +
                        ALLELE1_HTT + I(ALLELE1_HTT^2) + 
                        ALLELE2_HTT + I(ALLELE2_HTT^2) + 
                        ALLELE1_HTT:ALLELE2_HTT, 
                      data = SCHLM)
summary(model_onset_SCH)
step(model_onset_SCH, direction = "backward", trace = 0)
mod <- lm(ONSET_AGE ~ SEX + COFFEE + SMOKER, 
          data = SCH)
summary(mod)
tidy(mod, exponentiate = TRUE, conf.int = TRUE)

# ATXN1 and Age-of-Onset ----

#BD
BDLM <- BD %>%
  filter(!is.na(ALLELE1_SCA1), !is.na(ALLELE2_SCA1), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(ONSET_AGE))
model_onset_BD_ATXN1 <- lm(ONSET_AGE ~ SEX + COFFEE + SMOKER + 
                             ALLELE1_SCA1 + I(ALLELE1_SCA1^2) + 
                             ALLELE2_SCA1 + I(ALLELE2_SCA1^2) + 
                             ALLELE1_SCA1:ALLELE2_SCA1,
                           data = BDLM)
summary(model_onset_BD_ATXN1)
step(model_onset_BD_ATXN1, direction = "backward", trace = 0)
mod <- lm(ONSET_AGE ~ ALLELE1_SCA1 + I(ALLELE1_SCA1^2) + ALLELE2_SCA1 + 
            I(ALLELE2_SCA1^2) + ALLELE1_SCA1:ALLELE2_SCA1, 
          data = BD)
summary(mod)
tidy(mod, exponentiate = TRUE, conf.int = TRUE)

#SCH
SCHLM <- SCH %>%
  filter(!is.na(ALLELE1_SCA1), !is.na(ALLELE2_SCA1), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(ONSET_AGE))
model_onset_SCH_ATXN1 <- lm(ONSET_AGE ~ SEX + COFFEE + SMOKER + 
                              ALLELE1_SCA1 + I(ALLELE1_SCA1^2) + 
                              ALLELE2_SCA1 + I(ALLELE2_SCA1^2) + 
                              ALLELE1_SCA1:ALLELE2_SCA1,
                            data = SCHLM)
summary(model_onset_SCH_ATXN1)
step(model_onset_SCH_ATXN1, direction = "backward", trace = 0)
mod <- lm(ONSET_AGE ~ SEX + SMOKER + ALLELE1_SCA1 + ALLELE2_SCA1 + 
            ALLELE1_SCA1:ALLELE2_SCA1, 
          data = SCH)
summary(mod)
tidy(mod, exponentiate = TRUE, conf.int = TRUE)

# ATXN2 and Age-of-Onset ----

#BD
BDLM <- BD %>%
  filter(!is.na(ALLELE1_SCA2), !is.na(ALLELE2_SCA2), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(ONSET_AGE))
model_onset_BD_ATXN2 <- lm(ONSET_AGE ~ SEX + COFFEE + SMOKER + 
                             ALLELE1_SCA2 + I(ALLELE1_SCA2^2) + 
                             ALLELE2_SCA2 + I(ALLELE2_SCA2^2) + 
                             ALLELE1_SCA2:ALLELE2_SCA2,
                           data = BDLM)
summary(model_onset_BD_ATXN2)
step(model_onset_BD_ATXN2, direction = "backward", trace = 0)
mod <- lm(ONSET_AGE ~ ALLELE1_SCA2, 
          data = BD)
summary(mod)
tidy(mod, exponentiate = TRUE, conf.int = TRUE)

#SCH
SCHLM <- SCH %>%
  filter(!is.na(ALLELE1_SCA2), !is.na(ALLELE2_SCA2), !is.na(SEX), 
         !is.na(COFFEE), !is.na(SMOKER), !is.na(ONSET_AGE))
model_onset_SCH_ATXN2 <- lm(ONSET_AGE ~ SEX + COFFEE + SMOKER + 
                              ALLELE1_SCA2 + I(ALLELE1_SCA2^2) + 
                              ALLELE2_SCA2 + I(ALLELE2_SCA2^2) + 
                              ALLELE1_SCA2:ALLELE2_SCA2,
                            data = SCHLM)
summary(model_onset_SCH_ATXN2)
step(model_onset_SCH_ATXN2, direction = "backward", trace = 0)
mod <- lm(ONSET_AGE ~ SEX + COFFEE + SMOKER, 
          data = SCH)
summary(mod)
tidy(mod, exponentiate = TRUE, conf.int = TRUE)

#Session info ----
sessionInfo()