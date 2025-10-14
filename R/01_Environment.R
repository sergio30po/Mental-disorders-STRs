# Script name: 01_Environment.R
# ==============================================================================
# Title: Project environment.

# Author: Sergio Pérez Oliveira

# Purpose: Load and preprocess datasets for analysis of
#          frequencies and correlations in mental disorder cohorts.

# Dependencies: tidyverse, readxl, dplyr, ggplot2, etc.

# Description:
# This script loads mental patients and control datasets,
# converts relevant variables to factors or numeric types,
# subsets data by pathology groups (Bipolar Disorder, Schizophrenia),
# and combines datasets for further analyses.
# ==============================================================================

# Required packages ----
required_packages <- c(
  "pwr", "rcompanion","writexl", "FSA", "installr", "xlsx", "dplyr", "readxl", "ggplot2",
  "corrplot", "devtools", "ggpubr", "lsr", "Rcmdr", "survival", "KMsurv",
  "survMisc", "survminer", "ggfortify", "flexsurv", "actuar", "nortest",
  "rstatix", "gtsummary", "car", "DataExplorer", "effectsize", "lmtest", 
  "nnet","igraph","tidyverse","ggraph","tidygraph","visNetwork"
)

install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

install_and_load(required_packages)

# Input section ----
# Select working directory interactively (select the main branch as directory)
if (interactive()) {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    wd <- rstudioapi::selectDirectory(caption = "Select your working directory")
  } else if (.Platform$OS.type == "windows") {
    wd <- choose.dir(caption = "Select your working directory")
  } else {
    cat("Please set your working directory manually.\n")
    wd <- NULL
  }
  
  if (!is.null(wd)) {
    setwd(wd)
    cat("Working directory set to:", getwd(), "\n")
  } else {
    stop("No working directory selected. Please restart and select a directory.")
  }
} else {
  stop("Non-interactive session. Please set working directory manually.")
}

cat("Please select Mental disorders dataset (Excel):\n")
MENTAL <- read_excel(file.choose()) |> 
  select(-1) |> 
  select(-ADMISSION_AGE)

cat("Please select Control dataset (Excel):\n")
CONTROLS <- read_excel(file.choose())[-1, -1] |> 
  rename(AGE = DEATH_AGE)

# Functions ----

# Function to calculate mean and SD by category
mean_sd <- function(cohort, category, numeric_var, result_name) {
  categories <- unique(cohort[[category]])
  results <- data.frame()
  for (c in categories) {
    values <- na.omit(cohort[[numeric_var]][cohort[[category]] == c])
    results <- rbind(results, data.frame(
      CATEGORY = c,
      MEAN = round(mean(values), 2),
      SD = round(sd(values), 2)
    ))
  }
  colnames(results) <- c(category, paste(result_name, "MEAN"), paste(result_name, "SD"))
  return(results)
}

# Function to perform chi-square or Fisher's exact test with pairwise comparisons
my_function <- function(matrix) {
  chi_test <- chisq.test(matrix, correct = FALSE)
  perc_expected <- sum(chi_test$expected < 5) / length(chi_test$expected)
  
  if (perc_expected > 0.2) {
    fisher_test <- fisher.test(matrix)
    print("Fisher's test results:")
    print(fisher_test)
    test_result <- fisher_test$p.value
  } else {
    print("Chi-square test results:")
    print(chi_test)
    test_result <- chi_test$p.value
  }
  
  cat("\nMultiple comparisons:\n")
  print(pairwiseNominalIndependence(matrix, fisher = TRUE, chisq = TRUE, method = "holm"))
}

# Frequency table printing function
TABLE <- function(data, var1, var2, tableName) {
  cat(paste0("\nFrequency table: ", tableName, "\n"))
  table_mat <- xtabs(as.formula(paste("~", var1, "+", var2)), data = data)
  table_mat <- as.matrix(table_mat)
  table_mat[is.nan(table_mat)] <- 0
  print(table_mat)
}

# Data preprocessing ----

# Convert categorical variables to factors with labels
convert_factors <- function(df, list_of_vars_and_labels) {
  for (varname in names(list_of_vars_and_labels)) {
    df[[varname]] <- factor(df[[varname]], labels = list_of_vars_and_labels[[varname]])
  }
  return(df)
}

MENTAL <- convert_factors(MENTAL, list(
  PATHOLOGY = c("BD", "SCH"),
  SEX = c("Male", "Female"),
  HTT_CODE = c("NORMAL", "IA", "EXPANDED"),
  SCA1_CODE = c("NORMAL", "IA"),
  SCA2_CODE = c("NORMAL", "IA", "EXPANDED"),
  COFFEE = c("Non-coffee", "Coffee"),
  SMOKER = c("Smoking", "Non-smoking")
))

CONTROLS <- convert_factors(CONTROLS, list(
  PATHOLOGY = "CONTROL",
  SEX = c("Male", "Female"),
  HTT_CODE = c("NORMAL", "IA", "EXPANDED"),
  SCA1_CODE = c("NORMAL", "IA"),
  SCA2_CODE = c("NORMAL", "IA"),
  SMOKER = c("Smoking", "Non-smoking")
))

# Convert numeric variables
to_numeric <- function(df, vars) {
  df[vars] <- lapply(df[vars], as.numeric)
  return(df)
}

MENTAL <- to_numeric(MENTAL, c("ALLELE1_HTT", "ALLELE2_HTT", "ALLELE1_SCA1", "ALLELE2_SCA1",
                               "ALLELE1_SCA2", "ALLELE2_SCA2", "AGE", "ONSET_AGE", "DURATION"))

CONTROLS <- to_numeric(CONTROLS, c("ALLELE1_HTT", "ALLELE2_HTT", "ALLELE1_SCA1", "ALLELE2_SCA1",
                                   "ALLELE1_SCA2", "ALLELE2_SCA2", "AGE"))

# Subsetting by pathology groups ----

# Bipolar Disorder subset and transformations
BD <- subset(MENTAL, PATHOLOGY == "BD")
BD <- BD |> mutate(
  CD = factor(CD, labels = c("NO", "MILD", "MODERATE", "SEVERE", "VERY-SEVERE")),
  PATHOLOGY_TYPE = factor(PATHOLOGY_TYPE, labels = c("BD-I", "BD-II", "Cyclothymia", "Substance-related", "Other-BD")),
  CD_BINARY = factor(ifelse(CD == "NO", "NO", "CD"), levels = c("NO", "CD")),
  PATHOLOGY_TYPE_BINARY = factor(ifelse(PATHOLOGY_TYPE == "BD-I", "BD-I", "Other"), levels = c("BD-I", "Other"))
)

# Schizophrenia subset and transformations
SCH <- subset(MENTAL, PATHOLOGY == "SCH")
SCH <- SCH |> mutate(
  PATHOLOGY_TYPE_BINARY = factor(ifelse(PATHOLOGY_TYPE == 1, "SCH", "Other"), levels = c("SCH", "Other")),
  CD = factor(CD, labels = c("NO", "MILD", "MODERATE", "SEVERE")),
  CD_BINARY = factor(ifelse(CD == "NO", "NO", "CD"), levels = c("NO", "CD"))
)

# Add binary pathology columns to controls
CONTROLS$PATHOLOGY_TYPE_BINARY <- "CONTROL"
CONTROLS$CD_BINARY <- "CONTROL"

# Combine datasets ----

DT <- bind_rows(
  CONTROLS |> select(-CD_BINARY, -PATHOLOGY_TYPE_BINARY),
  MENTAL |> select(-COFFEE, -CD, -PATHOLOGY_TYPE, -ONSET_AGE, -DURATION)
)

DT$PATHOLOGY <- as.factor(DT$PATHOLOGY)

BD_CONTROLS <- bind_rows(
  BD |> select(-COFFEE, -CD, -ONSET_AGE, -DURATION, -PATHOLOGY_TYPE),
  CONTROLS
)

SCH_CONTROLS <- bind_rows(
  SCH |> select(-COFFEE,-CD, -ONSET_AGE, -DURATION, -PATHOLOGY_TYPE),
  CONTROLS
)

# Output section ----
# Save processed datasets for downstream analyses
if (!dir.exists("results")) {
  dir.create("results")
}

saveRDS(BD, file = "results/BD.rds")
saveRDS(SCH, file = "results/SCH.rds")
saveRDS(DT, file = "results/DT.rds")

# Save as XLSX
write_xlsx(BD, path = "results/BD.xlsx")
write_xlsx(SCH, path = "results/SCH.xlsx")
write_xlsx(DT, path = "results/DT.xlsx")

# Session info ----
sessionInfo()