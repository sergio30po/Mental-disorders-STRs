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
#   - Dataframes: DT, MENTAL, BD, SCZ

# Outputs:
#   - Logistic regression models for disorder risk
#   - Linear regression models for age-of-onset
#   - Model summaries, stepwise selection results, and effect sizes
#   - Figure 2

# ==============================================================================

# Load environment  ----
Env_path <- file.choose()
source(Env_path)
rm(Env_path)

# 1. MENTAL DISORDER RISK MODEL (Multinomial) ======
# Final: One model per gene (short + long + quadratic + interaction), adjusted for SEX + AGE + APOE_E4
# Plus: corrected block-deletion error for ATXN2.

# Data prep ----
DT <- DT %>%
  mutate(
    PATHOLOGY = relevel(factor(PATHOLOGY), ref = "CONTROL"),
    SEX = factor(SEX)
  ) %>%
  filter(
    ATXN1_CODE != "EXPANDED",
    ATXN2_CODE != "EXPANDED",
    HTT_CODE  != "EXPANDED"
  ) %>%
  mutate(
    ATXN1_CODE = droplevels(factor(ATXN1_CODE)),
    ATXN2_CODE = droplevels(factor(ATXN2_CODE)),
    HTT_CODE  = droplevels(factor(HTT_CODE))
  )

drop_na_vars <- function(df, vars) {
  dplyr::filter(df, dplyr::if_all(dplyr::all_of(vars), ~ !is.na(.)))
}

vars_allele <- c(
  "PATHOLOGY", "SEX", "AGE", "APOE_E4",
  "ALLELE1_ATXN1","ALLELE2_ATXN1","ALLELE1_ATXN2","ALLELE2_ATXN2","ALLELE1_HTT","ALLELE2_HTT"
)
DT_allele <- drop_na_vars(DT, vars_allele)

# Full poly model (used only for block-deletion AIC) ----
model_full_poly <- nnet::multinom(
  PATHOLOGY ~
    SEX + AGE + APOE_E4 +
    # ATXN1
    ALLELE1_ATXN1 + I(ALLELE1_ATXN1^2) +
    ALLELE2_ATXN1 + I(ALLELE2_ATXN1^2) +
    ALLELE1_ATXN1:ALLELE2_ATXN1 +
    # ATXN2
    ALLELE1_ATXN2 + I(ALLELE1_ATXN2^2) +
    ALLELE2_ATXN2 + I(ALLELE2_ATXN2^2) +
    ALLELE1_ATXN2:ALLELE2_ATXN2 +
    # HTT
    ALLELE1_HTT + I(ALLELE1_HTT^2) +
    ALLELE2_HTT + I(ALLELE2_HTT^2) +
    ALLELE1_HTT:ALLELE2_HTT,
  data = DT_allele,
  trace = FALSE
)

# Block-deletion models ----
model_no_HTT <- update(
  model_full_poly,
  . ~ . -
    ALLELE1_HTT - I(ALLELE1_HTT^2) -
    ALLELE2_HTT - I(ALLELE2_HTT^2) -
    ALLELE1_HTT:ALLELE2_HTT
)

model_no_ATXN1 <- update(
  model_full_poly,
  . ~ . -
    ALLELE1_ATXN1 - I(ALLELE1_ATXN1^2) -
    ALLELE2_ATXN1 - I(ALLELE2_ATXN1^2) -
    ALLELE1_ATXN1:ALLELE2_ATXN1
)

model_no_ATXN2 <- update(
  model_full_poly,
  . ~ . -
    ALLELE1_ATXN2 - I(ALLELE1_ATXN2^2) -
    ALLELE2_ATXN2 - I(ALLELE2_ATXN2^2) -
    ALLELE1_ATXN2:ALLELE2_ATXN2
)

aic_tbl <- tibble::tibble(
  model = c("Full (HTT+ATXN1+ATXN2)", "No HTT", "No ATXN1", "No ATXN2"),
  AIC = c(AIC(model_full_poly), AIC(model_no_HTT), AIC(model_no_ATXN1), AIC(model_no_ATXN2))
) %>%
  dplyr::mutate(delta_AIC = AIC - min(AIC)) %>%
  dplyr::arrange(AIC)

print(aic_tbl)

# Final models: ONE per gene ======
# Each model: covariates + (short linear+quadratic) + (long linear+quadratic) + interaction

fit_gene_multinom <- function(df, outcome = "PATHOLOGY",
                              covars = c("SEX","AGE","APOE_E4"),
                              short, long) {
  f <- as.formula(paste(
    outcome, "~",
    paste(covars, collapse = " + "), "+",
    paste0(short, " + I(", short, "^2) + ",
           long,  " + I(", long,  "^2) + ",
           short, ":", long)
  ))
  nnet::multinom(f, data = df, trace = FALSE)
}

tidy_multinom_or <- function(m) {
  broom::tidy(m, exponentiate = TRUE, conf.int = TRUE) %>%
    dplyr::filter(term != "(Intercept)")
}

# 1) HTT final
model_HTT_final <- fit_gene_multinom(
  df = DT_allele,
  short = "ALLELE1_HTT",
  long  = "ALLELE2_HTT"
)
summary(model_HTT_final)
tt_HTT_final <- tidy_multinom_or(model_HTT_final)
tt_HTT_final
lmtest::lrtest(model_HTT_final, model_full_poly)
lmtest::lrtest(model_no_HTT, model_full_poly)

# 2) ATXN1 final
model_ATXN1_final <- fit_gene_multinom(
  df = DT_allele,
  short = "ALLELE1_ATXN1",
  long  = "ALLELE2_ATXN1"
)
summary(model_ATXN1_final)
tt_ATXN1_final <- tidy_multinom_or(model_ATXN1_final)
tt_ATXN1_final
lmtest::lrtest(model_no_ATXN1, model_full_poly)

# 3) ATXN2 final
model_ATXN2_final <- fit_gene_multinom(
  df = DT_allele,
  short = "ALLELE1_ATXN2",
  long  = "ALLELE2_ATXN2"
)
summary(model_ATXN2_final)
tt_ATXN2_final <- tidy_multinom_or(model_ATXN2_final)
tt_ATXN2_final
lmtest::lrtest(model_no_ATXN2, model_full_poly)

# FIGURE 2.1 ATXN2 multinomial figure ----
# --- Settings
fig_dir <- "figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Ensure model and tidy data exist
model_plot <- model_ATXN2_final 
tt_plot <- broom::tidy(model_plot, exponentiate = TRUE, conf.int = TRUE)

# --- Settings Colors & Limits
cols_outcome <- c("BD" = "#8CBDE6", "SCZ" = "#F5A04D")
x_min <- 1e-3
x_max <- 1e4

# --- Data Preparation
df_fp <- tt_plot %>%
  filter(term != "(Intercept)") %>%
  mutate(
    outcome = factor(str_trim(as.character(y.level)), levels = c("BD", "SCZ")),
    
    # MODIFICATION: Using atop() to split labels into two lines
    term_label = case_when(
      term == "SEX[T.Female]" ~ "'Female sex'",
      term == "AGE" ~ "'Age'",
      
      term == "APOE_E4[T.E4+]" ~ "italic(APOE)~epsilon[4]~carrier",
      
      # Split ATXN2 labels:
      term == "ALLELE1_ATXN2" ~ "atop(italic(ATXN2), 'short allele (linear)')",
      stringr::str_detect(term, "^I\\(ALLELE1_ATXN2\\^2\\)") ~ "atop(italic(ATXN2), 'short allele (quadratic)')",
      
      term == "ALLELE2_ATXN2" ~ "atop(italic(ATXN2), 'long allele (linear)')",
      stringr::str_detect(term, "^I\\(ALLELE2_ATXN2\\^2\\)") ~ "atop(italic(ATXN2), 'long allele (quadratic)')",
      
      term == "ALLELE1_ATXN2:ALLELE2_ATXN2" ~ "atop(italic(ATXN2), 'short x long allele (interaction)')",
      TRUE ~ term
    ),
    
    # Stats processing
    p_plot = pmax(p.value, 1e-300),
    sig = -log10(p_plot),
    est_p = pmin(pmax(estimate,  x_min), x_max),
    lo_p  = pmin(pmax(conf.low,  x_min), x_max),
    hi_p  = pmin(pmax(conf.high, x_min), x_max),
    cut_left  = conf.low  < x_min,
    cut_right = conf.high > x_max
  )

# --- Order terms (Must match the strings in case_when EXACTLY)
order_terms <- c(
  "'Female sex'", 
  "'Age'", 
  "italic(APOE)~epsilon[4]~carrier",
  "atop(italic(ATXN2), 'short allele (linear)')",
  "atop(italic(ATXN2), 'short allele (quadratic)')",
  "atop(italic(ATXN2), 'long allele (linear)')",
  "atop(italic(ATXN2), 'long allele (quadratic)')",
  "atop(italic(ATXN2), 'short x long allele (interaction)')"
)

df_fp <- df_fp %>%
  mutate(term_label = factor(term_label, levels = rev(order_terms)))

# --- Plotting
pd <- position_dodge(width = 0.55)

g_forest <- ggplot(df_fp, aes(x = est_p, y = term_label)) +
  coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
  
  geom_vline(xintercept = 1, linetype = "dashed",
             linewidth = 0.6, color = "grey45") +
  
  # Error bars
  geom_errorbarh(
    aes(xmin = lo_p, xmax = hi_p, color = outcome),
    position = pd, height = 0.18, linewidth = 0.9
  ) +
  
  # Points
  geom_point(
    aes(fill = outcome, size = sig),
    shape = 21, color = "black", stroke = 0.35,
    position = pd, alpha = 0.95
  ) +
  
  # CI truncation arrows (Left)
  geom_segment(
    data = df_fp %>% filter(cut_left),
    aes(x = x_min * 1.35, xend = x_min * 1.08,
        y = term_label, yend = term_label, color = outcome),
    inherit.aes = FALSE,
    arrow = arrow(type = "closed", length = unit(0.16, "cm")),
    linewidth = 0.9
  ) +
  
  # CI truncation arrows (Right)
  geom_segment(
    data = df_fp %>% filter(cut_right),
    aes(x = x_max / 1.35, xend = x_max / 1.08,
        y = term_label, yend = term_label, color = outcome),
    inherit.aes = FALSE,
    arrow = arrow(type = "closed", length = unit(0.16, "cm")),
    linewidth = 0.9
  ) +
  
  scale_x_log10(
    limits = c(x_min, x_max), 
    expand = expansion(mult = c(0.05, 0)), 
    name = "Odds ratio (log scale)"
  ) +
  
  # parse=TRUE renders the atop() logic
  scale_y_discrete(labels = function(x) parse(text = x)) +
  
  scale_fill_manual(
    values = cols_outcome, 
    name = "Diagnosis"
  ) +
  scale_color_manual(
    values = cols_outcome, 
    name = "Diagnosis"
  ) +
  
  scale_size_continuous(
    name = expression(atop(-log[10](italic(p)), "(with 95% CI)")),
    range = c(2.6, 6.8) 
  ) +
  
  guides(
    color = "none",
    fill  = guide_legend(order = 1, override.aes = list(size = 5)),
    size  = guide_legend(order = 2)
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(lineheight = 0.8),
    legend.position = c(0.99, 0.99),
    legend.justification = c(1, 1),
    legend.box = "horizontal",
    legend.box.just = "top",
    legend.spacing.x = unit(0.4, "cm"),
    legend.background = element_rect(fill = "white", colour = "grey70", linewidth = 0.4),
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9),
    
    plot.margin = margin(8, 20, 8, 8)
  )

print(g_forest)

ggsave(
  filename = file.path(fig_dir, "ATXN2_forest_plot.tiff"),
  plot = g_forest,
  device = "tiff",
  width = 500, height = 160, units = "mm",
  dpi = 600, compression = "lzw"
)
# FIGURE 2.2/3: Predicted probabilities vs allele size----
# --- Model & Data Setup
DT_plot <- DT_allele 
model_plot <- model_ATXN2_final

# 1. Define SEX levels
sex_levels <- levels(DT_plot$SEX)
if (length(sex_levels) < 2) stop("SEX must have at least 2 levels.")

# 2. Define AGE quantiles (P25, P50, P75)
age_q <- quantile(DT_plot$AGE, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
age_df <- tibble(
  AGE = as.numeric(age_q),
  age_group = factor(c("Age P25", "Age P50", "Age P75"),
                     levels = c("Age P25", "Age P50", "Age P75"))
)

# 3. Define allele medians (to hold the non-varying allele constant)
med_short <- median(DT_plot$ALLELE1_ATXN2, na.rm = TRUE)
med_long  <- median(DT_plot$ALLELE2_ATXN2, na.rm = TRUE)

# 4. Define APOE reference
if(is.factor(DT_plot$APOE_E4)) {
  ref_apoe <- levels(DT_plot$APOE_E4)[1] 
} else {
  ref_apoe <- "E4-"
}

# 5. Prediction Grids
grid_short <- seq(min(DT_plot$ALLELE1_ATXN2, na.rm = TRUE), 
                  max(DT_plot$ALLELE1_ATXN2, na.rm = TRUE), by = 0.05)
grid_long  <- seq(min(DT_plot$ALLELE2_ATXN2, na.rm = TRUE), 
                  max(DT_plot$ALLELE2_ATXN2, na.rm = TRUE), by = 0.05)

# --- Helper Function
predict_probs_long <- function(model, newdata) {
  pr <- as.data.frame(predict(model, newdata = newdata, type = "probs"))
  bind_cols(newdata, pr) %>%
    pivot_longer(
      cols = all_of(colnames(pr)),
      names_to = "outcome",
      values_to = "p"
    )
}

# FIGURE 2.2: VARY LONG ALLELE ----

# Create newdata
nd_long <- expand_grid(
  ALLELE2_ATXN2 = grid_long,
  SEX = factor(sex_levels, levels = sex_levels),
  age_df,
  APOE_E4 = ref_apoe
) %>%
  mutate(ALLELE1_ATXN2 = med_short)

# Predict
pred_long <- predict_probs_long(model_plot, nd_long) %>%
  filter(outcome %in% c("BD", "SCZ")) %>%
  mutate(outcome = factor(outcome, levels = c("BD", "SCZ")))

# Plot
p_long <- ggplot(pred_long, aes(x = ALLELE2_ATXN2, y = p, color = outcome, linetype = SEX)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ age_group, nrow = 1) +
  
  scale_color_manual(values = cols_outcome, name = "Diagnosis") +
  scale_linetype_manual(values = c("solid", "dashed")[seq_along(sex_levels)], name = "Sex") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  
  labs(
    title = NULL,
    subtitle = NULL,
    x = expression(italic("ATXN2") ~ "long allele (CAG repeats)"), 
    y = "Predicted probability"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top", 
    legend.box = "horizontal", 
    legend.key.width = unit(1.5, "cm"), 
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold")
  )

print(p_long)

ggsave(
  filename = file.path(fig_dir, "Predicted_probs_long_byAge.tiff"),
  plot = p_long,
  device = "tiff",
  width = 10, height = 5, units = "in", 
  dpi = 600, compression = "lzw"
)


# FIGURE 2.3: VARY SHORT ALLELE ----

# Create newdata
nd_short <- expand_grid(
  ALLELE1_ATXN2 = grid_short,
  SEX = factor(sex_levels, levels = sex_levels),
  age_df,
  APOE_E4 = ref_apoe
) %>%
  mutate(ALLELE2_ATXN2 = med_long)

# Predict
pred_short <- predict_probs_long(model_plot, nd_short) %>%
  filter(outcome %in% c("BD", "SCZ")) %>%
  mutate(outcome = factor(outcome, levels = c("BD", "SCZ")))

# Plot
p_short <- ggplot(pred_short, aes(x = ALLELE1_ATXN2, y = p, color = outcome, linetype = SEX)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ age_group, nrow = 1) +
  
  scale_color_manual(values = cols_outcome, name = "Diagnosis") +
  scale_linetype_manual(values = c("solid", "dashed")[seq_along(sex_levels)], name = "Sex") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  
  labs(
    title = NULL,
    subtitle = NULL,
    x = expression(italic("ATXN2") ~ "short allele (CAG repeats)"),
    y = "Predicted probability"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top", 
    legend.box = "horizontal",
    legend.key.width = unit(1.5, "cm"),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold")
  )

print(p_short)

ggsave(
  filename = file.path(fig_dir, "Predicted_probs_short_byAge.tiff"),
  plot = p_short,
  device = "tiff",
  width = 10, height = 5, units = "in",
  dpi = 600, compression = "lzw"
)
# FIGURE 2: Layout Configuration ----
Figure_Composite <- 
  wrap_elements(g_forest) / 
  (p_short | p_long) +
  
  plot_layout(heights = c(1.3, 1)) + 
  
  plot_annotation(tag_levels = "A") & 
  
  theme(
    plot.tag = element_text(face = "bold", size = 20),
    plot.tag.position = c(0, 1), 
    plot.tag.padding = unit(5, "pt")
  )

# Preview
print(Figure_Composite)

# --- Save to File ---
ggsave(
  filename = file.path(fig_dir, "Figure_Composite_ATXN2.tiff"),
  plot = Figure_Composite,
  device = "tiff",
  width = 450, height = 400, units = "mm", 
  dpi = 600, compression = "lzw"
)

# 2. BINOMIAL MODELS: block-AIC selection + final per-gene models ======
# 1) Fit FULL model: covariates + gene blocks (each block: short + short^2 + long + long^2 + short:long)
# 2) Block-deletion AIC: remove each gene block from FULL and compare AIC (select best set by AIC)
# 3) (Optional) Greedy step-by-blocks: iteratively remove the gene whose removal improves AIC most
# 4) Fit FINAL per-gene models: covariates + one gene block (for each gene retained by selection)

# Helpers ----
drop_na_vars <- function(df, vars) {
  dplyr::filter(df, dplyr::if_all(dplyr::all_of(vars), ~ !is.na(.)))
}

ensure_factor_ref <- function(df, var, ref) {
  df[[var]] <- relevel(factor(df[[var]]), ref = ref)
  df
}

gene_block_string <- function(short, long) {
  paste0(
    short, " + I(", short, "^2) + ",
    long,  " + I(", long,  "^2) + ",
    short, ":", long
  )
}

make_formula <- function(outcome, covars, gene_blocks = character()) {
  rhs <- paste(covars, collapse = " + ")
  if (length(gene_blocks) > 0) rhs <- paste(rhs, paste(gene_blocks, collapse = " + "), sep = " + ")
  stats::as.formula(paste0(outcome, " ~ ", rhs))
}

fit_glm_binom <- function(df, formula) {
  stats::glm(formula, data = df, family = stats::binomial())
}

# Core function ----
fit_binomial_blockAIC <- function(df,
                                  outcome,
                                  covars,
                                  genes,
                                  outcome_ref = NULL,     # e.g., "BD"
                                  do_greedy = TRUE,       # greedy elimination of gene blocks
                                  keep_within_2AIC = FALSE,
                                  verbose = TRUE) {
  stopifnot(is.data.frame(df))
  stopifnot(all(c(outcome, covars) %in% names(df)))
  stopifnot(is.list(genes), length(genes) >= 1)
  
  # outcome as factor with chosen reference (important for interpretation)
  if (!is.null(outcome_ref)) df <- ensure_factor_ref(df, outcome, outcome_ref)
  df[[outcome]] <- factor(df[[outcome]])
  
  # common covariates as factors when present (keeps consistent with your pipeline)
  for (v in intersect(c("SEX", "COFFEE", "SMOKER", "APOE_E4"), c(outcome, covars))) {
    df[[v]] <- factor(df[[v]])
  }
  
  # Build gene block strings
  gene_names <- names(genes)
  if (is.null(gene_names) || any(gene_names == "")) {
    stop("genes must be a *named* list, e.g. list(HTT=c('ALLELE1_HTT','ALLELE2_HTT'), ...)")
  }
  gene_blocks <- purrr::imap_chr(genes, ~ gene_block_string(.x[1], .x[2]))
  
  # Complete-case for all vars used by FULL model
  allele_vars <- unique(unlist(genes))
  needed <- unique(c(outcome, covars, allele_vars))
  d <- drop_na_vars(df, needed)
  
  # FULL + NULL
  f_null <- make_formula(outcome, covars, gene_blocks = character())
  f_full <- make_formula(outcome, covars, gene_blocks = gene_blocks)
  
  m_null <- fit_glm_binom(d, f_null)
  m_full <- fit_glm_binom(d, f_full)
  
  # One-step block deletion (FULL minus each gene)
  del_tbl <- purrr::imap_dfr(gene_blocks, function(block, gname) {
    kept <- gene_blocks[setdiff(names(gene_blocks), gname)]
    f <- make_formula(outcome, covars, gene_blocks = kept)
    m <- fit_glm_binom(d, f)
    tibble::tibble(
      candidate = paste0("Full - ", gname),
      removed_gene = gname,
      k_genes = length(kept),
      AIC = stats::AIC(m)
    )
  })
  
  base_tbl <- tibble::tibble(
    candidate = c("Null (covars only)", "Full (all genes)"),
    removed_gene = c(NA_character_, NA_character_),
    k_genes = c(0L, length(gene_blocks)),
    AIC = c(stats::AIC(m_null), stats::AIC(m_full))
  )
  
  aic_once <- dplyr::bind_rows(base_tbl, del_tbl) %>%
    dplyr::mutate(delta_AIC = AIC - min(AIC)) %>%
    dplyr::arrange(AIC)
  
  # Greedy elimination: iteratively drop the single gene whose removal reduces AIC the most
  greedy_path <- NULL
  best_genes <- names(gene_blocks)
  best_model <- m_full
  best_aic <- stats::AIC(m_full)
  
  if (isTRUE(do_greedy)) {
    current_genes <- names(gene_blocks)
    current_blocks <- gene_blocks
    current_model <- m_full
    current_aic <- best_aic
    
    step_i <- 0L
    greedy_path <- tibble::tibble(
      step = step_i,
      action = "start",
      genes = paste(current_genes, collapse = ","),
      k_genes = length(current_genes),
      AIC = current_aic
    )
    
    repeat {
      if (length(current_genes) == 0) break
      
      # evaluate dropping each remaining gene
      candidates <- purrr::map_dfr(current_genes, function(g) {
        kept_genes <- setdiff(current_genes, g)
        kept_blocks <- current_blocks[kept_genes]
        f <- make_formula(outcome, covars, gene_blocks = kept_blocks)
        m <- fit_glm_binom(d, f)
        tibble::tibble(drop = g, AIC = stats::AIC(m))
      }) %>% dplyr::arrange(AIC)
      
      # pick best drop
      best_drop <- candidates$drop[1]
      best_drop_aic <- candidates$AIC[1]
      
      # stop if not improving
      if (!(best_drop_aic + 1e-8 < current_aic)) break
      
      # apply drop
      current_genes <- setdiff(current_genes, best_drop)
      current_blocks <- current_blocks[current_genes]
      current_aic <- best_drop_aic
      current_model <- fit_glm_binom(d, make_formula(outcome, covars, current_blocks))
      
      step_i <- step_i + 1L
      greedy_path <- dplyr::bind_rows(
        greedy_path,
        tibble::tibble(
          step = step_i,
          action = paste0("drop ", best_drop),
          genes = if (length(current_genes) > 0) paste(current_genes, collapse = ",") else "",
          k_genes = length(current_genes),
          AIC = current_aic
        )
      )
      
      # record best
      if (current_aic < best_aic) {
        best_aic <- current_aic
        best_genes <- current_genes
        best_model <- current_model
      }
    }
  }
  
  # Option: keep simplest within 2 AIC of best (from greedy end-models if available)
  if (isTRUE(keep_within_2AIC) && !is.null(greedy_path)) {
    minA <- min(greedy_path$AIC, na.rm = TRUE)
    candidates <- greedy_path %>% dplyr::filter(AIC <= minA + 2) %>% dplyr::arrange(k_genes, AIC)
    chosen <- candidates[1, , drop = FALSE]
    best_genes <- if (chosen$genes == "") character() else strsplit(chosen$genes, ",", fixed = TRUE)[[1]]
    best_model <- fit_glm_binom(d, make_formula(outcome, covars, gene_blocks[best_genes]))
    best_aic <- stats::AIC(best_model)
  }
  
  # Final per-gene models (covars + single gene block) for genes retained
  final_gene_models <- NULL
  final_gene_tidy <- NULL
  
  if (length(best_genes) > 0) {
    final_gene_models <- purrr::imap(genes[best_genes], function(pair, gname) {
      f <- make_formula(outcome, covars, gene_blocks = gene_block_string(pair[1], pair[2]))
      fit_glm_binom(d, f)
    })
    final_gene_tidy <- purrr::imap(final_gene_models, ~ broom::tidy(.x, exponentiate = TRUE, conf.int = TRUE))
  } else {
    final_gene_models <- list()
    final_gene_tidy <- list()
  }
  
  if (isTRUE(verbose)) {
    cat("\n# Block-deletion AIC (one-step):\n")
    print(aic_once)
    if (!is.null(greedy_path)) {
      cat("\n# Greedy block-step path:\n")
      print(greedy_path)
    }
    cat("\n# Selected genes:\n")
    print(best_genes)
    cat("\n# Selected model AIC:\n")
    print(best_aic)
  }
  
  list(
    data = d,
    formulas = list(null = f_null, full = f_full),
    models = list(null = m_null, full = m_full, selected = best_model),
    aic_once = aic_once,
    greedy_path = greedy_path,
    selected_genes = best_genes,
    final_gene_models = final_gene_models,
    final_gene_tidy = final_gene_tidy
  )
}

# 2.1 (BD vs SCZ) ======
genes <- list(
  HTT   = c("ALLELE1_HTT",  "ALLELE2_HTT"),
  ATXN1 = c("ALLELE1_ATXN1", "ALLELE2_ATXN1"),
  ATXN2 = c("ALLELE1_ATXN2", "ALLELE2_ATXN2")
)

covars_dx <- c("SEX","AGE","COFFEE","SMOKER","APOE_E4")

res_dx <- fit_binomial_blockAIC(
  df = MENTAL,
  outcome = "PATHOLOGY",
  outcome_ref = "BD",
  covars = covars_dx,
  genes = genes,
  do_greedy = TRUE,
  keep_within_2AIC = TRUE,
  verbose = TRUE
)

# # Final per-gene models (only for selected genes)
res_dx$final_gene_models
res_dx$final_gene_tidy


# 2.2 (CD in BD) ======
res_dx <- fit_binomial_blockAIC(
  df = BD,
  outcome = "CD_BINARY",
  outcome_ref = "CD",
  covars = covars_dx,
  genes = genes,
  do_greedy = TRUE,
  keep_within_2AIC = TRUE,
  verbose = TRUE
)

# # Final per-gene models (only for selected genes)
res_dx$final_gene_models
res_dx$final_gene_tidy

# 2.3 (CD in SCZ) ======
res_dx <- fit_binomial_blockAIC(
  df = SCZ,
  outcome = "CD_BINARY",
  outcome_ref = "CD",
  covars = covars_dx,
  genes = genes,
  do_greedy = TRUE,
  keep_within_2AIC = TRUE,
  verbose = TRUE
)

# # Final per-gene models (only for selected genes)
res_dx$final_gene_models
res_dx$final_gene_tidy

# 2.4 (TYPE in BD) ======
res_dx <- fit_binomial_blockAIC(
  df = BD,
  outcome = "PATHOLOGY_TYPE_BINARY",
  outcome_ref = "BD-I",
  covars = covars_dx,
  genes = genes,
  do_greedy = TRUE,
  keep_within_2AIC = TRUE,
  verbose = TRUE
)

# # Final per-gene models (only for selected genes)
res_dx$final_gene_models
res_dx$final_gene_tidy

# 2.5 (TYPE in SCZ) ======
res_dx <- fit_binomial_blockAIC(
  df = SCZ,
  outcome = "PATHOLOGY_TYPE_BINARY",
  outcome_ref = "SCZ",
  covars = covars_dx,
  genes = genes,
  do_greedy = TRUE,
  keep_within_2AIC = TRUE,
  verbose = TRUE
)

# # Final per-gene models (only for selected genes)
res_dx$final_gene_models
res_dx$final_gene_tidy

# Session info ----
sessionInfo()