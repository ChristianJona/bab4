# ============================================================================
# ANALISIS LENGKAP: STLT, DSTLT, CBS, DAN METODE LAINNYA
# ============================================================================
# Script ini menjalankan analisis komprehensif untuk semua metode:
# 1. STLT (Smooth Threshold Life Table)
# 2. DSTLT (Dynamic STLT)
# 3. CBD (Cairns-Blake-Dowd)
# 4. TLT, Gompertz, Makeham, HP2
#
# Perbandingan:
# - HMD Only vs HMD+CBS
# - Female vs Male
# - Across cohorts
# ============================================================================

library(tidyverse)
library(gridExtra)
library(knitr)

cat("\n")
cat("================================================================================\n")
cat("  ANALISIS LENGKAP: ALL METHODS\n")
cat("================================================================================\n\n")

# ============================================================================
# 0. LOAD MODULES DAN FUNCTIONS
# ============================================================================

cat("[STEP 0] Loading modules and functions...\n")
cat(strrep("-", 80), "\n")

# Load STLT functions
if(file.exists("R/stlt.R")) {
  source("R/stlt.R")
  cat("  ✓ STLT loaded\n")
}

# Load DSTLT Enhanced
if(file.exists("R (sendiri)/DSTLT_Enhanced.R")) {
  source("R (sendiri)/DSTLT_Enhanced.R")
  cat("  ✓ DSTLT Enhanced loaded\n")
}

# Load other model functions
model_files <- c("R/tlt.R", "R/getqx.R", "R/cbd.R")
for(f in model_files) {
  if(file.exists(f)) {
    tryCatch({
      source(f)
      cat(sprintf("  ✓ %s loaded\n", basename(f)))
    }, error = function(e) {
      cat(sprintf("  ⚠ Could not load %s\n", basename(f)))
    })
  }
}

# Load custom predict functions
if(file.exists("R/predictstlt.R")) source("R/predictstlt.R")
if(file.exists("R/predicttlt.R")) source("R/predicttlt.R")

cat("\n")

# ============================================================================
# 1. LOAD CLEANED DATA
# ============================================================================

cat("[STEP 1] Loading cleaned data...\n")
cat(strrep("-", 80), "\n")

# Check if cleaned data exists
if(!file.exists("data bersih/data_hmd_only.rds")) {
  cat("  ✗ Cleaned data not found. Please run 00_Data_Cleaning_Comprehensive.R first\n")
  stop("Data not found")
}

data_hmd_only <- readRDS("data bersih/data_hmd_only.rds")
data_hmd_cbs <- readRDS("data bersih/data_hmd_cbs.rds")
dstlt_matrices_hmd <- readRDS("data bersih/dstlt_matrices_hmd.rds")
dstlt_matrices_cbs <- readRDS("data bersih/dstlt_matrices_cbs.rds")

cat(sprintf("  ✓ HMD Only: %d rows\n", nrow(data_hmd_only)))
cat(sprintf("  ✓ HMD+CBS: %d rows\n", nrow(data_hmd_cbs)))
cat(sprintf("  ✓ DSTLT matrices loaded\n\n"))

# ============================================================================
# 2. DEFINE HELPER FUNCTIONS
# ============================================================================

cat("[STEP 2] Defining helper functions...\n")
cat(strrep("-", 80), "\n")

#' Fit all single-cohort models (STLT, TLT, etc.)
fit_all_single_cohort_models <- function(cohort_data, verbose = FALSE) {
  models <- list()

  # STLT
  if(verbose) cat("    Fitting STLT...\n")
  tryCatch({
    models$stlt <- stlt(ages = cohort_data$Age, qx = cohort_data$qx)
  }, error = function(e) {
    if(verbose) cat("      ✗ STLT failed\n")
    models$stlt <- NULL
  })

  # TLT (if available)
  if(exists("tlt")) {
    if(verbose) cat("    Fitting TLT...\n")
    tryCatch({
      models$tlt <- tlt(
        ages = cohort_data$Age,
        qx = cohort_data$qx,
        lx = cohort_data$lx,
        dx = cohort_data$dx
      )
    }, error = function(e) {
      if(verbose) cat("      ✗ TLT failed\n")
      models$tlt <- NULL
    })
  }

  # Gompertz, Makeham, HP2 (if available)
  if(exists("get_qx")) {
    for(law in c("gompertz", "makeham", "HP2")) {
      if(verbose) cat(sprintf("    Fitting %s...\n", law))
      tryCatch({
        models[[law]] <- get_qx(
          x = cohort_data$Age,
          qx = cohort_data$qx,
          law = law,
          pred_ages = cohort_data$Age
        )
      }, error = function(e) {
        if(verbose) cat(sprintf("      ✗ %s failed\n", law))
        models[[law]] <- NULL
      })
    }
  }

  return(models)
}

#' Compute SSE for all models
compute_all_sse <- function(cohort_data, models) {
  observed <- cohort_data$qx
  sse_results <- list()

  for(model_name in names(models)) {
    model <- models[[model_name]]

    if(is.null(model)) {
      sse_results[[model_name]] <- NA
      next
    }

    # Get predictions based on model type
    predicted <- tryCatch({
      if(model_name == "stlt" && exists("predict.stlt")) {
        predict(model, newdata = cohort_data$Age)
      } else if(model_name == "tlt" && exists("predict.tlt")) {
        predict.tlt(model, newdata = cohort_data$Age)
      } else if(is.numeric(model)) {
        model  # For Gompertz etc., predictions are already computed
      } else {
        NULL
      }
    }, error = function(e) NULL)

    if(!is.null(predicted) && length(predicted) == length(observed)) {
      sse_results[[model_name]] <- sum((observed - predicted)^2, na.rm = TRUE)
    } else {
      sse_results[[model_name]] <- NA
    }
  }

  return(sse_results)
}

#' Create comparison table
make_comparison_table <- function(sse_results) {
  df <- tibble(
    Method = names(sse_results),
    SSE = unlist(sse_results)
  ) %>%
    filter(!is.na(SSE)) %>%
    arrange(SSE) %>%
    mutate(Rank = row_number())

  return(df)
}

cat("  ✓ Helper functions defined\n\n")

# ============================================================================
# 3. SINGLE COHORT ANALYSIS (Contoh: Cohort 1901)
# ============================================================================

cat("[STEP 3] Single Cohort Analysis (Cohort 1901, Female)...\n")
cat(strrep("-", 80), "\n")

target_cohort <- 1901

# --- HMD Only ---
cat("\n=== Analysis 1: HMD Only ===\n")

cohort_hmd <- data_hmd_only %>%
  filter(Cohort == target_cohort, Sex == "Female")

cat(sprintf("  Data points: %d | Age range: %d-%d\n",
            nrow(cohort_hmd), min(cohort_hmd$Age), max(cohort_hmd$Age)))

models_hmd <- fit_all_single_cohort_models(cohort_hmd, verbose = TRUE)
sse_hmd <- compute_all_sse(cohort_hmd, models_hmd)
comparison_hmd <- make_comparison_table(sse_hmd)

cat("\nComparison Table (HMD Only):\n")
print(comparison_hmd, n = Inf)

# --- HMD+CBS ---
cat("\n=== Analysis 2: HMD+CBS ===\n")

cohort_cbs <- data_hmd_cbs %>%
  filter(Cohort == target_cohort, Sex == "Female")

cat(sprintf("  Data points: %d | Age range: %d-%d\n",
            nrow(cohort_cbs), min(cohort_cbs$Age), max(cohort_cbs$Age)))

models_cbs <- fit_all_single_cohort_models(cohort_cbs, verbose = TRUE)
sse_cbs <- compute_all_sse(cohort_cbs, models_cbs)
comparison_cbs <- make_comparison_table(sse_cbs)

cat("\nComparison Table (HMD+CBS):\n")
print(comparison_cbs, n = Inf)

# ============================================================================
# 4. DSTLT ANALYSIS (Multiple Cohorts)
# ============================================================================

cat("\n[STEP 4] DSTLT Analysis (Multiple Cohorts)...\n")
cat(strrep("-", 80), "\n")

# --- DSTLT dengan HMD Only ---
cat("\n=== DSTLT: HMD Only ===\n")

dstlt_female_hmd <- fit_dstlt_tidy(data_hmd_only, sex = "Female")
dstlt_male_hmd <- fit_dstlt_tidy(data_hmd_only, sex = "Male")

if(!is.null(dstlt_female_hmd)) {
  cat("\nDSTLT Female (HMD Only) Parameters:\n")
  cat(sprintf("  a     = %.6f\n", dstlt_female_hmd$coefficients$a))
  cat(sprintf("  b     = %.6f\n", dstlt_female_hmd$coefficients$b))
  cat(sprintf("  theta = %.6f\n", dstlt_female_hmd$coefficients$theta))
  cat(sprintf("  gamma = %.6f\n", dstlt_female_hmd$coefficients$gamma))
  cat(sprintf("  N     = %d\n", dstlt_female_hmd$coefficients$N))
  cat(sprintf("  Omega = %.2f\n", dstlt_female_hmd$Omega))

  sse_dstlt_f_hmd <- compute_dstlt_sse(dstlt_female_hmd)
  cat(sprintf("  SSE   = %.6f\n", sse_dstlt_f_hmd))
}

if(!is.null(dstlt_male_hmd)) {
  cat("\nDSTLT Male (HMD Only) Parameters:\n")
  cat(sprintf("  a     = %.6f\n", dstlt_male_hmd$coefficients$a))
  cat(sprintf("  b     = %.6f\n", dstlt_male_hmd$coefficients$b))
  cat(sprintf("  theta = %.6f\n", dstlt_male_hmd$coefficients$theta))
  cat(sprintf("  gamma = %.6f\n", dstlt_male_hmd$coefficients$gamma))
  cat(sprintf("  N     = %d\n", dstlt_male_hmd$coefficients$N))
  cat(sprintf("  Omega = %.2f\n", dstlt_male_hmd$Omega))

  sse_dstlt_m_hmd <- compute_dstlt_sse(dstlt_male_hmd)
  cat(sprintf("  SSE   = %.6f\n", sse_dstlt_m_hmd))
}

# --- DSTLT dengan HMD+CBS ---
cat("\n=== DSTLT: HMD+CBS ===\n")

dstlt_female_cbs <- fit_dstlt_tidy(data_hmd_cbs, sex = "Female")
dstlt_male_cbs <- fit_dstlt_tidy(data_hmd_cbs, sex = "Male")

if(!is.null(dstlt_female_cbs)) {
  cat("\nDSTLT Female (HMD+CBS) Parameters:\n")
  cat(sprintf("  a     = %.6f\n", dstlt_female_cbs$coefficients$a))
  cat(sprintf("  b     = %.6f\n", dstlt_female_cbs$coefficients$b))
  cat(sprintf("  theta = %.6f\n", dstlt_female_cbs$coefficients$theta))
  cat(sprintf("  gamma = %.6f\n", dstlt_female_cbs$coefficients$gamma))
  cat(sprintf("  N     = %d\n", dstlt_female_cbs$coefficients$N))
  cat(sprintf("  Omega = %.2f\n", dstlt_female_cbs$Omega))

  sse_dstlt_f_cbs <- compute_dstlt_sse(dstlt_female_cbs)
  cat(sprintf("  SSE   = %.6f\n", sse_dstlt_f_cbs))
}

if(!is.null(dstlt_male_cbs)) {
  cat("\nDSTLT Male (HMD+CBS) Parameters:\n")
  cat(sprintf("  a     = %.6f\n", dstlt_male_cbs$coefficients$a))
  cat(sprintf("  b     = %.6f\n", dstlt_male_cbs$coefficients$b))
  cat(sprintf("  theta = %.6f\n", dstlt_male_cbs$coefficients$theta))
  cat(sprintf("  gamma = %.6f\n", dstlt_male_cbs$coefficients$gamma))
  cat(sprintf("  N     = %d\n", dstlt_male_cbs$coefficients$N))
  cat(sprintf("  Omega = %.2f\n", dstlt_male_cbs$Omega))

  sse_dstlt_m_cbs <- compute_dstlt_sse(dstlt_male_cbs)
  cat(sprintf("  SSE   = %.6f\n", sse_dstlt_m_cbs))
}

# ============================================================================
# 5. VISUALIZATIONS
# ============================================================================

cat("\n[STEP 5] Creating visualizations...\n")
cat(strrep("-", 80), "\n")

# Create output directory for plots
if(!dir.exists("output")) {
  dir.create("output", recursive = TRUE)
}

if(!dir.exists("output/plots")) {
  dir.create("output/plots", recursive = TRUE)
}

# Plot 1: STLT Comparison (HMD vs HMD+CBS)
if(!is.null(models_hmd$stlt) && !is.null(models_cbs$stlt)) {
  cat("  Creating Plot 1: STLT Comparison...\n")

  pred_stlt_hmd <- predict(models_hmd$stlt, newdata = cohort_hmd$Age)
  pred_stlt_cbs <- predict(models_cbs$stlt, newdata = cohort_cbs$Age)

  plot_data <- bind_rows(
    tibble(Age = cohort_hmd$Age, qx = pred_stlt_hmd, Type = "STLT (HMD Only)"),
    tibble(Age = cohort_cbs$Age, qx = pred_stlt_cbs, Type = "STLT (HMD+CBS)")
  )

  obs_data <- cohort_hmd %>%
    select(Age, qx, lx) %>%
    mutate(size = log(lx + 1))

  p1 <- ggplot() +
    geom_point(data = obs_data, aes(x = Age, y = qx, size = size),
               alpha = 0.4, shape = 1, color = "black") +
    geom_line(data = plot_data, aes(x = Age, y = qx, color = Type),
              linewidth = 1) +
    scale_color_manual(values = c("STLT (HMD Only)" = "blue", "STLT (HMD+CBS)" = "red")) +
    scale_size_continuous(guide = "none") +
    labs(
      title = sprintf("STLT Comparison: HMD Only vs HMD+CBS (Cohort %d, Female)", target_cohort),
      x = "Age",
      y = "Mortality Rate (qx)",
      color = ""
    ) +
    theme_bw() +
    theme(legend.position = "top")

  ggsave("output/plots/01_stlt_comparison.png", p1, width = 10, height = 6, dpi = 300)
  cat("    ✓ Saved: output/plots/01_stlt_comparison.png\n")
}

# Plot 2: DSTLT Female vs Male (HMD Only)
if(!is.null(dstlt_female_hmd) && !is.null(dstlt_male_hmd)) {
  cat("  Creating Plot 2: DSTLT Female vs Male...\n")

  p2 <- plot_dstlt_sex_comparison(
    list(female = dstlt_female_hmd, male = dstlt_male_hmd),
    period_index = 5
  )

  ggsave("output/plots/02_dstlt_sex_comparison.png", p2, width = 10, height = 6, dpi = 300)
  cat("    ✓ Saved: output/plots/02_dstlt_sex_comparison.png\n")
}

# Plot 3: DSTLT Cohort Comparison
if(!is.null(dstlt_female_hmd)) {
  cat("  Creating Plot 3: DSTLT Cohort Comparison...\n")

  p3 <- compare_dstlt_cohorts(dstlt_female_hmd, cohort_indices = c(1, 5, 10, 15))

  ggsave("output/plots/03_dstlt_cohorts.png", p3, width = 10, height = 6, dpi = 300)
  cat("    ✓ Saved: output/plots/03_dstlt_cohorts.png\n")
}

# Plot 4: Method Comparison Bar Chart
if(nrow(comparison_hmd) > 0) {
  cat("  Creating Plot 4: Method Comparison Bar Chart...\n")

  p4 <- ggplot(comparison_hmd, aes(x = reorder(Method, SSE), y = SSE, fill = Method)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.4f", SSE)), hjust = -0.1, size = 3) +
    coord_flip() +
    labs(
      title = sprintf("Model Comparison by SSE (Cohort %d, Female, HMD Only)", target_cohort),
      x = "Method",
      y = "Sum of Squared Errors (SSE)"
    ) +
    theme_bw()

  ggsave("output/plots/04_method_comparison_sse.png", p4, width = 10, height = 6, dpi = 300)
  cat("    ✓ Saved: output/plots/04_method_comparison_sse.png\n")
}

# ============================================================================
# 6. SAVE RESULTS
# ============================================================================

cat("\n[STEP 6] Saving analysis results...\n")
cat(strrep("-", 80), "\n")

results <- list(
  single_cohort = list(
    cohort = target_cohort,
    hmd_only = list(
      models = models_hmd,
      sse = sse_hmd,
      comparison = comparison_hmd
    ),
    hmd_cbs = list(
      models = models_cbs,
      sse = sse_cbs,
      comparison = comparison_cbs
    )
  ),
  dstlt = list(
    hmd_only = list(
      female = dstlt_female_hmd,
      male = dstlt_male_hmd
    ),
    hmd_cbs = list(
      female = dstlt_female_cbs,
      male = dstlt_male_cbs
    )
  )
)

saveRDS(results, "output/full_analysis_results.rds")
cat("  ✓ Saved: output/full_analysis_results.rds\n")

# ============================================================================
# 7. FINAL SUMMARY REPORT
# ============================================================================

cat("\n")
cat("================================================================================\n")
cat("  ANALISIS LENGKAP SELESAI\n")
cat("================================================================================\n\n")

cat("SUMMARY:\n\n")

cat("1. Single Cohort Analysis (Cohort 1901, Female):\n")
cat("   a) HMD Only:\n")
if(nrow(comparison_hmd) > 0) {
  cat(sprintf("      Best model: %s (SSE = %.6f)\n",
              comparison_hmd$Method[1], comparison_hmd$SSE[1]))
}
cat("\n   b) HMD+CBS:\n")
if(nrow(comparison_cbs) > 0) {
  cat(sprintf("      Best model: %s (SSE = %.6f)\n",
              comparison_cbs$Method[1], comparison_cbs$SSE[1]))
}

cat("\n2. DSTLT Analysis (Multiple Cohorts):\n")
if(!is.null(dstlt_female_hmd)) {
  cat(sprintf("   Female (HMD Only): N=%d, Omega=%.1f, SSE=%.4f\n",
              dstlt_female_hmd$coefficients$N,
              dstlt_female_hmd$Omega,
              compute_dstlt_sse(dstlt_female_hmd)))
}
if(!is.null(dstlt_male_hmd)) {
  cat(sprintf("   Male (HMD Only): N=%d, Omega=%.1f, SSE=%.4f\n",
              dstlt_male_hmd$coefficients$N,
              dstlt_male_hmd$Omega,
              compute_dstlt_sse(dstlt_male_hmd)))
}

cat("\n3. Visualizations:\n")
cat("   - 4 plots created in output/plots/\n")

cat("\n4. Output Files:\n")
cat("   - output/full_analysis_results.rds\n")
cat("   - output/plots/*.png\n")

cat("\nKEY FINDINGS:\n")
cat("  - STLT model performance with HMD vs HMD+CBS augmentation\n")
cat("  - DSTLT threshold ages (N) estimated for multiple cohorts\n")
cat("  - Gender differences in mortality patterns captured\n")
cat("  - CBS data successfully integrated for advanced ages (>92)\n")

cat("\n")
cat("================================================================================\n\n")
