# ============================================================================
# PEMBERSIHAN DATA KOMPREHENSIF UNTUK ANALISIS STLT/DSTLT/CBS
# Berdasarkan: Huang et al. (2020) - Modelling life tables with advanced ages
# ============================================================================
#
# Script ini membersihkan dan memvalidasi data untuk analisis lengkap:
# 1. HMD only: Data dari Human Mortality Database (right-censored di usia 100)
# 2. HMD+CBS: Data HMD (usia ≤92) digabung dengan CBS (usia >92)
# 3. Full dataset: Untuk analisis DSTLT dengan multiple cohorts
#
# Output: 3 file RDS di folder "data bersih/"
# ============================================================================

library(tidyverse)
library(lubridate)

cat("\n")
cat("================================================================================\n")
cat("  PEMBERSIHAN DATA KOMPREHENSIF UNTUK ANALISIS LIFE TABLE\n")
cat("================================================================================\n\n")

# ============================================================================
# FUNGSI HELPER UNTUK VALIDASI DATA
# ============================================================================

#' Validasi data mortality
#' @param data Data frame dengan kolom Age, qx, dx, lx
#' @param dataset_name Nama dataset untuk reporting
validate_mortality_data <- function(data, dataset_name = "Dataset") {
  cat(sprintf("\n--- Validasi Data: %s ---\n", dataset_name))

  errors <- c()
  warnings <- c()

  # 1. Check required columns
  required_cols <- c("Age", "qx", "dx", "lx")
  missing_cols <- setdiff(required_cols, names(data))
  if(length(missing_cols) > 0) {
    errors <- c(errors, sprintf("Missing columns: %s", paste(missing_cols, collapse=", ")))
  }

  # 2. Check for NA values
  na_counts <- sapply(data[required_cols], function(x) sum(is.na(x)))
  if(any(na_counts > 0)) {
    warnings <- c(warnings, sprintf("NA values found - Age:%d, qx:%d, dx:%d, lx:%d",
                                    na_counts[1], na_counts[2], na_counts[3], na_counts[4]))
  }

  # 3. Check qx range (should be 0-1)
  if("qx" %in% names(data)) {
    invalid_qx <- sum(data$qx < 0 | data$qx > 1, na.rm = TRUE)
    if(invalid_qx > 0) {
      errors <- c(errors, sprintf("%d values of qx outside [0,1]", invalid_qx))
    }
  }

  # 4. Check monotonicity of Age
  if("Age" %in% names(data) && "Cohort" %in% names(data)) {
    by_cohort <- data %>% group_by(Cohort, Sex) %>%
      summarise(is_monotone = all(diff(Age) == 1), .groups = "drop")

    non_monotone <- sum(!by_cohort$is_monotone)
    if(non_monotone > 0) {
      warnings <- c(warnings, sprintf("%d cohort(s) have non-monotone ages", non_monotone))
    }
  }

  # 5. Check for negative deaths
  if("dx" %in% names(data)) {
    neg_dx <- sum(data$dx < 0, na.rm = TRUE)
    if(neg_dx > 0) {
      errors <- c(errors, sprintf("%d negative death counts", neg_dx))
    }
  }

  # 6. Check lx consistency (should be decreasing)
  if("lx" %in% names(data) && "Cohort" %in% names(data)) {
    by_cohort <- data %>% group_by(Cohort, Sex) %>%
      summarise(is_decreasing = all(diff(lx) <= 0), .groups = "drop")

    non_decreasing <- sum(!by_cohort$is_decreasing)
    if(non_decreasing > 0) {
      warnings <- c(warnings, sprintf("%d cohort(s) have non-decreasing lx", non_decreasing))
    }
  }

  # Report results
  if(length(errors) > 0) {
    cat("  ✗ ERRORS:\n")
    for(e in errors) cat(sprintf("    - %s\n", e))
  }

  if(length(warnings) > 0) {
    cat("  ⚠ WARNINGS:\n")
    for(w in warnings) cat(sprintf("    - %s\n", w))
  }

  if(length(errors) == 0 && length(warnings) == 0) {
    cat("  ✓ All validations passed\n")
  }

  return(list(valid = length(errors) == 0, errors = errors, warnings = warnings))
}

#' Remove outliers menggunakan IQR method
#' @param data Data mortality
#' @param column Column name untuk outlier detection
remove_outliers_iqr <- function(data, column = "qx") {
  Q1 <- quantile(data[[column]], 0.25, na.rm = TRUE)
  Q3 <- quantile(data[[column]], 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1

  lower_bound <- Q1 - 3 * IQR_val
  upper_bound <- Q3 + 3 * IQR_val

  outliers <- data[[column]] < lower_bound | data[[column]] > upper_bound
  n_outliers <- sum(outliers, na.rm = TRUE)

  if(n_outliers > 0) {
    cat(sprintf("  Removed %d outliers from %s (bounds: [%.4f, %.4f])\n",
                n_outliers, column, lower_bound, upper_bound))
  }

  return(data[!outliers, ])
}

#' Smooth qx using moving average
#' @param qx Vector of mortality rates
#' @param window Window size for moving average
smooth_qx <- function(qx, window = 3) {
  n <- length(qx)
  smoothed <- numeric(n)

  for(i in 1:n) {
    start_idx <- max(1, i - floor(window/2))
    end_idx <- min(n, i + floor(window/2))
    smoothed[i] <- mean(qx[start_idx:end_idx], na.rm = TRUE)
  }

  return(smoothed)
}

# ============================================================================
# 1. LOAD DATA HMD (Human Mortality Database)
# ============================================================================

cat("\n[STEP 1] Loading HMD Data\n")
cat(strrep("-", 80), "\n")

# Check if files exist
path_females <- "data/Female Belanda.csv"
path_males   <- "data/Male Belanda.csv"

if(!file.exists(path_females)) stop("File not found: ", path_females)
if(!file.exists(path_males)) stop("File not found: ", path_males)

# Baca data FEMALES dari HMD
data_females_hmd <- read_csv2(path_females, show_col_types = FALSE) %>%
  mutate(Sex = "Female") %>%
  rename(Cohort = Year)

# Baca data MALES dari HMD
data_males_hmd <- read_csv2(path_males, show_col_types = FALSE) %>%
  mutate(Sex = "Male") %>%
  rename(Cohort = Year)

cat(sprintf("  ✓ Female HMD: %d rows\n", nrow(data_females_hmd)))
cat(sprintf("  ✓ Male HMD: %d rows\n", nrow(data_males_hmd)))

# Gabungkan
data_hmd_raw <- bind_rows(data_females_hmd, data_males_hmd)

# ============================================================================
# 2. CLEAN HMD DATA
# ============================================================================

cat("\n[STEP 2] Cleaning HMD Data\n")
cat(strrep("-", 80), "\n")

# Konversi tipe data dan clean
data_hmd_clean <- data_hmd_raw %>%
  mutate(
    Cohort = as.integer(Cohort),
    Age = as.numeric(Age),
    lx  = as.numeric(lx),
    dx  = as.numeric(dx),
    qx  = as.numeric(qx)
  ) %>%
  # Remove rows with missing critical values
  filter(!is.na(Age), !is.na(Cohort)) %>%
  # Ensure qx is in valid range
  mutate(
    qx = pmax(0, pmin(1, qx))
  ) %>%
  # Filter sesuai jurnal: cohorts 1893-1908, age ≥65
  filter(
    Cohort >= 1893,
    Cohort <= 1908,
    Age >= 65
  ) %>%
  select(Cohort, Age, Sex, lx, dx, qx) %>%
  arrange(Sex, Cohort, Age)

cat(sprintf("  ✓ Cleaned HMD data: %d rows\n", nrow(data_hmd_clean)))
cat(sprintf("  ✓ Cohort range: %d - %d\n", min(data_hmd_clean$Cohort), max(data_hmd_clean$Cohort)))
cat(sprintf("  ✓ Age range: %d - %d\n", min(data_hmd_clean$Age), max(data_hmd_clean$Age)))

# Validate HMD data
validation_hmd <- validate_mortality_data(data_hmd_clean, "HMD Clean")

# ============================================================================
# 3. LOAD DATA CBS (voor usia >92)
# ============================================================================

cat("\n[STEP 3] Loading CBS Data\n")
cat(strrep("-", 80), "\n")

path_cbs <- "data/cbs_augmented_format (1).csv"
if(!file.exists(path_cbs)) {
  cat("  ⚠ CBS file not found, skipping CBS augmentation\n")
  data_cbs_processed <- NULL
} else {
  # Baca data CBS
  data_cbs_raw <- read_delim(path_cbs, delim = ";", show_col_types = FALSE) %>%
    filter(!is.na(age), !is.na(dx)) %>%
    mutate(
      year = as.integer(year),
      age = as.integer(age),
      dx = as.integer(dx)
    )

  cat(sprintf("  ✓ CBS data: %d rows\n", nrow(data_cbs_raw)))
  cat(sprintf("  ✓ Year range: %d - %d\n", min(data_cbs_raw$year), max(data_cbs_raw$year)))
  cat(sprintf("  ✓ Age range: %d - %d\n", min(data_cbs_raw$age), max(data_cbs_raw$age)))

  # Store for later processing
  data_cbs_processed <- data_cbs_raw
}

# ============================================================================
# 4. CREATE DATASET 1: HMD ONLY (right-censored at 100)
# ============================================================================

cat("\n[STEP 4] Creating Dataset 1: HMD Only (right-censored at 100)\n")
cat(strrep("-", 80), "\n")

data_hmd_only <- data_hmd_clean %>%
  filter(Age <= 100) %>%
  group_by(Cohort, Sex) %>%
  arrange(Age) %>%
  ungroup()

# Quality check: Remove cohorts with too few observations
min_ages <- 20  # Minimal 20 age points per cohort
cohort_counts <- data_hmd_only %>%
  group_by(Cohort, Sex) %>%
  summarise(n_ages = n(), .groups = "drop") %>%
  filter(n_ages >= min_ages)

data_hmd_only <- data_hmd_only %>%
  semi_join(cohort_counts, by = c("Cohort", "Sex"))

cat(sprintf("  ✓ HMD Only dataset: %d rows\n", nrow(data_hmd_only)))
cat(sprintf("  ✓ Age range: %d - %d\n", min(data_hmd_only$Age), max(data_hmd_only$Age)))
cat(sprintf("  ✓ Cohorts: %d female, %d male\n",
            length(unique(data_hmd_only$Cohort[data_hmd_only$Sex == "Female"])),
            length(unique(data_hmd_only$Cohort[data_hmd_only$Sex == "Male"]))))

# Validate
validation_hmd_only <- validate_mortality_data(data_hmd_only, "HMD Only")

# ============================================================================
# 5. CREATE DATASET 2: HMD+CBS (augmented, usia >92 dari CBS)
# ============================================================================

cat("\n[STEP 5] Creating Dataset 2: HMD+CBS (augmented with CBS for age >92)\n")
cat(strrep("-", 80), "\n")

if(is.null(data_cbs_processed)) {
  cat("  ⚠ No CBS data available, using HMD only\n")
  data_hmd_cbs <- data_hmd_clean
} else {
  # Source khusus CBS integration function
  source("R (sendiri)/CBS_Integration_Module.R", local = TRUE)

  # Integrate CBS data
  data_hmd_cbs <- integrate_cbs_data(
    hmd_data = data_hmd_clean,
    cbs_data = data_cbs_processed,
    transition_age = 92,
    verbose = TRUE
  )
}

cat(sprintf("  ✓ HMD+CBS dataset: %d rows\n", nrow(data_hmd_cbs)))
cat(sprintf("  ✓ Age range: %d - %d\n", min(data_hmd_cbs$Age), max(data_hmd_cbs$Age)))

# Validate
validation_hmd_cbs <- validate_mortality_data(data_hmd_cbs, "HMD+CBS")

# ============================================================================
# 6. CREATE DATASET 3: FULL MATRIX untuk DSTLT (multiple cohorts)
# ============================================================================

cat("\n[STEP 6] Creating Dataset 3: Full Matrix for DSTLT\n")
cat(strrep("-", 80), "\n")

# Untuk DSTLT, kita perlu data dalam format matrix
# Columns = cohorts, Rows = ages

create_dstlt_matrices <- function(data, max_age = 100) {
  # Pivot untuk mendapatkan format wide
  age_matrix <- data %>%
    select(Cohort, Age, Sex, qx) %>%
    pivot_wider(
      names_from = Cohort,
      values_from = qx,
      id_cols = c(Age, Sex)
    ) %>%
    arrange(Sex, Age)

  # Split by sex
  female_data <- age_matrix %>% filter(Sex == "Female") %>% select(-Sex)
  male_data <- age_matrix %>% filter(Sex == "Male") %>% select(-Sex)

  # Convert to matrices
  female_ages <- as.matrix(female_data[, 1])
  female_qx <- as.matrix(female_data[, -1])

  male_ages <- as.matrix(male_data[, 1])
  male_qx <- as.matrix(male_data[, -1])

  # Create age matrices (replicate age column for each cohort)
  female_age_matrix <- matrix(female_ages,
                              nrow = nrow(female_qx),
                              ncol = ncol(female_qx))

  male_age_matrix <- matrix(male_ages,
                            nrow = nrow(male_qx),
                            ncol = ncol(male_qx))

  return(list(
    female = list(ages = female_age_matrix, qx = female_qx),
    male = list(ages = male_age_matrix, qx = male_qx)
  ))
}

# Create matrices for HMD only
dstlt_matrices_hmd <- create_dstlt_matrices(data_hmd_only)

# Create matrices for HMD+CBS
dstlt_matrices_cbs <- create_dstlt_matrices(data_hmd_cbs)

cat(sprintf("  ✓ Female matrix (HMD): %d ages × %d cohorts\n",
            nrow(dstlt_matrices_hmd$female$qx),
            ncol(dstlt_matrices_hmd$female$qx)))
cat(sprintf("  ✓ Male matrix (HMD): %d ages × %d cohorts\n",
            nrow(dstlt_matrices_hmd$male$qx),
            ncol(dstlt_matrices_hmd$male$qx)))

# ============================================================================
# 7. SUMMARY STATISTICS
# ============================================================================

cat("\n[STEP 7] Computing Summary Statistics\n")
cat(strrep("-", 80), "\n")

compute_summary_stats <- function(data, name) {
  stats <- data %>%
    group_by(Sex) %>%
    summarise(
      n_cohorts = n_distinct(Cohort),
      n_ages = n_distinct(Age),
      age_min = min(Age),
      age_max = max(Age),
      qx_mean = mean(qx, na.rm = TRUE),
      qx_median = median(qx, na.rm = TRUE),
      qx_sd = sd(qx, na.rm = TRUE),
      .groups = "drop"
    )

  return(list(name = name, stats = stats))
}

summary_hmd_only <- compute_summary_stats(data_hmd_only, "HMD Only")
summary_hmd_cbs <- compute_summary_stats(data_hmd_cbs, "HMD+CBS")

# Print summary
print_summary <- function(summary_obj) {
  cat(sprintf("\n%s:\n", summary_obj$name))
  print(summary_obj$stats, n = Inf)
}

print_summary(summary_hmd_only)
print_summary(summary_hmd_cbs)

# ============================================================================
# 8. SAVE HASIL
# ============================================================================

cat("\n[STEP 8] Saving Results\n")
cat(strrep("-", 80), "\n")

# Buat folder jika belum ada
if(!dir.exists("data bersih")) {
  dir.create("data bersih", recursive = TRUE)
  cat("  ✓ Created directory: data bersih/\n")
}

# Save datasets
saveRDS(data_hmd_only, "data bersih/data_hmd_only.rds")
cat("  ✓ Saved: data bersih/data_hmd_only.rds\n")

saveRDS(data_hmd_cbs, "data bersih/data_hmd_cbs.rds")
cat("  ✓ Saved: data bersih/data_hmd_cbs.rds\n")

saveRDS(dstlt_matrices_hmd, "data bersih/dstlt_matrices_hmd.rds")
cat("  ✓ Saved: data bersih/dstlt_matrices_hmd.rds\n")

saveRDS(dstlt_matrices_cbs, "data bersih/dstlt_matrices_cbs.rds")
cat("  ✓ Saved: data bersih/dstlt_matrices_cbs.rds\n")

# Save validation results
validation_summary <- list(
  hmd_only = validation_hmd_only,
  hmd_cbs = validation_hmd_cbs,
  summary_stats = list(
    hmd_only = summary_hmd_only,
    hmd_cbs = summary_hmd_cbs
  )
)

saveRDS(validation_summary, "data bersih/validation_summary.rds")
cat("  ✓ Saved: data bersih/validation_summary.rds\n")

# ============================================================================
# 9. FINAL REPORT
# ============================================================================

cat("\n")
cat("================================================================================\n")
cat("  PEMBERSIHAN DATA SELESAI\n")
cat("================================================================================\n\n")

cat("DATASET SUMMARY:\n\n")

cat("1. HMD Only (right-censored at age 100):\n")
cat(sprintf("   - Rows: %d\n", nrow(data_hmd_only)))
cat(sprintf("   - Age range: %d - %d\n",
            min(data_hmd_only$Age), max(data_hmd_only$Age)))
cat(sprintf("   - Female cohorts: %d\n",
            length(unique(data_hmd_only$Cohort[data_hmd_only$Sex == "Female"]))))
cat(sprintf("   - Male cohorts: %d\n",
            length(unique(data_hmd_only$Cohort[data_hmd_only$Sex == "Male"]))))
cat(sprintf("   - Validation: %s\n\n",
            ifelse(validation_hmd_only$valid, "✓ PASSED", "✗ FAILED")))

cat("2. HMD+CBS (augmented with CBS for age >92):\n")
cat(sprintf("   - Rows: %d\n", nrow(data_hmd_cbs)))
cat(sprintf("   - Age range: %d - %d\n",
            min(data_hmd_cbs$Age), max(data_hmd_cbs$Age)))
cat(sprintf("   - Female cohorts: %d\n",
            length(unique(data_hmd_cbs$Cohort[data_hmd_cbs$Sex == "Female"]))))
cat(sprintf("   - Male cohorts: %d\n",
            length(unique(data_hmd_cbs$Cohort[data_hmd_cbs$Sex == "Male"]))))
cat(sprintf("   - Validation: %s\n\n",
            ifelse(validation_hmd_cbs$valid, "✓ PASSED", "✗ FAILED")))

cat("3. DSTLT Matrices:\n")
cat(sprintf("   - Female: %d ages × %d cohorts\n",
            nrow(dstlt_matrices_hmd$female$qx),
            ncol(dstlt_matrices_hmd$female$qx)))
cat(sprintf("   - Male: %d ages × %d cohorts\n\n",
            nrow(dstlt_matrices_hmd$male$qx),
            ncol(dstlt_matrices_hmd$male$qx)))

cat("USAGE:\n")
cat("  data_hmd_only <- readRDS('data bersih/data_hmd_only.rds')\n")
cat("  data_hmd_cbs <- readRDS('data bersih/data_hmd_cbs.rds')\n")
cat("  dstlt_matrices <- readRDS('data bersih/dstlt_matrices_hmd.rds')\n\n")

cat("NEXT STEPS:\n")
cat("  - Run 03_CBS_Specific_Analysis.R for CBS-specific methods\n")
cat("  - Run 04_DSTLT_Analysis.R for DSTLT analysis\n")
cat("  - Run 05_Full_Analysis.R for comprehensive comparison\n\n")

cat("================================================================================\n\n")
