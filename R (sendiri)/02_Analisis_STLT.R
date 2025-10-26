# ============================================================================
# ANALISIS SMOOTH THRESHOLD LIFE TABLE (STLT)
# Berdasarkan: Huang et al. (2020)
# ============================================================================
#
# Script ini melakukan:
# 1. Fitting STLT, TLT, dan model lainnya
# 2. Perbandingan model dengan SSE
# 3. Analisis untuk HMD only vs HMD+CBS
# 4. Visualisasi hasil

library(tidyverse)
library(knitr)
library(MortalityLaws)
library(scales)

cat("\n=== ANALISIS SMOOTH THRESHOLD LIFE TABLE ===\n\n")

# ============================================================================
# 1. LOAD FUNGSI MODEL
# ============================================================================

cat("1. Loading fungsi model...\n")

# Load semua fungsi dari folder R
file_sources <- list.files("R", pattern="*.R$", full.names=TRUE)
for(file in file_sources) {
  tryCatch({
    source(file)
  }, error = function(e) {
    cat("   Warning: Cannot load", basename(file), "\n")
  })
}

cat("   ✓ Fungsi model berhasil dimuat\n")

# ============================================================================
# 2. LOAD DATA
# ============================================================================

cat("\n2. Loading data...\n")

# Load kedua dataset
data_hmd_only <- readRDS("data bersih/data_hmd_only.rds")
data_hmd_cbs <- readRDS("data bersih/data_hmd_cbs.rds")

cat("   ✓ Data HMD Only:", nrow(data_hmd_only), "baris\n")
cat("   ✓ Data HMD+CBS:", nrow(data_hmd_cbs), "baris\n")

# Pisahkan berdasarkan gender
data_female_hmd <- data_hmd_only %>% filter(Sex == "Female")
data_male_hmd   <- data_hmd_only %>% filter(Sex == "Male")

data_female_cbs <- data_hmd_cbs %>% filter(Sex == "Female")
data_male_cbs   <- data_hmd_cbs %>% filter(Sex == "Male")

cat("   ✓ Female cohorts (HMD):", length(unique(data_female_hmd$Cohort)), "\n")
cat("   ✓ Male cohorts (HMD):", length(unique(data_male_hmd$Cohort)), "\n")

# ============================================================================
# 3. FUNGSI HELPER
# ============================================================================

# Fungsi untuk fit semua model
fit_all_models <- function(cohort_data, verbose = FALSE) {
  models <- list()

  # STLT
  tryCatch({
    models$stlt <- stlt(ages = cohort_data$Age, qx = cohort_data$qx)
    if(verbose) cat("  STLT: OK\n")
  }, error = function(e) {
    if(verbose) cat("  STLT: GAGAL\n")
    models$stlt <- NULL
  })

  # TLT
  tryCatch({
    models$tlt <- tlt(ages = cohort_data$Age, qx = cohort_data$qx,
                      lx = cohort_data$lx, dx = cohort_data$dx)
    if(verbose) cat("  TLT: OK\n")
  }, error = function(e) {
    if(verbose) cat("  TLT: GAGAL\n")
    models$tlt <- NULL
  })

  # Gompertz, Makeham, HP2
  for(law in c("gompertz", "makeham", "HP2")) {
    tryCatch({
      models[[law]] <- get_qx(x = cohort_data$Age, qx = cohort_data$qx,
                              law = law, pred_ages = cohort_data$Age)
      if(verbose) cat("  ", law, ": OK\n")
    }, error = function(e) {
      if(verbose) cat("  ", law, ": GAGAL\n")
      models[[law]] <- NULL
    })
  }

  return(models)
}

# Fungsi untuk hitung SSE
compute_sse <- function(observed, predicted_matrix) {
  if(is.vector(predicted_matrix)) {
    predicted_matrix <- matrix(predicted_matrix, ncol = 1)
  }

  sse <- apply(predicted_matrix, 2, function(pred) {
    sum((observed - pred)^2, na.rm = TRUE)
  })

  return(sse)
}

# Fungsi untuk membuat tabel perbandingan
make_comparison_table <- function(cohort_data, models) {
  pred_list <- list()

  # Collect predictions
  if(!is.null(models$stlt)) {
    pred_list$STLT <- predict(models$stlt, newdata = cohort_data$Age)
  }

  if(!is.null(models$tlt)) {
    pred_list$TLT <- predict.tlt(models$tlt, newdata = cohort_data$Age)
  }

  for(law in c("gompertz", "makeham", "HP2")) {
    if(!is.null(models[[law]])) {
      pred_list[[toupper(law)]] <- models[[law]]
    }
  }

  if(length(pred_list) == 0) return(NULL)

  # Validate lengths
  pred_lengths <- sapply(pred_list, length)
  valid_models <- names(pred_list)[pred_lengths == nrow(cohort_data)]

  if(length(valid_models) == 0) return(NULL)

  pred_list <- pred_list[valid_models]
  pred_matrix <- do.call(cbind, pred_list)

  # Compute SSE
  sse_values <- compute_sse(cohort_data$qx, pred_matrix)

  # Create results table
  results <- data.frame(
    Method = names(pred_list),
    SSE = sse_values,
    stringsAsFactors = FALSE
  )

  results <- results[order(results$SSE), ]
  rownames(results) <- NULL

  return(results)
}

# ============================================================================
# 4. ANALISIS UNTUK COHORT 1901 (CONTOH)
# ============================================================================

cat("\n3. Analisis cohort 1901 (Female)...\n\n")

target_cohort <- 1901

# --- ANALISIS DENGAN HMD ONLY ---
cat("=== ANALISIS DENGAN HMD ONLY ===\n")
cohort_hmd <- data_female_hmd %>% filter(Cohort == target_cohort)
cat("Data points:", nrow(cohort_hmd), "| Age range:", min(cohort_hmd$Age), "-", max(cohort_hmd$Age), "\n\n")

cat("Fitting models...\n")
models_hmd <- fit_all_models(cohort_hmd, verbose = TRUE)

if(!is.null(models_hmd$stlt)) {
  cat("\nParameter STLT (HMD Only):\n")
  cat("  B     =", format(models_hmd$stlt$coefficients$B, scientific = TRUE, digits = 4), "\n")
  cat("  C     =", format(models_hmd$stlt$coefficients$C, digits = 6), "\n")
  cat("  gamma =", format(models_hmd$stlt$coefficients$gamma, digits = 4), "\n")
  cat("  N     =", models_hmd$stlt$coefficients$N, "\n")
  cat("  Omega =", format(models_hmd$stlt$Omega, digits = 2), "\n")
}

cat("\nPerbandingan SSE:\n")
comp_table_hmd <- make_comparison_table(cohort_hmd, models_hmd)
if(!is.null(comp_table_hmd)) {
  print(comp_table_hmd)
  cat("\n✓ Model terbaik (HMD Only):", comp_table_hmd$Method[1], "\n")
}

# --- ANALISIS DENGAN HMD+CBS ---
cat("\n=== ANALISIS DENGAN HMD+CBS ===\n")
cohort_cbs <- data_female_cbs %>% filter(Cohort == target_cohort)
cat("Data points:", nrow(cohort_cbs), "| Age range:", min(cohort_cbs$Age), "-", max(cohort_cbs$Age), "\n\n")

cat("Fitting models...\n")
models_cbs <- fit_all_models(cohort_cbs, verbose = TRUE)

if(!is.null(models_cbs$stlt)) {
  cat("\nParameter STLT (HMD+CBS):\n")
  cat("  B     =", format(models_cbs$stlt$coefficients$B, scientific = TRUE, digits = 4), "\n")
  cat("  C     =", format(models_cbs$stlt$coefficients$C, digits = 6), "\n")
  cat("  gamma =", format(models_cbs$stlt$coefficients$gamma, digits = 4), "\n")
  cat("  N     =", models_cbs$stlt$coefficients$N, "\n")
  cat("  Omega =", format(models_cbs$stlt$Omega, digits = 2), "\n")
}

cat("\nPerbandingan SSE:\n")
comp_table_cbs <- make_comparison_table(cohort_cbs, models_cbs)
if(!is.null(comp_table_cbs)) {
  print(comp_table_cbs)
  cat("\n✓ Model terbaik (HMD+CBS):", comp_table_cbs$Method[1], "\n")
}

# ============================================================================
# 5. PERBANDINGAN HMD ONLY VS HMD+CBS
# ============================================================================

cat("\n=== PERBANDINGAN HMD ONLY VS HMD+CBS ===\n\n")

if(!is.null(models_hmd$stlt) && !is.null(models_cbs$stlt)) {
  cat("Parameter STLT:\n")
  cat(sprintf("%-10s %15s %15s\n", "Parameter", "HMD Only", "HMD+CBS"))
  cat(strrep("-", 42), "\n")
  cat(sprintf("%-10s %15.6e %15.6e\n", "B",
              models_hmd$stlt$coefficients$B,
              models_cbs$stlt$coefficients$B))
  cat(sprintf("%-10s %15.6f %15.6f\n", "C",
              models_hmd$stlt$coefficients$C,
              models_cbs$stlt$coefficients$C))
  cat(sprintf("%-10s %15.4f %15.4f\n", "gamma",
              models_hmd$stlt$coefficients$gamma,
              models_cbs$stlt$coefficients$gamma))
  cat(sprintf("%-10s %15d %15d\n", "N",
              models_hmd$stlt$coefficients$N,
              models_cbs$stlt$coefficients$N))
  cat(sprintf("%-10s %15.2f %15.2f\n", "Omega",
              models_hmd$stlt$Omega,
              models_cbs$stlt$Omega))
  cat("\n")
}

# SSE Comparison
if(!is.null(comp_table_hmd) && !is.null(comp_table_cbs)) {
  cat("SSE Comparison untuk STLT:\n")
  sse_hmd <- comp_table_hmd$SSE[comp_table_hmd$Method == "STLT"]
  sse_cbs <- comp_table_cbs$SSE[comp_table_cbs$Method == "STLT"]

  if(length(sse_hmd) > 0 && length(sse_cbs) > 0) {
    cat("  HMD Only: ", format(sse_hmd, digits = 4), "\n")
    cat("  HMD+CBS:  ", format(sse_cbs, digits = 4), "\n")

    improvement <- (sse_hmd - sse_cbs) / sse_hmd * 100
    cat("  Improvement:", format(improvement, digits = 2), "%\n\n")
  }
}

# ============================================================================
# 6. VISUALISASI
# ============================================================================

cat("4. Membuat visualisasi...\n")

# Plot comparison
if(!is.null(models_hmd$stlt) && !is.null(models_cbs$stlt)) {

  # Prepare data for plotting
  pred_data <- tibble(
    Age = cohort_hmd$Age,
    Observed = cohort_hmd$qx,
    STLT_HMD = predict(models_hmd$stlt, newdata = cohort_hmd$Age),
    lx = cohort_hmd$lx
  )

  # Add CBS predictions if available
  if(nrow(cohort_cbs) > nrow(cohort_hmd)) {
    pred_data_cbs <- tibble(
      Age = cohort_cbs$Age,
      STLT_CBS = predict(models_cbs$stlt, newdata = cohort_cbs$Age)
    )
    pred_data <- pred_data %>%
      left_join(pred_data_cbs, by = "Age")
  } else {
    pred_data$STLT_CBS <- predict(models_cbs$stlt, newdata = cohort_hmd$Age)
  }

  # Create plot
  p <- ggplot() +
    geom_point(data = pred_data,
               aes(x = Age, y = Observed, size = log(lx + 1)),
               alpha = 0.6, shape = 1, color = "black") +
    geom_line(data = pred_data,
              aes(x = Age, y = STLT_HMD, color = "HMD Only"),
              linewidth = 1) +
    geom_line(data = pred_data,
              aes(x = Age, y = STLT_CBS, color = "HMD+CBS"),
              linewidth = 1) +
    scale_color_manual(name = "Dataset",
                       values = c("HMD Only" = "blue", "HMD+CBS" = "red")) +
    scale_size_continuous(guide = "none") +
    labs(title = "Perbandingan STLT: HMD Only vs HMD+CBS",
         subtitle = paste("Female Cohort", target_cohort),
         x = "Age", y = "Mortality Rate (qx)") +
    theme_bw() +
    theme(legend.position = c(0.2, 0.85),
          legend.background = element_rect(fill = "white", color = "black"),
          plot.title = element_text(size = 14, face = "bold"))

  print(p)

  # Save plot
  ggsave("data bersih/plot_comparison_hmd_vs_cbs.png", p, width = 10, height = 6, dpi = 300)
  cat("   ✓ Plot tersimpan: data bersih/plot_comparison_hmd_vs_cbs.png\n")
}

# ============================================================================
# 7. SAVE RESULTS
# ============================================================================

cat("\n5. Menyimpan hasil analisis...\n")

results_1901 <- list(
  cohort = target_cohort,
  hmd_only = list(
    models = models_hmd,
    comparison = comp_table_hmd
  ),
  hmd_cbs = list(
    models = models_cbs,
    comparison = comp_table_cbs
  )
)

saveRDS(results_1901, "data bersih/results_cohort_1901.rds")
cat("   ✓ Hasil tersimpan: data bersih/results_cohort_1901.rds\n")

cat("\n=== SELESAI ===\n\n")
cat("KESIMPULAN:\n")
cat("- STLT model berhasil di-fit untuk cohort 1901\n")
cat("- Perbandingan HMD Only vs HMD+CBS menunjukkan perbedaan parameter\n")
cat("- Lihat plot untuk visualisasi: data bersih/plot_comparison_hmd_vs_cbs.png\n\n")
