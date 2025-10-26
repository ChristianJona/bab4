# ============================================================================
# DSTLT ENHANCED MODULE
# ============================================================================
# Modul lengkap untuk Dynamic Smooth Threshold Life Table (DSTLT)
#
# Berdasarkan: Huang et al. (2020), Section 2.3 - DSTLT Model
#
# Features:
# - Wrapper functions yang mudah digunakan
# - Plot functions yang sudah diaktifkan
# - Predict functions dengan berbagai options
# - Comparison tools
# ============================================================================

library(tidyverse)

# Source original DSTLT implementation
source("R/dstlt.R", local = TRUE)

#' Fit DSTLT model dengan data tidyverse format
#' @param data Data frame dengan kolom: Cohort, Age, Sex, qx
#' @param sex Gender untuk analisis ("Female" atau "Male")
#' @param startN Starting threshold age untuk testing
#' @param endN Ending threshold age untuk testing
#' @param verbose Print progress
#' @return DSTLT model object
fit_dstlt_tidy <- function(data, sex = "Female", startN = 80, endN = 105, verbose = TRUE) {

  if(verbose) cat(sprintf("\n=== Fitting DSTLT Model (%s) ===\n", sex))

  # Filter by sex
  data_sex <- data %>% filter(Sex == sex)

  # Get cohorts
  cohorts <- unique(data_sex$Cohort) %>% sort()
  n_cohorts <- length(cohorts)

  if(verbose) {
    cat(sprintf("  Cohorts: %d (range %d-%d)\n",
                n_cohorts, min(cohorts), max(cohorts)))
  }

  # Prepare matrices
  # Each column = one cohort
  # Rows = ages

  # Find common age range
  age_ranges <- data_sex %>%
    group_by(Cohort) %>%
    summarise(
      min_age = min(Age),
      max_age = max(Age),
      .groups = "drop"
    )

  common_min_age <- max(age_ranges$min_age)
  common_max_age <- min(age_ranges$max_age)

  if(verbose) {
    cat(sprintf("  Common age range: %d - %d\n", common_min_age, common_max_age))
  }

  # Create matrices
  ages_matrix <- matrix(NA, nrow = 200, ncol = n_cohorts)
  qx_matrix <- matrix(NA, nrow = 200, ncol = n_cohorts)

  for(j in 1:n_cohorts) {
    cohort_data <- data_sex %>%
      filter(Cohort == cohorts[j]) %>%
      arrange(Age)

    n_ages <- nrow(cohort_data)
    ages_matrix[1:n_ages, j] <- cohort_data$Age
    qx_matrix[1:n_ages, j] <- cohort_data$qx
  }

  # Fit DSTLT
  if(verbose) cat("  Fitting model...\n")

  tryCatch({
    model <- dstlt(
      ages = ages_matrix,
      qxs = qx_matrix,
      startN = startN,
      endN = endN,
      hessian = FALSE
    )

    if(verbose) {
      cat("  ✓ Model fitted successfully\n")
      cat(sprintf("    N (threshold age) = %d\n", model$coefficients$N))
      cat(sprintf("    Omega (max age) = %.2f\n", model$Omega))
    }

    # Add metadata
    model$sex <- sex
    model$cohorts <- cohorts
    model$n_cohorts <- n_cohorts

    return(model)

  }, error = function(e) {
    cat("  ✗ Error fitting DSTLT:\n")
    cat(sprintf("    %s\n", e$message))
    return(NULL)
  })
}

#' Plot DSTLT model
#' @param dstlt_model DSTLT model object from fit_dstlt_tidy
#' @param period_index Which period to plot (1 = first cohort, etc.)
#' @param age_range Age range untuk plot
#' @return ggplot object
plot_dstlt_model <- function(dstlt_model, period_index = 1, age_range = c(65, 110)) {

  if(is.null(dstlt_model)) {
    stop("DSTLT model is NULL")
  }

  # Extract parameters
  a <- dstlt_model$coefficients$a
  b <- dstlt_model$coefficients$b
  theta <- dstlt_model$coefficients$theta
  gamma <- dstlt_model$coefficients$gamma
  N <- dstlt_model$coefficients$N
  start <- dstlt_model$Start
  taus <- dstlt_model$Taus
  omega <- dstlt_model$Omega

  # Generate prediction ages
  ages <- seq(age_range[1], age_range[2], by = 0.5)

  # Predict for the specified period
  pred_qx <- predict_dstlt_model(dstlt_model, newdata = ages, t = period_index)

  # Get observed data for this period
  obs_qx <- dstlt_model$qxs[, period_index]
  obs_ages <- dstlt_model$ages[, period_index]

  # Remove NAs
  valid_obs <- !is.na(obs_qx) & !is.na(obs_ages)
  obs_qx <- obs_qx[valid_obs]
  obs_ages <- obs_ages[valid_obs]

  # Create plot data
  plot_data <- tibble(
    Age = ages,
    qx_pred = pred_qx,
    Type = "DSTLT Prediction"
  )

  obs_data <- tibble(
    Age = obs_ages,
    qx_obs = obs_qx
  )

  # Create plot
  p <- ggplot() +
    geom_line(data = plot_data,
              aes(x = Age, y = qx_pred, color = Type),
              linewidth = 1) +
    geom_point(data = obs_data,
               aes(x = Age, y = qx_obs),
               alpha = 0.6, size = 2, shape = 1) +
    geom_vline(xintercept = N, linetype = "dashed", color = "red", alpha = 0.5) +
    annotate("text", x = N, y = 0.1, label = sprintf("N=%d", N),
             angle = 90, vjust = -0.5, color = "red") +
    scale_color_manual(values = c("DSTLT Prediction" = "blue")) +
    labs(
      title = sprintf("DSTLT Model - %s, Period %d", dstlt_model$sex, period_index),
      subtitle = sprintf("N=%d, Omega=%.1f, γ=%.3f",
                         N, omega, gamma),
      x = "Age",
      y = "Mortality Rate (qx)",
      color = ""
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      plot.title = element_text(size = 14, face = "bold")
    )

  return(p)
}

#' Predict from DSTLT model
#' @param dstlt_model DSTLT model object
#' @param newdata Vector of ages untuk prediction
#' @param t Period index (1, 2, 3, ...)
#' @return Vector of predicted qx
predict_dstlt_model <- function(dstlt_model, newdata, t = 1) {

  # Extract parameters
  a <- dstlt_model$coefficients$a
  b <- dstlt_model$coefficients$b
  theta <- dstlt_model$coefficients$theta
  gamma <- dstlt_model$coefficients$gamma
  N <- dstlt_model$coefficients$N
  start <- dstlt_model$Start

  # Compute B and C for this period
  B <- exp(a + b * t)
  C <- 1 / (theta * B)^(1/N)

  # Predict qx for each age
  qx <- numeric(length(newdata))

  for(i in 1:length(newdata)) {
    x <- newdata[i]

    if(x < N) {
      # Gompertz component
      numerator <- exp(-B/log(C) * (C^x - 1)) - exp(-B/log(C) * (C^(x+1) - 1))
      denominator <- exp(-B/log(C) * (C^x - 1))

      if(denominator > 0) {
        qx[i] <- numerator / denominator
      } else {
        qx[i] <- 0
      }

    } else if(x < N + 1) {
      # Transition region
      term1 <- -1 + exp(-B/log(C) * (C^x - 1))
      term2 <- 1 - (exp(-B/log(C) * (C^N - 1))) * (1 + gamma * ((x+1) - N) / theta)^(-1/gamma)
      denominator <- exp(-B/log(C) * (C^x - 1))

      if(denominator > 0) {
        qx[i] <- (term1 + term2) / denominator
      } else {
        qx[i] <- 0
      }

    } else {
      # Pareto component
      omega <- N - theta / gamma

      if(x < (omega - 1)) {
        numerator <- (1 + gamma * (x - N) / theta)^(-1/gamma) -
                     (1 + gamma * ((x+1) - N) / theta)^(-1/gamma)
        denominator <- (1 + gamma * (x - N) / theta)^(-1/gamma)

        if(denominator > 0) {
          qx[i] <- numerator / denominator
        } else {
          qx[i] <- 0
        }

      } else if(x < omega) {
        qx[i] <- 1
      } else {
        qx[i] <- NA
      }
    }
  }

  # Bound to [0, 1]
  qx <- pmax(0, pmin(1, qx, na.rm = TRUE))

  return(qx)
}

#' Compare DSTLT across multiple cohorts
#' @param dstlt_model DSTLT model object
#' @param cohort_indices Which cohorts to compare
#' @param age_range Age range for comparison
#' @return ggplot object
compare_dstlt_cohorts <- function(dstlt_model, cohort_indices = 1:3, age_range = c(65, 110)) {

  if(is.null(dstlt_model)) {
    stop("DSTLT model is NULL")
  }

  # Prepare data
  plot_data <- list()

  for(t in cohort_indices) {
    if(t > dstlt_model$n_cohorts) next

    ages <- seq(age_range[1], age_range[2], by = 0.5)
    pred_qx <- predict_dstlt_model(dstlt_model, newdata = ages, t = t)

    cohort_year <- dstlt_model$cohorts[t]

    plot_data[[t]] <- tibble(
      Age = ages,
      qx = pred_qx,
      Cohort = as.character(cohort_year)
    )
  }

  combined_data <- bind_rows(plot_data)

  # Create plot
  p <- ggplot(combined_data, aes(x = Age, y = qx, color = Cohort)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = dstlt_model$coefficients$N,
               linetype = "dashed", alpha = 0.5) +
    labs(
      title = sprintf("DSTLT Comparison Across Cohorts (%s)", dstlt_model$sex),
      subtitle = sprintf("Threshold age N = %d", dstlt_model$coefficients$N),
      x = "Age",
      y = "Mortality Rate (qx)"
    ) +
    theme_bw() +
    theme(legend.position = "right")

  return(p)
}

#' Compute SSE untuk DSTLT model
#' @param dstlt_model DSTLT model object
#' @return SSE value
compute_dstlt_sse <- function(dstlt_model) {

  n_periods <- dstlt_model$n_cohorts
  sse_total <- 0

  for(t in 1:n_periods) {
    # Get observed data
    obs_qx <- dstlt_model$qxs[, t]
    obs_ages <- dstlt_model$ages[, t]

    # Remove NAs
    valid <- !is.na(obs_qx) & !is.na(obs_ages)
    obs_qx <- obs_qx[valid]
    obs_ages <- obs_ages[valid]

    # Predict
    pred_qx <- predict_dstlt_model(dstlt_model, newdata = obs_ages, t = t)

    # Compute SSE for this period
    sse_t <- sum((obs_qx - pred_qx)^2, na.rm = TRUE)
    sse_total <- sse_total + sse_t
  }

  return(sse_total)
}

#' Create summary table untuk DSTLT model
#' @param dstlt_model DSTLT model object
#' @return Data frame dengan summary statistics
summarise_dstlt <- function(dstlt_model) {

  if(is.null(dstlt_model)) {
    return(NULL)
  }

  # Extract parameters
  coefs <- dstlt_model$coefficients

  summary_df <- tibble(
    Parameter = c("a", "b", "theta", "gamma", "N", "Omega"),
    Value = c(
      coefs$a,
      coefs$b,
      coefs$theta,
      coefs$gamma,
      coefs$N,
      dstlt_model$Omega
    )
  )

  # Add SSE
  sse <- compute_dstlt_sse(dstlt_model)
  summary_df <- bind_rows(
    summary_df,
    tibble(Parameter = "SSE", Value = sse)
  )

  # Add model info
  summary_df <- bind_rows(
    summary_df,
    tibble(
      Parameter = c("Sex", "N_Cohorts", "Age_Start"),
      Value = c(NA, NA, NA)
    )
  )

  return(summary_df)
}

# ============================================================================
# WRAPPER FUNCTIONS UNTUK KEMUDAHAN PENGGUNAAN
# ============================================================================

#' Fit DSTLT untuk kedua gender
#' @param data Data frame dengan HMD/CBS data
#' @param startN Starting threshold age
#' @param endN Ending threshold age
#' @return List dengan model female dan male
fit_dstlt_both_sexes <- function(data, startN = 80, endN = 105) {

  cat("\n=== Fitting DSTLT for Both Sexes ===\n")

  # Fit female
  model_female <- fit_dstlt_tidy(data, sex = "Female", startN = startN, endN = endN)

  # Fit male
  model_male <- fit_dstlt_tidy(data, sex = "Male", startN = startN, endN = endN)

  return(list(
    female = model_female,
    male = model_male
  ))
}

#' Plot comparison antara DSTLT female dan male
#' @param models List dari fit_dstlt_both_sexes
#' @param period_index Which period to compare
#' @param age_range Age range
#' @return ggplot object
plot_dstlt_sex_comparison <- function(models, period_index = 1, age_range = c(65, 110)) {

  # Generate predictions
  ages <- seq(age_range[1], age_range[2], by = 0.5)

  pred_female <- predict_dstlt_model(models$female, newdata = ages, t = period_index)
  pred_male <- predict_dstlt_model(models$male, newdata = ages, t = period_index)

  plot_data <- bind_rows(
    tibble(Age = ages, qx = pred_female, Sex = "Female"),
    tibble(Age = ages, qx = pred_male, Sex = "Male")
  )

  # Create plot
  p <- ggplot(plot_data, aes(x = Age, y = qx, color = Sex)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = models$female$coefficients$N,
               linetype = "dashed", color = "red", alpha = 0.3) +
    geom_vline(xintercept = models$male$coefficients$N,
               linetype = "dashed", color = "blue", alpha = 0.3) +
    scale_color_manual(values = c("Female" = "red", "Male" = "blue")) +
    labs(
      title = sprintf("DSTLT Comparison: Female vs Male (Period %d)", period_index),
      subtitle = sprintf("N_female=%d, N_male=%d",
                         models$female$coefficients$N,
                         models$male$coefficients$N),
      x = "Age",
      y = "Mortality Rate (qx)"
    ) +
    theme_bw()

  return(p)
}

cat("✓ DSTLT Enhanced Module loaded\n")
cat("  Available functions:\n")
cat("  - fit_dstlt_tidy(): Fit DSTLT dengan tidy data\n")
cat("  - plot_dstlt_model(): Plot DSTLT model\n")
cat("  - predict_dstlt_model(): Predict dari DSTLT\n")
cat("  - compare_dstlt_cohorts(): Compare multiple cohorts\n")
cat("  - fit_dstlt_both_sexes(): Fit untuk female dan male\n")
cat("  - plot_dstlt_sex_comparison(): Plot comparison by sex\n\n")
