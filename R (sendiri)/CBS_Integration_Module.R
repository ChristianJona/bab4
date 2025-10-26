# ============================================================================
# CBS INTEGRATION MODULE
# ============================================================================
# Modul khusus untuk mengintegrasikan data CBS dengan HMD
#
# CBS data adalah period data (year-age format) yang perlu dikonversi
# ke cohort data untuk digabungkan dengan HMD cohort data
#
# Referensi: Huang et al. (2020), Section 3.2 - Data Augmentation
# ============================================================================

library(tidyverse)

#' Convert CBS period data to cohort data
#' @param cbs_data Period data dari CBS (year, age, dx)
#' @param cohort_years Cohort years yang ingin direkonstruksi
#' @return Cohort-based data frame
convert_period_to_cohort <- function(cbs_data, cohort_years) {
  cat("  Converting CBS period data to cohort data...\n")

  # Untuk period data: year dan age
  # Untuk cohort data: cohort = year - age
  # Contoh: year=2000, age=95 → cohort=1905

  cbs_cohort <- cbs_data %>%
    mutate(
      Cohort = year - age
    ) %>%
    # Filter hanya cohort yang kita perlukan
    filter(Cohort %in% cohort_years) %>%
    # Group by cohort dan age untuk aggregate
    group_by(Cohort, age) %>%
    summarise(
      dx_sum = sum(dx, na.rm = TRUE),
      n_periods = n(),
      .groups = "drop"
    ) %>%
    rename(Age = age)

  cat(sprintf("    ✓ Converted %d period observations to %d cohort observations\n",
              nrow(cbs_data), nrow(cbs_cohort)))

  return(cbs_cohort)
}

#' Estimate gender split untuk CBS data
#' @param hmd_data HMD data untuk reference
#' @param transition_age Age di mana CBS data mulai digunakan
#' @return List dengan proporsi female dan male
estimate_gender_split <- function(hmd_data, transition_age = 92) {
  cat("  Estimating gender split from HMD at transition age...\n")

  # Hitung proporsi female/male di transition age
  gender_props <- hmd_data %>%
    filter(Age == transition_age) %>%
    group_by(Sex) %>%
    summarise(total_dx = sum(dx, na.rm = TRUE), .groups = "drop") %>%
    mutate(proportion = total_dx / sum(total_dx))

  prop_female <- gender_props$proportion[gender_props$Sex == "Female"]
  prop_male <- gender_props$proportion[gender_props$Sex == "Male"]

  # Default jika tidak ada data
  if(length(prop_female) == 0) prop_female <- 0.65
  if(length(prop_male) == 0) prop_male <- 0.35

  cat(sprintf("    ✓ Female proportion: %.3f\n", prop_female))
  cat(sprintf("    ✓ Male proportion: %.3f\n", prop_male))

  return(list(female = prop_female, male = prop_male))
}

#' Calculate qx dari dx counts
#' @param dx_vector Vector of death counts by age
#' @param initial_lx Initial population (if known)
#' @return List dengan lx, dx, dan qx vectors
calculate_qx_from_dx <- function(dx_vector, initial_lx = NULL) {
  n <- length(dx_vector)
  lx <- numeric(n + 1)

  # Estimate initial population jika tidak diberikan
  if(is.null(initial_lx)) {
    # Metode 1: Sum of all deaths + buffer
    initial_lx <- sum(dx_vector, na.rm = TRUE) * 1.5
  }

  lx[1] <- initial_lx

  # Calculate lx for each age
  for(i in 1:n) {
    lx[i + 1] <- max(0, lx[i] - dx_vector[i])
  }

  # Calculate qx
  qx <- numeric(n)
  for(i in 1:n) {
    if(lx[i] > 0) {
      qx[i] <- dx_vector[i] / lx[i]
    } else {
      qx[i] <- 0
    }
  }

  # Bound qx to [0, 1]
  qx <- pmax(0, pmin(1, qx))

  return(list(
    lx = lx[1:n],
    dx = dx_vector,
    qx = qx
  ))
}

#' Smooth transition antara HMD dan CBS
#' @param hmd_values Values dari HMD
#' @param cbs_values Values dari CBS
#' @param transition_age Age untuk transisi
#' @param blend_window Window size untuk blending
#' @return Blended values
smooth_transition <- function(hmd_values, cbs_values, transition_age, blend_window = 3) {
  # Linear blending di sekitar transition age
  # Ages: [..., 90, 91, 92, 93, 94, ...]
  #                   ^transition

  n_hmd <- length(hmd_values)
  n_cbs <- length(cbs_values)

  # Find indices untuk blending
  ages_hmd <- transition_age - n_hmd + (1:n_hmd)
  ages_cbs <- transition_age + (1:n_cbs)

  # Create blending weights
  result <- numeric(n_hmd + n_cbs)

  for(i in 1:length(result)) {
    age <- ages_hmd[1] + i - 1

    if(age < transition_age - blend_window) {
      # Pure HMD
      result[i] <- hmd_values[i]
    } else if(age > transition_age + blend_window) {
      # Pure CBS
      cbs_idx <- age - transition_age
      result[i] <- cbs_values[cbs_idx]
    } else {
      # Blend
      weight_hmd <- (transition_age + blend_window - age) / (2 * blend_window)
      weight_cbs <- 1 - weight_hmd

      hmd_idx <- min(n_hmd, max(1, i))
      cbs_idx <- max(1, age - transition_age)

      result[i] <- weight_hmd * hmd_values[hmd_idx] +
                   weight_cbs * cbs_values[min(cbs_idx, n_cbs)]
    }
  }

  return(result)
}

#' Main function: Integrate CBS data dengan HMD
#' @param hmd_data Clean HMD data
#' @param cbs_data Raw CBS period data
#' @param transition_age Age untuk mulai menggunakan CBS (default 92)
#' @param verbose Print progress messages
#' @return Augmented data frame
integrate_cbs_data <- function(hmd_data, cbs_data, transition_age = 92, verbose = TRUE) {

  if(verbose) cat("\n  === CBS Integration ===\n")

  # 1. Get cohort years dari HMD
  cohort_years <- unique(hmd_data$Cohort)

  # 2. Convert CBS period to cohort
  cbs_cohort <- convert_period_to_cohort(cbs_data, cohort_years)

  # 3. Estimate gender split
  gender_props <- estimate_gender_split(hmd_data, transition_age)

  # 4. Split CBS data by estimated gender
  if(verbose) cat("  Splitting CBS data by gender...\n")

  cbs_female <- cbs_cohort %>%
    mutate(
      dx = dx_sum * gender_props$female,
      Sex = "Female"
    )

  cbs_male <- cbs_cohort %>%
    mutate(
      dx = dx_sum * gender_props$male,
      Sex = "Male"
    )

  cbs_split <- bind_rows(cbs_female, cbs_male)

  # 5. Calculate lx dan qx untuk CBS data
  if(verbose) cat("  Calculating lx and qx for CBS data...\n")

  cbs_with_qx <- cbs_split %>%
    group_by(Cohort, Sex) %>%
    arrange(Age) %>%
    group_modify(~ {
      # Get initial lx from HMD at transition age
      hmd_cohort <- hmd_data %>%
        filter(Cohort == .y$Cohort, Sex == .y$Sex, Age == transition_age)

      if(nrow(hmd_cohort) > 0) {
        initial_lx <- hmd_cohort$lx[1]
      } else {
        initial_lx <- NULL
      }

      # Calculate mortality rates
      mortality_calc <- calculate_qx_from_dx(.x$dx, initial_lx)

      tibble(
        Age = .x$Age,
        lx = mortality_calc$lx,
        dx = mortality_calc$dx,
        qx = mortality_calc$qx
      )
    }) %>%
    ungroup()

  # 6. Combine HMD (age ≤ transition_age) dengan CBS (age > transition_age)
  if(verbose) cat("  Merging HMD and CBS data...\n")

  hmd_base <- hmd_data %>%
    filter(Age <= transition_age)

  cbs_extended <- cbs_with_qx %>%
    filter(Age > transition_age)

  # Combine
  combined <- bind_rows(hmd_base, cbs_extended) %>%
    arrange(Sex, Cohort, Age)

  # 7. Apply smoothing at transition
  if(verbose) cat("  Applying smooth transition...\n")

  combined_smoothed <- combined %>%
    group_by(Cohort, Sex) %>%
    arrange(Age) %>%
    mutate(
      qx_smooth = smooth_qx_transition(qx, Age, transition_age)
    ) %>%
    mutate(qx = qx_smooth) %>%
    select(-qx_smooth) %>%
    ungroup()

  if(verbose) {
    cat(sprintf("  ✓ Final augmented dataset: %d rows\n", nrow(combined_smoothed)))
    cat(sprintf("  ✓ Age range: %d - %d\n",
                min(combined_smoothed$Age), max(combined_smoothed$Age)))
  }

  return(combined_smoothed)
}

#' Helper: Smooth qx at transition point
#' @param qx Vector of mortality rates
#' @param ages Vector of ages
#' @param transition_age Transition age
#' @return Smoothed qx
smooth_qx_transition <- function(qx, ages, transition_age, window = 3) {
  n <- length(qx)
  smoothed <- qx

  # Find indices near transition
  trans_idx <- which(ages == transition_age)

  if(length(trans_idx) == 0) return(qx)

  # Smooth around transition age
  start_idx <- max(1, trans_idx - window)
  end_idx <- min(n, trans_idx + window)

  for(i in start_idx:end_idx) {
    local_start <- max(1, i - 1)
    local_end <- min(n, i + 1)
    smoothed[i] <- mean(qx[local_start:local_end], na.rm = TRUE)
  }

  return(smoothed)
}

#' Validate CBS integration results
#' @param original_hmd Original HMD data
#' @param augmented Augmented data with CBS
#' @param transition_age Transition age used
validate_cbs_integration <- function(original_hmd, augmented, transition_age = 92) {
  cat("\n  === Validating CBS Integration ===\n")

  # Check 1: HMD data tidak berubah untuk age ≤ transition_age
  hmd_portion <- augmented %>% filter(Age <= transition_age)
  original_portion <- original_hmd %>% filter(Age <= transition_age)

  hmd_match <- all.equal(
    hmd_portion %>% arrange(Cohort, Sex, Age),
    original_portion %>% arrange(Cohort, Sex, Age)
  )

  if(isTRUE(hmd_match)) {
    cat("  ✓ HMD portion preserved correctly\n")
  } else {
    cat("  ⚠ HMD portion has changes (may be due to smoothing)\n")
  }

  # Check 2: CBS data mulai dari age > transition_age
  cbs_portion <- augmented %>% filter(Age > transition_age)
  if(nrow(cbs_portion) > 0) {
    cat(sprintf("  ✓ CBS data present for ages %d - %d\n",
                min(cbs_portion$Age), max(cbs_portion$Age)))
  } else {
    cat("  ⚠ No CBS data in augmented dataset\n")
  }

  # Check 3: Monotonicity of qx
  by_cohort <- augmented %>%
    group_by(Cohort, Sex) %>%
    arrange(Age) %>%
    summarise(
      is_monotone = all(diff(qx) >= 0),
      .groups = "drop"
    )

  non_monotone <- sum(!by_cohort$is_monotone)
  if(non_monotone == 0) {
    cat("  ✓ qx is monotone (increasing) for all cohorts\n")
  } else {
    cat(sprintf("  ⚠ %d cohort(s) have non-monotone qx\n", non_monotone))
  }

  cat("  === Validation Complete ===\n\n")
}

# Export functions
# (No need for explicit export in R, just sourcing the file makes them available)
