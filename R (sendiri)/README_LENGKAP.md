# Analisis Life Table dengan Metode STLT, DSTLT, dan CBS

## Overview

Repositori ini berisi implementasi lengkap untuk analisis life table dengan metode advanced, termasuk:

- **STLT** (Smooth Threshold Life Table) - untuk single cohort
- **DSTLT** (Dynamic Smooth Threshold Life Table) - untuk multiple cohorts
- **CBS Integration** - augmentasi data HMD dengan CBS untuk usia >92
- **Metode pembanding**: TLT, Gompertz, Makeham, HP2, CBD

Berdasarkan paper:
> Huang, H., Browne, C., & Xu, Y. (2020). Modelling life tables with advanced ages: An extreme value theory approach. *Insurance: Mathematics and Economics*.

## Struktur Folder

```
bab4/
├── data/                          # Data mentah
│   ├── Female Belanda.csv         # HMD data female
│   ├── Male Belanda.csv           # HMD data male
│   └── cbs_augmented_format (1).csv  # CBS data usia >92
│
├── data bersih/                   # Data yang sudah dibersihkan (output)
│   ├── data_hmd_only.rds         # Dataset 1: HMD saja
│   ├── data_hmd_cbs.rds          # Dataset 2: HMD+CBS augmented
│   ├── dstlt_matrices_hmd.rds    # Matrices untuk DSTLT (HMD)
│   ├── dstlt_matrices_cbs.rds    # Matrices untuk DSTLT (CBS)
│   └── validation_summary.rds    # Hasil validasi
│
├── R/                             # Fungsi-fungsi dari jurnal
│   ├── stlt.R                    # STLT implementation
│   ├── dstlt.R                   # DSTLT implementation
│   ├── tlt.R                     # TLT implementation
│   ├── cbd.R                     # CBD implementation
│   └── ... (fungsi lainnya)
│
├── R (sendiri)/                   # Script analisis utama (JALANKAN INI!)
│   ├── 00_Data_Cleaning_Comprehensive.R   # [1] Pembersihan data
│   ├── CBS_Integration_Module.R           # [Module] Integrasi CBS
│   ├── DSTLT_Enhanced.R                   # [Module] DSTLT enhanced
│   ├── 03_CBS_Specific_Analysis.R         # [2] Analisis CBS
│   ├── 05_Full_Analysis_All_Methods.R     # [3] Analisis lengkap
│   └── README_LENGKAP.md                  # Dokumentasi ini
│
└── output/                        # Hasil analisis (dibuat otomatis)
    ├── plots/                     # Visualisasi
    ├── cbs_analysis/             # Hasil analisis CBS
    └── full_analysis_results.rds # Hasil lengkap
```

## Prerequisites

### R Packages yang Diperlukan

```r
install.packages(c(
  "tidyverse",      # Data manipulation & visualization
  "lubridate",      # Date handling
  "gridExtra",      # Multiple plots
  "knitr",          # Tables
  "MortalityLaws",  # Mortality modeling (optional)
  "scales"          # Plotting scales
))
```

### Data Files

Pastikan file-file berikut ada di folder `data/`:
1. `Female Belanda.csv` - HMD data untuk perempuan Belanda
2. `Male Belanda.csv` - HMD data untuk laki-laki Belanda
3. `cbs_augmented_format (1).csv` - CBS data untuk usia lanjut

## Cara Menggunakan

### Step 1: Pembersihan Data

**Script**: `00_Data_Cleaning_Comprehensive.R`

Script ini membersihkan data mentah dan membuat 3 dataset:
- HMD Only (right-censored at age 100)
- HMD+CBS (augmented with CBS for age >92)
- DSTLT Matrices (format matrix untuk analisis DSTLT)

```r
# Jalankan di R Console atau RStudio
source("R (sendiri)/00_Data_Cleaning_Comprehensive.R")
```

**Output**:
- `data bersih/data_hmd_only.rds`
- `data bersih/data_hmd_cbs.rds`
- `data bersih/dstlt_matrices_hmd.rds`
- `data bersih/dstlt_matrices_cbs.rds`
- `data bersih/validation_summary.rds`

**Features**:
- ✅ Validasi data otomatis (check NA, outliers, monotonicity)
- ✅ Quality checks untuk setiap dataset
- ✅ Filtering cohorts (1893-1908) dan ages (≥65)
- ✅ Summary statistics

### Step 2: Analisis Khusus CBS (Opsional)

**Script**: `03_CBS_Specific_Analysis.R`

Script ini menganalisis data CBS secara mendalam:
- Period to cohort reconstruction
- Comparison HMD vs CBS at transition age
- Mortality patterns di usia lanjut
- Visualisasi lengkap

```r
source("R (sendiri)/03_CBS_Specific_Analysis.R")
```

**Output**:
- `output/cbs_analysis/` folder dengan 5 plots
- `output/cbs_analysis/cbs_analysis_results.rds`

**Visualizations**:
1. CBS Deaths by Age
2. CBS Deaths Over Time
3. Heatmap (Year × Age)
4. HMD vs CBS Comparison
5. Mortality Pattern Comparison

### Step 3: Analisis Lengkap Semua Metode

**Script**: `05_Full_Analysis_All_Methods.R`

Script utama yang menjalankan analisis lengkap:
- Single cohort analysis (STLT, TLT, Gompertz, dll)
- Multi-cohort analysis (DSTLT)
- Female vs Male comparison
- HMD Only vs HMD+CBS comparison

```r
source("R (sendiri)/05_Full_Analysis_All_Methods.R")
```

**Output**:
- `output/plots/` folder dengan 4 plots utama
- `output/full_analysis_results.rds`

**Visualizations**:
1. STLT Comparison (HMD vs HMD+CBS)
2. DSTLT Female vs Male
3. DSTLT Cohort Comparison
4. Method Comparison Bar Chart (SSE)

## Modules

### CBS_Integration_Module.R

Modul khusus untuk integrasi data CBS dengan HMD.

**Fungsi utama**:
- `integrate_cbs_data()` - Main integration function
- `convert_period_to_cohort()` - Convert CBS period to cohort
- `estimate_gender_split()` - Split CBS by gender
- `calculate_qx_from_dx()` - Calculate mortality rates
- `smooth_transition()` - Smooth HMD-CBS transition

**Usage**:
```r
source("R (sendiri)/CBS_Integration_Module.R")

augmented_data <- integrate_cbs_data(
  hmd_data = data_hmd,
  cbs_data = data_cbs,
  transition_age = 92,
  verbose = TRUE
)
```

### DSTLT_Enhanced.R

Modul enhanced untuk DSTLT dengan plot dan predict functions.

**Fungsi utama**:
- `fit_dstlt_tidy()` - Fit DSTLT dengan tidy data format
- `plot_dstlt_model()` - Plot DSTLT fitted model
- `predict_dstlt_model()` - Predict mortality rates
- `compare_dstlt_cohorts()` - Compare multiple cohorts
- `fit_dstlt_both_sexes()` - Fit untuk female & male
- `plot_dstlt_sex_comparison()` - Plot female vs male

**Usage**:
```r
source("R (sendiri)/DSTLT_Enhanced.R")

# Fit model
model_female <- fit_dstlt_tidy(data, sex = "Female")

# Plot
plot_dstlt_model(model_female, period_index = 1)

# Predict
predictions <- predict_dstlt_model(model_female, newdata = 65:110, t = 1)
```

## Metodologi

### 1. STLT (Smooth Threshold Life Table)

Model untuk single cohort yang membagi mortality pattern menjadi:
- **Gompertz component** (age < N): Exponential increase
- **Pareto component** (age ≥ N): Power law decay

**Parameters**:
- B: Gompertz baseline mortality
- C: Gompertz aging rate
- γ (gamma): Shape parameter Pareto
- N: Threshold age
- Ω (omega): Maximum lifespan

**Formula**:
```
For x < N:  S(x) = exp(-B/log(C) * (C^x - 1))
For x ≥ N:  S(x) = (1 + γ(x-N)/θ)^(-1/γ)
```

### 2. DSTLT (Dynamic Smooth Threshold Life Table)

Extends STLT untuk multiple cohorts dengan time-varying parameters:

**Time-varying component**:
- B(t) = exp(a + b*t)
- Parameters a, b capture mortality improvements over time

**Advantages**:
- Captures cohort effects
- Models mortality improvements
- Pooled estimation (more stable)

### 3. CBS Integration

**Why CBS data?**
- HMD data kurang reliable di usia >100
- CBS Netherlands has high-quality data for ages 92-107
- Augmentation improves parameter estimation

**Integration approach**:
1. Use HMD for ages ≤ 92
2. Use CBS for ages > 92
3. Smooth transition at age 92
4. Gender split based on HMD patterns at age 92

## Results Interpretation

### SSE (Sum of Squared Errors)

Lower SSE = better fit

**Typical SSE values**:
- STLT: 0.0001 - 0.001
- DSTLT: 0.0005 - 0.005 (higher due to multiple cohorts)
- Gompertz: 0.001 - 0.01
- CBD: varies widely

### Threshold Age (N)

**Expected ranges**:
- Female: 92-96
- Male: 89-93
- With CBS augmentation: may be slightly higher

**Interpretation**:
- N = age di mana mortality pattern beralih dari Gompertz ke Pareto
- Higher N = longer "normal" aging period
- Lower N = earlier onset of "extreme age" dynamics

### Omega (Ω) - Maximum Lifespan

**Formula**: Ω = N - θ/γ

**Expected ranges**:
- Female: 108-112
- Male: 105-109

**Interpretation**:
- Upper bound of possible lifespan
- Should be finite (realistic)
- CBS augmentation typically increases Ω estimates

## Troubleshooting

### Error: "File not found"

**Solution**:
```r
# Check working directory
getwd()

# Set working directory ke folder bab4
setwd("path/to/bab4")
```

### Error: "Package not found"

**Solution**:
```r
# Install missing packages
install.packages("tidyverse")
```

### Warning: "DSTLT failed to converge"

**Possible causes**:
- Insufficient data (< 3 cohorts)
- Too few age points per cohort
- Starting values inappropriate

**Solution**:
- Adjust startN and endN in fit_dstlt_tidy()
- Check data quality
- Try different cohort selections

### Error: "CBS file not found"

Script akan skip CBS augmentation dan menggunakan HMD saja.

**To fix**:
- Download CBS data
- Place in `data/` folder
- Re-run data cleaning script

## Advanced Usage

### Analisis Custom Cohort

```r
# Load data
data <- readRDS("data bersih/data_hmd_only.rds")

# Filter specific cohort
cohort_1900 <- data %>% filter(Cohort == 1900, Sex == "Female")

# Fit STLT
model <- stlt(ages = cohort_1900$Age, qx = cohort_1900$qx)

# Inspect parameters
model$coefficients

# Predict
predictions <- predict(model, newdata = 65:110)
```

### Comparison Multiple Cohorts

```r
# Fit DSTLT
models <- fit_dstlt_both_sexes(data)

# Compare cohorts
plot_dstlt_cohorts(models$female, cohort_indices = c(1, 5, 10, 15))

# Extract SSE
sse_female <- compute_dstlt_sse(models$female)
sse_male <- compute_dstlt_sse(models$male)
```

### Custom CBS Integration

```r
source("R (sendiri)/CBS_Integration_Module.R")

# Custom transition age
augmented <- integrate_cbs_data(
  hmd_data = data_hmd,
  cbs_data = data_cbs,
  transition_age = 95,  # Use 95 instead of 92
  verbose = TRUE
)
```

## FAQ

**Q: Perbedaan STLT dan DSTLT?**

A: STLT untuk single cohort, DSTLT untuk multiple cohorts dengan time-varying parameters.

**Q: Kapan menggunakan HMD Only vs HMD+CBS?**

A:
- HMD Only: Kalau fokus pada ages ≤100
- HMD+CBS: Kalau ingin analisis hingga age 105-108

**Q: Berapa lama running time?**

A:
- Data cleaning: ~30 detik
- CBS analysis: ~1 menit
- Full analysis: ~2-5 menit (tergantung DSTLT convergence)

**Q: Bagaimana cara mengubah cohort range?**

A: Edit di script `00_Data_Cleaning_Comprehensive.R`, baris 71-76:
```r
filter(
  Cohort >= 1893,  # Ubah ini
  Cohort <= 1908,  # Dan ini
  Age >= 65
)
```

## Referensi

1. Huang, H., Browne, C., & Xu, Y. (2020). Modelling life tables with advanced ages: An extreme value theory approach. *Insurance: Mathematics and Economics*, 96, 1-12.

2. Cairns, A. J., Blake, D., & Dowd, K. (2006). A two-factor model for stochastic mortality with parameter uncertainty: theory and calibration. *Journal of Risk and Insurance*, 73(4), 687-718.

3. Human Mortality Database. University of California, Berkeley (USA), and Max Planck Institute for Demographic Research (Germany). Available at www.mortality.org.

4. Statistics Netherlands (CBS). Available at www.cbs.nl.

## Contact

Untuk pertanyaan atau issues:
- Buat issue di GitHub repository
- Email: [your-email]

---

**Last Updated**: 2025-10-26

**Version**: 1.0.0
