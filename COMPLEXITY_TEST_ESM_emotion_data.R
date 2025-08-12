
# ---- CONTENTS ----

# 1 PREPARATION
## 1.1 Working Directory & Packages
## 1.2 Global Parameters & Variables
## 1.3 Load Helper Functions
## 1.4 Define Tests
## 1.5 Import Data

# 2 ANALYSIS I – BRINGMANN ET AL. (2013)
## 2.1 Data Preparation & Descriptives
## 2.2 Step 1 – Significance Tests
## 2.3 Step 2 – ACF Analysis
## 2.4 Step 3 – TV-AR
## 2.5 Step 4 – Change Point Analysis
## 2.6 Step 5 – Forecast Skill
## 2.7 Step 6 – Multifractality
## 2.8 Summary Statistics

# 3 ANALYSIS II – FRIED ET AL. (2021)
## 3.1 Data Preparation & Descriptives
## 3.2 Step 1 – Significance Tests
## 3.3 Step 2 – ACF Analysis
## 3.4 Step 3 – TV-AR
## 3.5 Step 4 – Change Point Analysis
## 3.6 Step 5 – Forecast Skill
## 3.7 Step 6 – Multifractality
## 3.8 Summary Statistics

# 4 BETWEEN-DATASET COMPARISONS
## 4.1 Proportion Test
## 4.2 Change Point Comparison
## 4.3 Prediction Decay Comparison
## 4.4 Partial Autocorrelation Comparison
## 4.5 Interpretation

# 5 SAVE WORKSPACE

# -----------------

# ==========================================================
# 1 CONTENT PREPARATION
# ==========================================================

# --- 1.1 Working Directory & Packages --- 

rm(list = ls())
setwd('C:/Users/bbe0870/Desktop/PhD1/complexity_in_emotion_ESM_data')
setwd('C:/Users/Sebastian Küppers/Desktop/Formal Theory of Co-Occuring Emotions (DFG project)/_PhD/_PhD_Study_1/complexity_in_emotion_ESM_data')

library(nonlinearTseries)
library(randtests)
library(tseries)
library(mice)
library(pwr)
library(stats)
library(mgcv)
library(ecp)
library(rEDM)


# --- 1.2 Global Parameters & Variables ---

set.seed(4168)

# GLOBAL VARS
MIN_NUM_OBS <- 20
ALPHA_LEVEL <- 0.05


# --- 1.3 Load Helper Functions ---

# Import helper functions for multifractality 
# (see https://osf.io/nm9b4/, Kelty-Stephen et al., 2023)
source('CJspectrum.R')
source('IAAFT.R')


# --- 1.4 Define Tests ---

tests <- c("Bartels", 
           "PAFC.n",
           "PAFC.max_lag",
           "TV-AR.EDF",
           "TV-AR.p",
           "KPSS.level", 
           "CP.n",
           "KPSS.trend", 
           "Keenan", 
           "Tsay", 
           "pred_decay",
           "multifractality")

significance_tests <- c("Bartels", 
                        "TV-AR.p", 
                        "KPSS.level",
                        "KPSS.trend",
                        "Keenan",
                        "Tsay",
                        "multifractality")


# --- 1.5 Import Data ---

# (see https://osf.io/nca2u/, Olthof et al., 2020)
data_Bringmann2013 <- readRDS("data_Bringmann2013.RDS")
data_Fried2021 <- readRDS("data_Fried2021.RDS")


# ==========================================================
# 2 ANALYSIS I – BRINGMANN ET AL. (2013)
# ==========================================================

# --- 2.1 Data Preparation & Descriptives ---

cols.Bringmann2013 <- c("Worried", "Sad", "Anxious", "Relaxed", "Excited")
subjects.Bringmann2013 <- unique(data_Bringmann2013$subj_id)

res_Bringmann2013 <- list()

for (subj in subjects.Bringmann2013) {
  res_df <- data.frame(matrix(NA,
                              nrow = length(tests),
                              ncol = length(cols.Bringmann2013),
                              dimnames = list(tests, cols.Bringmann2013)))
  res_Bringmann2013[[toString(subj)]] <- res_df
}

# Descriptives
ts_lengths <- c()
for (subj in subjects.Bringmann2013) {
  data <- data_Bringmann2013[data_Bringmann2013$subj_id == subj, cols.Bringmann2013]
  for (col in cols.Bringmann2013) {
    ts <- na.omit(data[[col]])
    N <- length(ts)
    
    # append to ts_lengths
    ts_lengths <- c(ts_lengths, N)
  }
}

mean(ts_lengths) #90.17538


# --- 2.2 Step 1 – Significance Tests

for (subj in subjects.Bringmann2013) {
  data <- data_Bringmann2013[data_Bringmann2013$subj_id == subj, cols.Bringmann2013]
  for (col in cols.Bringmann2013) {
    ts <- na.omit(data[[col]])
    N <- sum(!is.na(ts)) 
    if (N < MIN_NUM_OBS) next
    
    res_df <- res_Bringmann2013[[toString(subj)]]
    
    bartels_p <- tryCatch(bartels.rank.test(ts, alternative = "two.sided")$p.value, error = function(e) NA)
    res_df["Bartels", col] <- if (is.null(bartels_p) || length(bartels_p) == 0) NA else bartels_p
    
    kpss_level_p <- tryCatch(kpss.test(ts, null = "Level", lshort = TRUE)$p.value, error = function(e) NA)
    res_df["KPSS.level", col] <- if (is.null(kpss_level_p) || length(kpss_level_p) == 0) NA else kpss_level_p
    
    kpss_trend_p <- tryCatch(kpss.test(ts, null = "Trend", lshort = TRUE)$p.value, error = function(e) NA)
    res_df["KPSS.trend", col] <- if (is.null(kpss_trend_p) || length(kpss_trend_p) == 0) NA else kpss_trend_p
    
    keenan_p <- tryCatch(keenanTest(ts)$p.value, error = function(e) NA)
    res_df["Keenan", col] <- if (is.null(keenan_p) || length(keenan_p) == 0) NA else keenan_p
    
    tsay_p <- tryCatch(tsayTest(ts)$p.value, error = function(e) NA)
    res_df["Tsay", col] <- if (is.null(tsay_p) || length(tsay_p) == 0) NA else tsay_p
    
    res_Bringmann2013[[toString(subj)]] <- res_df
  }
}


# --- 2.3 Step 2 – ACF Analysis ---

for (subj in subjects.Bringmann2013) {
  data <- data_Bringmann2013[data_Bringmann2013$subj_id == subj, cols.Bringmann2013]
  for (col in cols.Bringmann2013) {
    ts <- na.omit(data[[col]])
    N <- sum(!is.na(ts)) 
    if (N < MIN_NUM_OBS) next
    
    res_df <- res_Bringmann2013[[toString(subj)]]
    
    acf_vals <- tryCatch(pacf(ts)$acf, error = function(e) NA)
    if (is.numeric(acf_vals) && length(acf_vals) > 1 && !all(is.na(acf_vals))) {
      thresh <- 2 / sqrt(N)
      sig_lags <- which(abs(acf_vals) > thresh)
      res_df["PAFC.n", col] <- length(sig_lags)
      res_df["PAFC.max_lag", col] <- if (length(sig_lags)>0) max(sig_lags) else NA
    }
    
    res_Bringmann2013[[toString(subj)]] <- res_df
  }
}


# --- 2.4 Step 3 – TV-AR ---

for (subj in subjects.Bringmann2013) {
  data <- data_Bringmann2013[data_Bringmann2013$subj_id == subj, cols.Bringmann2013]
  for (col in cols.Bringmann2013) {
    ts <- na.omit(data[[col]])
    N <- sum(!is.na(ts)) 
    if (N < MIN_NUM_OBS) next
    
    res_df <- res_Bringmann2013[[toString(subj)]]
    
    tt <- 1:(N-1)
    tv <- tryCatch(gam(ts[2:N] ~ s(tt, by = ts[1:(N-1)], k = 10, bs = "tp")), error = function(e) NULL)
    if (!is.null(tv)) {
      stvar <- summary(tv)
      res_df["TV-AR.EDF", col] <- round(stvar$edf, 3)
      res_df["TV-AR.p", col] <- round(stvar$s.table[4], 4)
    }
    
    res_Bringmann2013[[toString(subj)]] <- res_df
  }
}


# --- 2.5 Step 4 – Change Point Analysis ---

for (subj in subjects.Bringmann2013) {
  data <- data_Bringmann2013[data_Bringmann2013$subj_id == subj, cols.Bringmann2013]
  for (col in cols.Bringmann2013) {
    ts <- na.omit(data[[col]])
    N <- sum(!is.na(ts)) 
    if (N < MIN_NUM_OBS) next
    
    res_df <- res_Bringmann2013[[toString(subj)]]
    
    R <- ceiling(20 / ALPHA_LEVEL)
    cp.out <- tryCatch(e.divisive(matrix(ts), R = R, sig.lvl = ALPHA_LEVEL), error = function(e) NULL)
    if (!is.null(cp.out)) {
      cp.n <- length(which(cp.out$p.values < ALPHA_LEVEL))
      res_df["CP.n", col] <- cp.n
    }
    
    res_Bringmann2013[[toString(subj)]] <- res_df
  }
}


# --- 2.6 Step 5 – Forecast Skill ---

for (subj in subjects.Bringmann2013) {
  data <- data_Bringmann2013[data_Bringmann2013$subj_id == subj, cols.Bringmann2013]
  for (col in cols.Bringmann2013) {
    ts <- na.omit(data[[col]])
    N <- sum(!is.na(ts)) 
    
    # check if ts has minimum number of observations
    if (N < MIN_NUM_OBS) next
    
    # If ts has no variance, the calculation of timelag tau crashes the entire R
    # session, so check for that
    if (var(ts) == 0) next
    
    res_df <- res_Bringmann2013[[toString(subj)]]
    
    tau <- tryCatch(suppressWarnings(
      timeLag(unlist(ts), technique = "ami", selection.method = "first.minimum", lag.max = 10, do.plot = FALSE)
    ), error = function(e) NA)
    
    if (!is.na(tau)) {
      mid <- floor(N / 2)
      lib <- c(1, mid)
      pred <- c(mid + 1, N)
      
      e <- tryCatch({
        emb_out <- suppressWarnings(EmbedDimension(dataFrame = as.data.frame(ts), lib = lib, pred = pred, maxE = 10, Tp = 1, tau = tau, columns = "ts", target = "ts", noTime = TRUE, showPlot = FALSE))
        emb_out$E[which.max(emb_out$rho)]
      }, error = function(e) NA)
      
      if (!is.na(e)) {
        pi <- tryCatch(suppressWarnings(
          PredictInterval(dataFrame = as.data.frame(ts), noTime = TRUE, columns = "ts", target = "ts", 
                          E = e, tau = tau, lib = lib, pred = pred, maxTp = 20, showPlot = FALSE)
        ), error = function(e) NULL)
        
        if (!is.null(pi) && nrow(pi) >= 5) {
          decay <- tryCatch({
            fit <- lm(rho ~ Tp, data = pi[1:5,])
            coef(fit)[2] * 5
          }, error = function(e) NA)
          res_df["pred_decay", col] <- decay
        }
      }
    }
    
    res_Bringmann2013[[toString(subj)]] <- res_df
  }
}


# --- 2.7 Step 6 – Multifractality ---

for (subj in subjects.Bringmann2013) {
  print(subj)
  data <- data_Bringmann2013[data_Bringmann2013$subj_id == subj, cols.Bringmann2013]
  for (col in cols.Bringmann2013) {
    ts <- na.omit(data[[col]])
    N <- sum(!is.na(ts)) 
    if (N < MIN_NUM_OBS) next
    
    res_df <- res_Bringmann2013[[toString(subj)]]
    
    spectrum <- tryCatch(CJspectrum(ts, 1, N, 1), error = function(e) NULL)
    n_surr <- ceiling(pwr.t.test(d = 0.5, sig.level = 0.05, power = 0.80, type = "one.sample")$n)
    surrogates <- tryCatch(IAAFT(as.numeric(ts), n_surr, 1000), error = function(e) NULL)
    
    spectra <- c()
    if (!is.null(surrogates)) {
      for (i in 1:ncol(surrogates)) {
        surr <- surrogates[,i]
        spectrum_i <- tryCatch(CJspectrum(surr, 1, length(surr), 1), error = function(e) NULL)
        spectra <- c(spectra, spectrum_i)
      }
    }
    
    if (!is.null(spectrum)) {
      pval.t <- tryCatch(t.test(spectra, mu = spectrum)$p.value, error = function(e) NA)
      res_df["multifractality", col] <- pval.t
    }
    
    res_Bringmann2013[[toString(subj)]] <- res_df
  }
}


# --- 2.8 Summary Statistics ---

test_vectors.Bringmann2013 <- setNames(vector("list", length(tests)), tests)

for (res in res_Bringmann2013) {
  for (test in tests) {
    values <- as.numeric(res[test, , drop = TRUE])
    test_vectors.Bringmann2013[[test]] <- c(test_vectors.Bringmann2013[[test]], values)
  }
}

prop.sig.Bringmann2013 <- sum(unlist(test_vectors.Bringmann2013[significance_tests]) < ALPHA_LEVEL, na.rm = TRUE) /
  length(unlist(test_vectors.Bringmann2013[significance_tests])) # 0.4358242
prop.sig.Bringmann2013
# --> Approx. 43.5% of all test x participant x emotion time series display 
# complexity marker

mean(test_vectors.Bringmann2013$CP.n, na.rm = TRUE) # 0.6217054
mean(test_vectors.Bringmann2013$pred_decay, na.rm = TRUE) # -0.4609777
mean(test_vectors.Bringmann2013$PAFC.n, na.rm = TRUE) # 1.474178
mean(test_vectors.Bringmann2013$PAFC.max_lag, na.rm = TRUE) # 6.223684


# ==========================================================
# 3 ANALYSIS II – FRIED ET AL. (2021)
# ==========================================================

# --- 3.1 Data Preparation & Descriptives ---

cols.Fried2021 <- c("Relax", "Irritable", "Worry", "Nervous", "Angry")
subjects.Fried2021 <- unique(data_Fried2021$subj_id)

res_Fried2021 <- list()

for (subj in subjects.Fried2021) {
  res_df <- data.frame(matrix(NA,
                              nrow = length(tests),
                              ncol = length(cols.Fried2021),
                              dimnames = list(tests, cols.Fried2021)))
  res_Fried2021[[toString(subj)]] <- res_df
}

# Descriptives
ts_lengths.Fried2021 <- c()
for (subj in subjects.Fried2021) {
  data <- data_Fried2021[data_Fried2021$subj_id == subj, cols.Fried2021]
  for (col in cols.Fried2021) {
    ts <- na.omit(data[[col]])
    N <- length(ts)
    
    # append to ts_lengths
    ts_lengths.Fried2021 <- c(ts_lengths.Fried2021, N)
  }
}

mean(ts_lengths.Fried2021) # 49.95696


start <- Sys.time()

# --- 3.2 Step 1 – Significance Tests ---

for (subj in subjects.Fried2021) {
  data <- data_Fried2021[data_Fried2021$subj_id == subj, cols.Fried2021]
  for (col in cols.Fried2021) {
    ts <- na.omit(data[[col]])
    N <- sum(!is.na(ts)) 
    if (N < MIN_NUM_OBS) next
    
    res_df <- res_Fried2021[[toString(subj)]]
    
    bartels_p <- tryCatch(bartels.rank.test(ts, alternative = "two.sided")$p.value, error = function(e) NA)
    res_df["Bartels", col] <- if (is.null(bartels_p) || length(bartels_p) == 0) NA else bartels_p
    
    kpss_level_p <- tryCatch(kpss.test(ts, null = "Level", lshort = TRUE)$p.value, error = function(e) NA)
    res_df["KPSS.level", col] <- if (is.null(kpss_level_p) || length(kpss_level_p) == 0) NA else kpss_level_p
    
    kpss_trend_p <- tryCatch(kpss.test(ts, null = "Trend", lshort = TRUE)$p.value, error = function(e) NA)
    res_df["KPSS.trend", col] <- if (is.null(kpss_trend_p) || length(kpss_trend_p) == 0) NA else kpss_trend_p
    
    keenan_p <- tryCatch(keenanTest(ts)$p.value, error = function(e) NA)
    res_df["Keenan", col] <- if (is.null(keenan_p) || length(keenan_p) == 0) NA else keenan_p
    
    tsay_p <- tryCatch(tsayTest(ts)$p.value, error = function(e) NA)
    res_df["Tsay", col] <- if (is.null(tsay_p) || length(tsay_p) == 0) NA else tsay_p
    
    res_Fried2021[[toString(subj)]] <- res_df
  }
}


# --- 3.3 Step 2 – ACF Analysis ---

for (subj in subjects.Fried2021) {
  data <- data_Fried2021[data_Fried2021$subj_id == subj, cols.Fried2021]
  for (col in cols.Fried2021) {
    ts <- na.omit(data[[col]])
    N <- sum(!is.na(ts)) 
    if (N < MIN_NUM_OBS) next
    
    res_df <- res_Fried2021[[toString(subj)]]
    
    acf_vals <- tryCatch(pacf(ts)$acf, error = function(e) NA)
    if (is.numeric(acf_vals) && length(acf_vals) > 1 && !all(is.na(acf_vals))) {
      thresh <- 2 / sqrt(N)
      sig_lags <- which(abs(acf_vals) > thresh)
      res_df["PAFC.n", col] <- length(sig_lags)
      res_df["PAFC.max_lag", col] <- if (length(sig_lags)>0) max(sig_lags) else NA
    }
    
    res_Fried2021[[toString(subj)]] <- res_df
  }
}


# --- 3.4 Step 3 – TV-AR ---

for (subj in subjects.Fried2021) {
  data <- data_Fried2021[data_Fried2021$subj_id == subj, cols.Fried2021]
  for (col in cols.Fried2021) {
    ts <- na.omit(data[[col]])
    N <- sum(!is.na(ts)) 
    if (N < MIN_NUM_OBS) next
    
    res_df <- res_Fried2021[[toString(subj)]]
    
    tt <- 1:(N-1)
    tv <- tryCatch(gam(ts[2:N] ~ s(tt, by = ts[1:(N-1)], k = 10, bs = "tp")), error = function(e) NULL)
    if (!is.null(tv)) {
      stvar <- summary(tv)
      res_df["TV-AR.EDF", col] <- round(stvar$edf, 3)
      res_df["TV-AR.p", col] <- round(stvar$s.table[4], 4)
    }
    
    res_Fried2021[[toString(subj)]] <- res_df
  }
}


# --- 3.5 Step 4 – Change Point Analysis ---

for (subj in subjects.Fried2021) {
  data <- data_Fried2021[data_Fried2021$subj_id == subj, cols.Fried2021]
  for (col in cols.Fried2021) {
    ts <- na.omit(data[[col]])
    N <- sum(!is.na(ts)) 
    if (N < MIN_NUM_OBS) next
    
    res_df <- res_Fried2021[[toString(subj)]]
    
    R <- ceiling(20 / ALPHA_LEVEL)
    cp.out <- tryCatch(e.divisive(matrix(ts), R = R, sig.lvl = ALPHA_LEVEL), error = function(e) NULL)
    if (!is.null(cp.out)) {
      cp.n <- length(which(cp.out$p.values < ALPHA_LEVEL))
      res_df["CP.n", col] <- cp.n
    }
    
    res_Fried2021[[toString(subj)]] <- res_df
  }
}


# --- 3.6 Step 5 – Forecast Skill ---

for (subj in subjects.Fried2021) {
  data <- data_Fried2021[data_Fried2021$subj_id == subj, cols.Fried2021]
  for (col in cols.Fried2021) {
    ts <- na.omit(data[[col]])
    N <- sum(!is.na(ts)) 
    
    # check if ts has minimum number of observations
    if (N < MIN_NUM_OBS) next
    
    # If ts has no variance, the calculation of timelag tau crashes the entire R
    # session, so check for that
    if (var(ts) == 0) next
    
    res_df <- res_Fried2021[[toString(subj)]]
    
    tau <- tryCatch(suppressWarnings(
      timeLag(unlist(ts), technique = "ami", selection.method = "first.minimum", lag.max = 10, do.plot = FALSE)
    ), error = function(e) NA)
    
    if (!is.na(tau)) {
      mid <- floor(N / 2)
      lib <- c(1, mid)
      pred <- c(mid + 1, N)
      
      e <- tryCatch({
        emb_out <- suppressWarnings(EmbedDimension(dataFrame = as.data.frame(ts), lib = lib, pred = pred, maxE = 10, Tp = 1, tau = tau, columns = "ts", target = "ts", noTime = TRUE, showPlot = FALSE))
        emb_out$E[which.max(emb_out$rho)]
      }, error = function(e) NA)
      
      if (!is.na(e)) {
        pi <- tryCatch(suppressWarnings(
          PredictInterval(dataFrame = as.data.frame(ts), noTime = TRUE, columns = "ts", target = "ts", 
                          E = e, tau = tau, lib = lib, pred = pred, maxTp = 20, showPlot = FALSE)
        ), error = function(e) NULL)
        
        if (!is.null(pi) && nrow(pi) >= 5) {
          decay <- tryCatch({
            fit <- lm(rho ~ Tp, data = pi[1:5,])
            coef(fit)[2] * 5
          }, error = function(e) NA)
          res_df["pred_decay", col] <- decay
        }
      }
    }
    
    res_Fried2021[[toString(subj)]] <- res_df
  }
}


# --- 3.7 Step 6 – Multifractality ---

for (subj in subjects.Fried2021) {
  print(subj)
  data <- data_Fried2021[data_Fried2021$subj_id == subj, cols.Fried2021]
  for (col in cols.Fried2021) {
    ts <- na.omit(data[[col]])
    N <- sum(!is.na(ts)) 
    if (N < MIN_NUM_OBS) next
    
    res_df <- res_Fried2021[[toString(subj)]]
    
    spectrum <- tryCatch(CJspectrum(ts, 1, N, 1), error = function(e) NULL)
    n_surr <- ceiling(pwr.t.test(d = 0.5, sig.level = 0.05, power = 0.80, type = "one.sample")$n)
    surrogates <- tryCatch(IAAFT(as.numeric(ts), n_surr, 1000), error = function(e) NULL)
    
    spectra <- c()
    if (!is.null(surrogates)) {
      for (i in 1:ncol(surrogates)) {
        surr <- surrogates[,i]
        spectrum_i <- tryCatch(CJspectrum(surr, 1, length(surr), 1), error = function(e) NULL)
        spectra <- c(spectra, spectrum_i)
      }
    }
    
    if (!is.null(spectrum)) {
      pval.t <- tryCatch(t.test(spectra, mu = spectrum)$p.value, error = function(e) NA)
      res_df["multifractality", col] <- pval.t
    }
    
    res_Fried2021[[toString(subj)]] <- res_df
  }
}


# --- 3.8 Summary Statistics ---

test_vectors.Fried2021 <- setNames(vector("list", length(tests)), tests)

for (res in res_Fried2021) {
  for (test in tests) {
    values <- as.numeric(res[test, , drop = TRUE])
    test_vectors.Fried2021[[test]] <- c(test_vectors.Fried2021[[test]], values)
  }
}

prop.sig.Fried2021 <- sum(unlist(test_vectors.Fried2021[significance_tests]) < ALPHA_LEVEL, na.rm = TRUE) /
  length(unlist(test_vectors.Fried2021[significance_tests])) # 0.277396
prop.sig.Fried2021
# ---> approx. 27.7% of all test x participant x emotion time series display 
# complexity marker

mean(test_vectors.Fried2021$CP.n, na.rm = TRUE) # 0
mean(test_vectors.Fried2021$pred_decay, na.rm = TRUE) # -0.6590549
mean(test_vectors.Fried2021$PAFC.n, na.rm = TRUE) # 0.6910569
mean(test_vectors.Fried2021$PAFC.max_lag, na.rm = TRUE) # 5.13089


# ==========================================================
# 4 BETWEEN-DATASET COMPARISONS
# ==========================================================

# --- 4.1 Proportion Test ---

# test for significant difference between proportion of significant tests
sig_tests <- c(
  "Bartels",
  "TV-AR.p",
  "KPSS.level",
  "KPSS.trend",
  "Keenan",
  "Tsay",
  "multifractality")

n.sig_test <- length(sig_tests)

n_tests.Bringmann2013 <- n.sig_test * length(subjects.Bringmann2013) * length(cols.Bringmann2013)
n_tests.Fried2021 <- n.sig_test * length(subjects.Fried2021) * length(cols.Fried2021)

x.Bringmann2013 <- round(prop.sig.Bringmann2013 * n_tests.Bringmann2013)
x.Fried2021 <- round(prop.sig.Fried2021 * n_tests.Fried2021)

prop.test(x =  c(x.Bringmann2013, x.Fried2021), n = c(n_tests.Bringmann2013, n_tests.Fried2021), 
          alternative = "greater", correct = FALSE)

# 2-sample test for equality of proportions without continuity correction
# 
# data:  c(x.Bringmann2013, x.Fried2021) out of c(n_tests.Bringmann2013, n_tests.Fried2021)
# X-squared = 184, df = 1, p-value < 2.2e-16
# alternative hypothesis: greater
# 95 percent confidence interval:
#   0.1399256 1.0000000
# sample estimates:
#   prop 1    prop 2 
# 0.4358242 0.2773960 

# --> Bringmann et al. (2013) has more significant tests.

# --- 4.2 Change Point Comparison ---

t.test(test_vectors.Bringmann2013$CP.n,
       test_vectors.Fried2021$CP.n,
       alternative = "greater",
       var.equal = FALSE)

# Welch Two Sample t-test
# 
# data:  test_vectors.Bringmann2013$CP.n and test_vectors.Fried2021$CP.n
# t = 25.256, df = 644, p-value < 2.2e-16
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   0.5811576       Inf
# sample estimates:
#   mean of x mean of y 
# 0.6217054 0.0000000 

# --> Bringmann et al. (2013) has more change points.


# --- 4.3 Prediction Decay Comparison ---

t.test(abs(test_vectors.Fried2021$pred_decay),
       abs(test_vectors.Bringmann2013$pred_decay),
       alternative = "greater",
       var.equal = FALSE)

# Welch Two Sample t-test
# 
# data:  abs(test_vectors.Fried2021$pred_decay) and abs(test_vectors.Bringmann2013$pred_decay)
# t = 4.8987, df = 481.79, p-value = 6.599e-07
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   0.1109932       Inf
# sample estimates:
#   mean of x mean of y 
# 0.6659069 0.4986419 

# --> Fried et al. (2021) has higher abs. prediction decay.


# --- 4.4 Partial Autocorrelation Comparison ---

t.test(test_vectors.Bringmann2013$PAFC.n,
       test_vectors.Fried2021$PAFC.n,
       alternative = "greater",
       var.equal = FALSE)

# Welch Two Sample t-test
# 
# data:  test_vectors.Bringmann2013$PAFC.n and test_vectors.Fried2021$PAFC.n
# t = 13.608, df = 933.71, p-value < 2.2e-16
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   0.6883687       Inf
# sample estimates:
#   mean of x mean of y 
# 1.4741784 0.6910569 

# --> Bringmann et al. (2013) has more significant partial autocorrelations.


# --- 4.5 Interpretation ---

# ---> Yes, time series in Bringmann et al. (2013) displays more complexity.
# ---> I argue this is because time series are longer in this dataset.


# ==========================================================
# 5 SAVE WORKSPACE
# ==========================================================

save.image("COMPLEXITY_TESTS.RData")
