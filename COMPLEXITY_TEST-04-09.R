rm(list = ls())
setwd('C:/Users/bbe0870/Desktop/PhD1/complexity_in_emotion_ESM_data')

library(nonlinearTseries)
library(randtests)
library(tseries)
library(mice)
library(pwr)
library(stats)
library(mgcv)
library(ecp)
library(rEDM)

set.seed(4168)

# GLOBAL VARS
MIN_NUM_OBS <- 20
ALPHA_LEVEL <- 0.05

# Import helper functions for multifractality 
# (see https://osf.io/nm9b4/, Kelty-Stephen et al., 2023)
source('CJspectrum.R')
source('IAAFT.R')

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

# Datasets, preprocessed by Olthof et al.
data_Bringmann2013 <- readRDS("C:/Users/bbe0870/Desktop/R/Complexity in Emotion Data(1)/DataClean/Bringmann2013/data_Bringmann2013.RDS")
data_Bringmann2016 <- readRDS("C:/Users/bbe0870/Desktop/R/Complexity in Emotion Data(1)/DataClean/Bringmann2016/data_Bringmann2016.RDS")
data_Fisher2017 <- readRDS("C:/Users/bbe0870/Desktop/R/Complexity in Emotion Data(1)/DataClean/Fisher2017/data_Fisher2017.RDS")
data_Fried2022 <- readRDS("C:/Users/bbe0870/Desktop/R/Complexity in Emotion Data(1)/DataClean/Fried2021/data_Fried2021.RDS")
data_Rowland2020 <- readRDS("C:/Users/bbe0870/Desktop/R/Complexity in Emotion Data(1)/DataClean/Rowland2020/data_Rowland2020.RDS")
data_Wright2017 <- readRDS("C:/Users/bbe0870/Desktop/R/Complexity in Emotion Data(1)/DataClean/Wright2017/data_Wright2017.RDS")


# Bringmann et al. 2013

cols.Bringmann_2013 <- c("Worried", "Sad", "Anxious", "Relaxed", "Excited")
subjects.Bringmann_2013 <- unique(data_Bringmann2013$subj_id)

res_Bringmann_2013 <- list()

for (subj in subjects.Bringmann_2013) {
  res_df <- data.frame(matrix(NA,
                              nrow = length(tests),
                              ncol = length(cols.Bringmann_2013),
                              dimnames = list(tests, cols.Bringmann_2013)))
  res_Bringmann_2013[[toString(subj)]] <- res_df
}

# Descriptives
ts_lengths <- c()
for (subj in subjects.Bringmann_2013) {
  data <- data_Bringmann2013[data_Bringmann2013$subj_id == subj, cols.Bringmann_2013]
  for (col in cols.Bringmann_2013) {
    ts <- na.omit(data[[col]])
    N <- length(ts)
    
    # append to ts_lengths
    ts_lengths <- c(ts_lengths, N)
  }
}

mean(ts_lengths) #90.17538


start <- Sys.time()

# Step 1 – Fast Statistical Tests
for (subj in subjects.Bringmann_2013) {
  data <- data_Bringmann2013[data_Bringmann2013$subj_id == subj, cols.Bringmann_2013]
  for (col in cols.Bringmann_2013) {
    ts <- na.omit(data[[col]])
    N <- sum(!is.na(ts)) 
    if (N < MIN_NUM_OBS) next
    
    res_df <- res_Bringmann_2013[[toString(subj)]]
    
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
    
    res_Bringmann_2013[[toString(subj)]] <- res_df
  }
}

# Step 2 – ACF Analysis
for (subj in subjects.Bringmann_2013) {
  data <- data_Bringmann2013[data_Bringmann2013$subj_id == subj, cols.Bringmann_2013]
  for (col in cols.Bringmann_2013) {
    ts <- na.omit(data[[col]])
    N <- sum(!is.na(ts)) 
    if (N < MIN_NUM_OBS) next
    
    res_df <- res_Bringmann_2013[[toString(subj)]]
    
    acf_vals <- tryCatch(pacf(ts)$acf, error = function(e) NA)
    if (is.numeric(acf_vals) && length(acf_vals) > 1 && !all(is.na(acf_vals))) {
      thresh <- 2 / sqrt(N)
      sig_lags <- which(abs(acf_vals) > thresh)
      res_df["PAFC.n", col] <- length(sig_lags)
      res_df["PAFC.max_lag", col] <- if (length(sig_lags)>0) max(sig_lags) else NA
    }
    
    res_Bringmann_2013[[toString(subj)]] <- res_df
  }
}

# Step 3 – TV-AR
for (subj in subjects.Bringmann_2013) {
  data <- data_Bringmann2013[data_Bringmann2013$subj_id == subj, cols.Bringmann_2013]
  for (col in cols.Bringmann_2013) {
    ts <- na.omit(data[[col]])
    N <- sum(!is.na(ts)) 
    if (N < MIN_NUM_OBS) next
    
    res_df <- res_Bringmann_2013[[toString(subj)]]
    
    tt <- 1:(N-1)
    tv <- tryCatch(gam(ts[2:N] ~ s(tt, by = ts[1:(N-1)], k = 10, bs = "tp")), error = function(e) NULL)
    if (!is.null(tv)) {
      stvar <- summary(tv)
      res_df["TV-AR.EDF", col] <- round(stvar$edf, 3)
      res_df["TV-AR.p", col] <- round(stvar$s.table[4], 4)
    }
    
    res_Bringmann_2013[[toString(subj)]] <- res_df
  }
}

# Step 4 – Change Point
for (subj in subjects.Bringmann_2013) {
  data <- data_Bringmann2013[data_Bringmann2013$subj_id == subj, cols.Bringmann_2013]
  for (col in cols.Bringmann_2013) {
    ts <- na.omit(data[[col]])
    N <- sum(!is.na(ts)) 
    if (N < MIN_NUM_OBS) next
    
    res_df <- res_Bringmann_2013[[toString(subj)]]
    
    R <- ceiling(20 / ALPHA_LEVEL)
    cp.out <- tryCatch(e.divisive(matrix(ts), R = R, sig.lvl = ALPHA_LEVEL), error = function(e) NULL)
    if (!is.null(cp.out)) {
      cp.n <- length(which(cp.out$p.values < ALPHA_LEVEL))
      res_df["CP.n", col] <- cp.n
    }
    
    res_Bringmann_2013[[toString(subj)]] <- res_df
  }
}

# Step 5 – Forecast Skill
for (subj in subjects.Bringmann_2013) {
  data <- data_Bringmann2013[data_Bringmann2013$subj_id == subj, cols.Bringmann_2013]
  for (col in cols.Bringmann_2013) {
    ts <- na.omit(data[[col]])
    N <- sum(!is.na(ts)) 
    
    # check if ts has minimum number of observations
    if (N < MIN_NUM_OBS) next
    
    # If ts has no variance, the calculation of timelag tau crashes the entire R
    # session, so check for that
    if (var(ts) == 0) next
    
    res_df <- res_Bringmann_2013[[toString(subj)]]
    
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
    
    res_Bringmann_2013[[toString(subj)]] <- res_df
  }
}

# Step 6 – Multifractality
for (subj in subjects.Bringmann_2013) {
  print(subj)
  data <- data_Bringmann2013[data_Bringmann2013$subj_id == subj, cols.Bringmann_2013]
  for (col in cols.Bringmann_2013) {
    ts <- na.omit(data[[col]])
    N <- sum(!is.na(ts)) 
    if (N < MIN_NUM_OBS) next
    
    res_df <- res_Bringmann_2013[[toString(subj)]]
    
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
    
    res_Bringmann_2013[[toString(subj)]] <- res_df
  }
}

# get overview over all tests
test_vectors.Bringmann2013 <- setNames(vector("list", length(tests)), tests)

for (res in res_Bringmann_2013) {
  for (test in tests) {
    values <- as.numeric(res[test, , drop = TRUE])
    test_vectors.Bringmann2013[[test]] <- c(test_vectors.Bringmann2013[[test]], values)
  }
}

prop.sig.Bringmann2013 <- sum(unlist(test_vectors.Bringmann2013[significance_tests]) < ALPHA_LEVEL, na.rm = TRUE) /
  length(unlist(test_vectors.Bringmann2013[significance_tests])) # 0.4342857
prop.sig.Bringmann2013
# approx. 43.4 of all test x participant x emotion time series display 
# complexity marker

mean(test_vectors.Bringmann2013$CP.n, na.rm = TRUE) # 0.6265625
mean(test_vectors.Bringmann2013$pred_decay, na.rm = TRUE) # -0.4609777
mean(test_vectors.Bringmann2013$PAFC.n, na.rm = TRUE) # 1.47874
mean(test_vectors.Bringmann2013$PAFC.max_lag, na.rm = TRUE) # 6.215094

# write.csv(res_Bringmann_2013, file = "res_Bringmann_2013.csv")


# Fried et al. 2022
cols.Fried_2022 <- c("Relax", "Irritable", "Worry", "Nervous", "Angry")
subjects.Fried_2022 <- unique(data_Fried2022$subj_id)

res_Fried_2022 <- list()

for (subj in subjects.Fried_2022) {
  res_df <- data.frame(matrix(NA,
                              nrow = length(tests),
                              ncol = length(cols.Fried_2022),
                              dimnames = list(tests, cols.Fried_2022)))
  res_Fried_2022[[toString(subj)]] <- res_df
}

# Descriptives
ts_lengths.Fried_2022 <- c()
for (subj in subjects.Fried_2022) {
  data <- data_Fried2022[data_Fried2022$subj_id == subj, cols.Fried_2022]
  for (col in cols.Fried_2022) {
    ts <- na.omit(data[[col]])
    N <- length(ts)
    
    # append to ts_lengths
    ts_lengths.Fried_2022 <- c(ts_lengths.Fried_2022, N)
  }
}

mean(ts_lengths.Fried_2022) # 49.95696


start <- Sys.time()

# Step 1 – Fast Statistical Tests
for (subj in subjects.Fried_2022) {
  data <- data_Fried2022[data_Fried2022$subj_id == subj, cols.Fried_2022]
  for (col in cols.Fried_2022) {
    ts <- na.omit(data[[col]])
    N <- sum(!is.na(ts)) 
    if (N < MIN_NUM_OBS) next
    
    res_df <- res_Fried_2022[[toString(subj)]]
    
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
    
    res_Fried_2022[[toString(subj)]] <- res_df
  }
}

# Step 2 – ACF Analysis
for (subj in subjects.Fried_2022) {
  data <- data_Fried2022[data_Fried2022$subj_id == subj, cols.Fried_2022]
  for (col in cols.Fried_2022) {
    ts <- na.omit(data[[col]])
    N <- sum(!is.na(ts)) 
    if (N < MIN_NUM_OBS) next
    
    res_df <- res_Fried_2022[[toString(subj)]]
    
    acf_vals <- tryCatch(pacf(ts)$acf, error = function(e) NA)
    if (is.numeric(acf_vals) && length(acf_vals) > 1 && !all(is.na(acf_vals))) {
      thresh <- 2 / sqrt(N)
      sig_lags <- which(abs(acf_vals) > thresh)
      res_df["PAFC.n", col] <- length(sig_lags)
      res_df["PAFC.max_lag", col] <- if (length(sig_lags)>0) max(sig_lags) else NA
    }
    
    res_Fried_2022[[toString(subj)]] <- res_df
  }
}

# Step 3 – TV-AR
for (subj in subjects.Fried_2022) {
  data <- data_Fried2022[data_Fried2022$subj_id == subj, cols.Fried_2022]
  for (col in cols.Fried_2022) {
    ts <- na.omit(data[[col]])
    N <- sum(!is.na(ts)) 
    if (N < MIN_NUM_OBS) next
    
    res_df <- res_Fried_2022[[toString(subj)]]
    
    tt <- 1:(N-1)
    tv <- tryCatch(gam(ts[2:N] ~ s(tt, by = ts[1:(N-1)], k = 10, bs = "tp")), error = function(e) NULL)
    if (!is.null(tv)) {
      stvar <- summary(tv)
      res_df["TV-AR.EDF", col] <- round(stvar$edf, 3)
      res_df["TV-AR.p", col] <- round(stvar$s.table[4], 4)
    }
    
    res_Fried_2022[[toString(subj)]] <- res_df
  }
}

# Step 4 – Change Point
for (subj in subjects.Fried_2022) {
  data <- data_Fried2022[data_Fried2022$subj_id == subj, cols.Fried_2022]
  for (col in cols.Fried_2022) {
    ts <- na.omit(data[[col]])
    N <- sum(!is.na(ts)) 
    if (N < MIN_NUM_OBS) next
    
    res_df <- res_Fried_2022[[toString(subj)]]
    
    R <- ceiling(20 / ALPHA_LEVEL)
    cp.out <- tryCatch(e.divisive(matrix(ts), R = R, sig.lvl = ALPHA_LEVEL), error = function(e) NULL)
    if (!is.null(cp.out)) {
      cp.n <- length(which(cp.out$p.values < ALPHA_LEVEL))
      res_df["CP.n", col] <- cp.n
    }
    
    res_Fried_2022[[toString(subj)]] <- res_df
  }
}

# Step 5 – Forecast Skill
for (subj in subjects.Fried_2022) {
  data <- data_Fried2022[data_Fried2022$subj_id == subj, cols.Fried_2022]
  for (col in cols.Fried_2022) {
    ts <- na.omit(data[[col]])
    N <- sum(!is.na(ts)) 
    
    # check if ts has minimum number of observations
    if (N < MIN_NUM_OBS) next
    
    # If ts has no variance, the calculation of timelag tau crashes the entire R
    # session, so check for that
    if (var(ts) == 0) next
    
    res_df <- res_Fried_2022[[toString(subj)]]
    
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
    
    res_Fried_2022[[toString(subj)]] <- res_df
  }
}

# Step 6 – Multifractality
for (subj in subjects.Fried_2022) {
  print(subj)
  data <- data_Fried2022[data_Fried2022$subj_id == subj, cols.Fried_2022]
  for (col in cols.Fried_2022) {
    ts <- na.omit(data[[col]])
    N <- sum(!is.na(ts)) 
    if (N < MIN_NUM_OBS) next
    
    res_df <- res_Fried_2022[[toString(subj)]]
    
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
    
    res_Fried_2022[[toString(subj)]] <- res_df
  }
}

# get overview over all tests
test_vectors.Fried2022 <- setNames(vector("list", length(tests)), tests)

for (res in res_Fried_2022) {
  for (test in tests) {
    values <- as.numeric(res[test, , drop = TRUE])
    test_vectors.Fried2022[[test]] <- c(test_vectors.Fried2022[[test]], values)
  }
}

prop.sig.Fried2022 <- sum(unlist(test_vectors.Fried2022[significance_tests]) < ALPHA_LEVEL, na.rm = TRUE) /
  length(unlist(test_vectors.Fried2022[significance_tests])) # 0.2802893
prop.sig.Fried2022
# approx. 28% of all test x participant x emotion time series display 
# complexity marker

mean(test_vectors.Fried2022$CP.n, na.rm = TRUE) # 0
mean(test_vectors.Fried2022$pred_decay, na.rm = TRUE) # -0.6590549
mean(test_vectors.Fried2022$PAFC.n, na.rm = TRUE) # 0.6910569
mean(test_vectors.Fried2022$PAFC.max_lag, na.rm = TRUE) # 5.13089



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

n_tests.Bringmann2013 <- n.sig_test * length(subjects.Bringmann_2013) * length(cols.Bringmann_2013)
n_tests.Fried2022 <- n.sig_test * length(subjects.Fried_2022) * length(cols.Fried_2022)

x.Bringmann2013 <- round(prop.sig.Bringmann2013 * n_tests.Bringmann2013)
x.Fried2022 <- round(prop.sig.Fried2022 * n_tests.Fried2022)

prop.test(x = c(test_vectors.Bringmann2013$CP.n, mean(test_vectors.Fried2022$CP.n, na.rm = TRUE)), n = c(n_tests.Bringmann2013, n_tests.Fried2022), 
          alternative = "greater", correct = FALSE)

# 2-sample test for equality of proportions without continuity correction
# 
# data:  c(x.Bringmann2013, x.Fried2022) out of c(n_tests.Bringmann2013, n_tests.Fried2022)
# X-squared = 173.82, df = 1, p-value < 2.2e-16
# alternative hypothesis: greater
# 95 percent confidence interval:
#   0.1354632 1.0000000
# sample estimates:
#   prop 1    prop 2 
# 0.4342857 0.2802893 

# Change point analysis
t.test(test_vectors.Bringmann2013$CP.n,
       test_vectors.Fried2022$CP.n,
       alternative = "greater",
       var.equal = FALSE)
# 
# data:  test_vectors.Bringmann2013$CP.n and test_vectors.Fried2022$CP.n
# t = 25.354, df = 639, p-value < 2.2e-16
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   0.5858556       Inf
# sample estimates:
#   mean of x mean of y 
# 0.6265625 0.0000000 

# --> Bringmann et al. has more change points

# prediction decay on absolute values
t.test(abs(test_vectors.Fried2022$pred_decay),
       abs(test_vectors.Bringmann2013$pred_decay),
       alternative = "greater",
       var.equal = FALSE)

# data:  abs(test_vectors.Fried2022$pred_decay) and abs(test_vectors.Bringmann2013$pred_decay)
# t = 4.8987, df = 481.79, p-value = 6.599e-07
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   0.1109932       Inf
# sample estimates:
#   mean of x mean of y 
# 0.6659069 0.4986419 

# --> Bringmann et al. has actually lower prediction decay

# number of significant partial autocorrelations
t.test(test_vectors.Bringmann2013$PAFC.n,
       test_vectors.Fried2022$PAFC.n,
       alternative = "greater",
       var.equal = FALSE)

# data:  test_vectors.Bringmann2013$PAFC.n and test_vectors.Fried2022$PAFC.n
# t = 13.669, df = 932.87, p-value < 2.2e-16
# alternative hypothesis: true difference in means is greater than 0
# 95 percent confidence interval:
#   0.6928034       Inf
# sample estimates:
#   mean of x mean of y 
# 1.4787402 0.6910569 


# ---> Yes, time series in Bringmann et al. (2013) displays more complexity.
# ---> I argue this is because time series are longer in this dataset.

save.image("C:/Users/bbe0870/Desktop/R/Complexity in Emotion Data(1)/20-05-2025_20obs.RData")
