library(MASS)

set.seed(123)

library(MASS)
library(nonlinearTseries)
library(randtests)
library(tseries)
library(mice)
library(pwr)
library(stats)
library(mgcv)
library(ecp)
library(rEDM)


### 1 --- SETTINGS --- ###

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

ALPHA_LEVEL <- 0.05

# coef matrix
A <- matrix(c(
  0.4,  0.1, -0.05, -0.05, -0.05,
  0.1,  0.4, -0.05, -0.05, -0.05,
  -0.05, -0.05,  0.4,  0.1,  0.1,
  -0.05, -0.05,  0.1,  0.4,  0.1,
  -0.05, -0.05,  0.1,  0.1,  0.4
), 5, 5, byrow = TRUE)


Sigma <- diag(5)

T_list <- c(25, 50, 75, 100, 150, 200)


### 2 --- MIN-MAX SCALING FUNCTION --- ###

min_max_to_1_7 <- function(x) {
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  x_scaled <- 1 + (x - x_min) * (6 / (x_max - x_min))
  x_round <- pmin(7, pmax(1, round(x_scaled)))
  return(as.integer(x_round))
}

### 3 --- MAIN LOOP: SIMULATION FOR MULTIPLE T --- ###

result_list <- list()

for (T in T_list) {
  
  # raw VAR data
  Y <- matrix(0, T, 5)
  Y[1, ] <- rnorm(5)
  
  for (t in 2:T) {
    eps <- mvrnorm(1, rep(0, 5), Sigma)
    Y[t, ] <- A %*% Y[t - 1, ] + eps
  }
  
  ESM_raw <- as.data.frame(Y)
  colnames(ESM_raw) <- paste0("V", 1:5)
  
  # Likert-converted
  ESM_lik <- as.data.frame(lapply(ESM_raw, min_max_to_1_7))
  colnames(ESM_lik) <- paste0("V", 1:5, "_lik")
  
  # store in list
  result_list[[paste0("T", T)]] <- list(
    raw = ESM_raw,
    likert = ESM_lik
  )
}


### 4 --- APPLY TESTS --- ###

# Progress bars use base R only
n_blocks <- length(result_list)
pb_outer <- txtProgressBar(min = 0, max = n_blocks, style = 3)

test_results <- list()

for (i in seq_along(result_list)) {
  
  # name of the T block (e.g. "T50")
  T_name <- names(result_list)[i]
  
  # Likert data for this T
  data <- result_list[[i]]$likert
  cols <- colnames(data)
  
  # initialize empty results dataframe
  res_df <- matrix(NA, nrow = length(tests), ncol = length(cols),
                   dimnames = list(
                     tests,
                     cols
                   ))
  
  # inner progress bar: each column
  pb_inner <- txtProgressBar(min = 0, max = length(cols), style = 3)
  
  for (j in seq_along(cols)) {
    col <- cols[j]
    ts <- data[[col]]
    N <- length(ts)
    
    # Bartels rank test
    bartels_p <- tryCatch(bartels.rank.test(ts, alternative = "two.sided")$p.value, error = function(e) NA)
    res_df["Bartels", col] <- if (is.null(bartels_p) || length(bartels_p) == 0) NA else bartels_p
    
    # KPSS test level stationarity
    kpss_level_p <- tryCatch(kpss.test(ts, null = "Level", lshort = TRUE)$p.value, error = function(e) NA)
    res_df["KPSS.level", col] <- if (is.null(kpss_level_p) || length(kpss_level_p) == 0) NA else kpss_level_p
    
    # KPSS test trend
    kpss_trend_p <- tryCatch(kpss.test(ts, null = "Trend", lshort = TRUE)$p.value, error = function(e) NA)
    res_df["KPSS.trend", col] <- if (is.null(kpss_trend_p) || length(kpss_trend_p) == 0) NA else kpss_trend_p
    
    # Keenan's nonlinearity test
    keenan_p <- tryCatch(keenanTest(ts)$p.value, error = function(e) NA)
    res_df["Keenan", col] <- if (is.null(keenan_p) || length(keenan_p) == 0) NA else keenan_p
    
    # Tsay's nonlinearity test
    tsay_p <- tryCatch(tsayTest(ts)$p.value, error = function(e) NA)
    res_df["Tsay", col] <- if (is.null(tsay_p) || length(tsay_p) == 0) NA else tsay_p
    
    # Autocorrelation function analysis
    acf_vals <- tryCatch(pacf(ts, plot = FALSE)$acf, error = function(e) NA)
    if (is.numeric(acf_vals) && length(acf_vals) > 1 && !all(is.na(acf_vals))) {
      thresh <- 2 / sqrt(N)
      sig_lags <- which(abs(acf_vals) > thresh)
      res_df["PAFC.n", col] <- length(sig_lags)
      res_df["PAFC.max_lag", col] <- if (length(sig_lags)>0) max(sig_lags) else NA
    }
    
    # Time-varying autoregressive model
    tt <- 1:(N-1)
    tv <- tryCatch(gam(ts[2:N] ~ s(tt, by = ts[1:(N-1)], k = 10, bs = "tp")), error = function(e) NULL)
    if (!is.null(tv)) {
      stvar <- summary(tv)
      res_df["TV-AR.EDF", col] <- round(stvar$edf, 3)
      res_df["TV-AR.p", col] <- round(stvar$s.table[4], 4)
    }
    
    # Change point analysis
    R <- ceiling(20 / ALPHA_LEVEL)
    cp.out <- tryCatch(e.divisive(matrix(ts), R = R, sig.lvl = ALPHA_LEVEL), error = function(e) NULL)
    if (!is.null(cp.out)) {
      cp.n <- length(which(cp.out$p.values < ALPHA_LEVEL))
      res_df["CP.n", col] <- cp.n
    }
    
    # Forecast skill
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
    
    # Multifractality
    spectrum <- tryCatch(CJspectrum(ts, 1, N, 1), error = function(e) NULL)
    n_surr <- ceiling(pwr.t.test(d = 0.5, sig.level = 0.05, power = 0.80, type = "one.sample")$n)
    surrogates <- tryCatch(IAAFT(as.numeric(ts), n_surr, 1000), error = function(e) NULL)
    
    spectra <- c()
    if (!is.null(surrogates)) {
      for (i_surr in 1:ncol(surrogates)) {
        surr <- surrogates[,i_surr]
        spectrum_i <- tryCatch(CJspectrum(surr, 1, length(surr), 1), error = function(e) NULL)
        spectra <- c(spectra, spectrum_i)
      }
    }
    
    if (!is.null(spectrum)) {
      pval.t <- tryCatch(t.test(spectra, mu = spectrum)$p.value, error = function(e) NA)
      res_df["multifractality", col] <- pval.t
    }
    
    # update inner progress bar
    setTxtProgressBar(pb_inner, j)
  }
  
  close(pb_inner)
  
  # store results for this T
  test_results[[T_name]] <- as.data.frame(res_df)
  
  # update outer progress bar
  setTxtProgressBar(pb_outer, i)
}

close(pb_outer)



### 5 --- FILTER TESTS --- ###

filtered_results <- lapply(test_results, function(df) {
  df[rownames(df) %in% significance_tests, , drop = F]
})

logical_results <- lapply(filtered_results, function(df) df < ALPHA_LEVEL)

sig_ratios <- lapply(logical_results, function(x) mean(as.logical(x), na.rm = T))


### 6 --- COMPARE WITH BRINGMANN AND FRIED --- ###

filtered_results.Bringmann2013 <- lapply(res_Bringmann2013, function(df) {
  df[rownames(df) %in% significance_tests, , drop = F]
})

filtered_results.Fried2022 <- lapply(res_Fried2022, function(df) {
  df[rownames(df) %in% significance_tests, , drop = F]
})

logical_results.Fried2022 <- lapply(filtered_results.Fried2022, function(df) df < ALPHA_LEVEL)

for (i in 1:length(logical_results.Bringmann2013)) {
  print(logical_results.Bringmann2013[[i]])
}

for (i in 1:length(logical_results.Fried2022)) {
  print(logical_results.Fried2022[[i]])
}

