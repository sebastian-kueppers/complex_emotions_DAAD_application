
### 0 --- PACKAGES AND HELPER FUNCTIONS --- ###

setwd("C:/Users/Sebastian Küppers/Desktop/Formal Theory of Co-Occuring Emotions (DFG project)/_PhD/_PhD_Study_1/complexity_in_emotion_ESM_data/SIM_ESM")

# Import helper functions for data generation
source('data_generation.R')

# Import helper functions for multifractality 
# (see https://osf.io/nm9b4/, Kelty-Stephen et al., 2023)
source('../CJspectrum.R')
source('../IAAFT.R')

# get packages
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
library(dplyr)
library(tidyr)
library(ggplot2)
# library(future)
# library(future.apply)
library(furrr)

# set up parallelization
n_cores <- parallel::detectCores() - 1
plan(multisession, workers = n_cores)


### 1 --- SETTINGS --- ###

tests <- c(
  "Bartels",
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

significance_tests <- c(
  "Bartels", 
  "TV-AR.p", 
  "KPSS.level",
  "KPSS.trend",
  "Keenan",
  "Tsay",
  "multifractality")

non_significance_tests <- c(
  "TV-AR.EDF",
  "PAFC.n",
  "PAFC.max_lag",
  "CP.n")

edm_tests <- c(
  "E_opt", 
  "pred_decay", 
  "rho_theta0", 
  "rho_theta_opt", 
  "theta_opt", 
  "delta_rho")

ALPHA_LEVEL <- 0.05 / length(significance_tests)


### 2 --- DATA GENERATION PARAMETERS --- ###

T_list <- c(25, 50, 75, 100, 150, 200)
N_list <- c(100)

# 2.1 - VAR(1) MODEL - ###

# coef matrix
A <- matrix(c(
  0.4,  0.1, -0.05, -0.05, 
  0.1,  0.4, -0.05, -0.05, 
  -0.05, -0.05,  0.4,  0.1, 
  -0.05, -0.05,  0.1,  0.4
), 4, 4, byrow = TRUE)


Sigma <- diag(4)

# Parameter settings lead to model where approx. 31% of all variance is 
# explained (check by fitting a VAR)


# 2.2 - BISTABLE COMPLEX MODEL - ###
C <- matrix(c(-.2, .04, -.2, -.2,
              .04, -.2, -.2, -.2,
              -.2, -.2, -.2, .04,
              -.2, -.2, .04, -.2), 4, 4, byrow = TRUE)

# Growth rates
r <- c(1, 1, 1, 1) 

# Total simulation time
no_minutes <- 14 * 24 * 60  # 2 weeks in minutes

### 3 --- SIMULATION LOOP USING HELPER FUNCTIONS --- ###

# 3.1 - VAR(1) MODEL - ###

data_list.VAR <- list()

set.seed(2310)

for (N in N_list) {
  data_list.VAR[[paste0("N", N)]] <- list()
  
  for (T in T_list) {
    
    # --- Simulate raw VAR1 data for N subjects ---
    raw_data <- simulate_var_subjects(N = N, T = T, A = A, Sigma = Sigma)
    
    # --- Likert 1-7 scaling ---
    likert_1_7 <- lapply(raw_data, function(df) {
      df_lik <- as.data.frame(lapply(df, min_max_to_scale, new_min = 1, new_max = 7))
      colnames(df_lik) <- paste0(colnames(df_lik), "_lik7")
      df_lik
    })
    
    # --- Likert 1-100 scaling ---
    likert_1_100 <- lapply(raw_data, function(df) {
      df_lik <- as.data.frame(lapply(df, min_max_to_scale, new_min = 1, new_max = 100))
      colnames(df_lik) <- paste0(colnames(df_lik), "_lik100")
      df_lik
    })
    
    # --- Save to data_list.VAR ---
    data_list.VAR[[paste0("N", N)]][[paste0("T", T)]] <- list(
      raw = raw_data,
      likert_1_7 = likert_1_7,
      likert_1_100 = likert_1_100
    )
    
    cat("Simulated: N =", N, "T =", T, "\n")
  }
}


# 3.2 - BISTABLE COMPLEX MODEL - ###

data_list.bistable <- list()

for (N in N_list) {
  data_list.bistable[[paste0("N", N)]] <- list()
  
  # --- Progress bar for participants ---
  pb_part <- txtProgressBar(min = 0, max = N, style = 3)
  
  # --- Simulate once per participant ---
  raw_data_list <- vector("list", N)
  
  for (i in 1:N) {
    # Initialize starting state
    init <- c(rnorm(1, 1.3, 1), 
              rnorm(1, 1.3, 1), 
              rnorm(1, 4.8, 1), 
              rnorm(1, 4.8, 1))
    
    # Simulate full bistable trajectory
    data_full <- simulate_bistable_trajectory(time = no_minutes,
                                              timestep = 0.01,
                                              p = 4,
                                              C = C,
                                              mu = 1.6,
                                              r = r,
                                              init = init,
                                              noiseSD = 4.5,
                                              noise = TRUE,
                                              pbar = FALSE)
    
    colnames(data_full) <- c("time",
                             "V1",
                             "V2",
                             "V3",
                             "V4")
    
    data_full <- subset(data_full, select = -c(time))
    
    raw_data_list[[i]] <- data_full
    
    # Update progress bar inside participant loop
    setTxtProgressBar(pb_part, i)
  }
  
  close(pb_part)  # close after all participants are simulated
  
  # --- Loop over T (number of observations per series) ---
  for (T_ in T_list) {
    
    # Downsample each participant
    raw_thinned <- lapply(raw_data_list, function(data) {
      thinner <- seq(1, nrow(data), length.out = T_)
      data_thinned <- data[thinner, , drop = FALSE]
      rownames(data_thinned) <- NULL
      data_thinned
    })
    
    # --- Likert 1-7 scaling ---
    likert_1_7 <- lapply(raw_thinned, function(df) {
      df_lik <- as.data.frame(lapply(df, min_max_to_scale, new_min = 1, new_max = 7))
      colnames(df_lik) <- paste0(colnames(df_lik), "_lik7")
      df_lik
    })
    
    # --- Likert 1-100 scaling ---
    likert_1_100 <- lapply(raw_thinned, function(df) {
      df_lik <- as.data.frame(lapply(df, min_max_to_scale, new_min = 1, new_max = 100))
      colnames(df_lik) <- paste0(colnames(df_lik), "_lik100")
      df_lik
    })
    
    # --- Save in data_list.bistable ---
    data_list.bistable[[paste0("N", N)]][[paste0("T", T_)]] <- list(
      raw = raw_thinned,
      likert_1_7 = likert_1_7,
      likert_1_100 = likert_1_100
    )
    
    cat("Simulated bistable: N =", N, "T =", T_, "\n")
  }
}

save(data_list.VAR, data_list.bistable, file = "simulated_data.RData")


### 4 --- APPLY TESTS --- ###

# define data sources
data_sources <- list(
  "VAR" = data_list.VAR,
  "bistable" = data_list.bistable
)

# Progress bars 
# n_blocks <- length(data_list)
# pb_outer <- txtProgressBar(min = 0, max = n_blocks, style = 3)

# Initialize two lists to store results for each scale
test_results_1_7 <- list()
test_results_1_7.bistable <- list()
test_results_1_100.VAR <- list()
test_results_1_100.bistable <- list()

# due to runtime, test only for likert 1_100
scales <- list(
  # "1_7" = "likert_1_7", 
  "1_100" = "likert_1_100"
)

### 4.1 - SIGNIFICANCE TESTS (WITHOUT MULTIFRACTALITY) - ###

# sig_test_results_1_7.VAR <- list()
# sig_test_results_1_7.bistable <- list()
sig_test_results_1_100.VAR <- list()
sig_test_results_1_100.bistable <- list()

# Define function to be run in parallel

run_tests_for_person <- function(ts_data, cols, significance_tests_main, ALPHA_LEVEL) {
  
  res_mat <- matrix(
    NA,
    nrow = length(significance_tests_main),
    ncol = length(cols),
    dimnames = list(significance_tests_main, cols)
  )
  
  for (j in seq_along(cols)) {
    col <- cols[j]
    ts  <- ts_data[[col]]
    
    # --- Bartels ---
    bartels_p <- tryCatch(bartels.rank.test(ts)$p.value, error = function(e) NA)
    if (length(bartels_p) == 0) bartels_p <- NA
    res_mat["Bartels", col] <- bartels_p

    # --- KPSS level ---
    kpss_level_p <- tryCatch(kpss.test(ts, null = "Level", lshort = TRUE)$p.value, error=function(e) NA)
    if (length(kpss_level_p) == 0) kpss_level_p <- NA
    res_mat["KPSS.level", col] <- kpss_level_p

    # --- KPSS trend ---
    kpss_trend_p <- tryCatch(kpss.test(ts, null = "Trend", lshort = TRUE)$p.value, error=function(e) NA)
    if (length(kpss_trend_p) == 0) kpss_trend_p <- NA
    res_mat["KPSS.trend", col] <- kpss_trend_p

    # --- Keenan ---
    keenan_p <- tryCatch(keenanTest(ts)$p.value, error=function(e) NA)
    if (length(keenan_p) == 0) keenan_p <- NA
    res_mat["Keenan", col] <- keenan_p

    # --- Tsay ---
    tsay_p <- tryCatch(tsayTest(ts)$p.value, error=function(e) NA)
    if (length(tsay_p) == 0) tsay_p <- NA
    res_mat["Tsay", col] <- tsay_p
    
    # --- TV-AR GAM ---
    N_ts <- length(ts)
    if (N_ts > 3) {
      tt <- 1:(N_ts - 1)
      tv <- tryCatch(
        gam(
          ts[2:N_ts] ~ s(tt, by = ts[1:(N_ts - 1)], 
                         k = min(10, N_ts - 2), bs = "tp")
        ),
        error = function(e) NULL
      )
      
      if (!is.null(tv)) {
        stvar <- summary(tv)
        if (!is.null(stvar$s.table)) {
          res_mat["TV-AR.p", col] <- stvar$s.table[1, "p-value"]
        }
      }
    }
  }
  
  res_mat
}


# Loop over data sources
for (source_name in names(data_sources)) {
  
  data_list <- data_sources[[source_name]]
  
  # Loop over scales (only 1_100 for runtime)
  for (scale_name in names(scales)) {
    
    test_results <- list()  # temporary storage for this scale
    
    # Loop over sample sizes N
    for (N_name in names(data_list)) {
      
      test_results[[N_name]] <- list()
      
      # Loop over time series lengths T
      for (T_name in names(data_list[[N_name]])) {
        
        lik_data <- data_list[[N_name]][[T_name]][[scales[[scale_name]]]]  # List of persons
        N_persons <- length(lik_data)
        if (N_persons == 0) next
        
        cols <- colnames(lik_data[[1]])
        significance_tests_main <- setdiff(significance_tests, "multifractality")
        
        res_array <- array(
          NA,
          dim = c(length(significance_tests_main), length(cols), N_persons),
          dimnames = list(significance_tests_main, cols, paste0("P", 1:N_persons))
        )
        
        cat("\nRunning:", source_name, scale_name, N_name, T_name, "\n")
        
        # Loop over persons
        
        person_results <- future_map(
          lik_data,
          ~ run_tests_for_person(
            ts_data = .x,
            cols = cols,
            significance_tests_main = significance_tests_main,
            ALPHA_LEVEL = ALPHA_LEVEL
          ),
          .options = furrr_options(seed = TRUE)
        )
        
        for (p in seq_len(N_persons)) {
          res_array[ , , p] <- person_results[[p]]
        }
        
        # Aggregate across persons
        sig_rates <- apply(res_array, c(1,2), function(x) mean(x < ALPHA_LEVEL, na.rm=TRUE))
        
        test_results[[N_name]][[T_name]] <- list(
          raw_array = res_array,
          sig_rate = sig_rates
        )
        
        cat("Finished:", source_name, scale_name, N_name, T_name, "\n")
      }
    }
    
    # Save in proper results object
    if (scale_name == "1_7") {
      if (source_name == "VAR") sig_test_results_1_7.VAR <- test_results
      if (source_name == "bistable") sig_test_results_1_7.bistable <- test_results
    } else if (scale_name == "1_100") {
      if (source_name == "VAR") sig_test_results_1_100.VAR <- test_results
      if (source_name == "bistable") sig_test_results_1_100.bistable <- test_results
    }
  }
}


### 4.2 - MULTIFRACTALITY - ###

multifractality_pval <- function(ts) {
  
  ts <- na.omit(ts)
  N  <- length(ts)
  
  spectrum <- tryCatch(
    CJspectrum(ts, 1, N, 1),
    error = function(e) NULL
  )
  if (is.null(spectrum)) return(NA_real_)
  
  n_surr <- ceiling(
    pwr.t.test(
      d = 0.5,
      sig.level = 0.05,
      power = 0.80,
      type = "one.sample"
    )$n
  )
  
  surrogates <- tryCatch(
    IAAFT(as.numeric(ts), n_surr, 1000),
    error = function(e) NULL
  )
  if (is.null(surrogates)) return(NA_real_)
  
  spectra <- numeric(0)
  for (i in seq_len(ncol(surrogates))) {
    spec_i <- tryCatch(
      CJspectrum(surrogates[, i], 1, N, 1),
      error = function(e) NA
    )
    if (!is.na(spec_i))
      spectra <- c(spectra, spec_i)
  }
  
  if (length(spectra) <= 1) return(NA_real_)
  
  tryCatch(
    t.test(spectra, mu = spectrum)$p.value,
    error = function(e) NA_real_
  )
}

run_mf_for_person <- function(ts_data, cols) {
  
  res <- numeric(length(cols))
  names(res) <- cols
  
  for (col in cols) {
    res[col] <- multifractality_pval(ts_data[[col]])
  }
  
  res
}


# set.seed(2310)

for (source_name in names(data_sources)) {
  
  data_list <- data_sources[[source_name]]
  
  for (scale_name in names(scales)) {
    
    if (scale_name == "1_7") {
      if (source_name == "VAR") existing_results <- sig_test_results_1_7.VAR
      if (source_name == "bistable") existing_results <- sig_test_results_1_7.bistable
    } else if (scale_name == "1_100") {
      if (source_name == "VAR") existing_results <- sig_test_results_1_100.VAR
      if (source_name == "bistable") existing_results <- sig_test_results_1_100.bistable
    }
    
    # Choose the proper sig_test_results object
    sig_test_results <- list()
    
    for (N_name in names(data_list)) {
      
      sig_test_results[[N_name]] <- list()
      
      for (T_name in names(data_list[[N_name]])) {
        
        lik_data <- data_list[[N_name]][[T_name]][[scales[[scale_name]]]]
        N_persons <- length(lik_data)
        if (N_persons == 0) next
        
        cols <- colnames(lik_data[[1]])
        
        # Raw p-values matrix: Variables × Persons
        mf_array <- matrix(
          NA,
          nrow = length(cols),
          ncol = N_persons,
          dimnames = list(cols, paste0("P", 1:N_persons))
        )
        
        cat(
          "\nMF:", source_name,
          "| scale:", scale_name,
          "|", N_name,
          "|", T_name,
          "| persons:", N_persons, "\n"
        )
        
        # ---- PARALLEL EXECUTION ----
        # person_results <- future_lapply(
        #   lik_data,
        #   run_mf_for_person,
        #   cols = cols,
        #   future.seed = TRUE
        # )
        person_results <- future_map(
          lik_data,
          ~ run_mf_for_person(
            ts_data = .x,
            cols = cols
          ),
          .options = furrr_options(seed = TRUE)  # preserves reproducibility
        )
        
        # Combine into matrix: vars × persons
        mf_array <- do.call(cbind, person_results)
        colnames(mf_array) <- paste0("P", seq_len(ncol(mf_array)))
        rownames(mf_array) <- cols
        
        # Aggregate across persons
        mf_sig_rate <- apply(mf_array, 1, function(x) mean(x < ALPHA_LEVEL, na.rm = TRUE))
        
        # Merge with existing sig_rate and raw_array if they exist
        existing_sig_rate  <- existing_results[[N_name]][[T_name]]$sig_rate
        existing_raw_array <- existing_results[[N_name]][[T_name]]$raw_array
        
        if (!is.null(existing_sig_rate)) {
          new_sig_rate <- rbind(existing_sig_rate, multifractality = mf_sig_rate)
        } else {
          new_sig_rate <- mf_sig_rate
        }
        
        if (!is.null(existing_raw_array)) {
          new_raw_array <- abind::abind(existing_raw_array, multifractality = mf_array, along = 1)
        } else {
          new_raw_array <- mf_array
        }
        
        # Save in sig_test_results
        sig_test_results[[N_name]][[T_name]] <- list(
          sig_rate = new_sig_rate,
          raw_array = new_raw_array
        )
        
        cat("MF DONE:", source_name, scale_name, N_name, T_name, "\n")
        
      } # end T loop
    } # end N loop
    
    # Assign results to the proper output object
    if (scale_name == "1_7") {
      if (source_name == "VAR") sig_test_results_1_7.VAR <- sig_test_results
      if (source_name == "bistable") sig_test_results_1_7.bistable <- sig_test_results
    } else if (scale_name == "1_100") {
      if (source_name == "VAR") sig_test_results_1_100.VAR <- sig_test_results
      if (source_name == "bistable") sig_test_results_1_100.bistable <- sig_test_results
    }
    
  } # end scale loop
} # end source loop


### 4.3 - OTHER PARAMETRIC TESTS - ###

run_tv_cp_pacf_for_person <- function(ts_data, cols, ALPHA_LEVEL) {
  res <- matrix(NA, nrow = length(non_significance_tests), ncol = length(cols),
                dimnames = list(non_significance_tests, cols))
  
  for (col in cols) {
    ts <- ts_data[[col]]
    N_ts <- length(ts)
    
    # --- TV-AR GAM ---
    if (N_ts > 2) {
      tt <- 1:(N_ts-1)
      tv <- tryCatch(
        gam(ts[2:N_ts] ~ s(tt, by=ts[1:(N_ts-1)], k=min(10, N_ts-2), bs="tp")),
        error = function(e) NULL
      )
      if (!is.null(tv)) {
        stvar <- summary(tv)
        if (!is.null(stvar$s.table) && nrow(stvar$s.table) >= 1) {
          res["TV-AR.EDF", col] <- round(stvar$s.table[1, "edf"], 3)
        }
      }
    }
    
    # --- Change Point Analysis ---
    R <- ceiling(20 / ALPHA_LEVEL)
    cp.out <- tryCatch(e.divisive(matrix(ts), R=R, sig.lvl=ALPHA_LEVEL), error=function(e) NULL)
    if (!is.null(cp.out)) {
      res["CP.n", col] <- length(which(cp.out$p.values < ALPHA_LEVEL))
    }
    
    # --- PACF analysis ---
    pacf_vals <- tryCatch(pacf(ts, plot=FALSE)$acf, error=function(e) NA)
    if (is.numeric(pacf_vals) && length(pacf_vals) > 1 && !all(is.na(pacf_vals))) {
      thresh <- 2 / sqrt(N_ts)
      sig_lags <- which(abs(pacf_vals) > thresh)
      res["PAFC.n", col] <- length(sig_lags)
      res["PAFC.max_lag", col] <- if(length(sig_lags) > 0) max(sig_lags) else 0
    }
  }
  
  res
}


# Loop over data sources

for (source_name in names(data_sources)) {
  data_list <- data_sources[[source_name]]
  
  for (scale_name in names(scales)) {
    results_list <- list()
    
    for (N_name in names(data_list)) {
      results_list[[N_name]] <- list()
      
      for (T_name in names(data_list[[N_name]])) {
        lik_data <- data_list[[N_name]][[T_name]][[scales[[scale_name]]]]
        N_persons <- length(lik_data)
        if (N_persons == 0) next
        
        cols <- colnames(lik_data[[1]])
        
        # ---- PARALLEL EXECUTION ----
        person_results <- future_map(
          lik_data,
          ~ run_tv_cp_pacf_for_person(
            ts_data = .x,
            cols = cols,
            ALPHA_LEVEL = ALPHA_LEVEL
          ),
          .options = furrr_options(seed = TRUE)  # ensures reproducibility
        )
        
        # Combine results
        res_array <- array(
          NA, 
          dim = c(length(non_significance_tests), length(cols), N_persons),
          dimnames = list(non_significance_tests, cols, paste0("P", 1:N_persons))
        )
        
        for (p in seq_len(N_persons)) {
          res_array[,,p] <- person_results[[p]]
        }
        
        results_list[[N_name]][[T_name]] <- res_array
        
        cat("TV/CP/PACF DONE:", source_name, scale_name, N_name, T_name, "\n")
      }
    }
    
    # Assign results to proper output object
    if (scale_name == "1_7") {
      if (source_name == "VAR") tv_cp_pacf_results_1_7.VAR <- results_list
      if (source_name == "bistable") tv_cp_pacf_results_1_7.bistable <- results_list
    } else if (scale_name == "1_100") {
      if (source_name == "VAR") tv_cp_pacf_results_1_100.VAR <- results_list
      if (source_name == "bistable") tv_cp_pacf_results_1_100.bistable <- results_list
    }
    
  }
}

save(tv_cp_pacf_results_1_100.VAR,
     tv_cp_pacf_results_1_100.bistable,
     file = "non-sig_test_results.RData")

### 4.4 - EDM TESTS - ###

for (source_name in names(data_sources)) {
  
  data_list <- data_sources[[source_name]]
  
  for (scale_name in names(scales)) {
    
    results_list <- list()
    
    for (N_name in names(data_list)) {
      
      results_list[[N_name]] <- list()
      
      for (T_name in names(data_list[[N_name]])) {
        
        lik_data <- data_list[[N_name]][[T_name]][[scales[[scale_name]]]]
        N_persons <- length(lik_data)
        if (N_persons == 0) next
        
        cols <- colnames(lik_data[[1]])
        
        # Initialize arrays to store results
        res_array <- array(
          NA,
          dim = c(length(edm_tests), length(cols), N_persons),
          dimnames = list(edm_tests, cols, paste0("P", 1:N_persons))
        )

        cat("\nSource:", source_name, "| Scale:", scale_name, "| N:", N_name, "| T:", T_name, "| Persons:", N_persons, "\n")
        
        # Parallelized person loop
        person_results <- future_map(
          lik_data,
          function(ts_data) {
            res_person <- matrix(NA, nrow = length(edm_tests), ncol = length(cols),
                                 dimnames = list(edm_tests, cols))
            
            for (col in cols) {
              ts <- ts_data[[col]]
              # print(ts)
              N <- length(ts)
              
              tau <- tryCatch(
                suppressWarnings(timeLag(unlist(ts), technique = "ami",
                                         selection.method = "first.minimum", lag.max = 10, do.plot = FALSE)),
                error = function(e) NA
              )
              
              # FALLBACK
              # if (is.na(tau)) tau <- 1 
              
              if (!is.na(tau)) {
                mid <- floor(N / 2)
                lib <- c(1, mid)
                pred <- c(mid + 1, N)
                
                # maximum embed dimension to test depends on data length
                Tp <- 1
                maxE_lib <- floor((length(1:mid) - (Tp-1))/tau)     
                maxE_pred <- floor((length((mid+1):N) - (Tp-1))/tau)
                maxE_possible <- min(10, maxE_lib, maxE_pred)

                e <- tryCatch({
                  emb_out <- suppressWarnings(EmbedDimension(
                    dataFrame = as.data.frame(ts),
                    lib = lib,
                    pred = pred,
                    maxE = maxE_possible,
                    Tp = 1,
                    tau = tau,
                    columns = 'ts',
                    target = 'ts',
                    noTime = TRUE,
                    showPlot = FALSE
                  ))
                  emb_out$E[which.max(emb_out$rho)]
                }, error = function(e) NA)

                # print(e)
                if (!is.na(e)) {
                  res_person["E_opt", col] <- e

                  thetas <- seq(0, 8, by = 0.5)
                  rho_theta <- rep(NA, length(thetas))

                  for (k in seq_along(thetas)) {
                    sm <- tryCatch(
                      suppressWarnings(SMap(
                        dataFrame = as.data.frame(ts),
                        lib = lib,
                        pred = pred,
                        columns = "ts",
                        target = "ts",
                        E = e,
                        tau = tau,
                        Tp = 1,
                        theta = thetas[k],
                        noTime = TRUE
                      )),
                      error = function(e) NULL
                    )
                    
                    if (!is.null(sm)) {
                      sm.obs <- sm$predictions$Observations
                      sm.pred <- sm$predictions$Predictions
                      valid <- is.finite(sm.obs) & is.finite(sm.pred)
                      
                      rho <- cor(sm.obs[valid], sm.pred[valid])
                      rho_theta[k] <- rho
                    }
                  }

                  if (any(!is.na(rho_theta))) {
                    res_person["rho_theta0", col] <- rho_theta[thetas == 0]
                    res_person["rho_theta_opt", col] <- max(rho_theta, na.rm = TRUE)
                    res_person["theta_opt", col] <- thetas[which.max(rho_theta)]
                    res_person["delta_rho", col] <- res_person["rho_theta_opt", col] - res_person["rho_theta0", col]
                  }
                  
                  # maximum prediction horizon again depends on pred
                  Tp_max <- N - min(pred) - 1
                  Tp_max <- max(Tp_max, 1) 

                  pi <- tryCatch(suppressWarnings(PredictInterval(
                    dataFrame = as.data.frame(ts),
                    noTime = TRUE,
                    columns = "ts",
                    target = "ts",
                    E = e,
                    tau = tau,
                    lib = lib,
                    pred = pred,
                    maxTp = Tp_max,
                    showPlot = FALSE
                  )), error = function(e) NULL)

                  if (!is.null(pi) && nrow(pi) >= 5) {
                    decay <- tryCatch({
                      fit <- lm(rho ~ Tp, data = pi[1:5, ])
                      coef(fit)[2] * 5
                    }, error = function(e) NA)
                    res_person["pred_decay", col] <- decay
                  }
                }
              }
            }
            
            res_person
          },
          .options = furrr_options(seed = TRUE)
        )
        
        
        # Aggregate across persons if needed
        
        results_list[[N_name]][[T_name]] <- person_results
        
        cat("EDM Tests DONE:", source_name, scale_name, N_name, T_name, "\n")
        
      } # end T loop
    } # end N loop
    
    # Assign results to proper output object
    if (scale_name == "1_7") {
      if (source_name == "VAR") edm_results_1_7.VAR <- results_list
      if (source_name == "bistable") edm_results_1_7.bistable <- results_list
    } else if (scale_name == "1_100") {
      if (source_name == "VAR") edm_results_1_100.VAR <- results_list
      if (source_name == "bistable") edm_results_1_100.bistable <- results_list
    }
    
  } # end scale loop
} # end data source loop

save(edm_results_1_100.VAR, 
     edm_results_1_100.bistable,
     file = "EDM_test_results.RData")


### 5 --- POWERPOINT ---------- ###

# plot time series for illustrative purposes
ts.VAR <- ts(data_list.VAR$N100$T100$likert_1_100[[1]]$V1_lik100)
ts.bistable <- ts(data_list.bistable$N100$T100$likert_1_100[[1]]$V1_lik100)

plot(ts.VAR, col = "red")
plot(ts.bistable, col = "blue")



### 5 --- AGGREGATE RESULTS --- ###

### 5.X EDM -------------------- ###

# bistable and VAR lists
lst_bi  <- edm_results_1_100.bistable$N100$T200
lst_var <- edm_results_1_100.VAR$N100$T200

# extract all theta_opt values
theta_bi <- unlist(lapply(lst_bi, function(df) df["theta_opt", ]))
theta_var <- unlist(lapply(lst_var, function(df) df["theta_opt", ]))

# remove NA
theta_bi  <- theta_bi[!is.na(theta_bi)]
theta_var <- theta_var[!is.na(theta_var)]

summary(theta_bi)
summary(theta_var)

quantile(theta_bi, probs = c(.25, .5, .75))
quantile(theta_var, probs = c(.25, .5, .75))

mean(theta_bi == 0)
mean(theta_var == 0)



lst.VAR <- edm_results_1_100.VAR$N100$T200

# get row names (assumed identical across datasets)
rows <- rownames(lst.VAR[[1]])

# compute mean per row across all datasets and columns
row_means <- sapply(rows, function(r) {
  vals <- unlist(lapply(lst.VAR, function(df) df[r, ]))
  mean(vals, na.rm = TRUE)
})

# convert to data frame for readability
result.VAR <- data.frame(
  row = names(row_means),
  mean_value = row_means,
  row.names = NULL
)

result.VAR



lst.bistable <- edm_results_1_100.bistable$N100$T200

# get row names (assumed identical across datasets)
rows <- rownames(lst.bistable[[1]])

# compute mean per row across all datasets and columns
row_means <- sapply(rows, function(r) {
  vals <- unlist(lapply(lst.bistable, function(df) df[r, ]))
  mean(vals, na.rm = TRUE)
})

# convert to data frame for readability
result.bistable <- data.frame(
  row = names(row_means),
  mean_value = row_means,
  row.names = NULL
)

result.bistable


### --- AGGREGATE RESULTS ------ ###

### --- SIGNIFICANCE TESTS --- ###

# Initialize a list to store aggregated results
agg_list.VAR <- list()

# Loop over sample sizes (N)
for (N_name in names(sig_test_results_1_100.VAR)) {
  
  # Loop over time series lengths (T)
  for (T_name in names(sig_test_results_1_100.VAR[[N_name]])) {
    
    res_array <- sig_test_results_1_100.VAR[[N_name]][[T_name]]$raw_array  # get raw p-values array

    # Initialize vector to store proportion of significant p-values per test
    ratio_sig <- numeric(dim(res_array)[1])
    names(ratio_sig) <- dimnames(res_array)[[1]]

    # Aggregate over all variables × persons
    for (test_name in dimnames(res_array)[[1]]) {
      vals <- as.vector(res_array[test_name,,])  # flatten variable × person values
      ratio_sig[test_name] <- mean(vals < ALPHA_LEVEL, na.rm = TRUE)  # ratio (0-1)
    }

    # Store as a data frame
    agg_list.VAR[[paste0(N_name, "_", T_name)]] <- data.frame(
      Test = names(ratio_sig),
      model = "VAR",
      N = N_name,
      T = T_name,
      ratioSig = ratio_sig,  # ratio is the last column
      row.names = NULL
    )
    
  } # end T loop
} # end N loop

# Combine all data frames into one tidy data frame
agg_df.VAR <- bind_rows(agg_list.VAR)


# Initialize a list to store aggregated results
agg_list.bistable <- list()

# Loop over sample sizes (N)
for (N_name in names(sig_test_results_1_100.bistable)) {
  
  # Loop over time series lengths (T)
  for (T_name in names(sig_test_results_1_100.bistable[[N_name]])) {
    
    res_array <- sig_test_results_1_100.bistable[[N_name]][[T_name]]$raw_array  # get raw p-values array
    
    # Initialize vector to store proportion of significant p-values per test
    ratio_sig <- numeric(dim(res_array)[1])
    names(ratio_sig) <- dimnames(res_array)[[1]]
    
    # Aggregate over all variables × persons
    for (test_name in dimnames(res_array)[[1]]) {
      vals <- as.vector(res_array[test_name,,])  # flatten variable × person values
      ratio_sig[test_name] <- mean(vals < ALPHA_LEVEL, na.rm = TRUE)  # ratio (0-1)
    }
    
    # Store as a data frame
    agg_list.bistable[[paste0(N_name, "_", T_name)]] <- data.frame(
      Test = names(ratio_sig),
      model = "bistable",
      N = N_name,
      T = T_name,
      ratioSig = ratio_sig,  # ratio is the last column
      row.names = NULL
    )
    
  } # end T loop
} # end N loop

# Combine all data frames into one tidy data frame
agg_df.bistable <- bind_rows(agg_list.bistable)


# Combine into a single df for significance tests
agg_df.sig <- rbind(agg_df.VAR,
                    agg_df.bistable)

### --- OTHER PARAMETRIC TESTS --- ###

# Initialize a list to store aggregated results
agg_list.VAR.other <- list()

# Loop over sample sizes (N)
for (N_name in names(tv_cp_pacf_results_1_100.VAR)) {
  
  # Loop over time series lengths (T)
  for (T_name in names(tv_cp_pacf_results_1_100.VAR[[N_name]])) {
    
    res_array <- tv_cp_pacf_results_1_100.VAR[[N_name]][[T_name]]  
    res_mean <- apply(res_array, 1, mean, na.rm = TRUE)
    
    # Store as a data frame
    agg_list.VAR.other[[paste0(N_name, "_", T_name)]] <- data.frame(
      Test = names(res_mean),
      model = "VAR",
      N = N_name,
      T = T_name,
      value = res_mean,  # ratio is the last column
      row.names = NULL
    )
    
  } # end T loop
} # end N loop

# Combine all data frames into one tidy data frame
agg_df.VAR.other <- bind_rows(agg_list.VAR.other)


# Initialize a list to store aggregated results
agg_list.bistable.other <- list()

# Loop over sample sizes (N)
for (N_name in names(tv_cp_pacf_results_1_100.bistable)) {
  
  # Loop over time series lengths (T)
  for (T_name in names(tv_cp_pacf_results_1_100.bistable[[N_name]])) {
    
    res_array <- tv_cp_pacf_results_1_100.bistable[[N_name]][[T_name]]  
    res_mean <- apply(res_array, 1, mean, na.rm = TRUE)
    
    # Store as a data frame
    agg_list.bistable.other[[paste0(N_name, "_", T_name)]] <- data.frame(
      Test = names(res_mean),
      model = "bistable",
      N = N_name,
      T = T_name,
      value = res_mean,  # ratio is the last column
      row.names = NULL
    )
    
  } # end T loop
} # end N loop

# Combine all data frames into one tidy data frame
agg_df.bistable.other <- bind_rows(agg_list.bistable.other)

# Combine both dfs together
agg_df.other <- rbind(agg_df.VAR.other,
                      agg_df.bistable.other)


### --- EDM TESTS --- ###

# Initialize a list to store aggregated results
agg_list.VAR.EDM <- list()

# Loop over sample sizes (N)
for (N_name in names(edm_results_1_100.VAR)) {
  
  # Loop over time series lengths (T)
  for (T_name in names(edm_results_1_100.VAR[[N_name]])) {
    
    res_array <- simplify2array(edm_results_1_100.VAR[[N_name]][[T_name]])  
    res_mean <- apply(res_array, 1, mean, na.rm = TRUE)
    
    # Store as a data frame
    agg_list.VAR.EDM[[paste0(N_name, "_", T_name)]] <- data.frame(
      Test = names(res_mean),
      model = "VAR",
      N = N_name,
      T = T_name,
      value = res_mean,  # ratio is the last column
      row.names = NULL
    )
    
  } # end T loop
} # end N loop

# Combine all data frames into one tidy data frame
agg_df.VAR.EDM <- bind_rows(agg_list.VAR.EDM)


# Initialize a list to store aggregated results
agg_list.bistable.EDM <- list()

# Loop over sample sizes (N)
for (N_name in names(edm_results_1_100.bistable)) {
  
  # Loop over time series lengths (T)
  for (T_name in names(edm_results_1_100.bistable[[N_name]])) {
    
    res_array <- simplify2array(edm_results_1_100.bistable[[N_name]][[T_name]])  
    res_mean <- apply(res_array, 1, mean, na.rm = TRUE)
    
    # Store as a data frame
    agg_list.bistable.EDM[[paste0(N_name, "_", T_name)]] <- data.frame(
      Test = names(res_mean),
      model = "bistable",
      N = N_name,
      T = T_name,
      value = res_mean,  # ratio is the last column
      row.names = NULL
    )
    
  } # end T loop
} # end N loop

# Combine all data frames into one tidy data frame
agg_df.bistable.EDM <- bind_rows(agg_list.bistable.EDM)

# Combine into one df
agg_df.EDM <- rbind(agg_df.VAR.EDM,
                    agg_df.bistable.EDM)

### 6 --- PLOT --- ###

### 6.1 --- SIGNIFICANCE TESTS ---

plot_data.VAR <- agg_df.sig %>%
  filter(
    N == "N100",
    model == "VAR"
  ) %>%
  mutate(T = factor(T, levels = paste0("T", sort(T_list))))

plot_data.bistable <- agg_df.sig %>%
  filter(
    N == "N100",
    model == "bistable"    
  ) %>%
  mutate(T = factor(T, levels = paste0("T", sort(T_list))))

# Create the line plot
ggplot(plot_data.VAR, aes(x = T, y = ratioSig, color = Test, group = Test)) +
  geom_line(size = 1) +             # lines for each Test
  geom_point(size = 2) +            # points at each T
  labs(
    color = "Test"
  ) +
  theme_minimal(base_size = 14) +   # clean theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

# Create the line plot
ggplot(plot_data.bistable, aes(x = T, y = ratioSig, color = Test, group = Test)) +
  geom_line(size = 1) +             # lines for each Test
  geom_point(size = 2) +            # points at each T
  labs(
    color = "Test"
  ) +
  theme_minimal(base_size = 14) +   # clean theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )



### 6.1 --- OTHER PARAMETRIC TESTS ---

plot_data.other <- agg_df.other %>%
  filter(N == "N100") %>%
  mutate(
    T = factor(T, levels = paste0("T", sort(T_list))),
    model = factor(model, levels = c("VAR", "bistable"))
  )

test_colors <- c("TV-AR.EDF" = "steelblue",
                 "PAFC.n" = "darkgreen",
                 "PAFC.max_lag" = "orange",
                 "CP.n" = "purple")

# TV-AR.EDF plot
ggplot(
  filter(plot_data.other, Test == "TV-AR.EDF"),
  aes(
    x = T,
    y = value,
    group = model,
    linetype = model
  )
) +
  geom_line(aes(color = Test), size = 1.2) +
  geom_point(aes(color = Test), size = 2) +
  scale_color_manual(values = test_colors) +
  scale_linetype_manual(values = c("VAR" = "solid", "bistable" = "dashed")) +
  scale_y_continuous(limits = c(1, NA)) +   
  labs(
    title = "TV-AR.EDF",
    x = "T (time series length)",
    y = NULL,
    linetype = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.minor = element_line(color = "gray90", linetype = "dotted")
  ) + 
  guides(color = "none")

# PAFC.n
ggplot(
  filter(plot_data.other, Test == "PAFC.n"),
  aes(
    x = T,
    y = value,
    group = model,
    linetype = model
  )
) +
  geom_line(aes(color = Test), size = 1.2) +
  geom_point(aes(color = Test), size = 2) +
  scale_color_manual(values = test_colors) +
  scale_linetype_manual(values = c("VAR" = "solid", "bistable" = "dashed")) +
  scale_y_continuous(limits = c(0, NA)) +   
  labs(
    title = "PAFC.n",
    x = "T (time series length)",
    y = NULL,
    linetype = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.minor = element_line(color = "gray90", linetype = "dotted")
  ) + 
  guides(color = "none")

# PAFC.n
ggplot(
  filter(plot_data.other, Test == "PAFC.n"),
  aes(
    x = T,
    y = value,
    group = model,
    linetype = model
  )
) +
  geom_line(aes(color = Test), size = 1.2) +
  geom_point(aes(color = Test), size = 2) +
  scale_color_manual(values = test_colors) +
  scale_linetype_manual(values = c("VAR" = "solid", "bistable" = "dashed")) +
  scale_y_continuous(limits = c(0, NA)) +   
  labs(
    title = "PAFC.n",
    x = "T (time series length)",
    y = NULL,
    linetype = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.minor = element_line(color = "gray90", linetype = "dotted")
  ) + 
  guides(color = "none")

# PAFC.max_lag
ggplot(
  filter(plot_data.other, Test == "PAFC.max_lag"),
  aes(
    x = T,
    y = value,
    group = model,
    linetype = model
  )
) +
  geom_line(aes(color = Test), size = 1.2) +
  geom_point(aes(color = Test), size = 2) +
  scale_color_manual(values = test_colors) +
  scale_linetype_manual(values = c("VAR" = "solid", "bistable" = "dashed")) +
  scale_y_continuous(limits = c(0, NA)) +   
  labs(
    title = "PAFC.max_lag",
    x = "T (time series length)",
    y = NULL,
    linetype = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.minor = element_line(color = "gray90", linetype = "dotted")
  ) + 
  guides(color = "none")

# CP.n
ggplot(
  filter(plot_data.other, Test == "CP.n"),
  aes(
    x = T,
    y = value,
    group = model,
    linetype = model
  )
) +
  geom_line(aes(color = Test), size = 1.2) +
  geom_point(aes(color = Test), size = 2) +
  scale_color_manual(values = test_colors) +
  scale_linetype_manual(values = c("VAR" = "solid", "bistable" = "dashed")) +
  scale_y_continuous(limits = c(0, NA)) +   
  labs(
    title = "CP.n",
    x = "T (time series length)",
    y = NULL,
    linetype = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.minor = element_line(color = "gray90", linetype = "dotted")
  ) + 
  guides(color = "none")

### 6.3 --- EDM TESTS ---

plot_data.EDM <- agg_df.EDM %>%
  filter(N == "N100") %>%
  mutate(
    T = factor(T, levels = paste0("T", sort(T_list))),
    model = factor(model, levels = c("VAR", "bistable"))
  )

test_colors <- c("E_opt" = "steelblue",
                 "pred_decay" = "darkgreen",
                 "rho_theta0" = "orange",
                 "rho_theta_opt" = "purple",
                 "theta_opt" = "red",
                 "delta_rho" = "darkblue")
# E_opt
ggplot(
  filter(plot_data.EDM, Test == "E_opt"),
  aes(
    x = T,
    y = value,
    group = model,
    linetype = model
  )
) +
  geom_line(aes(color = Test), size = 1.2) +
  geom_point(aes(color = Test), size = 2) +
  scale_color_manual(values = test_colors) +
  scale_linetype_manual(values = c("VAR" = "solid", "bistable" = "dashed")) +
  scale_y_continuous(limits = c(0, NA)) +   
  labs(
    title = "E_opt",
    x = "T (time series length)",
    y = NULL,
    linetype = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.minor = element_line(color = "gray90", linetype = "dotted")
  ) +
  guides(color = "none")

# pred_decay
ggplot(
  filter(plot_data.EDM, Test == "pred_decay"),
  aes(
    x = T,
    y = value,
    group = model,
    linetype = model
  )
) +
  geom_line(aes(color = Test), size = 1.2) +
  geom_point(aes(color = Test), size = 2) +
  scale_color_manual(values = test_colors) +
  scale_linetype_manual(values = c("VAR" = "solid", "bistable" = "dashed")) +
  scale_y_continuous(limits = c(NA, 0)) +   
  labs(
    title = "pred_decay",
    x = "T (time series length)",
    y = NULL,
    linetype = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.minor = element_line(color = "gray90", linetype = "dotted")
  ) +
  guides(color = "none")

# rho_theta0
ggplot(
  filter(plot_data.EDM, Test == "rho_theta0"),
  aes(
    x = T,
    y = value,
    group = model,
    linetype = model
  )
) +
  geom_line(aes(color = Test), size = 1.2) +
  geom_point(aes(color = Test), size = 2) +
  scale_color_manual(values = test_colors) +
  scale_linetype_manual(values = c("VAR" = "solid", "bistable" = "dashed")) +
  scale_y_continuous(limits = c(0, NA)) +   
  labs(
    title = "rho_theta0",
    x = "T (time series length)",
    y = NULL,
    linetype = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.minor = element_line(color = "gray90", linetype = "dotted")
  ) +
  guides(color = "none")


# rho_theta_opt
ggplot(
  filter(plot_data.EDM, Test == "rho_theta_opt"),
  aes(
    x = T,
    y = value,
    group = model,
    linetype = model
  )
) +
  geom_line(aes(color = Test), size = 1.2) +
  geom_point(aes(color = Test), size = 2) +
  scale_color_manual(values = test_colors) +
  scale_linetype_manual(values = c("VAR" = "solid", "bistable" = "dashed")) +
  scale_y_continuous(limits = c(0, NA)) +   
  labs(
    title = "rho_theta_opt",
    x = "T (time series length)",
    y = NULL,
    linetype = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.minor = element_line(color = "gray90", linetype = "dotted")
  ) +
  guides(color = "none")

# theta_opt
ggplot(
  filter(plot_data.EDM, Test == "theta_opt"),
  aes(
    x = T,
    y = value,
    group = model,
    linetype = model
  )
) +
  geom_line(aes(color = Test), size = 1.2) +
  geom_point(aes(color = Test), size = 2) +
  scale_color_manual(values = test_colors) +
  scale_linetype_manual(values = c("VAR" = "solid", "bistable" = "dashed")) +
  scale_y_continuous(limits = c(0, NA)) +   
  labs(
    title = "theta_opt",
    x = "T (time series length)",
    y = NULL,
    linetype = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.minor = element_line(color = "gray90", linetype = "dotted")
  ) +
  guides(color = "none")

# delta_rho
ggplot(
  filter(plot_data.EDM, Test == "delta_rho"),
  aes(
    x = T,
    y = value,
    group = model,
    linetype = model
  )
) +
  geom_line(aes(color = Test), size = 1.2) +
  geom_point(aes(color = Test), size = 2) +
  scale_color_manual(values = test_colors) +
  scale_linetype_manual(values = c("VAR" = "solid", "bistable" = "dashed")) +
  scale_y_continuous(limits = c(0, NA)) +   
  labs(
    title = "delta_rho",
    x = "T (time series length)",
    y = NULL,
    linetype = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.minor = element_line(color = "gray90", linetype = "dotted")
  ) +
  guides(color = "none")



# 7 --- SOME MORE DESCRIPTIVES ---

mean(agg_df$ratioSig[agg_df$Test == "KPSS.level"])
