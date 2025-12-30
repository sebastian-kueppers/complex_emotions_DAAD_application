
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
        # pb_person <- txtProgressBar(min = 0, max = N_persons, style = 3)
        
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
        
        # for (p in 1:N_persons) {
        #   ts_data <- lik_data[[p]]
        #   
        #   for (j in seq_along(cols)) {
        #     col <- cols[j]
        #     ts <- ts_data[[col]]
        #     
        #     # --- Bartels ---
        #     bartels_p <- tryCatch(bartels.rank.test(ts)$p.value, error = function(e) NA)
        #     res_array["Bartels", col, p] <- bartels_p
        #     
        #     # --- KPSS level ---
        #     kpss_level_p <- tryCatch(kpss.test(ts, null = "Level", lshort = TRUE)$p.value, error=function(e) NA)
        #     res_array["KPSS.level", col, p] <- kpss_level_p
        #     
        #     # --- KPSS trend ---
        #     kpss_trend_p <- tryCatch(kpss.test(ts, null = "Trend", lshort = TRUE)$p.value, error=function(e) NA)
        #     res_array["KPSS.trend", col, p] <- kpss_trend_p
        #     
        #     # --- Keenan ---
        #     keenan_p <- tryCatch(keenanTest(ts)$p.value, error=function(e) NA)
        #     res_array["Keenan", col, p] <- keenan_p
        #     
        #     # --- Tsay ---
        #     tsay_p <- tryCatch(tsayTest(ts)$p.value, error=function(e) NA)
        #     res_array["Tsay", col, p] <- tsay_p
        #     
        #     # --- TV-AR GAM ---
        #     N_ts <- length(ts)
        #     tt <- 1:(N_ts-1)
        #     tv <- tryCatch(
        #       gam(ts[2:N_ts] ~ s(tt, by=ts[1:(N_ts-1)], k=min(10,N_ts-2), bs="tp")),
        #       error=function(e) NULL
        #     )
        #     if (!is.null(tv)) {
        #       stvar <- summary(tv)
        #       if (!is.null(stvar$s.table) && nrow(stvar$s.table) >= 1) {
        #         res_array["TV-AR.p", col, p] <- round(stvar$s.table[1,"p-value"],4)
        #       }
        #     }
        #   }
        #   setTxtProgressBar(pb_person, p)
        # } # end person loop
        
        # close(pb_person)
        
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


set.seed(2310)

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
        
        # TEST_MODE <- TRUE   # <<< switch this later
        # 
        # if (TEST_MODE) {
        #   lik_data <- lik_data[seq_len(min(2, length(lik_data)))]
        #   lik_data <- lapply(lik_data, function(df) df[, 1, drop = FALSE])
        #   cols <- colnames(lik_data[[1]])
        #   N_persons <- length(lik_data)
        #   
        #   cat("⚠️ MF TEST MODE:",
        #       "persons =", N_persons,
        #       "| vars =", length(cols), "\n")
        # }
        
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
        person_results <- future_lapply(
          lik_data,
          run_mf_for_person,
          cols = cols,
          future.seed = TRUE
        )
        
        # Combine into matrix: vars × persons
        mf_array <- do.call(cbind, person_results)
        colnames(mf_array) <- paste0("P", seq_len(ncol(mf_array)))
        rownames(mf_array) <- cols
        
        # cat(
        #   "\nMF:", source_name,
        #   "| scale:", scale_name,
        #   "|", N_name,
        #   "|", T_name,
        #   "| persons:", N_persons, "\n"
        # )
        # pb <- txtProgressBar(min = 0, max = N_persons, style = 3)
        # 
        # for (p in seq_len(N_persons)) {
        #   ts_data <- lik_data[[p]]
        #   
        #   for (col in cols) {
        #     mf_array[col, p] <- multifractality_pval(ts_data[[col]])
        #   }
        #   
        #   setTxtProgressBar(pb, p)
        # }
        # close(pb)
        
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
  
  tv_edf <- numeric(length(cols))
  names(tv_edf) <- cols
  
  print(head(ts_data))
  
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
          tv_edf[col] <- round(stvar$s.table[1, "edf"], 3)
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
      res["PACF.n", col] <- length(sig_lags)
      res["PACF.max_lag", col] <- if(length(sig_lags) > 0) max(sig_lags) else NA
    }
  }
  
  list(raw = res, tv_edf = tv_edf)
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
        
        # ---- TEST MODE ----
        if (TEST_MODE) {
          lik_data <- lik_data[1:min(2, N_persons)]
          # lik_data <- lapply(lik_data, function(df) df[, 1, drop=FALSE])
          cols <- colnames(lik_data[[1]])
          N_persons <- length(lik_data)
          cat("⚠️ TEST MODE:", source_name, scale_name, N_name, T_name,
              "| Persons:", N_persons, "| Vars:", length(cols), "\n")
        }
        
        # ---- PARALLEL EXECUTION ----
        person_results <- future_lapply(
          lik_data,
          run_tv_cp_pacf_for_person,
          cols = cols,
          ALPHA_LEVEL = ALPHA_LEVEL,
          future.seed = TRUE
        )
        
        # Combine results
        res_array <- array(
          NA, 
          dim = c(length(non_significance_tests), length(cols), N_persons),
          dimnames = list(non_significance_tests, cols, paste0("P", 1:N_persons))
        )
        
        res_df <- matrix(NA, nrow = 1, ncol = length(cols), dimnames = list("TV-AR.EDF", cols))
        
        for (p in seq_len(N_persons)) {
          res_array[,,p] <- person_results[[p]]$raw
          res_df[1,] <- person_results[[p]]$tv_edf  # will overwrite, but can aggregate if needed
        }
        
        # Aggregate across persons
        sig_rate <- apply(res_array, c(1,2), function(x) mean(x < ALPHA_LEVEL, na.rm=TRUE))
        
        results_list[[N_name]][[T_name]] <- list(
          raw_array = res_array,
          sig_rate = sig_rate,
          tv_edf = res_df
        )
        
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

# for (source_name in names(data_sources)) {
#   
#   data_list <- data_sources[[source_name]]
#   
#   for (scale_name in names(scales)) {
#     
#     results_list <- list()
#     
#     for (N_name in names(data_list)) {
#       
#       results_list[[N_name]] <- list()
#       
#       for (T_name in names(data_list[[N_name]])) {
#         
#         lik_data <- data_list[[N_name]][[T_name]][[scales[[scale_name]]]]
#         N_persons <- length(lik_data)
#         if (N_persons == 0) next
#         
#         cols <- colnames(lik_data[[1]])
#         
#         # Initialize arrays to store results
#         res_array <- array(
#           NA,
#           dim = c(length(non_significance_tests), length(cols), N_persons),
#           dimnames = list(non_significance_tests, cols, paste0("P", 1:N_persons))
#         )
#         
#         res_df <- matrix(NA, nrow = 1, ncol = length(cols), dimnames = list("TV-AR.EDF", cols))
#         
#         cat("\nSource:", source_name, "| Scale:", scale_name, "| N:", N_name, "| T:", T_name, "| Persons:", N_persons, "\n")
#         pb <- txtProgressBar(min = 0, max = N_persons, style = 3)
#         
#         for (p in seq_len(N_persons)) {
#           ts_data <- lik_data[[p]]
#           
#           for (col in cols) {
#             ts <- ts_data[[col]]
#             N_ts <- length(ts)
#             
#             # --- TV-AR GAM ---
#             tt <- 1:(N_ts-1)
#             tv <- tryCatch(
#               gam(ts[2:N_ts] ~ s(tt, by=ts[1:(N_ts-1)], k=min(10,N_ts-2), bs="tp")),
#               error = function(e) NULL
#             )
#             if (!is.null(tv)) {
#               stvar <- summary(tv)
#               if (!is.null(stvar$s.table) && nrow(stvar$s.table) >= 1) {
#                 res_array["TV-AR.EDF", col, p] <- round(stvar$s.table[1, "edf"], 3)
#               }
#             }
#             
#             # --- Change Point Analysis ---
#             R <- ceiling(20 / ALPHA_LEVEL)
#             cp.out <- tryCatch(e.divisive(matrix(ts), R=R, sig.lvl=ALPHA_LEVEL), error=function(e) NULL)
#             if (!is.null(cp.out)) {
#               cp.n <- length(which(cp.out$p.values < ALPHA_LEVEL))
#               res_array["CP.n", col, p] <- cp.n
#             }
#             
#             # --- PACF analysis ---
#             pacf_vals <- tryCatch(pacf(ts, plot=FALSE)$acf, error=function(e) NA)
#             if (is.numeric(pacf_vals) && length(pacf_vals) > 1 && !all(is.na(pacf_vals))) {
#               thresh <- 2 / sqrt(N_ts)
#               sig_lags <- which(abs(pacf_vals) > thresh)
#               res_array["PACF.n", col, p] <- length(sig_lags)
#               res_array["PACF.max_lag", col, p] <- if(length(sig_lags) > 0) max(sig_lags) else NA
#             }
#             
#           } # end variable loop
#           
#           setTxtProgressBar(pb, p)
#         } # end person loop
#         
#         close(pb)
#         
#         # Aggregate across persons if needed
#         sig_rate <- apply(res_array, c(1,2), function(x) mean(x < ALPHA_LEVEL, na.rm=TRUE))
#         
#         results_list[[N_name]][[T_name]] <- list(
#           raw_array = res_array,
#           sig_rate = sig_rate,
#           tv_edf = res_df
#         )
#         
#         cat("TV/CP/PACF DONE:", source_name, scale_name, N_name, T_name, "\n")
#         
#       } # end T loop
#     } # end N loop
#     
#     # Assign results to proper output object
#     if (scale_name == "1_7") {
#       if (source_name == "VAR") tv_cp_pacf_results_1_7.VAR <- results_list
#       if (source_name == "bistable") tv_cp_pacf_results_1_7.bistable <- results_list
#     } else if (scale_name == "1_100") {
#       if (source_name == "VAR") tv_cp_pacf_results_1_100.VAR <- results_list
#       if (source_name == "bistable") tv_cp_pacf_results_1_100.bistable <- results_list
#     }
#     
#   } # end scale loop
# } # end data source loop


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
        pb <- txtProgressBar(min = 0, max = N_persons, style = 3)
        
        for (p in seq_len(N_persons)) {
          ts_data <- lik_data[[p]]
          
          for (col in cols) {
            ts <- ts_data[[col]]
            N <- length(ts)
            
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
                res_array["E_opt", col, p] <- e
                
                ## --- S-map nonlinearity test ---
                thetas <- seq(0, 8, by = 0.5)
                rho_theta <- rep(NA, length(thetas))
                
                for (k in seq_along(thetas)) {
                  sm <- tryCatch(
                    suppressWarnings(
                      SMap(
                        dataFrame = as.data.frame(ts),
                        lib = lib,
                        pred = pred,
                        columns = "ts",
                        target = "ts",
                        E = e,
                        tau = tau,
                        Tp = 1,
                        theta = thetas[k],
                        noTime = TRUE,
                        silent = TRUE
                      )
                    ),
                    error = function(e) NULL
                  )
                  
                  if (!is.null(sm)) {
                    rho_theta[k] <- sm$rho
                  }
                }
                
                if (any(!is.na(rho_theta))) {
                  res_array["rho_theta0", col, p] <- rho_theta[thetas == 0]
                  res_array["rho_theta_opt", col, p] <- max(rho_theta, na.rm = TRUE)
                  res_array["theta_opt", col, p] <- thetas[which.max(rho_theta)]
                  res_array["delta_rho", col, p] <- 
                    res_array["rho_theta_opt", col, p] - 
                    res_array["rho_theta0", col, p]
                }
                
                pi <- tryCatch(suppressWarnings(
                  PredictInterval(dataFrame = as.data.frame(ts), noTime = TRUE, columns = "ts", target = "ts", 
                                  E = e, tau = tau, lib = lib, pred = pred, maxTp = 20, showPlot = FALSE)
                ), error = function(e) NULL)
                
                if (!is.null(pi) && nrow(pi) >= 5) {
                  decay <- tryCatch({
                    fit <- lm(rho ~ Tp, data = pi[1:5,])
                    coef(fit)[2] * 5
                  }, error = function(e) NA)
                  res_array["pred_decay", col, p] <- decay
                }
              }
            }
            
          } # end variable loop
          
          setTxtProgressBar(pb, p)
        } # end person loop
        
        close(pb)
        
        # Aggregate across persons if needed
        sig_rate <- apply(res_array, c(1,2), function(x) mean(x < ALPHA_LEVEL, na.rm=TRUE))
        
        results_list[[N_name]][[T_name]] <- list(
          raw_array = res_array,
          sig_rate = sig_rate,
          tv_edf = res_df
        )
        
        cat("TV/CP/PACF DONE:", source_name, scale_name, N_name, T_name, "\n")
        
      } # end T loop
    } # end N loop
    
    # Assign results to proper output object
    if (scale_name == "1_7") {
      if (source_name == "VAR") tv_cp_pacf_results_1_7.VAR <- results_list
      if (source_name == "bistable") tv_cp_pacf_results_1_7.bistable <- results_list
    } else if (scale_name == "1_100") {
      if (source_name == "VAR") tv_cp_pacf_results_1_100.VAR <- results_list
      if (source_name == "bistable") tv_cp_pacf_results_1_100.bistable <- results_list
    }
    
  } # end scale loop
} # end data source loop

### 5 --- AGGREGATE RESULTS --- ###

# Initialize a list to store aggregated results
agg_list <- list()

# Loop over scales: 1-7 and 1-100
for (scale_name in c("1_7", "1_100")) {
  
  # Select the corresponding results list
  results_scale <- if (scale_name == "1_7") test_results_1_7 else test_results_1_100
  
  # Loop over sample sizes (N)
  for (N_name in names(results_scale)) {
    
    # Loop over time series lengths (T)
    for (T_name in names(results_scale[[N_name]])) {
      
      res_array <- results_scale[[N_name]][[T_name]]$raw_array  # get raw p-values array
      
      # Initialize vector to store proportion of significant p-values per test
      ratio_sig <- numeric(dim(res_array)[1])
      names(ratio_sig) <- dimnames(res_array)[[1]]
      
      # Aggregate over all variables × persons
      for (test_name in dimnames(res_array)[[1]]) {
        vals <- as.vector(res_array[test_name,,])  # flatten variable × person values
        ratio_sig[test_name] <- mean(vals < ALPHA_LEVEL, na.rm = TRUE)  # ratio (0-1)
      }
      
      # Store as a data frame
      agg_list[[paste0(scale_name, "_", N_name, "_", T_name)]] <- data.frame(
        Test = names(ratio_sig),
        Likert = scale_name,
        N = N_name,
        T = T_name,
        ratioSig = ratio_sig,  # ratio is the last column
        row.names = NULL
      )
      
    } # end T loop
  } # end N loop
} # end scale loop

# Combine all data frames into one tidy data frame
agg_df <- bind_rows(agg_list)

# SHow significance_tests of agg_df
agg_df[agg_df$Test %in% significance_tests,]


### 6 --- PLOT --- ###
plot_data.N200.1_100 <- agg_df %>%
  filter(
    N == "N200",
    Likert == "1_100",
    Test %in% significance_tests,
    Test != "multifractality"      # ← remove multifractality
  ) %>%
  mutate(T = factor(T, levels = paste0("T", sort(T_list))))

plot_data.N200.1_7 <- agg_df %>%
  filter(
    N == "N200",
    Likert == "1_7",
    Test %in% significance_tests,
    Test != "multifractality"      # ← remove multifractality
  ) %>%
  mutate(T = factor(T, levels = paste0("T", sort(T_list))))

# Create the line plot
ggplot(plot_data.N200.1_100, aes(x = T, y = ratioSig, color = Test, group = Test)) +
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
ggplot(plot_data.N200.1_7, aes(x = T, y = ratioSig, color = Test, group = Test)) +
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


# 7 --- SOME MORE DESCRIPTIVES ---

mean(agg_df$ratioSig[agg_df$Test == "KPSS.level"])
