
### 0 --- PACKAGES AND HELPER FUNCTIONS --- ###

setwd("C:/Users/Sebastian Küppers/Desktop/Formal Theory of Co-Occuring Emotions (DFG project)/_PhD/_PhD_Study_1/complexity_in_emotion_ESM_data/SIM_ESM")

# To gain access to bistable model data, import previously simulated data
load("C:/Users/Sebastian Küppers/Desktop/Formal Theory of Co-Occuring Emotions (DFG project)/_PhD/_PhD_Study_1/complexity_in_emotion_ESM_data/SIM_ESM/simulated_data.RData")

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
library(vars)
library(tibble)

# set up parallelization
n_cores <- parallel::detectCores() - 1
plan(multisession, workers = 2)


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
A_list <- list(
  A1 = matrix(c(
    0.2,  0.1, 0, 0,
    0.1,  0.2, 0, 0,
    0,    0,  0.2, 0.1,
    0,    0,  0.1, 0.2
  ), 4, 4, byrow = TRUE),
  
  A2 = matrix(c(
    0.4,  0.1, -0.05, -0.05,
    0.1,  0.4, -0.05, -0.05,
    -0.05, -0.05, 0.4,  0.1,
    -0.05, -0.05, 0.1,  0.4
  ), 4, 4, byrow = TRUE),
  
  A3 = matrix(c(
    0.6,  0.1, 0, 0,
    0.1,  0.6, 0, 0,
    0,    0,  0.6, 0.1,
    0,    0,  0.1, 0.6
  ), 4, 4, byrow = TRUE),
  
  A4 = matrix(c(
    0.8,  0.1, 0, 0,
    0.1,  0.8, 0, 0,
    0,    0,  0.8, 0.1,
    0,    0,  0.1, 0.8
  ), 4, 4, byrow = TRUE)
)


Sigmas <- c(0.1,0.5,1,1.5,2)

# Parameter settings lead to model where approx. 31% of all variance is 
# explained (check by fitting a VAR)

### 3 --- SIMULATION LOOP USING HELPER FUNCTIONS --- ###

# 3.1 - VAR(1) MODEL - ###

data_list.VAR <- list()
set.seed(2310)

for (N in N_list) {
  data_list.VAR[[paste0("N", N)]] <- list()
  
  for (T in T_list) {
    data_list.VAR[[paste0("N", N)]][[paste0("T", T)]] <- list()
    
    for (A_name in names(A_list)) {          
      A <- A_list[[A_name]]
      
      for (Sig in Sigmas) {
        Sigma <- diag(4)
        diag(Sigma) <- Sig
        
        # --- Simulate raw VAR1 data ---
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
        data_list.VAR[[paste0("N", N)]][[paste0("T", T)]][[paste0("A", A_name)]][[paste0("Sigma", Sig)]] <- list(
          raw = raw_data,
          likert_1_7 = likert_1_7,
          likert_1_100 = likert_1_100
        )
        
        cat("Simulated: N =", N, "T =", T, "AR =", A_name, "Sigma =", Sig, "\n")
      }
    }
  }
}

# plot examples for PPT
Sigma_levels <- c(0.1, 1, 2)
A_levels <- names(A_list) 

# 3 rows (Sigma), 4 columns (A)
par(mfrow = c(3, length(A_levels)), mar = c(2, 2, 2, 1))

for (sigma in Sigma_levels) {
  sigma_name <- paste0("Sigma", sigma)
  
  for (A in A_levels) {
    A_name <- paste0("A", A)
    # pick first simulated subject for simplicity
    ts_data <- as.data.frame(
      data_list.VAR$N100$T100[[A_name]][[sigma_name]]$raw[[1]]
    )
    
    # convert to ts
    ts_var <- ts(ts_data[, 1])
    
    # plot
    plot(ts_var, type = "l", ylim = range(-4.5,4.5), 
         # main = paste0(A_name, ", ", sigma_name),
         ylab = "", xlab = "")
  }
}

# FIT VARS - LIKERT

# Initialize results
R2_results <- list(
  likert_1_100 = list(),
  raw          = list()
)

set.seed(2310)

for (A in names(A_list)) {
  A_name <- paste0("A", A)
  cat("=== Starting AR matrix:", A_name, "===\n")
  
  for (s in Sigmas) {
    sigma_name <- paste0("Sigma", s)
    cat("  Sigma:", sigma_name, "\n")
    
    # Access the simulated data
    sim_list_likert <- data_list.VAR$N100$T100[[A_name]][[sigma_name]]$likert_1_100
    sim_list_raw    <- data_list.VAR$N100$T100[[A_name]][[sigma_name]]$raw
    
    # Prepare matrices for R²
    R2_mat_likert <- matrix(NA, nrow = 100, ncol = 4)
    R2_mat_raw    <- matrix(NA, nrow = 100, ncol = 4)
    colnames(R2_mat_likert) <- colnames(R2_mat_raw) <- c("V1","V2","V3","V4")
    
    for (i in 1:100) {
      # --- Likert ---
      dat_i <- sim_list_likert[[i]]
      fit <- tryCatch(VAR(dat_i, p = 1, type = "const"), error = function(e) NULL)
      if (!is.null(fit)) R2_mat_likert[i, ] <- sapply(summary(fit)$varresult, function(x) x$r.squared)
      
      # --- Raw ---
      dat_i <- sim_list_raw[[i]]
      fit <- tryCatch(VAR(dat_i, p = 1, type = "const"), error = function(e) NULL)
      if (!is.null(fit)) R2_mat_raw[i, ] <- sapply(summary(fit)$varresult, function(x) x$r.squared)
    }
    
    # Save results
    if (is.null(R2_results$likert_1_100[[sigma_name]])) R2_results$likert_1_100[[sigma_name]] <- list()
    if (is.null(R2_results$raw[[sigma_name]])) R2_results$raw[[sigma_name]] <- list()
    
    R2_results$likert_1_100[[sigma_name]][[A_name]] <- R2_mat_likert
    R2_results$raw[[sigma_name]][[A_name]]          <- R2_mat_raw
  }
  
  cat("=== Finished AR matrix:", A_name, "===\n\n")
}

# summarize all R^2
summary_table_R2.raw <- matrix(NA, nrow = length(Sigmas), ncol = length(A_list),
                               dimnames = list(paste0("Sigma", Sigmas), paste0("A", names(A_list))))

summary_table_R2.lik <- matrix(NA, nrow = length(Sigmas), ncol = length(A_list),
                        dimnames = list(paste0("Sigma", Sigmas), paste0("A", names(A_list))))

for (A in names(A_list)) {
  A_name <- paste0("A",A)
  for (s in Sigmas) {
    sigma_name <- paste0("Sigma", s)
    
    R2_mat.raw <- R2_results$raw[[sigma_name]][[A_name]]
    summary_table_R2.raw[sigma_name, A_name] <- mean(colMeans(R2_mat.raw, na.rm = TRUE))
    
    R2_mat.lik <- R2_results$likert_1_100[[sigma_name]][[A_name]]
    summary_table_R2.lik[sigma_name, A_name] <- mean(colMeans(R2_mat.lik, na.rm = TRUE))
  }
}

summary_table_R2.raw
summary_table_R2.lik


# APPLY EDM COMPLEXITY TESTS TO VAR MODELS

VAR_EDM_results.1_100 <- list()

for (N in N_list) {
  N_name <- paste0("N", N)
  VAR_EDM_results.1_100[[N_name]] <- list()
  
  for (T in T_list) {
    T_name <- paste0("T", T)
    VAR_EDM_results.1_100[[N_name]][[T_name]] <- list()
    
    for (A in names(A_list)) {
      A_name <- paste0("A", A)
      VAR_EDM_results.1_100[[N_name]][[T_name]][[A_name]] <- list()
      
      for (S in Sigmas) {
        sigma_name <- paste0("Sigma", S)
        
        cat("Running EDM (parallel) for:",
            N_name, T_name, A_name, sigma_name, "\n")
        
        data_all <- data_list.VAR[[N_name]][[T_name]][[A_name]][[sigma_name]][["likert_1_100"]]
        
        cols <- colnames(data_all[[1]])
        
        ## -----------------------------
        ## PARALLEL over persons
        ## -----------------------------
        
        # Parallelized person loop
        person_results <- future_map(
          data_all,
          function(ts_data) {
            res_person <- matrix(NA, nrow = length(edm_tests), ncol = length(cols),
                                 dimnames = list(edm_tests, cols))
            
            for (col in cols) {
              ts <- ts_data[[col]]
              N <- length(ts)
              
              # tau <- tryCatch(
              #   suppressWarnings(timeLag(unlist(ts), technique = "ami",
              #                            selection.method = "first.minimum", lag.max = 10, do.plot = FALSE)),
              #   error = function(e) NA
              # )
              
              tau <- -1
              
              # FALLBACK
              # if (is.na(tau)) tau <- 1 
              
              if (!is.na(tau)) {
                mid <- floor(N * 0.7)
                lib <- c(1, mid)
                pred <- c(mid + 1, N)
                
                # maximum embed dimension to test depends on data length
                Tp <- 1
                
                e <- tryCatch({
                  emb_out <- suppressWarnings(EmbedDimension(
                    dataFrame = as.data.frame(ts),
                    lib = lib,
                    pred = pred,
                    maxE = 15,
                    Tp = 1,
                    # tau = tau,
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
                  # e <- 1

                  for (k in seq_along(thetas)) {
                    sm <- tryCatch(
                      suppressWarnings(SMap(
                        dataFrame = as.data.frame(ts),
                        lib = lib,
                        pred = pred,
                        columns = "ts",
                        target = "ts",
                        E = e,
                        # tau = -1,
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
                    # tau = tau,
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
        
        VAR_EDM_results.1_100[[N_name]][[T_name]][[A_name]][[sigma_name]] <- person_results
        
        cat("DONE:", N_name, T_name, A_name, sigma_name, "\n")
      }
    }
  }
}

save(VAR_EDM_results.1_100,
     file = "VAR_EDM_results_02-02.RData")


## AGGREGATE INTO DATAFRAME

out <- list()

for (N_name in names(VAR_EDM_results.1_100)) {
  N_list <- VAR_EDM_results.1_100[[N_name]]
  
  for (T_name in names(N_list)) {
    T_list <- N_list[[T_name]]
    
    for (A_name in names(T_list)) {
      A_list <- T_list[[A_name]]
      
      for (Sigma_name in names(A_list)) {
        Sigma_list <- A_list[[Sigma_name]]
        
        # Collect all 100 dfs into long format
        long_df <- do.call(rbind, lapply(Sigma_list, function(df) {
          df <- as.data.frame(df)
          df %>%
            rownames_to_column("Metric") %>%
            pivot_longer(
              cols = -Metric,
              names_to = "Variable",
              values_to = "Value"
            )
        }))
        
        summary_df <- long_df %>%
          group_by(Metric) %>%
          summarise(
            Value = mean(Value, na.rm = TRUE),
            .groups = "drop"
          )
        
        summary_df$N     <- N_name
        summary_df$T     <- T_name
        summary_df$A     <- A_name
        summary_df$Sigma <- Sigma_name
        
        out[[length(out) + 1]] <- summary_df
      }
    }
  }
}

VAR_EDM.final_df <- bind_rows(out)
VAR_EDM.final_df <- VAR_EDM.final_df[, c("N", "T", "A", "Sigma", "Metric", "Value")]

# Rename A (Autoregressive effects)
VAR_EDM.final_df$A[VAR_EDM.final_df$A == "AA1"] <- "A0.2"
VAR_EDM.final_df$A[VAR_EDM.final_df$A == "AA2"] <- "A0.4"
VAR_EDM.final_df$A[VAR_EDM.final_df$A == "AA3"] <- "A0.6"
VAR_EDM.final_df$A[VAR_EDM.final_df$A == "AA4"] <- "A0.8"

# load latest EDM results from server simulation
load("C:/Users/Sebastian Küppers/Desktop/Formal Theory of Co-Occuring Emotions (DFG project)/_PhD/_PhD_Study_1/complexity_in_emotion_ESM_data/SIM_ESM/EDM_test_results_02-02.RData")

# Aggregate bistable results

out_bistable <- list()

for (N_name in names(edm_results_1_100.bistable)) {
  N_list <- edm_results_1_100.bistable[[N_name]]
  
  for (T_name in names(N_list)) {
    T_list <- N_list[[T_name]]
    
    long_df <- do.call(rbind, lapply(T_list, function(df) {
      df <- as.data.frame(df)
      df %>%
        rownames_to_column("Metric") %>%
        pivot_longer(
          cols = -Metric,
          names_to = "Variable",
          values_to = "Value"
        )
    }))
    
    summary_df <- long_df %>%
      group_by(Metric) %>%
      summarise(
        Value = mean(Value, na.rm = TRUE),
        .groups = "drop"
      )
    
    summary_df$N <- N_name
    summary_df$T <- T_name
    
    out_bistable[[length(out_bistable) + 1]] <- summary_df
  }
}

BISTABLE_EDM_0202.final_df <- bind_rows(out_bistable)
BISTABLE_EDM_0202.final_df <- BISTABLE_EDM_0202.final_df[, c("N", "T", "Metric", "Value")]


# Aggregate VAR results from Feb 02

out_var <- list()

for (N_name in names(edm_results_1_100.VAR)) {
  N_list <- edm_results_1_100.VAR[[N_name]]
  
  for (T_name in names(N_list)) {
    T_list <- N_list[[T_name]]
    
    long_df <- do.call(rbind, lapply(T_list, function(df) {
      df <- as.data.frame(df)
      df %>%
        rownames_to_column("Metric") %>%
        pivot_longer(
          cols = -Metric,
          names_to = "Variable",
          values_to = "Value"
        )
    }))
    
    summary_df <- long_df %>%
      group_by(Metric) %>%
      summarise(
        Value = mean(Value, na.rm = TRUE),
        .groups = "drop"
      )
    
    summary_df$N <- N_name
    summary_df$T <- T_name
    
    out_var[[length(out_var) + 1]] <- summary_df
  }
}

VAR_EDM_0202.final_df <- bind_rows(out_var)
VAR_EDM_0202.final_df <- VAR_EDM_0202.final_df[, c("N", "T", "Metric", "Value")]


# PLOT

T_levels <- VAR_EDM.final_df %>%
  distinct(T) %>%
  mutate(T_num = as.numeric(sub("T", "", T))) %>%
  arrange(T_num) %>%
  pull(T)

