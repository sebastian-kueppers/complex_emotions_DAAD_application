##############################
## - Meeting #04 -------------
##############################

setwd("C:/Users/Sebastian Küppers/Desktop/Formal Theory of Co-Occuring Emotions (DFG project)/_PhD/_PhD_Study_1/complexity_in_emotion_ESM_data/SIM_ESM")

# Import helper functions for data generation
source('HELPER_EDM_tests.R')
source('data_generation.R')

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


## --------------------------
# OBJECTIVE: Scale VAR time series relative to each other.
## --------------------------

# -- Simulation Parameters --

T_list <- c(25, 50, 75, 100, 150, 200)
N_list <- c(100)

# coef matrix
A_list <- list(
  A0.2 = matrix(c(
    0.2,  0.1, -0.05, -0.05,
    0.1,  0.2, -0.05, -0.05,
    -0.05, -0.05,  0.2, 0.1,
    -0.05, -0.05,  0.1, 0.2
  ), 4, 4, byrow = TRUE),
  
  A0.4 = matrix(c(
    0.4,  0.1, -0.05, -0.05,
    0.1,  0.4, -0.05, -0.05,
    -0.05, -0.05, 0.4,  0.1,
    -0.05, -0.05, 0.1,  0.4
  ), 4, 4, byrow = TRUE),
  
  A0.6 = matrix(c(
    0.6,  0.1, -0.05, -0.05,
    0.1,  0.6, -0.05, -0.05,
    -0.05, -0.05,  0.6, 0.1,
    -0.05, -0.05,  0.1, 0.6
  ), 4, 4, byrow = TRUE),
  
  A0.8 = matrix(c(
    0.8,  0.1, -0.05, -0.05,
    0.1,  0.8, -0.05, -0.05,
    -0.05, -0.05,  0.8, 0.1,
    -0.05, -0.05,  0.1, 0.8
  ), 4, 4, byrow = TRUE)
)

Sigmas <- c(0.1,0.5,1,1.5,2)


# Generate VAR data - ###

data_list.VAR.scaled.04 <- list()
set.seed(2310)

for (N in N_list) {
  data_list.VAR.scaled.04[[paste0("N", N)]] <- list()
  
  for (T in T_list) {
    data_list.VAR.scaled.04[[paste0("N", N)]][[paste0("T", T)]] <- list()
    
    for (A_name in names(A_list)) {          
      A <- A_list[[A_name]]
      
      for (Sig in Sigmas) {
        Sigma <- diag(4)
        diag(Sigma) <- Sig
        
        # --- Simulate raw VAR1 data ---
        raw_data <- simulate_var_subjects(N = N, T = T, A = A, Sigma = Sigma)
        
        # --- Save to data_list.VAR.scaled.04 ---
        data_list.VAR.scaled.04[[paste0("N", N)]][[paste0("T", T)]][[paste0("A", A_name)]][[paste0("Sigma", Sig)]] <- list(
          raw = raw_data
        )
        
        cat("Simulated: N =", N, "T =", T, "AR =", A_name, "Sigma =", Sig, "\n")
      }
    }
  }
}



## -----------------------------
# get overall min and max ------
# ------------------------------

get_all_numeric_vals <- function(x) {
  if (is.data.frame(x)) {
    # Extract numeric columns only
    unlist(x[sapply(x, is.numeric)], use.names = FALSE)
    
  } else if (is.list(x)) {
    # Recurse into lists
    unlist(lapply(x, get_all_numeric_vals), use.names = FALSE)
    
  } else {
    # Ignore everything else
    NULL
  }
}

vals <- get_all_numeric_vals(data_list.VAR.scaled.04)

# sanity checks
length(vals) # 4800000 - seems legit

lower <- quantile(vals, 0.01, na.rm = TRUE)  # 1st percentile
upper <- quantile(vals, 0.99, na.rm = TRUE)  # 99th percentile


## ----------------------------
# Scale data globally ---------
## ----------------------------

data_list.VAR.scaled.04 <- add_likert_next_to_raw_quantile(
  data_list.VAR.scaled.04,
  lower = lower,
  upper = upper,
  new_min = 1,
  new_max = 100
)

## ----------------------------
# Plot examples for PPT -------
## ----------------------------

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
      data_list.VAR.scaled.04$N100$T100[[A_name]][[sigma_name]]$likert_1_100[[10]]$V1
    )
    
    # convert to ts
    ts_var <- ts(ts_data[, 1])
    
    # plot
    plot(ts_var, type = "l", ylim = range(0,100), 
         # main = paste0(A_name, ", ", sigma_name),
         ylab = "", xlab = "")
  }
}

save(data_list.VAR.scaled.04,
     file = "04/data_list-VAR-scaled-04.RData")

# APPLY EDM COMPLEXITY TESTS TO VAR MODELS

VAR_EDM_results.1_100.04 <- list()

for (N in N_list) {
  N_name <- paste0("N", N)
  VAR_EDM_results.1_100.04[[N_name]] <- list()
  
  for (T in T_list) {
    T_name <- paste0("T", T)
    VAR_EDM_results.1_100.04[[N_name]][[T_name]] <- list()
    
    for (A in names(A_list)) {
      A_name <- paste0("A", A)
      VAR_EDM_results.1_100.04[[N_name]][[T_name]][[A_name]] <- list()
      
      for (S in Sigmas) {
        sigma_name <- paste0("Sigma", S)
        
        cat("Running EDM (parallel) for:",
            N_name, T_name, A_name, sigma_name, "\n")
        
        data_all <- data_list.VAR.scaled.04[[N_name]][[T_name]][[A_name]][[sigma_name]][["likert_1_100"]]
        
        cols <- colnames(data_all[[1]])
        
        ## -----------------------------
        ## PARALLEL over persons
        ## -----------------------------
        
        # Parallelized person loop
        person_results <- future_map(
          data_all,
          function(ts_data) {
            EDM_tests(ts_data,
                      EDM.include = c("E_opt", "SMap"))
          },
          .options = furrr_options(seed = TRUE)
        )
        
        VAR_EDM_results.1_100.04[[N_name]][[T_name]][[A_name]][[sigma_name]] <- person_results
      }  
    }
  }
}

save(VAR_EDM_results.1_100.04,
     file = "04/VAR_EDM_results-04.RData")


