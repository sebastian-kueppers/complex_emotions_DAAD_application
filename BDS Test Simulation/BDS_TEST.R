
#####################################
### BDS Test Simulation #############
#####################################

setwd("C:/Users/Sebastian Küppers/Desktop/Formal Theory of Co-Occuring Emotions (DFG project)/_PhD/_PhD_Study_1/Research Exchange/PROJECT/Preregistration/BDS Test Simulation")

# import packages
library(dplyr)
library(tseries)
library(furrr)
library(future)

# import helpers
source("data_generation.R")

## ----------------------------------
# SIMULATE DATA.AR ------------------
## ----------------------------------

## ----------------------------------
# Simulation Parameters -------------

C <- matrix(c(-.2,  .04, -.2, -.2,
              .04, -.2, -.2, -.2,
              -.2,  -.2, -.2,  .04,
              -.2,  -.2,  .04, -.2),
            nrow = 4, byrow = TRUE)

r <- c(1, 1, 1, 1)
mu_bistable <- 0
p <- 4

# Total simulation time: 2 weeks in minutes
no_minutes <- 14 * 24 * 60  # 20160 minutes


## ----------------------------------
# Parallelization Setup -------------

# Detect available cores and leave one free for the OS
n_cores <- parallel::detectCores() - 1
plan(multisession, workers = n_cores)
cat("Running on", n_cores, "cores\n")

## ----------------------------------
# Simulation ------------------------

simulate_ar1_grid_full <- function(
    N_list,
    T_list,
    phi_list  = c(0.2, 0.4, 0.6, 0.8),
    sigma_list = c(0.5, 1.0, 1.5, 2.0),
    burnin = 0,
    likert_scales = list(
      "likert_1_7"   = c(1, 7),
      "likert_1_100" = c(1, 100)
    ),
    seed    = NULL,
    verbose = TRUE
) {
  if (!is.null(seed)) set.seed(seed)
  
  data_list <- list()
  
  for (N in N_list) {
    n_key <- paste0("N", N)
    data_list[[n_key]] <- list()
    
    for (T in T_list) {
      t_key <- paste0("T", T)
      data_list[[n_key]][[t_key]] <- list()
      
      for (phi in phi_list) {
        phi_key <- paste0("Phi", phi)
        data_list[[n_key]][[t_key]][[phi_key]] <- list()
        
        for (sigma in sigma_list) {
          sigma_key <- paste0("Sigma", sigma)
          
          # Simulate raw data
          raw_data <- simulate_ar1_subjects(
            N      = N,
            T      = T,
            phi    = phi,
            sigma  = sigma,
            burnin = burnin
          )
          
          # Discretise to Likert scales
          likert_data <- list()
          for (nm in names(likert_scales)) {
            sc <- likert_scales[[nm]]
            likert_data[[nm]] <- discretize_likert(
              raw_data,
              scale_min = sc[1],
              scale_max = sc[2],
              suffix    = paste0("_", nm)
            )
          }
          
          # Store in nested list
          data_list[[n_key]][[t_key]][[phi_key]][[sigma_key]] <- c(
            list(raw = raw_data),
            likert_data
          )
          
          if (verbose) {
            cat("Simulated: N =", N, "| T =", T,
                "| Phi =", phi, "| Sigma =", sigma, "\n")
          }
        }
      }
    }
  }
  
  data_list
}

simulate_bistable_grid <- function(
    N         = 100,
    T_list    = seq(25, 150, by = 25),
    timestep  = 1,
    p         = 4,
    C         = C,
    mu        = mu_bistable,
    r         = r,
    noiseSD   = 1,
    seed      = NULL,
    verbose   = TRUE
) {
  if (!is.null(seed)) set.seed(seed)
  
  total_iterations <- length(T_list)
  pb   <- txtProgressBar(min = 0, max = total_iterations, style = 3)
  iter <- 0
  
  data_bistable <- list()
  
  for (T in T_list) {
    t_key <- paste0("T", T)
    
    subject_list <- vector("list", N)
    
    for (i in seq_len(N)) {
      
      traj <- simulate_bistable_trajectory(
        time     = T,
        timestep = timestep,
        p        = p,
        C        = C,
        mu       = mu,
        r        = r,
        noiseSD  = noiseSD,
        noise    = TRUE,
        pbar     = FALSE
      )
      
      subject_list[[i]] <- traj
    }
    
    data_bistable[[t_key]] <- subject_list
    
    iter <- iter + 1
    setTxtProgressBar(pb, iter)
  }
  
  close(pb)
  data_bistable
}

DATA.AR <- simulate_ar1_grid_full(
  N_list     = c(100),
  T_list     = c(25, 50, 75, 100, 125, 150, 500),
  phi_list   = c(0.2, 0.4, 0.6, 0.8),
  sigma_list = c(0.5, 1.0, 1.5, 2.0),
  burnin     = 50,
  seed       = 42
)

DATA.BISTABLE <- simulate_bistable_grid(
  N        = 100,
  T_list   = seq(25, 150, by = 25),
  timestep = 1,
  p        = p,
  C        = C,
  mu       = mu_bistable,
  r        = r,
  noiseSD  = 1,
  seed     = 42,
  verbose  = TRUE
)

saveRDS(DATA.AR, file = "DATA_AR.rds")
saveRDS(DATA.BISTABLE, file = "DATA_BISTABLE.rds")

## ----------------------------------
# APPLY PROCEDURE -------------------
## ----------------------------------

## ----------------------------------
# Create bootstrapping function -----

bds_bootstrap <- function(
    ar.res_z,
    m     = 2,
    eps   = 1,
    B     = 499,
    alpha = 0.05,
    seed  = NULL      
) {
  
  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)
  
  # Observed BDS statistic on actual residuals
  obs <- tryCatch(
    tseries::bds.test(ar.res_z, m = m, eps = eps)$statistic[[1]],
    error = function(e) NA
  )
  
  if (is.na(obs)) return(NA)
  
  # Bootstrap: shuffle residuals B times and recompute statistic
  boot_stats <- replicate(B, {
    shuffled <- sample(ar.res_z, replace = FALSE)
    tryCatch(
      tseries::bds.test(shuffled, m = m, eps = eps)$statistic[[1]],
      error = function(e) NA
    )
  })
  
  boot_stats <- na.omit(boot_stats)
  
  # Bootstrap p-value: proportion of null statistics
  # as or more extreme than observed (two-tailed)
  p_boot <- mean(abs(boot_stats) >= abs(obs))
  
  p_boot < alpha
}

## ----------------------------------
# Parallelization Helpers -----------

process_subject_ar <- function(subject_df) {
  ts       <- na.omit(subject_df$V1_likert_1_100)
  ar.fit   <- arima(ts, order = c(1, 0, 0))
  ar.res   <- na.omit(residuals(ar.fit))
  ar.res_z <- scale(ar.res)[, 1]
  bds_bootstrap(ar.res_z, m = 2, eps = 1, B = 499)
}

process_subject_bistable <- function(traj) {
  ts       <- na.omit(traj$x1)
  ar.fit   <- arima(ts, order = c(1, 0, 0))
  ar.res   <- na.omit(residuals(ar.fit))
  ar.res_z <- scale(ar.res)[, 1]
  bds_bootstrap(ar.res_z, m = 2, eps = 1, B = 499)
}

## ----------------------------------
# Apply BDS test: AR ----------------

# Count total iterations for progress bar
total_iterations <- length(names(DATA.AR)) *
  length(names(DATA.AR[[1]])) *
  length(names(DATA.AR[[1]][[1]])) *
  length(names(DATA.AR[[1]][[1]][[1]]))

# Initialise progress bar
pb <- txtProgressBar(
  min   = 0,
  max   = total_iterations,
  style = 3
)

iter <- 0

# Initialize timestamp
time_start_ar <- Sys.time()
cat("\n[", format(time_start_ar, "%Y-%m-%d %H:%M:%S"), "]",
    "Starting AR loop...\n")

RESULTS.AR <- list()

for (n_key in names(DATA.AR)) {
  RESULTS.AR[[n_key]] <- list()
  
  for (t_key in names(DATA.AR[[n_key]])) {
    RESULTS.AR[[n_key]][[t_key]] <- list()
    
    if (t_key == "T500") {
      iter <- iter + 1
      setTxtProgressBar(pb, iter)
      next
    }
    
    for (phi_key in names(DATA.AR[[n_key]][[t_key]])) {
      RESULTS.AR[[n_key]][[t_key]][[phi_key]] <- list()
      
      for (sigma_key in names(DATA.AR[[n_key]][[t_key]][[phi_key]])) {
        
        raw_data <- DATA.AR[[n_key]][[t_key]][[phi_key]][[sigma_key]]$likert_1_100
        
        # Parallelise over subjects within each condition
        condition_results <- furrr::future_map(
          raw_data,
          process_subject_ar,
          .options = furrr::furrr_options(seed = 2310)
        )
        
        RESULTS.AR[[n_key]][[t_key]][[phi_key]][[sigma_key]] <- condition_results
        
        iter <- iter + 1
        setTxtProgressBar(pb, iter)
      }
    }
  }
}

close(pb)

time_end_ar <- Sys.time()

cat("\n[", format(time_end_ar, "%Y-%m-%d %H:%M:%S"), "]",
    "AR loop complete. Time elapsed:",
    round(difftime(time_end_ar, time_start_ar, units = "mins"), 2),
    "minutes\n")

# Save
saveRDS(RESULTS.AR, file = "RESULTS_AR.rds")



## ----------------------------------
# Apply BDS test: BISTABLE-----------

# Count total iterations for progress bar
total_iterations <- length(names(DATA.BISTABLE)) * 100

# Initialise progress bar
pb <- txtProgressBar(
  min   = 0,
  max   = total_iterations,
  style = 3
)

iter <- 0

# Initialize timestamp
time_start_bistable <- Sys.time()
cat("\n[", format(time_start_bistable, "%Y-%m-%d %H:%M:%S"), "]",
    "Starting bistable loop...\n")

RESULTS.BISTABLE <- list()

for (t_key in names(DATA.BISTABLE)) {
  RESULTS.BISTABLE[[t_key]] <- list()
  
  subject_list <- DATA.BISTABLE[[t_key]]
  
  # Parallelise over subjects within each T condition
  condition_results <- furrr::future_map(
    subject_list,
    process_subject_bistable,
    .options = furrr::furrr_options(seed = 2310)
  )
  
  RESULTS.BISTABLE[[t_key]] <- condition_results
  
  iter <- iter + 1
  setTxtProgressBar(pb, iter)
}

close(pb)

time_end_bistable <- Sys.time()
cat("\n[", format(time_end_bistable, "%Y-%m-%d %H:%M:%S"), "]",
    "Bistable loop complete. Time elapsed:",
    round(difftime(time_end_bistable, time_start_bistable, units = "mins"), 2),
    "minutes\n")

saveRDS(RESULTS.BISTABLE, file = "RESULTS_BISTABLE.rds")


## ----------------------------------
# Reset to Sequential Processing ----
plan(sequential)

## ----------------------------------
# SUMMARIZE RESULTS -----------------
## ----------------------------------

summary_df <- expand.grid(
  N     = names(RESULTS.AR),
  T     = names(RESULTS.AR[[1]])[names(RESULTS.AR[[1]]) != "T500"],
  Phi   = names(RESULTS.AR[[1]][[1]]),
  Sigma = names(RESULTS.AR[[1]][[1]][[1]]),
  stringsAsFactors = FALSE
) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    rejection_rate = mean(
      unlist(RESULTS.AR[[N]][[T]][[Phi]][[Sigma]]),
      na.rm = TRUE
    )
  ) |>
  dplyr::ungroup()

print(summary_df, n = 100)

