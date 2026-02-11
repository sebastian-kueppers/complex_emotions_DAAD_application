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

saveRDS(data_list.VAR.scaled.04,
     file = "04/data/data_list-VAR-scaled-04.rds")

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

# VAR

out <- list()

for (N_name in names(VAR_EDM_results.1_100.04)) {
  N_list <- VAR_EDM_results.1_100.04[[N_name]]
  
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
            mean  = mean(Value, na.rm = TRUE),
            sd    = sd(Value, na.rm = TRUE),
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
VAR_EDM.final_df <- VAR_EDM.final_df[, c("N", "T", "A", "Sigma", "Metric", "mean", "sd")]

VAR_EDM.final_df$A[VAR_EDM.final_df$A == "AA1"] <- "A0.2"
VAR_EDM.final_df$A[VAR_EDM.final_df$A == "AA2"] <- "A0.4"
VAR_EDM.final_df$A[VAR_EDM.final_df$A == "AA3"] <- "A0.6"
VAR_EDM.final_df$A[VAR_EDM.final_df$A == "AA4"] <- "A0.8"

saveRDS(VAR_EDM.final_df,
     file = "04/results/VAR_EDM-df.rds")


## -----------------------------
# DATA VISUALIZATION -----------
## -----------------------------

# Load bistable and Lorentz data
results.BISTABLE <- readRDS("04/results/BISTABLE_EDM-df.rds")
results.LORENZ <- readRDS("04/results/LORENZ_EDM-df.rds")
results.VAR <- readRDS("04/results/VAR_EDM-df.rds")

# -- Organize into nice dataframes
sigma_panels <- c("Sigma0.1", "Sigma1", "Sigma2")

# Make T into levels
T_levels <- results.VAR %>%
  distinct(T) %>%
  mutate(T_num = as.numeric(sub("T", "", T))) %>%
  arrange(T_num) %>%
  pull(T)

## ------------------------------
# E_OPT -------------------------

# Bistable
df.bistable.E_opt <- results.BISTABLE %>%
  filter(
    Metric == "E_opt"
  ) %>%
  mutate(
    T = factor(T, levels = T_levels),
    A = "BISTABLE",   
    Sigma = "BISTABLE",
    Source = "BISTABLE"
  )

df.bistable.E_opt.bySigma <- df.bistable.E_opt %>%
  select(-Sigma) %>%              # drop BISTABLE sigma
  tidyr::crossing(
    Sigma = sigma_panels          # replicate for each panel
  )


# Lorenz 

df.Lorenz.E_opt <- results.LORENZ %>%
  filter(
    Metric == "E_opt",
  ) %>%
  mutate(
    T = factor(T, levels = T_levels))

df.Lorenz.E_opt <- df.Lorenz.E_opt %>%
  as_tibble() %>%                    # convert to tibble
  mutate(
    T = factor(T, levels = T_levels),  # make T a factor
    N = as.character(N),
    Sigma = as.character(Sigma),
    A = as.character(A),
    Source = as.character(Source)
  )

df.Lorenz.E_opt.bySigma <- df.Lorenz.E_opt %>%
  select(-Sigma) %>%              # drop BISTABLE sigma
  tidyr::crossing(
    Sigma = sigma_panels          # replicate for each panel
  )

# VAR data
df.var.E_opt.bySigma <- results.VAR %>%
  filter(
    Metric == "E_opt",
    Sigma %in% sigma_panels
  ) %>%
  mutate(
    T = factor(T, levels = T_levels),
    Source = "VAR"
  )

# bind into one df
df.all.E_opt.bySigma <- bind_rows(
  df.var.E_opt.bySigma,
  df.bistable.E_opt.bySigma,
  df.Lorenz.E_opt.bySigma
)

df.all.E_opt.bySigma$A[df.all.E_opt.bySigma$A == "AA0.2"] <- "A0.2"
df.all.E_opt.bySigma$A[df.all.E_opt.bySigma$A == "AA0.4"] <- "A0.4"
df.all.E_opt.bySigma$A[df.all.E_opt.bySigma$A == "AA0.6"] <- "A0.6"
df.all.E_opt.bySigma$A[df.all.E_opt.bySigma$A == "AA0.8"] <- "A0.8"


# -- PLOT --

# PLOT
ggplot(
  df.all.E_opt.bySigma,
  aes(x = T, y = mean, group = A)
) +
  
  geom_line(
    data = subset(df.all.E_opt.bySigma, Source == "VAR"),
    aes(color = A, group = A),
    size = 1
  ) +
  
  geom_line(
    data = subset(df.all.E_opt.bySigma, Source == "BISTABLE"),
    aes(color = "BISTABLE", group = Source),
    size = 1.2
  ) +
  
  geom_line(
    data = subset(df.all.E_opt.bySigma, Source == "LORENZ"),
    aes(color = "LORENZ", group = Source),
    size = 1.2
  ) +
  
  facet_wrap(~ Sigma, nrow = 1) +
  
  scale_color_manual(
    name = "Dataset",
    values = c(
      "A0.2"     = "#c6dbef",
      "A0.4"     = "#9ecae1",
      "A0.6"     = "#6baed6",
      "A0.8"     = "#3182bd",
      "BISTABLE" = "orange",
      "LORENZ"   = "green"
    ),
    breaks = c("A0.2", "A0.4", "A0.6", "A0.8", "BISTABLE", "LORENZ"),
    labels = c("A = 0.2", "A = 0.4", "A = 0.6", "A = 0.8", "BISTABLE", "LORENZ")
  ) +
  
  coord_cartesian(ylim = c(0, 8)) +
  labs(
    y = "Optimal Embedding Dimension",
    x = "Time Series Length"
  ) +
  theme_minimal()


## ----------------------------------
# THETA_OPT -------------------------

# Bistable
df.bistable.theta_opt <- results.BISTABLE %>%
  filter(
    Metric == "theta_opt"
  ) %>%
  mutate(
    T = factor(T, levels = T_levels),
    A = "BISTABLE",   
    Sigma = "BISTABLE",
    Source = "BISTABLE"
  )

df.bistable.theta_opt.bySigma <- df.bistable.theta_opt %>%
  select(-Sigma) %>%              # drop BISTABLE sigma
  tidyr::crossing(
    Sigma = sigma_panels          # replicate for each panel
  )


# Lorenz 

df.Lorenz.theta_opt <- results.LORENZ %>%
  filter(
    Metric == "theta_opt"
  ) %>%
  mutate(
    T = factor(T, levels = T_levels))

df.Lorenz.theta_opt <- df.Lorenz.theta_opt %>%
  as_tibble() %>%                    # convert to tibble
  mutate(
    T = factor(T, levels = T_levels),  # make T a factor
    N = as.character(N),
    Sigma = as.character(Sigma),
    A = as.character(A),
    Source = as.character(Source)
  )

df.Lorenz.theta_opt.bySigma <- df.Lorenz.theta_opt %>%
  select(-Sigma) %>%              # drop BISTABLE sigma
  tidyr::crossing(
    Sigma = sigma_panels          # replicate for each panel
  )

# VAR data
df.var.theta_opt.bySigma <- results.VAR %>%
  filter(
    Metric == "theta_opt",
    Sigma %in% sigma_panels
  ) %>%
  mutate(
    T = factor(T, levels = T_levels),
    Source = "VAR"
  )

# bind into one df
df.all.theta_opt.bySigma <- bind_rows(
  df.var.theta_opt.bySigma,
  df.bistable.theta_opt.bySigma,
  df.Lorenz.theta_opt.bySigma
)

df.all.theta_opt.bySigma$A[df.all.theta_opt.bySigma$A == "AA0.2"] <- "A0.2"
df.all.theta_opt.bySigma$A[df.all.theta_opt.bySigma$A == "AA0.4"] <- "A0.4"
df.all.theta_opt.bySigma$A[df.all.theta_opt.bySigma$A == "AA0.6"] <- "A0.6"
df.all.theta_opt.bySigma$A[df.all.theta_opt.bySigma$A == "AA0.8"] <- "A0.8"

# PLOT
ggplot(
  df.all.theta_opt.bySigma,
  aes(x = T, y = mean, group = A)
) +
  
  geom_line(
    data = subset(df.all.theta_opt.bySigma, Source == "VAR"),
    aes(color = A, group = A),
    size = 1
  ) +
  
  geom_line(
    data = subset(df.all.theta_opt.bySigma, Source == "BISTABLE"),
    aes(color = "BISTABLE", group = Source),
    size = 1.2
  ) +
  
  geom_line(
    data = subset(df.all.theta_opt.bySigma, Source == "LORENZ"),
    aes(color = "LORENZ", group = Source),
    size = 1.2
  ) +
  
  facet_wrap(~ Sigma, nrow = 1) +
  
  scale_color_manual(
    name = "Dataset",
    values = c(
      "A0.2"     = "#c6dbef",
      "A0.4"     = "#9ecae1",
      "A0.6"     = "#6baed6",
      "A0.8"     = "#3182bd",
      "BISTABLE" = "orange",
      "LORENZ"   = "green"
    ),
    breaks = c("A0.2", "A0.4", "A0.6", "A0.8", "BISTABLE", "LORENZ"),
    labels = c("A = 0.2", "A = 0.4", "A = 0.6", "A = 0.8", "BISTABLE", "LORENZ")
  ) +
  
  coord_cartesian(ylim = c(0, 8)) +
  labs(
    y = "Optimal Theta (= Nonlinearity)",
    x = "Time Series Length"
  ) +
  theme_minimal()



## ------------------------
# EXTRA LONG TIME SERIES --
## ------------------------

# Question: Maybe 200 observations are still to few to see difference between 
# the different data generation models. 
# Idea: Simulate very long time series to be sure.

# Set simulatiob parameters
T_list.long <- c(500, 1000, 2000)
N_list.long <- c(100)

## ------------------------
# VAR MODEL ---------------

# coef matrix
A_list.long <- list(
  A0.2 = matrix(c(
    0.2,  0.1, 0, 0,
    0.1,  0.2, 0, 0,
    0,    0,  0.2, 0.1,
    0,    0,  0.1, 0.2
  ), 4, 4, byrow = TRUE),
  
  A0.4 = matrix(c(
    0.4,  0.1, -0.05, -0.05,
    0.1,  0.4, -0.05, -0.05,
    -0.05, -0.05, 0.4,  0.1,
    -0.05, -0.05, 0.1,  0.4
  ), 4, 4, byrow = TRUE),
  
  A0.6 = matrix(c(
    0.6,  0.1, 0, 0,
    0.1,  0.6, 0, 0,
    0,    0,  0.6, 0.1,
    0,    0,  0.1, 0.6
  ), 4, 4, byrow = TRUE),
  
  A0.8 = matrix(c(
    0.8,  0.1, 0, 0,
    0.1,  0.8, 0, 0,
    0,    0,  0.8, 0.1,
    0,    0,  0.1, 0.8
  ), 4, 4, byrow = TRUE)
)

Sigmas.long <- c(1)

data_list.VAR.scaled.long.04 <- list()

set.seed(2310)

for (N in N_list.long) {
  data_list.VAR.scaled.long.04[[paste0("N", N)]] <- list()
  
  for (T in T_list.long) {
    data_list.VAR.scaled.long.04[[paste0("N", N)]][[paste0("T", T)]] <- list()
    
    for (A_name in names(A_list.long)) {          
      A <- A_list[[A_name]]
      
      for (Sig in Sigmas.long) {
        Sigma <- diag(4)
        diag(Sigma) <- Sig
        
        # --- Simulate raw VAR1 data ---
        raw_data <- simulate_var_subjects(N = N, T = T, A = A, Sigma = Sigma)
        
        # --- Save to data_list.VAR.scaled.04 ---
        data_list.VAR.scaled.long.04[[paste0("N", N)]][[paste0("T", T)]][[paste0("A", A_name)]][[paste0("Sigma", Sig)]] <- list(
          raw = raw_data
        )
        
        cat("Simulated: N =", N, "T =", T, "AR =", A_name, "Sigma =", Sig, "\n")
      }
    }
  }
}

## ----------------------------
# Transform -------------------

vals <- get_all_numeric_vals(data_list.VAR.scaled.long.04)

# sanity checks
length(vals) # 4800000 - seems legit

lower <- quantile(vals, 0.01, na.rm = TRUE)  # 1st percentile
upper <- quantile(vals, 0.99, na.rm = TRUE)  # 99th percentile


## ----------------------------
# Scale data globally ---------

data_list.VAR.scaled.long.04 <- add_likert_next_to_raw_quantile(
  data_list.VAR.scaled.long.04,
  lower = lower,
  upper = upper,
  new_min = 1,
  new_max = 100
)


saveRDS(data_list.VAR.scaled.long.04,
     file = "04/data/data_list-VAR-scaled-long-04.rds")

# data_list.VAR.scaled.long.04 <- readRDS("04/data/data_list-VAR-scaled-long-04.rds")

## ------------------------
# BISTABLE MODEL ----------

data_list.bistable.long <- readRDS("04/data/data_list_bistable_long.rds")


## ------------------------
# APPLY EDM TESTS ON LONG -
## ------------------------

# APPLY EDM COMPLEXITY TESTS TO VAR MODELS

VAR_EDM_results.1_100.long.04 <- list()

for (N in N_list.long) {
  N_name <- paste0("N", N)
  VAR_EDM_results.1_100.long.04[[N_name]] <- list()
  
  for (T in T_list.long) {
    T_name <- paste0("T", T)
    VAR_EDM_results.1_100.long.04[[N_name]][[T_name]] <- list()
    
    for (A in names(A_list.long)) {
      A_name <- paste0("A", A)
      VAR_EDM_results.1_100.long.04[[N_name]][[T_name]][[A_name]] <- list()
      
      for (S in Sigmas.long) {
        sigma_name <- paste0("Sigma", S)
        
        cat("Running EDM (parallel) for:",
            N_name, T_name, A_name, sigma_name, "\n")
        
        data_all <- data_list.VAR.scaled.long.04[[N_name]][[T_name]][[A_name]][[sigma_name]][["likert_1_100"]]
        
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
        
        VAR_EDM_results.1_100.long.04[[N_name]][[T_name]][[A_name]][[sigma_name]] <- person_results
      }  
    }
  }
}

saveRDS(VAR_EDM_results.1_100.long.04,
        file = "04/results/VAR_EDM-long-list.rds")

out <- list()

for (N_name in names(VAR_EDM_results.1_100.long.04)) {
  N_list <- VAR_EDM_results.1_100.long.04[[N_name]]
  
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
            mean  = mean(Value, na.rm = TRUE),
            sd    = sd(Value, na.rm = TRUE),
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
VAR_EDM.final_df <- VAR_EDM.final_df[, c("N", "T", "A", "Sigma", "Metric", "mean", "sd")]

# VAR_EDM.final_df$A[VAR_EDM.final_df$A == "AA1"] <- "A0.2"
# VAR_EDM.final_df$A[VAR_EDM.final_df$A == "AA2"] <- "A0.4"
# VAR_EDM.final_df$A[VAR_EDM.final_df$A == "AA3"] <- "A0.6"
# VAR_EDM.final_df$A[VAR_EDM.final_df$A == "AA4"] <- "A0.8"

saveRDS(VAR_EDM.final_df,
        file = "04/results/VAR_EDM-long-df.rds")


# APPLY EDM TESTS TO BISTABLE MODEL

BISTABLE_EDM_results.1_100.long <- list()

for (N_ in N_list.long) {
  N_name <- paste0("N", N_)
  BISTABLE_EDM_results.1_100.long[[N_name]] <- list()
  
  for (T in T_list.long) {
    T_name <- paste0("T", T)
    BISTABLE_EDM_results.1_100.long[[N_name]][[T_name]] <- list()
    
    data_all <- data_list.bistable.long[[N_name]][[T_name]][["likert_1_100"]]
    
    cols <- colnames(data_all[[1]])
    
    ## -----------------------------
    ## PARALLEL over persons
    ## -----------------------------
    
    person_results <- future_map(
      data_all,
      function(ts_data) {
        EDM_tests(ts_data,
                  EDM.include = c("E_opt", "SMap"))
      },
      .options = furrr_options(seed = TRUE)
    )
    
    BISTABLE_EDM_results.1_100.long[[N_name]][[T_name]] <- person_results
    
    cat("DONE:", N_name, T_name, A_name, sigma_name, "\n")
    
  }
}


out_bistable <- list()

for (N_name in names(BISTABLE_EDM_results.1_100.long)) {
  N_list <- BISTABLE_EDM_results.1_100.long[[N_name]]
  
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
        mean = mean(Value, na.rm = TRUE),
        sd = sd(Value, na.rm = TRUE),
        .groups = "drop"
      )
    
    summary_df$N <- N_name
    summary_df$T <- T_name
    
    out_bistable[[length(out_bistable) + 1]] <- summary_df
  }
}

BISTABLE_EDM.final_df <- bind_rows(out_bistable)
BISTABLE_EDM.final_df <- BISTABLE_EDM.final_df[, c("N", "T", "Metric", "mean", "sd")]


saveRDS(BISTABLE_EDM.final_df,
        file = "04/results/BISTABLE_EDM-long-df.rds")


## ---------------------------
# LORENZ SYSTEM LONG ---------

library(deSolve)

# from https://www.sixhat.net/lorenz-attractor-in-r.html
parameters <- c(s = 10, r = 28, b = 8/3)
state <- c(X = 0, Y = 1, Z = 1)

Lorenz <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX <- s * (Y - X)
    dY <- X * (r - Z) - Y
    dZ <- X * Y - b * Z
    list(c(dX, dY, dZ))
  })
}

# Simulate Lorenz data for different T
data_list.Lorenz.long <- list()

for (i in 1:length(T_list.long)) {
  T <- T_list.long[[i]]
  T_name <- paste0("T", T)
  
  times <- seq(0, 50, by = 50 / T)[1:T]
  out <- ode(y = state, times = times, func = Lorenz, parms = parameters)
  
  out.df <- as.data.frame(out)
  
  out.df[, 2:4] <- lapply(out.df[, 2:4], 
                          min_max_to_scale, 
                          new_min = 1,
                          new_max = 100)
  
  data_list.Lorenz.long[[T_name]] <- out.df[,2:4]
}

EDM_results.lorenz.long <- lapply(data_list.Lorenz.long, EDM_tests)

# Aggregate Lorenz into df
# Initialize empty data.frame
Lorenz.X.final_df.long <- data.frame(
  N = numeric(),
  T = numeric(),
  Metric = character(),
  mean = numeric(),
  sd = numeric(),
  Sigma = numeric(),
  A = numeric(),
  Source = character(),
  stringsAsFactors = FALSE
)

for (T_name in names(EDM_results.lorenz.long)) {
  T_num <- as.numeric(sub("T", "", T_name))
  print(T_num)
  
  df <- EDM_results.lorenz.long[[T_name]][, "X", drop = FALSE]
  
  temp <- data.frame(
    Metric = rownames(df),
    N = "N100",
    T = T_name,
    mean = df,
    sd = 0,
    Sigma = "LORENZ",
    A = "LORENZ",
    Source = "LORENZ",
    stringsAsFactors = FALSE
  )
  
  Lorenz.X.final_df.long <- rbind(Lorenz.X.final_df.long, temp)
}

colnames(Lorenz.X.final_df.long) <- c("Metric", "N", "T", "mean","sd","Sigma","A", "Source")
rownames(Lorenz.X.final_df.long) <- c()

saveRDS(Lorenz.X.final_df.long,
        file = "04/results/LORENZ_EDM-long-df.rds")

