############################################################
## VAR DATA SIMULATION UTILITIES
## Author: ---
## Purpose: Simulate multivariate VAR(1) data and
##          discretize to Likert-type scales
############################################################

## =========================================================
## 1. Min–max scaling to Likert scale
## =========================================================

min_max_to_scale <- function(x, new_min = 1, new_max = 7) {
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  
  # Edge case: constant time series
  if (x_max == x_min) {
    mid <- round((new_min + new_max) / 2)
    return(as.integer(rep(mid, length(x))))
  }
  
  # Linear rescaling
  x_scaled <- new_min + (x - x_min) *
    ((new_max - new_min) / (x_max - x_min))
  
  # Round and clip
  x_round <- pmin(new_max, pmax(new_min, round(x_scaled)))
  
  as.integer(x_round)
}


## =========================================================
## 2. Simulate a single time series
## =========================================================

simulate_ar1 <- function(
    T,
    phi,        # AR(1) coefficient (scalar)
    sigma = 1,  # innovation SD
    burnin = 0,
    init = NULL
) {
  stopifnot(
    is.numeric(phi), length(phi) == 1,
    abs(phi) < 1,   # stationarity condition
    sigma > 0
  )
  
  TT <- T + burnin
  Y  <- numeric(TT)
  
  Y[1] <- if (is.null(init)) rnorm(1) else init
  
  for (t in 2:TT) {
    Y[t] <- phi * Y[t - 1] + rnorm(1, mean = 0, sd = sigma)
  }
  
  if (burnin > 0) Y <- Y[(burnin + 1):TT]
  
  data.frame(V1 = Y)
}


simulate_var1 <- function(
    T,
    A,
    Sigma,
    burnin = 0,
    init = NULL
) {
  stopifnot(
    is.matrix(A),
    nrow(A) == ncol(A),
    is.matrix(Sigma),
    nrow(Sigma) == ncol(Sigma),
    nrow(A) == nrow(Sigma)
  )
  
  K <- nrow(A)
  TT <- T + burnin
  
  Y <- matrix(0, TT, K)
  
  if (is.null(init)) {
    Y[1, ] <- rnorm(K)
  } else {
    stopifnot(length(init) == K)
    Y[1, ] <- init
  }
  
  for (t in 2:TT) {
    eps <- MASS::mvrnorm(1, rep(0, K), Sigma)
    Y[t, ] <- A %*% Y[t - 1, ] + eps
  }
  
  if (burnin > 0) {
    Y <- Y[(burnin + 1):TT, , drop = FALSE]
  }
  
  as.data.frame(Y)
}


## =========================================================
## 3. Simulate data for N subjects
## =========================================================

simulate_ar1_subjects <- function(
    N,
    T,
    phi,
    sigma = 1,
    burnin = 0
) {
  lapply(seq_len(N), function(i) {
    simulate_ar1(
      T      = T,
      phi    = phi,
      sigma  = sigma,
      burnin = burnin
    )
  })
}


simulate_var_subjects <- function(
    N,
    T,
    A,
    Sigma,
    burnin = 0
) {
  K <- nrow(A)
  
  raw_data <- vector("list", N)
  
  for (i in seq_len(N)) {
    df <- simulate_var1(
      T = T,
      A = A,
      Sigma = Sigma,
      burnin = burnin
    )
    colnames(df) <- paste0("V", seq_len(K))
    raw_data[[i]] <- df
  }
  
  raw_data
}


## =========================================================
## 4. Discretize subject data to Likert scales
## =========================================================

discretize_likert <- function(
    raw_data,
    scale_min = 1,
    scale_max = 7,
    suffix = NULL
) {
  likert_data <- lapply(raw_data, function(df) {
    out <- as.data.frame(
      lapply(df, min_max_to_scale,
             new_min = scale_min,
             new_max = scale_max)
    )
    if (!is.null(suffix)) {
      colnames(out) <- paste0(colnames(out), suffix)
    }
    out
  })
  
  likert_data
}


## =========================================================
## 5. Full simulation grid: N × T
## =========================================================

simulate_ar1_grid <- function(
    N_list,
    T_list,
    phi,
    sigma = 1,
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
    data_list[[paste0("N", N)]] <- list()
    
    for (T in T_list) {
      
      raw_data <- simulate_ar1_subjects(
        N      = N,
        T      = T,
        phi    = phi,
        sigma  = sigma,
        burnin = burnin
      )
      
      # Discretise to Likert scales (reuses your existing function)
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
      
      data_list[[paste0("N", N)]][[paste0("T", T)]] <- c(
        list(raw = raw_data),
        likert_data
      )
      
      if (verbose) cat("Simulated: N =", N, "| T =", T, "\n")
    }
  }
  
  data_list
}


simulate_var_grid <- function(
    N_list,
    T_list,
    A,
    Sigma,
    burnin = 0,
    likert_scales = list(
      "likert_1_7"   = c(1, 7),
      "likert_1_100" = c(1, 100)
    ),
    seed = NULL,
    verbose = TRUE
) {
  if (!is.null(seed)) set.seed(seed)
  
  data_list <- list()
  
  for (N in N_list) {
    data_list[[paste0("N", N)]] <- list()
    
    for (T in T_list) {
      
      raw_data <- simulate_var_subjects(
        N = N,
        T = T,
        A = A,
        Sigma = Sigma,
        burnin = burnin
      )
      
      likert_data <- list()
      
      for (nm in names(likert_scales)) {
        sc <- likert_scales[[nm]]
        likert_data[[nm]] <- discretize_likert(
          raw_data,
          scale_min = sc[1],
          scale_max = sc[2],
          suffix = paste0("_", nm)
        )
      }
      
      data_list[[paste0("N", N)]][[paste0("T", T)]] <- c(
        list(raw = raw_data),
        likert_data
      )
      
      if (verbose) {
        cat("Simulated: N =", N, "T =", T, "\n")
      }
    }
  }
  
  data_list
}


## =========================================================
## 6. Differential Base Equation
## =========================================================

# from Haslbeck & Ryan (2022) - found here:
# https://github.com/jmbh/RecoveringWithinPersonDynamics/tree/master

dXdt <- function(i,      # variable index
                 X,      # current state vector
                 r,      # growth rates vector
                 C,      # coupling matrix
                 mu,     # additive constant
                 epsilon, # noise term
                 timestep = 1) {
  
  # Compute derivative for variable i
  dX <- r[i] * X[i] + X[i] * sum(C[i, ] * X) + mu + epsilon
  dX * timestep
}



## =========================================================
## 7. Simulate a single bistable trajectory
## =========================================================

# from Haslbeck & Ryan (2022) - found here:
# https://github.com/jmbh/RecoveringWithinPersonDynamics/tree/master

simulate_bistable_trajectory <- function(time,        # total simulation time
                                         timestep,    # integration timestep
                                         p,           # number of variables
                                         C,           # coupling matrix
                                         mu,          # additive constant
                                         r = rep(1, p), # growth rates
                                         init = NULL, # initial state
                                         noiseSD = 1,
                                         noise = TRUE,
                                         pbar = TRUE) {
  
  nIter <- ceiling(time / timestep)  # number of integration steps
  
  # Initialize trajectory matrix
  X <- matrix(NA, nIter, p)
  
  # Initial state
  if (is.null(init)) {
    # use same starting variance as original
    X[1, ] <- c(rnorm(1, 1.3, 1),
                rnorm(1, 1.3, 1),
                rnorm(1, 4.8, 1),
                rnorm(1, 4.8, 1))
  } else {
    stopifnot(length(init) == p)
    X[1, ] <- init
  }
  
  # Progress bar
  if (pbar) pb <- txtProgressBar(min = 2, max = nIter, initial = 0, style = 3)
  
  # Euler integration
  for (j in 2:nIter) {
    for (i in 1:p) {
      epsilon <- if (noise) rnorm(1, 0, noiseSD) else 0
      X[j, i] <- X[j-1, i] + dXdt(i = i, X = X[j-1, ], r = r, C = C, mu = mu, epsilon = epsilon, timestep = timestep)
    }
    if (pbar) setTxtProgressBar(pb, j)
  }
  
  if (pbar) close(pb)
  
  # Return time + variables as data.frame
  v_time <- seq(0, time, length.out = nIter)
  out <- cbind(time = v_time, X)
  colnames(out) <- c("time", paste0("x", 1:p))
  as.data.frame(out)
}

