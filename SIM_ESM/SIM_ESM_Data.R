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

library(dplyr)
library(tidyr)

library(ggplot2)

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

# Parameter settings lead to model where approx. 22.5% of all variance is 
# explained (check by fitting a VAR)

T_list <- c(25, 50, 75, 100, 150, 200)
N_list <- c(100, 200)


### 2 --- MIN-MAX SCALING FUNCTION --- ###

min_max_to_scale <- function(x, new_min = 1, new_max = 7) {
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  
  # Sonderfall: konstante Zeitreihe
  if (x_max == x_min) {
    mid <- round((new_min + new_max) / 2)
    return(rep(mid, length(x)))
  }
  
  # Lineare Skalierung
  x_scaled <- new_min + (x - x_min) * ((new_max - new_min) / (x_max - x_min))
  
  # Rundung auf ganze Werte + Begrenzung
  x_round <- pmin(new_max, pmax(new_min, round(x_scaled)))
  
  return(as.integer(x_round))
}


### 3 --- SIMULATION LOOP: N × T --- ###

data_list <- list()

set.seed(2310)

for (N in N_list) {
  
  data_list[[paste0("N", N)]] <- list()
  
  for (T in T_list) {
    
    raw_data_list <- vector("list", N)
    
    for (i in 1:N) {
      Y <- matrix(0, T, 5)
      Y[1, ] <- rnorm(5)
      
      for (t in 2:T) {
        eps <- mvrnorm(1, rep(0,5), Sigma)
        Y[t, ] <- A %*% Y[t-1, ] + eps
      }
      
      raw_data_list[[i]] <- as.data.frame(Y)
      colnames(raw_data_list[[i]]) <- paste0("V", 1:5)
    }
    
    raw_data <- raw_data_list
    
    # Likert 1-7
    likert_1_7 <- lapply(raw_data, function(df) {
      df <- as.data.frame(lapply(df, min_max_to_scale, new_min = 1, new_max = 7))
      colnames(df) <- paste0("V", 1:5, "_lik7")  # set names for every person
      return(df)
    })
    
    # Likert 1-100
    likert_1_100 <- lapply(raw_data, function(df) {
      df <- as.data.frame(lapply(df, min_max_to_scale, new_min = 1, new_max = 100))
      colnames(df) <- paste0("V", 1:5, "_lik100")
      return(df)
    })
    
    # save to data_list
    data_list[[paste0("N", N)]][[paste0("T", T)]] <- list(
      raw = raw_data,
      likert_1_7 = likert_1_7,
      likert_1_100 = likert_1_100
    )
    
    cat("Simulated: N =", N, "T =", T, "\n")
  }
}

### 4 --- APPLY TESTS --- ###

# Progress bars use base R only
n_blocks <- length(data_list)
pb_outer <- txtProgressBar(min = 0, max = n_blocks, style = 3)

# Initialize two lists to store results for each scale
test_results_1_7 <- list()
test_results_1_100 <- list()

scales <- list(
  "1_7" = "likert_1_7",
  "1_100" = "likert_1_100"
)

# Loop over scales
for (scale_name in names(scales)) {
  
  test_results <- list()  # temporary storage for this scale
  
  # Loop over sample sizes N
  for (N_name in names(data_list)) {
    
    test_results[[N_name]] <- list()  # storage per N
    
    # Loop over time series lengths T
    for (T_name in names(data_list[[N_name]])) {
      
      lik_data <- data_list[[N_name]][[T_name]][[scales[[scale_name]]]]  # List of persons
      N_persons <- length(lik_data)
      
      # Skip if no data
      if (N_persons == 0) next
      
      cols <- colnames(lik_data[[1]])
      
      # Robust initialization of array: Tests × Variables × Persons
      res_array <- array(NA, dim = c(length(tests), length(cols), N_persons),
                         dimnames = list(tests, cols, paste0("P", 1:N_persons)))
      
      pb_outer <- txtProgressBar(min = 0, max = N_persons, style = 3)
      
      # Loop over persons
      for (p in 1:N_persons) {
        ts_data <- lik_data[[p]]  # Dataframe of variables for this person
        
        # Loop over variables
        for (j in seq_along(cols)) {
          col <- cols[j]
          ts <- ts_data[[col]]
          N <- length(ts)
          
          # --- Bartels rank test ---
          bartels_p <- tryCatch(bartels.rank.test(ts, alternative = "two.sided")$p.value, error = function(e) NA)
          res_array["Bartels", col, p] <- if (is.null(bartels_p) || length(bartels_p) == 0) NA else bartels_p
          
          # --- KPSS test level stationarity ---
          kpss_level_p <- tryCatch(kpss.test(ts, null = "Level", lshort = TRUE)$p.value, error = function(e) NA)
          res_array["KPSS.level", col, p] <- if (is.null(kpss_level_p) || length(kpss_level_p) == 0) NA else kpss_level_p

          # --- KPSS test trend ---
          kpss_trend_p <- tryCatch(kpss.test(ts, null = "Trend", lshort = TRUE)$p.value, error = function(e) NA)
          res_array["KPSS.trend", col, p] <- if (is.null(kpss_trend_p) || length(kpss_trend_p) == 0) NA else kpss_trend_p

          # --- Keenan's nonlinearity test ---
          keenan_p <- tryCatch(keenanTest(ts)$p.value, error = function(e) NA)
          res_array["Keenan", col, p] <- if (is.null(keenan_p) || length(keenan_p) == 0) NA else keenan_p

          # --- Tsay's nonlinearity test ---
          tsay_p <- tryCatch(tsayTest(ts)$p.value, error = function(e) NA)
          res_array["Tsay", col, p] <- if (is.null(tsay_p) || length(tsay_p) == 0) NA else tsay_p

          # --- PACF analysis ---
          # pacf_vals <- tryCatch(pacf(ts, plot=FALSE)$acf, error=function(e) NA)
          # if (is.numeric(pacf_vals) && length(pacf_vals) > 1 && !all(is.na(pacf_vals))) {
          #   thresh <- 2 / sqrt(N)
          #   sig_lags <- which(abs(pacf_vals) > thresh)
          #   res_array["PAFC.n", col, p] <- length(sig_lags)
          #   res_array["PAFC.max_lag", col, p] <- if(length(sig_lags) > 0) max(sig_lags) else NA
          # }
          
          # --- TV-AR GAM ---
          tt <- 1:(N-1)
          tv <- tryCatch(
            gam(ts[2:N] ~ s(tt, by=ts[1:(N-1)], k=min(10,N-2), bs="tp")),
            error=function(e) NULL
          )
          if (!is.null(tv)) {
            stvar <- summary(tv)
            if (!is.null(stvar$s.table) && nrow(stvar$s.table) >= 1) {
              # res_array["TV-AR.EDF", col, p] <- round(stvar$s.table[1,"edf"],3)
              res_array["TV-AR.p", col, p] <- round(stvar$s.table[1,"p-value"],4)
            }
          }
          
          # --- Change Point Analysis ---
          # R <- ceiling(20 / ALPHA_LEVEL)
          # cp.out <- tryCatch(e.divisive(matrix(ts), R=R, sig.lvl=ALPHA_LEVEL), error=function(e) NULL)
          # if (!is.null(cp.out)) {
          #   cp.n <- length(which(cp.out$p.values < ALPHA_LEVEL))
          #   res_array["CP.n", col, p] <- cp.n
          # }
          
          
        } # end variable loop
        
        setTxtProgressBar(pb_outer, p)
      } # end person loop
      
      close(pb_outer)
      
      # --- Aggregate across persons ---
      sig_rates <- apply(res_array, c(1,2), function(x) mean(x < ALPHA_LEVEL, na.rm=TRUE))
      
      test_results[[N_name]][[T_name]] <- list(
        raw_array = res_array,
        sig_rate = sig_rates
      )
      
      cat("Finished:", scale_name, N_name, T_name, "\n")
      
    } # end T loop
  } # end N loop
  
  # Save results in the proper scale list
  if (scale_name == "1_7") test_results_1_7 <- test_results
  if (scale_name == "1_100") test_results_1_100 <- test_results
}



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
