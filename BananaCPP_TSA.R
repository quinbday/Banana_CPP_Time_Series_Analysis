library(quantmod)
library(ggplot2)
library(forecast)
library(zoo)
library(gridExtra)
library(tseries)   
library(lmtest)    

getSymbols("APU0000711211", src = "FRED")

bananas <- data.frame(
  Date  = zoo::index(APU0000711211),
  Price = as.numeric(APU0000711211)
)

bananas <- subset(bananas, Date >= "2010-01-01" & Date <= "2026-01-01")
bananas$Price[is.na(bananas$Price)] <- 0.667

start_year  <- as.integer(format(min(bananas$Date), "%Y"))
start_month <- as.integer(format(min(bananas$Date), "%m"))
bananas_ts  <- ts(bananas$Price, start = c(start_year, start_month), frequency = 12)

banana_df <- data.frame(
  Date  = as.numeric(time(bananas_ts)),
  Price = as.numeric(bananas_ts)
)

bimonthly_breaks <- seq(
  from = min(banana_df$Date), 
  to   = max(banana_df$Date), 
  by   = 2/12
)

ggplot(banana_df, aes(x = Date, y = Price)) +
  geom_vline(xintercept = bimonthly_breaks, color = "red", linetype = "dashed", alpha = 0.3) +
  geom_line(color = "orange", linewidth = 0.9) +
  labs(
    title = "Historical Banana Prices: Bimonthly Grid",
    x = "Year",
    y = "Price per Pound (USD)"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  )


profile_distribution <- function(x, label) {
  x   <- as.numeric(x)
  n   <- length(x)
  sw  <- shapiro.test(x)
  cat(sprintf("  %-40s  n=%d\n", label, n))
  cat(sprintf("    Mean=%.4f  SD=%.4f\n",
              mean(x), sd(x)))
  cat(sprintf("    Shapiro-Wilk  W=%.5f  p=%.6f  -> %s\n",
              sw$statistic, sw$p.value,
              ifelse(sw$p.value > 0.05, "NORMAL",
                     ifelse(sw$statistic > 0.97, "minor tails (W>0.97, practically normal)",
                            ifelse(sw$statistic > 0.95, "mild tails (W>0.95, acceptable)",
                                   "non-normal (W<0.95, investigate)")))))
  
}

# Raw prices
prof_raw   <- profile_distribution(bananas_ts, "Raw prices (bananas_ts)")
# Log prices
prof_log   <- profile_distribution(log(bananas_ts), "Log prices  log(bananas_ts)")
# First differences of log (log returns)
prof_ret   <- profile_distribution(diff(log(bananas_ts)), "Log returns  diff(log(bananas_ts))")
# STL remainder of raw prices
stl_raw    <- stl(bananas_ts, s.window = "periodic")
prof_rem   <- profile_distribution(stl_raw$time.series[,"remainder"],
                                   "STL remainder of raw prices")
# STL remainder of log prices
stl_log    <- stl(log(bananas_ts), s.window = "periodic")
prof_remlog <- profile_distribution(stl_log$time.series[,"remainder"],
                                    "STL remainder of log prices")

lambda <- BoxCox.lambda(bananas_ts)
print(lambda)
lambda     <- 1
bananas_bc <- BoxCox(bananas_ts, lambda)

p_full <- ggplot(bananas, aes(x=Date, y=Price)) +
  geom_line(color="goldenrod", linewidth=0.8) +
  labs(title    = "STEP 4: Historical Banana Prices",
       subtitle = "Source: FRED (APU0000711211)",
       x="Year", y="Price per Pound (USD)") +
  theme_minimal()
print(p_full)

bananas_stl <- stl(bananas_bc, s.window = "periodic")

# Extract the seasonal component as a time series
seasonal_ts <- bananas_stl$time.series[, "seasonal"]
# For a monthly series (frequency = 12)
seasonal_cycle <- ts(seasonal_ts, frequency = frequency(seasonal_ts))
boxplot(seasonal_cycle ~ cycle(seasonal_cycle),
        main = "Monthly Seasonal Effects",
        xlab = "Month", ylab = "Seasonal Effect")
plot(seasonal_ts, ylab = "Seasonal Effect", main = "STL Seasonal Breakdown")

p_stl <- autoplot(bananas_stl) +
  labs(title = "STEP 5: STL Decomposition")
print(p_stl)

remainder_stl <- bananas_stl$time.series[, "remainder"]

print(ggAcf(remainder_stl, lag.max=48) +
        labs(title    = "ACF of STL Remainder from 2010-2026"))

print(ggPacf(remainder_stl, lag.max=48) +
        labs(title    = "PACF of STL Remainder from 1990-2026"))

run_stationarity <- function(series, label) {
  adf_r  <- adf.test(series)
  kpss_r <- kpss.test(series, null="Level")
  
  adf_p_exact  <- as.numeric(adf_r$p.value)
  kpss_p_exact <- as.numeric(kpss_r$p.value)
  
  cat(sprintf("ADF  p-value = %.6f  ->  %s\n", adf_p_exact,
              ifelse(adf_p_exact > 0.05, "NON-STATIONARY (difference needed)", "STATIONARY")))
  cat(sprintf("KPSS p-value = %.6f  ->  %s\n\n", kpss_p_exact,
              ifelse(kpss_p_exact < 0.05, "NON-STATIONARY (difference needed)", "STATIONARY")))
  
  invisible(list(adf=adf_r, kpss=kpss_r))
}

bananas_diff1       <- diff(bananas_bc, differences=1)
bananas_diff_annual <- diff(bananas_bc, lag=12, differences=1)
bananas_diff12      <- diff(bananas_diff1, lag=12, differences=1)

run_stationarity(bananas_bc,         "Raw transformed series")
run_stationarity(bananas_diff1,       "After d=1 only")
run_stationarity(bananas_diff_annual, "After D=1 only (lag=12)")
run_stationarity(bananas_diff12,      "After d=1 AND D=1  <- target")

print(ggAcf(bananas_diff_annual, lag.max=48) +
        labs(title="ACF after Annual Diff Only (D=1)"))

print(ggAcf(bananas_diff1, lag.max=48) +
        labs(title="ACF after d=1 Only"))

print(ggPacf(bananas_diff1, lag.max=48) +
        labs(title="PACF after d=1 Only"))

print(ggAcf(bananas_diff12, lag.max=48) +
        labs(title="ACF after d=1, D=1"))

print(ggPacf(bananas_diff12, lag.max=48) +
        labs(title="PACF after d=1, D=1"))

fit_mean   <- meanf(bananas_ts,  h=24)
fit_naive  <- naive(bananas_ts,  h=24)
fit_drift  <- rwf(bananas_ts,    h=24, drift=TRUE)
fit_snaive <- snaive(bananas_ts, h=24)
fit_ar1    <- Arima(bananas_ts,  order=c(1,0,0), lambda=lambda)

if (TRUE) {
  
  results_list <- list()
  counter      <- 1
  
  p_range <- 0:1
  d_range <- 1
  q_range <- 0:2
  P_range <- 0:1
  D_range <- 0:1
  Q_range <- 0:1
  
  
  total_iters <- length(p_range) * length(d_range) * length(q_range) * 
    length(P_range) * length(D_range) * length(Q_range)
  
  
  cat(sprintf("%d combinations", total_iters))
  for (p in p_range) {
    for (d in d_range) {
      for (q in q_range) {
        for (P in P_range) {
          for (D in D_range) {
            for (Q in Q_range) {
              
              model_label <- sprintf("SARIMA(%d,%d,%d)(%d,%d,%d)", p, d, q, P, D, Q)
              
              fit <- tryCatch({
                Arima(bananas_ts, 
                      order    = c(p, d, q), 
                      seasonal = list(order = c(P, D, Q), period = 12), 
                      lambda   = lambda)
              }, error = function(e) { return(NULL) })
              
              if (!is.null(fit)) {
                acc_metrics <- accuracy(fit)
                
                resids <- as.numeric(residuals(fit))
                resids <- resids[!is.na(resids)]
                sw_metric <- tryCatch({ shapiro.test(resids) }, error = function(e) { list(statistic=0, p.value=0) })
                
                results_list[[counter]] <- data.frame(
                  Model        = model_label,
                  RMSE         = acc_metrics[1, "RMSE"],
                  MASE         = acc_metrics[1, "MASE"],
                  AICc         = fit$aicc,
                  Shapiro_p    = sw_metric$p.value,
                  stringsAsFactors = FALSE
                )
                counter <- counter + 1
              }
            }
          }
        }
      }
    }
  }
  
  grid_results <- do.call(rbind, results_list)
  grid_results <- grid_results[order(grid_results$MASE), ]
  
  rownames(grid_results) <- NULL
  print(head(grid_results, 15))
}

fit_temp <- auto.arima(bananas_ts, max.p = 1, d=1, max.D=1, max.q=2, max.P=1,max.Q=1)
fit_temp
autoplot(fit_temp)
fit_target <- Arima(bananas_ts, order=c(0,1,2), seasonal = list(order = c(0,0,0), period=12), lambda = lambda)
fit_alt_Q <- Arima(bananas_ts, order=c(1,1,2), seasonal = list(order = c(0,1,1), period=12), lambda = lambda)


extract_diagnostics <- function(model_obj, label, fit_df) {
  if (is.null(model_obj)) return(NULL)
  
  acc  <- accuracy(model_obj)[1, ]
  lb   <- Box.test(residuals(model_obj), lag = 36, type = "Ljung-Box", fitdf = fit_df)
  
  data.frame(
    Model         = label,
    LjungBox_Stat = lb$statistic,
    LjungBox_p    = lb$p.value,
    RMSE          = acc["RMSE"],
    MASE          = acc["MASE"],
    stringsAsFactors = FALSE
  )
}

diagnostic_grid <- do.call(rbind, list(
  extract_diagnostics(fit_target, "target", fit_df = 2),
  extract_diagnostics(fit_alt_Q, "alt_Q", fit_df = 2)
))

rownames(diagnostic_grid) <- NULL
print(diagnostic_grid)

fit_norm_winner <- fit_target
print(fit_auto)
print(fit_target)

get_in_sample_metrics <- function(model_obj, label) {
  if (is.null(model_obj)) return(NULL)
  
  acc <- accuracy(model_obj)
  if (nrow(acc) > 1) {
    metrics <- acc["Training set", ]
  } else {
    metrics <- acc[1, ]
  }
  
  data.frame(
    Model = label,
    RMSE  = metrics["RMSE"],
    MASE  = metrics["MASE"],
    stringsAsFactors = FALSE
  )
}

performance_matrix <- do.call(rbind, list(
  get_in_sample_metrics(fit_ar1,    "AR(1)"),
  get_in_sample_metrics(fit_mean,    "Mean Model"),
  get_in_sample_metrics(fit_naive,   "Naive Model"),
  get_in_sample_metrics(fit_drift,   "Random Walk w/ Drift"),
  get_in_sample_metrics(fit_snaive,  "Seasonal Naive Model"),
  get_in_sample_metrics(fit_target,  "0 1 2"),
  get_in_sample_metrics(fit_alt_Q,   "0 1 2 0 1 1")
))

rownames(performance_matrix) <- NULL
print(performance_matrix)

best_model_row <- performance_matrix[order(performance_matrix$RMSE), ][1, ]
best_name      <- best_model_row$Model

best_fit <- switch(best_name,
                   "Mean Model"                  = fit_mean,
                   "Naive Model"                 = fit_naive,
                   "Random Walk w/ Drift"        = fit_drift,
                   "Seasonal Naive Model"        = fit_snaive,
                   "AR(1) Model Structural"      = fit_target,
                   "SARIMA(1,1,1)(0,0,1)[12]"    = fit_alt_Q,
                   "Auto-ARIMA Optimized Matrix" = fit_auto)







fit_target
fit_alt_Q

fc_mean  <- forecast(fit_mean, h=12)
fc_naive <- forecast(fit_naive, h=12)
fc_drift <- forecast(fit_drift, h=12)
fc_snaive <- forecast(fit_snaive, h=12)
fc_target <- forecast(fit_target, h=12)
fc_alt_Q <- forecast(fit_alt_Q, h=12)
fc_ar1 <- forecast(fit_ar1, h=12)


forecast_list <- list(
  Target = fc_target, 
  AltQ = fc_alt_Q,
  Mean = fc_mean,
  Naive = fc_naive,
  AR1 = fc_ar1,
  SNaive = fc_snaive,
  Drift = fc_drift 
  
)

all_forecasts <- do.call(rbind, mapply(function(fc_obj, name) {
  data.frame(
    Month = 1:12,
    Model = name,
    Mean = round(as.numeric(fc_obj$mean), 3),
    Low = round(as.numeric(fc_obj$lower[, 2]), 3),
    High = round(as.numeric(fc_obj$upper[, 2]), 3)
  )
}, forecast_list, names(forecast_list), SIMPLIFY = FALSE))

print(all_forecasts)
