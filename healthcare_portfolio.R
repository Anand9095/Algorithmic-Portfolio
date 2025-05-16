
# Purpose: Institutional Healthcare Long-Only Portfolio (Min Vol + Max Sharpe)

# ================================
# 1. Load Libraries and Define Setup
# ================================
packages <- c("tidyquant", "openxlsx", "readxl", "xts", "dplyr", "PerformanceAnalytics",
              "quadprog", "ggplot2", "tidyr", "gganimate", "av", "zoo")
installed <- rownames(installed.packages())
for (pkg in packages) if (!(pkg %in% installed)) install.packages(pkg)
invisible(lapply(packages, library, character.only = TRUE))

# ================================
# 2. Download Healthcare Stock Prices from Yahoo
# ================================
healthcare_tickers <- c("LLY", "UNH", "JNJ", "ABBV", "ABT", "MRK", "TMO", "ISRG", "AMGN", "BSX")
start_date <- as.Date("2014-01-01")
end_date <- Sys.Date()

wb <- createWorkbook()
for (ticker in healthcare_tickers) {
  df <- tq_get(ticker, from = start_date, to = end_date)
  df <- df %>% select(date, open, high, low, close, volume, adjusted)
  addWorksheet(wb, ticker)
  writeData(wb, sheet = ticker, x = df)
}
saveWorkbook(wb, "C:/r/healthcare.xlsx", overwrite = TRUE)

# ================================
# 3. Read Excel Data and Build Return Series
# ================================
file_path <- "C:/r/healthcare.xlsx"
adjusted_prices_list <- list()
for (ticker in healthcare_tickers) {
  df <- read_excel(file_path, sheet = ticker)
  if ("adjusted" %in% names(df) && "date" %in% names(df)) {
    adjusted_prices_list[[ticker]] <- xts(df$adjusted, order.by = as.Date(df$date))
  }
}
prices_xts <- do.call(merge, adjusted_prices_list)
returns_xts <- na.omit(Return.calculate(prices_xts, method = "log"))
saveRDS(returns_xts, "data/healthcare_log_returns.rds")

# ================================
# 4. Rolling Optimization: Min Volatility and Max Sharpe
# ================================
window_size <- 120
rebalance_every <- 30
Rf_daily <- 0.045 / 252

weights_list <- list()
weights_sharpe <- list()
rebalance_points <- seq(window_size, nrow(returns_xts) - rebalance_every, by = rebalance_every)

for (i in rebalance_points) {
  window_returns <- returns_xts[(i - window_size + 1):i, ]
  mu_hat <- colMeans(window_returns)
  cov_matrix <- cov(window_returns)
  N <- ncol(window_returns)
  Amat <- cbind(rep(1, N), diag(N))
  bvec <- c(1, rep(0, N))
  meq <- 1
  
  # Min Vol
  minvol <- solve.QP(2 * cov_matrix, rep(0, N), Amat, bvec, meq)
  weights_list[[as.character(index(returns_xts)[i])]] <- minvol$solution
  
  # Max Sharpe
  excess_returns <- mu_hat - Rf_daily
  result <- tryCatch({
    solve.QP(2 * cov_matrix, dvec = -excess_returns, Amat, bvec, meq)
  }, error = function(e) NULL)
  if (!is.null(result)) weights_sharpe[[as.character(index(returns_xts)[i])]] <- result$solution
}

# ================================
# 5. Simulate Portfolio Performance (30-Day Holding)
# ================================
portfolio_ret_minvol <- xts(rep(NA, nrow(returns_xts)), order.by = index(returns_xts))
portfolio_ret_sharpe <- xts(rep(NA, nrow(returns_xts)), order.by = index(returns_xts))

for (dt_str in names(weights_list)) {
  dt <- as.Date(dt_str)
  i <- which.min(abs(as.numeric(index(returns_xts)) - as.numeric(dt)))
  if ((i + 30) > nrow(returns_xts)) next
  next_returns <- returns_xts[(i + 1):(i + 30), ]
  r_matrix <- as.matrix(next_returns)
  
  # Min Vol
  w_minvol <- weights_list[[dt_str]]
  if (length(w_minvol) == ncol(r_matrix))
    portfolio_ret_minvol[(i + 1):(i + 30)] <- as.numeric(r_matrix %*% w_minvol)
  
  # Max Sharpe
  if (!is.null(weights_sharpe[[dt_str]])) {
    w_sharpe <- weights_sharpe[[dt_str]]
    if (length(w_sharpe) == ncol(r_matrix))
      portfolio_ret_sharpe[(i + 1):(i + 30)] <- as.numeric(r_matrix %*% w_sharpe)
  }
}

# ================================
# 6. Performance Summary and Visualization
# ================================
combined_returns <- na.omit(merge(portfolio_ret_minvol, portfolio_ret_sharpe))
colnames(combined_returns) <- c("Min Volatility", "Max Sharpe")
charts.PerformanceSummary(combined_returns, main = "Strategy Performance Comparison")
print(table.AnnualizedReturns(combined_returns))

# ================================
# 7. Animate Cumulative Return and Weight Allocation
# ================================
cumulative_minvol <- na.omit(cumsum(portfolio_ret_minvol))
weights_df <- do.call(rbind, weights_list) %>% as.data.frame()
weights_df$date <- as.Date(names(weights_list))
weights_df <- weights_df %>% relocate(date)
colnames(weights_df)[-1] <- healthcare_tickers

weights_long <- weights_df %>%
  pivot_longer(-date, names_to = "Ticker", values_to = "Weight")

# Weights animation
p_weights <- ggplot(weights_long, aes(x = Ticker, y = Weight, fill = Ticker)) +
  geom_bar(stat = "identity") +
  labs(title = "Portfolio Weights on {frame_time}", x = "Stock", y = "Weight") +
  transition_time(date) +
  ease_aes('linear') +
  theme_minimal()
animate(p_weights, renderer = av_renderer("weights_simulation.mp4"), width = 800, height = 500, fps = 10, duration = 15)

# Cumulative return animation
cumulative_df <- data.frame(date = index(cumulative_minvol), CumulativeReturn = as.numeric(cumulative_minvol))
p_return <- ggplot(cumulative_df, aes(x = date, y = CumulativeReturn)) +
  geom_line(color = "steelblue", size = 1.2) +
  labs(title = "Cumulative Return Progression", x = "Date", y = "Cumulative Return") +
  transition_reveal(date) +
  theme_minimal()
animate(p_return, renderer = av_renderer("returns_simulation.mp4"), width = 800, height = 500, fps = 10, duration = 15)

# ================================
# 8. Rolling VaR & ES Analysis
# ================================
var_es_results <- data.frame()
for (dt_str in names(weights_list)) {
  dt <- as.Date(dt_str)
  i <- which.min(abs(as.numeric(index(returns_xts)) - as.numeric(dt)))
  if ((i + 30) > nrow(returns_xts)) next
  window_returns <- returns_xts[(i + 1):(i + 30), ]
  r_matrix <- as.matrix(window_returns)
  
  w_minvol <- weights_list[[dt_str]]
  if (!is.null(w_minvol) && length(w_minvol) == ncol(r_matrix)) {
    port_minvol <- as.numeric(r_matrix %*% w_minvol)
    var_min <- VaR(port_minvol, p = 0.99, method = "historical")
    es_min <- ES(port_minvol, p = 0.99, method = "historical")
    var_es_results <- rbind(var_es_results, data.frame(Date = dt, Strategy = "MinVol", VaR_1pct = var_min, ES_1pct = es_min))
  }
  
  w_sharpe <- weights_sharpe[[dt_str]]
  if (!is.null(w_sharpe) && length(w_sharpe) == ncol(r_matrix)) {
    port_sharpe <- as.numeric(r_matrix %*% w_sharpe)
    var_sharp <- VaR(port_sharpe, p = 0.99, method = "historical")
    es_sharp <- ES(port_sharpe, p = 0.99, method = "historical")
    var_es_results <- rbind(var_es_results, data.frame(Date = dt, Strategy = "MaxSharpe", VaR_1pct = var_sharp, ES_1pct = es_sharp))
  }
}

# Plot VaR & ES
ggplot(var_es_results, aes(x = Date, y = VaR_1pct, color = Strategy)) +
  geom_line() +
  labs(title = "Rolling 1% VaR", y = "VaR", x = "Date") +
  theme_minimal()

ggplot(var_es_results, aes(x = Date, y = ES_1pct, color = Strategy)) +
  geom_line() +
  labs(title = "Rolling 1% Expected Shortfall", y = "ES", x = "Date") +
  theme_minimal()


code portfolio , healthcare
