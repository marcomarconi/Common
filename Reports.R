
strategy_performance <- function(returns, dates = NULL, period = 252) {
  # if(any(is.na(returns)))
  #   stop("strategy_performance: NAs are not allowed, replace them with zeros")
  if (!is.null(dates) & class(dates) != "Date") {
    stop(paste("strategy_performance: dates should be of class Date"))
  }
  if (class(returns) != "numeric") {
    stop(paste("strategy_performance: returns should be of class numeric"))
  }
  if (is.null(dates)) {
    dates <- seq(as.Date("1970/01/01"), by = "day", length.out = length(returns))
  }
  df <- data.frame(Dates = dates, Returns = returns) %>% arrange(Dates)
  annual_returns <- group_by(df, year(Dates)) %>%
    summarise(Dates = first(year(Dates)), Returns = sum(Returns, na.rm = TRUE)) %>%
    pull(Returns)
  monthly_returns <- group_by(df, yearmonth(Dates)) %>%
    summarise(Dates = first(yearmonth(Dates)), Returns = sum(Returns, na.rm = TRUE)) %>%
    pull(Returns)
  weekly_returns <- group_by(df, yearweek(Dates)) %>%
    summarise(Dates = first(yearweek(Dates)), Returns = sum(Returns, na.rm = TRUE)) %>%
    pull(Returns)
  mean_ann_ret <- mean(annual_returns, na.rm = TRUE) * 100
  ann_sd <- sd(df$Returns, na.rm = TRUE) * sqrt(period) * 100
  sr <- mean(df$Returns, na.rm = TRUE) / sd(df$Returns, na.rm = TRUE) * sqrt(period)
  skew_ <- skew(weekly_returns)
  kurtosis_ <- kurtosis(weekly_returns)
  q <- quantile(df$Returns[df$Returns != 0], probs = c(0.01, 0.3, 0.7, 0.99), na.rm = TRUE)
  lower_tail <- as.numeric(q[1] / q[2] / 4.43)
  upper_tail <- as.numeric(q[4] / q[3] / 4.43)
  cum_returns <- cumsum(replace(df$Returns, is.na(df$Returns), 0))
  peak <- cummax(cum_returns)
  drawdown <- peak - cum_returns
  max_drawdown <- -(exp(drawdown[which.max(drawdown)]) - 1) * 100
  avg_drawdown <- -(exp(mean(drawdown)) - 1) * 100
  gsr <- sr * (1 + skew_/6*sr - (kurtosis_ - 3)/24*sr^2 )
  cum_annual_returns <- cumsum(replace(annual_returns, is.na(annual_returns), 0))
  r2 <- summary(lm(1:length(cum_annual_returns) ~ 0 + cum_annual_returns))$adj.r.squared
  # turnover <- round(length(rle(as.vector(na.omit(sign(returns))))$length) / (length(returns) / period), 1)
  results <- list(
    "Mean annual return" = mean_ann_ret,
    "Annualized SD" = ann_sd,
    "Sharpe ratio" = sr,
    "Skew" = skew_,
    "Lower tail" = lower_tail,
    "Upper tail" = upper_tail,
    "Max DD" = max_drawdown,
    "Avg DD" = avg_drawdown,
    "Adj Avg DD" = avg_drawdown / ann_sd,
    "GSR" = gsr,
    "R2" = r2
  )
  return(lapply(results, round, 2))
}
portfolio_summary <- function(portfolio, dates = NULL, period = 252, benchmark.dates = NULL, benchmark.returns = NULL, plot_stats = FALSE, symbol_wise = FALSE) {
  # if(any(is.na(portfolio)))
  #   stop("portfolio_summary: NAs are not allowed, replace them with zeros")
  if (!is.null(dates) & class(dates) != "Date") {
    stop(paste("portfolio_summary: dates should be of class Date"))
  }
  if (class(portfolio) != "matrix") {
    stop(paste("portfolio_summary: portfolio should be of class matrix"))
  }
  if (is.null(dates)) {
    dates <- seq(as.Date("1970/01/01"), by = "day", length.out = length(returns))
  }
  returns <- rowSums(portfolio, na.rm = TRUE)
  results <- strategy_performance(returns, dates = dates, period = period)
  SRs <- apply(portfolio, 2, function(x) mean(x, na.rm = TRUE) / sd(x, na.rm = TRUE) * sqrt(period))
  results[["SR avg"]] <- round(mean(SRs), 2)
  results[["SR sd"]] <- round(sd(SRs), 2)
  # Benchmark comparison
  alpha <- NA
  beta <- NA
  if (!is.null(benchmark.dates)) {
    benchmark <- data.frame(Dates = benchmark.dates, Returns = benchmark.returns)
    df <- data.frame(Dates = dates, Returns = returns)
    z <- merge(benchmark, df, by = "Dates") %>% na.omit()
    z <- group_by(z, yearmonth(Dates)) %>% summarise(X = sum(Returns.x, na.rm = TRUE), Y = sum(Returns.y, na.rm = TRUE))
    fit <- (lm(Y ~ X, z))
    alpha <- as.numeric(coef(fit)[1]) * 100 * 12
    beta <- as.numeric(coef(fit)[2])
    results[["alpha"]] <- round(alpha, 2)
    results[["beta"]] <- round(beta, 2)
    results[["correlation"]] <- round(cor(z$X, z$Y, use = "pairwise.complete.obs"), 2)
  }
  if (plot_stats) {
    layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
    returns <- replace(returns, is.na(returns), 0) * 100
    cum_returns <- cumsum(returns)
    plot(sort(as.Date(dates)), cum_returns, ylab = "Equity curve log%")
    chunks <- group_by(data.frame(Date = as.Date(dates), ret = returns), year(Date)) %>%
      summarise(sum = round(sum(ret), 1), first = first(Date), .groups = "drop")
    abline(v = chunks$first, lty = 2, lwd = 0.5)
    text(x = chunks$first + period / 2, y = max(cum_returns), labels = chunks$sum, cex = 0.75)
    matplot(apply(portfolio, 2, function(x) cumsum(replace(x, is.na(x), 0))), type = "l", lwd = 2, lty = 1, ylab = "Assets curves log%", xaxt = "n")
    abline(h = 0, lwd = 2)
    axis(side = 1, labels = dates, at = seq(1, length(dates)), tick = FALSE)
    peak <- cummax(cum_returns)
    drawdown <- peak - cum_returns
    plot(sort(as.Date(dates)), -drawdown, ylab = "Drawdown log%", type="l")
  }
  symbol_result <- NULL
  if (symbol_wise) {
    symbol_result <- apply(portfolio, 2, function(x) unlist(strategy_performance(x, dates = dates))) %>% t() %>% as.data.frame
  }
  return(list(Aggregate = results, Symbols = symbol_result))
}
merge_portfolio_list <- function(portfolio_list) {
  full_df <- Reduce(function(...) full_join(..., by = "Date", all = TRUE, incomparables = NA), portfolio_list) %>% arrange(Date)
  colnames(full_df) <- c("Date", names(portfolio_list))
  # full_df[is.na(full_df)] <- 0
  return(full_df)
}
# in percentages
calculate_volatility <- function(returns, long_span=252, short_span=35,  weights=c(0.3, 0.7), period=252){
  vol_short <- sqrt(EMA(replace(returns, is.na(returns), 0)^2, short_span))
  vol_long <- runMean(vol_short, long_span)
  vol <-  (weights[1] * vol_long + weights[2] * vol_short) * sqrt(period) # one year instead of ten
  return(vol)
}
