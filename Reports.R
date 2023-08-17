library(tidyverse)
library(Rfast)

# functions to calculare report statistics

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
  daily_returns <- df$Returns
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
  ann_sd <- sd(daily_returns, na.rm = TRUE) * sqrt(period) * 100
  sr <- mean(daily_returns, na.rm = TRUE) / sd(daily_returns, na.rm = TRUE) * sqrt(period)
  skew_ <- skew(weekly_returns[weekly_returns!=0])
  kurtosis_ <- kurt(weekly_returns[weekly_returns!=0])
  demeaned_returns <- daily_returns[daily_returns!=0] - mean(daily_returns[daily_returns!=0])
  q <- quantile(demeaned_returns, probs = c(0.01, 0.3, 0.7, 0.99), na.rm = TRUE)
  lower_tail <- as.numeric(q[1] / q[2] / 4.43)
  upper_tail <- as.numeric(q[4] / q[3] / 4.43)
  cum_returns <- cumsum(replace(daily_returns, is.na(daily_returns), 0))
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


# Get statistics out of a single series of log returns and trades
get_trades_statistics <- function(log_returns, trades, ID, period=252){
    stats = list()  
    if(sd(log_returns) > 0.1)
        warning("log_returns std > 0.1, maybe they are in percentages?")
    if(!all(trades == 1 | trades == -1 | trades == 0))
        stop("trades should be one of these value: -1 0 1")
    # tradesID <- rep(1:length(rle(trades)$length), times=rle(trades)$length)
    # tradesID[trades==0] <- NA
    all_t <- group_by(data.frame(ID, log_returns, trades), ID) %>% na.omit %>% summarise(type = mean(trades), days = n(), return = sum(log_returns), .groups = 'drop') %>% na.omit()
    win_t <- dplyr::filter(all_t, return > 0)
    lose_t <- dplyr::filter(all_t, return < 0)
    stats[['total_trades']] <- nrow(all_t)
    stats[['trades_freq']] <- nrow(all_t) / (length(trades) / period)
    stats[['sharpe_ratio']] <- mean(log_returns) / sd(log_returns) * sqrt(period)
    stats[['mean_return_%']] <- mean(log_returns) * 100
    stats[['volatility']] <- sd(log_returns) * 100
    stats[['profit_factor']]  <- sum(win_t$return) / abs(sum(lose_t$return))
    stats[['time_in_trade_%']] <- (sum(trades != 0) / length(trades)) * 100
    stats[['avg_trades_len']] <- mean(all_t$days)
    stats[['avg_win_len']]  <- mean(win_t$days)
    stats[['avg_lost_len']]  <- mean(lose_t$days)
    stats[['win_trades_%']]  <- nrow(win_t) / nrow(all_t) * 100
    stats[['win_long_%']]  <- nrow(dplyr::filter(all_t, type > 0 & return > 0)) / nrow(dplyr::filter(all_t, type > 0))    * 100
    stats[['win_short_%']]  <- nrow(dplyr::filter(all_t, type < 0 & return > 0)) / nrow(dplyr::filter(all_t, type < 0))    * 100
    stats[['avg_trade_return_%']]  <- mean(all_t$return) * 100
    stats[['avg_trade_win_return_%']]  <- mean(win_t$return) * 100
    stats[['avg_trade_lost_return_%']]  <- mean(lose_t$return) * 100
    stats[['avg_long_return_%']]  <- mean(pull(dplyr::filter(all_t, type > 0), return)) * 100
    stats[['avg_short_return_%']]  <- mean(pull(dplyr::filter(all_t, type < 0), return)) * 100
    stats[['R2']] <- summary(lm(1:length(log_returns) ~ 0+cumsum(log_returns)))$adj.r.squared
    # fit <- bayesglm((log_returns*100) ~ 1, data = data.frame(log_returns), family = "gaussian")
    # stats[['est_mean_return_%']]  <- as.numeric(fit$coefficients) * 100
    # stats[['est_stderr_return_%']]  <- as.numeric(sqrt(diag(vcov(fit)))) * 100
    return(stats)
}

# Get statistics out of a series of log returns
get_returns_statistics <- function(portfolio, dates=NULL, risk_free_rate = 0.0, trading_cost = NULL, holding_cost = NULL, period=252,  plot=FALSE) {
    if(any(is.na(portfolio)))
        stop("Na are not allowed, replace them with zeros")
    log_returns <- rowSums(portfolio)
    statistics = list()  
    if(sd(log_returns) > 0.1)
        warning("log_returns std > 0.1, maybe they portfolio is in percentages?")
    sds <- apply(portfolio, 2, sd)
    statistics[['portfolio_variance_%']] <- sqrt(t(sds) %*% cor(portfolio) %*% sds) * 100
    statistics[['cumulative_returns_%']] <- (exp(sum(log_returns)) - 1) * 100
    statistics[['annualized_returns_%']] <- mean(log_returns ) * period * 100
    statistics[['annualized_volatility_%']] <- sd(log_returns ) * sqrt(period) * 100
    downside <- sd(log_returns[log_returns<0]) 
    statistics[['adj_sortino_ratio']] <- (mean(log_returns) - risk_free_rate) / (downside * sqrt(2)) * sqrt(period) 
    statistics[['sharpe_ratio']] <- (mean(log_returns) - risk_free_rate) / sd(log_returns) * sqrt(period)
    # if(!is.null(dates)) { # monthly version of GPR
    #   monthly <- group_by(data.frame(log_returns=log_returns,  Month=month(as_date(dates)), Year=year(as_date(dates)))  %>% unite(col=Date, Month, Year), Date) %>% summarize(Sum=sum(log_returns)) 
    #   statistics[['GPR']] <- sum(monthly$Sum) / sum(abs(monthly$Sum[monthly$Sum<0]))
    # }
    statistics[['GPR']] <- sum(log_returns) / sum(abs(log_returns[log_returns<0]))
    q <- quantile(log_returns, probs=c(0.1, 0.9))
    statistics[['tail_ratio']] <-  mean(log_returns[log_returns>q[2]]) / abs(mean(log_returns[log_returns<q[1]]))
    cum_returns = cumsum(log_returns) 
    peak = cummax(cum_returns)
    drawdown = peak - cum_returns
    max_idx = which.max(drawdown)
    statistics[['max_drawdown_%']] <- as.numeric((1 - exp(cum_returns[max_idx]) / exp(peak[max_idx])) * 100)
    statistics[['R2']] <- summary(lm(1:length(cum_returns) ~ 0+cum_returns))$adj.r.squared
    if(plot){
        par(mfrow=c(2,2), mar=c(2,4,1,0.5))
        if(is.null(dates)) dates <- seq(as.Date("1970/01/01"), by = "day", length.out = length(log_returns))
        plot(sort(as.Date(dates)), cumsum(log_returns)*100, ylab="Equity curve %", type="o", pch=16)
        chunks <- group_by(data.frame(Date=as.Date(dates), ret=log_returns), year(Date)) %>% 
            summarise(sum=round(sum(ret)*100,1), first=first(Date), .groups = 'drop')
        abline(v=chunks$first, lty=2, lwd=0.5)
        text(x=chunks$first+period/2, y=max(cumsum(log_returns*100)), labels=chunks$sum)
        matplot(apply(portfolio, 2, cumsum)*100, type = "l", lwd=2, lty=1, ylab="Assets curves %", xaxt="n")
        axis(side = 1, labels=dates, at=seq(1, length(dates)), tick = FALSE)
        sr <- apply(portfolio, 2, mean, na.rm=TRUE) / apply(portfolio, 2, sd, na.rm=TRUE) * sqrt(period)
        smallPlot({plot(x=jitter(rep(0, length(sr))), sr, xlim=c(-0.1,0.1), pch=16, xaxt="n", xlab="", ylab="", main="Sharpe Rs")}, x1=0.15, x2=0.2, y1=0.6, y2=0.9, mar=0)
        r2s <- apply(portfolio, 2, function(y){fit <- lm(1:length(y) ~ 0+cumsum(y)); sign(coef(fit)[1]) * summary(fit)$adj.r.squared})
        smallPlot({plot(x=jitter(rep(0, length(r2s))), r2s, xlim=c(-0.1,0.1),ylim=c(-1,1), pch=16, xaxt="n", xlab="", ylab="", main="R2s")}, x1=0.25, x2=0.30, y1=0.6, y2=0.9, mar=0)
        plot(sort(as.Date(dates)), -drawdown*100, type="l", ylab="Drawdown %")
        plot(sort(as.Date(dates)), runSD(log_returns*100, 20), type="l", ylab="Volatility %")
        # log_returns_zeroless <- log_returns[log_returns!=0]
        # hist(log_returns_zeroless, 100, freq = FALSE); 
        # mtext(paste("Skewness:", round(skewness(log_returns_zeroless), 2)), side = 3, line = -2.5, adj = 0.1, cex = 1.5)
        # mtext(paste("Kurtosis:", round(kurtosis(log_returns_zeroless), 2)), side = 3, line = -4.5, adj = 0.1, cex = 1.5)
    }
    return(statistics)
}

