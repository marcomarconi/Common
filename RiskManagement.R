

# Calculate volatility from log returns
calculate_volatility <- function(returns, long_span=252, short_span=32,  weights=c(0.3, 0.7), period=252){
    vol_short <- sqrt(EMA(replace(returns, is.na(returns), 0)^2, short_span))
    vol_long <- runMean(vol_short, long_span)
    vol <-  (weights[1] * vol_long + weights[2] * vol_short) * sqrt(period) # one year instead of ten
    return(vol)
}


# Rob Carver's ATFS book, strategy 17.
normalize_price <- function(adjclose, close, volatility, period=252) {
    np <- rep(NA, length(close))
    np[1] <- 0
    for(i in 2:length(close)) {
        np[i] <-  (100 * (adjclose[i] - adjclose[i-1]) / (close[i] * volatility[i] / sqrt(period))) + np[i-1]
        if(is.na(np[i]))
            np[i] <- np[i-1]
    }
    return(np)
}

cap_forecast <- function(x, cap=20) {
    return(ifelse(x > cap, cap, ifelse(x < -cap, -cap, x ) ))
}

decimalplaces <- function(x) {
    if ((x %% 1) != 0) {
        nchar(strsplit(sub('0+$', '', format(x, scientific = FALSE)), ".", fixed=TRUE)[[1]][[2]])
    } else {
        return(0)
    }
}

# round_position <- function(position, min_position, position_tick) {
#     return(ifelse(abs(position) <= min_position,  0, round(position, sapply(position_tick, decimalplaces))))
# }

round_position <- function(position, min_position, position_tick) {
    rounded <- round(position, sapply(position_tick, decimalplaces))
    return(ifelse(abs(rounded) < min_position, 0, rounded))
}

# Running correlation matrix, usually to be run with weekly data.
runCorMatrix <- function(M, n=25, absolute_correlation=FALSE) {
    nas <- lapply(1:(n-1), function(i)matrix(NA, nrow = ncol(M), ncol = ncol(M)))
    run_corr_matrices <- lapply(n:nrow(M), function(i)
        c(identity, abs)[[absolute_correlation+1]](cor(M[(i-n):i,], use="pairwise.complete.obs")))
    run_corr_matrices <- c(nas, run_corr_matrices)
    run_corr_vectors <- lapply(run_corr_matrices, as.vector) 
    corr_by_date <- do.call(cbind, run_corr_vectors)
    ema_corr_by_date <- apply(corr_by_date, 1, EMA, n) %>% t 
    Q <- lapply(1:ncol(ema_corr_by_date), 
                function(i) matrix(ema_corr_by_date[,i], ncol=ncol(M), dimnames = list(colnames(M), colnames(M))))
    return(Q)
}

# Rob Carver's ATFS book, Appendix B. I use absolute correlation values instead of capping negatives to zero.
calculateIDM <- function(returns, weights=NULL, absolute_correlation=TRUE) {
    corr <- cor(returns, use = "pairwise.complete.obs")
    if(absolute_correlation) 
        corr <- abs(corr)
    if(is.null(weights))
        weights <- rep(1/ncol(returns), ncol(returns))
    idm <- 1/(weights %*% corr %*% weights)^0.5
    
    return(as.numeric(idm))
}


## Greedy algorithm to find a set of positions closest to the optimal provided. Adapted from "Advanced futures trading strategies (2022)".
## refer to the book for details. 
# capital : your money in your account currency
# optimal_positions : vector of the best un-rounded positions (in contracts) corresponding to your forecast. 
#                     Forecasts are usually from -20 (max short position) to 20. Optimal positions in contracts are calculated as following:
#                     (capital * instrument_weight * IDM  * target_volatility / instrument_volatility  * FX_rate * Forecast/10) / (contract_size * price) 
# notional_exposures : vector of the values in account currency of one contract (usually contract_size * price / FX_rate)
# cov_matrix : covariance matrix of instruments returns, usually calculated from the last 6 months of (daily or weekly) returns.
# previous_position : vector of the previous optimal positions. All zeroes if not provided.
# max_positions : vector of the max allowed positions (in absolute contracts), usually corresponding to a forecast of 20 (see above formula). Ignored if NULL
# min_positions : vector of the min allowed positions (in absolute contracts). If NULL, it is set to the minimum incremental step (1 contract for futures).
# costs_per_contract : vector of the costs to trade one contract, in price scale. 
# trade_shadow_cost : a factor multiplier of the cost per contracts.
# fractional : TRUE is your broker allow fractional contracts, like for CFDs. The algorithm will use the decimal part of the positions as incremental step 
#              in the greedy algorithm. If you are trading futures where all contracts are 1, just set it to FALSE.
# max_factor : maximum multiple of optimal position allowed (e.g. if optimal position == 2 and max_factor == 2, the optimized position will be <= 4). 
#
# returned value: a vector of optimized positions according to the dynamic portfolio algorithm.
dynamic_portfolio <- function(capital, optimal_positions, notional_exposures, cov_matrix, 
                              previous_position = NULL, min_positions=NULL, max_positions=NULL, costs_per_contract = NULL, trade_shadow_cost = 1, fractional=TRUE, max_factor=2) {
    # Calculate the cost of making trades. trade_shadow_cost represents the number of expected trades in year
    calculate_costs <- function(weights) {
        trade_gap <- abs(weights_previous - weights)
        trade_costs <- trade_shadow_cost * sum(trade_gap * costs_per_trade_in_weight)
        return(trade_costs)
    }
    # Calculate the error of given weights from the optimal weights considering instruments correlations, plus optional costs
    evaluate <- function(weights_optimal, weights, cov_matrix) {
        solution_gap <- weights_optimal - weights
        track_error <- as.numeric(sqrt(t(solution_gap) %*% cov_matrix %*% solution_gap))
        trade_costs <- calculate_costs(weights)
        return(track_error + trade_costs)
    }
    # The greedy algorithm (see https://qoppac.blogspot.com/2021/10/mr-greedy-and-tale-of-minimum-tracking.html)
    find_possible_new_best <- function(weights_optimal, weights_max, weights_per_contract, direction, best_solution, best_value, cov_matrix, max_factor, buffer) {
        new_best_value <- best_value
        new_solution <- best_solution
        count_assets <- length(best_solution)
        for (i in sample(1:count_assets)) {
            temp_step <- best_solution
            if(temp_step[i] == 0) {
                temp_step[i] <- temp_step[i] + weights_min[i] * direction[i]
            } else {
                temp_step[i] <- temp_step[i] + weights_per_contract[i] * fractional[i] * direction[i]
            }
            if(abs(temp_step[i]) > weights_max[i])
                temp_step[i] <- weights_max[i] * sign(temp_step[i])
            else if (abs(temp_step[i]) > max_factor * abs(weights_optimal[i]))
                temp_step[i] <- max_factor * weights_optimal[i]
            temp_objective_value <- evaluate(weights_optimal, temp_step, cov_matrix)
            if (temp_objective_value < new_best_value) {
                new_best_value <- temp_objective_value
                new_solution <- temp_step
            }
        }
        return(list(new_best_value, new_solution))
    }
    
    # Number os instruments
    n <- nrow(cov_matrix)
    # Set previous positions as zero if not specified
    if (is.null(previous_position)) {
        previous_position <- rep(0, n)
    }
    # Set trading costs to zero if not specified
    if (is.null(costs_per_contract)) {
        costs_per_contract <- rep(0, n)
    }
    # Find a fractional increments from positions (e.g. if position == 1.2 then the increment is 0.1)
    if (!fractional) {
        fractional <- rep(1, n)
    } else {
        fractional <-  10^(floor(log10(abs(optimal_positions)))-1)
    }
    weights_per_contract <- notional_exposures / capital
    weights_optimal <- optimal_positions * weights_per_contract 
    weights_max <- if(!is.null(max_positions)) max_positions * weights_per_contract else rep(Inf, n)
    weights_min <- if(!is.null(min_positions)) min_positions * weights_per_contract else weights_per_contract * fractional
    weights_previous <- previous_position * weights_per_contract
    costs_per_trade_in_weight <- (costs_per_contract  / capital) / weights_per_contract
    best_solution <- rep(0,n)
    best_value <- evaluate(weights_optimal, best_solution, cov_matrix)
    while (1) {
        res <- find_possible_new_best(weights_optimal, weights_max, weights_per_contract, sign(weights_optimal), best_solution, best_value, cov_matrix, max_factor)
        new_best_value <- res[[1]]
        new_solution <- res[[2]]
        if (new_best_value < best_value) {
            best_value <- new_best_value
            best_solution <- new_solution
        } else {
            break
        }
    }
    return(best_solution / weights_per_contract)
}

## Dynamic portfolio buffering.  Adapted from "Advanced futures trading strategies (2022)".
# capital : your money in your account currency
# optimized_positions : vector of optimized positions returned from the function "dynamic_portfolio"
# previous_position : vector of the previous optimal positions. All zeroes if not provided.
# notional_exposures : vector of the values in account currency of one contract (usually contract_size * price / FX_rate)
# cov_matrix : covariance matrix of instruments returns, usually calculated from the last 6 months of (daily or weekly) returns.
# target_volatility : your portfolio volatility target (e.g. 0.25)
# portfolio_buffering_level : the deviance representing the edges of the buffering. The highest this number the less frequent the portfolio updates.
#
# returned value: a list of: a vector of required positions updates to take (all zero if the adjustment factor is negative), 
#                            the tracking error of the portfolio 
#                            the adjustment factor (the percentage of position to adjust from the current to the optimized position)  
buffering_portfolio <- function(capital, optimized_positions, previous_positions, notional_exposures, cov_matrix, target_volatility, portfolio_buffering_level=0.1) {
    weights_per_contract <- notional_exposures / capital
    optimized_portfolio_weight <- optimized_positions * weights_per_contract 
    previous_portfolio_weight <- previous_positions * weights_per_contract 
    tracking_error_current_weight <- optimized_portfolio_weight - previous_portfolio_weight
    tracking_error <- as.numeric(sqrt(t(tracking_error_current_weight) %*% cov_matrix %*% tracking_error_current_weight))
    adjustment_factor <- max((tracking_error - portfolio_buffering_level/2 * target_volatility) / tracking_error, 0)
    required_trade <- adjustment_factor * (optimized_positions - previous_positions) 
    return(list(required_trade, tracking_error, adjustment_factor))
}



## Pardo's HRP
{
# Helper function to calculate the HRP
HRP <- function(returns, absolute_correlation=FALSE) {
    if(any(is.na(returns)))
        stop("HRP: returns matrix cannot contain NAs")
    cov_matrix <-  cov(returns, use = "pairwise.complete.obs")
    corr_matrix <-  cor(returns, use = "pairwise.complete.obs")
    if(absolute_correlation)
        corr_matrix <-  abs(corr_matrix)
    return(setNames(getHRP(cov_matrix, corr_matrix)$w, colnames(returns)))
}

# Helper function to calculate running  HRP
runHRP <- function(returns, lookback=100, absolute_correlation=FALSE) {
    weights <- matrix(0, nrow(returns), ncol(returns)); 
    failed <- 0
    for(i in (lookback+1):nrow(returns)) {  
        sl <-   returns[(i-lookback):i,]
        if(any(colSums(sl)==0) | any(is.na(sl))) {
            weights[i,] <- 1/ncol(returns)
            failed <- failed + 1
        } else {
            weights[i,] <- HRP(sl, absolute_correlation = absolute_correlation)
        }
    }    
    warning(paste("runHRP: HRP calculation failed for",failed,"rows." ))
    return(weights)
}

# HRP algorithm as described by Pardo
getHRP <- function(cov, corr, max = NULL, min = NULL, return_raw = NULL, robust_cov = F) {
    # Construct a hierarchical portfolio
    if (robust_cov == T) {
        if(is.null(return_raw))
            stop("You must provide a raw returns if you want to shrink the covariance matrix")
        cov <- cov.shrink(return_raw, verbose = FALSE)
        corr <- cov2cor(cov)
    }
    
    # Set the constraint matrix
    if (is.null(max)) max <- rep(1,ncol(cov))
    if (is.null(min)) min <- rep(0,ncol(cov))
    if (length(max)==1) max <- rep(max,ncol(cov)) else if (length(max)<ncol(cov)) stop("Provide correct weights")
    if (length(min)==1) min <- rep(min,ncol(cov)) else if (length(min)<ncol(cov)) stop("Provide correct weights")
    const <- rbind(max, min)
    
    # check constraints
    if (sum(const[1,]) < 1 | sum(const[2,]) > 1) stop("Incompatible weights")
    
    distmat <- ((1 - corr) / 2)^0.5
    clustOrder <- hclust(dist(distmat), method = 'mcquitty')$order
    out <- getRecBipart(cov, clustOrder, const)
    return(out)
}

getClusterVar <- function(cov, cItems) {
    # compute cluster variance from the inverse variance portfolio above
    covSlice <- cov[cItems, cItems]
    weights <- getIVP(covSlice)
    cVar <- t(weights) %*% as.matrix(covSlice) %*% weights
    return(cVar)
}

getRecBipart <- function(cov, sortIx, const) {
    
    w <- rep(1, ncol(cov))
    
    # create recursion function within parent function to avoid use of globalenv
    recurFun <- function(cov, sortIx, const) {
        # get first half of sortIx which is a cluster order
        subIdx <- 1:trunc(length(sortIx)/2)
        
        # subdivide ordering into first half and second half
        cItems0 <- sortIx[subIdx]
        cItems1 <- sortIx[-subIdx]
        
        # compute cluster variances of covariance matrices indexed
        # on first half and second half of ordering
        cVar0 <- getClusterVar(cov, cItems0)
        cVar1 <- getClusterVar(cov, cItems1)
        alpha <- 1 - cVar0/(cVar0 + cVar1)
        
        # determining whether weight constraint binds
        alpha <- min(sum(const[1,cItems0]) / w[cItems0[1]],
                     max(sum(const[2,cItems0]) / w[cItems0[1]], 
                         alpha))
        alpha <- 1 - min(sum(const[1,cItems1]) / w[cItems1[1]], 
                         max(sum(const[2,cItems1]) / w[cItems1[1]], 
                             1 - alpha))
        
        w[cItems0] <<- w[cItems0] * rep(alpha, length(cItems0))
        w[cItems1] <<- w[cItems1] * rep((1-alpha), length(cItems1))
        
        # rerun the function on a half if the length of that half is greater than 1
        if(length(cItems0) > 1) {
            recurFun(cov, cItems0, const)
        }
        if(length(cItems1) > 1) {
            recurFun(cov, cItems1, const)
        }
        
    }
    
    # run recursion function
    recurFun(cov, sortIx, const)
    return(list(w=w,sortIx=sortIx))
}

getIVP <- function(covMat) {
    invDiag <- 1/diag(as.matrix(covMat))
    weights <- invDiag/sum(invDiag)
    return(weights)
}
}


# Get returns weight by considering SD and correlations TRY TO IMPLEMENT THE FORMULA IN LEVERAGED TRADER
get_portofolio_weights <- function(returns, SD = TRUE, CORR = FALSE,  HRP = FALSE) {
  if(any(is.na(returns)))
    stop("NAs in the returns matrix")
  if(any(is.numeric(returns)))
    stop("Non-numbers in the returns matrix (have you removed the Date column?)")
  if(SD) {
    sds <- apply(returns, 2, sd, na.rm=TRUE)
    sds <- 1/sds
    sds_weights <- t(replicate(nrow(returns), sds/sum(sds)))
  } else
    sds_weights <- t(replicate(nrow(returns), rep(1, ncol(returns))))
  if(CORR) {
    if(!HRP) {
      corr_weights <- 1-(abs(cor(returns)) %>% colMeans())
      corr_weights <- corr_weights/sum(corr_weights)
      corr_weights <- matrix(corr_weights, nrow(returns), ncol(returns), byrow = TRUE) 
    } else
      corr_weights <- runHRP(returns, lookback = 50)
  } else
    corr_weights <- t(replicate(nrow(returns), rep(1, ncol(returns))))
  weights <- sds_weights * corr_weights
  weights[is.na(weights)] <- 0
  weights <- t(apply(weights, 1, function(x) if(sum(x)==0) setNames(rep(0, length(x)), colnames(weights)) else x/sum(x)))
  return(weights)
}





combineStrategies <- function(Equities) {
  dates <- sort(as.Date(unique(unlist((lapply(Equities, function(x) lapply(x, function(y) y$Date)))))))
  symbols <- sort(unique(unlist((lapply(Equities, function(x) lapply(x, function(y) y$Symbol))))))
  trades <- matrix(0, ncol=length(symbols), nrow=length(dates))
  colnames(trades) <- symbols
  rownames(trades) <- as.character(dates)
  for(n in names(Equities)) {
    for(s in names(Equities[[n]])){
      df <- Equities[[n]][[s]]
      for(i in 1:nrow(df)){
        trades[df$Date[i], s] <- trades[df$Date[i], s]  + df$Trade[i]
      }
    }
  }
  joint <- reshape2::melt(trades) %>% setnames(old = colnames(.), new = c("Date", "Symbol", "Trade"))
  joint$Date <- as.character(joint$Date )
  joint$Symbol <- as.character(joint$Symbol )
  return(joint)

}


