

# Calculate volatility from log returns
calculate_volatility <- function(returns, long_span=252, short_span=32,  weights=c(0.3, 0.7), period=252){
    if(length(returns) < long_span+short_span)
        return(rep(NA, length(returns)))
    vol_short <- sqrt(EMA(replace(returns, is.na(returns), 0)^2, short_span))
    vol_long <- runMean(vol_short, long_span)
    vol <-  (weights[1] * vol_long + weights[2] * vol_short) * sqrt(period) # one year instead of ten
    return(vol)
}


# ?
normalize_price_ <- function(adjclose, close, volatility, period=252, mult=100) {
    np <- rep(NA, length(close))
    np[1] <- 0
    for(i in 2:length(close)) {
        np[i] <-  (mult * (adjclose[i] - adjclose[i-1]) / (close[i] * volatility[i] / sqrt(period))) + np[i-1]
        if(is.na(np[i]))
            np[i] <- np[i-1]
    }
    return(np)
}
# Rob Carver's ATFS book, strategy 17.
normalize_price <- function(close, volatility, period=1, mult=1) {
    np <- rep(NA, length(close))
    np[1] <- 0
    for(i in 2:length(close)) {
        np[i] <-  (mult * (close[i] - close[i-1]) / (volatility[i] / sqrt(period))) + np[i-1]
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


round_position <- function(position, min_position, position_tick) {
    rounded <- round(position, sapply(position_tick, decimalplaces))
    return(ifelse(abs(rounded) < min_position, 0, rounded))
}

# Running correlation matrix, usually to be run with weekly data.
runCorMatrix <- function(M, n=25, maType=TTR::EMA, absolute_correlation=FALSE) {
    nas <- lapply(1:(n-1), function(i)matrix(NA, nrow = ncol(M), ncol = ncol(M)))
    run_corr_matrices <- lapply(n:nrow(M), function(i)
        c(identity, abs)[[absolute_correlation+1]](cor(M[(i-n):i,], use="pairwise.complete.obs")))
    run_corr_matrices <- c(nas, run_corr_matrices)
    run_corr_vectors <- lapply(run_corr_matrices, as.vector) 
    corr_by_date <- do.call(cbind, run_corr_vectors)
    ema_corr_by_date <- apply(corr_by_date, 1, maType, n) %>% t 
    Q <- lapply(1:ncol(ema_corr_by_date), 
                function(i) matrix(ema_corr_by_date[,i], ncol=ncol(M), dimnames = list(colnames(M), colnames(M))))
    return(Q)
}

# Rob Carver's ATFS book, Appendix B. 
calculate_IDM <- function(returns, weights=NULL, floor_correlation=TRUE) {
    corr <- cor(returns, use = "pairwise.complete.obs")
    if(floor_correlation) 
        corr[corr<0] <- 0
    if(is.null(weights))
        weights <- rep(1/ncol(returns), ncol(returns))
    idm <- 1/(weights %*% corr %*% weights)^0.5
    
    return(as.numeric(idm))
}

calculate_portfolio_volatility <- function(capital, positions, contractsizes, Prices, FXs, cov_matrix) {
    w <- as.numeric(positions * contractsizes * Prices / FXs / capital) 
    return(as.numeric(sqrt(w %*% cov_matrix %*% w)))
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
# min_positions : vector of the min allowed positions (in absolute contracts). Ignored if NULL
# costs_per_contract : vector of the costs to trade one contract, in price scale. 
# trade_shadow_cost : a factor multiplier of the cost per contracts.
# fractional : If your broker allow fractional contracts, like for CFDs, this is the minimal position increment. If you are trading futures where all contracts are 1, just ignore it.
# returned value: a vector of optimized positions according to the dynamic portfolio algorithm.
dynamic_portfolio <- function(capital, optimal_positions, notional_exposures, cov_matrix, 
                              previous_position = NULL, min_positions=NULL, max_positions=NULL, costs_per_contract = NULL, 
                              trade_shadow_cost = 1, fractional=NULL) {
    # Calculate the cost of making trades. trade_shadow_cost represents the number of expected trades in year 
    calculate_costs <- function(weights) {
        trade_gap <- abs(weights_previous - weights)
        trade_costs <- trade_shadow_cost * sum(trade_gap * costs_per_trade_in_weight)
        return(trade_costs)
    }
    # Calculate the error of given weights from the optimal weights considering instruments correlations, plus optional costs
    evaluate <- function(weights_optimal, weights, cov_matrix) {
        solution_gap <- weights - weights_optimal 
        track_error <- as.numeric(sqrt(t(solution_gap) %*% cov_matrix %*% solution_gap))
        if(any(is.nan(track_error)))           
            stop(paste("dynamic_portfolio: NAs is the tracking error"))
        trade_costs <- calculate_costs(weights)
        return(track_error + trade_costs)
    }
    # The greedy algorithm (see https://qoppac.blogspot.com/2021/10/mr-greedy-and-tale-of-minimum-tracking.html)
    find_possible_new_best <- function(weights_optimal, weights_max, weights_min, weights_per_contract, direction, best_solution, best_value, cov_matrix) {
        new_best_value <- best_value
        new_solution <- best_solution
        count_assets <- length(best_solution)
        for (i in sample(1:count_assets)) {
            temp_step <- best_solution
            # The first increment will be the minimum position (as portfolio weight), if provided
            if(temp_step[i] == 0 & weights_min[i] > 0) {
                temp_step[i] <- temp_step[i] + weights_min[i] * direction[i]
            # Otherwise, do the normal greedy increments
            } else {
                temp_step[i] <- temp_step[i] + weights_per_contract[i] * fractional[i] * direction[i]
            }
            # Check if we have exceeded the maximum position (as portfolio weight), if provided
            if(abs(temp_step[i]) > weights_max[i])
                temp_step[i] <- weights_max[i] * sign(temp_step[i])
            # Check we haven't exceed the portfolio capital (weight = 1)
            #if(sum(abs(temp_step)) > 1)
            #    temp_step[i] <- best_solution[i]
            # Evaluate this solution and update the current best
            temp_objective_value <- evaluate(weights_optimal, temp_step, cov_matrix)
            if (temp_objective_value < new_best_value) {
                new_best_value <- temp_objective_value
                new_solution <- temp_step
            }
        }
        return(list(new_best_value, new_solution))
    }
    
    # Number of instruments
    n <- nrow(cov_matrix)
    # Set previous positions as zero if not specified
    if (is.null(previous_position)) 
        previous_position <- rep(0, n)
    # Set trading costs to zero if not specified
    if (is.null(costs_per_contract)) 
        costs_per_contract <- rep(0, n)
    # If fractional is not provided we assume an increment of 1
    if (is.null(fractional)) 
        fractional <- rep(1, n)
    # Calculate contracts weights
    weights_per_contract <- notional_exposures / capital
    weights_optimal <- optimal_positions * weights_per_contract 
    weights_max <- if(!is.null(max_positions)) max_positions * weights_per_contract else rep(Inf, n)
    weights_min <- if(!is.null(min_positions)) min_positions * weights_per_contract else rep(0, n)
    weights_previous <- previous_position * weights_per_contract
    costs_per_trade_in_weight <- (costs_per_contract  / capital) / weights_per_contract
    best_solution <- rep(0,n)
    best_value <- evaluate(weights_optimal, best_solution, cov_matrix)
    while (1) {
        res <- find_possible_new_best(weights_optimal, weights_max, weights_min, weights_per_contract, sign(weights_optimal), best_solution, best_value, cov_matrix)
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



