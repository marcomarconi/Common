

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



instrumentVolatility <- function(returns, period=252) { 
  return(apply(returns, 2, sd, na.rm=TRUE)*sqrt(period))
}


calculateIDM <- function(returns, weights=NULL, absolute_correlation=TRUE) {
  corr <- cor(returns, use = "pairwise.complete.obs")
  if(absolute_correlation) 
    corr <- abs(corr)
  if(is.null(weights))
    weights <- rep(1/ncol(returns), ncol(returns))
  idm <- 1/(weights %*% corr %*% weights)^0.5
  
  return(as.numeric(idm))
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


