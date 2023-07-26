# some helper functions we'll use throughout

# more stable than log(sum(exp(x))) 
log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

# more stable than log(mean(exp(x)))
log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
}

# compute log of raw importance ratios
# sums over observations *not* over posterior samples
sum_log_ratios <- function(loglik, ids = NULL) {
  if (!is.null(ids)) loglik <- loglik[, ids, drop = FALSE]
  rowSums(loglik)
}

# for printing comparisons later
rbind_print <- function(...) {
  round(rbind(...), digits = 2)
}


plot_ks <- function(ks, ids, thres = 0.6) {
  dat_ks <- data.frame(ks = ks, ids = ids)
  ggplot(dat_ks, aes(x = ids, y = ks)) +
    geom_point(aes(color = ks > thres), shape = 3, show.legend = FALSE) +
    geom_hline(yintercept = thres, linetype = 2, color = "red2") +
    scale_color_manual(values = c("cornflowerblue", "darkblue")) +
    labs(x = "Data point", y = "Pareto k") +
    ylim(-0.5, 1.5)
}


k_thres <- 0.7




LOO_LFO <- function(fit, df, L, M) {
  N <- nrow(df)
  approx_elpds_Msap <- rep(NA, N)
  
  # initialize the process for i = L
  past <- 1:L
  oos <- (L + 1):(L + M)
  df_past <- df[past, , drop = FALSE]
  df_oos <- df[c(past, oos), , drop = FALSE]
  fit_past <- update(fit, newdata = df_past, recompile = FALSE, cores=4)
  loglik <- log_lik(fit_past, newdata = df_oos, oos = oos)
  loglikm <- rowSums(loglik[, oos, drop=F])
  approx_elpds_1sap[L + 1] <- log_mean_exp(loglikm)
  
  # iterate over i > L
  i_refit <- L
  refits <- L
  ks <- NULL
  for (i in (L + 1):(N - M)) {
    past <- 1:i
    oos <- (i + 1):(i + M)
    df_past <- df[past, , drop = FALSE]
    df_oos <- df[c(past, oos), , drop = FALSE]
    loglik <- log_lik(fit_past, newdata = df_oos, oos = oos)
    
    logratio <- sum_log_ratios(loglik, (i_refit + 1):i)
    psis_obj <- suppressWarnings(psis(logratio))
    k <- pareto_k_values(psis_obj)
    ks <- c(ks, k)
    if (k > k_thres) {
      # refit the model based on the first i observations
      i_refit <- i
      refits <- c(refits, i)
      fit_past <- update(fit_past, newdata = df_past, recompile = FALSE)
      loglik <- log_lik(fit_past, newdata = df_oos, oos = oos)
      loglikm <- rowSums(loglik[, oos, drop=F])
      approx_elpds_Msap[i + 1] <- log_mean_exp(loglikm)
    } else {
      lw <- weights(psis_obj, normalize = TRUE)[, 1]
      loglikm <- rowSums(loglik[, oos, drop=F])
      approx_elpds_Msap[i + 1] <- log_sum_exp(lw + loglikm)
    }
  }
  
  return(c("ks"=ks, "approx_elpds_Msap"=approx_elpds_Msap))
}


LOO_LFO_1 <- function(model, y, L, model_data=list() , chains=1, samples=250, SEED=1983) {
  N = length(y)
  approx_elpds_1sap <- rep(NA, N)
  fits <- list()
  m_pred <- rep(NA, N)
  # initialize the process for i = L
  fit_past <- model$sample(data=c(list(N = length(y), y = y, trainset=L ), model_data) , chains=chains, parallel_chains=chains, iter_warmup = samples,iter_sampling = samples, seed = SEED)
  fits[[as.character(L)]] <- fit_past
  mp <- fit_past$draws("m_pred") %>% merge_chains %>% colMeans() %>% as.vector()
  loglik <- fit_past$draws("log_lik") %>% merge_chains()
  loglik <- matrix(loglik, dim(loglik)[1], dim(loglik)[3])
  approx_elpds_1sap[L + 1] <- log_mean_exp(loglik[, oos])

  # iterate over i > L
  i_refit <- L
  refits <- L
  ks <- NULL
  for (i in (L + 1):(N - 1)) {
    print(i)
    logratio <- sum_log_ratios(loglik, (i_refit + 1):i)
    psis_obj <- suppressWarnings(psis(logratio))
    k <- pareto_k_values(psis_obj)
    ks <- c(ks, k)
    if (k > k_thres) {
      # refit the model based on the first i observations
      i_refit <- i
      refits <- c(refits, i)
      fit_past <- model$sample(data=c(list(N = length(y), y = y,  trainset=i  ), model_data), 
                               chains=chains, parallel_chains=chains, iter_warmup = samples, iter_sampling = samples, seed = SEED)
      fits[[as.character(i)]] <- fit_past
      mp <- fit_past$draws("m_pred") %>% merge_chains %>% colMeans() %>% as.vector()
      loglik <- fit_past$draws("log_lik") %>% merge_chains()
      loglik <- matrix(loglik, dim(loglik)[1], dim(loglik)[3])
      approx_elpds_1sap[L + 1] <- log_mean_exp(loglik[, oos])
    } else {
      lw <- weights(psis_obj, normalize = TRUE)[, 1]
      approx_elpds_1sap[i + 1] <- log_sum_exp(lw + loglik[, oos])
    }
    m_pred[i] <- mp[i]
  }
  return(list(approx_elpds_1sap=approx_elpds_1sap, ks=ks, refits=refits, fits=fits, m_pred=m_pred))
}



AR_loglik <- function(y, pars, K){
  N <- length(y)
  S <- nrow(pars$lp__)
  log_lik <- matrix(nrow = S, ncol = N)
  alpha <- pars$alpha
  beta <- pars$beta
  sigma <- pars$sigma
  for (t in 1:K) 
    log_lik[,t] = dnorm(y[t], alpha, sigma, log = TRUE);
  for (t in (K+1):N) {
    mu <- alpha;
    for (k in 1:K)
      mu <- mu + beta[,k] * y[t-k];
    log_lik[,t] = dnorm(y[t], mu, sigma, log = TRUE);
  }
  return (log_lik)
}

LOO_LFO_norefit <- function(fit, y, L, model_data=list() ,SEED=1983) {
  N = length(y)
  pars <- rstan::extract(fit)
  S <- nrow(pars$lp__)
  approx_elpds_1sap <- rep(NA, N)
  past <- 1:L
  oos <- L + 1
  y_past <- y[past]
  y_oos <- y[oos:length(y)]
  loglik <- AR_loglik(y, pars, model_data$K)
  approx_elpds_1sap[L + 1] <- log_mean_exp(loglik[, oos])
  ks <- NULL
  for (i in (L + 1):(N - 1)) {
    past <- 1:i
    oos <- i + 1
    y_past <- y[past]
    y_oos <- y[oos:length(y)]
    logratio <- sum_log_ratios(loglik, (L+1):i)
    psis_obj <- suppressWarnings(psis(logratio))
    k <- pareto_k_values(psis_obj)
    ks <- c(ks, k)
    lw <- weights(psis_obj, normalize = TRUE)[, 1]
    approx_elpds_1sap[i + 1] <- log_sum_exp(lw + loglik[, oos])
  }
  return(list(approx_elpds_1sap=approx_elpds_1sap, ks=ks))
}

