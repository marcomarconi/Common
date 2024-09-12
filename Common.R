library(tidyverse)
library(lubridate)
library(berryFunctions)
library(TTR)
options(stringsAsFactors = FALSE)
Sys.setlocale("LC_TIME", "en_US.UTF-8")


matplot2 <- function(df, ...){
  matplot(df, type="l", lty=1, ...=...)
}


runZscore <- function(x, n=10) {
  return((x-runMean(x, n))/runSD(x, n))
}



run_ntile <- function(x, k=length(x), n) {
  rv <- runquantile(x, k, probs = seq(0, 1, length.out=n+1), align = "right")
  q <- sapply(1:(ncol(rv)-1), function(j) ifelse(x >= rv[,j] & x < rv[,j+1] , j, 0))
  q[,(ncol(rv)-1)] <- ifelse(x == rv[ ,ncol(rv)], ncol(rv)-1, q[,(ncol(rv)-1)])
  fc <-  as.integer(rowSums(q))
  return(fc)
  
}

montecarlo_resampler <- function(x, n, f, m=1, ...) {
    return(replicate(n, f(sample(x, max(m, sample.int(length(x), 1)), replace = T), ...)))
}



random_ohlc <- function(n=252, m=1440, mu=0, sigma=0.001, lambda=1000) {
    x <- rnorm(n*m, mu, sigma) %>% cumsum %>% exp
    df <- data.frame(Time=rep(1:n, each=m), X=x) %>% group_by(Time) %>% summarize(Open=first(X), High=max(X), Low=min(X), Close=last(X))
    df <- mutate(df, Ticks=rpois(n, lambda))
    return(df)
}

gbm_vec <- function(nsim = 1, t = 365, mu = 0, sigma = 0.3, S0 = 100, dt = 1./365) {
    # matrix of random draws - one for each day for each simulation
    epsilon <- matrix(rnorm((t-1)*nsim), ncol = nsim, nrow = t-1)  
    # get GBM and convert to price paths
    gbm <- exp((mu - sigma * sigma / 2) * dt + sigma * epsilon * sqrt(dt))
    gbm <- apply(rbind(rep(S0, nsim), gbm), 2, cumprod)
    return(gbm)
}

# Load some prices OLD
load_all_data <- function(dir="/home/marco/trading/Historical Data/Yahoo/Scraping/", asset_file="/home/marco/trading/Historical Data/Yahoo/Scraping/Assets.txt", scraper_file="/home/marco/trading/Historical Data/Yahoo/retrieve.sh", download_investing=TRUE, download_yahoo=TRUE) {
  
  if(download_investing) {
    print("Loading investing commodities data")
    setwd("/home/marco/trading/Historical Data/Investing//Commodities/")
    unlink("*.csv")
    system("python3 ../investing.py commodities")
    {
      Symbols <- list()
      dir <- "/home/marco/trading/Historical Data/Investing/Commodities/"
      for(f in list.files(dir, pattern = ".csv")){
        print(f)
        a <- read.csv(paste0(dir, f)) 
        if(length(unique(a$Date)) != nrow(a))
          stop(paste0("Duplicate Dates in ", f))
        a <- zoo::na.locf(a, na.rm = FALSE)
        f <- sub("\\.csv", "", f)
        a$Close <- ifelse(a$Close > 0, a$Close, 0)
        a$Return <- c(0, diff(log(a$Close)))
        a$Symbol <- f
        Symbols[[f]] <-  dplyr::filter(a, !(lubridate::wday(lubridate::as_date(Date)) %in% c(1, 7))) 
      }
      # put the merged returns by date on a matrix
      df = Reduce(function(...) full_join(..., by="Date"), Symbols)
      Returns <- df[,grep("^Date|Return", colnames(df))]
      colnames(Returns) <- c("Date", names(Symbols))
      Returns$Date <- as.Date(Returns$Date)
      Returns <- arrange(Returns, Date)
      Closes <- df[,grep("^Date|^Close", colnames(df))]
      colnames(Closes) <- c("Date", names(Symbols))
      Closes$Date <- as.Date(Closes$Date)
      Closes <- arrange(Closes, Date)
      Returns_investing_commodities <<- Returns
      Closes_investing_commodities <<- Closes
    }
    
    
    print("Loading investing CMC data")
    # Load some prices (Investing CMC)
    dir <- "/home/marco/trading/Historical Data/Investing/CMC/"
    setwd(dir)
    unlink("*.csv")
    if(download_investing)
      system("python3 ../investing.py CMC")
    {
      Symbols <- list()
      for(f in list.files(dir, pattern = ".csv")){
        print(f)
        a <- read.csv(paste0(dir, f)) 
        if(length(unique(a$Date)) != nrow(a))
          stop(paste0("Duplicate Dates in ", f))
        a <- zoo::na.locf(a, na.rm = FALSE)
        f <- sub("\\.csv", "", f)
        a$Return <- c(0, diff(log(a$Close)))
        a$Symbol <- f
        Symbols[[f]] <-  dplyr::filter(a, !(lubridate::wday(lubridate::as_date(Date)) %in% c(1, 7))) 
      }
      # put the merged returns by date on a matrix
      df = Reduce(function(...) full_join(..., by="Date"), Symbols)
      Returns <- df[,grep("^Date|Return", colnames(df))]
      colnames(Returns) <- c("Date", names(Symbols))
      Returns$Date <- as.Date(Returns$Date)
      Returns <- arrange(Returns, Date)
      Closes <- df[,grep("^Date|^Close", colnames(df))]
      colnames(Closes) <- c("Date", names(Symbols))
      Closes$Date <- as.Date(Closes$Date)
      Closes <- arrange(Closes, Date)
      Returns_investing_CMC <<- Returns
      Closes_investing_CMC <<- Closes
      }
    
    Returns_investing <<- merge(Returns_investing_CMC, Returns_investing_commodities, by="Date")
    Closes_investing <<- merge(Closes_investing_CMC, Closes_investing_commodities, by="Date")
  }
  
  
  if(download_yahoo) {
    print("Loading yahoo data")
    # Load some prices (Yahoo)
    setwd(dir)
    unlink("*.csv")
    system(paste0("bash ", "\"", scraper_file, "\" ", "\"", asset_file,"\""))
    {
      Symbols <- list()
      for(f in list.files(dir, pattern = ".csv")){
        print(f)
        a <- read.csv(paste0(dir, f)) 
        if(length(unique(a$Date)) != nrow(a))
          stop(paste0("Duplicate Dates in ", f))
        a <- zoo::na.locf(a, na.rm = FALSE)
        f <- sub("\\.csv", "", f)
        a$Close <- ifelse(a$Close > 0, a$Close, 0)
        a$Return <- c(0, diff(log(a$Close)))
        a$Symbol <- f
        Symbols[[f]] <-  dplyr::filter(a, !(lubridate::wday(lubridate::as_date(Date)) %in% c(1, 7))) %>% arrange(as.Date(Date))
      }
      # put the merged returns by date on a matrix
      df = Reduce(function(...) full_join(..., by="Date"), Symbols)
      Returns <- df[,grep("^Date|Return", colnames(df))]
      colnames(Returns) <- c("Date", names(Symbols))
      Returns$Date <- as.Date(Returns$Date)
      Returns <- arrange(Returns, Date)
      Closes <- df[,grep("^Date|^Close", colnames(df))]
      colnames(Closes) <- c("Date", names(Symbols))
      Closes$Date <- as.Date(Closes$Date)
      Closes <- arrange(Closes, Date)
      Returns_yahoo <<- Returns
      Closes_yahoo <<- Closes
    }
  }
}






