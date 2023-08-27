library(tidyverse)
library(lubridate)
library(berryFunctions)
library(TTR)
options(stringsAsFactors = FALSE)
Sys.setlocale("LC_TIME", "en_US.UTF-8")


matplot2 <- function(df, ...){
  matplot(df, type="l", lty=1, ...=...)
}


# in percentages
calculate_volatility <- function(returns, long_span=252, short_span=35,  weights=c(0.3, 0.7), period=252){
    vol_short <- sqrt(EMA(replace(returns, is.na(returns), 0)^2, short_span))
    vol_long <- runMean(vol_short, long_span)
    vol <-  (weights[1] * vol_long + weights[2] * vol_short) * sqrt(period) # one year instead of ten
    return(vol)
}


random_ohlc <- function(n=100, m=24, mu=0, sigma=1, lambda=1000) {
    x <- rnorm(n*m, mu, sigma) %>% cumsum
    df <- data.frame(N=rep(1:n, each=m), M=rep(1:m, n), X=x) %>% group_by(N) %>% summarize(Open=first(X), High=max(X), Low=min(X), Close=last(X))
    df <- select(df, -N) %>% mutate(Ticks=rpois(n, lambda))
    return(df)
}



# Load some prices
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


montecarlo_resampler <- function(x, n, f, ...) {
    res <- rep(NA, n)
    for(i in 1:n){
        r <- sample(x, length(x), replace = T)
        res[i] <- f(r, ...)
    }
    return(res)
}






