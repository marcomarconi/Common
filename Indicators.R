# Technical indicators

GAS <- function(y) {
  h <- c(0,0)
  err <- c()
  err[1] <- y[1]
  err[2] <- y[2] - y[1]
  for(t in 3:length(y)) {
    err[t] = y[t] - 2*(LaplacesDemon::invlogit(h[t-1])-0.5) * y[t-1];
    h[t] <-  (err[t] * y[t-1]) 
  }
  return(2*(LaplacesDemon::invlogit(h)-0.5))
}


BearsBullsImpulse <- function(df, n = 20, ma=TTR::SMA) {
  mov_avg <- ma(df$Close)  
  Bulls = df$High - mov_avg;
  Bears = df$Low - mov_avg;
  avg = Bears + Bulls;
  buffer1 <- ifelse(avg >= 0, 1, -1)
  buffer2 <- ifelse(avg < 0, 1, -1)
  return((buffer1 - buffer2) / 2)
}

SSL <- function(HLC, n = 10, type="simple") {
  High = HLC[,1]
  Low = HLC[,2]
  Close = HLC[,3]
  hlv <- 0; 
  if(type=="simple"){
    MAHigh <- SMA(High, n = n)
    MALow <- SMA(Low, n = n)
  } else if(type=="exponential") {
    MAHigh <- EMA(High, n = n)
    MALow <- EMA(Low, n = n)
  }  
  SSLd <- c()
  SSLu <- c()
  for(i in 2:nrow(HLC)) { 
    if(is.na(MAHigh[i])) {
      SSLd[i] <- NA; 
      SSLu[i] <- NA; 
      next;
    }
    if(Close[i] > MAHigh[i]) 
      hlv <- 1;
    if(Close[i] < MALow[i]) 
      hlv <- -1;  
    if(hlv==-1) {
      SSLd[i] <- MAHigh[i]; 
      SSLu[i] <- MALow[i]
    }
    else {
      SSLd[i] <- MALow[i]
      SSLu[i] <- MAHigh[i]
    }
  }
  return(SSLu - SSLd)
}

DiDi <- function(x, slow = 20, avg = 8, ma=TTR::SMA) {
  hlv <- 0; 
  MAslow <- ma(x, n = slow)
  MAavg <- ma(x, n = avg)
  DiDi <- 1.0 - MAslow / MAavg
  return(DiDi)
}

Braid <- function(df, MAperiod1 = 3, MAperiod2 = 7, MAperiod3 = 14, ATRPeriod = 14, PipsMinSepPercent = 40, ma=TTR::SMA) {
  ma1 <- ma(df$Close, MAperiod1)
  ma2 <- ma(df$Open, MAperiod2)
  ma3 <- ma(df$Close, MAperiod3)
  maxs <- apply(cbind(ma1, ma2, ma3), 1, max) 
  mins <- apply(cbind(ma1, ma2, ma3), 1, min)
  dif <- maxs - mins
  atr <- TTR::ATR(df[,3:5], ATRPeriod*2-1)[,2]
  filter <- atr * PipsMinSepPercent / 100
  buffer <- sapply(1:length(filter), function(i) tryCatch({if(ma1[i] > ma2[i] & dif[i] > filter[i]) 1 else if (ma1[i] < ma2[i] & dif[i] > filter[i]) -1 else 0}, error=function(cond) 0))
  return(buffer)
}

TII <- function(x, P = 60, ma=TTR::EMA) {
  ma_p <- ma(x, P)
  diff <- x - ma_p
  pos_count <- runSum(diff>0, floor(P/2))
  return(400 * (pos_count) / P - 100)
}

AbsoluteStrength <- function(x, n = 20, pre = 1, post = 1, sma=TTR::EMA) {
  x <- sma(x, pre)
  bulls <- abs(c(0,diff(x))) + x - lag(x)
  bears <- abs(c(0,diff(x))) - x + lag(x)
  avgbulls <- sma(bulls, n)
  avgbears <- sma(bears, n)
  return(sma(avgbulls, post) - sma(avgbears, post))
}

KalmanFilterIndicator <- function(x, sharpness=1, K=1) {
  n <- length(x)
  value <- rep(NA, n)
  distance <- rep(NA, n)
  velocity <- rep(NA, n)
  error <- rep(NA, n)
  value[1] <- x[1]
  velocity[1] <- 0
  distance[1] <- 0
  error[1] <- 0
  for(i in 2:length(x)){
    distance[i] <- x[i] - value[i-1]
    error[i] <- value[i-1] + distance[i] * sqrt(sharpness * K / 100)
    velocity[i] <- velocity[i-1] + distance[i]*K/100
    value[i] <- error[i]+velocity[i]
  }
  return(cbind(value=value, velocity=velocity, distance=distance, error=error))
}

McGinleyDynamicIndicator <- function(x, n=20) {
  md <- rep(NA, length(x))
  md[1] <- x[1]
  for(i in 2:length(x)) {
    md[i] <-  md[i-1] + (x[i] - md[i-1]) / (n * (x[i] / md[i-1])^4)
  }
  return(md)
}

VolatilityRatio <- function(x, n=20, sma=TTR::SMA) {
  return(runSD(x, n) / sma(runSD(x, n),n) - 1 )
}

WaddahAttarExplosion <- function(x, fast=20, slow=40, channel=20, mult=2) {
  macd <- MACD(x, fast, slow)
  bbU <- runMean(x, channel) + mult*runSD(x, channel)
  bbL <- runMean(x, channel) - mult*runSD(x, channel)
  t1 = abs(c(NA,diff(macd[,1])))
  e1 = bbU - bbL
  return(t1 - e1)
}

NormalizedVolume <- function(x, n=20, ma=TTR::SMA) {
  return(x / ma(x,n) - 1 )
}

RVI <- function(df, n=10) {
  value1 <- ((df$Close - df$Open)  + 2*lag(df$Close - df$Open, 1)  + 2*lag(df$Close - df$Open, 2)  + 2*lag(df$Close - df$Open, 3))  / 6
  value2 <- ((df$High - df$Low)  + 2*lag(df$High - df$Low, 1)  + 2*lag(df$High - df$Low, 2)  + 2*lag(df$High - df$Low, 3))  / 6
  rvi <- runSum(value1, n) / runSum(value2, 10)
  rvi_sig <- (rvi + 2 * lag(rvi) + 2*lag(rvi,2) + lag(rvi,2))/6
  return(rvi_sig)
}

SDL <- function(x, n1=20, n2=20, ma=TTR::SMA) {
  return(ROC(ma(x, n1), n2))
}

ER <- function(x, n = 20) {
  return(abs(c(rep(0, n), diff(x, n))) / runSum(abs(c(0,diff(x))), n))
}

ER_indicator <- function(x, n = 20, ma=TTR::SMA ) {
  return(ER(x, n) - ma(ER(x, n), n))
}


FD <- function(x, n = 20) {
  dx2 <- (1/n)^2
  L <- runSum(sqrt(dx2 + abs(c(0, diff(x))) / (runMax(x) - runMin(x))), n)
  return(1 - (log(L)+log(2))/(log(2*n)))
}

FD_indicator <- function(x, n = 20, ma=TTR::SMA ) {
  return(FD(x, n) - ma(FD(x, n), n))
}

# BEWARE we always use the "Close" here, not the adjclose
indicators_tester <- function(df, indi, ...) {
  if(is.null(df$Ticks))
    volume = df$Volume
  else
    volume = df$Ticks
  if(is.null(volume))
    stop("Neither Volume nor Ticks in the data.frame?")
  # Directional
  if(indi == "ADX_dir"){
    a <- ADX(df[,3:5], n = 100, ...=...)
    a <- a[,1]-a[,2]
  } else if (indi == "McGinleyDynamicIndicator"){
    a <- McGinleyDynamicIndicator(df$Close)
  } else if (indi == "KalmanFilterIndicator"){
    a <- KalmanFilterIndicator(df$Close)[,1]
  } else if (indi == "KalmanFilterIndicator_velocity"){
    a <- KalmanFilterIndicator(df$Close, sharpness = 1, K=1)[,2]
  }  else if (indi == "TII"){
    a <- TII(df$Close)
  } else if (indi == "BearsBullsImpulse"){
    a <- BearsBullsImpulse(df)
  } else if (indi == "DiDi"){
    a <- DiDi(df$Close)
  }  else if (indi == "Braid"){
    a <- Braid(df)
  }  else if (indi == "SSL10"){
    df_ <- SSL(df, n = 10)
    a <- df_[,2] - df_[,1]
  }  else if (indi == "SSL5"){
    df_ <- SSL(df, n = 5)
    a <- df_[,2] - df_[,1]
  }  else if (indi == "ALMA"){
    a <- ALMA(df$Close, n=20)
  }  else if (indi == "HMA"){
    a <- HMA(df$Close, n=20)
  }else if (indi == "EVWMA"){
    a <- EVWMA(price = df$Close, volume = volume, n=20)
  } else if (indi == "SMA"){
    a <- SMA(df$Close, n=20)
  }   else if (indi == "SDL"){
    a <- SDL(df$Close, n1 = 10, n2 = 10)
  }  
  else if (indi == "aroon"){
    a <- aroon(df[,3:4],n = 20, ...=...)
    a <- a[,3]
  }else if (indi == "RSI"){
    a <- RSI(df$Close)
    a <- a - 50
  } else if (indi == "RVI"){
    a <- RVI(df, n = 10)
  } 
  else if (indi == "CCI"){
    a <- CCI(df[,3:5],n = 20, ...=...) - 50
  }else if (indi == "CLV"){
    a <- CLV(df[,3:5],  ...=...)
  }else if (indi == "MFI"){
    a <- MFI(df[,3:5],n = 20, volume, ...=...) - 50
  }else if (indi == "ROC"){
    a <- ROC(df[,5],n = 20, ...=...)
  }else if (indi == "ROC_ma"){
    a <- SMA(ROC(df[,5],n = 20, ...=...), 3)
  }else if (indi == "CMF"){
    a <- CMF(df[,3:5],volume = volume, ...=..., n = 20)
  }else if (indi == "CTI"){
    a <- c(rep(NA, 19), CTI(df[,5], ...=...))
  }else if (indi == "CMO"){
    a <- CMO(df[,5], n = 20, ...=...)
  } else if (indi == "DVI"){
    a <- DVI(df[,5],n = 20, ...=...)[,3]-0.5
  }else if (indi == "EMV"){
    a <- EMV(df[,3:4], n=20, volume = volume+1, ...=...)[,2]
  } else if (indi == "KST"){
    a <- KST(df$Close, n = 20, ...=...)[,1]
  } else if (indi == "MACD"){
    a <- MACD(df$Close, ...=...)[,1]
  } else if (indi == "runPercentRank"){
    a <- runPercentRank(df$Close, n = 20, ...=...)-0.5
  }else if (indi == "SAR"){
    a <- (df$Close - SAR(df[,3:4], ...=...)[,1])
  }else if (indi == "stoch"){
    a <- stoch(df[,3:5], ...=...)
    return(a[,2]-a[,3])
  }else if (indi == "TDI"){
    a <- TDI(df$Close, n = 20, ...=...)[,1]
  }else if (indi == "TRIX"){
    a <- TRIX(df$Close, n = 20, ...=...)[,1]
  }else if (indi == "ultimateOscillator"){
    a <- ultimateOscillator(df[,3:5],  ...=...)-50
  }else if(indi == "WPR") {
    a <- WPR(df[,3:5], n = 20, ...=...)-0.5
  }
  # non-Direction
  else if (indi == "WaddahAttarExplosion"){
    a <- WaddahAttarExplosion(df$Close)
  } else if (indi == "VolatilityRatio"){
    a <- VolatilityRatio(df$Close)
  }else if (indi == "NormalizedVolume"){
    a <- NormalizedVolume(volume, 20)
  }else if(indi == "ADX_str"){
    a <- ADX(df[,3:5], n = 100, ...=...)
    #return(ifelse(a[,1] > a[,2], 1, ifelse(a[,1] < a[,2], -1, 0 )  ))
    a <- a[,3]-a[,4]
  } else if (indi == "OBV"){
    a <- OBV(df[,5], volume, ...=...)
    a <- a - SMA(a,20)
  } else if (indi == "chaikinAD"){
    a <- chaikinAD(df[,3:5], volume, ...=...)
    a <- a - SMA(a,20)
  }else if (indi == "williamsAD"){
    a <- williamsAD(df[,3:5], ...=...)
  }else if (indi == "chaikinVolatility"){
    a <- TTR::chaikinVolatility(df[,3:4], ...=...)
  }else if (indi == "volatility"){
    a <- volatility(df[,2:5],  ...=...) 
    a <- a / SMA(a, 20) - 1
  }else if(indi == "REX") {
    TVB <- 3*df$Close-(df$Low+df$Open+df$High)
    rex <- SMA(TVB,n = 5, ...=...)
    signal <- SMA(rex,n = 5, ...=...)
    a <- signal
  }  else if (indi == "SNR"){
    a <- SNR(df[,3:5], n=10, ...=...)
    a <- a - SMA(a, 20) - 1
  }else if (indi == "VHF"){
    a <- VHF(volume, n = 20, ...=...)
  }else if (indi == "VHF_ma"){
    a <- VHF(volume, n = 20, ...=...)
    a <- a-SMA(a, 20)
  }else if (indi == "ATR_ma"){
    a <- ATR(df[,3:5], n = 20)[,2]
    a <- a-SMA(a, 20)
  }else if (indi == "ER"){
    a <- ER_indicator(df$Close, n = 20)
  }else if (indi == "FD"){
    a <- FD_indicator(df$Close, n = 20)
  }else if (indi == "ones"){
    a <- rep(1, length(df$Close))
  }else if (indi == "zeros"){
    a <- rep(0, length(df$Close))
  }else if(indi == "GAS") {
    a <- TTR::SMA(GAS(df$Return*100), 10) / runSD(df$Return, 100)
  }
  
  
  else
    stop("Indicator name not found.")
  return(a)
}
