# Inner functions for time series simulation




#' \code{genpoll} Simulates pollution data with error
#' 
#' @param N Number of days of data
#' @param errtype Error type.  One of berkson, classical, amplitude, or shift
#' @param mean Mean of pollution data on the log scale 
#' @param sd Standard deviation of pollution data on the log scale 
#' @param magerr Magnitude of error.  For Berkson or classical error, this is the standard deviation.  For amplitude, this is the multiplying factor.  For shift, this is the days of shift.
#' @param timetrend Time trends per year (2 for cold/warm, 12 for monthly)
#' @export
genpoll <- function(N = 1000, errtype = "berkson", mean = 1, sd = .1,
                    magerr = 0.01, timetrend = 0, check = F, day = NULL) {

  # sequence of dates
  dates <- seq(1, N)
  
  # seasonal, monthly, and daily trends for true time series
  seas <- ttrend(2, dates, 1)
  month <-  ttrend(12, dates, 1)
  if(is.null(day)) day <- rnorm(N, sd= 1)
  x1 <- day + seas + month
  
  # variance = 1 to add error
  x  <- x1 / (sd(x1))
  # true time series
  xt <- x

  # by error type, create error time series, xe
  if(errtype %in% c("amplitude", "shift")) {

    # seasonal or monthly error
    if(timetrend == 2) {
      change <- seas / sd(x1)
    } else if(timetrend == 12) {
      change <- month / sd(x1)
    }
    # non-error component of xt
    xother <- x - change
    
    # if error type is amplitude change, multiply by magerr for designated frequency
    if(errtype == "amplitude") {

      xe <- xother + change * magerr
      
    # if error type is shift, shift the time series at designated frequency
    } else if(errtype == "shift") {
      xe <- xother + Lag(change, magerr)
    }

  # Berskon error added to designated frequency
  } else if(errtype == "berkson") {
    # simulate random error
    err <- rnorm(N, sd = magerr)
    
    # create error time series
    xe <- (day /sd(x1) + err) / (1 + magerr^2) + (month + seas)/ sd(x1)

  # Classical error added to designated frequency
  } else if(errtype == "classical") {
    # simulate random error
    err <- rnorm(N, sd = magerr)
    
    # create error time series
    xe <- xt + err
  } 
  
  # denormalize:
  zt <- xt * sd + mean
  ze <- xe * sd + mean
  zt <- exp(zt)
  ze <- exp(ze)

  # format output
  zs <- data.frame(truth = zt, error = ze)
  zs <- mutate(zs, date = seq(1, nrow(zs))) %>%
    gather(., type, value, -date)

  # return truth and estimated
  out <- zs
  # return full xt (for checking)
  if(check) out <- list(zs = zs, seas = seas, month = month, day = day, xe = xe, xt = xt)
  return(out)
}



#' \code{ttrend} Creates time trend in simulated pollution data
#' 
#' @title ttrend
#' @param timetrend Time trends per year (2 for cold/warm, 12 for monthly)
#' @param dates Dates (e.g., sequence)
#' @param magerr Standard deviation for time series
#' @export
ttrend <- function(timetrend, dates, magerr) {
  # pi instead of 2 pi so that timetrend corresponds to seasons, etc.
  err <- (1 * cos( pi/ (365 / timetrend)  * dates) )
  
  # standardize (for adding error relative to truth)
  (err - mean(err)) / sd(err) * magerr
}








#' \code{wraptsdecomp} Wrapper function to run tsdecomp
#' 
#' @param dat Dataset with column value indicating time series
#' @param breaks Vector of breaks for decomposition (if null, specify to see seasonal, monthly, acute in simulated data)
#' @export
wraptsdecomp <- function(dat, breaks = NULL) {
  # remove missing data
  dat <- dat[complete.cases(dat), ]
  # arrange by date
  dat <- arrange(dat, date)
  
  # find breaks for time series decomposition
  if(is.null(breaks)) {
    N <- nrow(dat)
    # number of years
    nyears <- ceiling(N/ 365)
    # number of month trends
    nmonth <- nyears * 6
    # maximum
    Ntot <- nyears * 365
    nyears2 <- round(Ntot / 2) + 100
    # breaks
    breaks <- c(1, (nmonth - nyears) / 2 + nyears, (nyears2 - nmonth)/ 2 + nmonth, nyears2)
  }
   
  # get time series decomposition
  ts1 <- tsModel::tsdecomp(dat$value, breaks = breaks)
  ts1 <- as.matrix(ts1)
  
  # format output
  cn1 <-  breaks[-length(breaks)]
  cn2 <- breaks[-1]
  colnames(ts1) <- paste0("b", cn1, ".", cn2)
  dat <- select(dat, -value)
  ts1 <- data.frame(dat, ts1)
  
  # fix range
  ts1 <- gather(ts1, cut, value, -date) %>%
    mutate(cut = gsub("b", "[", cut),
           cut = gsub("\\.", ",", cut),
           cut = paste0(cut, ")")) 
  
  return(ts1)
}


#' \code{lvr} Log variance ratio
#' 
#' @param truth True time series
#' @param error Error time series
#' @export
lvr <- function(truth, error) {
  log(var(error, na.rm = T) / var(truth, na.rm = T))
}







#' \code{healthsim} Simulate health outcome data
#' 
#' @title healthsim
#' @param truth True pollutant and health matrix
#' @param beta0 Base rate of health outcome (intercept)
#' @param beta1 Association (log relative risk)
#' @export
healthsim <- function(truth, beta0, beta1) {
  
  # remove time trends
  lt <- length(truth)
  nyears <- lt / 365
  # 6 df / year, residuals
  seas <- ttrend(2, seq(1, lt), 1)
  month <-  ttrend(12, seq(1, lt), 1)  
  zt <- truth - seas - month
 
  # get mean
  #detrended z:
  df <- ceiling(lt / 365) * 12
  dates <- data.frame(seq(1, lt))
  temp <- lm(zt ~ ns(seq(1, lt), df))
  # residual
  acute <- temp$resid
  #predicted
  long <- predict(temp, dates)
  mu <- exp(beta0 + beta1 * acute + .1 * long )
  # simulate health data
  y <- rpois(length(mu), mu)
  # return results
  y
}







#' \code{fithealth} Fit health models
#' 
#' @title fithealth
#' @param dat Pollutant and health matrix
#' @param rmtimetrend Indicator of whether to remove timetrend
#' @param df Degrees of freedom for natural spline
#' @export
fithealth <- function(dat, rmtimetrend = F, df = 12) {
  
  nyears <- ceiling(nrow(dat) / 365)
  # remove time trend
  if(rmtimetrend) dat$value <- lm(dat$value ~ ns(dat$date, df = df * nyears))$resid
  
  # fit health model
  glm1 <- glm(y ~ value + ns(date * nyears, df = df * nyears), family = "quasipoisson", data = dat) %>%
    tidy() %>%
    filter(., term == "value") %>%
    select(., -term)
  
  # return results
  return(glm1)
}






#' \code{getmeansd} Summarize results from simulation study
#' 
#' @param dat Dataset with columns for meas, val, correlation (cor), lvr
#' @export
getmeansd <- function(dat) {
  gdat <- gather(dat, meas, val, cor, lvr) %>%
    group_by(., cut, error, timetrend, meas, magerr) %>%
    summarize(., mean1 = mean(val), sd1 = sd(val), n1 = n(),
              med1 = median(val), q1 = quantile(val, probs = 0.25),
              q3 = quantile(val, probs = 0.75)) %>%
    mutate(., se= sd1 / sqrt(n1), lb = mean1- 1.96 * se, ub = 1.96 * se)
  gdat
}


#' \code{spfn} Split text
#' 
#' @param x String
#' @export
spfn <- function(x) {
  x1 <- substring(x, 2)
  x1 <- strsplit(x1, ",") %>% sapply(., function(x) x[[1]])
  
  as.numeric(x1)
}



