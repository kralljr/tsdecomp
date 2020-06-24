# Outer functions for time series simulation




#' \code{wrapsim} Wrap simulation
#' 
#' @title wrapsim
#' 
#' @param nsim Number of simulated datasets
#' @param N Number of days of data
#' @param errtype Error type.  One of berkson, classical, amplitude, or shift
#' @param mean Mean of pollution data on the log scale 
#' @param sd Standard deviation of pollution data on the log scale 
#' @param magerr Magnitude of error.  For Berkson or classical error, this is the standard deviation.  For amplitude, this is the multiplying factor.  For shift, this is the days of shift.
#' @param timetrend Time trends per year (2 for cold/warm, 12 for monthly)
#' @param health Indicator of whether health simulation should be conducted
#' @param beta0 Base rate for health outcome
#' @param rr Relative risk per IQR increase in pollutant
#' @param iqrt IQR for truth.  
#' @export
wrapsim <- function(nsim = 10, N = 365, errtype = "berkson", mean = 1, sd = 2,
                    magerr = 1, zs = NULL,
                    beta0 = 3, rr = 1.01, iqrt = 10, timetrend = NULL,
                    health = F) {
  # get beta1 from RR and IQR
  beta1 <- log(rr) / iqrt
  
  # run pollution simulation and format results
  out <- rerun(nsim, innersim(N = N, errtype = errtype, mean = mean, sd = sd, magerr = magerr, timetrend = timetrend,
                               health = health, beta0 = beta0, rr = rr, iqrt = iqrt))
  zout <- map_dfr(lapply(out, function(x) x[["zout"]]), ~ as.data.frame(.)) 
  
  # format health
  if(health) {
    # add confidence intervals
    
    healthres <- map_dfr(lapply(out, function(x) x[["healthout"]]), ~ as.data.frame(.)) %>%
      mutate(., lb = estimate - 1.96 * std.error, 
             ub = estimate + 1.96 * std.error,
             # find coverage
             coverage = ifelse(lb < beta1 & ub > beta1, 1, 0)) %>%
      group_by(type) %>%
      # find RMSE, bias, coverage
      summarize(rmse = sqrt(mean((estimate - beta1)^2, na.rm = T)),
                meanbias =  mean(estimate - beta1, na.rm = T), 
                coverage = mean(coverage, na.rm = T))
    
    out <- list(zout = zout, healthres = healthres)
  } else {
    out <- list(zout = zout)
  }
  
  # return results
  return(out)
}





#' \code{outersim} Simulation for multiple paramters
#' 
#' @title outersim
#' @param errs Vector of length L of error types.  One of berkson, classical, amplitude, or shift
#' @param mns Vector of length L of mean of pollution data on the log scale 
#' @param sds Vector of length L of standard deviation of pollution data on the log scale 
#' @param magerrs Vector of length L of magnitude of error.  For Berkson or classical error, this is the standard deviation.  For amplitude, this is the multiplying factor.  For shift, this is the days of shift.
#' @param timetrends Vector of length L of time trends per year (2 for cold/warm, 12 for monthly)
#' @param nsim Number of simulated datasets
#' @param N Number of days of data
#' @param health Indicator of whether health simulation should be conducted
#' @param beta0 Base rate for health outcome
#' @param rrs Vector of length L of relative risk per IQR increase in pollutant
#' @param iqrt IQR for truth.  
#' @export
outersim <- function(errs, mns, sds, magerrs, timetrends, nsim = 10, N = 365, 
                     health = F, beta0 = 3, rrs = 1.01, iqrt = 10) {
  
  for(i in 1 : length(errs)) {
    
    # if no health (i.e., rrs is default)
    if(length(rrs) == length(errs)) {
      rrs1 <- rrs[i]
    } else {
      rrs1 <- rrs
    }

    if(length(mns) == 1) {
      mns <- rep(mns, length = length(errs))
      sds <- rep(sds, length = length(errs))
    }
    
    # Run simulation nsim times for this configuration
    res1 <- (wrapsim(nsim, N, errtype = errs[i], mean = mns[i], sd = sds[i],
                     magerr = magerrs[i], beta0 = beta0, 
                     rr = rrs1, iqrt = iqrt, timetrend = timetrends[i], health = health))
    
    # format dataset
    zout1 <- res1$zout %>%
      mutate(., error = errs[i], mean = mns[i], sd = sds[i],
             magerr = magerrs[i], timetrend = timetrends[i])
    
    
    if(i == 1) {
      zout <- zout1
    } else {
      zout <- suppressMessages(full_join(zout, zout1))
    }
    
    # Health only
    if(health) {
      health1 <- mutate(res1$healthres, error = errs[i],
                        mean = mns[i], sd = sds[i],
                        magerr = magerrs[i], rr = rrs1, timetrend = timetrends[i])
      if(i == 1) {
        healthall <- health1
      } else {
        healthall <- full_join(healthall, health1)
      }
    
    # end health
    }  
    

  # end loop  
  }
  
  # format output
  zout <- mutate(zout, order = spfn(cut), 
               cut = reorder(cut, order),
               error = ifelse(error == "amplitude" & magerr > 1, "classical-like",
                              ifelse(error == "amplitude" & magerr < 1, "Berkson-like",
                                     ifelse(error == "berkson", "Berkson-like",
                                            ifelse(error == "classical", "classical-like", error)))),
               cut2 = factor(cut, levels = levels(cut), labels = c("seasonal", "monthly", "acute", "overall")),
               timetrend = factor(timetrend, levels = c(2, 12, 0), labels = c("seasonal", "monthly", "acute"))) %>%
    select(., -order)
  
  
  
  
  # what to return
  if(health ) {
    healthall <- mutate(healthall, error = ifelse(error == "amplitude" & magerr > 1, "classical-like",
                                       ifelse(error == "amplitude" & magerr < 1, "Berkson-like",
                                              ifelse(error == "berkson", "Berkson-like",
                                                     ifelse(error == "classical", "classical-like", error)))),
                        timetrend = factor(timetrend, levels = c(2, 12, 0), labels = c("seasonal", "monthly", "acute")))
    out <- list(zout = zout, healthall = healthall)
  } else {
    out <- list(zout = zout)
  }

  
  return(out)  
}






#' \code{innersim} Run simulation comparing error and true time series
#' 
#' @title innersim
#' @param N Number of days of data
#' @param errtype Error type.  One of berkson, classical, amplitude, or shift
#' @param mean Mean of pollution data on the log scale 
#' @param sd Standard deviation of pollution data on the log scale 
#' @param magerr Magnitude of error.  For Berkson or classical error, this is the standard deviation.  For amplitude, this is the multiplying factor.  For shift, this is the days of shift.
#' @param timetrend Time trends per year (2 for cold/warm, 12 for monthly)
#' @param health Indicator of whether health simulation should be conducted
#' @param beta0 Base rate for health outcome
#' @param rr Relative risk per IQR increase in pollutant
#' @param iqrt IQR for truth.  
#' @export
innersim <- function(N = 1000, errtype = "berkson", mean = 1, sd = 0.1,
                     magerr = 0.01, timetrend = NULL, 
                     health = F, beta0 = 3, rr = 1.01, iqrt = 10) {
  
  # simulate pollution data
  zs <- genpoll(N, errtype = errtype, mean = mean, sd = sd, magerr = magerr, timetrend = timetrend) 

  # do time series decomposition, get summary statistics
  zsts <-  nest(zs, data = c(date, value)) %>%
    mutate(ts = map(data, ~ wraptsdecomp(dat = .))) %>%
    select(., -data) %>%
    unnest(ts) %>%
    spread(., type, value) %>%
    # frequency specific correlation/lvr
    group_by(., cut)  %>%
    summarize(., cor = cor(error, truth, use = "pair"), lvr = lvr(truth, error)) 
  
  # get overall correlation/lvr
  overall <- zs %>%
    spread(., type, value) %>%
    summarize(., cor = cor(truth, error, use = "pair"), lvr = lvr(truth, error),
              rmse = mean((truth - error)^2, na.rm = T)) %>%
    mutate(., cut = "overall")
  

  # add in overall correlation, lvr
  zout <- full_join(zsts, overall)
  
  # if health indicator, run health simulation
  if(health) {
    healthout <- spread(zs, type, value) %>% innerhealth(., zout, beta0, rr, iqrt)
    out <- list(zout = zout, healthout = healthout)
  } else {
    out <- list(zout = zout)
  }
  
  # return results
  return(out)
  
}






#' \code{innerhealth} Run simulaation comparing error and true time series
#' 
#' @title innerhealth
#' @param zs Dataframe of pollutant data
#' @param res1 Results from pollutant error simulation
#' @param truth True simulated pollutant data
#' @param beta0 Base rate for health outcome
#' @param rr Relative risk per IQR increase in pollutant
#' @param iqrt IQR for truth.  
#' @export
innerhealth <- function(zs, res1,
                        beta0 = 3, rr = 1.01, iqrt = 10) {
  
  # find beta1 association from RR and IQR
  beta1 <- log(rr)/ iqrt

  # simulate health data
  dat <- mutate(zs, y = healthsim(truth, beta0 = beta0, beta1 = beta1)) %>%
    gather(., type, value, truth, error)

  # find interquartile range by error type
  iqrs <- group_by(dat, type) %>%
    summarize(iqr = IQR(value, na.rm = T))
  

  # get association for different error types as IQR increase
  betares <- nest(dat, data = c(date, y, value)) %>%
    mutate(fit = map(data, ~ fithealth(dat = .))) %>%
    unnest(cols = c(fit)) %>%
    select(., -data, -statistic, -p.value) %>%
    full_join(., iqrs) %>%
    mutate(., rr = exp(estimate * iqr), cut = "overall")
  
  # results
  return(betares)
  
}