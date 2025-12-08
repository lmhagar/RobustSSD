## this code file is used to reproduce the results from Section 5.1

## load necessary packages
require(foreach)
require(doParallel)
require(doSNOW)
require(survival)
require(ggplot2)
require(ggpubr)
require(cowplot)

expit <- function(x){1/(1 + exp(-x))}
logit <- function(x){log(x) - log(1-x)}

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

m <- 10000
registerDoSNOW(cl)
pb <- txtProgressBar(max = m, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## this function generates data with a Weibull baseline hazard and implements censoring;
## the inputs are described below
## n: sample size (both groups combined)
## age_mean and age_sd: parameters to generate age covariate
## volume_prop: proportion of high-volume patients (binary)
## meta_prop0: proportion of zero-inflation Possion model for # of bone metastates
## meta_pois: Possion rate for non-zero component of ZIP model
## ALP_mean and ALP_sd: mean and sd for lognormal distribution of Alkaline Phosphatase
## primary_mean and primary_sd: mean and sd for lognormal distribution of # months since primary cancer diagnosed
## psa_mean and psa_sd: mean and sd for lognormal distribution of prostate-specific antigen
## prop_coef: coefficients for propensity score (logistic regression)
## fail_coef: regression coefficients for survival time model
## fail_shape and fail_scale: parameters for Weibull baseline hazard of survival model
## cens_coef: regression coefficients for censoring time model
## cens_shape and cens_scale: parameters for Weibull baseline hazard of censoring model
gen_weibull <- function(n, age_mean = 69, age_sd = 8, volume_prop = 0.55, meta_prop0 = 0.12,
                        meta_pois = 5.5, ALP_mean = 4.5, ALP_sd = sqrt(2*(log(190) - 4.5)),
                        primary_mean = 3.05, primary_sd = sqrt(2*(log(25) - 3.05)),
                        psa_mean = 2.65, psa_sd = sqrt(2*(log(15) - 2.65)),
                        prop_coef = c(2.6, -0.06, 1, 0.24, 0.34, -0.07, 0.01),
                        fail_coef = c(-0.475, 0.2, 0.01, 0.125, 0.03),
                        fail_shape = 1.5, fail_scale = 0.0003280816,
                        cens_coef = c(-0.05, 0.15, 0.01, 0.02),
                        cens_shape = 1.5, cens_scale = 0.00004){
  
  ## generate the covariates
  age <- rnorm(n, age_mean, age_sd)
  volume <- rbinom(n, 1, volume_prop)
  meta <- ifelse(runif(n) < meta_prop0, 0, rpois(n, meta_pois))
  ALP <- rlnorm(n, ALP_mean, ALP_sd)
  primary <- rlnorm(n, primary_mean, primary_sd)
  PSA <- rlnorm(n, psa_mean, psa_sd)
  
  ## generate treatment assignment probabilities via propensity score
  logit_ps <- prop_coef[1] + prop_coef[2]*age + prop_coef[3]*volume + prop_coef[4]*meta + 
    prop_coef[5]*log(ALP) + prop_coef[6]*primary + prop_coef[7]*volume*primary
  A <- rbinom(n, 1, plogis(logit_ps))
  
  ## generate survival and censoring times with Weibull baseline hazard, adjusting for covariates
  Ttime <- rweibull(n, shape = fail_shape, 
                    scale = fail_scale^(-1/fail_shape)*exp(-1*(fail_coef[1]*A + fail_coef[2]*volume +
                                                                 fail_coef[3]*primary + fail_coef[4]*log(PSA) +
                                                                 fail_coef[5]*meta)/fail_shape))
  
  Ctime <- rweibull(n, shape = cens_shape, 
                    scale = cens_scale^(-1/cens_shape)*exp(-1*(cens_coef[1]*A + cens_coef[2]*volume +
                                                                 cens_coef[3]*primary +
                                                                 cens_coef[4]*meta)/cens_shape))
  
  ## return outcomes and event indicators that account for censoring
  ftime <- pmin(Ttime, Ctime)
  ftype <- as.numeric(Ttime <= Ctime)
  
  ## return data frame
  return(data.frame(A = A, age = age, volume = volume, meta = meta, ALP = ALP, 
                    logALP = log(ALP), primary = primary, PSA = PSA, logPSA = log(PSA),
                    ftime, ftype))
}

## this function generates data with a lognormal baseline hazard and implements censoring;
## the inputs are described below
## n: sample size (both groups combined)
## age_mean and age_sd: parameters to generate age covariate
## volume_prop: proportion of high-volume patients (binary)
## meta_prop0: proportion of zero-inflation Possion model for # of bone metastates
## meta_pois: Possion rate for non-zero component of ZIP model
## ALP_mean and ALP_sd: mean and sd for lognormal distribution of Alkaline Phosphatase
## primary_mean and primary_sd: mean and sd for lognormal distribution of # months since primary cancer diagnosed
## psa_mean and psa_sd: mean and sd for lognormal distribution of prostate-specific antigen
## prop_coef: coefficients for propensity score (logistic regression)
## fail_coef: regression coefficients for survival time model
## fail_shape and fail_scale: parameters for lognormal baseline hazard of survival model
## cens_coef: regression coefficients for censoring time model
## cens_shape and cens_scale: parameters for lognormal baseline hazard of censoring model
gen_lnorm <- function(n, age_mean = 69, age_sd = 8, volume_prop = 0.55, meta_prop0 = 0.12,
                      meta_pois = 5.5, ALP_mean = 4.5, ALP_sd = sqrt(2*(log(190) - 4.5)),
                      primary_mean = 3.05, primary_sd = sqrt(2*(log(25) - 3.05)),
                      psa_mean = 2.65, psa_sd = sqrt(2*(log(15) - 2.65)),
                      prop_coef = c(2.6, -0.06, 1, 0.24, 0.34, -0.07, 0.01),
                      fail_coef = c(0.1735, -0.2, -0.01, -0.125, -0.03),
                      fail_shape = 5.211, fail_scale = 0.5,
                      cens_coef = c(0.035, -0.15, -0.01, -0.02),
                      cens_shape = 5.265, cens_scale = 0.4){
  
  ## generate the covariates
  age <- rnorm(n, age_mean, age_sd)
  volume <- rbinom(n, 1, volume_prop)
  meta <- ifelse(runif(n) < meta_prop0, 0, rpois(n, meta_pois))
  ALP <- rlnorm(n, ALP_mean, ALP_sd)
  primary <- rlnorm(n, primary_mean, primary_sd)
  PSA <- rlnorm(n, psa_mean, psa_sd)
  
  ## generate treatment assignment probabilities via propensity score
  logit_ps <- prop_coef[1] + prop_coef[2]*age + prop_coef[3]*volume + prop_coef[4]*meta + 
    prop_coef[5]*log(ALP) + prop_coef[6]*primary + prop_coef[7]*volume*primary
  A <- rbinom(n, 1, plogis(logit_ps))
  
  ## generate survival and censoring times with lognormal baseline hazard, adjusting for covariates
  Ttime <- exp(rnorm(n, fail_shape + fail_coef[1]*A + fail_coef[2]*volume +
                       fail_coef[3]*primary + fail_coef[4]*log(PSA) +
                       fail_coef[5]*meta, fail_scale))
  
  Ctime <- exp(rnorm(n, cens_shape + cens_coef[1]*A + cens_coef[2]*volume +
                       cens_coef[3]*primary +  cens_coef[4]*meta, cens_scale))
  
  ## return outcomes and event indicators that account for censoring
  ftime <- pmin(Ttime, Ctime)
  ftype <- as.numeric(Ttime <= Ctime) # 1=event, 0=censored
  
  ## return data frame
  return(data.frame(A = A, age = age, volume = volume, meta = meta, ALP = ALP, 
                    logALP = log(ALP), primary = primary, PSA = PSA, logPSA = log(PSA),
                    ftime, ftype))
}

## this function generates data with a piecewise constant baseline hazard and implements censoring;
## the inputs are described below
## n: sample size (both groups combined)
## age_mean and age_sd: parameters to generate age covariate
## volume_prop: proportion of high-volume patients (binary)
## meta_prop0: proportion of zero-inflation Possion model for # of bone metastates
## meta_pois: Possion rate for non-zero component of ZIP model
## ALP_mean and ALP_sd: mean and sd for lognormal distribution of Alkaline Phosphatase
## primary_mean and primary_sd: mean and sd for lognormal distribution of # months since primary cancer diagnosed
## psa_mean and psa_sd: mean and sd for lognormal distribution of prostate-specific antigen
## prop_coef: coefficients for propensity score (logistic regression)
## fail_coef: regression coefficients for survival time model
## fail_hazards: piecewise constant baseline hazards of survival model
## cens_coef: regression coefficients for censoring time model
## cens_hazards: piecewise constant baseline hazards of censoring model
## cutpoints: cutpoints for the baseline piecewise function
gen_piece <- function(n, age_mean = 69, age_sd = 8, volume_prop = 0.55, meta_prop0 = 0.12,
                      meta_pois = 5.5, ALP_mean = 4.5, ALP_sd = sqrt(2*(log(190) - 4.5)),
                      primary_mean = 3.05, primary_sd = sqrt(2*(log(25) - 3.05)),
                      psa_mean = 2.65, psa_sd = sqrt(2*(log(15) - 2.65)),
                      prop_coef = c(2.6, -0.06, 1, 0.24, 0.34, -0.07, 0.01),
                      fail_coef = c(-0.47, 0.2, 0.01, 0.125, 0.03),
                      fail_hazards = 0.845*c(0.001, 0.002, 0.003, 0.004, 0.005, 0.006),
                      cens_coef = c(-0.06, 0.15, 0.01, 0.02),
                      cens_hazards = 0.131*c(0.001, 0.002, 0.003, 0.004, 0.005, 0.006),
                      cutpoints = c(12, 24, 36, 48, 60, 600)) {
  
  ## Compute interval widths
  intervals <- c(0, cutpoints)
  widths <- diff(intervals)
  
  ## generate the covariates
  age <- rnorm(n, age_mean, age_sd)
  volume <- rbinom(n, 1, volume_prop)
  meta <- ifelse(runif(n) < meta_prop0, 0, rpois(n, meta_pois))
  ALP <- rlnorm(n, ALP_mean, ALP_sd)
  primary <- rlnorm(n, primary_mean, primary_sd)
  PSA <- rlnorm(n, psa_mean, psa_sd)
  
  ## generate treatment assignment probabilities via propensity score
  logit_ps <- prop_coef[1] + prop_coef[2]*age + prop_coef[3]*volume + prop_coef[4]*meta + 
    prop_coef[5]*log(ALP) + prop_coef[6]*primary + prop_coef[7]*volume*primary
  A <- rbinom(n, 1, plogis(logit_ps))
  
  ## Initialize results and generate pseudorandom sequence for CDF inversion
  t <- rep(0, n)
  event <- rep(0, n)
  u <- runif(n)
  
  ## get the incremental hazards for each segment and participant
  incremental_hazards <- NULL
  for (k in 1:length(fail_hazards)){
    incremental_hazards <- rbind(incremental_hazards, 
                                 widths[k]*fail_hazards[k]*exp(fail_coef[1]*A + fail_coef[2]*volume +
                                                                 fail_coef[3]*primary + fail_coef[4]*log(PSA) +
                                                                 fail_coef[5]*meta))
  }
  
  ## use incremental hazards to get cumulative hazards
  cumulative_hazards <- apply(incremental_hazards, 2, cumsum)
  
  ## get survival times via CDF inversion for each participant
  which.done <- NULL
  ## iterate over the piecewise segments
  for (k in 1:length(fail_hazards)){
    
    ## keep track of which observations have already failed and which have not
    which.above <- which(u > exp(-cumulative_hazards[k,]))
    which.current <- subset(which.above, !(which.above %in% which.done))
    which.done <- c(which.done, which.current)
    
    ## for the observations that fail in this interval, get their survival times
    if (length(which.current) > 0){
      t[which.current] <- intervals[k] - (log(u[which.current]) + cumulative_hazards[k, which.current]
                                          - incremental_hazards[k, which.current])/(fail_hazards[k]*exp(fail_coef[1]*A + fail_coef[2]*volume +
                                                                                                          fail_coef[3]*primary + fail_coef[4]*log(PSA) +
                                                                                                          fail_coef[5]*meta)[which.current])
      
      ## update the event indicator to denote failure
      event[which.current] <- 1
    }
  }
  
  ## adminstrative censoring for observations that have not failed at the end of the last segment
  which.admin <- which(event == 0)
  if (length(which.admin) > 0){
    t[which.admin] <- max(intervals)
  }
  
  ## record the survival times
  Ttime <- t
  
  ## repeat this process for the censoring times
  u <- runif(n)
  
  incremental_hazards <- NULL
  for (k in 1:length(cens_hazards)){
    incremental_hazards <- rbind(incremental_hazards, 
                                 widths[k]*cens_hazards[k]*exp(cens_coef[1]*A + cens_coef[2]*volume +
                                                                 cens_coef[3]*primary + cens_coef[4]*meta))
  }
  
  cumulative_hazards <- apply(incremental_hazards, 2, cumsum)
  
  which.done <- NULL
  for (k in 1:length(cens_hazards)){
    which.above <- which(u > exp(-cumulative_hazards[k,]))
    which.current <- subset(which.above, !(which.above %in% which.done))
    which.done <- c(which.done, which.current)
    
    if (length(which.current) > 0){
      t[which.current] <- intervals[k] - (log(u[which.current]) + cumulative_hazards[k, which.current]
                                          - incremental_hazards[k, which.current])/(cens_hazards[k]*exp(cens_coef[1]*A + cens_coef[2]*volume +
                                                                                                          cens_coef[3]*primary + cens_coef[4]*meta)[which.current])
      event[which.current] <- 1
    }
  }
  
  which.admin <- which(event == 0)
  if (length(which.admin) > 0){
    t[which.admin] <- max(intervals)
  }
  
  ## record the censoring times
  Ctime <- t
  
  ## record the censored outcomes and event indicators
  ftime <- pmin(Ttime, Ctime)
  ftype <- as.numeric(Ttime <= Ctime)
  
  ## account for administrative censoring
  ftype <- ifelse(ftime == tail(cutpoints, 1), 0, ftype)
  
  ## return the data frame
  return(data.frame(data.frame(A = A, age = age, volume = volume, meta = meta, ALP = ALP, 
                               logALP = log(ALP), primary = primary, PSA = PSA, logPSA = log(PSA),
                               ftime, ftype)))
}

## this is a wrapper function that computes a doubly robust
## estimate of the difference in overall survival at time t0
## the inputs are as follows:
## dat: a data frame (generated by a function like gen_weibull())
## t0: the time at which to compute the survival difference (t0 = 60 months by default)
## prop_fun: a function to estimate the propensity score
## surv_fun: a function to estimate survival times
## cens_fun: a function to estimate censoring times
dr_master <- function(dat, t0 = 60, prop_fun, surv_fun, cens_fun){
  
  ## extract columns from data set
  A <- dat$A
  ftime <- dat$ftime
  ftype <- dat$ftype
  
  ## get estimated treatment assignment probabilities from prop_fun()
  ps_hat <- prop_fun(dat)
  
  ## get the counterfactual survival times for each participant from surv_fun() 
  surv_temp <- surv_fun(dat, t0 = t0)
  
  ## extract the survival times under each treatment
  S_pred1 <- surv_temp[[1]]
  S_pred0 <- surv_temp[[2]]
  
  ## estimate the censoring probabilities for each participant from cens_fun()
  G_pred <- cens_fun(dat, t0 = t0)
  
  ## get the weights; incorporate treatment assignment and censoring
  w_treat <- ifelse(A == 1, 1 / ps_hat, 1 / (1 - ps_hat))
  w_censor <- 1 / G_pred
  
  ## only observations that have not been censored before t0
  ## contribute to the estimate
  Y <- 1*(ftime > t0)
  I_obs <- (ftime >= t0) | (ftime <= t0 & ftype==1)
  
  ## get the weights for treatments 0 and 1
  w1 <- ifelse(A == 1, 1, 0)*as.numeric(I_obs)*w_treat*w_censor
  w0 <- ifelse(A == 0, 1, 0)*as.numeric(I_obs)*w_treat*w_censor 
  
  ## estimate survival probability for each group
  S1_DR <- mean(S_pred1 + w1*((ftime >= t0) - S_pred1))
  S0_DR <- mean(S_pred0 + w0*((ftime >= t0) - S_pred0))
  
  ## get the risk difference
  RD_DR <- S1_DR - S0_DR
  return(RD_DR)
}

## this is a wrapper function that implements the nonparametric bootstrap
## to quantify uncertainty; the inputs are as follows
## dat: a data frame (generated by a function like gen_weibull())
## B: the number of bootstrap replicates
## prop_fun: a function to estimate the propensity score
## surv_fun: a function to estimate survival times
## cens_fun: a function to estimate censoring times
## t0: the time at which to compute the survival difference (t0 = 60 months by default)
boot_master <- function(dat, B, prop_fun, surv_fun, cens_fun, t0 = 60){
  
  n_row <- nrow(dat)
  boot_dr <- NULL
  for (b in 1:B) {
    ## resample the data
    idx <- sample(1:n_row, n_row, replace = TRUE)
    ## implement dr_master() with the resampled data
    boot_dr[b] <- dr_master(dat[idx, ], t0 = t0, prop_fun = prop_fun,
                            surv_fun = surv_fun, cens_fun = cens_fun)
  }
  
  ## approximate the standard deviation using the mean absolute deviation
  boot_se_dr <- mad(boot_dr)
  return(boot_se_dr)
  
}

## this function estimates the propensity score
## using the correct model; the input is as follows:
## dat: a data frame (generated by a function like gen_weibull())
prop_fun_correct <- function(dat){
  
  ## estimate probability of being assigned to new treatment
  ps_fit <- glm(A ~ age + volume + meta + logALP + primary + volume*primary, 
                family = binomial, data = dat)
  ps_hat <- predict(ps_fit, type = "response")
  
  ## truncate ps weights for stability
  ps_hat <- pmax(pmin(ps_hat, 1 - 0.01), 0.01)
}

## this function estimates the survival time
## using the correct model; the inputs are as follows:
## dat: a data frame (generated by a function like gen_weibull())
## t0: the time at which to compute the survival difference (t0 = 60 months by default)
surv_fun_correct <- function(dat, t0){
  
  ## fit survival model
  fail_fit <- coxph(Surv(ftime, ftype) ~ A + volume + meta + logPSA + primary, 
                    data = dat)
  
  ## baseline cumulative hazard for failure
  bh_fail <- basehaz(fail_fit, centered = FALSE)
  
  ## interpolate cumulative baseline hazard at t0
  H0_t0_fail <- approx(bh_fail$time, bh_fail$hazard, xout = t0, rule = 2)$y
  
  ## create copies of the data set for each treatment
  dat0 <- dat
  dat0$A <- 0
  
  dat1 <- dat
  dat1$A <- 1
  
  ## get linear predictor under A=1 and A=0
  lp1_fail <- predict(fail_fit, newdata = dat1, type = "lp")
  lp0_fail <- predict(fail_fit, newdata = dat0, type = "lp")
  
  ## get subject-specific survival predictions at t0 under A=1 and A=0
  S_pred1 <- exp(- H0_t0_fail * exp(lp1_fail))
  S_pred0 <- exp(- H0_t0_fail * exp(lp0_fail))
  
  ## return predictions for each treatment group as a list
  return(list(S_pred1, S_pred0))
}

## this function estimates the censoring time
## using the correct model; the inputs are as follows:
## dat: a data frame (generated by a function like gen_weibull())
## t0: the time at which to compute the survival difference (t0 = 60 months by default)
cens_fun_correct <- function(dat, t0){
  ## fit the censoring model on all data
  C_ind <- 1 - dat$ftype
  censor_fit <- coxph(Surv(ftime, C_ind) ~ A + volume + meta + primary,
                      data = dat)
  
  ## get the baseline cumulative hazard for censoring
  bh_censor <- basehaz(censor_fit, centered = FALSE)
  
  ## interpolate the baseline hazard at the time t0 and the observed event time
  H0_t0_censor <- approx(bh_censor$time, bh_censor$hazard, xout = t0, rule = 2)$y
  H0_ftime_censor <- approx(bh_censor$time, bh_censor$hazard, xout = dat$ftime, rule = 2)$y
  
  ## linear predictor for censoring (observed treatment)
  lp_censor <- predict(censor_fit, newdata = dat, type = "lp")
  
  ## probability of being uncensored at t0
  G_pred_t0 <- exp(- H0_t0_censor * exp(lp_censor))
  
  ## probability of being uncensored at ftime
  G_pred_ftime <- exp(- H0_ftime_censor * exp(lp_censor))
  
  ## get the censoring probability at min(t0, ftime)
  G_pred <- ifelse(dat$ftime <= t0, G_pred_ftime, G_pred_t0)
  
  ## truncate G_pred for stability and return results
  G_pred <- pmax(G_pred, 0.01)
  return(G_pred)
}

## consider the six scenarios
ns <- seq(500,1200,50)

## Psi_1: iterate over sample sizes
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('survival'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(m*(j-1) + k)
                       
                       ## generate data for sample size n
                       dat.temp <- gen_weibull(n)
                       
                       ## get DR estimate
                       dr.temp <- dr_master(dat.temp, t0 = 60, prop_fun = prop_fun_correct,
                                            surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       ## estimate the standard error
                       se.temp <- boot_master(dat.temp, B = 100, prop_fun = prop_fun_correct,
                                              surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       ## compute the p-value
                       test.temp <- dr.temp/se.temp
                       p.temp <- 2*pnorm(abs(test.temp), lower.tail = FALSE)
                       
                       ## return p-value and estimate (the latter was used to check the DR property)
                       c(dr.temp, p.temp)
                       
                       
                     }
  
  ## return results as a .csv file
  write.csv(sim.res, paste0("dr_correct_n_", n, ".csv"), row.names = FALSE)
}

## repeat this process for Psi_2
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('survival'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(150000 + m*(j-1) + k)
                       
                       ## the last interaction coefficient is 0 for the true propensity score model
                       dat.temp <- gen_weibull(n,
                                               prop_coef = c(1.25, -0.06, 1, 0.24, 0.34, 0, 0))
                       dr.temp <- dr_master(dat.temp, t0 = 60, prop_fun = prop_fun_correct,
                                            surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       se.temp <- boot_master(dat.temp, B = 100, prop_fun = prop_fun_correct,
                                              surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       test.temp <- dr.temp/se.temp
                       p.temp <- 2*pnorm(abs(test.temp), lower.tail = FALSE)
                       c(dr.temp, p.temp)
                       
                       
                     }
  write.csv(sim.res, paste0("dr_extra_prop_n_", n, ".csv"), row.names = FALSE)
}

## repeat process for Psi_3
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('survival'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(300000 + m*(j-1) + k)
                       
                       ## the coefficient for logPSA is equal to 0 for the true survival model
                       dat.temp <- gen_weibull(n,fail_coef = c(-0.475, 0.2, 0.01, 0, 0.03),
                                               fail_scale = 0.000458)
                       dr.temp <- dr_master(dat.temp, t0 = 60, prop_fun = prop_fun_correct,
                                            surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       se.temp <- boot_master(dat.temp, B = 100, prop_fun = prop_fun_correct,
                                              surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       test.temp <- dr.temp/se.temp
                       p.temp <- 2*pnorm(abs(test.temp), lower.tail = FALSE)
                       c(dr.temp, p.temp)
                       
                       
                     }
  write.csv(sim.res, paste0("dr_extra_out_n_", n, ".csv"), row.names = FALSE)
}

## repeat process for Psi_4
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('survival'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(450000 + m*(j-1) + k)
                       
                       ## change the distribution parameters for the age covariate
                       dat.temp <- gen_weibull(n, age_mean = 59, age_sd = 6)
                       dr.temp <- dr_master(dat.temp, t0 = 60, prop_fun = prop_fun_correct,
                                            surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       se.temp <- boot_master(dat.temp, B = 100, prop_fun = prop_fun_correct,
                                              surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       test.temp <- dr.temp/se.temp
                       p.temp <- 2*pnorm(abs(test.temp), lower.tail = FALSE)
                       c(dr.temp, p.temp)
                       
                       
                     }
  write.csv(sim.res, paste0("dr_age_n_", n, ".csv"), row.names = FALSE)
}

## repeat process for Psi_5
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('survival'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(600000 + m*(j-1) + k)
                       
                       ## generate data with a lognormal baseline hazard this time
                       dat.temp <- gen_lnorm(n)
                       dr.temp <- dr_master(dat.temp, t0 = 60, prop_fun = prop_fun_correct,
                                            surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       se.temp <- boot_master(dat.temp, B = 100, prop_fun = prop_fun_correct,
                                              surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       test.temp <- dr.temp/se.temp
                       p.temp <- 2*pnorm(abs(test.temp), lower.tail = FALSE)
                       c(dr.temp, p.temp)
                       
                       
                     }
  write.csv(sim.res, paste0("dr_lnorm_n_", n, ".csv"), row.names = FALSE)
}

## repeat process for Psi_6
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('survival'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(750000 + m*(j-1) + k)
                       
                       ## generate data with a piecewise constant baseline hazard this time
                       dat.temp <- gen_piece(n)
                       dr.temp <- dr_master(dat.temp, t0 = 60, prop_fun = prop_fun_correct,
                                            surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       se.temp <- boot_master(dat.temp, B = 100, prop_fun = prop_fun_correct,
                                              surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       test.temp <- dr.temp/se.temp
                       p.temp <- 2*pnorm(abs(test.temp), lower.tail = FALSE)
                       c(dr.temp, p.temp)
                       
                       
                     }
  write.csv(sim.res, paste0("dr_piece_n_", n, ".csv"), row.names = FALSE)
}

## get the plots for Figure 1

## this function takes the sampling distribution estimates of p-values
## and constructs linear approximations to estimate power for a point null hypotheses;
## the inputs are described as follows:
## m1 is the first sampling distribution estimate
## m2 is the second sampling distribution estimate
## n0 is first sample size
## n1 is second sample size
## lb is lower bound for linear approximations
## ub is upper bound for linear approximations
## by is the increment for linear approximations
## gam is the threshold (i.e., alpha for type I error rate)
pwr_point_null <- function(m1, m2, n0, n1, lb, ub, by, gam){
  
  ## get logits for the p-values in one tail of 
  ls <- logit(m1/2)
  li <- logit(m2/2)
  
  ## adjust any infinite logits
  ls <- ifelse(ls == -Inf, min(subset(ls, is.finite(ls))) - 1, ls)
  ls <- ifelse(ls == Inf, max(subset(ls, is.finite(ls))) + 1, ls)
  
  ## adjust any infinite logits
  li <- ifelse(li == -Inf, min(subset(li, is.finite(li))) - 1, li)
  li <- ifelse(li == Inf, max(subset(li, is.finite(li))) + 1, li)
  
  slopes <- NULL
  ints <- NULL
  
  ## construct the slopes using the order statistics of the sampling
  ## distribution estimates
  ls_s <- ls[order(ls)]
  li_s <- li[order(li)]
    
  slopes <- (li_s - ls_s)/(n1-n0)
  ints <- ls_s - slopes*n0
  
  ## create matrix to calculate power
  samps <- seq(lb,ub,by)
  res.vec <- NULL
  for (i in 1:length(samps)){
    ## check which probabilities are less than the threshold 
    stop.temp <- 2*expit(ints + slopes*samps[i]) <= gam
    res.vec[i] <- mean(stop.temp)
  }
  
  ## return matrix with samples sizes and estimated power
  return(cbind(samps, res.vec))
  
}

## get power curve using naive simulation
alpha <- 0.1
pwr <- NULL
## read in results from each saved .csv file and calculate power empirically
for (j in 1:length(n_vals)){
  n <- n_vals[j]
  temp <- read.csv(paste0("dr_correct_n_", n, ".csv"))$V2
  pwr[j] <- mean(temp <= alpha)
}
## save results
pwr1 <- pwr

## get power curve using Algorithm 2
temp1.df <- pwr_point_null(read.csv(paste0("dr_correct_n_600.csv"))$V2,
                           read.csv(paste0("dr_correct_n_1100.csv"))$V2,
                           600, 1100, 500, 1200,1, 0.1)

## combine different power curve estimates into one data frame
df1 <- data.frame(n = c(seq(500, 1200, 50), seq(500, 1200, 1)),
                  power = c(pwr1, temp1.df[,2]),
                  curve = c(rep("C_Simulation", length(seq(500, 1200, 50))),
                            rep("A_Algorithm 2", length(seq(500, 1200, 1)))))

## load in colour palette
cbb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## create subplot
plot1 <- ggplot(df1, aes(x=n)) + theme_bw() +
  geom_line(aes(y = power, color=as.factor(curve), linetype = as.factor(curve)), 
            alpha = 0.9, size = 1) +
  labs(title=bquote(psi[1]*': Weibull Baseline Hazard')) +
  labs(x= bquote(italic(n)), y= bquote('Power')) +
  theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                  margin = margin(t = 0, 0, 5, 0))) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  theme(legend.position="none") +
  scale_color_manual(name = " ", 
                     labels = c("Algorithm 2  ", "Simulation"),
                     values = c("steelblue1", "firebrick")) +
  scale_linetype_manual(name = " ", 
                        labels = c("Algorithm 2  ", "Simulation"),
                        values = c(1, 5)) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.key.size = unit(1, "cm")) +
  ylim(0.55,1) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  geom_hline(yintercept=0.8, lty = 2)

## repeat for the second model
pwr <- NULL
for (j in 1:length(n_vals)){
  n <- n_vals[j]
  temp <- read.csv(paste0("dr_extra_prop_n_", n, ".csv"))$V2
  pwr[j] <- mean(temp <= alpha)
}
pwr2 <- pwr

## get power curve using Algorithm 2
temp2.df <- pwr_point_null(read.csv(paste0("dr_extra_prop_n_600.csv"))$V2,
                           read.csv(paste0("dr_extra_prop_n_1100.csv"))$V2,
                           600, 1100, 500, 1200,1, 0.1)

## combine different power curve estimates into one data frame
df2 <- data.frame(n = c(seq(500, 1200, 50), seq(500, 1200, 1)),
                  power = c(pwr2, temp2.df[,2]),
                  curve = c(rep("C_Simulation", length(seq(500, 1200, 50))),
                            rep("A_Algorithm 2", length(seq(500, 1200, 1)))))

## create subplot
plot2 <- ggplot(df2, aes(x=n)) + theme_bw() +
  geom_line(aes(y = power, color=as.factor(curve), linetype = as.factor(curve)), 
            alpha = 0.9, size = 1) +
  labs(title=bquote(psi[2]*': Simpler Propensity Score')) +
  labs(x= bquote(italic(n)), y= bquote('Power')) +
  theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                  margin = margin(t = 0, 0, 5, 0))) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  theme(legend.position="none") +
  scale_color_manual(name = " ", 
                     labels = c("Algorithm 2  ", "Simulation"),
                     values = c("steelblue1", "firebrick")) +
  scale_linetype_manual(name = " ", 
                        labels = c("Algorithm 2  ", "Simulation"),
                        values = c(1, 5)) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.key.size = unit(1, "cm")) +
  ylim(0.55,1) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  geom_hline(yintercept=0.8, lty = 2)

## repeat for the third model
pwr <- NULL
for (j in 1:length(n_vals)){
  n <- n_vals[j]
  temp <- read.csv(paste0("dr_extra_out_n_", n, ".csv"))$V2
  pwr[j] <- mean(temp <= alpha)
}
pwr3 <- pwr

## get power curve using Algorithm 2
temp3.df <- pwr_point_null(read.csv(paste0("dr_extra_out_n_600.csv"))$V2,
                           read.csv(paste0("dr_extra_out_n_1100.csv"))$V2,
                           600, 1100, 500, 1200,1, 0.1)

## combine different power curve estimates into one data frame
df3 <- data.frame(n = c(seq(500, 1200, 50), seq(500, 1200, 1)),
                  power = c(pwr3, temp3.df[,2]),
                  curve = c(rep("C_Simulation", length(seq(500, 1200, 50))),
                            rep("A_Algorithm 2", length(seq(500, 1200, 1)))))

## create subplot
plot3 <- ggplot(df3, aes(x=n)) + theme_bw() +
  geom_line(aes(y = power, color=as.factor(curve), linetype = as.factor(curve)), 
            alpha = 0.9, size = 1) +
  labs(title=bquote(psi[3]*': Simpler Survival Model')) +
  labs(x= bquote(italic(n)), y= bquote('Power')) +
  theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                  margin = margin(t = 0, 0, 5, 0))) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  theme(legend.position="none") +
  scale_color_manual(name = " ", 
                     labels = c("Algorithm 2  ", "Simulation"),
                     values = c("steelblue1", "firebrick")) +
  scale_linetype_manual(name = " ", 
                        labels = c("Algorithm 2  ", "Simulation"),
                        values = c(1, 5)) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.key.size = unit(1, "cm")) +
  ylim(0.55,1) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  geom_hline(yintercept=0.8, lty = 2)

## repeat for the fourth model
pwr <- NULL
for (j in 1:length(n_vals)){
  n <- n_vals[j]
  temp <- read.csv(paste0("dr_age_n_", n, ".csv"))$V2
  pwr[j] <- mean(temp <= alpha)
}
pwr4 <- pwr

## get power curve using Algorithm 2
temp4.df <- pwr_point_null(read.csv(paste0("dr_age_n_600.csv"))$V2,
                           read.csv(paste0("dr_age_n_1050.csv"))$V2,
                           600, 1050, 500, 1200,1, 0.1)

## combine different power curve estimates into one data frame
df4 <- data.frame(n = c(seq(500, 1200, 50), seq(500, 1200, 1)),
                  power = c(pwr4, temp4.df[,2]),
                  curve = c(rep("C_Simulation", length(seq(500, 1200, 50))),
                            rep("A_Algorithm 2", length(seq(500, 1200, 1)))))

## create subplot
plot4 <- ggplot(df4, aes(x=n)) + theme_bw() +
  geom_line(aes(y = power, color=as.factor(curve), linetype = as.factor(curve)), 
            alpha = 0.9, size = 1) +
  labs(title=bquote(psi[4]*': New Covariate Distributions')) +
  labs(x= bquote(italic(n)), y= bquote('Power')) +
  theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                  margin = margin(t = 0, 0, 5, 0))) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  theme(legend.position="none") +
  scale_color_manual(name = " ", 
                     labels = c("Algorithm 2  ", "Simulation"),
                     values = c("steelblue1", "firebrick")) +
  scale_linetype_manual(name = " ", 
                        labels = c("Algorithm 2  ", "Simulation"),
                        values = c(1, 5)) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.key.size = unit(1, "cm")) +
  ylim(0.55,1) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  geom_hline(yintercept=0.8, lty = 2)

## repeat for the fifth model
pwr <- NULL
for (j in 1:length(n_vals)){
  n <- n_vals[j]
  temp <- read.csv(paste0("dr_lnorm_n_", n, ".csv"))$V2
  pwr[j] <- mean(temp <= alpha)
}
pwr5 <- pwr

## get power curve using Algorithm 2
temp5.df <- pwr_point_null(read.csv(paste0("dr_lnorm_n_600.csv"))$V2,
                           read.csv(paste0("dr_lnorm_n_1100.csv"))$V2,
                           600, 1100, 500, 1200,1, 0.1)

## combine different power curve estimates into one data frame
df5 <- data.frame(n = c(seq(500, 1200, 50), seq(500, 1200, 1)),
                  power = c(pwr5, temp5.df[,2]),
                  curve = c(rep("C_Simulation", length(seq(500, 1200, 50))),
                            rep("A_Algorithm 2", length(seq(500, 1200, 1)))))

## create subplot
plot5 <- ggplot(df5, aes(x=n)) + theme_bw() +
  geom_line(aes(y = power, color=as.factor(curve), linetype = as.factor(curve)), 
            alpha = 0.9, size = 1) +
  labs(title=bquote(psi[5]*': Lognormal Baseline Hazard')) +
  labs(x= bquote(italic(n)), y= bquote('Power')) +
  theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                  margin = margin(t = 0, 0, 5, 0))) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  theme(legend.position="none") +
  scale_color_manual(name = " ", 
                     labels = c("Algorithm 2  ", "Simulation"),
                     values = c("steelblue1", "firebrick")) +
  scale_linetype_manual(name = " ", 
                        labels = c("Algorithm 2  ", "Simulation"),
                        values = c(1, 5)) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.key.size = unit(1, "cm")) +
  ylim(0.55,1) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  geom_hline(yintercept=0.8, lty = 2)

## repeat for the sixth model
pwr <- NULL
for (j in 1:length(n_vals)){
  n <- n_vals[j]
  temp <- read.csv(paste0("dr_piece_n_", n, ".csv"))$V2
  pwr[j] <- mean(temp <= alpha)
}
pwr6 <- pwr

## get power curve using Algorithm 2
temp6.df <- pwr_point_null(read.csv(paste0("dr_piece_n_600.csv"))$V2,
                           read.csv(paste0("dr_piece_n_1100.csv"))$V2,
                           600, 1100, 500, 1200,1, 0.1)

## combine different power curve estimates into one data frame
df6 <- data.frame(n = c(seq(500, 1200, 50), seq(500, 1200, 1)),
                  power = c(pwr6, temp6.df[,2]),
                  curve = c(rep("C_Simulation", length(seq(500, 1200, 50))),
                            rep("A_Algorithm 2", length(seq(500, 1200, 1)))))

## create subplot
plot6 <- ggplot(df6, aes(x=n)) + theme_bw() +
  geom_line(aes(y = power, color=as.factor(curve), linetype = as.factor(curve)), 
            alpha = 0.9, size = 1) +
  labs(title=bquote(psi[6]*': Piecewise Baseline Hazard')) +
  labs(x= bquote(italic(n)), y= bquote('Power')) +
  theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                  margin = margin(t = 0, 0, 5, 0))) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  theme(legend.position="none") +
  scale_color_manual(name = " ", 
                     labels = c("Algorithm 2  ", "Simulation"),
                     values = c("steelblue1", "firebrick")) +
  scale_linetype_manual(name = " ", 
                        labels = c("Algorithm 2  ", "Simulation"),
                        values = c(1, 5)) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.key.size = unit(1, "cm")) +
  ylim(0.55,1) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  geom_hline(yintercept=0.8, lty = 2)

## arrange for the final plot
figp.row1 <- plot_grid(plot1 + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")), 
                       plot2 + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")),
                       rel_widths = c(1, 1))
figp.row2 <- plot_grid(plot3 + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")), 
                       plot4 + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")),
                       rel_widths = c(1, 1))
figp.row3 <- plot_grid(plot5 + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")), 
                       plot6 + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")),
                       rel_widths = c(1, 1))
figp <- plot_grid(figp.row1, figp.row2, figp.row3, nrow = 3)

## get a common legend for the larger figure
plot1.legend <- ggplot(df1, aes(x=n)) + theme_bw() +
  geom_line(aes(y = power, color=as.factor(curve), linetype = as.factor(curve)), 
            alpha = 0.9, size = 1) +
  labs(title=bquote(psi[1]*': Weibull Baseline Hazard')) +
  labs(x= bquote(italic(n)), y= bquote('Power')) +
  theme(plot.title = element_text(size=18,face="bold", hjust =  0.5,
                                  margin = margin(t = 0, 0, 5, 0))) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  theme(legend.position="bottom") +
  scale_color_manual(name = " ", 
                     labels = c("Algorithm 2  ", "Simulation"),
                     values = c("steelblue1", "firebrick")) +
  scale_linetype_manual(name = " ", 
                        labels = c("Algorithm 2  ", "Simulation"),
                        values = c(1, 5)) +
  ylim(0.55,1) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  geom_hline(yintercept=0.8, lty = 2) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.key.size = unit(1.5, "cm"))

## add legend to bottom of the plot
fig_final <- plot_grid(figp, ggpubr::get_legend(plot1.legend), ncol = 1, rel_heights = c(2.1, .1))

# output as .pdf file for the article
pdf(file = "Fig_DR.pdf",   # The directory you want to save the file in
    width = 10.5, # The width of the plot in inches
    height = 9.75) # The height of the plot in inches

fig_final

dev.off()

## now implement the confirmatory simulations for Appendix C

## get recommended sample size (robust)
ns <- 951

## confirmation for Psi_1 (power)
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('survival'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(900000 + m*(j-1) + k)
                       
                       dat.temp <- gen_weibull(n)
                       dr.temp <- dr_master(dat.temp, t0 = 60, prop_fun = prop_fun_correct,
                                            surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       se.temp <- boot_master(dat.temp, B = 100, prop_fun = prop_fun_correct,
                                              surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       test.temp <- dr.temp/se.temp
                       p.temp <- 2*pnorm(abs(test.temp), lower.tail = FALSE)
                       c(dr.temp, p.temp)
                       
                       
                     }
  write.csv(sim.res, paste0("conf_dr_correct_n_", n, ".csv"), row.names = FALSE)
}

## confirmation for Psi_2 (power)
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('survival'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(910000 + m*(j-1) + k)
                       
                       dat.temp <- gen_weibull(n,
                                               prop_coef = c(1.25, -0.06, 1, 0.24, 0.34, 0, 0))
                       dr.temp <- dr_master(dat.temp, t0 = 60, prop_fun = prop_fun_correct,
                                            surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       se.temp <- boot_master(dat.temp, B = 100, prop_fun = prop_fun_correct,
                                              surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       test.temp <- dr.temp/se.temp
                       p.temp <- 2*pnorm(abs(test.temp), lower.tail = FALSE)
                       c(dr.temp, p.temp)
                       
                       
                     }
  write.csv(sim.res, paste0("conf_dr_extra_prop_n_", n, ".csv"), row.names = FALSE)
}

## confirmation for Psi_3 (power)
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('survival'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(920000 + m*(j-1) + k)
                       
                       dat.temp <- gen_weibull(n,fail_coef = c(-0.475, 0.2, 0.01, 0, 0.03),
                                               fail_scale = 0.000458)
                       dr.temp <- dr_master(dat.temp, t0 = 60, prop_fun = prop_fun_correct,
                                            surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       se.temp <- boot_master(dat.temp, B = 100, prop_fun = prop_fun_correct,
                                              surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       test.temp <- dr.temp/se.temp
                       p.temp <- 2*pnorm(abs(test.temp), lower.tail = FALSE)
                       c(dr.temp, p.temp)
                       
                       
                     }
  write.csv(sim.res, paste0("conf_dr_extra_out_n_", n, ".csv"), row.names = FALSE)
}

## confirmation for Psi_4 (power)
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('survival'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(930000 + m*(j-1) + k)
                       
                       dat.temp <- gen_weibull(n, age_mean = 59, age_sd = 6)
                       dr.temp <- dr_master(dat.temp, t0 = 60, prop_fun = prop_fun_correct,
                                            surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       se.temp <- boot_master(dat.temp, B = 100, prop_fun = prop_fun_correct,
                                              surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       test.temp <- dr.temp/se.temp
                       p.temp <- 2*pnorm(abs(test.temp), lower.tail = FALSE)
                       c(dr.temp, p.temp)
                       
                       
                     }
  write.csv(sim.res, paste0("conf_dr_age_n_", n, ".csv"), row.names = FALSE)
}

## confirmation for Psi_5 (power)
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('survival'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(940000 + m*(j-1) + k)
                       
                       dat.temp <- gen_lnorm(n)
                       dr.temp <- dr_master(dat.temp, t0 = 60, prop_fun = prop_fun_correct,
                                            surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       se.temp <- boot_master(dat.temp, B = 100, prop_fun = prop_fun_correct,
                                              surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       test.temp <- dr.temp/se.temp
                       p.temp <- 2*pnorm(abs(test.temp), lower.tail = FALSE)
                       c(dr.temp, p.temp)
                       
                       
                     }
  write.csv(sim.res, paste0("conf_dr_lnorm_n_", n, ".csv"), row.names = FALSE)
}

## confirmation for Psi_6 (power)
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('survival'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(950000 + m*(j-1) + k)
                       
                       dat.temp <- gen_piece(n)
                       dr.temp <- dr_master(dat.temp, t0 = 60, prop_fun = prop_fun_correct,
                                            surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       se.temp <- boot_master(dat.temp, B = 100, prop_fun = prop_fun_correct,
                                              surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       test.temp <- dr.temp/se.temp
                       p.temp <- 2*pnorm(abs(test.temp), lower.tail = FALSE)
                       c(dr.temp, p.temp)
                       
                       
                     }
  write.csv(sim.res, paste0("conf_dr_piece_n_", n, ".csv"), row.names = FALSE)
}

## repeat simulations under H0 to check type I error

## confirmation for Psi_1 (type I error)
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('survival'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(960000 + m*(j-1) + k)
                       
                       dat.temp <- gen_weibull(n,
                                               fail_coef = c(0, 0.2, 0.01, 0.125, 0.03),
                                               cens_coef = c(0, 0.15, 0.01, 0.02))
                       dr.temp <- dr_master(dat.temp, t0 = 60, prop_fun = prop_fun_correct,
                                            surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       se.temp <- boot_master(dat.temp, B = 100, prop_fun = prop_fun_correct,
                                              surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       test.temp <- dr.temp/se.temp
                       p.temp <- 2*pnorm(abs(test.temp), lower.tail = FALSE)
                       c(dr.temp, p.temp)
                       
                       
                     }
  write.csv(sim.res, paste0("conf_H0_dr_correct_n_", n, ".csv"), row.names = FALSE)
}

## confirmation for Psi_2 (type I error)
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('survival'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(970000 + m*(j-1) + k)
                       
                       dat.temp <- gen_weibull(n,
                                               prop_coef = c(1.25, -0.06, 1, 0.24, 0.34, 0, 0),
                                               fail_coef = c(0, 0.2, 0.01, 0.125, 0.03),
                                               cens_coef = c(0, 0.15, 0.01, 0.02))
                       dr.temp <- dr_master(dat.temp, t0 = 60, prop_fun = prop_fun_correct,
                                            surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       se.temp <- boot_master(dat.temp, B = 100, prop_fun = prop_fun_correct,
                                              surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       test.temp <- dr.temp/se.temp
                       p.temp <- 2*pnorm(abs(test.temp), lower.tail = FALSE)
                       c(dr.temp, p.temp)
                       
                       
                     }
  write.csv(sim.res, paste0("conf_H0_dr_extra_prop_n_", n, ".csv"), row.names = FALSE)
}

## confirmation for Psi_3 (type I error)
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('survival'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(980000 + m*(j-1) + k)
                       
                       dat.temp <- gen_weibull(n,fail_coef = c(0, 0.2, 0.01, 0, 0.03),
                                               cens_coef = c(0, 0.15, 0.01, 0.02),
                                               fail_scale = 0.000458)
                       dr.temp <- dr_master(dat.temp, t0 = 60, prop_fun = prop_fun_correct,
                                            surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       se.temp <- boot_master(dat.temp, B = 100, prop_fun = prop_fun_correct,
                                              surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       test.temp <- dr.temp/se.temp
                       p.temp <- 2*pnorm(abs(test.temp), lower.tail = FALSE)
                       c(dr.temp, p.temp)
                       
                       
                     }
  write.csv(sim.res, paste0("conf_H0_dr_extra_out_n_", n, ".csv"), row.names = FALSE)
}

## confirmation for Psi_4 (type I error)
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('survival'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(990000 + m*(j-1) + k)
                       
                       dat.temp <- gen_weibull(n, age_mean = 59, age_sd = 6,
                                               fail_coef = c(0, 0.2, 0.01, 0.125, 0.03),
                                               cens_coef = c(0, 0.15, 0.01, 0.02))
                       dr.temp <- dr_master(dat.temp, t0 = 60, prop_fun = prop_fun_correct,
                                            surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       se.temp <- boot_master(dat.temp, B = 100, prop_fun = prop_fun_correct,
                                              surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       test.temp <- dr.temp/se.temp
                       p.temp <- 2*pnorm(abs(test.temp), lower.tail = FALSE)
                       c(dr.temp, p.temp)
                       
                       
                     }
  write.csv(sim.res, paste0("conf_H0_dr_age_n_", n, ".csv"), row.names = FALSE)
}

## confirmation for Psi_5 (type I error)
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('survival'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(1000000 + m*(j-1) + k)
                       
                       dat.temp <- gen_lnorm(n,
                                             fail_coef = c(0, -0.2, -0.01, -0.125, -0.03),
                                             cens_coef = c(0, -0.15, -0.01, -0.02))
                       dr.temp <- dr_master(dat.temp, t0 = 60, prop_fun = prop_fun_correct,
                                            surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       se.temp <- boot_master(dat.temp, B = 100, prop_fun = prop_fun_correct,
                                              surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       test.temp <- dr.temp/se.temp
                       p.temp <- 2*pnorm(abs(test.temp), lower.tail = FALSE)
                       c(dr.temp, p.temp)
                       
                       
                     }
  write.csv(sim.res, paste0("conf_H0_dr_lnorm_n_", n, ".csv"), row.names = FALSE)
}

## confirmation for Psi_6 (type I error)
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('survival'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(1010000 + m*(j-1) + k)
                       
                       dat.temp <- gen_piece(n,
                                             fail_coef = c(0, 0.2, 0.01, 0.125, 0.03),
                                             cens_coef = c(0, 0.15, 0.01, 0.02))
                       dr.temp <- dr_master(dat.temp, t0 = 60, prop_fun = prop_fun_correct,
                                            surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       se.temp <- boot_master(dat.temp, B = 100, prop_fun = prop_fun_correct,
                                              surv_fun = surv_fun_correct, cens_fun = cens_fun_correct)
                       
                       test.temp <- dr.temp/se.temp
                       p.temp <- 2*pnorm(abs(test.temp), lower.tail = FALSE)
                       c(dr.temp, p.temp)
                       
                       
                     }
  write.csv(sim.res, paste0("conf_H0_dr_piece_n_", n, ".csv"), row.names = FALSE)
}

## get the confirmatory estimates for power and the type I error rate

## get the file names for power
files <- c("conf_dr_correct_n_1001.csv", "conf_dr_extra_prop_n_1001.csv",
          "conf_dr_extra_out_n_1001.csv", "conf_dr_age_n_1001.csv",
          "conf_dr_lnorm_n_1001.csv", "conf_dr_piece_n_1001.csv")

## get a vector of power estimates for each scenario
pwr_conf <- NULL
for (j in 1:length(files)){
  sim.res <- read.csv(files[j])$V2
  pwr_conf[j] <- mean(sim.res <= 0.1)
}

## get file names for H0
filesH0 <- c("conf_H0_dr_correct_n_1001.csv", "conf_H0_dr_extra_prop_n_1001.csv",
             "conf_H0_dr_extra_out_n_1001.csv", "conf_H0_dr_age_n_1001.csv",
             "conf_H0_dr_lnorm_n_1001.csv", "conf_H0_dr_piece_n_1001.csv")

## get a vector of type I error estimates for each scenario
t1E_conf <- NULL
for (j in 1:length(files)){
  sim.res <- read.csv(filesH0[j])$V2
  t1E_conf[j] <- mean(sim.res <= 0.1)
}