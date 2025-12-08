## this code file is used to reproduce the results from Section 5.3

## load necessary packages
require(foreach)
require(doParallel)
require(doSNOW)
require(geepack)
require(MASS)
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

## this function generates data with a Poisson mean model, overdispersion, and
## potential latent correlation; the inputs are described below
## n: sample size (both groups combined)
## c: the number of observations per cluster
## beta: the regression coefficients for the Poisson mean model
## trt_prob: the probability of being assigned to the new treatment
## nb_size: overdispersion parameter for count data
## corr_struct: latent correlation structure (independent, exchangeable, ar1, or unstructured)
## corr_param: parameter(s) for correlation structure
## diagnostics: boolean to check induced marginal dependence structure and overdispersion
gen_cluster <- function(n, c = 5, beta = c(1.42, 0, -0.1, 0), trt_prob = 0.5, nb_size = 15, 
                        corr_struct = "independent", corr_param = 0.2,
                        diagnostics = FALSE){
  
  ## generate covariates
  id <- rep(1:n, each = c)
  x1 <- rep(rbinom(n, 1, trt_prob), each = c)
  x2 <- rep(c(0, 1, 1, 1, 1), n)
  x1x2 <- x1*x2
  t <- rep(c(8, 2, 2, 2, 2), n)
  
  ## create the marginal mean for each observation
  eta <- beta[1] + beta[2]*x1 + beta[3]*x2 + beta[4]*x1*x2
  mu_vec <- as.numeric(t*exp(eta))
  
  ## helper function to create correlation matrix
  make_R <- function(c, struct, param, distances = NULL) {
    if (struct == "independent"){
      return(diag(c))
    }
    else if (struct == "exchangeable") {
      rho <- param
      R <- matrix(rho, nrow = c, ncol = c)
      diag(R) <- 1
      return(R)
    }
    else if (struct == "ar1") {
      rho <- param
      idx <- 1:c
      R <- outer(idx, idx, function(i,j) rho^(abs(i-j)))
      return(R)
    }
    else if (struct == "unstructured") {
      R <- diag(c)
      R[lower.tri(R)] <- param
      R[upper.tri(R)] <- t(R)[upper.tri(R)]
      return(R)
    }
  }
  
  ## generate vectors for response and copula models
  Y <- numeric(n*c)
  Lambda <- numeric(n*c)
  Uvals <- numeric(n*c)
  R <- make_R(c, corr_struct, corr_param)
  for (i in 1:n) {
    idx <- ((i-1)*c + 1):(i*c)
    mu_block <- mu_vec[idx]
    ## draw multivariate normal with mean 0 and covariance R
    Z <- MASS::mvrnorm(1, mu = rep(0, c), Sigma = R)
    U <- pnorm(Z)  ## convert to Gaussian copula
    ## transform to NB marginal via the quantile function
    Ys <- qnbinom(U, size = nb_size, mu = mu_block)
    Y[idx] <- Ys
    Lambda[idx] <- mu_block
    Uvals[idx] <- U
  }
  
  ## create data frame
  df <- data.frame(id = id, x1 = x1, x2 = x2, x1x2 = x1x2, t = t, mu = mu_vec, y = Y)
  
  ## if diagnostics is TRUE, return information to calibrate parameters
  if (diagnostics){
    ## empirical overdispersion
    emp_mean <- mean(df$y)
    emp_var  <- var(df$y)
    emp_phi  <- emp_var/emp_mean
    
    ## construct wide cluster-by-position matrix
    Ymat <- matrix(df$y, nrow = n, ncol = c, byrow = TRUE)
    
    ## variable-by-variable correlation (across clusters)
    cm <- cor(Ymat)
    emp_within_rho <- mean(cm[lower.tri(cm)]) ## average correlation across pairwise components
    
    ## return empirical dispersion summaries
    return(list(emp_mean = emp_mean, emp_var = emp_var, emp_phi = emp_phi,
                cm = cm, emp_within_rho = emp_within_rho))
  }
  
  ## return data frame
  if (diagnostics == FALSE){
    return(df)
  }
  
}

## consider the scenarios
ns <- seq(30, 90, 5)

## Psi_1: iterate over sample sizes
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('geepack', 'MASS'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(m*(j-1) + k)
                       
                       ## generate data
                       temp.dat <- gen_cluster(n)
                       
                       fit <- geese(y~x1*x2+offset(log(t)), id=id, corstr="independence",
                                    family="poisson", data=temp.dat)
                       
                       ## p-values for each one-sided hypothesis
                       eff.temp <- as.numeric(fit$beta[4])
                       se.temp <- sqrt(fit$vbeta[4,4])
                       
                       p1.temp <- pnorm((eff.temp + log(4/3))/se.temp, lower.tail = FALSE)
                       p2.temp <- pnorm((eff.temp - log(4/3))/se.temp)
                       
                       ## return results
                       as.numeric(c(p1.temp, p2.temp))
                       
                     }
  
  ## save results as a .csv file
  write.csv(sim.res, paste0("gee_ind_n_", n, ".csv"), row.names = FALSE)
}

## Psi_2: iterate over sample sizes
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('geepack', 'MASS'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(m*(j-1) + k + 130000)
                       
                       ## generate data
                       temp.dat <- gen_cluster(n, corr_struct = "exchangeable",
                                               corr_param = 0.25)
                       
                       fit <- geese(y~x1*x2+offset(log(t)), id=id, corstr="independence",
                                    family="poisson", data=temp.dat)
                       
                       ## p-values for each one-sided hypothesis
                       eff.temp <- as.numeric(fit$beta[4])
                       se.temp <- sqrt(fit$vbeta[4,4])
                       
                       p1.temp <- pnorm((eff.temp + log(4/3))/se.temp, lower.tail = FALSE)
                       p2.temp <- pnorm((eff.temp - log(4/3))/se.temp)
                       
                       ## return results
                       as.numeric(c(p1.temp, p2.temp))
                       
                     }
  
  ## save results as a .csv file
  write.csv(sim.res, paste0("gee_exc_n_", n, ".csv"), row.names = FALSE)
}

## Psi_3: iterate over sample sizes
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('geepack', 'MASS'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(m*(j-1) + k + 260000)
                       
                       ## generate data
                       temp.dat <- gen_cluster(n, corr_struct = "ar1",
                                               corr_param = 0.5)
                       
                       fit <- geese(y~x1*x2+offset(log(t)), id=id, corstr="independence",
                                    family="poisson", data=temp.dat)
                       
                       ## p-values for each one-sided hypothesis
                       eff.temp <- as.numeric(fit$beta[4])
                       se.temp <- sqrt(fit$vbeta[4,4])
                       
                       p1.temp <- pnorm((eff.temp + log(4/3))/se.temp, lower.tail = FALSE)
                       p2.temp <- pnorm((eff.temp - log(4/3))/se.temp)
                       
                       ## return results
                       as.numeric(c(p1.temp, p2.temp))
                       
                     }
  
  ## save results as a .csv file
  write.csv(sim.res, paste0("gee_ar1_n_", n, ".csv"), row.names = FALSE)
}

## Psi_4: iterate over sample sizes
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('geepack', 'MASS'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(m*(j-1) + k + 390000)
                       
                       ## generate data
                       temp.dat <- gen_cluster(n, corr_struct = "unstructured",
                                               corr_param = c(0.05, 0.05, 0.05, 0.05, 0.3, 0.2, 0.1, 0.3,0.2,0.1))
                       
                       fit <- geese(y~x1*x2+offset(log(t)), id=id, corstr="independence",
                                    family="poisson", data=temp.dat)
                       
                       ## p-values for each one-sided hypothesis
                       eff.temp <- as.numeric(fit$beta[4])
                       se.temp <- sqrt(fit$vbeta[4,4])
                       
                       p1.temp <- pnorm((eff.temp + log(4/3))/se.temp, lower.tail = FALSE)
                       p2.temp <- pnorm((eff.temp - log(4/3))/se.temp)
                       
                       ## return results
                       as.numeric(c(p1.temp, p2.temp))
                       
                     }
  
  ## save results as a .csv file
  write.csv(sim.res, paste0("gee_uns_n_", n, ".csv"), row.names = FALSE)
}

## get the plots for Figure 3

## this function takes the sampling distribution estimates of p-values
## and constructs linear approximations to estimate power for equivalence tests;
## the inputs are described as follows:
## m1 is the first sampling distribution estimate (one column for each hypothesis)
## m2 is the second sampling distribution estimate (one column for each hypothesis)
## n0 is first sample size
## n1 is second sample size
## lb is lower bound for linear approximations
## ub is upper bound for linear approximations
## by is the increment for linear approximations
## gam is the threshold (i.e., alpha for type I error rate)

## get the power curve for the equivalence test
pwr_equiv <- function(m1, m2, n0, n1, lb, ub, by, gam){
  
  ## get logits for the p-values for both hypotheses 
  ls <- logit(m1)
  li <- logit(m2)
  
  ## adjust any infinite logits
  for (i in 1:2){
    ls[,i] <- ifelse(ls[,i] == -Inf, min(subset(ls[,i], is.finite(ls[,i]))) - 1, ls[,i])
    ls[,i] <- ifelse(ls[,i] == Inf, max(subset(ls[,i], is.finite(ls[,i]))) + 1, ls[,i])
    
    li[,i] <- ifelse(li[,i] == -Inf, min(subset(li[,i], is.finite(li[,i]))) - 1, li[,i])
    li[,i] <- ifelse(li[,i] == Inf, max(subset(li[,i], is.finite(li[,i]))) + 1, li[,i])
  }
  
  ## get indexes to combine individual logits later
  ls <- cbind(ls, seq(1, nrow(ls), 1))
  li <- cbind(li, seq(1, nrow(li), 1))
  
  slopes <- NULL
  ints <- NULL
  
  ## construct the slopes separately for each hypothesis
  for (j in 1:2){
    ls_s <- ls[order(ls[,j]),j]
    li_s <- li[order(li[,j]),j]
    
    l_slope <- (li_s - ls_s)/(n1-n0)
    l_int <- ls_s - l_slope*n0
    
    ## reorder according to smaller sample size
    l_slope[ls[order(ls[,j]),ncol(ls)]] <- l_slope 
    l_int[ls[order(ls[,j]),ncol(ls)]] <- l_int 
    
    slopes <- cbind(slopes, l_slope)
    ints <- cbind(ints, l_int)
  }
  
  ## create matrix to calculate power
  samps <- seq(lb,ub,by)
  res.vec <- NULL
  for (i in 1:length(samps)){
    ## check which probabilities are less than the threshold;
    ## we need to model the logit of the p-value for both hypotheses
    logit1.temp <- ints[,1] + slopes[,1]*samps[i] 
    logit2.temp <- ints[,2] + slopes[,2]*samps[i]
    stop.temp <- pmax(logit1.temp, logit2.temp) <= logit(gam) 
    res.vec[i] <- mean(stop.temp)
  }
  
  ## return matrix with samples sizes and estimated power
  return(cbind(samps, res.vec))
  
}

n_vals <- seq(30,90,5)
## get power curve using naive simulation
alpha <- 0.05
pwr <- NULL
## read in results from each saved .csv file and calculate power empirically
for (j in 1:length(n_vals)){
  n <- n_vals[j]
  temp <- read.csv(paste0("gee_ind_n_", n, ".csv"))
  pwr[j] <- mean(pmax(temp$V1, temp$V2) <= alpha)
}
## save results
pwr1 <- pwr

## get power curve using Algorithm 2
temp1.df <- pwr_equiv(read.csv(paste0("gee_ind_n_40.csv")), 
                      read.csv(paste0("gee_ind_n_80.csv")), 
                      40, 80, 30, 90, 1, 0.05)

## combine different power curve estimates into one data frame
df1 <- data.frame(n = c(seq(30, 90, 5), seq(30, 90, 1)),
                  power = c(pwr1, temp1.df[,2]),
                  curve = c(rep("C_Simulation", length(seq(30, 90, 5))),
                            rep("A_Algorithm 2", length(seq(30, 90, 1)))))

## load in colour palette
cbb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## create subplot
plot1 <- ggplot(df1, aes(x=n)) + theme_bw() +
  geom_line(aes(y = power, color=as.factor(curve), linetype = as.factor(curve)), 
            alpha = 0.9, size = 1) +
  labs(title=bquote(psi[1]*': Independence')) +
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
  ylim(0.2,1) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  geom_hline(yintercept=0.8, lty = 2)

## repeat for the second scenario
pwr <- NULL
## read in results from each saved .csv file and calculate power empirically
for (j in 1:length(n_vals)){
  n <- n_vals[j]
  temp <- read.csv(paste0("gee_exc_n_", n, ".csv"))
  pwr[j] <- mean(pmax(temp$V1, temp$V2) <= alpha)
}
## save results
pwr2 <- pwr

## get power curve using Algorithm 2
temp2.df <- pwr_equiv(read.csv(paste0("gee_exc_n_40.csv")), 
                      read.csv(paste0("gee_exc_n_80.csv")), 
                      40, 80, 30, 90, 1, 0.05)

## combine different power curve estimates into one data frame
df2 <- data.frame(n = c(seq(30, 90, 5), seq(30, 90, 1)),
                  power = c(pwr2, temp2.df[,2]),
                  curve = c(rep("C_Simulation", length(seq(30, 90, 5))),
                            rep("A_Algorithm 2", length(seq(30, 90, 1)))))

## create subplot
plot2 <- ggplot(df2, aes(x=n)) + theme_bw() +
  geom_line(aes(y = power, color=as.factor(curve), linetype = as.factor(curve)), 
            alpha = 0.9, size = 1) +
  labs(title=bquote(psi[2]*': Exchangeable')) +
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
  ylim(0.2,1) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  geom_hline(yintercept=0.8, lty = 2)

## repeat for the third scenario
pwr <- NULL
## read in results from each saved .csv file and calculate power empirically
for (j in 1:length(n_vals)){
  n <- n_vals[j]
  temp <- read.csv(paste0("gee_ar1_n_", n, ".csv"))
  pwr[j] <- mean(pmax(temp$V1, temp$V2) <= alpha)
}
## save results
pwr3 <- pwr

## get power curve using Algorithm 2
temp3.df <- pwr_equiv(read.csv(paste0("gee_ar1_n_40.csv")), 
                      read.csv(paste0("gee_ar1_n_80.csv")), 
                      40, 80, 30, 90, 1, 0.05)

## combine different power curve estimates into one data frame
df3 <- data.frame(n = c(seq(30, 90, 5), seq(30, 90, 1)),
                  power = c(pwr3, temp3.df[,2]),
                  curve = c(rep("C_Simulation", length(seq(30, 90, 5))),
                            rep("A_Algorithm 2", length(seq(30, 90, 1)))))

## create subplot
plot3 <- ggplot(df3, aes(x=n)) + theme_bw() +
  geom_line(aes(y = power, color=as.factor(curve), linetype = as.factor(curve)), 
            alpha = 0.9, size = 1) +
  labs(title=bquote(psi[3]*': AR(1)')) +
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
  ylim(0.2,1) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  geom_hline(yintercept=0.8, lty = 2)

## repeat for the fourth scenario
pwr <- NULL
## read in results from each saved .csv file and calculate power empirically
for (j in 1:length(n_vals)){
  n <- n_vals[j]
  temp <- read.csv(paste0("gee_uns_n_", n, ".csv"))
  pwr[j] <- mean(pmax(temp$V1, temp$V2) <= alpha)
}
## save results
pwr4 <- pwr

## get power curve using Algorithm 2
temp4.df <- pwr_equiv(read.csv(paste0("gee_uns_n_40.csv")), 
                      read.csv(paste0("gee_uns_n_80.csv")), 
                      40, 80, 30, 90, 1, 0.05)

## combine different power curve estimates into one data frame
df4 <- data.frame(n = c(seq(30, 90, 5), seq(30, 90, 1)),
                  power = c(pwr4, temp4.df[,2]),
                  curve = c(rep("C_Simulation", length(seq(30, 90, 5))),
                            rep("A_Algorithm 2", length(seq(30, 90, 1)))))

## create subplot
plot4 <- ggplot(df4, aes(x=n)) + theme_bw() +
  geom_line(aes(y = power, color=as.factor(curve), linetype = as.factor(curve)), 
            alpha = 0.9, size = 1) +
  labs(title=bquote(psi[4]*': Unstructured')) +
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
  ylim(0.2,1) + 
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
figp <- plot_grid(figp.row1, figp.row2, nrow = 2)

## get a common legend for the larger figure
plot1.legend <- ggplot(df1, aes(x=n)) + theme_bw() +
  geom_line(aes(y = power, color=as.factor(curve), linetype = as.factor(curve)), 
            alpha = 0.9, size = 1) +
  labs(title=bquote(psi[1]*': Independence')) +
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
  ylim(0.2,1) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  geom_hline(yintercept=0.8, lty = 2) +
  theme(legend.text=element_text(size=18)) +
  theme(legend.key.size = unit(1.5, "cm"))

## add legend to bottom of the plot
fig_final <- plot_grid(figp, ggpubr::get_legend(plot1.legend), ncol = 1, rel_heights = c(1.5, .1))

# output as .pdf file for the article
pdf(file = "Fig_GEE.pdf",   # The directory you want to save the file in
    width = 10.5, # The width of the plot in inches
    height = 6.75) # The height of the plot in inches

fig_final

dev.off()

## now implement the confirmatory simulations for Appendix E

## get recommended sample size (robust)
ns <- 69

## confirmation for Psi_1 (power)
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('geepack', 'MASS'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(m*(j-1) + k + 520000)
                       
                       ## generate data
                       temp.dat <- gen_cluster(n)
                       
                       fit <- geese(y~x1*x2+offset(log(t)), id=id, corstr="independence",
                                    family="poisson", data=temp.dat)
                       
                       ## p-values for each one-sided hypothesis
                       eff.temp <- as.numeric(fit$beta[4])
                       se.temp <- sqrt(fit$vbeta[4,4])
                       
                       p1.temp <- pnorm((eff.temp + log(4/3))/se.temp, lower.tail = FALSE)
                       p2.temp <- pnorm((eff.temp - log(4/3))/se.temp)
                       
                       ## return results
                       as.numeric(c(p1.temp, p2.temp))
                       
                     }
  
  ## save results as a .csv file
  write.csv(sim.res, paste0("conf_gee_ind_n_", n, ".csv"), row.names = FALSE)
}

## confirmation for Psi_2 (power)
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('geepack', 'MASS'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(m*(j-1) + k + 530000)
                       
                       ## generate data
                       temp.dat <- gen_cluster(n, corr_struct = "exchangeable",
                                               corr_param = 0.25)
                       
                       fit <- geese(y~x1*x2+offset(log(t)), id=id, corstr="independence",
                                    family="poisson", data=temp.dat)
                       
                       ## p-values for each one-sided hypothesis
                       eff.temp <- as.numeric(fit$beta[4])
                       se.temp <- sqrt(fit$vbeta[4,4])
                       
                       p1.temp <- pnorm((eff.temp + log(4/3))/se.temp, lower.tail = FALSE)
                       p2.temp <- pnorm((eff.temp - log(4/3))/se.temp)
                       
                       ## return results
                       as.numeric(c(p1.temp, p2.temp))
                       
                     }
  
  ## save results as a .csv file
  write.csv(sim.res, paste0("conf_gee_exc_n_", n, ".csv"), row.names = FALSE)
}

## confirmation for Psi_3 (power)
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('geepack', 'MASS'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(m*(j-1) + k + 540000)
                       
                       ## generate data
                       temp.dat <- gen_cluster(n, corr_struct = "ar1",
                                               corr_param = 0.5)
                       
                       fit <- geese(y~x1*x2+offset(log(t)), id=id, corstr="independence",
                                    family="poisson", data=temp.dat)
                       
                       ## p-values for each one-sided hypothesis
                       eff.temp <- as.numeric(fit$beta[4])
                       se.temp <- sqrt(fit$vbeta[4,4])
                       
                       p1.temp <- pnorm((eff.temp + log(4/3))/se.temp, lower.tail = FALSE)
                       p2.temp <- pnorm((eff.temp - log(4/3))/se.temp)
                       
                       ## return results
                       as.numeric(c(p1.temp, p2.temp))
                       
                     }
  
  ## save results as a .csv file
  write.csv(sim.res, paste0("conf_gee_ar1_n_", n, ".csv"), row.names = FALSE)
}

## confirmation for Psi_4 (power)
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('geepack', 'MASS'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(m*(j-1) + k + 550000)
                       
                       ## generate data
                       temp.dat <- gen_cluster(n, corr_struct = "unstructured",
                                               corr_param = c(0.05, 0.05, 0.05, 0.05, 0.3, 0.2, 0.1, 0.3,0.2,0.1))
                       
                       fit <- geese(y~x1*x2+offset(log(t)), id=id, corstr="independence",
                                    family="poisson", data=temp.dat)
                       
                       ## p-values for each one-sided hypothesis
                       eff.temp <- as.numeric(fit$beta[4])
                       se.temp <- sqrt(fit$vbeta[4,4])
                       
                       p1.temp <- pnorm((eff.temp + log(4/3))/se.temp, lower.tail = FALSE)
                       p2.temp <- pnorm((eff.temp - log(4/3))/se.temp)
                       
                       ## return results
                       as.numeric(c(p1.temp, p2.temp))
                       
                     }
  
  ## save results as a .csv file
  write.csv(sim.res, paste0("conf_gee_uns_n_", n, ".csv"), row.names = FALSE)
}

## confirmation for Psi_1 (type I error)
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('geepack', 'MASS'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(m*(j-1) + k + 560000)
                       
                       ## generate data
                       temp.dat <- gen_cluster(n,
                                               beta = c(1.42, 0, -0.1, ifelse(k%%2 == 0, log(4/3), -log(4/3))))
                       
                       fit <- geese(y~x1*x2+offset(log(t)), id=id, corstr="independence",
                                    family="poisson", data=temp.dat)
                       
                       ## p-values for each one-sided hypothesis
                       eff.temp <- as.numeric(fit$beta[4])
                       se.temp <- sqrt(fit$vbeta[4,4])
                       
                       p1.temp <- pnorm((eff.temp + log(4/3))/se.temp, lower.tail = FALSE)
                       p2.temp <- pnorm((eff.temp - log(4/3))/se.temp)
                       
                       ## return results
                       as.numeric(c(p1.temp, p2.temp))
                       
                     }
  
  ## save results as a .csv file
  write.csv(sim.res, paste0("conf_H0_gee_ind_n_", n, ".csv"), row.names = FALSE)
}

## confirmation for Psi_2 (type I error)
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('geepack', 'MASS'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(m*(j-1) + k + 570000)
                       
                       ## generate data
                       temp.dat <- gen_cluster(n, corr_struct = "exchangeable",
                                               corr_param = 0.25,
                                               beta = c(1.42, 0, -0.1, ifelse(k%%2 == 0, log(4/3), -log(4/3))))
                       
                       fit <- geese(y~x1*x2+offset(log(t)), id=id, corstr="independence",
                                    family="poisson", data=temp.dat)
                       
                       ## p-values for each one-sided hypothesis
                       eff.temp <- as.numeric(fit$beta[4])
                       se.temp <- sqrt(fit$vbeta[4,4])
                       
                       p1.temp <- pnorm((eff.temp + log(4/3))/se.temp, lower.tail = FALSE)
                       p2.temp <- pnorm((eff.temp - log(4/3))/se.temp)
                       
                       ## return results
                       as.numeric(c(p1.temp, p2.temp))
                       
                     }
  
  ## save results as a .csv file
  write.csv(sim.res, paste0("conf_H0_gee_exc_n_", n, ".csv"), row.names = FALSE)
}

## confirmation for Psi_3 (type I error)
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('geepack', 'MASS'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(m*(j-1) + k + 580000)
                       
                       ## generate data
                       temp.dat <- gen_cluster(n, corr_struct = "ar1",
                                               corr_param = 0.5,
                                               beta = c(1.42, 0, -0.1, ifelse(k%%2 == 0, log(4/3), -log(4/3))))
                       
                       fit <- geese(y~x1*x2+offset(log(t)), id=id, corstr="independence",
                                    family="poisson", data=temp.dat)
                       
                       ## p-values for each one-sided hypothesis
                       eff.temp <- as.numeric(fit$beta[4])
                       se.temp <- sqrt(fit$vbeta[4,4])
                       
                       p1.temp <- pnorm((eff.temp + log(4/3))/se.temp, lower.tail = FALSE)
                       p2.temp <- pnorm((eff.temp - log(4/3))/se.temp)
                       
                       ## return results
                       as.numeric(c(p1.temp, p2.temp))
                       
                     }
  
  ## save results as a .csv file
  write.csv(sim.res, paste0("conf_H0_gee_ar1_n_", n, ".csv"), row.names = FALSE)
}

## confirmation for Psi_4 (type I error)
for (j in 1:length(ns)){
  n <- ns[j]
  
  sim.res <- foreach(k=1:10000, .combine=rbind, .packages=c('geepack', 'MASS'),
                     .options.snow=opts) %dopar% {
                       
                       set.seed(m*(j-1) + k + 590000)
                       
                       ## generate data
                       temp.dat <- gen_cluster(n, corr_struct = "unstructured",
                                               corr_param = c(0.05, 0.05, 0.05, 0.05, 0.3, 0.2, 0.1, 0.3,0.2,0.1),
                                               beta = c(1.42, 0, -0.1, ifelse(k%%2 == 0, log(4/3), -log(4/3))))
                       
                       fit <- geese(y~x1*x2+offset(log(t)), id=id, corstr="independence",
                                    family="poisson", data=temp.dat)
                       
                       ## p-values for each one-sided hypothesis
                       eff.temp <- as.numeric(fit$beta[4])
                       se.temp <- sqrt(fit$vbeta[4,4])
                       
                       p1.temp <- pnorm((eff.temp + log(4/3))/se.temp, lower.tail = FALSE)
                       p2.temp <- pnorm((eff.temp - log(4/3))/se.temp)
                       
                       ## return results
                       as.numeric(c(p1.temp, p2.temp))
                       
                     }
  
  ## save results as a .csv file
  write.csv(sim.res, paste0("conf_H0_gee_uns_n_", n, ".csv"), row.names = FALSE)
}

## get the confirmatory estimates for power and the type I error rate

## get the file names for power
files <- c("conf_gee_ind_n_69.csv", "conf_gee_exc_n_69.csv",
           "conf_gee_ar1_n_69.csv", "conf_gee_uns_n_69.csv")

## get a vector of power estimates for each scenario
pwr_conf <- NULL
for (j in 1:length(files)){
  sim.res <- read.csv(files[j])
  pwr_conf[j] <- mean(pmax(sim.res$V1, sim.res$V2) <= 0.05)
}

## get file names for H0
filesH0 <- c("conf_H0_gee_ind_n_69.csv", "conf_H0_gee_exc_n_69.csv",
             "conf_H0_gee_ar1_n_69.csv", "conf_H0_gee_uns_n_69.csv")

## get a vector of type I error estimates for each scenario
t1E_conf <- NULL
for (j in 1:length(files)){
  sim.res <- read.csv(filesH0[j])
  t1E_conf[j] <- mean(pmax(sim.res$V1, sim.res$V2) <= 0.05)
}

## get the sample size recommendation based on a weighted average
## of data generation processes

## this function is similar to pwr_equiv(), but now the slopes
## and intercepts are returned instead of the power estimates;
## these inputs are a subset of the inputs form pwr_equiv()
pwr_lines <- function(m1, m2, n0, n1){
  
  ## get logits for the p-values for both hypotheses 
  ls <- logit(m1)
  li <- logit(m2)
  
  ## adjust any infinite logits
  for (i in 1:2){
    ls[,i] <- ifelse(ls[,i] == -Inf, min(subset(ls[,i], is.finite(ls[,i]))) - 1, ls[,i])
    ls[,i] <- ifelse(ls[,i] == Inf, max(subset(ls[,i], is.finite(ls[,i]))) + 1, ls[,i])
    
    li[,i] <- ifelse(li[,i] == -Inf, min(subset(li[,i], is.finite(li[,i]))) - 1, li[,i])
    li[,i] <- ifelse(li[,i] == Inf, max(subset(li[,i], is.finite(li[,i]))) + 1, li[,i])
  }
  
  ## get indexes to combine individual logits later
  ls <- cbind(ls, seq(1, nrow(ls), 1))
  li <- cbind(li, seq(1, nrow(li), 1))
  
  slopes <- NULL
  ints <- NULL
  
  ## construct the slopes separately for each hypothesis
  for (j in 1:2){
    ls_s <- ls[order(ls[,j]),j]
    li_s <- li[order(li[,j]),j]
    
    l_slope <- (li_s - ls_s)/(n1-n0)
    l_int <- ls_s - l_slope*n0
    
    ## reorder according to smaller sample size
    l_slope[ls[order(ls[,j]),ncol(ls)]] <- l_slope 
    l_int[ls[order(ls[,j]),ncol(ls)]] <- l_int 
    
    slopes <- cbind(slopes, l_slope)
    ints <- cbind(ints, l_int)
  }
  
  ## return matrix with samples sizes and estimated power
  return(list(ints, slopes))
  
}

## get the slopes and intercepts for all submodels
lines1 <- pwr_lines(read.csv(paste0("gee_ind_n_40.csv")), 
                    read.csv(paste0("gee_ind_n_80.csv")), 
                    40, 80)

lines2 <- pwr_lines(read.csv(paste0("gee_exc_n_40.csv")), 
                    read.csv(paste0("gee_exc_n_80.csv")), 
                    40, 80)

lines3 <- pwr_lines(read.csv(paste0("gee_ar1_n_40.csv")), 
                    read.csv(paste0("gee_ar1_n_80.csv")), 
                    40, 80)

lines4 <- pwr_lines(read.csv(paste0("gee_uns_n_40.csv")), 
                    read.csv(paste0("gee_uns_n_80.csv")), 
                    40, 80)

## this function takes the intercepts (ints) and slopes above
## and computes the power curve based on a weighted average
## of data-generation processes
weighted_pwr <- function(ints, slopes, lb, ub, by, gam){
  
  ## create matrix to calculate power
  samps <- seq(lb,ub,by)
  res.vec <- NULL
  for (i in 1:length(samps)){
    ## check which probabilities are less than the threshold;
    ## we need to model the logit of the p-value for both hypotheses
    logit1.temp <- ints[,1] + slopes[,1]*samps[i] 
    logit2.temp <- ints[,2] + slopes[,2]*samps[i]
    stop.temp <- pmax(logit1.temp, logit2.temp) <= logit(gam) 
    res.vec[i] <- mean(stop.temp)
  }
  
  ## return matrix with samples sizes and estimated power
  return(cbind(samps, res.vec))
  
}

ints.all <- rbind(lines1[[1]], lines2[[1]], lines3[[1]], lines4[[1]])
slopes.all <- rbind(lines1[[2]], lines2[[2]], lines3[[2]], lines4[[2]])

w.pwr <- weighted_pwr(ints.all, slopes.all, 30, 90, 1, 0.05)