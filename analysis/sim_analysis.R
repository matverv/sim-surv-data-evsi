rm(list=ls(all=TRUE))	

#####################################################################################
# load functions
#####################################################################################
source(here::here("Analysis", "sim_ipd_functions.R"))
source(here::here("Analysis", "evsi_functions.R"))
source(here::here("Analysis", "data_gen_functions.R"))

#####################################################################################
# Analysis settings 
#####################################################################################
set.seed(2)   # set the seed

n_sim <- 10e3               # number of simulations in the probabilistic analysis
cycle_year <- 12            # number of model cycles per year
t_h <- 30*cycle_year        # model time horizon
d_r <- 0.035                # annual discount rate for costs and effects

waning <- FALSE             # apply waning of treatment effect (TRUE/FALSE)
wane_start <- 2*cycle_year  # mean point at which waning of treatment effect starts
wane_end <- 4*cycle_year    # mean point at which waning of treatment effect ends (expected hazards are equal between arms from this point)

threshold <- 8e4            # opportunity cost/willingness to pay

#####################################################################################
# Simulate IPD
#####################################################################################

# number of patients in each treatment arm
npat1 <- 100
npat2 <- 100

# observed follow-up time
t1 <- 24

# Weibull shapes and scales for simulating the hypothetical datasets
shape_w1 <- c(exp(0.30), exp(0.15))
shape_w2 <- c(exp(0.35), exp(0.20))
scale_w1 <- c(exp(4.10), exp(3.60))
scale_w2 <- c(exp(3.85), exp(3.30))

os_pfs_times1 <- mapply(gen_times_q_fun, npat1, shape_w1, scale_w1)
os_pfs_times2 <- mapply(gen_times_q_fun, npat2, shape_w2, scale_w2)

ipd1 <- apply(os_pfs_times1, 2, ipd_cens_fun,  tmin = 0, tmax = t1, treat = 1)
ipd2 <- apply(os_pfs_times2, 2, ipd_cens_fun,  tmin = 0, tmax = t1, treat = 2)
names(ipd1) <- names(ipd2) <- c("OS", "PFS")


#####################################################################################
# Fit survival models and compute survival curves
#####################################################################################

# Fit survival models using MLE and obtain summaries
dists <- "weibull"
surv_summ1 <- survfit_fun(dists, ipd1, t_h, c("OS", "PFS"))
surv_summ2 <- survfit_fun(dists, ipd2, t_h, c("OS", "PFS"))

## Weibull model parameters
theta1_os <- MASS::mvrnorm(n_sim,
                           surv_summ1$OS$fit$weibull[,1],
                           surv_summ1$OS$fit$weibull[,2:3])

theta2_os <- MASS::mvrnorm(n_sim,
                           surv_summ2$OS$fit$weibull[,1],
                           surv_summ2$OS$fit$weibull[,2:3])

theta1_pfs <- MASS::mvrnorm(n_sim,
                            surv_summ1$PFS$fit$weibull[,1],
                            surv_summ1$PFS$fit$weibull[,2:3])

theta2_pfs <- MASS::mvrnorm(n_sim,
                            surv_summ2$PFS$fit$weibull[,1],
                            surv_summ2$PFS$fit$weibull[,2:3])

# compute the proportion surviving at each model cycle
S_os1 <- apply(theta1_os, 1, function (x) pweibull(0:t_h, exp(x[1]), exp(x[2]), lower.tail = F))
S_os2 <- apply(theta2_os, 1, function (x) pweibull(0:t_h, exp(x[1]), exp(x[2]), lower.tail = F))
S_pfs1 <- apply(theta1_pfs, 1, function (x) pweibull(0:t_h, exp(x[1]), exp(x[2]), lower.tail = F))
S_pfs2 <- apply(theta2_pfs, 1, function (x) pweibull(0:t_h, exp(x[1]), exp(x[2]), lower.tail = F))


#####################################################################################
# Waning of the treatment effect
#####################################################################################
if(waning == TRUE) {
  
  # specify lognormal distributions for start and end times of the waning period
  wane_start_par <- lnorm_mom(wane_start, wane_start) # start of waning period
  wane_end_par <- lnorm_mom(wane_end-wane_start, wane_end-wane_start) # duration of waning period
  
  # sample waning start times
  wane_start_runif <- runif(n_sim,
                            0,
                            plnorm(t_h, wane_start_par$mu, sdlog = sqrt(wane_start_par$sigma2))
  )
  wane_start_times <- round(qlnorm(wane_start_runif, wane_start_par$mu, sqrt(wane_start_par$sigma2)))
  
  # sample waning end times
  wane_end_runif <- runif(n_sim,
                          0,
                          plnorm(t_h, wane_end_par$mu, sdlog = sqrt(wane_end_par$sigma2))
  )
  wane_end_times <- wane_start_times + round(qlnorm(wane_end_runif, wane_end_par$mu, sqrt(wane_end_par$sigma2)))
  wane_end_times[wane_end_times>t_h] <- t_h
  
  ### compute the hazard
  # Overall survival
  h_os1 <- diff(-log(rowMeans(S_os1))) # hweibull(1:t_h, exp(mean(theta1_os[,1])), exp(mean(theta1_os[,2]))) #
  h_os2 <- diff(-log(rowMeans(S_os2))) # hweibull(1:t_h, exp(mean(theta2_os[,1])), exp(mean(theta2_os[,2]))) #
  
  # Progression-free survival
  h_pfs1 <- diff(-log(rowMeans(S_pfs1))) # hweibull(1:t_h, exp(mean(theta1_pfs[,1])), exp(mean(theta1_pfs[,2]))) 
  h_pfs2 <- diff(-log(rowMeans(S_pfs2))) # hweibull(1:t_h, exp(mean(theta2_pfs[,1])), exp(mean(theta2_pfs[,2]))) 
  
  for (i in 1:n_sim) {
    
    #### overall survival
    # compute additive waning hazards
    h_wane_os1 <- wane_fun(h_os1, h_os2, wane_start_times[i], wane_end_times[i])
    add_h_wane_os1 <- h_wane_os1 - h_os1
    
    # add the "waning" hazards to the unadjusted hazards
    h_temp <- diff(-log(S_os1[,i])) + add_h_wane_os1
    
    # compute adjusted survival curves
    S_os1[,i] <- cumprod(c(1, exp(-h_temp)))
    
    #### progression-free survival
    # compute additive waning hazards
    h_wane_pfs1 <- wane_fun(h_pfs1, h_pfs2, wane_start_times[i], wane_end_times[i])
    add_h_wane_pfs1 <- h_wane_pfs1 - h_pfs1
    
    # add the "waning" hazards to the unadjusted hazards
    h_temp <- diff(-log(S_pfs1[,i])) + add_h_wane_pfs1
    
    # compute adjusted survival curves
    S_pfs1[,i] <- cumprod(c(1, exp(-h_temp)))
    
  }
  
  S_os1[is.nan(S_os1)] <- S_pfs1[is.nan(S_pfs1)] <- 0 # to avoid NaN due to extremely small numbers
  
}    

# adjust PFS curves to ensure PFS doesn't exceed OS
S_pfs1 <- matrix(mapply(min, S_pfs1, S_os1), nrow(S_pfs1), ncol(S_pfs1))
S_pfs2 <- matrix(mapply(min, S_pfs2, S_os2), nrow(S_pfs2), ncol(S_pfs2))

# plot survival curves
plot(x = c(0:(t_h)), y = rowMeans(S_os1), 
     type = "l", col = "2", lty = 1, xlim = c(0,t_h), ylim = c(0,1), 
     xlab="Months", ylab="Survival probability",  
     yaxt="n",  bty = "l", yaxs = "i", xaxs="i", yaxp = c(0,1,4), las = 1,  xaxp = c(0, t_h, 15)) 
axis(side = 2, las = 2, at = c(0,0.25,0.50,0.75,1), labels = c("0","0.25","0.50","0.75","1"))
lines(rowMeans(S_os2), x = 0:t_h, type = "l", col = "4", lty = 1)
lines(x = c(0:t_h), rowMeans(S_pfs1), type = "l", col = "2", lty = 2)
lines(x = c(0:t_h), rowMeans(S_pfs2), type = "l", col = "4", lty = 2)
legend(x="topright", legend = c("OS new treatment", "OS standard care", "PFS new treatment", "PFS standard care"),
       lty = c(1,1,2,2), col = c(2, 4, 2, 4), bty = "n", cex = 1, y.intersp=1)


#####################################################################################
# Utilities and costs
#####################################################################################

# quality of life weights 
u_pfs <- rbeta(n_sim, 80, 20)
u_pps <- rbeta(n_sim, 50, 50)

# drug costs (per cycle)
c_drug1 <- rep(14400, n_sim)/cycle_year #22000

# medical costs (per cycle)
c_med1 <- rgamma(n_sim, 500^2/250^2, 1/(250^2/500)) #rgamma(n_sim, 6000^2/3000^2, 1/(3000^2/6000))/cycle_year
c_med2 <- rgamma(n_sim, 500^2/250^2, 1/(250^2/500)) #rgamma(n_sim, 6000^2/3000^2, 1/(3000^2/6000))/cycle_year


#####################################################################################
# Compute prior net benefit
#####################################################################################

# apply discounting to the survival curves
v_d_r <- 1/(1+d_r)^seq.int(0, t_h/cycle_year, (1/cycle_year)) # vector of monthly discounting factors

S_os1_d <- S_os1*v_d_r
S_os2_d <- S_os2*v_d_r
S_pfs1_d <- S_pfs1*v_d_r
S_pfs2_d <- S_pfs2*v_d_r

# compute life-years
ly1 <- colSums(S_os1_d)/cycle_year
ly2 <- colSums(S_os2_d)/cycle_year

# compute QALYs
qalys_pfs1 <- colSums(mapply("*", as.data.frame(S_pfs1_d), u_pfs))/cycle_year
qalys_pfs2 <- colSums(mapply("*", as.data.frame(S_pfs2_d), u_pfs))/cycle_year

qalys_pps1 <- colSums(mapply("*", as.data.frame(S_os1_d-S_pfs1_d), u_pps))/cycle_year
qalys_pps2 <- colSums(mapply("*", as.data.frame(S_os2_d-S_pfs2_d), u_pps))/cycle_year

qalys_tot1 <- qalys_pfs1 + qalys_pps1
qalys_tot2 <- qalys_pfs2 + qalys_pps2

# compute costs
if(waning == TRUE) {
  costs_drug1 <- colSums(sapply(1:n_sim, function (x) { # drug costs applied only until start of waning period 
    S_pfs1_d[,x] * c(rep(c_drug1[x], wane_start), rep(0, t_h-wane_start+1))
  })
  ) 
} else {
  costs_drug1 <- colSums(mapply("*", as.data.frame(S_pfs1_d), c_drug1))
}

costs_med1 <- colSums(mapply("*", as.data.frame(S_os1_d), c_med1))
costs_med2 <- colSums(mapply("*", as.data.frame(S_os2_d), c_med2))

costs_tot1 <- costs_med1 + costs_drug1
costs_tot2 <- costs_med2 # + costs_drug2

# cost-effectiveness results
mean(costs_tot1) # mean total costs for intervention 1
mean(costs_tot2) # mean total costs for intervention 2
mean(qalys_tot1) # mean total QALYs for intervention 1
mean(qalys_tot2) # mean total QALYs for intervention 2

paste("ICER = ", round(mean(costs_tot1-costs_tot2)/mean(qalys_tot1-qalys_tot2),0), sep = "")


#####################################################################################
# (Partial Expected Value of Perfect Information
#####################################################################################
os_curves <- list(S_os1, S_os2)           
pfs_curves <- list(S_pfs1, S_pfs2)        
theta_os <- list(theta1_os, theta2_os)    
theta_pfs <- list(theta1_pfs, theta2_pfs) 

inb <- threshold * (qalys_tot1-qalys_tot2) - (costs_tot1 - costs_tot2) # incremental net monetary benefit

##### EVPI #####
################
evpi <- mean(pmax(inb,0)) - mean(inb)
evpi

##### partial EVPI #####
########################
evppi_os <- surv_evppi_fun(os_curves)
evppi_os

evppi_pfs <- surv_evppi_fun(pfs_curves)
evppi_pfs

evppi_os_pfs <- surv_evppi_fun(os_curves, pfs_curves)
evppi_os_pfs


#####################################################################################
# Expected Value of Sample Information
#####################################################################################
start_times_os <- list(rep(0,100),        # start times for patients at risk for OS
                       rep(0,100)) 
start_times_pfs <- list(rep(0,100),       # start times for patients at risk for PFS
                        rep(0,100)) 
add_fu <- c(12, 24, 36, 48, 60, 120)      # additional follow-up times


##### Overall survival #####
evsi_os_int <- surv_evsi_fun(
  inb = inb,
  os_curves = os_curves,
  start_times_os = start_times_os,
  add_fu = add_fu,
  method = "interpolation",
  seed = 1
)

evsi_os_disc <- surv_evsi_fun(
  inb = inb,
  os_curves = os_curves,
  start_times_os = start_times_os,
  add_fu = add_fu,
  method = "discrete",
  fast = FALSE,
  seed = 1
)

evsi_os_std <- surv_evsi_fun(
  inb = inb,
  os_curves = os_curves,
  theta_os = theta_os,
  start_times_os = start_times_os,
  add_fu = add_fu,
  method = "standard",
  seed = 1
)

##### Overall survival + progression-free survival #####
evsi_os_pfs_int <- surv_evsi_fun(
  inb = inb,
  os_curves = os_curves,
  pfs_curves = pfs_curves,
  start_times_os = start_times_os,
  start_times_pfs = start_times_pfs,
  add_fu = add_fu,
  method = "interpolation",
  seed = 1
)

evsi_os_pfs_disc <- surv_evsi_fun(
  inb,
  os_curves = os_curves,
  pfs_curves = pfs_curves,
  start_times_os = start_times_os,
  start_times_pfs = start_times_pfs,
  add_fu = add_fu,
  method = "discrete",
  fast = F,
  seed = 1
)

evsi_os_pfs_std <- surv_evsi_fun(
  inb,
  os_curves = os_curves,
  pfs_curves = pfs_curves,
  theta_os = theta_os,
  theta_pfs = theta_pfs,
  start_times_os = start_times_os,
  start_times_pfs = start_times_pfs,
  add_fu = add_fu,
  method = "standard",
  seed = 1
)

# store results
evsi_no_wane <- list("evpi" = evpi, "evppi_os" = evppi_os, "evppi_pfs" = evppi_pfs, "evppi_os_pfs" = evppi_os_pfs, 
                     "evsi_os_int"=evsi_os_int, "evsi_os_disc"=evsi_os_disc, "evsi_os_std"=evsi_os_std,
                     "evsi_os_pfs_int"=evsi_os_pfs_int, "evsi_os_pfs_disc"=evsi_os_pfs_disc, "evsi_os_pfs_std"=evsi_os_pfs_std) 

surv_curves_no_wane <- list("S_os1"=S_os1, "S_os2"=S_os2, "S_pfs1"=S_pfs1, "S_pfs2"=S_pfs2, 
                            "theta_os1"=theta1_os,"theta2_os" =theta2_os, "theta1_pfs" = theta1_pfs, "theta2_pfs"=theta2_pfs)

ce_results_no_wane <- list("inb"=inb,
                           "costs_tot1"=costs_tot1, "costs_tot2"=costs_tot2, 
                           "qalys_tot1"=qalys_tot1, "qalys_tot2"=qalys_tot2, 
                           "ly1"=ly1, "ly2"=ly2)  
  

  

  
  
  
  
  
#####################################################################################
# load functions
#####################################################################################
source(here::here("Analysis", "sim_ipd_functions.R"))
source(here::here("Analysis", "evsi_functions.R"))
source(here::here("Analysis", "data_gen_functions.R"))

#####################################################################################
# Analysis settings 
#####################################################################################
set.seed(2)   # set the seed

n_sim <- 10e3               # number of simulations in the probabilistic analysis
cycle_year <- 12            # number of model cycles per year
t_h <- 30*cycle_year        # model time horizon
d_r <- 0.035                # annual discount rate for costs and effects

waning <- TRUE              # apply waning of treatment effect (TRUE/FALSE)
wane_start <- 2*cycle_year  # mean point at which waning of treatment effect starts
wane_end <- 4*cycle_year    # mean point at which waning of treatment effect ends (expected hazards are equal between arms from this point)

threshold <- 8e4            # opportunity cost/willingness to pay

#####################################################################################
# Simulate IPD
#####################################################################################

# number of patients in each treatment arm
npat1 <- 100
npat2 <- 100

# observed follow-up time
t1 <- 24

# Weibull shapes and scales for simulating the hypothetical datasets
shape_w1 <- c(exp(0.30), exp(0.15))
shape_w2 <- c(exp(0.35), exp(0.20))
scale_w1 <- c(exp(4.10), exp(3.60))
scale_w2 <- c(exp(3.85), exp(3.30))

os_pfs_times1 <- mapply(gen_times_q_fun, npat1, shape_w1, scale_w1)
os_pfs_times2 <- mapply(gen_times_q_fun, npat2, shape_w2, scale_w2)

ipd1 <- apply(os_pfs_times1, 2, ipd_cens_fun,  tmin = 0, tmax = t1, treat = 1)
ipd2 <- apply(os_pfs_times2, 2, ipd_cens_fun,  tmin = 0, tmax = t1, treat = 2)
names(ipd1) <- names(ipd2) <- c("OS", "PFS")


#####################################################################################
# Fit survival models and compute survival curves
#####################################################################################

# Fit survival models using MLE and obtain summaries
dists <- "weibull"
surv_summ1 <- survfit_fun(dists, ipd1, t_h, c("OS", "PFS"))
surv_summ2 <- survfit_fun(dists, ipd2, t_h, c("OS", "PFS"))

## Weibull model parameters
theta1_os <- MASS::mvrnorm(n_sim,
                           surv_summ1$OS$fit$weibull[,1],
                           surv_summ1$OS$fit$weibull[,2:3])

theta2_os <- MASS::mvrnorm(n_sim,
                           surv_summ2$OS$fit$weibull[,1],
                           surv_summ2$OS$fit$weibull[,2:3])

theta1_pfs <- MASS::mvrnorm(n_sim,
                            surv_summ1$PFS$fit$weibull[,1],
                            surv_summ1$PFS$fit$weibull[,2:3])

theta2_pfs <- MASS::mvrnorm(n_sim,
                            surv_summ2$PFS$fit$weibull[,1],
                            surv_summ2$PFS$fit$weibull[,2:3])

# compute the proportion surviving at each model cycle
S_os1 <- apply(theta1_os, 1, function (x) pweibull(0:t_h, exp(x[1]), exp(x[2]), lower.tail = F))
S_os2 <- apply(theta2_os, 1, function (x) pweibull(0:t_h, exp(x[1]), exp(x[2]), lower.tail = F))
S_pfs1 <- apply(theta1_pfs, 1, function (x) pweibull(0:t_h, exp(x[1]), exp(x[2]), lower.tail = F))
S_pfs2 <- apply(theta2_pfs, 1, function (x) pweibull(0:t_h, exp(x[1]), exp(x[2]), lower.tail = F))


#####################################################################################
# Waning of the treatment effect
#####################################################################################
if(waning == TRUE) {
  
  # specify lognormal distributions for start and end times of the waning period
  wane_start_par <- lnorm_mom(wane_start, wane_start) # start of waning period
  wane_end_par <- lnorm_mom(wane_end-wane_start, wane_end-wane_start) # duration of waning period
  
  # sample waning start times
  wane_start_runif <- runif(n_sim,
                            0,
                            plnorm(t_h, wane_start_par$mu, sdlog = sqrt(wane_start_par$sigma2))
  )
  wane_start_times <- round(qlnorm(wane_start_runif, wane_start_par$mu, sqrt(wane_start_par$sigma2)))
  
  # sample waning end times
  wane_end_runif <- runif(n_sim,
                          0,
                          plnorm(t_h, wane_end_par$mu, sdlog = sqrt(wane_end_par$sigma2))
  )
  wane_end_times <- wane_start_times + round(qlnorm(wane_end_runif, wane_end_par$mu, sqrt(wane_end_par$sigma2)))
  wane_end_times[wane_end_times>t_h] <- t_h
  
  ### compute the hazard
  # Overall survival
  h_os1 <- diff(-log(rowMeans(S_os1))) # hweibull(1:t_h, exp(mean(theta1_os[,1])), exp(mean(theta1_os[,2]))) #
  h_os2 <- diff(-log(rowMeans(S_os2))) # hweibull(1:t_h, exp(mean(theta2_os[,1])), exp(mean(theta2_os[,2]))) #
  
  # Progression-free survival
  h_pfs1 <- diff(-log(rowMeans(S_pfs1))) # hweibull(1:t_h, exp(mean(theta1_pfs[,1])), exp(mean(theta1_pfs[,2]))) 
  h_pfs2 <- diff(-log(rowMeans(S_pfs2))) # hweibull(1:t_h, exp(mean(theta2_pfs[,1])), exp(mean(theta2_pfs[,2]))) 
  
  for (i in 1:n_sim) {
    
    #### overall survival
    # compute additive waning hazards
    h_wane_os1 <- wane_fun(h_os1, h_os2, wane_start_times[i], wane_end_times[i])
    add_h_wane_os1 <- h_wane_os1 - h_os1
    
    # add the "waning" hazards to the unadjusted hazards
    h_temp <- diff(-log(S_os1[,i])) + add_h_wane_os1
    
    # compute adjusted survival curves
    S_os1[,i] <- cumprod(c(1, exp(-h_temp)))
    
    #### progression-free survival
    # compute additive waning hazards
    h_wane_pfs1 <- wane_fun(h_pfs1, h_pfs2, wane_start_times[i], wane_end_times[i])
    add_h_wane_pfs1 <- h_wane_pfs1 - h_pfs1
    
    # add the "waning" hazards to the unadjusted hazards
    h_temp <- diff(-log(S_pfs1[,i])) + add_h_wane_pfs1
    
    # compute adjusted survival curves
    S_pfs1[,i] <- cumprod(c(1, exp(-h_temp)))
    
  }
  
  S_os1[is.nan(S_os1)] <- S_pfs1[is.nan(S_pfs1)] <- 0 # to avoid NaN due to extremely small numbers
  
}    

# adjust PFS curves to ensure PFS doesn't exceed OS
S_pfs1 <- matrix(mapply(min, S_pfs1, S_os1), nrow(S_pfs1), ncol(S_pfs1))
S_pfs2 <- matrix(mapply(min, S_pfs2, S_os2), nrow(S_pfs2), ncol(S_pfs2))

# plot survival curves
plot(x = c(0:(t_h)), y = rowMeans(S_os1), 
     type = "l", col = "2", lty = 1, xlim = c(0,t_h), ylim = c(0,1), 
     xlab="Months", ylab="Survival probability",  
     yaxt="n",  bty = "l", yaxs = "i", xaxs="i", yaxp = c(0,1,4), las = 1,  xaxp = c(0, t_h, 15)) 
axis(side = 2, las = 2, at = c(0,0.25,0.50,0.75,1), labels = c("0","0.25","0.50","0.75","1"))
lines(rowMeans(S_os2), x = 0:t_h, type = "l", col = "4", lty = 1)
lines(x = c(0:t_h), rowMeans(S_pfs1), type = "l", col = "2", lty = 2)
lines(x = c(0:t_h), rowMeans(S_pfs2), type = "l", col = "4", lty = 2)
legend(x="topright", legend = c("OS new treatment", "OS standard care", "PFS new treatment", "PFS standard care"),
       lty = c(1,1,2,2), col = c(2, 4, 2, 4), bty = "n", cex = 1, y.intersp=1)


#####################################################################################
# Utilities and costs
#####################################################################################

# quality of life weights 
u_pfs <- rbeta(n_sim, 80, 20)
u_pps <- rbeta(n_sim, 50, 50)

# drug costs (per cycle)
c_drug1 <- rep(14400, n_sim)/cycle_year #22000

# medical costs (per cycle)
c_med1 <- rgamma(n_sim, 500^2/250^2, 1/(250^2/500)) #rgamma(n_sim, 6000^2/3000^2, 1/(3000^2/6000))/cycle_year
c_med2 <- rgamma(n_sim, 500^2/250^2, 1/(250^2/500)) #rgamma(n_sim, 6000^2/3000^2, 1/(3000^2/6000))/cycle_year


#####################################################################################
# Compute prior net benefit
#####################################################################################

# apply discounting to the survival curves
v_d_r <- 1/(1+d_r)^seq.int(0, t_h/cycle_year, (1/cycle_year)) # vector of monthly discounting factors

S_os1_d <- S_os1*v_d_r
S_os2_d <- S_os2*v_d_r
S_pfs1_d <- S_pfs1*v_d_r
S_pfs2_d <- S_pfs2*v_d_r

# compute life-years
ly1 <- colSums(S_os1_d)/cycle_year
ly2 <- colSums(S_os2_d)/cycle_year

# compute QALYs
qalys_pfs1 <- colSums(mapply("*", as.data.frame(S_pfs1_d), u_pfs))/cycle_year
qalys_pfs2 <- colSums(mapply("*", as.data.frame(S_pfs2_d), u_pfs))/cycle_year

qalys_pps1 <- colSums(mapply("*", as.data.frame(S_os1_d-S_pfs1_d), u_pps))/cycle_year
qalys_pps2 <- colSums(mapply("*", as.data.frame(S_os2_d-S_pfs2_d), u_pps))/cycle_year

qalys_tot1 <- qalys_pfs1 + qalys_pps1
qalys_tot2 <- qalys_pfs2 + qalys_pps2

# compute costs
if(waning == TRUE) {
  costs_drug1 <- colSums(sapply(1:n_sim, function (x) { # drug costs applied only until start of waning period 
    S_pfs1_d[,x] * c(rep(c_drug1[x], wane_start), rep(0, t_h-wane_start+1))
  })
  ) 
} else {
  costs_drug1 <- colSums(mapply("*", as.data.frame(S_pfs1_d), c_drug1))
}

costs_med1 <- colSums(mapply("*", as.data.frame(S_os1_d), c_med1))
costs_med2 <- colSums(mapply("*", as.data.frame(S_os2_d), c_med2))

costs_tot1 <- costs_med1 + costs_drug1
costs_tot2 <- costs_med2 # + costs_drug2

# cost-effectiveness results
mean(costs_tot1) # mean total costs for intervention 1
mean(costs_tot2) # mean total costs for intervention 2
mean(qalys_tot1) # mean total QALYs for intervention 1
mean(qalys_tot2) # mean total QALYs for intervention 2

paste("ICER = ", round(mean(costs_tot1-costs_tot2)/mean(qalys_tot1-qalys_tot2),0), sep = "")


#####################################################################################
# (Partial Expected Value of Perfect Information
#####################################################################################
os_curves <- list(S_os1, S_os2)           
pfs_curves <- list(S_pfs1, S_pfs2)        
theta_os <- list(theta1_os, theta2_os)    
theta_pfs <- list(theta1_pfs, theta2_pfs) 

inb <- threshold * (qalys_tot1-qalys_tot2) - (costs_tot1 - costs_tot2) # incremental net monetary benefit

##### EVPI #####
################
evpi <- mean(pmax(inb,0)) - mean(inb)
evpi

##### partial EVPI #####
########################
evppi_os <- surv_evppi_fun(os_curves)
evppi_os

evppi_pfs <- surv_evppi_fun(pfs_curves)
evppi_pfs

evppi_os_pfs <- surv_evppi_fun(os_curves, pfs_curves)
evppi_os_pfs


#####################################################################################
# Expected Value of Sample Information
#####################################################################################
start_times_os <- list(rep(0,100),        # start times for patients at risk for OS
                       rep(0,100)) 
start_times_pfs <- list(rep(0,100),       # start times for patients at risk for PFS
                        rep(0,100)) 
add_fu <- c(12, 24, 36, 48, 60, 120)      # additional follow-up times


##### Overall survival #####
evsi_os_int <- surv_evsi_fun(
  inb = inb,
  os_curves = os_curves,
  start_times_os = start_times_os,
  add_fu = add_fu,
  method = "interpolation",
  seed = 1
)

evsi_os_disc <- surv_evsi_fun(
  inb = inb,
  os_curves = os_curves,
  start_times_os = start_times_os,
  add_fu = add_fu,
  method = "discrete",
  fast = FALSE,
  seed = 1
)

evsi_os_std <- surv_evsi_fun(
  inb = inb,
  os_curves = os_curves,
  theta_os = theta_os,
  start_times_os = start_times_os,
  add_fu = add_fu,
  method = "standard",
  seed = 1
)

##### Overall survival + progression-free survival #####
evsi_os_pfs_int <- surv_evsi_fun(
  inb = inb,
  os_curves = os_curves,
  pfs_curves = pfs_curves,
  start_times_os = start_times_os,
  start_times_pfs = start_times_pfs,
  add_fu = add_fu,
  method = "interpolation",
  seed = 1
)

evsi_os_pfs_disc <- surv_evsi_fun(
  inb,
  os_curves = os_curves,
  pfs_curves = pfs_curves,
  start_times_os = start_times_os,
  start_times_pfs = start_times_pfs,
  add_fu = add_fu,
  method = "discrete",
  fast = F,
  seed = 1
)

evsi_os_pfs_std <- surv_evsi_fun(
  inb,
  os_curves = os_curves,
  pfs_curves = pfs_curves,
  theta_os = theta_os,
  theta_pfs = theta_pfs,
  start_times_os = start_times_os,
  start_times_pfs = start_times_pfs,
  add_fu = add_fu,
  method = "standard",
  seed = 1
)


evsi_wane <- list("evpi" = evpi, "evppi_os" = evppi_os, "evppi_pfs" = evppi_pfs, "evppi_os_pfs" = evppi_os_pfs, 
                  "evsi_os_int"=evsi_os_int, "evsi_os_disc"=evsi_os_disc, "evsi_os_std"=evsi_os_std,
                  "evsi_os_pfs_int"=evsi_os_pfs_int, "evsi_os_pfs_disc"=evsi_os_pfs_disc, "evsi_os_pfs_std"=evsi_os_pfs_std)
                                 


surv_curves_wane <- list("S_os1"=S_os1, "S_os2"=S_os2, "S_pfs1"=S_pfs1, "S_pfs2"=S_pfs2, 
                         "theta_os1"=theta1_os,"theta2_os" =theta2_os, "theta1_pfs" = theta1_pfs, "theta2_pfs"=theta2_pfs)

ce_results_wane <- list("inb"=inb,
                        "costs_tot1"=costs_tot1, "costs_tot2"=costs_tot2, 
                        "qalys_tot1"=qalys_tot1, "qalys_tot2"=qalys_tot2, 
                        "ly1"=ly1, "ly2"=ly2)  






#####################################################################################
# Save results
#####################################################################################
save.image(file = here::here("output", paste0("output_", Sys.Date(), ".RData")))

# if(waning == TRUE) {
#   save.image(file = here::here("output", paste0("output_wane_", Sys.Date(), ".RData")))
# } else {
#   save.image(file = here::here("output", paste0("output_", Sys.Date(), ".RData")))
# }













#####################################################################################
# Plot the simulation results
#####################################################################################
rm(list=ls(all=TRUE))	

library(here)
load(here::here("output", "output_2022-05-05.RData"))

library(dplyr)

# create df with EVSI results
df_evsi_no_wane <- lapply(evsi_no_wane[5:10], t)
df_evsi_no_wane <- as.data.frame(do.call(rbind, df_evsi_no_wane))
df_evsi_no_wane$method <- rep(c("int", "disc", "std"), 2 , each = 6)
df_evsi_no_wane$data <- rep(c("os", "os_pfs"), 1, each = 18)
df_evsi_no_wane$wane <- "no"

names(df_evsi_no_wane) <- c("add. follow-up", "evsi", "se", "lower.2.5%", "upper.97.5%", "method", "data", "wane") 


# plot settings
cols=colorspace::rainbow_hcl(3)
inits <- seq(0.6,5.6, by=1)
space <- c(0,0.2,0.4)#c(0,0.25,0.5)
bottom <- 0
top <- ceiling(max(df_evsi_no_wane$`upper.97.5%`))
K_ <- n_sim
t_1 <- add_fu
labels <- c("Overall survival", "Overall- and progression-free survival")
s_data <- c("os", "os_pfs")

# create subsets of the data
df_sub <- vector(mode = "list", length = 2)
for (i in 1:2) {
    df_sub[[i]] <- df_evsi_no_wane %>% dplyr::filter(data == unique(df_evsi_no_wane$data)[i])
}





lwd_cust <- 2
pch_cust <- 8

dev.off()
par(tcl = -0.25) # The length of tick marks as a fraction of the height of a line of text. 
par(mgp = c(2, 0.6, 0)) # axis label locations
#par(mar = c(1.5,1.5,1.5,1.5), oma = c(4, 4, 0.5, 0.5))

par(mar = c(1.5,3.5,1.5,1.5), oma = c(4, 4, 0.5, 0.5)) # oma = outer margins
#m_layout <- rbind(c(1,1,1), do.call("rbind", replicate(8, c(2,3), simplify = FALSE)))
m_layout <- rbind(c(1,1),
                  c(2,3),
                  c(2,3),
                  c(2,3),
                  c(2,3),
                  c(2,3),
                  c(2,3),
                  c(2,3))
layout(mat = m_layout, heights = c(1, rep(0.5, 2) ))


plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("bottom",c("Standard","Interpolated","Discrete", "EVPPI"),
       pch=c(1,2,3, NA),
       lwd=2,
       col=c(cols[1],cols[2],cols[3], "black"),
       lty=c(1,1,1,2),
       bty = "n",
       horiz = T,
       pt.cex = 1,
       cex = 1,
       x.intersp=0.5,
       text.width=c(0.04,0.04,0.04, 0.03))


for (i in 1:length(df_sub)) {
  
  if (i %in% c(1,3)) {
    plot(1, main = c("Overall survival", "Overall- and progression-free survival")[i], cex.lab = 1.5, yaxt='n', ylab = "Follow-up", xlab = "EVSI", type = "n", xlim=c(bottom,12000),ylim=c(0.5,6.5)) # xlim=c(bottom,top)
  } else {
    plot(1, main = c("Overall survival", "Overall- and progression-free survival")[i], yaxt='n', ylab = "", xlab = "EVSI", type = "n", xlim=c(bottom,12000),ylim=c(0.5,6.5)) # xlim=c(bottom,top)
  }
  

abline(v = c(evsi_no_wane$evppi_os$evppi, evsi_no_wane$evppi_os_pfs$evppi)[i], lty = 2, lwd = lwd_cust)
  
    for(j in 1:length(add_fu)){
    
points(subset(df_sub[[i]]$evsi, df_sub[[i]]$method == "std")[j],rep(inits[j]+space[3],1),pch=1,lwd=lwd_cust,col=cols[1])
points(c(subset(df_sub[[i]]$`lower.2.5%`, df_sub[[i]]$method == "std")[j],subset(df_sub[[i]]$`upper.97.5%`, df_sub[[i]]$method == "std")[j]) ,rep(inits[j]+space[3],2),type="l",lwd=lwd_cust-1,col=cols[1])

points(subset(df_sub[[i]]$evsi, df_sub[[i]]$method == "int")[j],rep(inits[j]+space[2],1),pch=2,lwd=lwd_cust,col=cols[2], lty = 1)
points(c(subset(df_sub[[i]]$`lower.2.5%`, df_sub[[i]]$method == "int")[j],subset(df_sub[[i]]$`upper.97.5%`, df_sub[[i]]$method == "std")[j]) ,rep(inits[j]+space[2],2),type="l",lwd=lwd_cust-1,col=cols[2], lty= 1)

points(subset(df_sub[[i]]$evsi, df_sub[[i]]$method == "disc")[j],rep(inits[j]+space[1],1),pch=3,lwd=lwd_cust,col=cols[3], lty = 1)
points(c(subset(df_sub[[i]]$`lower.2.5%`, df_sub[[i]]$method == "disc")[j],subset(df_sub[[i]]$`upper.97.5%`, df_sub[[i]]$method == "std")[j]) ,rep(inits[j]+space[1],2),type="l",lwd=lwd_cust-1,col=cols[3], lty= 1)

#axis(2,inits[j]+mean(space),labels=paste("N=",N[j]))

  if (i %in% c(1,3)) {
    axis(2,inits[j]+mean(space),labels=paste(add_fu[j]/12), cex.axis = 1)
  }
}
}


mtext("EVSI", side = 1, outer = TRUE, cex = 1, line = 2.2,
      col = "black")


# 
# #if (i %in% c(7,8,9))
# axis(1,  at = seq(0, 12000, 12), cex.axis = 1) #col = "grey40", col.axis = "grey20",






dev.off()
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
#par(mar = c(1.5,1.5,1.5,1.5), oma = c(4, 4, 0.5, 0.5))

par(mar = c(1.5,1.5,1.5,1.5), oma = c(4, 4, 0.5, 0.5))
#m_layout <- rbind(c(1,1,1), do.call("rbind", replicate(8, c(2,3), simplify = FALSE)))
m_layout <- rbind(c(1,1),
                  c(2,3),
                  c(2,3),
                  c(2,3),
                  c(2,3))
layout(mat = m_layout, heights = c(1, rep(0.5, 2) ))


plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("bottom",c("Standard","Interpolated","Discrete"),pch=c(1,2,3),lwd=2,
       col=c(cols[1],cols[2],cols[3]),
       lty=c(1,1,1),
       bty = "n",
       horiz = T,
       pt.cex = 1,
       cex = 1,
       x.intersp=0.5,
       text.width=c(0.08,0.08,0.02))

for (i in 1:length(df_sub)) {
  
  plot(1, axes = FALSE, type = "n", xlim=c(bottom,12000),ylim=c(0.5,6.5)) # xlim=c(bottom,top)
  
  for(j in 1:length(add_fu)){
    
    points(subset(df_sub[[i]]$evsi, df_sub[[i]]$method == "std")[j],rep(inits[j]+space[3],1),pch=1,lwd=1,col=cols[1])
    points(c(subset(df_sub[[i]]$`lower.2.5%`, df_sub[[i]]$method == "std")[j],subset(df_sub[[i]]$`upper.97.5%`, df_sub[[i]]$method == "std")[j]) ,rep(inits[j]+space[3],2),type="l",lwd=1,col=cols[1])
    
    points(subset(df_sub[[i]]$evsi, df_sub[[i]]$method == "int")[j],rep(inits[j]+space[2],1),pch=2,lwd=1,col=cols[2], lty = 1)
    points(c(subset(df_sub[[i]]$`lower.2.5%`, df_sub[[i]]$method == "int")[j],subset(df_sub[[i]]$`upper.97.5%`, df_sub[[i]]$method == "std")[j]) ,rep(inits[j]+space[2],2),type="l",lwd=1,col=cols[2], lty= 1)
    
    points(subset(df_sub[[i]]$evsi, df_sub[[i]]$method == "disc")[j],rep(inits[j]+space[1],1),pch=3,lwd=1,col=cols[3], lty = 1)
    points(c(subset(df_sub[[i]]$`lower.2.5%`, df_sub[[i]]$method == "disc")[j],subset(df_sub[[i]]$`upper.97.5%`, df_sub[[i]]$method == "std")[j]) ,rep(inits[j]+space[1],2),type="l",lwd=1,col=cols[3], lty= 1)
    
    #axis(2,inits[j]+mean(space),labels=paste("N=",N[j]))
    
    if (i %in% c(1,3))
      axis(2,inits[j]+mean(space),labels=paste(add_fu[j]/12), cex.axis = 1)
    
    #if (i %in% c(7,8,9))
    axis(1,  at = seq(0, 12000, 12), cex.axis = 1) #col = "grey40", col.axis = "grey20",
    
    #grid(nx = NULL, ny = NA,
    #     lty = 3,      # Grid line type
    #     col = "lightgray", # Grid line color
    #     lwd = 0.5)  
    
    
  }
  mtext(labels[i], side = 3, line = -1.5, adj = 0.1, cex = 0.7,  col = "black") #letters[i]
  
  box(col = "grey20")
}

mtext("EVSI", side = 1, outer = TRUE, cex = 1, line = 2.2,
      col = "black")
# mtext("y axis", side = 2, outer = TRUE, cex = 1, line = 2.2,
#       col = "grey20")

















#####################################################################################
# Plot the simulation results
#####################################################################################
rm(list=ls(all=TRUE))	

library(here)
load(here::here("output", "output_2022-05-05.RData"))

library(dplyr)

# create df with EVSI results
df_evsi_no_wane <- lapply(evsi_no_wane[5:10], t)
df_evsi_no_wane <- as.data.frame(do.call(rbind, df_evsi_no_wane))
df_evsi_no_wane$method <- rep(c("int", "disc", "std"), 2 , each = 6)
df_evsi_no_wane$data <- rep(c("os", "os_pfs"), 1, each = 18)
df_evsi_no_wane$wane <- "no"

df_evsi_wane <- lapply(evsi_wane[5:10], t)
df_evsi_wane <- as.data.frame(do.call(rbind, df_evsi_wane))
df_evsi_wane$method <- rep(c("int", "disc", "std"), 2 , each = 6)
df_evsi_wane$data <- rep(c("os", "os_pfs"), 1, each = 18)
df_evsi_wane$wane <- "yes"

df_evsi <- rbind(df_evsi_wane, df_evsi_no_wane)
names(df_evsi) <- c("add. follow-up", "evsi", "se", "lower.2.5%", "upper.97.5%", "method", "data", "wane") 


# plot settings
cols=colorspace::rainbow_hcl(3)
inits <- seq(0.6,5.6, by=1)
space <- c(0,0.2,0.4)#c(0,0.25,0.5)
bottom <- 0
top <- ceiling(max(df_evsi$`upper.97.5%`))
K_ <- n_sim
t_1 <- add_fu
labels <- c("Overall survival, no waning", "Overall- and progression-free survival, no waning", 
            "Overall survival, no waning", "Overall- and progression-free survival, no waning")
s_wane <- c("no", "yes")
s_data <- c("os", "os_pfs")

# create subsets of the data
df_sub <- rep(list(vector(mode = "list", length = 2)), 2)
for (i in 1:2) {
  for (j in 1:2) {
    df_sub[[i]][[j]] <- df_evsi %>% dplyr::filter(wane == s_wane[i] & data == s_data[j])
  }
}
df_sub <- unlist(df_sub, recursive=FALSE)


dev.off()
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
#par(mar = c(1.5,1.5,1.5,1.5), oma = c(4, 4, 0.5, 0.5))

par(mar = c(1.5,1.5,1.5,1.5), oma = c(4, 4, 0.5, 0.5))
#m_layout <- rbind(c(1,1,1), do.call("rbind", replicate(8, c(2,3), simplify = FALSE)))
m_layout <- rbind(c(1,1),
                  c(2,3),
                  c(2,3),
                  c(2,3),
                  c(4,5),
                  c(4,5),
                  c(4,5))

layout(mat = m_layout, heights = c(1, rep(0.5,9)))


plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("bottom",c("Standard","Interpolated","Discrete"),pch=c(1,2,3),lwd=2,
       col=c(cols[1],cols[2],cols[3]),
       lty=c(1,1,1),
       bty = "n",
       horiz = T,
       pt.cex = 1,
       cex = 1,
       x.intersp=0.5,
       text.width=c(0.08,0.08,0.02))


for (i in 1:length(df_sub)) {
  
  plot(1, axes = FALSE, type = "n", xlim=c(bottom,top),ylim=c(0.5,5.5))
  
  for(j in 1:length(add_fu)){
    
    points(subset(df_sub[[i]]$evsi, df_sub[[i]]$method == "std")[j],rep(inits[j]+space[3],1),pch=1,lwd=1,col=cols[1])
    points(c(subset(df_sub[[i]]$`lower.2.5%`, df_sub[[i]]$method == "std")[j],subset(df_sub[[i]]$`upper.97.5%`, df_sub[[i]]$method == "std")[j]) ,rep(inits[j]+space[3],2),type="l",lwd=1,col=cols[1])
    
    points(subset(df_sub[[i]]$evsi, df_sub[[i]]$method == "int")[j],rep(inits[j]+space[2],1),pch=2,lwd=1,col=cols[2], lty = 1)
    points(c(subset(df_sub[[i]]$`lower.2.5%`, df_sub[[i]]$method == "int")[j],subset(df_sub[[i]]$`upper.97.5%`, df_sub[[i]]$method == "std")[j]) ,rep(inits[j]+space[2],2),type="l",lwd=1,col=cols[2], lty= 1)
    
    points(subset(df_sub[[i]]$evsi, df_sub[[i]]$method == "disc")[j],rep(inits[j]+space[1],1),pch=3,lwd=1,col=cols[3], lty = 1)
    points(c(subset(df_sub[[i]]$`lower.2.5%`, df_sub[[i]]$method == "disc")[j],subset(df_sub[[i]]$`upper.97.5%`, df_sub[[i]]$method == "std")[j]) ,rep(inits[j]+space[1],2),type="l",lwd=1,col=cols[3], lty= 1)
    
    #axis(2,inits[j]+mean(space),labels=paste("N=",N[j]))
    
    if (i %in% c(1,3))
      axis(2,inits[j]+mean(space),labels=paste(add_fu[j]/12), cex.axis = 1)
    
    #if (i %in% c(7,8,9))
    axis(1, col = "grey40", col.axis = "grey20", at = seq(0, top, 10), cex.axis = 1)
    
    #grid(nx = NULL, ny = NA,
    #     lty = 3,      # Grid line type
    #     col = "lightgray", # Grid line color
    #     lwd = 0.5)  
    
    
  }
  mtext(labels[i], side = 3, line = -1.5, adj = 0.1, cex = 0.7,  col = "black") #letters[i]
  
  box(col = "grey20")
}

mtext("EVSI", side = 1, outer = TRUE, cex = 1, line = 2.2,
      col = "black")
# mtext("y axis", side = 2, outer = TRUE, cex = 1, line = 2.2,
#       col = "grey20")

