rm(list=ls(all=TRUE))	

#####################################################################################
# load functions
#####################################################################################
source(here::here("Analysis", "gen_functions.R"))
source(here::here("Analysis", "voi_functions.R"))
source(here::here("Analysis", "data_sim_functions.R"))

#####################################################################################
# Analysis settings 
#####################################################################################
set.seed(2)   # set the seed

n_sim <- 2e3               # number of simulations in the probabilistic analysis
cycle_year <- 12            # number of model cycles per year
t_h <- 30*cycle_year        # model time horizon
d_r <- 0.035                # annual discount rate for costs and effects

waning <- TRUE             # apply waning of treatment effect (TRUE/FALSE)
wane_start <- 2*cycle_year  # mean point at which waning of treatment effect starts
wane_end <- 4*cycle_year    # mean point at which waning of treatment effect ends (expected hazards are equal between arms from this point)

threshold <- 8e4            # opportunity cost/willingness to pay threshold

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
  
  start_time <- Sys.time()
  
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
  
  ### Overall survival
  
  # compute mean survival
  mean_os1 <- mean_surv_fun(0:t_h, theta1_os)
  mean_os2 <- mean_surv_fun(0:t_h, theta2_os)

  # compute the difference in hazards for mean survival between the treatment arms
  haz_mean_os1 <- -log(1-(1-(mean_os1[2:(t_h+1)]/mean_os1[1:t_h])))
  haz_mean_os2 <- -log(1-(1-(mean_os2[2:(t_h+1)]/mean_os2[1:t_h])))
  diff_haz_mean_os <- haz_mean_os2 - haz_mean_os1
  
  # compute the waning-adjusted survival curves
  S_os1 <- parallel::mclapply(1:n_sim, function (k) {
    sapply(1:t_h, function (i) {
      S_wane_os_fun(i, logshape = theta1_os[k,1], logscale=theta1_os[k,2],
                                 t_start=wane_start_times[k], t_end=wane_end_times[k], 
                                 diff_haz_mean=diff_haz_mean_os) 
    })
  })
  S_os1 <- rbind(rep(1, n_sim), matrix(unlist(S_os1), nrow = t_h, ncol = n_sim))
  
  ### Progression-free survival
  
  # compute mean survival
  mean_pfs1 <- mean_surv_fun(0:t_h, theta1_pfs)
  mean_pfs2 <- mean_surv_fun(0:t_h, theta2_pfs)
  
  # compute the difference in hazards for mean survival between the treatment arms
  haz_mean_pfs1 <- -log(1-(1-(mean_pfs1[2:(t_h+1)]/mean_pfs1[1:t_h])))
  haz_mean_pfs2 <- -log(1-(1-(mean_pfs2[2:(t_h+1)]/mean_pfs2[1:t_h])))
  diff_haz_mean_pfs <- haz_mean_pfs2 - haz_mean_pfs1
  
  
  # compute the waning-adjusted survival curves
  S_pfs1 <- parallel::mclapply(1:n_sim, function (k) {
    sapply(1:t_h, function (i) {
      S_wane_pfs_fun(i, logshape = theta1_pfs[k,1], logscale=theta1_pfs[k,2],
                    t_start=wane_start_times[k], t_end=wane_end_times[k], 
                    diff_haz_mean=diff_haz_mean_pfs) 
    })
  })
  S_pfs1 <- rbind(rep(1, n_sim), matrix(unlist(S_pfs1), nrow = t_h, ncol = n_sim))
  
  psa_time <- Sys.time() - start_time
  psa_time
  
} # end waning function

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
c_med1 <- rgamma(n_sim, 500^2/250^2, 1/(250^2/500)) # rgamma(n_sim, 6000^2/3000^2, 1/(3000^2/6000))/cycle_year
c_med2 <- rgamma(n_sim, 500^2/250^2, 1/(250^2/500)) # rgamma(n_sim, 6000^2/3000^2, 1/(3000^2/6000))/cycle_year


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
evppi_os$evppi

evppi_pfs <- surv_evppi_fun(pfs_curves)
evppi_pfs$evppi

evppi_os_pfs <- surv_evppi_fun(os_curves, pfs_curves)
evppi_os_pfs$evppi


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
  fast = TRUE, # set to false in order to use the same random numbers as for the other methods
  seed = 1
)

if(waning == TRUE) {
  evsi_os_std <- surv_evsi_fun(
    inb = inb,
    theta_os = theta_os,
    start_times_os = start_times_os,
    add_fu = add_fu,
    method = "standard_num",
    seed = 1,
    sdists = list("S_wane_os_fun", "Sweibull"), #S_wane_os_fun
    wane_start_times = wane_start_times, 
    wane_end_times = wane_end_times
  )
} else {
  evsi_os_std <- surv_evsi_fun(
    inb = inb,
    theta_os = theta_os,
    start_times_os = start_times_os,
    add_fu = add_fu,
    method = "standard_ana",
    seed = 1
  )
}

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
  fast = TRUE, # set to false in order to use the same random numbers as for the other methods
  seed = 1
)

if(waning == TRUE) {
  evsi_os_pfs_std<- surv_evsi_fun(
    inb = inb,
    os_curves = os_curves,
    pfs_curves = pfs_curves,
    theta_os = theta_os,
    theta_pfs = theta_pfs,
    start_times_os = start_times_os,
    start_times_pfs = start_times_pfs,
    add_fu = add_fu,
    method = "standard_num",
    seed = 1,
    sdists = list("S_wane_os_fun", "Sweibull"), # S_wane_os_fun
    wane_start_times = wane_start_times, 
    wane_end_times = wane_end_times
  )
} else {
  evsi_os_pfs_std <- surv_evsi_fun(
    inb,
    os_curves = os_curves,
    pfs_curves = pfs_curves,
    theta_os = theta_os,
    theta_pfs = theta_pfs,
    start_times_os = start_times_os,
    start_times_pfs = start_times_pfs,
    add_fu = add_fu,
    method = "standard_ana",
    seed = 1
  )
}


#####################################################################################
# Save results
#####################################################################################
# if(waning == TRUE) {
#   save.image(file = here::here("output", paste0("output_wane_", Sys.Date(), ".RData")))
# } else {
#   save.image(file = here::here("output", paste0("output_", Sys.Date(), ".RData")))
# }
