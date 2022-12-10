#####################################################################################
# Function for generating survival times evenly spaced from a Weibull distribution
#####################################################################################
gen_times_q_fun <- function(npat, 
                            shape_w, 
                            scale_w) {

  if(any(grepl("package:dplyr", search()))) detach("package:dplyr") #else message("dplyr not loaded")
  
  # moving average function to compute the values between the quantiles
  ma <- function(x, n = 2) as.numeric(filter(x, rep(1 / n, n), sides = 2))[-length(x)]
  
  quantile_values <- ma(seq(0, 1, length.out = npat + 1))
  
  tt <- qweibull(quantile_values, shape_w, scale_w)
  
  return(tt)
}



#####################################################################################
# Function for truncating and right-censoring a dataset
#####################################################################################
ipd_cens_fun <- function(times, 
                         tmin, 
                         tmax, 
                         treat) {
  
  tt <- subset(times, times > tmin) # select times larger than left-truncation point
  event <- ifelse(tt > tmax, 0, 1) # event indicator
  tt <- pmin(tt, tmax) # right-censor times at tmax
  trunc <- rep(tmin, length(tt)) # left-truncate times at tmin
  treat <- treat
  ipd <-  as.data.frame(cbind(tt, event, trunc, treat))
  return(ipd)
}


#####################################################################################
# Function to fit survival models
#####################################################################################
survfit_fun <- function (dists, 
                         ipd, 
                         t_horizon, 
                         names_scenarios) {  
  
  library(flexsurv)
  
  # create lists to store model summaries
  surv_summ <- vector(mode = "list", length = length(ipd))
  names(surv_summ) <- names_scenarios
  surv_summ <- lapply(surv_summ, function(x)
    x = vector(mode = "list", length = 4))
  for (i in 1:length(ipd)) {
    names(surv_summ[[i]]) <- c("fit", "aic", "Aw", "mean_survival")
  }
  
  for (i in 1:length(ipd)) {
    
    surv_summ[[i]]$fit <- vector(mode = "list", length = length(dists)) # list of models fits
    surv_summ[[i]]$aic <- matrix(NaN, length(dists), 1) # vectors of aic scores
    surv_summ[[i]]$Aw <- matrix(NaN, length(dists), 1) # vectors of aic weights
    surv_summ[[i]]$mean_survival <- matrix(NaN, length(dists), 1) # mean survival estimates
    surv_summ[[i]]$events <- matrix(NaN, length(dists), 1) # event numbers
    surv_summ[[i]]$trisk <- matrix(NaN, length(dists), 1) # time at risk

    names(surv_summ[[i]]$fit) <-
      names(surv_summ[[i]]$aic) <-
      names(surv_summ[[i]]$Aw) <-
      names(surv_summ[[i]]$mean_survival) <-
      names(surv_summ[[i]]$events) <-
      names(surv_summ[[i]]$trisk) <- c(dists)
  }  
  
  # fit the models
  for (j in 1:length(ipd)) {
    
    for (i in 1:length(dists)) {
      
      if(dists[i] == "exp" || dists[i] == "gengamma") {
        
        model <- flexsurvreg(Surv(tt, event) ~ 1, data = ipd[[j]], 
                             dist = dists[i], method = "BFGS")
        
      } else {
        
        model <- flexsurvreg(Surv(tt, event) ~ 1, data = ipd[[j]], 
                             dist = dists[i], method = "Nelder-Mead")
        
      }
      
      surv_summ[[j]]$fit[[i]] <- cbind(model$coefficients, model$cov)
      
      surv_summ[[j]]$aic[[i]] <- model$AIC
      
      surv_summ[[j]]$mean_survival[[i]] <- 
        eval(parse(text = paste0("rmst_", dists[i], "(", t_horizon, ",", 
                                 paste0(as.vector(model$res[,1]), collapse=","), ")")  )) 
      
      surv_summ[[j]]$events[[i]] <- model$events
      
      surv_summ[[j]]$trisk[[i]] <- model$trisk
      
    }
  }
  
  # compute AIC-based model weights
  aic_min <- lapply(surv_summ, function (x) {min(x$aic)} ) # identify the lowest AIC for each dataset
  
  # transform the AIC differences back to the scale of probabilities
  Ak <- lapply(1:length(ipd), 
               function (x) {exp(-0.5*(surv_summ[[x]]$aic-aic_min[[x]]))}) 
  # compute weights
  for (i in 1:length(ipd)) {
    surv_summ[[i]]$Aw <- apply(Ak[[i]], 2, function (y) {
      Ak[[i]] / sum(y)
    })
  }
  
  for (i in 1:length(ipd)) {
    names(surv_summ[[i]]$Aw) <- c(dists)
  }  
  
  # function output
  return(surv_summ)
}


#####################################################################################
# Function for sampling from a prior multivariate normal distribution
#####################################################################################
mvrnorm_fun <- function (surv_summ, 
                         dists, 
                         t_horizon, 
                         n, 
                         seed) {
  
  library(MASS)
  
  # set the seed
  set.seed(seed)
  
  # sample K survival distributions
  dist_select <- sample(x = dists, n, replace = T, prob = surv_summ$Aw)
  
  # matrix of transformed prior mean vectors
  mu <- lapply(dist_select, function (x) { 
    eval(parse(text = paste0("surv_summ$fit$", x, "[,1]"))) 
  })
  
  # matrix of transformed prior variance matrices
  cov_matrix <- lapply(dist_select, function (x) { 
    eval(parse(text = paste0("surv_summ$fit$", x, "[,2:","ncol(surv_summ$fit$", x,")]")))
  })
  
  # sample from priors on transformed scale
  prior_samples_trans <- lapply(1:n, function (x) {
    set.seed(seed + x)
    mvrnorm(1, mu = mu[[x]], Sigma = cov_matrix[[x]])
  })
  
  # transform parameter samples back to original scale
  prior_samples <- lapply(1:n, function (x) {
    sapply(1:length(prior_samples_trans[[x]]), function (y)  {
      eval(parse(text = paste0("flexsurv.dists$", dist_select[[x]],
                               "$inv.transforms[[y]](prior_samples_trans[[x]][y])")))
    })
  })
  
  prior_surv <- lapply(1:n, function (x) {
    eval(parse(text = paste0("rmst_", dist_select[[x]], "(", t_horizon, 
                             ",", paste(prior_samples[[x]], collapse = ",", sep = ","), ")")))
  })
  
  # return output
  output <- list(distribution = dist_select, theta = prior_samples, mean_surv = prior_surv)
  return(output)
}


##############################################################################################
# function computing meanlog and log-sd using method of moments
# https://devinincerti.com/2018/02/10/psa.html#gamma-and-lognormal-distributions
##############################################################################################
lnorm_mom <- function(mean, sd){
  if (mean > 0){
    sigma2 <- log((sd^2 + mean^2)/mean^2)
    mu <- log(mean) - 1/2 * sigma2
  } else{
    stop("Mean must be positive")
  }
  
  return(list(mu = mu, sigma2 = sigma2))
}


################################################################################################
# Waning-adjusted hazard and survivor functions
################################################################################################

# function for computing mean survival
mean_surv_fun <- function (x, theta) rowMeans(apply(theta, 1, function (y) pweibull(x, exp(y[1]), exp(y[2]), lower.tail = F ) ))

# function for spline interpolation of the difference in mean hazards between the treatment arms
spline_fun <- function (q, x, y) spline(x=x, y=y, xout = q)$y

# waning hazard function
whaz_fun <- function (q, logshape, logscale, t_start, t_end, diff_haz_mean) {
  hweibull(q, exp(logshape), exp(logscale)) +
    punif((q - t_start) / (t_end - t_start), 0, 1) *
    pmax(0,  spline_fun(q, x=1:t_h, y=diff_haz_mean))
}

# survivor function for OS
S_wane_os_fun <- function (q, logshape, logscale, t_start, t_end, diff_haz_mean = diff_haz_mean_os) {
  cum_haz_fun <- function (q) {integrate(whaz_fun, 0, q, logshape = logshape, logscale = logscale, t_start=t_start, t_end=t_end, diff_haz_mean = diff_haz_mean)$value}
  cum_haz_fun <- Vectorize(cum_haz_fun)
  cum_haz <- cum_haz_fun(q)
  surv <- exp(-cum_haz)
  return(surv)
}

# survivor function for PFS
S_wane_pfs_fun <- function (q, logshape, logscale, t_start, t_end, diff_haz_mean = diff_haz_mean_pfs) {
  cum_haz_fun <- function (q) {integrate(whaz_fun, 0, q, logshape = logshape, logscale = logscale, t_start=t_start, t_end=t_end, diff_haz_mean = diff_haz_mean)$value}
  cum_haz_fun <- Vectorize(cum_haz_fun)
  cum_haz <- cum_haz_fun(q)
  surv <- exp(-cum_haz)
  return(surv)
}

