#--------------------------------------------------------------------------------------------#
# Waning-adjusted hazard and survivor functions
#--------------------------------------------------------------------------------------------#

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

# weibull survivor function that ignores additional unused arguments
Sweibull <- function (q, logshape, logscale, ...) {
  pweibull(q, shape = exp(logshape), scale = exp(logscale), lower.tail = F)
} 


#--------------------------------------------------------------------------------------------#
# OS functions
#--------------------------------------------------------------------------------------------#

##############################################################################################
# function for generating overall survival datasets using spline interpolation
##############################################################################################
interpol_os_data_fun <- function (curves, start_times, t, seed, ...) {
  
  gc()
  
  print(paste("Generating datasets.."))
  
  atrisk <- sapply(start_times, function (x) length(x))           # number of patients at risk at t1
  cycles <-  0:(length(os_curves[[1]][,1])-1)                     # number of model cycles
  n_sim <- ncol(os_curves[[1]])                                   # number of PA simulations
  arm_indic <- rep(1:length(curves), 2, each = 2)                 # trial arm indicator
  arms <- unique(arm_indic)
    
  set.seed(seed) # set the seed
  
  # Sample values from a uniform distribution for each survival curve
  rand_num <- lapply(1:length(curves), function (y)  {
    apply(curves[[y]], 2, function (x) runif(atrisk[[y]], min(x), spline(cycles, x, xout=(start_times[[y]]), ties = max, method="hyman")$y)) 
  })
  
  # Look up the sampled values in the vector of survival probabilities and record the corresponding survival times (i.e. model cycles)
  surv_times <- lapply(1:length(curves), function (y) {
    sapply(1:n_sim, function (x) {
      spline(round(curves[[y]][,x], 100), cycles, xout = rand_num[[y]][, x], ties = max, method="hyman")$y
    })
  })

  # compute censoring indicator for each dataset
  cens <- lapply(1:length(surv_times), function (y) {
    apply(surv_times[[y]], 2, function (x)
      ifelse(x > (start_times[[y]] + t), 0, 1))
  })
  
  # censor the survival times
  surv_times <- lapply(1:length(surv_times), function (y) {
    apply(surv_times[[y]], 2, function (x) pmin(x, (start_times[[y]] + t)))    
  })
  
  # create datasets
  datasets <- lapply(1:length(surv_times), function (y) {
    lapply(1:n_sim, function (x) data.frame(cbind(surv_times[[y]][,x], cens[[y]][,x])))
  })
  
  # Compute a summary statistic using the number of observed events and time at risk for each dataset
  summ_stat <- lapply(1:length(datasets), function (y) {
    sapply(datasets[[y]], function (x) {
      par <- c(sum(x[, 2]), (sum(x[, 1]) - sum(start_times[[y]])))
    })
  })
  summ_stat <- as.data.frame(t(do.call(rbind, summ_stat)))
  
} 



##############################################################################################
# function for generating overall survival datasets using discrete sampling
##############################################################################################
discrete_os_data_fun <- function (curves, start_times, t, seed, fast = TRUE, ...) {
  
  gc()
  
  print(paste("Generating datasets.."))
  
  atrisk <- sapply(start_times, function (x) length(x))           # number of patients at risk at t1
  cycles <-  0:(length(os_curves[[1]][,1])-1)                     # number of model cycles
  n_sim <- ncol(curves[[1]])                                      # number of PA simulations
  arm_indic <- rep(1:length(os_curves), 2, each = 2)              # trial arm indicator
  arms <- unique(arm_indic)
    
  # function to find the half-cycle time corresponding to the uniform value rounded down to the nearest survival probability
  minpositive_fun = function(x) which.min(x[x > 0]) - 0.5 # minimum of positive values
  #minabs_fun = function(x) which.min(abs(x)) - 0.5 # minimum of absolute values
  
  if(fast == TRUE) {
    
    set.seed(seed) # set the seed
  
    # sample model cycles (i.e. survival times) using the densities
    surv_times <- lapply(1:length(curves), function (y)  {
      apply(curves[[y]], 2, function (x)
        sample(
          0.5:(max(cycles) - 0.5),
          size = atrisk[[y]],
          replace = T,
          prob = abs(diff(x))
        ))
    })
  } else {

    set.seed(seed) # set the seed
    
    # Sample values from a uniform distribution for each survival curve
    rand_num <- lapply(1:length(curves), function (y)  {
      apply(curves[[y]], 2, function (x) runif(atrisk[[y]], 
                                             min(x), 
                                             x[round(start_times[[y]])+1])) 
    })
    
    # Look up the sampled values in the vector of survival probabilities and record the corresponding survival times (i.e. model cycles)
    surv_times <- lapply(1:length(curves), function (y) {
      sapply(1:n_sim, function (x) {
        curves_temp <- curves[[y]][,x]
        sapply(1:atrisk[[y]], function (n) {
          minpositive_fun(curves_temp - rand_num[[y]][, x][n])
        })
      })
    })
    
  }
  
  # compute censoring indicator for each dataset
  cens <- lapply(1:length(surv_times), function (y) {
    apply(surv_times[[y]], 2, function (x)
      ifelse(x > (start_times[[y]] + t), 0, 1))
  })
  
  # censor the survival times
  surv_times <- lapply(1:length(surv_times), function (y) {
    apply(surv_times[[y]], 2, function (x) pmin(x, (start_times[[y]] + t)))    
  })
  
  # create datasets
  datasets <- lapply(1:length(surv_times), function (y) {
    lapply(1:n_sim, function (x) data.frame(cbind(surv_times[[y]][,x], cens[[y]][,x])))
  })
  
  # Compute a summary statistic using the number of observed events and time at risk for each dataset
  summ_stat <- lapply(1:length(datasets), function (y) {
    sapply(datasets[[y]], function (x) {
      par <- c(sum(x[, 2]), (sum(x[, 1]) - sum(start_times[[y]])))
    })
  })
  summ_stat <- as.data.frame(t(do.call(rbind, summ_stat)))
  
} 


##############################################################################################
# function for generating overall survival datasets using standard ITS
##############################################################################################
weib_os_data_fun <- function (theta, start_times, t, seed, ...) {
  
  gc()
  
  print(paste("Generating datasets.."))
  
  atrisk <- sapply(start_times, function (x) length(x))       # number of patients at risk at t1 for OS
  cycles <-  0:nrow(theta[[1]])
  n_sim <- nrow(theta[[1]])                                   # number of PA simulations
  arm_indic <- rep(1:length(theta), 2, each = 2)              # trial arm indicator
  arms <- unique(arm_indic)
  
  set.seed(seed) # set the seed
  
  # Sample values from a uniform distribution for each survival curve
  rand_num <- lapply(1:length(theta), function (y)  {
    sapply(1:n_sim, function (x) {
      runif(
        atrisk[[y]], 
        pweibull(max(cycles), exp(theta[[y]][x,1]), exp(theta[[y]][x,2]), lower.tail = F),
        pweibull(start_times[[y]], exp(theta[[y]][x,1]), exp(theta[[y]][x,2]), lower.tail = F)
      )
    })
  })
  
  # plug the random numbers into the corresponding inverse survivor function
  surv_times <- lapply(1:length(theta), function (y) {
    sapply(1:n_sim, function (x) {
      qweibull(rand_num[[y]][,x], exp(theta[[y]][x,1]), exp(theta[[y]][x,2]), lower.tail = F)
    })
  })
  
  # compute censoring indicator for each dataset
  cens <- lapply(1:length(surv_times), function (y) {
    apply(surv_times[[y]], 2, function (x)
      ifelse(x > (start_times[[y]] + t), 0, 1))
  })
  
  # censor the survival times
  surv_times <- lapply(1:length(surv_times), function (y) {
    apply(surv_times[[y]], 2, function (x) pmin(x, (start_times[[y]] + t)))    
  })
  
  # create datasets
  datasets <- lapply(1:length(surv_times), function (y) {
    lapply(1:n_sim, function (x) data.frame(cbind(surv_times[[y]][,x], cens[[y]][,x])))
  })
  
  # Compute a summary statistic using the number of observed events and time at risk for each dataset
  summ_stat <- lapply(1:length(datasets), function (y) {
    sapply(datasets[[y]], function (x) {
      par <- c(sum(x[, 2]), (sum(x[, 1]) - sum(start_times[[y]])))
    })
  })
  summ_stat <- as.data.frame(t(do.call(rbind, summ_stat)))
  
} # end data generation function


##############################################################################################
# function for generating overall survival datasets using ITS with numeric methods
##############################################################################################
wane_weib_os_data_fun <- function (theta, start_times, t, seed, sdists, wane_start_times, wane_end_times, ...) {
  
  gc()
  
  print(paste("Generating datasets.."))
  
  atrisk <- sapply(start_times, function (x) length(x))       # number of patients at risk at t1 for OS
  cycles <-  0:nrow(theta[[1]])
  n_sim <- nrow(theta[[1]])                                   # number of PA simulations
  arm_indic <- rep(1:length(theta), 2, each = 2)              # trial arm indicator
  arms <- unique(arm_indic)
  
  set.seed(seed) # set the seed
  
  # Sample values from a uniform distribution for each value for theta
  rand_num <- lapply(1:length(theta), function (y)  {
    sdist_temp <- eval(parse(text=sdists[[y]]))
    sapply(1:n_sim, function (x) {
      runif(
        atrisk[[y]], 
        sdist_temp(max(cycles), logshape = theta[[y]][x,1], logscale = theta[[y]][x,2], t_start=wane_start_times[x],t_end=wane_end_times[x]),
        sdist_temp(pmax(start_times[[y]], 10e-8), logshape = theta[[y]][x,1], logscale = theta[[y]][x,2], t_start=wane_start_times[x],t_end=wane_end_times[x])
      )
    })
  })
  
  # plug the random numbers into the corresponding inverse survivor function
  surv_times <- lapply(1:length(theta), function (y) {
    sdist_temp <- eval(parse(text=sdists[[y]]))
    sapply(1:n_sim, function (x) {
      #print(c(y,x))
      sapply(rand_num[[y]][,x], function (z) {
        qgeneric(sdist_temp, p = z,
                 logshape=theta[[y]][x,1], logscale=theta[[y]][x,2],
                 t_start=wane_start_times[x],t_end=wane_end_times[x])
      })
    })
  })
  
  # compute censoring indicator for each dataset
  cens <- lapply(1:length(surv_times), function (y) {
    apply(surv_times[[y]], 2, function (x)
      ifelse(x > (start_times[[y]] + t), 0, 1))
  })
  
  # censor the survival times
  surv_times <- lapply(1:length(surv_times), function (y) {
    apply(surv_times[[y]], 2, function (x) pmin(x, (start_times[[y]] + t)))    
  })
  
  # create datasets
  datasets <- lapply(1:length(surv_times), function (y) {
    lapply(1:n_sim, function (x) data.frame(cbind(surv_times[[y]][,x], cens[[y]][,x])))
  })
  
  # Compute a summary statistic using the number of observed events and time at risk for each dataset
  summ_stat <- lapply(1:length(datasets), function (y) {
    sapply(datasets[[y]], function (x) {
      par <- c(sum(x[, 2]), (sum(x[, 1]) - sum(start_times[[y]])))
    })
  })
  summ_stat <- as.data.frame(t(do.call(rbind, summ_stat)))
  
} # end data generation function



#--------------------------------------------------------------------------------------------#
# OS and PFS functions
#--------------------------------------------------------------------------------------------#



############################################################################################################
# function for generating overall survival and progression-free survival datasets using spline interpolation
############################################################################################################
interpol_os_pfs_data_fun <- function (os_curves, pfs_curves, start_times_os, start_times_pfs, t, seed, ...) {
  
  gc()
  
  print(paste("Generating datasets.."))
  
  atrisk_os <- sapply(start_times_os, function (x) length(x))     # number of patients at risk at t1 for OS
  atrisk_pfs <- sapply(start_times_pfs, function (x) length(x))   # number of patients at risk at t1 for OS
  cycles <-  0:(length(os_curves[[1]][,1])-1)                     # number of model cycles
  n_sim <- ncol(os_curves[[1]])                                   # number of PA simulations
  arm_indic <- rep(1:length(os_curves), 2, each = 2)              # trial arm indicator
  arms <- unique(arm_indic)
  
  set.seed(seed) # set the seed
  
  ### Overall survival
  
  # Sample values from a uniform distribution for each survival curve
  rand_num <- lapply(1:length(os_curves), function (y)  {
    apply(os_curves[[y]], 2, function (x) runif(atrisk_os[[y]], min(x), spline(cycles, x, xout=(start_times_os[[y]]), ties = max, method="hyman")$y)) 
  })
  
  # Look up the sampled values in the vector of survival probabilities and record the corresponding survival times (i.e. model cycles)
  surv_times <- lapply(1:length(os_curves), function (y) {
    sapply(1:n_sim, function (x) {
      spline(round(os_curves[[y]][,x], 100), cycles, xout = rand_num[[y]][, x], ties = max, method="hyman")$y
    })
  })
  
  # compute censoring indicator for each dataset
  cens <- lapply(1:length(surv_times), function (y) {
    apply(surv_times[[y]], 2, function (x)
      ifelse(x > (start_times_os[[y]] + t), 0, 1))
  })
  
  # censor the survival times
  surv_times <- lapply(1:length(surv_times), function (y) {
    apply(surv_times[[y]], 2, function (x) pmin(x, (start_times_os[[y]] + t)))    
  })
  
  # create datasets
  datasets <- lapply(1:length(surv_times), function (y) {
    lapply(1:n_sim, function (x) data.frame(cbind(surv_times[[y]][,x], cens[[y]][,x])))
  })
  
  # Compute a summary statistic using the number of observed events and time at risk for each dataset
  summ_stat_os <- lapply(1:length(datasets), function (y) {
    sapply(datasets[[y]], function (x) {
      par <- c(sum(x[, 2]), (sum(x[, 1]) - sum(start_times_os[[y]])))
    })
  })
  
  ### Progression-free survival
  
  set.seed(seed+1) # set the seed
  
  # Sample values from a uniform distribution for each survival curve
  rand_num <- lapply(1:length(pfs_curves), function (y)  {
    apply(pfs_curves[[y]], 2, function (x) runif(atrisk_pfs[[y]], min(x), spline(cycles, x, xout=(start_times_pfs[[y]]), ties = max, method="hyman")$y)) 
  })
  
  # Look up the sampled values in the vector of survival probabilities and record the corresponding survival times (i.e. model cycles)
  surv_times <- lapply(1:length(pfs_curves), function (y) {
    sapply(1:n_sim, function (x) {
      spline(round(pfs_curves[[y]][,x], 100), cycles, xout = rand_num[[y]][, x], ties = max, method="hyman")$y
    })
  })
  
  # select the times at which the PFS curve crosses the OS curve
  t_pfs_cross_os <- lapply(1:length(pfs_curves), function (x) {
    (pfs_curves[[x]] !=  os_curves[[x]])[-1,]
  })
  
  # select the first time at which the PFS curve and OS curve cross
  t_pfs_cross_os <- lapply(1:length(pfs_curves), function (x) {
    apply(t_pfs_cross_os[[x]], 2, function (y) {
      
      temp <- which(y == FALSE)
      
      if (length(temp) > 1) { # only select the time cycle if OS = PFS for at least two subsequent time cycles
        min(temp[temp > 1])  
      } else {
        Inf
      }
      
    })
  })
  
  # compute censoring indicator for each dataset
  cens <- lapply(1:length(surv_times), function (y) {
    sapply(1:n_sim, function (x) {
      ifelse(surv_times[[y]][, x] > min((start_times_pfs[[y]] + t), t_pfs_cross_os[[y]][x]), 0, 1)
    })
  })
  
  # censor the survival times
  surv_times <- lapply(1:length(surv_times), function (y) {
    apply(surv_times[[y]], 2, function (x) pmin(x, (start_times_pfs[[y]] + t)))    
  })
  
  # create datasets
  datasets <- lapply(1:length(surv_times), function (y) {
    lapply(1:n_sim, function (x) data.frame(cbind(surv_times[[y]][,x], cens[[y]][,x])))
  })
  
  # Compute a summary statistic using the number of observed events and time at risk for each dataset
  summ_stat_prog <- lapply(1:length(datasets), function (y) {
    sapply(datasets[[y]], function (x) {
      par <- c(sum(x[, 2]), (sum(x[, 1]) - sum(start_times_pfs[[y]])))
    })
  })
  
  # combine summary statistics for OS and for progression in one dataframe
  summ_stat <- as.data.frame(t(do.call(rbind,  c(summ_stat_os, summ_stat_prog))))

} # end data generation function


#########################################################################################################
# function for generating overall survival and progression-free survival datasets using discrete sampling
#########################################################################################################
discrete_os_pfs_data_fun <- function (os_curves, pfs_curves, start_times_os, start_times_pfs, t, seed, fast = TRUE, ...) {
  
  gc()
  
  print(paste("Generating datasets.."))
  
  atrisk_os <- sapply(start_times_os, function (x) length(x))     # number of patients at risk at t1 for OS
  atrisk_pfs <- sapply(start_times_pfs, function (x) length(x))   # number of patients at risk at t1 for OS
  cycles <-  0:(length(os_curves[[1]][,1])-1)                     # number of model cycles
  n_sim <- ncol(os_curves[[1]])                                   # number of PA simulations
  arm_indic <- rep(1:length(os_curves), 2, each = 2)              # trial arm indicator
  arms <- unique(arm_indic)
  
  # function to find the half-cycle time corresponding to the uniform value rounded down to the nearest survival probability
  minpositive_fun <- function(x) which.min(x[x > 0]) - 0.5 # minimum of positive values
  #minabs_fun <- function(x) which.min(abs(x)) - 0.5 # minimum of absolute values
  
  if(fast == TRUE) {
    
    set.seed(seed) # set the seed
    
    # sample model cycles (i.e. survival times) using the densities
    surv_times <- lapply(1:length(os_curves), function (y)  {
      apply(os_curves[[y]], 2, function (x)
        sample(
          0.5:(max(cycles) - 0.5),
          size = atrisk_os[[y]],
          replace = T,
          prob = abs(diff(x))
        ))
    })
    
  } else {
    
    set.seed(seed) # set the seed
    
    ### Overall survival
    
    # Sample values from a uniform distribution for each survival curve
    rand_num <- lapply(1:length(os_curves), function (y)  {
      apply(os_curves[[y]], 2, function (x) runif(atrisk_os[[y]], 
                                                  min(x), 
                                                  x[round(start_times_os[[y]])+1])) 
    })
    
    
    # Look up the sampled values in the vector of survival probabilities and record the corresponding survival times (i.e. model cycles)
    surv_times <- lapply(1:length(os_curves), function (y) {
      sapply(1:n_sim, function (x) {
        curves_temp <- os_curves[[y]][,x]
        sapply(1:atrisk_os[[y]], function (n) {
          minpositive_fun(curves_temp - rand_num[[y]][, x][n])
        })
      })
    })
    
  }
  
  # compute censoring indicator for each dataset
  cens <- lapply(1:length(surv_times), function (y) {
    apply(surv_times[[y]], 2, function (x)
      ifelse(x > (start_times_os[[y]] + t), 0, 1))
  })
  
  # censor the survival times
  surv_times <- lapply(1:length(surv_times), function (y) {
    apply(surv_times[[y]], 2, function (x) pmin(x, (start_times_os[[y]] + t)))    
  })
  
  # create datasets
  datasets <- lapply(1:length(surv_times), function (y) {
    lapply(1:n_sim, function (x) data.frame(cbind(surv_times[[y]][,x], cens[[y]][,x])))
  })
  
  # Compute a summary statistic using the number of observed events and time at risk for each dataset
  summ_stat_os <- lapply(1:length(datasets), function (y) {
    sapply(datasets[[y]], function (x) {
      par <- c(sum(x[, 2]), (sum(x[, 1]) - sum(start_times_os[[y]])))
    })
  })
  
  ### Progression-free survival
  if(fast == TRUE) {
    
    set.seed(seed+1) # set the seed
    
    # sample model cycles (i.e. survival times) using the densities
    surv_times <- lapply(1:length(pfs_curves), function (y)  {
      apply(pfs_curves[[y]], 2, function (x)
        sample(
          0.5:(max(cycles) - 0.5),
          size = atrisk_pfs[[y]],
          replace = T,
          prob = abs(diff(x))
        ))
    })
    
  } else {
    
    set.seed(seed+1) # set the seed
    
    # Sample values from a uniform distribution for each survival curve
    rand_num <- lapply(1:length(pfs_curves), function (y)  {
      apply(pfs_curves[[y]], 2, function (x) runif(atrisk_pfs[[y]], 
                                                   min(x), 
                                                   x[round(start_times_pfs[[y]])+1])) 
    })
    
    # Look up the sampled values in the vector of survival probabilities and record the corresponding survival times (i.e. model cycles)
    surv_times <- lapply(1:length(pfs_curves), function (y) {
      sapply(1:n_sim, function (x) {
        curves_temp <- pfs_curves[[y]][,x]
        sapply(1:atrisk_pfs[[y]], function (n) {
          minpositive_fun(curves_temp - rand_num[[y]][, x][n])
        })
      })
    })
    
  }
  
  # select the times at which the PFS curve crosses the OS curve
  t_pfs_cross_os <- lapply(1:length(pfs_curves), function (x) {
    (pfs_curves[[x]] !=  os_curves[[x]])[-1,]
  })
  
  # select the first time at which the PFS curve and OS curve cross
  t_pfs_cross_os <- lapply(1:length(pfs_curves), function (x) {
    apply(t_pfs_cross_os[[x]], 2, function (y) {
      
      temp <- which(y == FALSE)
      
      if (length(temp) > 1) { # only select the time cycle if OS = PFS for at least two subsequent time cycles
        min(temp[temp > 1])  
      } else {
        Inf
      }
      
    })
  })
  
  # compute censoring indicator for each dataset
  cens <- lapply(1:length(surv_times), function (y) {
    sapply(1:n_sim, function (x) {
      ifelse(surv_times[[y]][, x] > min((start_times_pfs[[y]] + t), t_pfs_cross_os[[y]][x]), 0, 1)
    })
  })
  
  # censor the survival times
  surv_times <- lapply(1:length(surv_times), function (y) {
    apply(surv_times[[y]], 2, function (x) pmin(x, (start_times_pfs[[y]] + t)))    
  })
  
  # create datasets
  datasets <- lapply(1:length(surv_times), function (y) {
    lapply(1:n_sim, function (x) data.frame(cbind(surv_times[[y]][,x], cens[[y]][,x])))
  })
  
  # Compute a summary statistic using the number of observed events and time at risk for each dataset
  summ_stat_prog <- lapply(1:length(datasets), function (y) {
    sapply(datasets[[y]], function (x) {
      par <- c(sum(x[, 2]), (sum(x[, 1]) - sum(start_times_pfs[[y]])))
    })
  })
  
  # combine summary statistics for OS and for progression in one dataframe
  summ_stat <- as.data.frame(t(do.call(rbind,  c(summ_stat_os, summ_stat_prog))))
  
} # end data generation function


####################################################################################################
# function for generating overall survival and progression-free survival datasets using standard ITS
####################################################################################################
weib_os_pfs_data_fun <- function (theta_os, theta_pfs, start_times_os, start_times_pfs, t, seed, ...) {
  
  gc()
  
  print(paste("Generating datasets.."))
  
  atrisk_os <- sapply(start_times_os, function (x) length(x))     # number of patients at risk at t1 for OS
  atrisk_pfs <- sapply(start_times_pfs, function (x) length(x))   # number of patients at risk at t1 for OS
  cycles <-  0:(length(os_curves[[1]][,1])-1)                     # number of model cycles
  arm_indic <- rep(1:length(os_curves), 2, each = 2)              # trial arm indicator
  arms <- unique(arm_indic)
  
  set.seed(seed) # set the seed
  
  ### Overall survival
  
  # Sample values from a uniform distribution for each survival curve
  rand_num <- lapply(1:length(theta_os), function (y)  {
    sapply(1:n_sim, function (x) {
      runif(
        atrisk_os[[y]], 
        pweibull(max(cycles), exp(theta_os[[y]][x,1]), exp(theta_os[[y]][x,2]), lower.tail = F),
        pweibull(start_times_os[[y]], exp(theta_os[[y]][x,1]), exp(theta_os[[y]][x,2]), lower.tail = F)
      )
    })
  })
  
  # Look up the sampled values in the vector of survival probabilities and record the corresponding survival times (i.e. model cycles)
  surv_times <- lapply(1:length(theta_os), function (y) {
    sapply(1:n_sim, function (x) {
      qweibull(rand_num[[y]][,x], exp(theta_os[[y]][x,1]), exp(theta_os[[y]][x,2]), lower.tail = F)
    })
  })
  
  # compute censoring indicator for each dataset
  cens <- lapply(1:length(surv_times), function (y) {
    apply(surv_times[[y]], 2, function (x)
      ifelse(x > (start_times_os[[y]] + t), 0, 1))
  })
  
  # censor the survival times
  surv_times <- lapply(1:length(surv_times), function (y) {
    apply(surv_times[[y]], 2, function (x) pmin(x, (start_times_os[[y]] + t)))    
  })
  
  # create datasets
  datasets <- lapply(1:length(surv_times), function (y) {
    lapply(1:n_sim, function (x) data.frame(cbind(surv_times[[y]][,x], cens[[y]][,x])))
  })
  
  # Compute a summary statistic using the number of observed events and time at risk for each dataset
  summ_stat_os <- lapply(1:length(datasets), function (y) {
    sapply(datasets[[y]], function (x) {
      par <- c(sum(x[, 2]), (sum(x[, 1]) - sum(start_times_os[[y]])))
    })
  })
  
  ### Progression-free survival
  
  set.seed(seed+1) # set the seed
  
  # Sample values from a uniform distribution for each survival curve
  rand_num <- lapply(1:length(theta_pfs), function (y)  {
    sapply(1:n_sim, function (x) {
      runif(
        atrisk_pfs[[y]], 
        pweibull(max(cycles), exp(theta_pfs[[y]][x,1]), exp(theta_pfs[[y]][x,2]), lower.tail = F),
        pweibull(start_times_os[[y]], exp(theta_pfs[[y]][x,1]), exp(theta_pfs[[y]][x,2]), lower.tail = F)
      )
    })
  })
  
  # Look up the sampled values in the vector of survival probabilities and record the corresponding survival times (i.e. model cycles)
  surv_times <- lapply(1:length(theta_pfs), function (y) {
    sapply(1:n_sim, function (x) {
      qweibull(rand_num[[y]][,x], exp(theta_pfs[[y]][x,1]), exp(theta_pfs[[y]][x,2]), lower.tail = F)
    })
  })
  
  # select the times at which the PFS curve crosses the OS curve
  t_pfs_cross_os <- lapply(1:length(pfs_curves), function (x) {
    (pfs_curves[[x]] !=  os_curves[[x]])[-1,]
  })
  
  # select the first time at which the PFS curve and OS curve cross
  t_pfs_cross_os <- lapply(1:length(pfs_curves), function (x) {
    apply(t_pfs_cross_os[[x]], 2, function (y) {
      
      temp <- which(y == FALSE)
      
      if (length(temp) > 1) { # only select the time cycle if OS = PFS for at least two subsequent time cycles
        min(temp[temp > 1])  
      } else {
        Inf
      }
      
    })
  })
  
  # compute censoring indicator for each dataset
  cens <- lapply(1:length(surv_times), function (y) {
    sapply(1:n_sim, function (x) {
      ifelse(surv_times[[y]][, x] > min((start_times_pfs[[y]] + t), t_pfs_cross_os[[y]][x]), 0, 1)
    })
  })
  
  # censor the survival times
  surv_times <- lapply(1:length(surv_times), function (y) {
    apply(surv_times[[y]], 2, function (x) pmin(x, (start_times_pfs[[y]] + t)))    
  })
  
  # create datasets
  datasets <- lapply(1:length(surv_times), function (y) {
    lapply(1:n_sim, function (x) data.frame(cbind(surv_times[[y]][,x], cens[[y]][,x])))
  })
  
  # Compute a summary statistic using the number of observed events and time at risk for each dataset
  summ_stat_prog <- lapply(1:length(datasets), function (y) {
    sapply(datasets[[y]], function (x) {
      par <- c(sum(x[, 2]), (sum(x[, 1]) - sum(start_times_pfs[[y]])))
    })
  })
  
  # combine summary statistics for OS and for progression in one dataframe
  summ_stat <- as.data.frame(t(do.call(rbind,  c(summ_stat_os, summ_stat_prog))))
  
} # end data generation function


###########################################################################################################################
# function for generating overall survival and progression-free survival datasets using standard ITS with numerical methods
###########################################################################################################################
wane_weib_os_pfs_data_fun <- function (theta_os, theta_pfs, start_times_os, start_times_pfs, t, seed, sdists, wane_start_times, wane_end_times, ...) {
  
  gc()
  
  print(paste("Generating datasets.."))
  
  atrisk_os <- sapply(start_times_os, function (x) length(x))     # number of patients at risk at t1 for OS
  atrisk_pfs <- sapply(start_times_pfs, function (x) length(x))   # number of patients at risk at t1 for OS
  cycles <-  0:(length(os_curves[[1]][,1])-1)                     # number of model cycles
  n_sim <- nrow(theta_os[[1]])                                    # number of simulations
  arm_indic <- rep(1:length(theta_os), 2, each = 2)               # trial arm indicator
  arms <- unique(arm_indic)
  
  set.seed(seed) # set the seed
  
  ### Overall survival
  
  # Sample values from a uniform distribution for each value for theta
  rand_num <- lapply(1:length(theta_os), function (y)  {
    sdist_temp <- eval(parse(text=sdists[[y]]))
    sapply(1:n_sim, function (x) {
      runif(
        atrisk_os[[y]], 
        sdist_temp(max(cycles), logshape = theta_os[[y]][x,1], logscale = theta_os[[y]][x,2], t_start=wane_start_times[x],t_end=wane_end_times[x]),
        sdist_temp(pmax(start_times_os[[y]], 10e-8), logshape = theta_os[[y]][x,1], logscale = theta_os[[y]][x,2], t_start=wane_start_times[x],t_end=wane_end_times[x])
      )
    })
  })
  
  # evaluate the inverse survivor functions at the random uniform numbers
  surv_times <- lapply(1:length(theta_os), function (y) {
    sdist_temp <- eval(parse(text=sdists[[y]]))
    sapply(1:n_sim, function (x) {
      #print(c(y,x))
      sapply(rand_num[[y]][,x], function (z) {
        qgeneric(sdist_temp, p = z,
                 logshape=theta_os[[y]][x,1], logscale=theta_os[[y]][x,2],
                 t_start=wane_start_times[x], t_end=wane_end_times[x])
      })
    })
  })
  
  # compute censoring indicator for each dataset
  cens <- lapply(1:length(surv_times), function (y) {
    apply(surv_times[[y]], 2, function (x)
      ifelse(x > (start_times_os[[y]] + t), 0, 1))
  })
  
  # censor the survival times
  surv_times <- lapply(1:length(surv_times), function (y) {
    apply(surv_times[[y]], 2, function (x) pmin(x, (start_times_os[[y]] + t)))    
  })
  
  # create datasets
  datasets <- lapply(1:length(surv_times), function (y) {
    lapply(1:n_sim, function (x) data.frame(cbind(surv_times[[y]][,x], cens[[y]][,x])))
  })
  
  # Compute a summary statistic using the number of observed events and time at risk for each dataset
  summ_stat_os <- lapply(1:length(datasets), function (y) {
    sapply(datasets[[y]], function (x) {
      par <- c(sum(x[, 2]), (sum(x[, 1]) - sum(start_times_os[[y]])))
    })
  })
  
  ### Progression-free survival
  
  set.seed(seed+1) # set the seed
  
  # Sample values from a uniform distribution for each survival curve
  rand_num <- lapply(1:length(theta_pfs), function (y)  {
    sdist_temp <- eval(parse(text=sdists[[y]]))
    sapply(1:n_sim, function (x) {
      runif(
        atrisk_pfs[[y]], 
        sdist_temp(max(cycles), logshape = theta_pfs[[y]][x,1], logscale = theta_pfs[[y]][x,2], t_start=wane_start_times[x],t_end=wane_end_times[x]),
        sdist_temp(pmax(start_times_pfs[[y]], 10e-8), logshape = theta_pfs[[y]][x,1], logscale = theta_pfs[[y]][x,2], t_start=wane_start_times[x],t_end=wane_end_times[x])
      )
    })
  })
  
  # plug the random numbers into the corresponding inverse survivor function
  surv_times <- lapply(1:length(theta_pfs), function (y) {
    sdist_temp <- eval(parse(text=sdists[[y]]))
    sapply(1:n_sim, function (x) {
      sapply(rand_num[[y]][,x], function (z) {
        qgeneric(sdist_temp, p = z,
                 logshape=theta_pfs[[y]][x,1], logscale=theta_pfs[[y]][x,2],
                 t_start=wane_start_times[x],t_end=wane_end_times[x])
      })
    })
  })
  
  # select the times at which the PFS curve crosses the OS curve
  t_pfs_cross_os <- lapply(1:length(pfs_curves), function (x) {
    (pfs_curves[[x]] !=  os_curves[[x]])[-1,]
  })
  
  # select the first time at which the PFS curve and OS curve cross
  t_pfs_cross_os <- lapply(1:length(pfs_curves), function (x) {
    apply(t_pfs_cross_os[[x]], 2, function (y) {
      
      temp <- which(y == FALSE)
      
      if (length(temp) > 1) { # only select the time cycle if OS = PFS for at least two subsequent time cycles
        min(temp[temp > 1])  
      } else {
        Inf
      }
      
    })
  })
  
  # compute censoring indicator for each dataset
  cens <- lapply(1:length(surv_times), function (y) {
    sapply(1:n_sim, function (x) {
      ifelse(surv_times[[y]][, x] > min((start_times_pfs[[y]] + t), t_pfs_cross_os[[y]][x]), 0, 1)
    })
  })
  
  # censor the survival times
  surv_times <- lapply(1:length(surv_times), function (y) {
    apply(surv_times[[y]], 2, function (x) pmin(x, (start_times_pfs[[y]] + t)))    
  })
  
  # create datasets
  datasets <- lapply(1:length(surv_times), function (y) {
    lapply(1:n_sim, function (x) data.frame(cbind(surv_times[[y]][,x], cens[[y]][,x])))
  })
  
  # Compute a summary statistic using the number of observed events and time at risk for each dataset
  summ_stat_prog <- lapply(1:length(datasets), function (y) {
    sapply(datasets[[y]], function (x) {
      par <- c(sum(x[, 2]), (sum(x[, 1]) - sum(start_times_pfs[[y]])))
    })
  })
  
  # combine summary statistics for OS and for progression in one dataframe
  summ_stat <- as.data.frame(t(do.call(rbind,  c(summ_stat_os, summ_stat_prog))))
  
} # end data generation function


##############################################################################################
# Function for identify implausible combinations of OS-PFS datasets
##############################################################################################
check_os_pfs_fun <- function (summ_stat, start_times_os, start_times_pfs) {
  
  trunc_times <- c(start_times_os, start_times_pfs)
  
  even_num_os <- seq(2, (ncol(summ_stat) / 2), 2)
  uneven_num_os <- seq(1, (ncol(summ_stat) / 2), 2)
  even_num_pfs <- seq(((ncol(summ_stat) / 2)) + 2, ncol(summ_stat) , 2)
  uneven_num_pfs <- seq(((ncol(summ_stat) / 2) + 1), ncol(summ_stat) , 2)
  
  # identify implausible datasets where the sum of observed progression times exceeds total time at risk for OS
  err_trisk <- lapply(1:length(even_num_os), function (i) {
    which(summ_stat[, even_num_os[i]] < summ_stat[, even_num_pfs[i]])
  }) 
  
  # merge error indices
  err_index <- lapply(1:length(err_trisk), function (i) err_trisk[[i]])
  err_index <- rep(err_index, 2, each = 1)
  
}

