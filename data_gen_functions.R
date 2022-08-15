
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
weib_os_data_fun <- function (theta, start_times, cycles, t, seed, ...) {
  
  gc()
  
  print(paste("Generating datasets.."))
  
  atrisk <- sapply(start_times, function (x) length(x))       # number of patients at risk at t1 for OS
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
weib_os_pfs_data_fun <- function (theta_os, theta_pfs, start_times_os, start_times_pfs, cycles, t, seed, ...) {
  
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
