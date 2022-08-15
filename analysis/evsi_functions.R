##############################################################################################
# Function for computing partial EVPI using GAM
##############################################################################################
surv_evppi_fun <- function (os_curves=NULL, pfs_curves=NULL) {
  
  ######### Overall survival only #########
  if(is.null(pfs_curves)) {
    summ_stat <- data.frame(sapply(os_curves, function (x) {
      colSums(x)
    }))
    arm_indic <- rep(1:length(os_curves), 1, each = 1)         # trial arm indicator
    arms <- unique(arm_indic)
  }
  
  ######### Porgression-free survival only #########
  if(is.null(os_curves)) {
    summ_stat <- data.frame(sapply(pfs_curves, function (x) {
      colSums(x)
    }))
    arm_indic <- rep(1:length(os_curves), 1, each = 1)         # trial arm indicator
    arms <- unique(arm_indic)
  }
  
  ######### Overall survival and progression-free survival #########
  if(is.null(os_curves) == FALSE & is.null(pfs_curves) == FALSE) {
    curves <- c(os_curves, pfs_curves)
    summ_stat <- data.frame(sapply(curves, function (x) {
      colSums(x)
    }))
    arm_indic <- rep(1:length(os_curves), 2, each = 1)         # trial arm indicator
    arms <- unique(arm_indic)
  }
  
  ######### Compute partial EVPI using GAM #########
  regr_model <- reg_mod_fun(summ_stat, arms, arm_indic) 
  evppi <- evsi_fun(inb, summ_stat, regr_model)
  names(evppi) <- c("evppi", "se", "lower", "upper")
  
  return(evppi)
}

##############################################################################################
# function for generating survival data and computing EVSI for OS and/or PFS
##############################################################################################
surv_evsi_fun <- function (inb,
                           os_curves=NULL, 
                           pfs_curves=NULL,
                           theta_os=NULL,
                           theta_pfs=NULL,
                           start_times_os=NULL, 
                           start_times_pfs=NULL,
                           add_fu,
                           pfs2os=NULL,
                           method = "interpolation",
                           fast = TRUE,
                           seed) {
  

  sapply(add_fu, function (t) { 
    
    ######### Overall survival only #########
    if(is.null(pfs_curves)) {

      arm_indic <- rep(1:length(os_curves), 1, each = 2)         # trial arm indicator for each survival curve
      arms <- unique(arm_indic)                                  # number of trial arms
      
      # generate survival data
      if(method == "interpolation") {
        summ_stat <- interpol_os_data_fun(os_curves, start_times_os, t, seed) 
      } 
      
      if(method == "discrete") {
        summ_stat <- discrete_os_data_fun(os_curves, start_times_os, t, seed, fast = fast)
      } 
      
      if(method == "standard") {
        cycles <-  0:(length(os_curves[[1]][,1])-1)                     # number of model cycles
        summ_stat <- weib_os_data_fun(theta_os, start_times_os, cycles, t, seed) 
      } 
      
    } 
    
    ######### Overall survival and progression-free survival #########
    if(is.null(os_curves) == FALSE & is.null(pfs_curves) == FALSE) { #& is.null(start_times_os) == FALSE
      
      arm_indic <- rep(1:length(os_curves), 2, each = 2)              # trial arm indicator for each survival curve
      arms <- unique(arm_indic)                                       # number of trial arms
      
      
      # generate survival data
      if(method == "interpolation") {
        summ_stat <- interpol_os_pfs_data_fun(os_curves, pfs_curves, start_times_os, start_times_pfs, t, seed) 
      } 
      
      if(method == "discrete") {
        summ_stat <- discrete_os_pfs_data_fun(os_curves, pfs_curves, start_times_os, start_times_pfs, t, seed, fast = fast)
      } 
      
      if(method == "standard") {
        cycles <-  0:(length(os_curves[[1]][,1])-1)                     # number of model cycles
        summ_stat <- weib_os_pfs_data_fun(theta_os, theta_pfs, start_times_os, start_times_pfs, cycles, t, seed) 
      } 

      # check for implausible OS-progression datasets (i.e. sum of observed progression times > time at risk for OS)
      err_index <- check_os_pfs_fun(summ_stat, start_times_os, start_times_pfs)
      uneven_num_os <- seq(1, (ncol(summ_stat) / 2), 2)
      
      for (i in 1:length(uneven_num_os)) {
        print(paste0("Implausible datasets intervention ", i, " = ", round(length(unlist(err_index[uneven_num_os[i]]))/n_sim,4)*100, "%"))
      }
      
    }

    ######### Compute EVSI using GAM #########
    regr_model <- reg_mod_fun(summ_stat, arms, arm_indic)
    evsi <- evsi_fun(inb, summ_stat, regr_model)
    
    print(c("add. follow-up" = round(t,2), "evsi" = round(evsi[[1]],2), "se" = round(evsi[[2]],2), "lower" = round(evsi[[3]],2), "upper" = round(evsi[[4]],2) )) 
    
    return(c("add. follow-up" = round(t,2), "evsi" = round(evsi[[1]],2), "se" = round(evsi[[2]],2), "lower" = round(evsi[[3]],2), "upper" = round(evsi[[4]],2) )) 
    
  })
  
}

##############################################################################################
# Function for computing EVSI using GAM
##############################################################################################
evsi_fun <- function (inb, summ_stat, regr_model) {
  
  library(mgcv)
  library(MASS)
  
  print(paste("Estimating posterior INB.."))
  
  f <- update(formula(inb ~ .), formula(paste(".~", regr_model)))
  #mod <- bam(f, data = data.frame(summ_stat), discrete=TRUE)
  mod <- gam(f, data = data.frame(summ_stat))
  
  # extract fitted values                  
  g_hat <- mod$fitted
  
  # compute EVSI
  evsi <- sum(abs(g_hat[g_hat<0]))/length(inb)
  
  ### compute standard errors
  
  # extract the basis function values
  Xstar <- model.matrix(mod)
  
  # extract coefficients
  beta <- mod$coef
  
  # covariance matrix
  v_cov <- vcov(mod)
  
  # sample from the parameter distributions
  set.seed(123)
  parameter_draws <- t(mvrnorm(2000, beta, v_cov))
  
  # from these draws, calculate draws from fitted values
  fitted_draws <- Xstar %*% parameter_draws
  
  # compute EVSI for each sample
  evsi_samples <- apply(fitted_draws, 2, function (x) {
    sum(abs(x[x<0]))/length(x)
  })
  se <- sd(evsi_samples)
  ci <- quantile(evsi_samples, c(0.025, 0.975))
  
  return(list("evsi"= round(evsi,2), "se"=round(se,2), "upper" = round(ci[1],2), "lower" = round(ci[2],2) ))
}


##############################################################################################
# Function to define a GAM regression model
##############################################################################################
reg_mod_fun <- function (summ_stat, arms, arm_indic) {
    
  var_names <- lapply(arms, function (x) {
    colnames(summ_stat[arm_indic == x])
  })
  regr_model <- lapply(var_names, function (x) {
    paste("te(", paste(x,collapse = ","), ",k=4)", sep = "")
  })
  
  regr_model <- paste(unlist(regr_model), collapse = "+")
  
}
