
#==================== - File for iterative full optimization routines - ========================================
# - This file includes the functions which can be run in order to carry out the 1st stage of the 
# - optimization process, where the neg-loglikelihood function is optimized over all parameters simultaneously.
#===============================================================================================================


# - Blackburn-Sherris model with independent factors
## - Use of the DNK parameter estimates as starting values set by default
it_f_opt_BSi <- function(mu_bar, x0=c(6.960591e-03, 9.017154e-03, 5.091784e-03), delta=c(-1.305830e-06, -5.220474e-02, -1.013210e-01), kappa=c(1.162624e-02, 6.787268e-02, 5.061539e-03), sigma=exp(c(-6.806310, -6.790270, -7.559145)), r=exp(c( -3.327060e+01, -6.086479e-01, -1.553156e+01)), max_iter=10, tol_lik=10, opt_met = 'Nelder-Mead'){
  par_est_table <- matrix(NA, nrow=max_iter, ncol=length(c(x0, delta, kappa, sigma, r)) + 2)
  n_factors <- length(kappa)
  colnames(par_est_table) <- c(sprintf("x0_%d", c(1:n_factors)), sprintf("delta_%d", c(1:n_factors)), sprintf("kappa_%d", c(1:n_factors)), sprintf("sigma_%d", c(1:n_factors)), c("r1", "r2", "rc"), "log_lik", "Code")
  
  l_sigma <- log(sigma)
  l_r <- log(r)
  
  par_opt_uKD0 <- optim(c(x0, delta, kappa, l_sigma, l_r), nLL_BSi_uKD, mu_bar=mu_bar, gr = NULL, method = opt_met, hessian = FALSE, control=list(trace=TRUE, maxit = 10000))  
  log_lik <- par_opt_uKD0$value
  par_est_table[1, (1:length(c(x0, delta, kappa, sigma, r)))] <- par_opt_uKD0$par
  par_est_table[1, length(c(x0, delta, kappa, sigma, r))+1] <- - 0.5 * par_opt_uKD0$value - 0.5 * ncol(mu_bar) * nrow(mu_bar)
  par_est_table[1, length(c(x0, delta, kappa, sigma, r))+2] <- par_opt_uKD0$convergence
  
  iter_count <- 2
  
  repeat{
    par_opt_uKD0 <- optim(par_opt_uKD0$par, nLL_BSi_uKD, mu_bar=mu_bar, gr = NULL, method = opt_met, hessian = FALSE, control=list(trace=TRUE, maxit = 10000)) 
    par_est_table[iter_count, (1:length(c(x0, delta, kappa, l_sigma, r)))] <- par_opt_uKD0$par
    par_est_table[iter_count, length(c(x0, delta, kappa, sigma, r))+1] <- -0.5 * par_opt_uKD0$value - 0.5 * ncol(mu_bar) * nrow(mu_bar)
    par_est_table[iter_count, length(c(x0, delta, kappa, sigma, r))+2] <- par_opt_uKD0$convergence
    if ((abs(par_opt_uKD0$value - log_lik) < tol_lik) | (iter_count==max_iter) ){
      break
    }else{
      # - Update log-likelihood
      log_lik <- par_opt_uKD0$value 
      iter_count <- iter_count + 1
    }
  }
  
  par_est_table[,(n_factors * 3 + 1):(n_factors * 5)] <- exp(par_est_table[,(n_factors * 3 + 1):(n_factors * 5)])
  
  # - Return unconstrained parameter estimates
  x0_est <- par_opt_uKD0$par[1:n_factors]
  delta_est <- par_opt_uKD0$par[(n_factors + 1):(n_factors * 2)]
  kappa_est <- par_opt_uKD0$par[(n_factors * 2 + 1):(n_factors * 3)]
  sigma_est <- exp(par_opt_uKD0$par[(n_factors * 3 + 1):(n_factors * 4)])
  r1_est <- exp(par_opt_uKD0$par[n_factors * 4 + 1])
  r2_est <- exp(par_opt_uKD0$par[n_factors * 4 + 2])
  rc_est <- exp(par_opt_uKD0$par[n_factors * 4 + 3])
  
  return(list(par_est = list(x0=x0_est, delta=delta_est, kappa=kappa_est, sigma=sigma_est, r1=r1_est, r2=r2_est, rc=rc_est), log_lik = (-0.5 * par_opt_uKD0$value - 0.5 * ncol(mu_bar) * nrow(mu_bar)), par_table = par_est_table[1:iter_count,]))
}

## - Improve parameter estimate local search
f_opt_BSi_LS <- function(mu_bar, x0=c(6.960591e-03, 9.017154e-03, 5.091784e-03), delta=c(-1.305830e-06, -5.220474e-02, -1.013210e-01), kappa=c(1.162624e-02, 6.787268e-02, 5.061539e-03), sigma=exp(c(-6.806310, -6.790270, -7.559145)), r=exp(c( -3.327060e+01, -6.086479e-01, -1.553156e+01))){
  par_sv <- c(x0, delta, kappa, log(sigma), log(r))  
  n_factors <- length(kappa)
  par_opt_LS <- optim(par_sv, nLL_BSi_uKD, mu_bar=mu_bar, gr = NULL, method = "BFGS", hessian = TRUE, control=list(trace=TRUE, maxit = 10000))
  
  x0_est <- par_opt_LS$par[1:n_factors]
  delta_est <- par_opt_LS$par[(n_factors + 1):(n_factors * 2)]
  kappa_est <- par_opt_LS$par[(n_factors * 2 + 1):(n_factors * 3)]
  sigma_est <- exp(par_opt_LS$par[(n_factors * 3 + 1):(n_factors * 4)])
  r1_est <- exp(par_opt_LS$par[n_factors * 4 + 1])
  r2_est <- exp(par_opt_LS$par[n_factors * 4 + 2])
  rc_est <- exp(par_opt_LS$par[n_factors * 4 + 3])
  
  return(list(par_est = list(x0=x0_est, delta=delta_est, kappa=kappa_est, sigma=sigma_est, r1=r1_est, r2=r2_est, rc=rc_est), log_lik = (-0.5 * par_opt_LS$value - 0.5 * nrow(mu_bar) * ncol(mu_bar))))
}


# - Blackburn-Sherris model with three dependent factors
## - Use of the DNK parameter estimates as starting values set by default
it_f_opt_BSd_3F <- function(mu_bar, x0=c(2.191140e-03, -8.855686e-03, 2.711990e-02), delta=c(-5.175933e-02, 4.578315e-01, -5.175921e-02, -2.299199e-01, 1.383445e-02, -6.310253e-02), kappa=c(3.455255e-02, 1.075876e-02, 1.000030e-02), sigma_dg=c(6.789215e-04, 2.036748e-03, 1.875928e-03), Sigma_cov=c(-1.260778e-06, 1.194974e-06, -3.718267e-06), r=exp(c(-3.345631e+01, -6.015438e-01, -1.557244e+01)), max_iter=10, tol_lik=10, opt_met = 'Nelder-Mead'){
  n_factors <- 3
  # - Table parameter estimation construction
  par_est_table <- matrix(NA, nrow=max_iter, ncol=length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r)) + 2)
  colnames(par_est_table) <- c(sprintf("x0_%d", c(1:n_factors)), "delta_11", 'delta_21', 'delta_22', 'delta_31', 'delta_32', 'delta_33', sprintf("kappa_%d", c(1:n_factors)), 'sigma_11', 'sigma_21', 'sigma_22', 'sigma_31', 'sigma_32', 'sigma_33', c("r1", "r2", "rc"), "log_lik", 'Code')

  # - Transformation of parameters
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  l_r <- log(r)
  
  # - Start estimating the model with the supplied starting values
  par_opt_uKD0 <- optim(c(x0, delta, kappa, dg_l_Sigma_chol, odg_Sigma_chol, l_r), nLL_BSd_3F_uKD, mu_bar=mu_bar, gr = NULL, method = opt_met, hessian = TRUE, control=list(trace=TRUE, maxit = 10000))  
  log_lik <- par_opt_uKD0$value
  par_est_table[1, (1:length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r)))] <- par_opt_uKD0$par
  par_est_table[1, length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1] <- -0.5 * par_opt_uKD0$value - 0.5 * ncol(mu_bar) * nrow(mu_bar)
  par_est_table[1, length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+2] <- par_opt_uKD0$convergence
  par_est_table[1,13:18] <- parest2cov(par_est_table[1,13:15], par_est_table[1,16:18])
  par_est_table[1,19:21] <- exp(par_est_table[1,19:21])
  
  iter_count <- 2
  
  repeat{
    par_opt_uKD0 <- optim(par_opt_uKD0$par, nLL_BSd_3F_uKD, mu_bar=mu_bar, gr = NULL, method = opt_met, hessian = FALSE, control=list(trace=TRUE, maxit = 10000))  
    par_est_table[iter_count, (1:length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r)))] <- par_opt_uKD0$par
    par_est_table[iter_count, length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1] <- -0.5 * par_opt_uKD0$value - 0.5 * ncol(mu_bar) * nrow(mu_bar)
    par_est_table[iter_count, length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+2] <- par_opt_uKD0$convergence
    
    # - Return table with backtransformed parameters
    par_est_table[iter_count,13:18] <- parest2cov(par_est_table[iter_count,13:15], par_est_table[iter_count,16:18])
    par_est_table[iter_count,19:21] <- exp(par_est_table[iter_count,19:21])
    
    if ((abs(par_opt_uKD0$value - log_lik) < tol_lik) | (iter_count==max_iter) ){
      break
    }else{
      # - Update log-likelihood
      log_lik <- par_opt_uKD0$value 
      iter_count <- iter_count + 1
    }
  }
  
  # - Return unconstrained parameter estimates
  x0_est <- par_opt_uKD0$par[1:3]
  delta_est <- par_opt_uKD0$par[4:9]
  kappa_est <- par_opt_uKD0$par[10:12]
  
  Sigma_est <- list(sigma_11 = par_est_table[iter_count,13], sigma_21 = par_est_table[iter_count,14], sigma_22 = par_est_table[iter_count,15], sigma_31 = par_est_table[iter_count,16], sigma_32 = par_est_table[iter_count,17], sigma_33 = par_est_table[iter_count,18])
  
  r1_est <- exp(par_opt_uKD0$par[19])
  r2_est <- exp(par_opt_uKD0$par[20])
  rc_est <- exp(par_opt_uKD0$par[21])
  
  return(list(par_est = list(x0=x0_est, delta=delta_est, kappa=kappa_est, Sigma=Sigma_est, r1=r1_est, r2=r2_est, rc=rc_est), log_lik = (-0.5 * par_opt_uKD0$value - 0.5 * ncol(mu_bar) * nrow(mu_bar)), par_table = par_est_table[1:iter_count,]))
}


## - Improve parameter estimate local search
f_opt_BSd_3F_LS <- function(mu_bar, x0=c(2.191140e-03, -8.855686e-03, 2.711990e-02), delta=c(-5.175933e-02, 4.578315e-01, -5.175921e-02, -2.299199e-01, 1.383445e-02, -6.310253e-02), kappa=c(3.455255e-02, 1.075876e-02, 1.000030e-02), sigma_dg=c(6.789215e-04, 2.036748e-03, 1.875928e-03), Sigma_cov=c(-1.260778e-06, 1.194974e-06, -3.718267e-06), r=exp(c(-3.345631e+01, -6.015438e-01, -1.557244e+01))){
  # - Transformation of parameters
  n_factors <- 3
  
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  l_r <- log(r)
  
  par_sv <- c(x0, delta, kappa, dg_l_Sigma_chol, odg_Sigma_chol, l_r)  
  par_opt_LS <- optim(par_sv, nLL_BSd_3F_uKD, mu_bar=mu_bar, gr = NULL, method = "BFGS", hessian = TRUE, control=list(trace=TRUE, maxit = 10000))
  
  # - Return unconstrained parameter estimates
  x0_est <- par_opt_LS$par[1:3]
  delta_est <- par_opt_LS$par[4:9]
  kappa_est <- par_opt_LS$par[10:12]
  
  Sigma_est <- list(sigma_11 = parest2cov(par_opt_LS$par[13:15],par_opt_LS$par[16:18])[1], sigma_21 = parest2cov(par_opt_LS$par[13:15],par_opt_LS$par[16:18])[2], sigma_22 = parest2cov(par_opt_LS$par[13:15],par_opt_LS$par[16:18])[3], sigma_31 = parest2cov(par_opt_LS$par[13:15],par_opt_LS$par[16:18])[4], sigma_32 = parest2cov(par_opt_LS$par[13:15],par_opt_LS$par[16:18])[5], sigma_33 = parest2cov(par_opt_LS$par[13:15],par_opt_LS$par[16:18])[6])
  
  r1_est <- exp(par_opt_LS$par[19])
  r2_est <- exp(par_opt_LS$par[20])
  rc_est <- exp(par_opt_LS$par[21])
  
  return(list(par_est = list(x0=x0_est, delta=delta_est, kappa=kappa_est, Sigma=Sigma_est, r1=r1_est, r2=r2_est, rc=rc_est), log_lik = (-0.5 * par_opt_LS$value - 0.5 * nrow(mu_bar) * ncol(mu_bar))))
}


# - AFNS model with independent factors
it_f_opt_AFNSi <- function(mu_bar, x0=c(1.091714e-02, 1.002960e-02, -5.990785e-04), delta=-8.304334e-02, kappa=c(9.154603e-03, 1.067658e-02, 7.439715e-03), sigma=exp(c(-7.318991, -7.535594, -8.456025)), r=exp(c(-3.371775e+01, -5.887962e-01, -1.548729e+01)), max_iter=10, tol_lik=10, opt_met = 'Nelder-Mead'){
  n_factors <- 3
  par_est_table <- matrix(NA, nrow=max_iter, ncol=length(c(x0, delta, kappa, sigma, r)) + 2)
  colnames(par_est_table) <- c('x0_L', 'x0_S', 'x0_C', "delta", 'kappa_L', 'kappa_S', 'kappa_C', 'sigma_L', 'sigma_S', 'sigma_C', c("r1", "r2", "rc"), "log_lik", "Code")
  
  l_sigma <- log(sigma)
  l_r <- log(r)
  
  # - Start estimating the model with the supplied starting values
  par_opt_uKD0 <- optim(c(x0, delta, kappa, l_sigma, l_r), nLL_AFNSi_uKD, mu_bar=mu_bar, gr = NULL, method = opt_met, hessian = FALSE, control=list(trace=TRUE, maxit = 10000))  
  log_lik <- par_opt_uKD0$value
  par_est_table[1, (1:length(c(x0, delta, kappa, sigma, r)))] <- par_opt_uKD0$par
  par_est_table[1, length(c(x0, delta, kappa, sigma, r))+1] <- -0.5 * par_opt_uKD0$value - 0.5 * ncol(mu_bar) * nrow(mu_bar)
  par_est_table[1, length(c(x0, delta, kappa, sigma, r))+2] <- par_opt_uKD0$convergence
  
  iter_count <- 2
  
  repeat{
    par_opt_uKD0 <- optim(par_opt_uKD0$par, nLL_AFNSi_uKD, mu_bar=mu_bar, gr = NULL, method = opt_met, hessian = FALSE, control=list(trace=TRUE, maxit = 10000)) 
    par_est_table[iter_count, (1:length(c(x0, delta, kappa, l_sigma, r)))] <- par_opt_uKD0$par
    par_est_table[iter_count, length(c(x0, delta, kappa, l_sigma, r))+1] <- -0.5 * par_opt_uKD0$value - 0.5 * ncol(mu_bar) * nrow(mu_bar)
    par_est_table[iter_count, length(c(x0, delta, kappa, l_sigma, r))+2] <- par_opt_uKD0$convergence
    if ((abs(par_opt_uKD0$value - log_lik) < tol_lik) | (iter_count==max_iter) ){
      break
    }else{
      # - Update log-likelihood
      log_lik <- par_opt_uKD0$value 
      iter_count <- iter_count + 1
    }
  }
  
  par_est_table[,8:13] <- exp(par_est_table[,8:13])
  
  # - Return unconstrained parameter estimates
  x0_est <- par_opt_uKD0$par[1:n_factors]
  delta_est <- par_opt_uKD0$par[4]
  kappa_est <- par_opt_uKD0$par[5:7]
  sigma_est <- exp(par_opt_uKD0$par[8:10])
  r1_est <- exp(par_opt_uKD0$par[11])
  r2_est <- exp(par_opt_uKD0$par[12])
  rc_est <- exp(par_opt_uKD0$par[13])
  
  return(list(par_est = list(x0=x0_est, delta=delta_est, kappa=kappa_est, sigma=sigma_est, r1=r1_est, r2=r2_est, rc=rc_est), log_lik = (-0.5 * par_opt_uKD0$value - 0.5 * ncol(mu_bar) * nrow(mu_bar)), par_table = par_est_table[1:iter_count,]))
}

f_opt_AFNSi_LS <- function(mu_bar, x0=c(1.091714e-02, 1.002960e-02, -5.990785e-04), delta=-8.304334e-02, kappa=c(9.154603e-03, 1.067658e-02, 7.439715e-03), sigma=exp(c(-7.318991, -7.535594, -8.456025)), r=exp(c(-3.371775e+01, -5.887962e-01, -1.548729e+01))){
  n_factors <- 3
  par_sv <- c(x0, delta, kappa, log(sigma), log(r))  
  par_opt_LS <- optim(par_sv, nLL_AFNSi_uKD, mu_bar=mu_bar, gr = NULL, method = "BFGS", hessian = TRUE, control=list(trace=TRUE, maxit = 10000))
  
  x0_est <- par_opt_LS$par[1:3]
  delta_est <- par_opt_LS$par[4]
  kappa_est <- par_opt_LS$par[5:7]
  sigma_est <- exp(par_opt_LS$par[8:10])
  r1_est <- exp(par_opt_LS$par[11])
  r2_est <- exp(par_opt_LS$par[12])
  rc_est <- exp(par_opt_LS$par[13])
  
  return(list(par_est = list(x0=x0_est, delta=delta_est, kappa=kappa_est, sigma=sigma_est, r1=r1_est, r2=r2_est, rc=rc_est), log_lik = (-0.5 * par_opt_LS$value - 0.5 * nrow(mu_bar) * ncol(mu_bar))))
}


# - AFNS model with dependent factors
it_f_opt_AFNSd <- function(mu_bar, x0=c(9.582516e-03, 1.094110e-02, -1.503155e-03), delta=-7.487697e-02, kappa=c(1.389363e-02, 3.525542e-03, 3.004883e-03), sigma_dg=c(3.215422e-03, 2.625474e-03, 1.164715e-03), Sigma_cov=c(-8.328978e-06, -3.685028e-06, 3.036376e-06), r=exp(c(-3.335725e+01, -6.066149e-01, -1.552061e+01)), max_iter=10, tol_lik=10, opt_met = 'Nelder-Mead'){
  # - Table parameter estimation construction
  n_factors <- 3
  par_est_table <- matrix(NA, nrow=max_iter, ncol=length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r)) + 2)
  colnames(par_est_table) <- c('x0_L', 'x0_S', 'x0_C', "delta", 'kappa_L', 'kappa_S', 'kappa_C', 'sigma_L', 'sigma_LS', 'sigma_S', 'sigma_LC', 'sigma_SC', 'sigma_C', c("r1", "r2", "rc"), "log_lik", 'Code')
  
  # - Transformation of parameters
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  l_r <- log(r)
  
  # - Start estimating the model with the supplied starting values
  par_opt_uKD0 <- optim(c(x0, delta, kappa, dg_l_Sigma_chol, odg_Sigma_chol, l_r), nLL_AFNSd_uKD, mu_bar=mu_bar, gr = NULL, method = opt_met, hessian = TRUE, control=list(trace=TRUE, maxit = 10000))  
  log_lik <- par_opt_uKD0$value
  par_est_table[1, (1:length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r)))] <- par_opt_uKD0$par
  par_est_table[1, length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1] <- -0.5 * par_opt_uKD0$value - 0.5 * ncol(mu_bar) * nrow(mu_bar)
  par_est_table[1, length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+2] <- par_opt_uKD0$convergence
  par_est_table[1,8:13] <- parest2cov(par_est_table[1,8:10], par_est_table[1,11:13])
  par_est_table[1,14:16] <- exp(par_est_table[1,14:16])
  
  iter_count <- 2
  
  repeat{
    par_opt_uKD0 <- optim(par_opt_uKD0$par, nLL_AFNSd_uKD, mu_bar=mu_bar, gr = NULL, method = opt_met, hessian = FALSE, control=list(trace=TRUE, maxit = 10000))  
    par_est_table[iter_count, (1:length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r)))] <- par_opt_uKD0$par
    par_est_table[iter_count, length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1] <- -0.5 * par_opt_uKD0$value - 0.5 * ncol(mu_bar) * nrow(mu_bar)
    par_est_table[iter_count, length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+2] <- par_opt_uKD0$convergence
    
    # - Return table with backtransformed parameters
    par_est_table[iter_count,8:13] <- parest2cov(par_est_table[iter_count,8:10], par_est_table[iter_count,11:13])
    par_est_table[iter_count,14:16] <- exp(par_est_table[iter_count,14:16])
    
    if ((abs(par_opt_uKD0$value - log_lik) < tol_lik) | (iter_count==max_iter) ){
      break
    }else{
      # - Update log-likelihood
      log_lik <- par_opt_uKD0$value 
      iter_count <- iter_count + 1
    }
  }
  
  # - Return unconstrained parameter estimates
  x0_est <- par_opt_uKD0$par[1:3]
  delta_est <- par_opt_uKD0$par[4]
  kappa_est <- par_opt_uKD0$par[5:7]
  
  Sigma_est <- list(sigma_L = par_est_table[iter_count,8], sigma_LS = par_est_table[iter_count,9], sigma_S = par_est_table[iter_count,10], sigma_LC = par_est_table[iter_count,11], sigma_SC = par_est_table[iter_count,12], sigma_C = par_est_table[iter_count,13])
  
  r1_est <- exp(par_opt_uKD0$par[14])
  r2_est <- exp(par_opt_uKD0$par[15])
  rc_est <- exp(par_opt_uKD0$par[16])
  
  return(list(par_est = list(x0=x0_est, delta=delta_est, kappa=kappa_est, Sigma=Sigma_est, r1=r1_est, r2=r2_est, rc=rc_est), log_lik = (-0.5 * par_opt_uKD0$value - 0.5 * ncol(mu_bar) * nrow(mu_bar)), par_table = par_est_table[1:iter_count,]))
}


## - Improve parameter estimate local search
f_opt_AFNSd_LS <- function(mu_bar, x0=c(9.582516e-03, 1.094110e-02, -1.503155e-03), delta=-7.487697e-02, kappa=c(1.389363e-02, 3.525542e-03, 3.004883e-03), sigma_dg=c(3.215422e-03, 2.625474e-03, 1.164715e-03), Sigma_cov=c(-8.328978e-06, -3.685028e-06, 3.036376e-06), r=exp(c(-3.335725e+01, -6.066149e-01, -1.552061e+01))){
  n_factors <- 3
  
  # - Transformation of parameters
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  l_r <- log(r)
  
  par_sv <- c(x0, delta, kappa, dg_l_Sigma_chol, odg_Sigma_chol, l_r)  
  par_opt_LS <- optim(par_sv, nLL_AFNSd_uKD, mu_bar=mu_bar, gr = NULL, method = "BFGS", hessian = TRUE, control=list(trace=TRUE, maxit = 10000))
  
  # - Return unconstrained parameter estimates
  x0_est <- par_opt_LS$par[1:3]
  delta_est <- par_opt_LS$par[4]
  kappa_est <- par_opt_LS$par[5:7]
  
  Sigma_est <- list(sigma_L = parest2cov(par_opt_LS$par[8:10],par_opt_LS$par[11:13])[1], sigma_LS = parest2cov(par_opt_LS$par[8:10],par_opt_LS$par[11:13])[2], sigma_S = parest2cov(par_opt_LS$par[8:10],par_opt_LS$par[11:13])[3], sigma_LC = parest2cov(par_opt_LS$par[8:10],par_opt_LS$par[11:13])[4], sigma_SC = parest2cov(par_opt_LS$par[8:10],par_opt_LS$par[11:13])[5], sigma_C = parest2cov(par_opt_LS$par[8:10],par_opt_LS$par[11:13])[6])
  
  r1_est <- exp(par_opt_LS$par[14])
  r2_est <- exp(par_opt_LS$par[15])
  rc_est <- exp(par_opt_LS$par[16])
  
  return(list(par_est = list(x0=x0_est, delta=delta_est, kappa=kappa_est, Sigma=Sigma_est, r1=r1_est, r2=r2_est, rc=rc_est), log_lik = (-0.5 * par_opt_LS$value - 0.5 * nrow(mu_bar) * ncol(mu_bar))))
}

# - CIR model
it_f_opt_CIR <- function(mu_bar, x0=c(1.611524e-03, 5.763081e-03, 1.208483e-02), delta=c(-0.12379389, -0.06208546, -0.08131285), kappa=c(1.665062e-16, 3.477558e-01, 4.619791e-02), sigma=c(4.143351e-03, 6.242207e-02, 1.797287e-02), theta_Q = c(9.322613e-10, 8.457568e-03, 4.661882e-10), theta_P=c(0.01, 3.792994e-03, 6.272185e-03), r=c(2.952881e-15, 5.445661e-01, 1.493218e-07), max_iter=10, tol_lik=10, opt_met = 'Nelder-Mead'){
  n_factors <- length(kappa)
  par_est_table <- matrix(NA, nrow=max_iter, ncol=length(c(x0, delta, kappa, sigma, theta_Q, theta_P, r)) + 2)
  colnames(par_est_table) <- c(sprintf("x0_%d", c(1:n_factors)), sprintf("delta_%d", c(1:n_factors)), sprintf("kappa_%d", c(1:n_factors)), sprintf("sigma_%d", c(1:n_factors)), sprintf("theta_Q_%d", c(1:n_factors)), sprintf("theta_P_%d", c(1:n_factors)),  c("r1", "r2", "rc"), "log_lik", "Code")
  
  l_x0 <- log(x0)
  l_kappa <- log(kappa)
  l_sigma <- log(sigma)
  l_theta_Q <- log(theta_Q)
  l_theta_P <- log(theta_P)
  l_r <- log(r)
  
  par_opt_uKD0 <- optim(c(l_x0, delta, l_kappa, l_sigma, l_theta_Q, l_theta_P, l_r), nLL_CIR_uKD_bd, mu_bar=mu_bar, gr = NULL, method = opt_met, hessian = FALSE, control=list(trace=TRUE, maxit = 10000))  
  log_lik <- par_opt_uKD0$value
  par_est_table[1, (1:length(c(x0, delta, kappa, sigma, theta_Q, theta_P, r)))] <- par_opt_uKD0$par
  par_est_table[1, length(c(x0, delta, kappa, sigma, theta_Q, theta_P, r))+1] <- -0.5 * par_opt_uKD0$value - 0.5 * ncol(mu_bar) * nrow(mu_bar)
  par_est_table[1, length(c(x0, delta, kappa, sigma, theta_Q, theta_P, r))+2] <- par_opt_uKD0$convergence
  
  iter_count <- 2
  
  repeat{
    par_opt_uKD0 <- optim(par_opt_uKD0$par, nLL_CIR_uKD_bd, mu_bar=mu_bar, gr = NULL, method = opt_met, hessian = FALSE, control=list(trace=TRUE, maxit = 10000)) 
    par_est_table[iter_count, (1:length(c(x0, delta, kappa, sigma, theta_Q, theta_P, r)))] <- par_opt_uKD0$par
    par_est_table[iter_count, length(c(x0, delta, kappa, sigma, theta_Q, theta_P, r))+1] <- -0.5 * par_opt_uKD0$value - 0.5 * ncol(mu_bar) * nrow(mu_bar)
    par_est_table[iter_count, length(c(x0, delta, kappa, sigma, theta_Q, theta_P, r))+2] <- par_opt_uKD0$convergence
    if ((abs(par_opt_uKD0$value - log_lik) < tol_lik) | (iter_count==max_iter) ){
      break
    }else{
      # - Update log-likelihood
      log_lik <- par_opt_uKD0$value 
      iter_count <- iter_count + 1
    }
  }
  
  par_est_table[,c(1:n_factors, (n_factors * 2 + 1):(n_factors * 7))] <- exp(par_est_table[,c(1:n_factors, (n_factors * 2 + 1):(n_factors * 7))])
  
  # - Return unconstrained parameter estimates
  x0_est <- exp(par_opt_uKD0$par[1:n_factors])
  delta_est <- par_opt_uKD0$par[(n_factors + 1):(n_factors * 2)]
  kappa_est <- exp(par_opt_uKD0$par[(n_factors * 2 + 1):(n_factors * 3)])
  sigma_est <- exp(par_opt_uKD0$par[(n_factors * 3 + 1):(n_factors * 4)])
  theta_Q_est <- exp(par_opt_uKD0$par[(n_factors * 4 + 1):(n_factors * 5)])
  theta_P_est <- exp(par_opt_uKD0$par[(n_factors * 5 + 1):(n_factors * 6)])
  r1_est <- exp(par_opt_uKD0$par[n_factors * 6 + 1])
  r2_est <- exp(par_opt_uKD0$par[n_factors * 6 + 2])
  rc_est <- exp(par_opt_uKD0$par[n_factors * 6 + 3])
  
  return(list(par_est = list(x0=x0_est, delta=delta_est, kappa=kappa_est, sigma=sigma_est, theta_Q=theta_Q_est, theta_P=theta_P_est, r1=r1_est, r2=r2_est, rc=rc_est), log_lik = (-0.5 * par_opt_uKD0$value - 0.5 * ncol(mu_bar) * nrow(mu_bar)), par_table = par_est_table[1:iter_count,]))
}


## - Improve parameter estimate local search
f_opt_CIR_LS <- function(mu_bar, x0=c(1.611524e-03, 5.763081e-03, 1.208483e-02), delta=c(-0.12379389, -0.06208546, -0.08131285), kappa=c(1.665062e-16, 3.477558e-01, 4.619791e-02), sigma=c(4.143351e-03, 6.242207e-02, 1.797287e-02), theta_Q = c(9.322613e-10, 8.457568e-03, 4.661882e-10), theta_P=c(0.01, 3.792994e-03, 6.272185e-03), r=c(2.952881e-15, 5.445661e-01, 1.493218e-07)){
  n_factors <- length(kappa)
  par_sv <- c(log(x0), delta, log(kappa), log(sigma), log(theta_Q), log(theta_P), log(r))  
  par_opt_LS <- optim(par_sv, nLL_CIR_uKD_bd, mu_bar=mu_bar, gr = NULL, method = "BFGS", hessian = FALSE, control=list(trace=TRUE, maxit = 10000))
  
  x0_est <- exp(par_opt_LS$par[1:n_factors])
  delta_est <- par_opt_LS$par[(n_factors + 1):(n_factors * 2)]
  kappa_est <- exp(par_opt_LS$par[(n_factors * 2 + 1):(n_factors * 3)])
  sigma_est <- exp(par_opt_LS$par[(n_factors * 3 + 1):(n_factors * 4)])
  theta_Q_est <- exp(par_opt_LS$par[(n_factors * 4 + 1):(n_factors * 5)])
  theta_P_est <- exp(par_opt_LS$par[(n_factors * 5 + 1):(n_factors * 6)])
  r1_est <- exp(par_opt_LS$par[n_factors * 6 + 1])
  r2_est <- exp(par_opt_LS$par[n_factors * 6 + 2])
  rc_est <- exp(par_opt_LS$par[n_factors * 6 + 3])
  
  return(list(par_est = list(x0=x0_est, delta=delta_est, kappa=kappa_est, sigma=sigma_est, theta_Q=theta_Q_est, theta_P=theta_P_est, r1=r1_est, r2=r2_est, rc=rc_est), log_lik = (-0.5 * par_opt_LS$value - 0.5 * nrow(mu_bar) * ncol(mu_bar))))
}













