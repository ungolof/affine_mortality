###############################################################################################
# - Code to obtain the standard errors of the parameter estimates.
# - This implementation simplifies the ordinary one, since it obtains x0 from the 
# - Smoothing distribution. In this way, we account for the variability in parameter
# - estimates due to unknown x0, but we avoid to optimize also over x0, this way simplifying
# - the optimization procedure
###############################################################################################

# - Functions to get the standardized residuals

# - Blackburn-Sherris independent factor model
std_res_BSi <- function(x0, delta, kappa, sigma, r, mu_bar){
  
  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])
  
  n_factors <- length(kappa)
  
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years
  
  v_ti <- matrix(0, nrow(mu_bar), ncol(mu_bar))
  F_ti <- matrix(1, nrow(mu_bar), ncol(mu_bar))
  
  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance
  
  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSi(age, sigma, delta)  
    B_tT[age,] <- B_BSi(age,delta)  
  }
  
  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  R <- diag(((sigma^2) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)
  
  for(t in 1:n_years){
    
    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1,1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))   
    
    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t] 
      
      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) 
      
      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i 
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i,1] - B_tT[i,] %*% x_ti
      
    }
  }
  
  e_hat <- v_ti / sqrt(F_ti)
  
  return(e_hat)
}

y_star_BSi <- function(x0, delta, kappa, sigma, r, e_star){
  
  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])
  
  n_factors <- length(kappa)
  
  n_ages <- nrow(e_star)   # - Number of ages
  n_years <- ncol(e_star)  # - Number of years
  
  y_star <- matrix(0, nrow(e_star), ncol(e_star))
  F_ti <- matrix(1, nrow(e_star), ncol(e_star))
  
  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance
  
  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSi(age, exp(l_sigma), delta)  
    B_tT[age,] <- B_BSi(age,delta)  
  }
  
  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  R <- diag(((sigma^2) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)
  
  for(t in 1:n_years){
    
    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))    
    y_star[1,t] <- A_tT[1,1] + B_tT[1,] %*% x_ti + sqrt(F_ti[1,t]) * e_star[1,t]
    
    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / sqrt(F_ti[i-1,t])) %*% e_star[i-1,t] 
      
      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) 
      
      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i 
      y_star[i,t] <- A_tT[i,1] + B_tT[i,] %*% x_ti + sqrt(F_ti[i,t]) * e_star[i,t]
      
    }
  }
  
  return(y_star)
}

nLL_BSi_uKD_BS <- function(vdParameters, mu_bar_BS){
  
  n_factors <- length(vdParameters) / 4
  
  delta <- vdParameters[1:n_factors]
  kappa <- vdParameters[(n_factors+1):(n_factors*2)]
  l_sigma <- vdParameters[(n_factors*2+1):(n_factors*3)] # - take logs to ensure positivity
  r_c <- vdParameters[n_factors*3 + 3]
  r_1 <- vdParameters[n_factors*3 + 1]
  r_2 <- vdParameters[n_factors*3 + 2]
  
  n_ages <- nrow(mu_bar_BS)   # - Number of ages
  n_years <- ncol(mu_bar_BS)  # - Number of years
  
  v_ti <- mu_bar_BS
  F_ti <- mu_bar_BS
  
  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance
  
  # - Initialize X and Sigma
  x_ti <- x0_s #init_X
  P_ti <- 1e-10 * diag(1, n_factors)
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSi(age, exp(l_sigma), delta)  
    B_tT[age,] <- B_BSi(age,delta)  
  }
  
  Phi <- diag(exp(-kappa), n_factors) 
  R <- diag(((exp(l_sigma*2)) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)
  
  for(t in 1:n_years){
    
    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar_BS[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))   
    
    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t] 
      
      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) 
      
      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i 
      v_ti[i,t] <- mu_bar_BS[i,t] - A_tT[i] - B_tT[i,] %*% x_ti
      
    }
  }
  
  log_lik <- sum(log(F_ti) + (v_ti^2) / F_ti)
  return(log_lik)
}


## - Boostrap estimation of the standard errors (Covariance estimation)

CovEst_BS_BSi <- function(x0, delta, kappa, sigma, r, mu_bar, n_BS=500, t_ex = 4){
  par_table_BSi <- matrix(NA, n_BS, length(c(delta, kappa, sigma, r)))
  colnames(par_table_BSi) <- c(sprintf("delta_%d", c(1:n_factors)), sprintf("kappa_%d", c(1:n_factors)), sprintf("sigma_%d", c(1:n_factors)), c("r1", "r2", "rc"))
  
  res_table <- matrix(NA, n_ages, n_years)
  
  # - 4) Parameter estimation
  ## - 4.1) Get filtered estimates
  Filtering <- KF_BSi_uKD(x0, delta, kappa, sigma, r, mu_bar)
  X_t_fil <- Filtering$X_t
  X_t_c_fil <- Filtering$X_t_c
  S_t_fil <- Filtering$S_t
  S_t_c_fil <- Filtering$S_t_c
  
  ## - 4.2) Get smoothed estimate of x0
  Smoothing <- RTS_sm_bas(X_t_fil, X_t_c_fil, S_t_fil, S_t_c_fil, kappa)
  X_t_sm <- Smoothing$X_t_sm[,1]
  S_t_sm <- Smoothing$S_t_sm[,1:n_factors]
  
  # - 1) Get standardized residuals
  std_r <- std_res_BSi(x0, delta, kappa, sigma, r, mu_bar)
  ## - Fill the first t_ex residuals, which stay the same
  res_table[,1:t_ex] <- std_r[,1:t_ex]
  
  # - Repeat bootstrapping samples
  for(i in 1:n_BS){
    
    # - 2) resample standardized residuals
    t_sample <- sample(c((t_ex+1):n_years), (n_years-t_ex), replace=TRUE)  # - starting from 5 to account for KF start up irregularity
    res_table[,(t_ex + 1):n_years] <- std_r[,t_sample]
    
    ## - 4.3) Sample x0 from the smoothing distribution
    x0_s <- mvrnorm(n = 1, mu=X_t_sm, Sigma = S_t_sm, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    
    # - 3) Get mu_bar_star
    mu_bar_BS <- y_star_BSi(x0_s, delta, kappa, sigma, r, res_table)
    
    ## - 4.4) Get bootstrapped parameter estimates
    par_table_BSi[i,] <- optim(c(delta, kappa, log(sigma), log(r)), nLL_BSi_uKD_BS, mu_bar_BS=mu_bar_BS, gr = NULL, method = "BFGS", hessian = FALSE, control=list(trace=TRUE, maxit = 10000))$par # or coordinate ascent
    ## - Backtransformation
    par_table_BSi[i,((2 * n_factors + 1):(n_factors*4))] <- exp(par_table_BSi[i,((2 * n_factors + 1):(n_factors*4))])
  }

  # 5) Get st. err. parameter estimates (unconstrained parameters)
  cov_pe <- cov(par_table_BSi)
  serr_pe <- sqrt(diag(cov_pe_BSi))
  
return(list(Cov = cov_pe, St.err = serr_pe))
}


# - Blackburn-Sherris dependent factor model
std_res_BSd_3F <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar){
  
  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])
  
  n_factors <- 3
  
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years
  
  delta_matrix <- low_trg_fill(delta)
  
  v_ti <- matrix(0, nrow(mu_bar), ncol(mu_bar))
  F_ti <- matrix(1, nrow(mu_bar), ncol(mu_bar))
  
  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance
  
  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)
  
  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)
  
  # - Get Sigma (covariance matrix of the diffusion process)
  Sigma_diff <- Low_chol %*% t(Low_chol) # matrix(0, n_factors, n_factors)
  
  # - Get R (covariance of the state variable)
  for(row in 1:n_factors){
    for(col in 1:n_factors){
      R[row,col] <- Sigma_diff[row,col] * (1 - exp(- kappa[row] - kappa[col])) / (kappa[row] + kappa[col])
    }
  }
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSd_3F(age, Low_chol, delta_matrix)####### A_ind(age, exp(l_sigma), delta)  
    B_tT[age,] <- B_BSd_3F(age, delta_matrix)  ###### B_ind(age,delta)  
  }
  
  for(t in 1:n_years){
    
    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))    #max(0.000001, r_c + r_1 * exp(r_2))#exp(r_1 + r_2)  In case of the error specification of Xu et al. (2019)
    
    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t] # * (v_ti[i-1,t] / F_ti[i-1,t])   ### - Same change from + to minus following Xu, Sherris and Ziveyi 2019SAJ               # - a_{t,i+1}
      
      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) #exp(r_1 + r_2 * (i-1)) %*% t(K_ti) ####(max(0.001, r_c + r_1 * sum(exp(r_2 * c(1:(i-1)))) / (i-1))) %*% t(K_ti)
      
      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- t(B_tT[i,]) %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i) 
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti
      
    }
  }
  
  e_hat <- v_ti / sqrt(F_ti)
  
  return(e_hat)
}

y_star_BSd_3F <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, e_star){
  
  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])
  
  n_factors <- 3
  
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years
  
  delta_matrix <- low_trg_fill(delta)
  
  y_star <- matrix(0, nrow(mu_bar), ncol(mu_bar))
  F_ti <- matrix(1, nrow(mu_bar), ncol(mu_bar))
  
  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance
  
  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)
  
  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)
  
  # - Get Sigma (covariance matrix of the diffusion process)
  Sigma_diff <- Low_chol %*% t(Low_chol) # matrix(0, n_factors, n_factors)
  
  # - Get R (covariance of the state variable)
  for(row in 1:n_factors){
    for(col in 1:n_factors){
      R[row,col] <- Sigma_diff[row,col] * (1 - exp(- kappa[row] - kappa[col])) / (kappa[row] + kappa[col])
    }
  }
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSd_3F(age, Low_chol, delta_matrix)####### A_ind(age, exp(l_sigma), delta)  
    B_tT[age,] <- B_BSd_3F(age, delta_matrix)  ###### B_ind(age,delta)  
  }
  
  for(t in 1:n_years){
    
    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))    #max(0.000001, r_c + r_1 * exp(r_2))#exp(r_1 + r_2)  In case of the error specification of Xu et al. (2019)
    y_star[1,t] <- A_tT[1,1] + B_tT[1,] %*% x_ti + sqrt(F_ti[1,t]) * e_star[1,t]
    
    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / sqrt(F_ti[i-1,t])) %*% e_star[i-1,t] # * (v_ti[i-1,t] / F_ti[i-1,t])   ### - Same change from + to minus following Xu, Sherris and Ziveyi 2019SAJ               # - a_{t,i+1}
      
      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) #exp(r_1 + r_2 * (i-1)) %*% t(K_ti) ####(max(0.001, r_c + r_1 * sum(exp(r_2 * c(1:(i-1)))) / (i-1))) %*% t(K_ti)
      
      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i) 
      y_star[i,t] <- A_tT[i,1] + B_tT[i,] %*% x_ti + sqrt(F_ti[i,t]) * e_star[i,t]
      
    }
  }
  
  return(y_star)
}

nLL_BSd_3F_uKD_BS <- function(vdParameters, mu_bar_BS){
  n_factors <- 3
  # - Parameters
  delta <- vdParameters[1:6]####### - TO BE CHANGED
  kappa <- vdParameters[7:9] # - take logs to ensure positivity
  dg_l_Sigma_chol <- vdParameters[10:12] # - log diagonal elements of the lower cholesky dec.
  odg_Sigma_chol <- vdParameters[13:15] # - off diagonal elements of the lower cholesky dec.
  r_1 <- vdParameters[16]
  r_2 <- vdParameters[17]
  r_c <- vdParameters[18]
  
  n_ages <- nrow(mu_bar_BS)   # - Number of ages
  n_years <- ncol(mu_bar_BS)  # - Number of years
  
  v_ti <- mu_bar_BS
  F_ti <- mu_bar_BS
  
  delta_matrix <- low_trg_fill(delta)
  
  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance
  
  # - Initialize X and Sigma
  x_ti <- x0_s #init_X
  P_ti <- 1e-10 * diag(1, n_factors)
  
  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)
  
  # - Get Sigma (covariance matrix of the diffusion process)
  Sigma_diff <- Low_chol %*% t(Low_chol) # matrix(0, n_factors, n_factors)
  
  # - Get R (covariance of the state variable)
  for(row in 1:n_factors){
    for(col in 1:n_factors){
      R[row,col] <- Sigma_diff[row,col] * (1 - exp(- kappa[row] - kappa[col])) / (kappa[row] + kappa[col])
    }
  }
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSd_3F(age, Low_chol, delta_matrix)####### A_ind(age, exp(l_sigma), delta)  
    B_tT[age,] <- B_BSd_3F(age, delta_matrix)  ###### B_ind(age,delta)  
  }
  
  for(t in 1:n_years){
    
    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar_BS[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))  
    
    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t] 
      
      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) 
      
      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i 
      v_ti[i,t] <- mu_bar_BS[i,t] - A_tT[i] - B_tT[i,] %*% x_ti
    }
  }
  
  log_lik <- sum(log(F_ti) + (v_ti^2) / F_ti)
  return(log_lik)
}

CovEst_BS_BSd_3F <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar, n_BS=500, t_ex = 4){
  
  par_table_BSd <- matrix(NA, n_BS, length(c(delta, kappa, sigma_dg, Sigma_cov, r)))
  colnames(par_table_BSd) <- c("delta_11", 'delta_21', 'delta_22', 'delta_31', 'delta_32', 'delta_33', sprintf("kappa_%d", c(1:n_factors)), 'sigma_11', 'sigma_21', 'sigma_22', 'sigma_31', 'sigma_32', 'sigma_33', c("r1", "r2", "rc"))
  res_table <- matrix(NA, n_ages, n_years)
  
  # - 4) Parameter estimation
  ## - 4.1) Get filtered estimates
  Filtering <- KF_BSd_3F_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r)
  X_t_fil <- Filtering$X_t
  X_t_c_fil <- Filtering$X_t_c
  S_t_fil <- Filtering$S_t
  S_t_c_fil <- Filtering$S_t_c
  
  ## - 4.2) Get smoothed estimate of x0
  Smoothing <- RTS_sm_bas(X_t_fil, X_t_c_fil, S_t_fil, S_t_c_fil, kappa)
  X_t_sm <- Smoothing$X_t_sm[,1]
  S_t_sm <- Smoothing$S_t_sm[,1:3]
  
  # - 1) Get standardized residuals
  std_r <- std_res_BSd_3F(x0, delta, kappa, sigma_dg, Sigma_cov, r)
  ## - Fill the first t_ex residuals, which stay the same
  res_table[,1:t_ex] <- std_r[,1:t_ex]
  
  # - for starting values
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  
  # - Repeat bootstrapping samples
  for(i in 1:n_BS){
    
    # - 2) resample standardized residuals
    t_sample <- sample(c((t_ex+1):n_years), (n_years-t_ex), replace=TRUE)  # - starting from 5 to account for KF start up irregularity
    res_table[,(t_ex + 1):n_years] <- std_r[,t_sample]
    
    # - 4) Parameter estimation
    
    ## - 4.3) Sample x0 from the smoothing distribution
    x0_s <- mvrnorm(n = 1, mu=X_t_sm, Sigma = S_t_sm, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    
    # - 3) Get mu_bar_star (small correction wrt the initial code)
    mu_bar_BS <- y_star_BSd_3F(x0_s, delta, kappa, sigma_dg, Sigma_cov, r, res_table)
    
    ## - 4.4) Get bootstrapped parameter estimates
    par_table_BSd[i,] <- optim(c(delta, kappa, dg_l_Sigma_chol, odg_Sigma_chol, log(r)), nLL_BSd_3F_uKD_BS, mu_bar_BS=mu_bar_BS, gr = NULL, method = "Nelder-Mead", hessian = FALSE, control=list(trace=TRUE, maxit = 10000))$par # or coordinate ascent
    
    ## - Backtransformation
    par_table_BSd[i,10:15] <- parest2cov(par_table_BSd[i,10:12], par_table_BSd[i,13:15])
    par_table_BSd[i,16:18] <- exp(par_table_BSd[i,16:18])
  }
  
  # 5) Get st. err. parameter estimates
  cov_pe <- cov(par_table_BSd)
  serr_pe <- sqrt(diag(cov_pe_BSd))
  
  return(list(Cov = cov_pe, St.err = serr_pe))
}



# - AFNS model with independent factors

std_res_AFNSi <- function(x0, delta, kappa, sigma, r, mu_bar){
  
  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])
  
  n_factors <- 3

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years
  
  v_ti <- matrix(0, nrow(mu_bar), ncol(mu_bar))
  F_ti <- matrix(1, nrow(mu_bar), ncol(mu_bar))
  
  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance
  
  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFNSi(age, sigma, delta)  
    B_tT[age,] <- B_AFNS(age,delta)  
  }
  
  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  R <- diag(((sigma^2) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)
  
  for(t in 1:n_years){
    
    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))    #max(0.000001, r_c + r_1 * exp(r_2))#exp(r_1 + r_2)  In case of the error specification of Xu et al. (2019)
    
    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t] # * (v_ti[i-1,t] / F_ti[i-1,t])   ### - Same change from + to minus following Xu, Sherris and Ziveyi 2019SAJ               # - a_{t,i+1}
      
      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) #exp(r_1 + r_2 * (i-1)) %*% t(K_ti) ####(max(0.001, r_c + r_1 * sum(exp(r_2 * c(1:(i-1)))) / (i-1))) %*% t(K_ti)
      
      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i) 
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti
      
    }
  }
  
  e_hat <- v_ti / sqrt(F_ti)
  
  return(e_hat)
}

y_star_AFNSi <- function(x0, delta, kappa, sigma, r, e_star){
  
  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])
  
  n_factors <- 3

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years
  
  y_star <- matrix(0, nrow(mu_bar), ncol(mu_bar))
  F_ti <- matrix(1, nrow(mu_bar), ncol(mu_bar))
  
  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance
  
  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFNSi(age, sigma, delta)  
    B_tT[age,] <- B_AFNS(age,delta)  
  }
  
  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  R <- diag(((sigma^2) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)
  
  for(t in 1:n_years){
    
    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))    #max(0.000001, r_c + r_1 * exp(r_2))#exp(r_1 + r_2)  In case of the error specification of Xu et al. (2019)
    y_star[1,t] <- A_tT[1,1] + B_tT[1,] %*% x_ti + sqrt(F_ti[1,t]) * e_star[1,t]
    
    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / sqrt(F_ti[i-1,t])) %*% e_star[i-1,t] # * (v_ti[i-1,t] / F_ti[i-1,t])   ### - Same change from + to minus following Xu, Sherris and Ziveyi 2019SAJ               # - a_{t,i+1}
      
      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) #exp(r_1 + r_2 * (i-1)) %*% t(K_ti) ####(max(0.001, r_c + r_1 * sum(exp(r_2 * c(1:(i-1)))) / (i-1))) %*% t(K_ti)
      
      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i) 
      y_star[i,t] <- A_tT[i,1] + B_tT[i,] %*% x_ti + sqrt(F_ti[i,t]) * e_star[i,t]
      
    }
  }
  
  return(y_star)
}

nLL_AFNSi_uKD_BS <- function(vdParameters, mu_bar_BS){
  
  n_factors <- 3
  
  delta <- vdParameters[1]
  kappa <- vdParameters[2:(n_factors + 1)]
  l_sigma <- vdParameters[(n_factors + 2):(n_factors*2 + 1)] # - take logs to ensure positivity
  r_c <- vdParameters[n_factors*2 + 4]
  r_1 <- vdParameters[n_factors*2 + 2]
  r_2 <- vdParameters[n_factors*2 + 3]
  
  n_ages <- nrow(mu_bar_BS)   # - Number of ages
  n_years <- ncol(mu_bar_BS)  # - Number of years
  
  v_ti <- mu_bar_BS
  F_ti <- mu_bar_BS
  
  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance
  
  # - Initialize X and Sigma
  x_ti <- x0_s #init_X
  P_ti <- 1e-10 * diag(1, n_factors)
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFNSi(age, exp(l_sigma), delta)  
    B_tT[age,] <- B_AFNS(age,delta)  
  }
  
  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  R <- diag(((exp(l_sigma*2)) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)
  
  for(t in 1:n_years){
    
    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar_BS[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))   
    
    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t] 
      
      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) 
      
      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i 
      v_ti[i,t] <- mu_bar_BS[i,t] - A_tT[i] - B_tT[i,] %*% x_ti
      
    }
  }
  
  log_lik <- sum(log(F_ti) + (v_ti^2) / F_ti)
  return(log_lik)
}


CovEst_BS_AFNSi <- function(x0, delta, kappa, sigma, r, mu_bar, n_BS=500, t_ex = 4){
  par_table_AFNSi <- matrix(NA, n_BS, length(c(delta, kappa, sigma, r)))
  colnames(par_table_AFNSi) <- c("delta", sprintf("kappa_%s", c("L", 'S', 'C')), sprintf("sigma_%s", c("L", 'S', 'C')), c("r1", "r2", "rc"))
  
  res_table <- matrix(NA, n_ages, n_years)
  
  # - 4) Parameter estimation
  ## - 4.1) Get filtered estimates
  Filtering <- KF_AFNSi_uKD(x0, delta, kappa, sigma, r)
  X_t_fil <- Filtering$X_t
  X_t_c_fil <- Filtering$X_t_c
  S_t_fil <- Filtering$S_t
  S_t_c_fil <- Filtering$S_t_c
  
  ## - 4.2) Get smoothed estimate of x0
  Smoothing <- RTS_sm_bas(X_t_fil, X_t_c_fil, S_t_fil, S_t_c_fil, kappa)
  X_t_sm <- Smoothing$X_t_sm[,1]
  S_t_sm <- Smoothing$S_t_sm[,1:n_factors]
  
  # - 1) Get standardized residuals
  std_r <- std_res_AFNSi(x0, delta, kappa, sigma, r)
  ## - Fill the first t_ex residuals, which stay the same
  res_table[,1:t_ex] <- std_r[,1:t_ex]
  
  # - Repeat bootstrapping samples
  for(i in 1:n_BS){
    
    # - 2) resample standardized residuals
    t_sample <- sample(c((t_ex+1):n_years), (n_years-t_ex), replace=TRUE)  # - starting from 5 to account for KF start up irregularity
    res_table[,(t_ex + 1):n_years] <- std_r[,t_sample]
    
    ## - 4.3) Sample x0 from the smoothing distribution
    x0_s <- mvrnorm(n = 1, mu=X_t_sm, Sigma = S_t_sm, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    
    # - 3) Get mu_bar_star
    mu_bar_BS <- y_star_AFNSi(x0_s, delta, kappa, sigma, r, res_table)
    
    ## - 4.4) Get bootstrapped parameter estimates
    par_table_AFNSi[i,] <- optim(c(delta, kappa, log(sigma), log(r)), nLL_AFNSi_uKD_BS, mu_bar_BS=mu_bar_BS, gr = NULL, method = "BFGS", hessian = FALSE, control=list(trace=TRUE, maxit = 10000))$par # or coordinate ascent
    ## - Backtransformation
    par_table_AFNSi[i,5:10] <- exp(par_table_AFNSi[i,5:10])
  }
  
  # 5) Get st. err. parameter estimates (unconstrained parameters)
  cov_pe <- cov(par_table_AFNSi)
  serr_pe <- sqrt(diag(cov_pe_AFNSi))
  
  return(list(Cov = cov_pe, St.err = serr_pe))
}


# - AFNS model with dependent factors

std_res_AFNSd <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar){
  
  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])
  
  n_factors <- 3
  
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years
  
  v_ti <- matrix(0, nrow(mu_bar), ncol(mu_bar))
  F_ti <- matrix(1, nrow(mu_bar), ncol(mu_bar))
  
  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance
  
  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)
  
  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  
  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)
  
  # - Get Sigma (covariance matrix of the diffusion process)
  Sigma_diff <- Low_chol %*% t(Low_chol) # matrix(0, n_factors, n_factors)
  
  # - Get R (covariance of the state variable)
  for(row in 1:n_factors){
    for(col in 1:n_factors){
      R[row,col] <- Sigma_diff[row,col] * (1 - exp(- kappa[row] - kappa[col])) / (kappa[row] + kappa[col])
    }
  }
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFNSg(age, Low_chol, delta)####### A_ind(age, exp(l_sigma), delta)  
    B_tT[age,] <- B_AFNS(age, delta)  ###### B_ind(age,delta)  
  }
  
  for(t in 1:n_years){
    
    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))   
    
    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t] 
      
      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) 
      
      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i) 
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti
      
    }
  }
  
  e_hat <- v_ti / sqrt(F_ti)
  return(e_hat)
}

y_star_AFNSd <- function(x0, delta, kappa, dg_l_Sigma_chol, odg_Sigma_chol, r, e_star){
  
  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])
  
  n_factors <- 3
  
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years
  
  y_star <- matrix(0, nrow(mu_bar), ncol(mu_bar))
  F_ti <- matrix(1, nrow(mu_bar), ncol(mu_bar))
  
  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance
  
  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)
  
  Phi <- diag(exp(-kappa), n_factors)
  
  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)
  
  # - Get Sigma (covariance matrix of the diffusion process)
  Sigma_diff <- Low_chol %*% t(Low_chol) # matrix(0, n_factors, n_factors)
  
  # - Get R (covariance of the state variable)
  for(row in 1:n_factors){
    for(col in 1:n_factors){
      R[row,col] <- Sigma_diff[row,col] * (1 - exp(- kappa[row] - kappa[col])) / (kappa[row] + kappa[col])
    }
  }
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFNSg(age, Low_chol, delta)####### A_ind(age, exp(l_sigma), delta)  
    B_tT[age,] <- B_AFNS(age, delta)  ###### B_ind(age,delta)  
  }
  
  for(t in 1:n_years){
    
    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))    #max(0.000001, r_c + r_1 * exp(r_2))#exp(r_1 + r_2)  In case of the error specification of Xu et al. (2019)
    y_star[1,t] <- A_tT[1,1] + B_tT[1,] %*% x_ti + sqrt(F_ti[1,t]) * e_star[1,t]
    
    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / sqrt(F_ti[i-1,t])) %*% e_star[i-1,t] # * (v_ti[i-1,t] / F_ti[i-1,t])   ### - Same change from + to minus following Xu, Sherris and Ziveyi 2019SAJ               # - a_{t,i+1}
      
      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) #exp(r_1 + r_2 * (i-1)) %*% t(K_ti) ####(max(0.001, r_c + r_1 * sum(exp(r_2 * c(1:(i-1)))) / (i-1))) %*% t(K_ti)
      
      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i) 
      y_star[i,t] <- A_tT[i,1] + B_tT[i,] %*% x_ti + sqrt(F_ti[i,t]) * e_star[i,t]
      
    }
  }
  
  return(y_star)
}

nLL_AFNSd_uKD_BS <- function(vdParameters, mu_bar_BS){
  n_factors <- 3
  # - Parameters
  delta <- vdParameters[1]
  kappa <- vdParameters[2:(n_factors + 1)] # - take logs to ensure positivity    n_factors * (3 + n_factors) / 2 + 1
  dg_l_Sigma_chol <- vdParameters[(n_factors + 2):(2 * n_factors + 1)] # - log diagonal elements of the lower cholesky dec.
  odg_Sigma_chol <- vdParameters[(2 * n_factors + 2):(n_factors * (3 + n_factors) / 2 + 1)] # - off diagonal elements of the lower cholesky dec.
  r_c <- vdParameters[n_factors * (3 + n_factors) / 2 + 4]
  r_1 <- vdParameters[n_factors * (3 + n_factors) / 2 + 2]
  r_2 <- vdParameters[n_factors * (3 + n_factors) / 2 + 3]
  
  n_ages <- nrow(mu_bar_BS)   # - Number of ages
  n_years <- ncol(mu_bar_BS)  # - Number of years
  
  v_ti <- mu_bar_BS
  F_ti <- mu_bar_BS
  
  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance
  
  # - Initialize X and Sigma
  x_ti <- x0_s #init_X
  P_ti <- 1e-10 * diag(1, n_factors)
  
  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  
  # - Build diffusion process
  ## - Build lower cholesky factor
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)
  
  # - Get Sigma (covariance matrix of the diffusion process)
  Sigma_diff <- Low_chol %*% t(Low_chol) # matrix(0, n_factors, n_factors)
  
  # - Get R (covariance of the state variable)
  for(row in 1:n_factors){
    for(col in 1:n_factors){
      R[row,col] <- Sigma_diff[row,col] * (1 - exp(- kappa[row] - kappa[col])) / (kappa[row] + kappa[col])
    }
  }
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFNSg(age, Low_chol, delta)####### A_ind(age, exp(l_sigma), delta)  
    B_tT[age,] <- B_AFNS(age, delta)  ###### B_ind(age,delta)  
  }
  
  for(t in 1:n_years){
    
    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar_BS[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))    
    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t] 
      
      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) 

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i 
      v_ti[i,t] <- mu_bar_BS[i,t] - A_tT[i] - B_tT[i,] %*% x_ti
      
    }
  }
  
  log_lik <- sum(log(F_ti) + (v_ti^2) / F_ti)
  return(log_lik)
}


CovEst_BS_AFNSd <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar, n_BS=500, t_ex = 4){
  
  par_table_AFNSd <- matrix(NA, n_BS, length(c(delta, kappa, sigma_dg, Sigma_cov, r)))
  colnames(par_table_AFNSd) <- c("delta", sprintf("kappa_%d", c(1:n_factors)), 'sigma_L', 'sigma_LS', 'sigma_S', 'sigma_LC', 'sigma_SC', 'sigma_C', c("r1", "r2", "rc"))
  res_table <- matrix(NA, n_ages, n_years)
  
  # - 4) Parameter estimation
  ## - 4.1) Get filtered estimates
  Filtering <- KF_AFNSd_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r)
  X_t_fil <- Filtering$X_t
  X_t_c_fil <- Filtering$X_t_c
  S_t_fil <- Filtering$S_t
  S_t_c_fil <- Filtering$S_t_c
  
  ## - 4.2) Get smoothed estimate of x0
  Smoothing <- RTS_sm_bas(X_t_fil, X_t_c_fil, S_t_fil, S_t_c_fil, kappa)
  X_t_sm <- Smoothing$X_t_sm[,1]
  S_t_sm <- Smoothing$S_t_sm[,1:3]
  
  # - 1) Get standardized residuals
  std_r <- std_res_AFNSd(x0, delta, kappa, sigma_dg, Sigma_cov, r)
  ## - Fill the first t_ex residuals, which stay the same
  res_table[,1:t_ex] <- std_r[,1:t_ex]
  
  # - for starting values
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  
  # - Repeat bootstrapping samples
  for(i in 1:n_BS){
    
    # - 2) resample standardized residuals
    t_sample <- sample(c((t_ex+1):n_years), (n_years-t_ex), replace=TRUE)  # - starting from 5 to account for KF start up irregularity
    res_table[,(t_ex + 1):n_years] <- std_r[,t_sample]
    
    # - 4) Parameter estimation
    
    ## - 4.3) Sample x0 from the smoothing distribution
    x0_s <- mvrnorm(n = 1, mu=X_t_sm, Sigma = S_t_sm, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    
    # - 3) Get mu_bar_star (small correction wrt the initial code)
    mu_bar_BS <- y_star_AFNSd(x0_s, delta, kappa, sigma_dg, Sigma_cov, r, res_table)
    
    ## - 4.4) Get bootstrapped parameter estimates
    par_table_AFNSd[i,] <- optim(c(delta, kappa, dg_l_Sigma_chol, odg_Sigma_chol, log(r)), nLL_AFNSd_uKD_BS, mu_bar_BS=mu_bar_BS, gr = NULL, method = "Nelder-Mead", hessian = FALSE, control=list(trace=TRUE, maxit = 10000))$par # or coordinate ascent
    
    ## - Backtransformation
    par_table_AFNSd[i,5:10] <- parest2cov(par_table_AFNSd[i,5:7], par_table_AFNSd[i,8:10])
    par_table_AFNSd[i,11:13] <- exp(par_table_BSd[i,11:13])
  }
  
  # 5) Get st. err. parameter estimates
  cov_pe <- cov(par_table_AFNSd)
  serr_pe <- sqrt(diag(cov_pe_AFNSd))
  
  return(list(Cov = cov_pe, St.err = serr_pe))
}


# - CIR

std_res_CIR <- function(x0, delta, kappa, sigma, theta_Q, theta_P, r, mu_bar){
  
  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])
  
  n_factors <- length(kappa)
  
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years
  
  v_ti <- matrix(0, nrow(mu_bar), ncol(mu_bar))
  F_ti <- matrix(1, nrow(mu_bar), ncol(mu_bar))
  
  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors * (n_years + 1)) # - Factor covariance
  
  # - Initialize X and Sigma
  #x_ti <- x_0 #init_X
  x_ti <- x0
  P_ti <- 1e-10 * diag(1, n_factors)
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_CIR(age, theta_Q, sigma, delta)  
    B_tT[age,] <- B_CIR(age, sigma, delta)  
  }
  
  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  
  for(t in 1:n_years){
    
    R[,(t * n_factors + 1):((t+1) * n_factors)] <- diag((sigma^2) * ((1 - exp(-kappa)) / kappa) * (0.5 * theta_P * (1 - exp(-kappa)) + exp(-kappa) * x_ti[1:n_factors]),n_factors)
    
    # - First observation
    x_ti <- Phi %*% x_ti + theta_P * (1 - exp(-kappa))
    x_ti <- l_bound(x_ti)
    P_ti <- Phi %*% P_ti %*% t(Phi) + R[,(t * n_factors + 1):((t+1) * n_factors)]     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))    #max(0.000001, r_c + r_1 * exp(r_2))#exp(r_1 + r_2)  In case of the error specification of Xu et al. (2019)
    
    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t]
      x_ti <- l_bound(x_ti)
      
      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) #exp(r_1 + r_2 * (i-1)) %*% t(K_ti) ####(max(0.001, r_c + r_1 * sum(exp(r_2 * c(1:(i-1)))) / (i-1))) %*% t(K_ti)
      
      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i) 
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti
      
    }
  }
  
  e_hat <- v_ti / sqrt(F_ti)
  return(e_hat)
}

y_star_CIR <- function(x0, delta, kappa, sigma, theta_Q, theta_P, r, e_star){
  
  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])
  
  n_factors <- length(kappa)
  
  n_ages <- nrow(e_star)   # - Number of ages
  n_years <- ncol(e_star)  # - Number of years
  
  y_star <- matrix(0, nrow(e_star), ncol(e_star))
  F_ti <- matrix(1, nrow(e_star), ncol(e_star))
  
  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors * (n_years + 1)) # - Factor covariance
  
  # - Initialize X and Sigma
  #x_ti <- x_0 #init_X
  x_ti <- x0
  P_ti <- 1e-10 * diag(1, n_factors)
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_CIR(age, theta_Q, sigma, delta)  
    B_tT[age,] <- B_CIR(age, sigma, delta)  
  }
  
  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  
  for(t in 1:n_years){
    
    R[,(t * n_factors + 1):((t+1) * n_factors)] <- diag((sigma^2) * ((1 - exp(-kappa)) / kappa) * (0.5 * theta_P * (1 - exp(-kappa)) + exp(-kappa) * x_ti[1:n_factors]),n_factors)
    
    # - First observation
    x_ti <- Phi %*% x_ti + theta_P * (1 - exp(-kappa))
    x_ti <- l_bound(x_ti)
    P_ti <- Phi %*% P_ti %*% t(Phi) + R[,(t * n_factors + 1):((t+1) * n_factors)]     # - P_{1,t}
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))    
    y_star[1,t] <- A_tT[1,1] + B_tT[1,] %*% x_ti + sqrt(F_ti[1,t]) * e_star[1,t]
    
    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / sqrt(F_ti[i-1,t])) %*% e_star[i-1,t]
      x_ti <- l_bound(x_ti)
      
      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) 
      
      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i) 
      y_star[i,t] <- A_tT[i,1] + B_tT[i,] %*% x_ti + sqrt(F_ti[i,t]) * e_star[i,t]
      
    }
  }
  
  return(y_star)
}

nLL_CIR_uKD_BS <- function(vdParameters, mu_bar_BS){
  n_factors <- length(vdParameters) / 6
  # - Parameters
  delta <- vdParameters[1:n_factors]
  l_kappa <- vdParameters[(n_factors+1):(n_factors*2)]
  l_sigma <- vdParameters[(n_factors*2+1):(n_factors*3)] # - take logs to ensure positivity
  l_theta_Q <- vdParameters[(n_factors*3+1):(n_factors*4)]
  l_theta_P <- vdParameters[(n_factors*4+1):(n_factors*5)]
  r_1 <- vdParameters[n_factors*5 + 1]
  r_2 <- vdParameters[n_factors*5 + 2]
  r_c <- vdParameters[n_factors*5 + 3]
  
  n_ages <- nrow(mu_bar_BS)   # - Number of ages
  n_years <- ncol(mu_bar_BS)  # - Number of years
  
  v_ti <- mu_bar_BS
  F_ti <- mu_bar_BS
  
  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors * (n_years + 1)) # - Factor covariance
  
  # - Initialize X and Sigma
  #x_ti <- x_0 #init_X
  x_ti <- x0_s
  P_ti <- 1e-10 * diag(1, n_factors)
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_CIR(age, exp(l_theta_Q), exp(l_sigma), delta)  
    B_tT[age,] <- B_CIR(age, exp(l_sigma), delta)  
  }
  
  Phi <- diag(exp(-exp(l_kappa)), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  
  for(t in 1:n_years){
    
    R[,(t * n_factors + 1):((t+1) * n_factors)] <- diag(exp(2 * l_sigma) * ((1 - exp(-exp(l_kappa))) / exp(l_kappa)) * (0.5 * exp(l_theta_P) * (1 - exp(-exp(l_kappa))) + exp(-exp(l_kappa)) * x_ti[1:n_factors]),n_factors)
    
    # - First observation
    x_ti <- Phi %*% x_ti + exp(l_theta_P) * (1 - exp(-exp(l_kappa)))
    x_ti <- l_bound(x_ti)
    P_ti <- Phi %*% P_ti %*% t(Phi) + R[,(t * n_factors + 1):((t+1) * n_factors)]     # - P_{1,t}
    v_ti[1,t] <- mu_bar_BS[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))  
    
    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t]
      x_ti <- l_bound(x_ti)
      
      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) 
      
      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i 
      v_ti[i,t] <- mu_bar_BS[i,t] - A_tT[i] - B_tT[i,] %*% x_ti
      
    }
  }
  
  log_lik <- sum(log(F_ti) + (v_ti^2) / F_ti)
  return(log_lik)
}


CovEst_BS_CIR <- function(x0, delta, kappa, sigma, theta_Q, theta_P, r, mu_bar, n_BS=500, t_ex = 4){
  par_table_CIR <- matrix(NA, n_BS, length(c(delta, kappa, sigma, theta_Q, theta_P, r)))
  colnames(par_table_CIR) <- c(sprintf("delta_%d", c(1:n_factors)), sprintf("kappa_%d", c(1:n_factors)), sprintf("sigma_%d", c(1:n_factors)), sprintf("theta_Q_%d", c(1:n_factors)), sprintf("theta_P_%d", c(1:n_factors)), c("r1", "r2", "rc"))
  
  res_table <- matrix(NA, n_ages, n_years)
  
  # - 4) Parameter estimation
  ## - 4.1) Get filtered estimates
  Filtering <- KF_CIR_uKD(x0, delta, kappa, sigma, theta_Q, theta_P, r)
  X_t_fil <- Filtering$X_t
  X_t_c_fil <- Filtering$X_t_c
  S_t_fil <- Filtering$S_t
  S_t_c_fil <- Filtering$S_t_c
  
  ## - 4.2) Get smoothed estimate of x0
  Smoothing <- RTS_sm_bas(X_t_fil, X_t_c_fil, S_t_fil, S_t_c_fil, kappa)
  X_t_sm <- Smoothing$X_t_sm[,1]
  S_t_sm <- Smoothing$S_t_sm[,1:n_factors]
  
  # - 1) Get standardized residuals
  std_r <- std_res_CIR(x0, delta, kappa, sigma, theta_Q, theta_P, r)
  ## - Fill the first t_ex residuals, which stay the same
  res_table[,1:t_ex] <- std_r[,1:t_ex]
  
  # - Repeat bootstrapping samples
  for(i in 1:n_BS){
    
    # - 2) resample standardized residuals
    t_sample <- sample(c((t_ex+1):n_years), (n_years-t_ex), replace=TRUE)  # - starting from 5 to account for KF start up irregularity
    res_table[,(t_ex + 1):n_years] <- std_r[,t_sample]
    
    ## - 4.3) Sample x0 from the smoothing distribution - adjusted to ensure a lower bound of zero
    x0_s <- l_bound(mvrnorm(n = 1, mu=X_t_sm, Sigma = S_t_sm, tol = 1e-6, empirical = FALSE, EISPACK = FALSE))
    
    # - 3) Get mu_bar_star
    mu_bar_BS <- y_star_CIR(x0_s, delta, kappa, sigma, theta_Q, theta_P, r, res_table)
    
    ## - 4.4) Get bootstrapped parameter estimates
    par_table_CIR[i,] <- optim(c(delta, kappa, log(sigma), log(theta_Q), log(theta_P), log(r)), nLL_CIR_uKD_BS, mu_bar_BS=mu_bar_BS, gr = NULL, method = "BFGS", hessian = FALSE, control=list(trace=TRUE, maxit = 10000))$par # or coordinate ascent
    ## - Backtransformation
    par_table_CIR[i,((n_factors + 1):(n_factors*6))] <- exp(par_table_CIR[i,((n_factors + 1):(n_factors*6))])
  }
  
  # 5) Get st. err. parameter estimates (unconstrained parameters)
  cov_pe <- cov(par_table_CIR)
  serr_pe <- sqrt(diag(cov_pe_CIR))
  
  return(list(Cov = cov_pe, St.err = serr_pe))
}



