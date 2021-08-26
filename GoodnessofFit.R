source('FilterSmoother.R')

#==========================================
# - Fitted rates (mu_bar_hat)
#==========================================

# - Blackburn-Sherris
## - Independent
mu_bar_hat_BSi <- function(x0, delta, kappa, sigma, r, mu_bar){
  n_factors <- length(kappa)

  KF_outcome <- KF_BSi_uKD(x0, delta, kappa, sigma, r, mu_bar)
  X_t <- KF_outcome$X_t
  
  mu_bar_hat <- matrix(NA, n_ages, n_years)
  rownames(mu_bar_hat) <- AgeRange
  colnames(mu_bar_hat) <- CohortRange
    
  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSi(age, sigma, delta)  
    B_tT[age,] <- B_BSi(age,delta)  
  }
  
  for(t in 1:n_years){
    mu_bar_hat[,t] <- A_tT + B_tT %*% X_t[,t+1]
  }
  return(mu_bar_hat)
}

## - Dependent
mu_bar_hat_BSd_2F <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar){
  n_factors <- 2

  delta_matrix <- low_trg_fill(delta) # - create lower diagonal matrix
  
  KF_outcome <- KF_BSd_2F_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)
  X_t <- KF_outcome$X_t
  
  mu_bar_hat <- matrix(NA, n_ages, n_years)
  rownames(mu_bar_hat) <- AgeRange
  colnames(mu_bar_hat) <- CohortRange
  
  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  
  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSd_2F(age, Low_chol, delta_matrix)####### A_ind(age, exp(l_sigma), delta)  
    B_tT[age,] <- B_BSd_2F(age, delta_matrix)  ###### B_ind(age,delta)  
  }
  
  for(t in 1:n_years){
    mu_bar_hat[,t] <- A_tT + B_tT %*% X_t[,t+1]
  }
  return(mu_bar_hat)
}

mu_bar_hat_BSd_3F <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar){
  n_factors <- 3

  delta_matrix <- low_trg_fill(delta) # - create lower diagonal matrix
  
  KF_outcome <- KF_BSd_3F_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)
  X_t <- KF_outcome$X_t
  
  mu_bar_hat <- matrix(NA, n_ages, n_years)
  rownames(mu_bar_hat) <- AgeRange
  colnames(mu_bar_hat) <- CohortRange
  
  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  
  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSd_3F(age, Low_chol, delta_matrix)####### A_ind(age, exp(l_sigma), delta)  
    B_tT[age,] <- B_BSd_3F(age, delta_matrix)  ###### B_ind(age,delta)  
  }
  
  for(t in 1:n_years){
    mu_bar_hat[,t] <- A_tT + B_tT %*% X_t[,t+1]
  }
  return(mu_bar_hat)
}

# - AFNS
## - Independent
mu_bar_hat_AFNSi <- function(x0, delta, kappa, sigma, r, mu_bar){
  n_factors <- 3
  
  KF_outcome <- KF_AFNSi_uKD(x0, delta, kappa, sigma, r, mu_bar)
  X_t <- KF_outcome$X_t
  
  mu_bar_hat <- matrix(NA, n_ages, n_years)
  rownames(mu_bar_hat) <- AgeRange
  colnames(mu_bar_hat) <- CohortRange
  
  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFNSi(age, sigma, delta)  
    B_tT[age,] <- B_AFNS(age,delta)  
  }
  
  for(t in 1:n_years){
    mu_bar_hat[,t] <- A_tT + B_tT %*% X_t[,t+1]
  }
  return(mu_bar_hat)
}

## - Dependent
mu_bar_hat_AFNSd <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar){
  n_factors <- 3

  KF_outcome <- KF_AFNSd_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)
  X_t <- KF_outcome$X_t
  
  mu_bar_hat <- matrix(NA, n_ages, n_years)
  rownames(mu_bar_hat) <- AgeRange
  colnames(mu_bar_hat) <- CohortRange
  
  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  
  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFNSg(age, Low_chol, delta)####### A_ind(age, exp(l_sigma), delta)  
    B_tT[age,] <- B_AFNS(age, delta)  ###### B_ind(age,delta)  
  }
  
  for(t in 1:n_years){
    mu_bar_hat[,t] <- A_tT + B_tT %*% X_t[,t+1]
  }
  return(mu_bar_hat)
}

## - CIR in the basic version (since the two log-likelihoods yield very different results)
mu_bar_hat_CIR <- function(x0, delta, kappa, sigma, theta_Q, theta_P, r, mu_bar){
  n_factors <- length(kappa)

  KF_outcome <- KF_CIR_uKD(x0, delta, kappa, sigma, theta_Q, theta_P, r, mu_bar)
  X_t <- KF_outcome$X_t
  
  mu_bar_hat <- matrix(NA, n_ages, n_years)
  rownames(mu_bar_hat) <- AgeRange
  colnames(mu_bar_hat) <- CohortRange
  
  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_CIR(age, theta_Q, sigma, delta)  
    B_tT[age,] <- B_CIR(age, sigma, delta)  
  }
  
  for(t in 1:n_years){
    mu_bar_hat[,t] <- A_tT + B_tT %*% X_t[,t+1]
  }
  return(mu_bar_hat)
}


#==========================================================
# - Model Selection (AIC and BIC) - Build cross-validation
#==========================================================

AIC_BIC <- function(log_likelihood, n_par, n_obs){
  
  AIC <- - 2 * log_likelihood + 2 * n_par
  BIC <- - 2 * log_likelihood + n_par * log(n_obs)
  
  return(list(AIC=AIC, BIC=BIC))
}

#==================================================
# - Root Mean Squared Error (RMSE)
#==================================================

# - First reproduce mu_bar_hat for all models

RMSE <- function(observed, estimated){
  RMSE <- sqrt(mean((observed - estimated)^2))
  return(RMSE)  
}


#=======================================================
# - Mean Absolute Percentage Error (MAPE)
#=======================================================

MAPE_row <- function(observed, estimated){
  MAPE <- rowMeans(abs((observed - estimated)/observed))
  return(MAPE)  
}


#=======================================================
# - Fitted rates (useful for residuals plot and RMSE)
#=======================================================

# - mu_bar hat

# - Residuals
residuals_f <- function(observed, estimated){
  residuals <- observed - estimated
  colnames(residuals) <- colnames(mu_bar)
  rownames(residuals) <- rownames(mu_bar)
  return(residuals)
}

## - 0/1 residuals (0 if negative, 1 if positive)
residuals_01 <- function(observed, estimated){
  residuals <- ifelse(observed - estimated > 0, 1, 0)
  colnames(residuals) <- colnames(observed)
  rownames(residuals) <- rownames(observed)
  return(residuals)
}


# - Standardized residuals
residuals_std_BSi <- function(observed, x0, delta, kappa, sigma, r){
  n_factors <- length(kappa)
  
  l_r1 <- log(r[1])
  l_r2 <- log(r[2])
  l_rc <- log(r[3])
  
  for(age in 1:n_ages){    # - scroll over the ages
    B_tT[age,] <- B_BSi(age,delta)  
  }
  
  residuals <- observed - mu_bar_hat_BSi(x0, delta, kappa, sigma, r, observed)
  H <- meas_err_BS(l_r1, l_r2, l_rc)
  Var_y <- H + B_tT %*% diag(((sigma^2) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors) %*% t(B_tT)
  #Chol_dec <- chol(Var_y)
  premultiplier <- solve(sqrtm(Var_y))#solve(t(Chol_dec))
  
  for(t in 1:n_years){
    residuals[,t] <- premultiplier %*% residuals[,t]
  }
  
  colnames(residuals) <- colnames(observed)
  rownames(residuals) <- rownames(observed)
  return(residuals)
}

residuals_std_BSd_3F <- function(observed, x0, delta, kappa, sigma_dg, Sigma_cov, r){
  n_factors <- 3
  
  l_r1 <- log(r[1])
  l_r2 <- log(r[2])
  l_rc <- log(r[3])
  
  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)
  
  delta_matrix <- low_trg_fill(delta) # - create lower diagonal matrix
  
  R_matrix <- matrix(0, n_factors, n_factors)
  # - Get Sigma (covariance matrix of the diffusion process)
  Sigma_diff <- Low_chol %*% t(Low_chol) # matrix(0, n_factors, n_factors)
  
  # - Get R (covariance of the state variable)
  for(row in 1:n_factors){
    for(col in 1:n_factors){
      R_matrix[row,col] <- Sigma_diff[row,col] * (1 - exp(- kappa[row] - kappa[col])) / (kappa[row] + kappa[col])
    }
  }
  
  for(age in 1:n_ages){    # - scroll over the ages
    B_tT[age,] <- B_BSd_3F(age, delta_matrix) 
  }
  
  residuals <- observed - mu_bar_hat_BSd_3F(x0, delta, kappa, sigma_dg, Sigma_cov, r, observed)
  H <- meas_err_BS(l_r1, l_r2, l_rc)
  
  Var_y <- H + B_tT %*% R_matrix %*% t(B_tT)
  premultiplier <- solve(sqrtm(Var_y))
  
  for(t in 1:n_years){
    residuals[,t] <- premultiplier %*% residuals[,t]
  }
  
  colnames(residuals) <- colnames(observed)
  rownames(residuals) <- rownames(observed)
  return(residuals)
}

residuals_std_AFNSi <- function(observed, x0, delta, kappa, sigma, r){
  n_factors <- 3
  
  l_r1 <- log(r[1])
  l_r2 <- log(r[2])
  l_rc <- log(r[3])
  
  for(age in 1:n_ages){    # - scroll over the ages
    B_tT[age,] <- B_AFNS(age,delta)  
  }
  
  residuals <- observed - mu_bar_hat_AFNSi(x0, delta, kappa, sigma, r)
  H <- meas_err_BS(l_r1, l_r2, l_rc)

  Var_y <- H + B_tT %*% diag(((sigma^2) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors) %*% t(B_tT)

  premultiplier <- solve(sqrtm(Var_y))
  
  for(t in 1:n_years){
    residuals[,t] <- premultiplier %*% residuals[,t]
  }
  
  colnames(residuals) <- colnames(observed)
  rownames(residuals) <- rownames(observed)
  return(residuals)
}

residuals_std_AFNSd <- function(observed, delta, kappa, sigma_dg, Sigma_cov, r){
  n_factors <- 3
  
  l_r1 <- log(r[1])
  l_r2 <- log(r[2])
  l_rc <- log(r[3])
  
  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)
  delta_matrix <- low_trg_fill(delta) # - create lower diagonal matrix
  
  R_matrix <- matrix(0, n_factors, n_factors)
  # - Get Sigma (covariance matrix of the diffusion process)
  Sigma_diff <- Low_chol %*% t(Low_chol) # matrix(0, n_factors, n_factors)
  
  # - Get R (covariance of the state variable)
  for(row in 1:n_factors){
    for(col in 1:n_factors){
      R_matrix[row,col] <- Sigma_diff[row,col] * (1 - exp(- kappa[row] - kappa[col])) / (kappa[row] + kappa[col])
    }
  }
  
  for(age in 1:n_ages){    # - scroll over the ages
    B_tT[age,] <- B_AFNS(age,delta)  
  }
  
  residuals <- observed - mu_bar_hat_AFNSd(x0, delta, kappa, sigma_dg, Sigma_cov, r, observed)
  H <- meas_err_BS(l_r1, l_r2, l_rc)
  Var_y <- H + B_tT %*% R_matrix %*% t(B_tT)

  premultiplier <- solve(sqrtm(Var_y))
  
  for(t in 1:n_years){
    residuals[,t] <- premultiplier %*% residuals[,t]
  }
  
  colnames(residuals) <- colnames(observed)
  rownames(residuals) <- rownames(observed)
  return(residuals)
}

residuals_std_CIR <- function(observed, x0, delta, kappa, sigma, theta_Q, theta_P, r){
  n_factors <- length(kappa)

  l_r1 <- log(r[1])
  l_r2 <- log(r[2])
  l_rc <- log(r[3])
  
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance
  
  for(age in 1:n_ages){    # - scroll over the ages
    B_tT[age,] <- B_CIR(age, sigma, delta)  
  }
  
  residuals <- observed - mu_bar_hat_CIR(x0, delta, kappa, sigma, theta_Q, theta_P, r, observed)
  H <- meas_err_BS(l_r1, l_r2, l_rc)
  
  X_t <- KF_CIR_uKD(x0, delta, kappa, sigma, theta_Q, theta_P, r)$X_t
  P_ti <- 1e-10 * diag(1, n_factors)
  
  Var_y <- matrix(0, n_years, n_years)
  
  
  for(t in 1:n_years){
    R <- diag((sigma^2) * ((1 - exp(-kappa)) / kappa) * (0.5 * theta_P * (1 - exp(-kappa)) + exp(-kappa) * X_t[,t+1]), n_factors)
    Var_y <- H + (B_tT %*% R %*% t(B_tT)) 
    premultiplier <- solve(sqrtm(Var_y))
    
    residuals[,t] <- premultiplier %*% residuals[,t]
  }
  
  colnames(residuals) <- colnames(observed)
  rownames(residuals) <- rownames(observed)
  return(residuals)
}


# - Poisson standardized residuals

poi_r_std <- function(mu_bar_hat, deaths, exposures){
  mu_hat <- mu_bar_hat
  for(i in 2:nrow(mu_hat)){
    mu_hat[i,] <- i * mu_bar_hat[i,] - (i-1) * mu_bar_hat[i-1,]
  }
  D_hat <- mu_hat * exposures
  
  poi_res_std <- (deaths - D_hat) / sqrt(D_hat)
  
  return(poi_res_std)
}


avg2rates <- function(mu_bar){
  mu <- mu_bar
  for(i in 2:nrow(mu_bar)){
    mu[i,] <- i * mu_bar[i,] - (i-1) * mu_bar[i-1,]
  }
  return(mu)
}
















