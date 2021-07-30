
#====================== - Projected survival curves - ================

# - Blackburn - Sherris independent factor model

S_t_BSi_proj <- function(x0, delta, kappa, sigma, r, proj_years){
  
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)
  
  B_tT <- matrix(NA, n_ages, n_factors)
  X_t_last <- KF_BSi_uKD(x0, delta, kappa, sigma, r)$X_t[,n_years+1]

  E_X_t1 <- exp(-kappa * proj_years) * X_t_last
  
  S_prj <- matrix(NA, n_ages, 1)
  mu_hat_prj <- matrix(NA, n_ages, 1)
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSi(age, sigma, delta)  
    B_tT[age,] <- B_BSi(age, delta)  
    S_prj[age, 1] <- exp(-age * (A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1)) 
    mu_hat_prj[age,1] <- A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1
  }
  return(S_prj)
}

# - Blackburn - Sherris dependent factor model
S_t_BSd_2F_proj <- function(x0, delta, kappa, sigma, r, proj_years){
  
  n_factors <- 2
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)
  
  B_tT <- matrix(NA, n_ages, n_factors)
  X_t_last <- KF_BSd_2F_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r)$X_t[,n_years+1]
  
  E_X_t1 <- exp(-kappa * proj_years) * X_t_last
  
  S_prj <- matrix(NA, n_ages, 1)
  mu_hat_prj <- matrix(NA, n_ages, 1)
  
  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSd_2F(age, Low_chol, delta_matrix)
    B_tT[age,] <- B_BSd_2F(age, delta_matrix) 
    S_prj[age, 1] <- exp(-age * (A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1)) 
    mu_hat_prj[age,1] <- A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1
  }
  return(S_prj)
}

S_t_BSd_3F_proj <- function(x0, delta, kappa, sigma, r, proj_years){
  
  n_factors <- 3
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)
  
  B_tT <- matrix(NA, n_ages, n_factors)
  X_t_last <- KF_BSd_3F_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r)$X_t[,n_years+1]
  
  E_X_t1 <- exp(-kappa * proj_years) * X_t_last
  
  S_prj <- matrix(NA, n_ages, 1)
  mu_hat_prj <- matrix(NA, n_ages, 1)
  
  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSd_3F(age, Low_chol, delta_matrix)
    B_tT[age,] <- B_BSd_3F(age, delta_matrix) 
    S_prj[age, 1] <- exp(-age * (A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1)) 
    mu_hat_prj[age,1] <- A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1
  }
  return(S_prj)
}

# - AFNS independent factor model
S_t_AFNSi_proj <- function(x0, delta, kappa, sigma, r, proj_years){
  
  n_factors <- 3
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)
  
  B_tT <- matrix(NA, n_ages, n_factors)
  X_t_last <- KF_AFNSi_uKD(x0, delta, kappa, sigma, r)$X_t[,n_years+1]
  
  E_X_t1 <- exp(-kappa * proj_years) * X_t_last
  
  S_prj <- matrix(NA, n_ages, 1)
  mu_hat_prj <- matrix(NA, n_ages, 1)
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFNSi(age, sigma, delta)  
    B_tT[age,] <- B_AFNS(age,delta)  
    S_prj[age, 1] <- exp(-age * (A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1)) 
    mu_hat_prj[age,1] <- A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1
  }
  return(S_prj)
}

# - AFNS dependent factor model
S_t_AFNSd_proj <- function(x0, delta, kappa, sigma, r, proj_years){
  
  n_factors <- 3
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)
  
  B_tT <- matrix(NA, n_ages, n_factors)
  X_t_last <- KF_AFNSd_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r)$X_t[,n_years+1]
  
  E_X_t1 <- exp(-kappa * proj_years) * X_t_last
  
  S_prj <- matrix(NA, n_ages, 1)
  mu_hat_prj <- matrix(NA, n_ages, 1)
  
  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFNSg(age, Low_chol, delta)
    B_tT[age,] <- B_AFNS(age,delta)  
    S_prj[age, 1] <- exp(-age * (A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1)) 
    mu_hat_prj[age,1] <- A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1
  }
  return(S_prj)
}

# - CIR model
S_t_CIR_proj <- function(x0, delta, kappa, sigma, theta_Q, theta_P, r, proj_years){
  
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)
  
  B_tT <- matrix(NA, n_ages, n_factors)
  X_t_last <- KF_CIR_uKD(x0, delta, kappa, sigma, theta_Q, theta_P, r)$X_t[,n_years+1]
  
  E_X_t1 <- exp(-kappa * proj_years) + theta_P * (1 - exp(-kappa)) * sum(exp(-kappa * c(0:(proj_years-1)))) #exp(-kappa) * X_t_last + theta_P * (1 - exp(-kappa)) # - FIX THE RECURSION TO GENERALIZE
  
  S_prj <- matrix(NA, n_ages, 1)
  mu_hat_prj <- matrix(NA, n_ages, 1)
  
  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_CIR(age, theta_Q, sigma, delta)  
    B_tT[age,] <- B_CIR(age, sigma, delta)  
    S_prj[age, 1] <- exp(-age * (A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1)) 
    mu_hat_prj[age,1] <- A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1
  }
  return(S_prj)
}


