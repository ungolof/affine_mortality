library(HMDHFDplus)
library(plot3D)
library(expm)

#========================================= Functions for the continuous time affine models ========================


### ============= Factor loadings functions (we consider a measurement of the type y_t = A + B*X_t) ===============

# - Blackburn and Sherris
## - Independent
### - B(t,T) function for independent model (already divided by T-t, and coded as if t=0)
B_BSi <- function(T,delta){
  return((1 - exp(-delta * T))/(delta * T))
}

### - A(t,T) function for independent model (already divided by T-t, and coded as if t=0)
A_BSi <- function(T, sigma, delta){
  value <- - 0.5 * sum(((sigma^2) / (delta^3)) * ( 0.5 * (1 - exp(- 2 * delta * T)) -  2 * (1 - exp( - delta * T)) + delta * T)) / T
  return(value)
}

## - Dependent
### - 2 factors
#### - B(t,T) function for dependent model (already divided by T-t, and coded as if t=0)
B_BSd_2F <- function(T,delta_mat){   # - function of a matrix delta which is lower diagonal
  D1 <- delta_mat[2,1] / (delta_mat[1,1] - delta_mat[2,2])
  
  B1 <- - (1 + D1) * (1 - exp(- delta_mat[1,1] * T)) / (delta_mat[1,1] * T) + D1 * (1 - exp(- delta_mat[2,2] * T)) / (delta_mat[2,2] * T) 
  B2 <- - (1 - exp(- delta_mat[2,2] * T)) / (delta_mat[2,2] * T)
  
  return(-c(B1, B2))
}

#### -  A(t,T) function for dependent model (already divided by T-t, and coded as if t=0)
A_BSd_2F <- function(T, sigma_mat, delta_mat){   # - function of matrices delta and sigma which are lower diagonal
  
  # D terms as in Huang et. al. 2020
  D1 <- delta_mat[2,1] / (delta_mat[1,1] - delta_mat[2,2])
  
  # F terms as in Huang et. al. 2020
  F1 <- (sigma_mat[1,1] ^ 2) * ((1 + D1) ^2)
  F2 <- (sigma_mat[1,1] ^ 2) * (D1 ^2) - 2 * sigma_mat[1,1] * sigma_mat[2,1] * D1 + ((sigma_mat[2,1] ^ 2) + (sigma_mat[2,2] ^ 2)) 
  F4 <- - 2 * (sigma_mat[1,1] ^ 2) * (1 + D1) * D1 + 2 * sigma_mat[1,1] * sigma_mat[2,1] * (1 + D1)
  
  value <- - 0.5 * ( (F1 / delta_mat[1,1]) * ( 0.5 * (1 - exp( - 2 * T * delta_mat[1,1])) - 2 * (1 - exp( - T * delta_mat[1,1])) + delta_mat[1,1] * T) + 
                       (F2 / delta_mat[2,2]) * ( 0.5 * (1 - exp( - 2 * T * delta_mat[2,2])) - 2 * (1 - exp( - T * delta_mat[2,2])) + delta_mat[2,2] * T) + 
                       (F4 / (delta_mat[1,1] * delta_mat[2,2])) * (T - (1 - exp(- T * delta_mat[1,1])) / delta_mat[1,1] - (1 - exp(- T * delta_mat[2,2])) / delta_mat[2,2] + (1 - exp(- T * (delta_mat[1,1] + delta_mat[2,2]))) / (delta_mat[1,1] + delta_mat[2,2]))  ) / T
  
  return(value)
}

### - 3 factors
#### - B(t,T) function for dependent model (already divided by T-t, and coded as if t=0)
B_BSd_3F <- function(T,delta_mat){   # - function of a matrix delta which is lower diagonal
  D1 <- delta_mat[2,1] / (delta_mat[1,1] - delta_mat[2,2])
  D2 <- delta_mat[3,2] / (delta_mat[2,2] - delta_mat[3,3])
  D3 <- delta_mat[2,1] / (delta_mat[1,1] - delta_mat[3,3])
  D4 <- delta_mat[3,1] / (delta_mat[1,1] - delta_mat[3,3])
  D5 <- delta_mat[3,2] / (delta_mat[1,1] - delta_mat[3,3])
  E1 <- 1 + D1 + D1 * D5 + D4
  E2 <- D1 * (1 + D2)
  E3 <- D2 * D3 - D4
  
  B1 <- - E1 * (1 - exp(- delta_mat[1,1] * T)) / (delta_mat[1,1] * T) + E2 * (1 - exp(- delta_mat[2,2] * T)) / (delta_mat[2,2] * T) - E3 * (1 - exp(- delta_mat[3,3] * T)) / (delta_mat[3,3] * T)
  B2 <- - (1 + D2) * (1 - exp(- delta_mat[2,2] * T)) / (delta_mat[2,2] * T) + D2 * (1 - exp(- delta_mat[3,3] * T)) / (delta_mat[3,3] * T)
  B3 <- - (1 - exp(- delta_mat[3,3] * T)) / (delta_mat[3,3] * T)
  
  return(-c(B1, B2, B3))
}

#### -  A(t,T) function for dependent model (already divided by T-t, and coded as if t=0)
A_BSd_3F <- function(T, sigma_mat, delta_mat){   # - function of matrices delta and sigma which are lower diagonal
  
  # D terms as in Huang et. al. 2020
  D1 <- delta_mat[2,1] / (delta_mat[1,1] - delta_mat[2,2])
  D2 <- delta_mat[3,2] / (delta_mat[2,2] - delta_mat[3,3])
  D3 <- delta_mat[2,1] / (delta_mat[1,1] - delta_mat[3,3])
  D4 <- delta_mat[3,1] / (delta_mat[1,1] - delta_mat[3,3])
  D5 <- delta_mat[3,2] / (delta_mat[1,1] - delta_mat[3,3])
  
  # E terms as in Huang et. al. 2020
  E1 <- 1 + D1 + D1 * D5 + D4
  E2 <- D1 * (1 + D2)
  E3 <- D2 * D3 - D4
  
  # F terms as in Huang et. al. 2020
  F1 <- (sigma_mat[1,1] ^ 2) * (E1 ^2)
  F2 <- (sigma_mat[1,1] ^ 2) * (E2 ^2) - 2 * sigma_mat[1,1] * sigma_mat[2,1] * E2 * (1 + D2) + ((sigma_mat[2,1] ^ 2) + (sigma_mat[2,2] ^ 2)) * ((1 + D2) ^ 2)  
  F3 <- (sigma_mat[1,1] ^ 2) * (E3 ^2) - 2 * sigma_mat[1,1] * sigma_mat[2,1] * E3 * D2 + ((sigma_mat[2,1] ^ 2) + (sigma_mat[2,2] ^ 2)) * (D2 ^ 2) - 2 * (sigma_mat[2,1] * sigma_mat[3,1] + sigma_mat[2,2] * sigma_mat[3,2]) * D2 + 2 * sigma_mat[1,1] * sigma_mat[3,1] * E3 + ((sigma_mat[3,1] ^ 2) + (sigma_mat[3,2] ^ 2) + (sigma_mat[3,3] ^ 2))
  F4 <- - 2 * (sigma_mat[1,1] ^ 2) * E1 * E2 + 2 * sigma_mat[1,1] * sigma_mat[2,1] * E1 * (1 + D2)
  F5 <-   2 * (sigma_mat[1,1] ^ 2) * E1 * E3 - 2 * sigma_mat[1,1] * sigma_mat[2,1] * E1 * D2 + 2 * sigma_mat[1,1] * sigma_mat[3,1] * E1
  F6 <- - 2 * (sigma_mat[1,1] ^ 2) * E2 * E3 + 2 * sigma_mat[1,1] * sigma_mat[2,1] * (E3 * (1 + D2) + E2 * D2) - 2 * sigma_mat[1,1] * sigma_mat[3,1] * E2 + 2 * (sigma_mat[2,1] * sigma_mat[3,1] + sigma_mat[2,2] * sigma_mat[3,2]) * (1 + D2) - 2 * ((sigma_mat[2,1] ^ 2) + (sigma_mat[2,2] ^ 2)) * D2 * (1 + D2)
  
  value <- - 0.5 * ( (F1 / (delta_mat[1,1]^3)) * ( 0.5 * (1 - exp( - 2 * T * delta_mat[1,1])) - 2 * (1 - exp( - T * delta_mat[1,1])) + delta_mat[1,1] * T) + 
                       (F2 / (delta_mat[2,2]^3)) * ( 0.5 * (1 - exp( - 2 * T * delta_mat[2,2])) - 2 * (1 - exp( - T * delta_mat[2,2])) + delta_mat[2,2] * T) + 
                       (F3 / (delta_mat[3,3]^3)) * ( 0.5 * (1 - exp( - 2 * T * delta_mat[3,3])) - 2 * (1 - exp( - T * delta_mat[3,3])) + delta_mat[3,3] * T) + 
                       (F4 / (delta_mat[1,1] * delta_mat[2,2])) * (T - (1 - exp(- T * delta_mat[1,1])) / delta_mat[1,1] - (1 - exp(- T * delta_mat[2,2])) / delta_mat[2,2] + (1 - exp(- T * (delta_mat[1,1] + delta_mat[2,2]))) / (delta_mat[1,1] + delta_mat[2,2])) + 
                       (F5 / (delta_mat[1,1] * delta_mat[3,3])) * (T - (1 - exp(- T * delta_mat[1,1])) / delta_mat[1,1] - (1 - exp(- T * delta_mat[3,3])) / delta_mat[3,3] + (1 - exp(- T * (delta_mat[1,1] + delta_mat[3,3]))) / (delta_mat[1,1] + delta_mat[3,3])) + 
                       (F6 / (delta_mat[2,2] * delta_mat[3,3])) * (T - (1 - exp(- T * delta_mat[2,2])) / delta_mat[2,2] - (1 - exp(- T * delta_mat[3,3])) / delta_mat[3,3] + (1 - exp(- T * (delta_mat[2,2] + delta_mat[3,3]))) / (delta_mat[2,2] + delta_mat[3,3])) ) / T
  
  return(value)
}

# - AFNS
## - Independent
### - B(t,T) function for AFNS (already divided by T-t, and coded as if t=0)
B_AFNS <- function(T, delta){
  B_1 <- -1
  B_2 <- -(1 - exp(-delta * T))/(delta * T)
  B_3 <- exp(- delta * T) - (1 - exp(-delta * T))/(delta * T)
  return(-c(B_1, B_2, B_3))
}

### - A(t,T) function for independent AFNS model (already divided by T-t, and coded as if t=0)
A_AFNSi <- function(T, sigma, delta){
  value <- (sigma[1]^2) * (T^2) / 6 + (sigma[2]^2) * (0.5 / (delta^2) - (1 - exp(-delta * T))/(T * (delta^3)) + 0.25 * (1 - exp(-delta * 2 * T))/(T * (delta^3)) ) + 
    (sigma[3]^2) * (0.5 / (delta^2) + exp(- delta * T)/(delta^2) - 0.25 * T * exp(- 2 * T * delta) / delta - 0.75 * exp(- 2 * T * delta) / (delta^2) - 2 * (1 - exp( - delta * T)) / (T * (delta^3)) + 0.625 * (1 - exp( - 2 * delta * T)) / (T * (delta^3)) )
  return(- value)
}

## - AFNS for dependent (and independent) factor models
### -  A(t,T) function for dependent model (already divided by T-t), general version applicable also to independent models
# - from Christensen et. al. 2011
A_AFNSg <- function(T, sigma_mat, delta){
  
  A_bar <- (sigma_mat[1,1] ^ 2) + (sigma_mat[1,2] ^ 2) + (sigma_mat[1,3] ^ 2)
  B_bar <- (sigma_mat[2,1] ^ 2) + (sigma_mat[2,2] ^ 2) + (sigma_mat[2,3] ^ 2)
  C_bar <- (sigma_mat[3,1] ^ 2) + (sigma_mat[3,2] ^ 2) + (sigma_mat[3,3] ^ 2)
  D_bar <- sigma_mat[1,1] * sigma_mat[2,1] + sigma_mat[1,2] * sigma_mat[2,2] + sigma_mat[1,3] * sigma_mat[2,3]
  E_bar <- sigma_mat[1,1] * sigma_mat[3,1] + sigma_mat[1,2] * sigma_mat[3,2] + sigma_mat[1,3] * sigma_mat[3,3]
  F_bar <- sigma_mat[2,1] * sigma_mat[3,1] + sigma_mat[2,2] * sigma_mat[3,2] + sigma_mat[2,3] * sigma_mat[3,3]
  
  value <- A_bar * (T^2) / 6 + B_bar * (0.5 / (delta^2) - (1 - exp(-delta * T))/(T * (delta^3)) + 0.25 * (1 - exp(-delta * 2 * T))/(T * (delta^3)) ) + 
    C_bar * (0.5 / (delta^2) + exp(- delta * T) / (delta^2) - 0.25 * T * exp(- 2 * T * delta) / delta - 0.75 * exp(- 2 * T * delta) / (delta^2) - 2 * (1 - exp( - delta * T)) / (T * (delta^3)) + 0.625 * (1 - exp( - 2 * delta * T)) / (T * (delta^3)) ) +
    D_bar * (0.5 * T / delta + exp(- delta * T) / (delta ^ 2) - (1 - exp(- T * delta)) / (T * (delta ^ 3)) ) + 
    E_bar * (3 * exp(- T * delta) / (delta ^ 2) + 0.5 * T / delta + T * exp(- T * delta) / delta - 3 * (1 - exp(- T * delta)) / (T * (delta ^ 3)) ) + 
    F_bar * ( 1 / (delta ^ 2) + exp(- delta * T) / (delta ^ 2) - 0.5 * exp(- 2 * delta * T) / (delta ^ 2) - 3 * (1 - exp(- T * delta)) / (T * (delta ^ 3)) + 0.75 * (1 - exp(- 2 * T * delta)) / (T * (delta ^ 3)) )
  
  return(- value)
}

# - AFGNS
## - Independent factors
A_AFGNSi <- function(T, sigma, delta1, delta2){
  value <- (sigma[1]^2) * (T^2) / 6 + (sigma[2]^2) * (0.5 / (delta1^2) - (1 - exp(-delta1 * T))/(T * (delta1^3)) + 0.25 * (1 - exp(-delta1 * 2 * T))/(T * (delta1^3)) ) + 
                                      (sigma[3]^2) * (0.5 / (delta2^2) - (1 - exp(-delta2 * T))/(T * (delta2^3)) + 0.25 * (1 - exp(-delta2 * 2 * T))/(T * (delta2^3)) ) + 
           (sigma[4]^2) * (0.5 / (delta1^2) + exp(- delta1 * T)/(delta1^2) - 0.25 * T * exp(- 2 * T * delta1) / delta1 - 0.75 * exp(- 2 * T * delta1) / (delta1^2) - 2 * (1 - exp( - delta1 * T)) / (T * (delta1^3)) + 0.625 * (1 - exp( - 2 * delta1 * T)) / (T * (delta1^3)) ) + 
           (sigma[5]^2) * (0.5 / (delta2^2) + exp(- delta2 * T)/(delta2^2) - 0.25 * T * exp(- 2 * T * delta2) / delta2 - 0.75 * exp(- 2 * T * delta2) / (delta2^2) - 2 * (1 - exp( - delta2 * T)) / (T * (delta2^3)) + 0.625 * (1 - exp( - 2 * delta2 * T)) / (T * (delta2^3)) )
  return(- value)
}

A_AFGNSg <- function(T, sigma_mat, delta1, delta2){
  A_bar <- (sigma_mat[1,1] ^ 2) + (sigma_mat[1,2] ^ 2) + (sigma_mat[1,3] ^ 2) + (sigma_mat[1,4] ^ 2) + (sigma_mat[1,5] ^ 2)
  B_bar <- (sigma_mat[2,1] ^ 2) + (sigma_mat[2,2] ^ 2) + (sigma_mat[2,3] ^ 2) + (sigma_mat[2,4] ^ 2) + (sigma_mat[2,5] ^ 2)
  C_bar <- (sigma_mat[3,1] ^ 2) + (sigma_mat[3,2] ^ 2) + (sigma_mat[3,3] ^ 2) + (sigma_mat[3,4] ^ 2) + (sigma_mat[3,5] ^ 2)
  D_bar <- (sigma_mat[4,1] ^ 2) + (sigma_mat[4,2] ^ 2) + (sigma_mat[4,3] ^ 2) + (sigma_mat[4,4] ^ 2) + (sigma_mat[4,5] ^ 2)
  E_bar <- (sigma_mat[5,1] ^ 2) + (sigma_mat[5,2] ^ 2) + (sigma_mat[5,3] ^ 2) + (sigma_mat[5,4] ^ 2) + (sigma_mat[5,5] ^ 2)
  
  F_bar <- sigma_mat[1,1] * sigma_mat[2,1] + sigma_mat[1,2] * sigma_mat[2,2] + sigma_mat[1,3] * sigma_mat[2,3] + sigma_mat[1,4] * sigma_mat[2,4] + sigma_mat[1,5] * sigma_mat[2,5]
  G_bar <- sigma_mat[1,1] * sigma_mat[3,1] + sigma_mat[1,2] * sigma_mat[3,2] + sigma_mat[1,3] * sigma_mat[3,3] + sigma_mat[1,4] * sigma_mat[3,4] + sigma_mat[1,5] * sigma_mat[3,5]
  H_bar <- sigma_mat[1,1] * sigma_mat[4,1] + sigma_mat[1,2] * sigma_mat[4,2] + sigma_mat[1,3] * sigma_mat[4,3] + sigma_mat[1,4] * sigma_mat[4,4] + sigma_mat[1,5] * sigma_mat[4,5]
  I_bar <- sigma_mat[1,1] * sigma_mat[5,1] + sigma_mat[1,2] * sigma_mat[5,2] + sigma_mat[1,3] * sigma_mat[5,3] + sigma_mat[1,4] * sigma_mat[5,4] + sigma_mat[1,5] * sigma_mat[5,5]
  J_bar <- sigma_mat[2,1] * sigma_mat[3,1] + sigma_mat[2,2] * sigma_mat[3,2] + sigma_mat[2,3] * sigma_mat[3,3] + sigma_mat[2,4] * sigma_mat[3,4] + sigma_mat[2,5] * sigma_mat[3,5]
  K_bar <- sigma_mat[2,1] * sigma_mat[4,1] + sigma_mat[2,2] * sigma_mat[4,2] + sigma_mat[2,3] * sigma_mat[4,3] + sigma_mat[2,4] * sigma_mat[4,4] + sigma_mat[2,5] * sigma_mat[4,5]
  L_bar <- sigma_mat[2,1] * sigma_mat[5,1] + sigma_mat[2,2] * sigma_mat[5,2] + sigma_mat[2,3] * sigma_mat[5,3] + sigma_mat[2,4] * sigma_mat[5,4] + sigma_mat[2,5] * sigma_mat[5,5]
  M_bar <- sigma_mat[3,1] * sigma_mat[4,1] + sigma_mat[3,2] * sigma_mat[4,2] + sigma_mat[3,3] * sigma_mat[4,3] + sigma_mat[3,4] * sigma_mat[4,4] + sigma_mat[3,5] * sigma_mat[4,5]
  N_bar <- sigma_mat[3,1] * sigma_mat[5,1] + sigma_mat[3,2] * sigma_mat[5,2] + sigma_mat[3,3] * sigma_mat[5,3] + sigma_mat[3,4] * sigma_mat[5,4] + sigma_mat[3,5] * sigma_mat[5,5]
  O_bar <- sigma_mat[4,1] * sigma_mat[5,1] + sigma_mat[4,2] * sigma_mat[5,2] + sigma_mat[4,3] * sigma_mat[5,3] + sigma_mat[4,4] * sigma_mat[5,4] + sigma_mat[4,5] * sigma_mat[5,5]
  
  value <- A_bar * (T^2) / 6 + B_bar * (0.5 / (delta1^2) - (1 - exp(- delta1 * T))/(T * (delta1^3)) + 0.25 * (1 - exp(-delta1 * 2 * T))/(T * (delta1^3)) ) + 
                               C_bar * (0.5 / (delta2^2) - (1 - exp(- delta2 * T))/(T * (delta2^3)) + 0.25 * (1 - exp(-delta2 * 2 * T))/(T * (delta2^3)) ) + 
           D_bar * (0.5 / (delta1^2) + exp(- delta1 * T)/(delta1^2) - 0.25 * T * exp(- 2 * T * delta1) / delta1 - 0.75 * exp(- 2 * T * delta1) / (delta1^2) - 2 * (1 - exp( - delta1 * T)) / (T * (delta1^3)) + 0.625 * (1 - exp( - 2 * delta1 * T)) / (T * (delta1^3)) ) + 
           E_bar * (0.5 / (delta2^2) + exp(- delta2 * T)/(delta2^2) - 0.25 * T * exp(- 2 * T * delta2) / delta2 - 0.75 * exp(- 2 * T * delta2) / (delta2^2) - 2 * (1 - exp( - delta2 * T)) / (T * (delta2^3)) + 0.625 * (1 - exp( - 2 * delta2 * T)) / (T * (delta2^3)) ) +
           F_bar * (0.5 * T / delta1 + exp(- delta1 * T) / (delta1 ^ 2) - (1 - exp(- T * delta1)) / (T * (delta1 ^ 3)) ) + 
           G_bar * (0.5 * T / delta2 + exp(- delta2 * T) / (delta2 ^ 2) - (1 - exp(- T * delta2)) / (T * (delta2 ^ 3)) ) + 
           H_bar * (3 * exp(- T * delta1) / (delta1 ^ 2) + 0.5 * T / delta1 + T * exp(- T * delta1) / delta1 - 3 * (1 - exp(- T * delta1)) / (T * (delta1 ^ 3)) ) + 
           I_bar * (3 * exp(- T * delta2) / (delta2 ^ 2) + 0.5 * T / delta2 + T * exp(- T * delta2) / delta2 - 3 * (1 - exp(- T * delta2)) / (T * (delta2 ^ 3)) ) + 
           J_bar * (1 / (delta1 * delta2) - (1 - exp(- T * delta1)) / (T * (delta1 ^ 2) * delta2) - (1 - exp(- T * delta2)) / (T * (delta2 ^ 2) * delta1) + (1 - exp(- T * (delta1 + delta2))) / (T * delta1 * delta2 * (delta1 + delta2) ) ) + 
           K_bar * (1 / (delta1 ^ 2) + exp(- delta1 * T) / (delta1 ^ 2) - 0.5 * exp(- 2 * delta1 * T) / (delta1 ^ 2) - 3 * (1 - exp(- T * delta1)) / (T * (delta1 ^ 3)) + 0.75 * (1 - exp(- 2 * T * delta1)) / (T * (delta1 ^ 3)) ) + 
           L_bar * (1 / (delta1 * delta2) + exp(- T * delta2) / (delta1 * delta2) - exp(- T * (delta1 + delta2)) / (delta1 * (delta1 + delta2)) - (1 - exp(- T * delta1)) / (T * (delta1 ^ 2) * delta2) - 2 * (1 - exp(- T * delta2)) / (T * (delta2 ^ 2) * delta1) + (delta1 + 2 * delta2) * (1 - exp(- T * (delta1 + delta2))) / (delta1 * delta2 * ((delta1 + delta2)^2)) ) + 
           M_bar * (1 / (delta1 * delta2) + exp(- T * delta1) / (delta1 * delta2) - exp(- T * (delta1 + delta2)) / (delta2 * (delta1 + delta2)) - (1 - exp(- T * delta2)) / (T * (delta2 ^ 2) * delta1) - 2 * (1 - exp(- T * delta1)) / (T * (delta1 ^ 2) * delta2) + (delta2 + 2 * delta1) * (1 - exp(- T * (delta1 + delta2))) / (delta1 * delta2 * ((delta1 + delta2)^2)) ) +
           N_bar * (1 / (delta2 ^ 2) + exp(- delta2 * T) / (delta2 ^ 2) - 0.5 * exp(- 2 * delta2 * T) / (delta2 ^ 2) - 3 * (1 - exp(- T * delta2)) / (T * (delta2 ^ 3)) + 0.75 * (1 - exp(- 2 * T * delta2)) / (T * (delta2 ^ 3)) ) + 
           O_bar * (1 / (delta1 * delta2) + (exp(- T * delta1) + exp(- T * delta2)) / (delta1 * delta2) - (1/delta1 + 1/delta2) * exp(- T * (delta1 * delta2)) / (delta1 + delta2) - 2 * exp(-T * (delta1 + delta2)) / ((delta1 + delta2)^2) - T * exp(-T * (delta1 + delta2)) / (delta1 + delta2) - 2 * (1 - exp(- T * delta1)) / (T * delta2 * (delta1^2)) - 2 * (1 - exp(- T * delta2)) / (T * delta1 * (delta2^2)) + 
                      2 * (1 - exp(- T * (delta1 + delta2))) / (T * ((delta1 + delta2)^3))  + (1 / delta1 + 1 / delta2) * (1 - exp(- T * (delta1 + delta2))) / (T * ((delta1 + delta2)^2)) + (1 - exp(- T * (delta1 + delta2))) / (T * delta1 * delta2 * (delta1 + delta2) ) )
    
    return(- value)
}



# - CIR model
### -  A(t,T) function (already divided by T-t)
A_CIR <- function(T, theta_Q, sigma, delta){
  value <- - 2 * sum((delta * theta_Q / (sigma^2)) * (log(2) + 0.5 * log((delta ^ 2) + 2 * (sigma ^ 2)) + 0.5 * T * (delta + sqrt((delta ^ 2) + 2 * (sigma ^ 2))) - log( (delta + sqrt((delta ^ 2) + 2 * (sigma ^ 2))) * (exp(T * sqrt((delta ^ 2) + 2 * (sigma ^ 2))) - 1) + 2 * sqrt((delta ^ 2) + 2 * (sigma ^ 2))  ))) / T
  return(value)
}

### - Same A_CIR as before, but using ^0.5 instead of sqrt() for computational reasons
A_CIR <- function(T, theta_Q, sigma, delta){
  value <- - 2 * sum((delta * theta_Q / (sigma^2)) * (log(2) + 0.5 * log((delta ^ 2) + 2 * (sigma ^ 2)) + 0.5 * T * (delta + ((delta ^ 2) + 2 * (sigma ^ 2))^0.5) - log( (delta + ((delta ^ 2) + 2 * (sigma ^ 2))^0.5) * (exp(T * ((delta ^ 2) + 2 * (sigma ^ 2))^0.5) - 1) + 2 * ((delta ^ 2) + 2 * (sigma ^ 2))^0.5  ))) / T
  return(value)
}

### - Same B_CIR as before, but using ^0.5 instead of sqrt() for computational reasons
B_CIR <- function(T, sigma, delta){
  B <- 2 * (exp(T * ((delta ^ 2) + 2 * (sigma ^ 2))^0.5) - 1) / ((delta + ((delta ^ 2) + 2 * (sigma ^ 2))^0.5) * (exp(T * ((delta ^ 2) + 2 * (sigma ^ 2))^0.5) - 1) + 2 * ((delta ^ 2) + 2 * (sigma ^ 2))^0.5) / T
  return(B)
}



#================================================ Measurement error functions ======================================

# - Measurement error as in Blackburn and Sherris with exponentiated parameters
meas_err_BS <- function(r1, r2, rc){
  n_ages <- nrow(mu_bar)
  H <- matrix(0, n_ages, n_ages)
  
  for(age in 1:n_ages){
    # - Original diag from Blackburn-Sherris (2013)
    H[age,age] <- exp(rc) + exp(r1) * sum(exp(exp(r2) * c(1:age))) / age
  }
  return(H)
}

# - Measurement error as in Xu, Sherris and Ziveyi (with exponentiated parameters to ensure positivity)
meas_err_XSZ <- function(r1, r2){
  n_ages <- nrow(mu_bar)
  for(age in 1:n_ages){
    # - Original diag from Blackburn-Sherris (2013)
    H[age,age] <- exp(r1) * exp(exp(r2) * age)
  }
  return(H)
}

#======================================== Accessory functions =================================

# - Fill lower triangular matrix (by row) from a vector: this turns out useful when coding delta and Sigma
# - in the dependent AFNS and Blackburn-Sherris model.
low_trg_fill <- function(par_vector){
  n_factors <- (-1 + sqrt(1 + 8 * length(par_vector))) * 0.5
  par_matrix <- matrix(0, n_factors, n_factors) # - lower triangular matrix of delta coefficients
  count_vec <- 1
  for(row in 1:n_factors){
    for(col in 1:row){
      par_matrix[row,col] <- par_vector[count_vec]
      count_vec <- count_vec + 1
    }
  }
  return(par_matrix)
}

# - Fill lower triangular matrix with 0s on the main diagonal (this is also useful when coding Sigma)
low_trg_fill_0diag <- function(par_vector){
  n_factors <- (1 + sqrt(1 + 8 * length(par_vector))) * 0.5
  par_matrix <- matrix(0, n_factors, n_factors) # - lower triangular matrix of delta coefficients
  count_vec <- 1
  for(row in 2:n_factors){
    for(col in 1:(row-1)){
      par_matrix[row,col] <- par_vector[count_vec]
      count_vec <- count_vec + 1
    }
  }
  return(par_matrix)
}

# - Function to get the parameters which are to be input in the optimization routine (as starting value), and viceversa
# - cov2par returns the elements of the parameter vector for fitting from the elements of the lower triangular part of the covariance matrix
# - The function yields two vectors: the first dg_l_Sigma_chol gives the logarithm of the elements on the diagonal of the
# - Cholesky decomposition of Sigma (which must be positive to ensure uniqueness), while odg_Sigma_chol provides the elements
# - below the diagonal of the Cholesky decomposition of Sigma, which do not need to be constrained
# - (example of 3x3 cov. Input: sigma^2_11, sigma^2_22, sigma^2_33, sigma_21, sigma_31, sigma_32; i.e. first input elements on the diagonal,
# - output: log(Chol(Sigma)[1,1]), log(Chol(Sigma)[2,2]), log(Chol(Sigma)[3,3]), Chol(Sigma)[2,1], Chol(Sigma)[2,2], Chol(Sigma)[3,3])

cov2par <- function(Sigma_el){
  n_factors <- (-1 + sqrt(1 + 8 * length(Sigma_el))) * 0.5
  
  Sigma_diffusion <- matrix(diag(Sigma_el[1:n_factors]), n_factors, n_factors)
  Sigma_diffusion <- Sigma_diffusion + low_trg_fill_0diag(Sigma_el[(n_factors+1):length(Sigma_el)]) + t(low_trg_fill_0diag(Sigma_el[(n_factors+1):length(Sigma_el)]))
  
  Low_Sigma <- t(chol(Sigma_diffusion))
  odg_Sigma_chol <- rep(0, length(Sigma_el) - n_factors)
  
  count_vec <- 1
  for(row in 2:n_factors){
    for(col in 1:(row-1)){
      odg_Sigma_chol[count_vec] <- Low_Sigma[row,col]
      count_vec <- count_vec + 1
    }
  }
  
  return(list(dg_l_Sigma_chol=log(diag(Low_Sigma)), odg_Sigma_chol=odg_Sigma_chol))
}

# - From parameter vector for fitting function to elements of the covariance matrix (input: log(diag(Chol(Low_chol))), 
# - off-diagonal elements of the Chol.). Output: st. deviation and covariances
parest2cov <- function(dg_l_Sigma_chol, odg_Sigma_chol){
  n_factors <- (-1 + sqrt(1 + 8 * (length(dg_l_Sigma_chol) + length(odg_Sigma_chol)))) * 0.5
  
  Low_Chol <- matrix(0, n_factors, n_factors)
  diag(Low_Chol) <- exp(dg_l_Sigma_chol)
  Low_Chol <- Low_Chol + low_trg_fill_0diag(odg_Sigma_chol)
  
  Sigma_diffusion <- matrix(NA, n_factors, n_factors)
  Sigma_diffusion <- Low_Chol %*% t(Low_Chol)
  diag(Sigma_diffusion) <- sqrt(diag(Sigma_diffusion))
  
  Sigma_el <- rep(0, (length(dg_l_Sigma_chol) + length(odg_Sigma_chol)))
  count_vec <- 1
  for(row in 1:n_factors){
    for(col in 1:row){
      Sigma_el[count_vec] <- Sigma_diffusion[row,col]
      count_vec <- count_vec + 1
    }
  }
  
  return(Sigma_el)
}


get_L <- function(sigma_dg, Sigma_cov){
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)
  return(Low_chol)
}

# - This function takes the maximum between 1e-10 and X elementwise for the vector X. This is useful to bound below the value of X
# - by construction in the CIR model
l_bound <- function(X){
  Xb <- X
  for(j in 1:length(X)){
    Xb[j] <- max(1e-10, X[j])
  }
  return(Xb)
}

l_bound_mod <- function(X,Y,low_bound=0){
  Xb <- X
  for(j in 1:length(X)){
    Xb[j] <- ifelse(X[j] > low_bound, X[j], Y[j])
  }
  return(Xb) 
}


#========= Objective functions for the negative log-likelihoods (Without constants) for direct minimization ==================

#============= Negative log-likelihoods (to be directly minimized) for the independent BS models =============================


# - Basic univariate version of Koopman and Durbin (2002)
# - r_1 and r_2 have been enforced to be positive in order to allow for an increasing
#   Poisson variability with respect to the age
nLL_BSi_uKD <- function(vdParameters, mu_bar){
  
  n_factors <- (length(vdParameters) - 3) / 4
  
  x0 <- vdParameters[1:n_factors]
  delta <- vdParameters[(n_factors+1):(n_factors*2)]
  kappa <- vdParameters[(n_factors*2+1):(n_factors*3)] # - take logs to ensure positivity
  l_sigma <- vdParameters[(n_factors*3+1):(n_factors*4)]
  r_c <- vdParameters[n_factors*4 + 3]
  r_1 <- vdParameters[n_factors*4 + 1]
  r_2 <- vdParameters[n_factors*4 + 2]
  
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years
  
  v_ti <- mu_bar
  F_ti <- mu_bar
  
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
  R <- diag(((exp(l_sigma*2)) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)
  
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
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i 
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti
      
    }
  }
  
  log_lik <- sum(log(F_ti) + (v_ti^2) / F_ti)
  return(log_lik)
}


# - Basic univariate version of Koopman and Durbin (2002)
nLL_BSd_2F_uKD <- function(vdParameters, mu_bar){
  n_factors <- 2
  # - Parameters
  x0 <- vdParameters[1:n_factors]
  delta <- vdParameters[(n_factors+1):(n_factors * (3 + n_factors) / 2)]####### - TO BE CHANGED
  kappa <- vdParameters[(n_factors * (3 + n_factors) / 2 + 1):(n_factors * (5 + n_factors) / 2)] # - take logs to ensure positivity
  dg_l_Sigma_chol <- vdParameters[(n_factors * (5 + n_factors) / 2 + 1):(n_factors * (7 + n_factors) / 2)] # - log diagonal elements of the lower cholesky dec.
  odg_Sigma_chol <- vdParameters[(n_factors * (7 + n_factors) / 2 + 1):(n_factors * (3 + n_factors))] # - off diagonal elements of the lower cholesky dec.
  r_c <- vdParameters[n_factors * (3 + n_factors) + 3]
  r_1 <- vdParameters[n_factors * (3 + n_factors) + 1]
  r_2 <- vdParameters[n_factors * (3 + n_factors) + 2]
  
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years
  
  delta_matrix <- low_trg_fill(delta)
  
  v_ti <- mu_bar
  F_ti <- mu_bar
  
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
    A_tT[age,1] <- A_BSd_2F(age, Low_chol, delta_matrix)####### A_ind(age, exp(l_sigma), delta)  
    B_tT[age,] <- B_BSd_2F(age, delta_matrix)  ###### B_ind(age,delta)  
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
      F_ti[i,t] <- t(B_tT[i,]) %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i) 
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti
      
    }
  }
  
  log_lik <- sum(log(F_ti) + (v_ti^2) / F_ti)
  return(log_lik)
}

nLL_BSd_3F_uKD <- function(vdParameters, mu_bar){
  n_factors <- 3
  # - Parameters
  x_0 <- vdParameters[1:n_factors]
  delta <- vdParameters[(n_factors+1):(n_factors * (3 + n_factors) / 2)]####### - TO BE CHANGED
  kappa <- vdParameters[(n_factors * (3 + n_factors) / 2 + 1):(n_factors * (5 + n_factors) / 2)] # - take logs to ensure positivity
  dg_l_Sigma_chol <- vdParameters[(n_factors * (5 + n_factors) / 2 + 1):(n_factors * (7 + n_factors) / 2)] # - log diagonal elements of the lower cholesky dec.
  odg_Sigma_chol <- vdParameters[(n_factors * (7 + n_factors) / 2 + 1):(n_factors * (3 + n_factors))] # - off diagonal elements of the lower cholesky dec.
  r_c <- vdParameters[n_factors * (3 + n_factors) + 3]
  r_1 <- vdParameters[n_factors * (3 + n_factors) + 1]
  r_2 <- vdParameters[n_factors * (3 + n_factors) + 2]
  
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years
  
  delta_matrix <- low_trg_fill(delta)
  
  v_ti <- mu_bar
  F_ti <- mu_bar
  
  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance
  
  # - Initialize X and Sigma
  x_ti <- x_0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)
  
  Phi <- diag(exp(-kappa), n_factors) 
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
    A_tT[age,1] <- A_BSd_3F(age, Low_chol, delta_matrix)####### A_ind(age, exp(l_sigma), delta)  
    B_tT[age,] <- B_BSd_3F(age, delta_matrix)  ###### B_ind(age,delta)  
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
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i 
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti
      
    }
  }
  
  log_lik <- sum(log(F_ti) + (v_ti^2) / F_ti)
  return(log_lik)
}


#================== Negative log-likelihoods for the independent AFNS models =================


# - Basic univariate version of Koopman and Durbin (2002)
# - Note: this function is exactly the same as for the Blackburn-Sherris, 
#   with the only change in the Section Parameters and calculation of A(t,T) and B(t,T)
nLL_AFNSi_uKD <- function(vdParameters, mu_bar){
  n_factors <- 3
  # - Parameters
  x0 <- vdParameters[1:n_factors]
  delta <- vdParameters[n_factors+1]
  kappa <- vdParameters[(n_factors+2):(n_factors*2+1)] # - take logs to ensure positivity
  l_sigma <- vdParameters[(n_factors*2+2):(n_factors*3+1)]
  r_c <- vdParameters[n_factors*3 + 4]
  r_1 <- vdParameters[n_factors*3 + 2]
  r_2 <- vdParameters[n_factors*3 + 3]
  
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years
  
  v_ti <- mu_bar
  F_ti <- mu_bar
  
  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance
  
  # - Initialize X and Sigma
  x_ti <- x0 #init_X
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
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))  

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t] 
      
      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) 

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i 
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti
      
    }
  }
  
  log_lik <- sum(log(F_ti) + (v_ti^2) / F_ti)
  return(log_lik)
}

#================== Negative log-likelihoods for the dependent AFNS models =================

# - Basic univariate version of Koopman and Durbin (2002)
nLL_AFNSd_uKD <- function(vdParameters, mu_bar){
  n_factors <- 3
  # - Parameters
  x0 <- vdParameters[1:n_factors]
  delta <- vdParameters[n_factors+1]
  kappa <- vdParameters[(n_factors+2):(2 * n_factors + 1)] # - take logs to ensure positivity    n_factors * (3 + n_factors) / 2 + 1
  dg_l_Sigma_chol <- vdParameters[(2 * n_factors + 2):(3 * n_factors + 1)] # - log diagonal elements of the lower cholesky dec.
  odg_Sigma_chol <- vdParameters[(3 * n_factors + 2):(n_factors * (5 + n_factors) / 2 + 1)] # - off diagonal elements of the lower cholesky dec.
  r_c <- vdParameters[n_factors * (5 + n_factors) / 2 + 4]
  r_1 <- vdParameters[n_factors * (5 + n_factors) / 2 + 2]
  r_2 <- vdParameters[n_factors * (5 + n_factors) / 2 + 3]
  
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years
  
  v_ti <- mu_bar
  F_ti <- mu_bar
  
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
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i 
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti
      
    }
  }
  
  log_lik <- sum(log(F_ti) + (v_ti^2) / F_ti)
  return(log_lik)
}

#======================= Negative log-likelihoods for the CIR model (independent factors) =====================================

# - Version of the CIR univariate filtering which includes a lower bound for the states, as well as exp(l_kappa) in order to avoid negative definite R
nLL_CIR_uKD_bd <- function(vdParameters, mu_bar){
  n_factors <- (length(vdParameters) - 3) / 6
  # - Parameters
  l_x0 <- vdParameters[1:n_factors]
  delta <- vdParameters[(n_factors+1):(n_factors*2)]
  l_kappa <- vdParameters[(n_factors*2+1):(n_factors*3)] # - take logs to ensure positivity
  l_sigma <- vdParameters[(n_factors*3+1):(n_factors*4)]
  l_theta_Q <- vdParameters[(n_factors*4+1):(n_factors*5)]
  l_theta_P <- vdParameters[(n_factors*5+1):(n_factors*6)]
  r_c <- vdParameters[n_factors*6 + 3]
  r_1 <- vdParameters[n_factors*6 + 1]
  r_2 <- vdParameters[n_factors*6 + 2]
  
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years
  
  v_ti <- mu_bar
  F_ti <- mu_bar
  
  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors * (n_years + 1)) # - Factor covariance
  
  # - Initialize X and Sigma
  #x_ti <- x_0 #init_X
  x_ti <- exp(l_x0)
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
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))  

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t]
      x_ti <- l_bound(x_ti)
      
      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) 

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i 
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti
      
    }
  }
  
  log_lik <- sum(log(F_ti) + (v_ti^2) / F_ti)
  return(log_lik)
}














