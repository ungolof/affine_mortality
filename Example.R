#============ - Sample code for implementation - ==================================================


#======================== - Call files for estimation - ================================
source('Est_fun.R') # - Read all functions to optimize the models and the packages needed
source('Full_opt.R') # - Read all functions which perform the full optimization
source('Coordinate_Ascent.R') # - Read all functions which perform optimization by coordinate ascent
source('GoodnessofFit.R') # - Read functions to compute the Goodness of fit measures and of residuals
source('FilterSmoother.R') # - Read functions to perform filtering and smoothing of the latent state variables
source('StErr_Bootstrap.R') # - Read functions to perform bootstrap estimation of the covariance of the parameter estimates
source('StErr_MultImp.R') # - Read functions to estimate the covariance of the parameter estimates by Multiple Imputation
source('Projection.R') # - Read functions yielding the projected cohort survival curves



#======================== - Variable definition - ================================
AgeRange <- c(50:99) # - Set of ages under analysis
CohortRange <- seq(1883, 1915)  # - Set of cohorts under analysis
PeriodRange <- CohortRange + AgeRange[1] # - Period data

d_xt <- matrix(NA, length(AgeRange),length(CohortRange)) # - Deaths
E_xt <- matrix(NA, length(AgeRange),length(CohortRange)) # - Central Exposures
m_xt <- matrix(NA, length(AgeRange),length(CohortRange)) # - central death rates matrix

S_xt <- matrix(NA, length(AgeRange),length(CohortRange)) # - Survival function for each age and cohort
mu_bar_USA <- matrix(NA, length(AgeRange),length(CohortRange)) # - mu_bar for each age (row) and cohort (column)

n_ages <- length(AgeRange) # - Number of ages (N)
n_years <- length(CohortRange) # - Number of cohort years (T)

colnames(m_xt) <- colnames(E_xt) <- colnames(d_xt) <- CohortRange
rownames(m_xt) <- rownames(E_xt) <- rownames(d_xt) <- AgeRange


#================== - Read file deaths and exposures directly from the HMD dataset=========
# - Read file deaths and exposures directly from the HMD dataset
deaths <- readHMDweb(CNTRY = "USA", item = "Deaths_1x1", username, password, fixup = TRUE)
exposures <- readHMDweb(CNTRY = "USA", item = "Exposures_1x1", username, password, fixup = TRUE)


# - Working dataset by cohort (considering we are using period data)
## - Data from HMD dataset to working datasets
col <- 1
for(t in PeriodRange){
  row <- 1
  cohort_count <- 0
  for(x in AgeRange){
    d_xt[row,col] <- deaths[((deaths$Age==x) & (deaths$Year==(t + cohort_count))),"Male"]
    E_xt[row,col] <- exposures[((exposures$Age==x) & (exposures$Year==t + cohort_count)),"Male"] # -  & (exposures$Cohort==coh)
    row <- row + 1
    cohort_count <- cohort_count + 1
  }
  col <- col + 1
}

## - Mortality rates by age (x) and cohort (c)
m_xt <- d_xt / E_xt # - central death rates

for(i in 1:ncol(S_xt)){
  for(j in 1:nrow(S_xt)){
    S_xt[j,i] <- exp(-sum(m_xt[1:j,i]))
    mu_bar_USA[j,i] <- mean(m_xt[1:j,i])
  }
}

#========================== - Parameter estimation - =================



# - iterative full optimization function
## - Set starting values
x0_sv <- c(0, 0, 0)
delta_sv <- c(-1.3058e-06, -0.05220, -0.10132) # c(delta_11, delta_22, delta_33)
kappa_sv <- c(0.01163, 0.06787, 0.00506) # c(kappa_1, kappa_2, kappa_3)
sigma_sv <- c(0.00111, 0.00112, 0.00052) # c(sigma_11, sigma_22, sigma_33)
r_sv <- c(3.5544e-15, 0.54409, 1.797e-07) # c(r_1, r_2, r_c)

## - Fit
fit_BSi_3F <- it_f_opt_BSi(mu_bar=mu_bar_USA, x0=x0_sv, delta=delta_sv, kappa=kappa_sv, sigma=sigma_sv, r=r_sv, max_iter=3, tol_lik=10)

# - Improve estimation by local search
## - Initialize parameters for local search
x0_sv <- fit_BSi_3F$par_est$x0
delta_sv <- fit_BSi_3F$par_est$delta # c(delta_11, delta_22, delta_33)
kappa_sv <- fit_BSi_3F$par_est$kappa # c(kappa_1, kappa_2, kappa_3)
sigma_sv <- fit_BSi_3F$par_est$sigma # c(sigma_11, sigma_22, sigma_33)
r_sv <- c(fit_BSi_3F$par_est$r1, fit_BSi_3F$par_est$r2, fit_BSi_3F$par_est$rc) # c(r_1, r_2, r_c)

par_est_LS <- f_opt_BSi_LS(mu_bar=mu_bar_USA, x0=x0_sv, delta=delta_sv, kappa=kappa_sv, sigma=sigma_sv, r=r_sv)

# - Coordinate Ascent optimization function
## - Set starting values
x0_sv <- par_est_LS$par_est$x0
delta_sv <- par_est_LS$par_est$delta # c(delta_11, delta_22, delta_33)
kappa_sv <- par_est_LS$par_est$kappa # c(kappa_1, kappa_2, kappa_3)
sigma_sv <- par_est_LS$par_est$sigma # c(sigma_11, sigma_22, sigma_33)
r_sv <- c(par_est_LS$par_est$r1, par_est_LS$par_est$r2, par_est_LS$par_est$rc) # c(r_1, r_2, r_c)

## - Fit
fit_BSi_CA_3F <- co_asc_BSi(mu_bar=mu_bar_USA, x0=x0_sv, delta=delta_sv, kappa=kappa_sv, sigma=sigma_sv, r=r_sv, max_iter=100, tol_lik=0.1)


# - Final improvement of estimation by local search
## - Initialize parameters for local search
x0_sv <- fit_BSi_CA_3F$par_est$x0
delta_sv <- fit_BSi_CA_3F$par_est$delta # c(delta_11, delta_22, delta_33)
kappa_sv <- fit_BSi_CA_3F$par_est$kappa # c(kappa_1, kappa_2, kappa_3)
sigma_sv <- fit_BSi_CA_3F$par_est$sigma # c(sigma_11, sigma_22, sigma_33)
r_sv <- c(fit_BSi_CA_3F$par_est$r1, fit_BSi_CA_3F$par_est$r2, fit_BSi_CA_3F$par_est$rc) # c(r_1, r_2, r_c)

par_est_LS_final <- f_opt_BSi_LS(mu_bar=mu_bar_USA, x0=x0_sv, delta=delta_sv, kappa=kappa_sv, sigma=sigma_sv, r=r_sv)

## - Get parameter estimates
pe_BSi_3F <- par_est_LS_final$par_est

#========================== - Goodness of fit - =================

x0_est <- pe_BSi_3F$x0
delta_est <- pe_BSi_3F$delta
kappa_est <- pe_BSi_3F$kappa
sigma_est <- pe_BSi_3F$sigma
r1_est <- pe_BSi_3F$r1
r2_est <- pe_BSi_3F$r2
rc_est <- pe_BSi_3F$rc










