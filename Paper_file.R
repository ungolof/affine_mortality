#==================================================================================================
#============ - Sample code for implementation - ==================================================
#==================================================================================================


#======================== - Call files for estimation - ================================
source('Est_fun.R') # - Read all functions to optimize the models and the packages needed
source('Full_opt.R') # - Read all functions which perform the full optimization
source('Coordinate_Ascent.R') # - Read all functions which perform optimization by coordinate ascent


#======================== - Variable definition - ================================
AgeRange <- c(50:99) # - Set of ages under analysis
CohortRange <- seq(1795, 1914)  # - Set of cohorts under analysis
PeriodRange <- CohortRange + AgeRange[1] # - Period data

D_xt <- matrix(NA, length(AgeRange),length(CohortRange)) # - Deaths
E_xt <- matrix(NA, length(AgeRange),length(CohortRange)) # - Central Exposures

S_xt <- matrix(NA, length(AgeRange),length(CohortRange)) # - Survival function for each age and cohort
mu_bar <- matrix(NA, length(AgeRange),length(CohortRange)) # - mu_bar for each age (row) and cohort (column)

n_ages <- length(AgeRange) # - Number of ages (N)
n_years <- length(CohortRange) # - Number of cohort years (T)

colnames(m_xt) <- colnames(E_xt) <- colnames(D_xt) <- CohortRange
rownames(m_xt) <- rownames(E_xt) <- rownames(D_xt) <- AgeRange


#================== - Read file deaths and exposures directly from the HMD dataset=========
# - Read file deaths and exposures directly from the HMD dataset
deaths <- readHMDweb(CNTRY = "GBRCENW", item = "Deaths_1x1", username, password, fixup = TRUE)
exposures <- readHMDweb(CNTRY = "GBRCENW", item = "Exposures_1x1", username, password, fixup = TRUE)


# - Working dataset by cohort (considering we are using period data)
## - Data from HMD dataset to working datasets
col <- 1
for(t in PeriodRange){
  row <- 1
  cohort_count <- 0
  for(x in AgeRange){
    D_xt[row,col] <- deaths[((deaths$Age==x) & (deaths$Year==(t + cohort_count))),"Male"]
    E_xt[row,col] <- exposures[((exposures$Age==x) & (exposures$Year==t + cohort_count)),"Male"] # -  & (exposures$Cohort==coh)
    row <- row + 1
    cohort_count <- cohort_count + 1
  }
  col <- col + 1
}

## - Mortality rates by age (x) and cohort (c)
m_xt <- D_xt / E_xt # - central death rates

for(i in 1:ncol(S_xt)){
  for(j in 1:nrow(S_xt)){
    S_xt[j,i] <- exp(-sum(m_xt[1:j,i]))
    mu_bar[j,i] <- mean(m_xt[1:j,i])
  }
}

#========================== - Parameter estimation - =================

# - Set starting values
x0_sv <- c(0, 0, 0)
delta_sv <- c(-1.3058e-06, -0.05220, -0.10132) # c(delta_11, delta_22, delta_33)
kappa_sv <- c(0.01163, 0.06787, 0.00506) # c(kappa_1, kappa_2, kappa_3)
sigma_sv <- c(0.00111, 0.00112, 0.00052) # c(sigma_11, sigma_22, sigma_33)
r_sv <- c(3.5544e-15, 0.54409, 1.797e-07) # c(r_1, r_2, r_c)

n_iter <- 10 # - maximum number of iteration for the full optimization process (default setting in the it_f_opt_BSi)
tol_lik <- 10 # - minimum increase to go to next iterations
optim_m <- "Nelder-Mead"  # - Choice of the optimization method within the function optim

# - iterative full optimization function
fit_BSi_3F_full_opt <- it_f_opt_BSi()

## - Get parameter estimates
par_est_BSi_3F <- fit_BSi_3F_full_opt$par_est
## - Negative log-likelihood value (without constants)
neg_log_lik_v <- fit_BSi_3F_full_opt$log_lik

# - Improvement with Local search from BFGS
x0_sv <- fit_BSi_3F_full_opt$par_est$x0
delta_sv <- fit_BSi_3F_full_opt$par_est$delta
kappa_sv <- fit_BSi_3F_full_opt$par_est$kappa
sigma_sv <-fit_BSi_3F_full_opt$par_est$sigma
r_sv <- c(fit_BSi_3F_full_opt$par_est$r1, fit_BSi_3F_full_opt$par_est$r2, fit_BSi_3F_full_opt$par_est$rc)
  
par_est_LS <- f_opt_BSi_LS(x0_sv, delta_sv, kappa_sv, sigma_sv, r_sv)

x0_sv <- par_est_LS$par_est$x0
delta_sv <- par_est_LS$par_est$delta
kappa_sv <- par_est_LS$par_est$kappa
sigma_sv <-par_est_LS$par_est$sigma
r_sv <- c(par_est_LS$par_est$r1, par_est_LS$par_est$r2, par_est_LS$par_est$rc)


fit_BSi_3F_CA <- co_asc_BSi(x0_sv, delta_sv, kappa_sv, sigma_sv, r_sv, 10, 0.1)

x0_sv <- fit_BSi_3F_CA$par_est$x0
delta_sv <- fit_BSi_3F_CA$par_est$delta
kappa_sv <- fit_BSi_3F_CA$par_est$kappa
sigma_sv <-fit_BSi_3F_CA$par_est$sigma
r_sv <- c(fit_BSi_3F_CA$par_est$r1, fit_BSi_3F_CA$par_est$r2, fit_BSi_3F_CA$par_est$rc)

par_est_LS <- f_opt_BSi_LS(x0_sv, delta_sv, kappa_sv, sigma_sv, r_sv)

# - AFNS dependent factor model
fit_AFNSd_full_opt <- it_f_opt_AFNSd()







