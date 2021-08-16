#============ - Sample code for implementation - ==================================================
# - We show the code for analysing the Blackburn-Sherris model with three independent factors
# - for the USA mortality dataset
# - The analysis of other models follow along similar lines
#==================================================================================================

#=============================== - Call files for estimation - ====================================
source('Est_fun.R') # - Read all functions to optimize the models and the packages needed
source('Full_opt.R') # - Read all functions which perform the full optimization
source('Coordinate_Ascent.R') # - Read all functions which perform optimization by coordinate ascent
source('GoodnessofFit.R') # - Read functions to compute the Goodness of fit measures and of residuals
source('FilterSmoother.R') # - Read functions to perform filtering and smoothing of the latent state variables
source('StErr_Bootstrap.R') # - Read functions to perform bootstrap estimation of the covariance of the parameter estimates
source('StErr_MultImp.R') # - Read functions to estimate the covariance of the parameter estimates by Multiple Imputation
source('Projection.R') # - Read functions yielding the projected cohort survival curves



#=============================== - Variable definition - ==========================================
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

rownames(mu_bar_USA) <- AgeRange
colnames(mu_bar_USA) <- CohortRange


#================== - Read file deaths and exposures directly from the HMD dataset - ==============
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

#=================================== - Parameter estimation - ===============================



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
log_lik_BSi_3F <- par_est_LS_final$log_lik

#====================================== - Goodness of fit - ==============================

# - Get estimated parameters
x0_est <- pe_BSi_3F$x0
delta_est <- pe_BSi_3F$delta
kappa_est <- pe_BSi_3F$kappa
sigma_est <- pe_BSi_3F$sigma
r1_est <- pe_BSi_3F$r1
r2_est <- pe_BSi_3F$r2
rc_est <- pe_BSi_3F$rc

# - Get fitted average force of mortality

mu_bar_hat_BSi_3F <- mu_bar_hat_BSi(x0=x0_est, delta=delta_est, kappa=kappa_est, sigma=sigma_est, r=c(r1_est, r2_est, rc_est), mu_bar=mu_bar_USA)

# - Get fitted force of mortality
mu_hat_BSi_3F <- avg2rates(mu_bar_hat_BSi_3F)

# - AIC
AIC_BSi_3F <- AIC_BIC(log_lik_BSi_3F, 12, nrow(mu_bar_USA)*ncol(mu_bar_USA))$AIC
# - BIC
BIC_BSi_3F <- AIC_BIC(log_lik_BSi_3F, 12, nrow(mu_bar_USA)*ncol(mu_bar_USA))$BIC

# - RMSE
RMSE_BSi_3F <- RMSE(mu_bar_USA, mu_bar_hat_BSi_3F)

# - MAPE (based on Survival probabilities)
## - Step 1 -  Get fitted survival probabilities
S_xt_BSi_3F <- matrix(NA, length(AgeRange),length(CohortRange))

for(i in 1:ncol(S_xt)){
  for(j in 1:nrow(S_xt)){
    S_xt_BSi_3F[j,i] <- exp(-mu_bar_hat_BSi_3F[j,i] * j)
  }
}

## - Step 2 - Get MAPE (for each age)
MAPE_BSi_3F <- MAPE_row(S_xt, S_xt_BSi_3F)

## - Step 3 - Plot

plot(c(50:85), MAPE_BSi_3F[1:36]*100, type="o", pch=1, cex=0.5, ylim=c(0, max(MAPE_BSi_3F[1:36]*100)), xlim=c(50,100), ylab="Percentage error", xlab="Age", lty=1)
par(new = TRUE)
plot(c(50:99), MAPE_BSi_3F[1:50]*100, type="o", pch=1, cex=0.5, lty=1, axes = FALSE, bty = "n", xlab = "", ylab = "", ylim=c(0,max(MAPE_BSi_3F[36:50]*100)),col='white')
lines(c(85:99), MAPE_BSi_3F[36:50]*100, type="o", pch=1, cex=0.5, lty=1)
abline(v=84.5)
axis(side=4)

#========================= - Residuals graphical checks - =============================

# - Libraries for heatmaps
library(heatmaply)
library(fields)

# - Residuals
res_BSi_3F <- residuals_f(mu_bar_USA, mu_bar_hat_BSi_3F)

## - 3D Plotting
persp3D(CohortRange, AgeRange, t(res_BSi), main="", zlab="", xlab="Cohort", ylab="Age", theta = 35, phi = 25, shade = 0.5, 
        axes = T, box=TRUE, nticks=5, ticktype="detailed", xlim=c(1880,1920), ylim = c(50,100), zlim=c(-0.001,0.01))

# - 0-1 Residuals
res01_BSi_3F <- residuals_01(mu_bar_USA, mu_bar_hat_BSi_3F)

## - Heatmap plotting example 1 (use package heatmaply)
heatmaply(res01_BSi_3F[c(50:1),], dendrogram = "none", xlab = "Cohort", ylab="Age", xaxis_font_size="6px", yaxis_font_size="6px")


# - Standardized residuals
res_std_BSi_3F <- residuals_std_BSi(observed=mu_bar_USA, x0=x0_est, delta=delta_est, kappa=kappa_est, sigma=sigma_est, r=c(r1_est, r2_est, rc_est))

## - Heatmap plotting example 2 (use package fields)

### - Choice of the color
col <- colorRampPalette(RColorBrewer::brewer.pal(9, "Greys"))(64)
col <- colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(64)

maxRes <- max(abs(res_std_BSi_3F), na.rm = TRUE)
reslim <- c(-maxRes, maxRes)

image.plot(CohortRange, AgeRange, t(res_std_BSi_3F), zlim = reslim, ylab = "Age", xlab = "Cohort", col = col)


#========================= - Projection - =============================

# - 10 year projection
Proj_S_BSi_3F <- S_t_BSi_proj(x0=x0_est, delta=delta_est, kappa=kappa_est, sigma=sigma_est, r=c(r1_est, r2_est, rc_est), mu_bar=mu_bar_USA, proj_years = 10)

## - Plotting example
plot(c(51:100), Proj_S_BSi_3F, type="l", xlab="Age", ylab="S(t)", lwd=1)


#======================= - Estimation of the standard errors of the parameter estimates

Par_unc_BSi_3F <- CovEst_BS_BSi(x0=x0_est, delta=delta_est, kappa=kappa_est, sigma=sigma_est, r=c(r1_est, r2_est, rc_est), mu_bar=mu_bar_USA, n_BS=100, t_ex=4)









