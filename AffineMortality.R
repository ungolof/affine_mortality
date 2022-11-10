# - Functions for AffineMortality package

#source("Coordinate_Ascent.R")
#source("Est_fun.R")
#source("FilterSmoother.R")
#source("GoodnessofFit.R")
#source("Projection.R")
#source("StErr_Bootstrap.R")
#source("StErr_MultImp.R")


## - Fit

### - default starting values
sv_default <- list(BSi=list(x0=c(6.960591e-03, 9.017154e-03, 5.091784e-03, 0), delta=c(-1.305830e-06, -5.220474e-02, -1.013210e-01, -0.02), kappa=c(1.162624e-02, 6.787268e-02, 5.061539e-03, 0.0001), sigma=exp(c(-6.806310, -6.790270, -7.559145, -8)), r1=exp(-3.327060e+01), r2=exp(-6.086479e-01), rc=exp(-1.553156e+01)),
                   BSd=list(x0=c(2.191140e-03, -8.855686e-03, 2.711990e-02), delta=c(-5.175933e-02, 4.578315e-01, -5.175921e-02, -2.299199e-01, 1.383445e-02, -6.310253e-02), kappa=c(3.455255e-02, 1.075876e-02, 1.000030e-02), sigma_dg=c(6.789215e-04, 2.036748e-03, 1.875928e-03), Sigma_cov=c(-1.260778e-06, 1.194974e-06, -3.718267e-06), r1=exp(-3.345631e+01), r2=exp(-6.015438e-01), rc=exp(-1.557244e+01)),
                   AFNSi=list(x0=c(1.091714e-02, 1.002960e-02, -5.990785e-04), delta=-8.304334e-02, kappa=c(9.154603e-03, 1.067658e-02, 7.439715e-03), sigma=exp(c(-7.318991, -7.535594, -8.456025)), r1=exp(-3.371775e+01), r2=exp(-5.887962e-01), rc=exp(-1.548729e+01)),
                   AFNSd=list(x0=c(9.582516e-03, 1.094110e-02, -1.503155e-03), delta=-7.487697e-02, kappa=c(1.389363e-02, 3.525542e-03, 3.004883e-03), sigma_dg=c(3.215422e-03, 2.625474e-03, 1.164715e-03), Sigma_cov=c(-8.328978e-06, -3.685028e-06, 3.036376e-06), r1=exp(-3.335725e+01), r2=exp(-6.066149e-01), rc=exp(-1.552061e+01)),
                   AFGNSi=list(x0=c(1.091714e-02, 1.002960e-02, 5.990785e-04, 0, 0), delta=c(8.304334e-02, delta2=-0.05), kappa=c(9.154603e-03, 1.067658e-02, 7.439715e-03, 0.003, 0.001), sigma=exp(c(-7.318991, -7.535594, -8.456025, -7, -7)), r1=exp(-3.371775e+01), r2=exp(-5.887962e-01), rc=exp(-1.548729e+01)),
                   AFGNSd=list(x0=c(1.091714e-02, 1.002960e-02, 5.990785e-04, 0, 0), delta=c(8.304334e-02, delta2=-0.05), kappa=c(9.154603e-03, 1.067658e-02, 7.439715e-03, 0.003, 0.001), sigma_dg=c(3.215422e-03, 2.625474e-03, 1.164715e-03, 0.0003, 0.0001), Sigma_cov=rep(0, 10), r1=exp(-3.335725e+01), r2=exp(-6.066149e-01), rc=exp(-1.552061e+01)),
                   CIR=list(x0=c(5.080033e-11, 1.535266e-02, 2.972783e-03, 5.091784e-03), delta=c(-0.2183696, 0.2865177, -0.1307629, -0.06), kappa=c(0.001292314, 0.485964737, 0.129874648, 0.001), sigma=c(0.002777236, 0.004947513, 0.021948997, 0.0136), theta_P=c(5.368400e-03, 7.074616e-03, 5.759856e-09, 0.00137), r1=2.776231e-22, r2=8.408350e-01, rc=1.668493e-07))




# - Starting values must be provided by the user in the format of a list, eg.: list(x0=c(...))
## - eventually add the option for using subplex optimization
#' @title affine_fit
#'
#' @description Estimation of affine mortality models
#'
#' @param model Specific model to be fit. E.g. for the Blackburn-Sherris model we have model="BS", and so on
#' @param fact_dep Boolean parameter indicating whether estimate models with factor dependence (fact_dep=TRUE) or independence
#' @param n_factors Number of factors. For some models, these are set by default
#' @param data Table with the average mortality rates
#' @param st_val Starting value for the parameters. If not set, then default values will be used
#' @param max_iter Maximum number of iterations for the Coordinate Ascent estimation algorithm in absence of convergence
#' @param tolerance Minimum value of improvement of the log-likelihood value before convergence
#' @param wd Working directory for saving the partial output of the estimation process
#'
#' @return A data frame object that contains a summary of a sample that
#'     can later be converted to a TeX output using \code{overview_print}
#' @examples
#' data(toydata)
#' output_table <- overview_tab(dat = toydata, id = ccode, time = year)
#' @export
affine_fit <- function(model="BS", fact_dep=FALSE, n_factors=3, data=data_default, st_val=0, max_iter=200, tolerance=0.1, wd=0){

  if(model=="AFNS"){
    if(fact_dep==TRUE){
      if(st_val==0){ # - then use default starting value
        st_val <- sv_default$AFNSd
      }
      fit <- co_asc_AFNSd(mu_bar = data, x0=st_val$x0, delta=st_val$delta, kappa=st_val$kappa, sigma_dg = st_val$sigma_dg, Sigma_cov = st_val$Sigma_cov, r=c(st_val$r1, st_val$r2, st_val$rc), max_iter=max_iter, tol_lik=tolerance, workdir=wd)
      return(list(model=model, fit=fit, n.parameters=16, AIC=AIC_BIC(fit$log_lik, 16, (nrow(data) * ncol(data)))$AIC, BIC=AIC_BIC(fit$log_lik, 16, (nrow(data) * ncol(data)))$BIC))
    }
    else{
      if(st_val==0){ # - then use default starting value
        st_val <- sv_default$AFNSi
      }
      fit <- co_asc_AFNSi(mu_bar = data, x0=st_val$x0, delta=st_val$delta, kappa=st_val$kappa, sigma = st_val$sigma, r=c(st_val$r1, st_val$r2, st_val$rc), max_iter=max_iter, tol_lik=tolerance, workdir=wd)
      return(list(model=model, fit=fit, n.parameters=13, AIC=AIC_BIC(fit$log_lik, 13, (nrow(data) * ncol(data)))$AIC, BIC=AIC_BIC(fit$log_lik, 13, (nrow(data) * ncol(data)))$BIC))
    }
  }else{
    if(model=="BS"){
      if(fact_dep==TRUE){
          if(n_factors==2){
            if(st_val==0){
              st_val <- list(x0=sv_default$BSd$x0[1:2], delta=sv_default$BSd$delta[1:3], kappa=sv_default$BSd$kappa[1:2], sigma_dg=sv_default$BSd$sigma_dg[1:2], Sigma_cov=sv_default$BSd$Sigma_cov[1], r=c(sv_default$BSd$r1, sv_default$BSd$r2, sv_default$BSd$rc)) # - change in order to take into account the number of factors
            }
          fit <- co_asc_BSd_2F(mu_bar = data, x0=st_val$x0, delta=st_val$delta, kappa=st_val$kappa, sigma_dg = st_val$sigma_dg, Sigma_cov = st_val$Sigma_cov, r=c(st_val$r1, st_val$r2, st_val$rc), max_iter=max_iter, tol_lik=tolerance, workdir=wd)
          return(list(model=model, fit=fit, n.parameters=13, AIC=AIC_BIC(fit$log_lik, 13, (nrow(data) * ncol(data)))$AIC, BIC=AIC_BIC(fit$log_lik, 13, (nrow(data) * ncol(data)))$BIC))
          } else{
            if(st_val==0){
              st_val <- sv_default$BSd # - change in order to take into account the number of factors
            }
            fit <- co_asc_BSd_3F(mu_bar = data, x0=st_val$x0, delta=st_val$delta, kappa=st_val$kappa, sigma_dg = st_val$sigma_dg, Sigma_cov = st_val$Sigma_cov, r=c(st_val$r1, st_val$r2, st_val$rc), max_iter=max_iter, tol_lik=tolerance, workdir=wd)
            return(list(model=model, fit=fit, n.parameters=21, AIC=AIC_BIC(fit$log_lik, 21, (nrow(data) * ncol(data)))$AIC, BIC=AIC_BIC(fit$log_lik, 21, (nrow(data) * ncol(data)))$BIC))
          }
        }else{ #i.e. factor independence
            if(st_val==0){
              st_val <- sv_default$BSi
            }
            fit <- co_asc_BSi(mu_bar = data, x0=st_val$x0[1:n_factors], delta=st_val$delta[1:n_factors], kappa=st_val$kappa[1:n_factors], sigma = st_val$sigma[1:n_factors], r=c(st_val$r1, st_val$r2, st_val$rc), max_iter=max_iter, tol_lik=tolerance, workdir=wd)
            return(list(model=model, fit=fit, n.parameters=(3 + n_factors*4), AIC=AIC_BIC(fit$log_lik, (3 + n_factors*4), (nrow(data) * ncol(data)))$AIC, BIC=AIC_BIC(fit$log_lik, (3 + n_factors*4), (nrow(data) * ncol(data)))$BIC))
        }}else{
       if(model=="CIR"){
         if(st_val==0){
           st_val <- sv_default$CIR
         }
         fit <- co_asc_CIR(mu_bar = mu_bar, x0=st_val$x0, delta=st_val$delta, kappa=st_val$kappa, sigma=st_val$sigma, theta_P=st_val$theta_P, r=c(st_val$r1, st_val$r2, st_val$rc), max_iter=max_iter, tol_lik=tolerance, workdir=wd)
         return(list(model=model, fit=fit, n.parameters=(3 + n_factors*5), AIC=AIC_BIC(fit$log_lik, (3 + n_factors*5), (nrow(data) * ncol(data)))$AIC, BIC=AIC_BIC(fit$log_lik, (3 + n_factors*5), (nrow(data) * ncol(data)))$BIC))
       }else{
         if(model=="AFUNS"){
           if(fact_dep==TRUE){
             if(st_val==0){
             st_val <- sv_default$AFUNSd
             }
             fit <- co_asc_AFUNSd(mu_bar = data, x0=st_val$x0, delta=st_val$delta, kappa=st_val$kappa, sigma_dg = st_val$sigma_dg, Sigma_cov = st_val$Sigma_cov, r=c(st_val$r1, st_val$r2, st_val$rc), max_iter=max_iter, tol_lik=tolerance, workdir=wd)
             return(list(model=model, fit=fit, n.parameters=18, AIC=AIC_BIC(fit$log_lik, 18, (nrow(data) * ncol(data)))$AIC, BIC=AIC_BIC(fit$log_lik, 18, (nrow(data) * ncol(data)))$BIC))
           }else{
             if(st_val==0){
               st_val <- sv_default$AFUNSi
             }
             fit <- co_asc_AFUNSi(mu_bar = data, x0=st_val$x0, delta=st_val$delta, kappa=st_val$kappa, sigma = st_val$sigma, r=c(st_val$r1, st_val$r2, st_val$rc), max_iter=max_iter, tol_lik=tolerance, workdir=wd)
             return(list(model=model, fit=fit, n.parameters=15, AIC=AIC_BIC(fit$log_lik, 15, (nrow(data) * ncol(data)))$AIC, BIC=AIC_BIC(fit$log_lik, 15, (nrow(data) * ncol(data)))$BIC))
           }}else{
             if(model=="AFRNS"){
               if(fact_dep==TRUE){
                 if(st_val==0){
                   st_val <- sv_default$AFRNSd
                 }
                 fit <- co_asc_AFRNSd(mu_bar = data, x0=st_val$x0, delta=st_val$delta, kappa=st_val$kappa, sigma_dg = st_val$sigma_dg, Sigma_cov = st_val$Sigma_cov, r=c(st_val$r1, st_val$r2, st_val$rc), max_iter=max_iter, tol_lik=tolerance, workdir=wd)
                 return(list(model=model, fit=fit, n.parameters=13, AIC=AIC_BIC(fit$log_lik, 13, (nrow(data) * ncol(data)))$AIC, BIC=AIC_BIC(fit$log_lik, 13, (nrow(data) * ncol(data)))$BIC))
               }else{
                 if(st_val==0){
                   st_val <- sv_default$AFRNSi
                 }
                 fit <- co_asc_AFRNSi(mu_bar = data, x0=st_val$x0, delta=st_val$delta, kappa=st_val$kappa, sigma = st_val$sigma, r=c(st_val$r1, st_val$r2, st_val$rc), max_iter=max_iter, tol_lik=tolerance, workdir=wd)
                 return(list(model=model, fit=fit, n.parameters=10, AIC=AIC_BIC(fit$log_lik, 10, (nrow(data) * ncol(data)))$AIC, BIC=AIC_BIC(fit$log_lik, 10, (nrow(data) * ncol(data)))$BIC))
               }
             }else{
             if(model=="GMk"){
               if(fact_dep==TRUE){
                 if(st_val==0){
                   st_val <- sv_default$GMkd
                 }
                 fit <- co_asc_GMkd(mu_bar = data, x0=st_val$x0, delta=st_val$delta, kappa=st_val$kappa, gamma=st_val$gamma, sigma_dg = st_val$sigma_dg, Sigma_cov = st_val$Sigma_cov, r=c(st_val$r1, st_val$r2, st_val$rc), max_iter=max_iter, tol_lik=tolerance, workdir=wd)
                 return(list(model=model, fit=fit, n.parameters=13, AIC=AIC_BIC(fit$log_lik, 13, (nrow(data) * ncol(data)))$AIC, BIC=AIC_BIC(fit$log_lik, 13, (nrow(data) * ncol(data)))$BIC))
               }else{
                 if(st_val==0){
                   st_val <- sv_default$GMki
                 }
                 fit <- co_asc_GMki(mu_bar = data, x0=st_val$x0, delta=st_val$delta, kappa=st_val$kappa, sigma = st_val$sigma, r=c(st_val$r1, st_val$r2, st_val$rc), max_iter=max_iter, tol_lik=tolerance, workdir=wd)
                 return(list(model=model, fit=fit, n.parameters=11, AIC=AIC_BIC(fit$log_lik, 11, (nrow(data) * ncol(data)))$AIC, BIC=AIC_BIC(fit$log_lik, 11, (nrow(data) * ncol(data)))$BIC))
               }
             }
           }


         }
         }
   }}}


#======================== - Filtering distribution - ===================================

# - Filtering distribution (xfilter)
#' @title xfilter
#'
#' @description Estimation of the moments of the latent states
#'
#' @param model Affine model. E.g. for the Blackburn-Sherris model we have model="BS", and so on
#' @param fact_dep Boolean parameter indicating whether estimate models with factor dependence (fact_dep=TRUE) or independence
#' @param n_factors Number of factors. For some models, these are set by default
#' @param parameters Starting value for the parameters. If not set, then default values will be used
#' @param data Table with the average mortality rates
#'
#' @return A data frame object that contains a summary of a sample that
#'     can later be converted to a TeX output using \code{overview_print}
#' @examples
#' data(toydata)
#' output_table <- overview_tab(dat = toydata, id = ccode, time = year)
#' @export
xfilter <- function(model="BS", fact_dep=FALSE, n_factors=3, parameters=0, data=data_default){
  if(model=="AFNS"){
    if(fact_dep==TRUE){
      if(parameters==0){ # - then use default starting value
        parameters <- sv_default$AFNSd
      }
      filter <- KF_AFNSd_uKD(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
      return(filter)
    } else{
      if(parameters==0){
        parameters <- sv_default$AFNSi
      }
      filter <- KF_AFNSi_uKD(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
      return(filter)
    }
  } else{
    if(model=="BS"){
      if(fact_dep==TRUE){
        if(parameters==0){
          if(n_factors==2){
            parameters <- list(x0=sv_default$BSd$x0[1:2], delta=sv_default$BSd$delta[1:3], kappa=sv_default$BSd$kappa[1:2], sigma_dg=sv_default$BSd$sigma_dg[1:2], Sigma_cov=sv_default$BSd$Sigma_cov[1], r=c(sv_default$BSd$r1, sv_default$BSd$r2, sv_default$BSd$rc)) # - change in order to take into account the number of factors
          }
          filter <- KF_BSd_2F_uKD(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
          return(filter)

        } else{
          if(n_factors==2){
            filter <- KF_BSd_2F_uKD(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
            return(filter)
          } else{
            filter <- KF_BSd_3F_uKD(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
            return(filter)
          }
        }
      } else{ # - it is the Blackburn-Sherris model with independent factors
        if(parameters==0){
          parameters <- sv_default$BSi
        }
        filter <- KF_BSi_uKD(x0=parameters$x0[1:n_factors], delta=parameters$delta[1:n_factors], kappa=parameters$kappa[1:n_factors], sigma = parameters$sigma[1:n_factors], r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
        return(filter)

      }
    } else{
      if(model=="CIR"){
        if(parameters==0){
          parameters <- sv_default$CIR
        }
        filter <- KF_CIR_uKD(x0=parameters$x0[1:n_factors], delta=parameters$delta[1:n_factors], kappa=parameters$kappa[1:n_factors], sigma=parameters$sigma[1:n_factors], theta_P=parameters$theta_P[1:n_factors], r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
        return(filter)
      } else{
        if(model=="AFGNS"){# - if none of the previous model was selected, then it is an AFGNS
        if(fact_dep==TRUE){
          if(parameters==0){
            parameters <- sv_default$AFGNSd
          }
          filter <- KF_AFGNSd_uKD(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
          return(filter)
        } else{
          if(parameters==0){
            parameters <- sv_default$AFGNSi
          }
          filter <- KF_AFGNSi_uKD(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
          return(filter)
        }} else{
          if(model=="AFUNS"){
            if(fact_dep==TRUE){
              if(parameters==0){
                parameters <- sv_default$AFUNSd
              }
              filter <- KF_AFUNSd_uKD(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
              return(filter)
            } else{
                if(parameters==0){
                  parameters <- sv_default$AFUNSi
                }
                filter <- KF_AFUNSi_uKD(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
                return(filter)
              }
          }else{
            if(model=="AFRNS"){
              if(fact_dep==TRUE){
                if(parameters==0){
                  parameters <- sv_default$AFRNSd
                }
                filter <- KF_AFRNSd_uKD(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
                return(filter)
              } else{
                if(parameters==0){
                  parameters <- sv_default$AFRNSi
                }
                filter <- KF_AFRNSi_uKD(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
                return(filter)
              }
            } else{
              if(model=="GMk"){
                if(fact_dep==TRUE){
                  if(parameters==0){
                    parameters <- sv_default$GMkd
                  }
                  filter <- KF_GMkd_uKD(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
                  return(filter)
                } else{
                  if(parameters==0){
                    parameters <- sv_default$GMki
                  }
                  filter <- KF_GMki_uKD(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
                  return(filter)
                }
              }
            }
          }
        }
      }
    }
  }
}

#======================== - Smoothing - ===================================
#' @title xsmooth
#'
#' @description Estimation of the moments of the latent states
#'
#' @param filterobject An object as of the output of xfilter
#' @param kappa parameter kappa from the real-world dynamics of the latent variables
#'
#' @return A data frame object that contains a summary of a sample that
#'     can later be converted to a TeX output using \code{overview_print}
#' @examples
#' data(toydata)
#' output_table <- overview_tab(dat = toydata, id = ccode, time = year)
#' @export
xsmooth <- function(filterobject, kappa){
  smooth <- RTS_sm_bas(filterobject$X_t, filterobject$X_t_c, filterobject$S_t, filterobject$S_t_c, kappa, (ncol(filterobject$X_t)-1))
  return(smooth)
}

#======================== - Goodness of fit - ===================================

## - Fitted rates
#' @title mubar_hat
#'
#' @description Function returning the fitted values of the average mortality rates
#'
#' @param model Affine model. E.g. for the Blackburn-Sherris model we have model="BS", and so on
#' @param fact_dep Boolean parameter indicating whether estimate models with factor dependence (fact_dep=TRUE) or independence
#' @param n_factors Number of factors. For some models, these are set by default
#' @param parameters Starting value for the parameters. If not set, then default values will be used
#' @param data Table with the average mortality rates
#'
#' @return A data frame object that contains a summary of a sample that
#'     can later be converted to a TeX output using \code{overview_print}
#' @examples
#' data(toydata)
#' output_table <- overview_tab(dat = toydata, id = ccode, time = year)
#' @export
mubar_hat <- function(model="BS", fact_dep=FALSE, n_factors=3, parameters=0, data=data_default){
  if(model=="AFNS"){
    if(fact_dep==TRUE){
      if(parameters==0){ # - then use default starting value
        parameters <- sv_default$AFNSd
      }
      fitted <- mu_bar_hat_AFNSd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
      return(fitted)
    } else{
      if(parameters==0){
        parameters <- sv_default$AFNSi
      }
      fitted <- mu_bar_hat_AFNSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
      return(fitted)
    }
  } else{
    if(model=="BS"){
      if(fact_dep==TRUE){
        if(parameters==0){
          if(n_factors==2){
            parameters <- list(x0=sv_default$BSd$x0[1:2], delta=sv_default$BSd$delta[1:3], kappa=sv_default$BSd$kappa[1:2], sigma_dg=sv_default$BSd$sigma_dg[1:2], Sigma_cov=sv_default$BSd$Sigma_cov[1], r=c(sv_default$BSd$r1, sv_default$BSd$r2, sv_default$BSd$rc)) # - change in order to take into account the number of factors
            fitted <- mu_bar_hat_BSd_2F(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
            return(fitted)
          } else{ # - the factors must be three

            parameters <- sv_default$BSd # - change in order to take into account the number of factors
            fitted <- mu_bar_hat_BSd_3F(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
            return(fitted)

          }
        } else{
          if(n_factors==2){
            fitted <- mu_bar_hat_BSd_2F(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
            return(fitted)
          } else{
            fitted <- mu_bar_hat_BSd_3F(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
            return(fitted)
          }
        }
      } else{ # - it is the Blackburn-Sherris model with independent factors
        if(parameters==0){
          parameters <- sv_default$BSi
        }
        fitted <- mu_bar_hat_BSi(x0=parameters$x0[1:n_factors], delta=parameters$delta[1:n_factors], kappa=parameters$kappa[1:n_factors], sigma = parameters$sigma[1:n_factors], r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
        return(fitted)
      }
    } else{
      if(model=="CIR"){
        if(parameters==0){
          parameters <- sv_default$CIR
        }
        fitted <- mu_bar_hat_CIR(x0=parameters$x0[1:n_factors], delta=parameters$delta[1:n_factors], kappa=parameters$kappa[1:n_factors], sigma=parameters$sigma[1:n_factors], theta_P=parameters$theta_P[1:n_factors], r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
        return(fitted)
        } else{
          if(model=="AFGNS"){    # - if none of the previous model was selected, then it is an AFGNS
        if(fact_dep==TRUE){
          if(parameters==0){
            parameters <- sv_default$AFGNSd
          }
          fitted <- mu_bar_hat_AFGNSd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
          return(fitted)
          } else{
          if(parameters==0){
            parameters <- sv_default$AFGNSi
          }
            fitted <- mu_bar_hat_AFGNSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
            return(fitted)
          }}else{
            if(model=="AFUNS"){
              if(fact_dep==TRUE){
                if(parameters==0){
                  parameters <- sv_default$AFUNSd
                }
                fitted <- mu_bar_hat_AFUNSd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
                return(fitted)
              } else{
                if(parameters==0){
                  parameters <- sv_default$AFGNSi
                }
                fitted <- mu_bar_hat_AFUNSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
                return(fitted)
              }
            }else{
              if(model=="AFRNS"){
                if(fact_dep==TRUE){
                  if(parameters==0){
                    parameters <- sv_default$AFRNSd
                  }
                  fitted <- mu_bar_hat_AFRNSd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
                  return(fitted)
                } else{
                  if(parameters==0){
                    parameters <- sv_default$AFRNSi
                  }
                  fitted <- mu_bar_hat_AFRNSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
                  return(fitted)
                }
              }else{
                if(model=="GMk"){
                  if(fact_dep==TRUE){
                    if(parameters==0){
                      parameters <- sv_default$GMkd
                    }
                    fitted <- mu_bar_hat_GMkd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
                    return(fitted)
                  } else{
                    if(parameters==0){
                      parameters <- sv_default$GMki
                    }
                    fitted <- mu_bar_hat_GMki(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
                    return(fitted)
                  }
                }
              }
            }
          }
        }
      }
    }
  }

## - Fitted rates
#' @title std_res
#'
#' @description Function returning the standardized residuals
#'
#' @param model Affine model. E.g. for the Blackburn-Sherris model we have model="BS", and so on
#' @param fact_dep Boolean parameter indicating whether estimate models with factor dependence (fact_dep=TRUE) or independence
#' @param n_factors Number of factors. For some models, these are set by default
#' @param parameters Starting value for the parameters. If not set, then default values will be used
#' @param data Table with the average mortality rates
#'
#' @return A data frame object that contains a summary of a sample that
#'     can later be converted to a TeX output using \code{overview_print}
#' @examples
#' data(toydata)
#' output_table <- overview_tab(dat = toydata, id = ccode, time = year)
#' @export
std_res <- function(model="BS", fact_dep=FALSE, n_factors=3, parameters=0, data=data_default){
  if(model=="AFNS"){
    if(fact_dep==TRUE){
      if(parameters==0){ # - then use default starting value
        parameters <- sv_default$AFNSd
      }
    } else{
      if(parameters==0){
        parameters <- sv_default$AFNSi
      }
        std_res <- residuals_std_AFNSi(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc))
        return(std_res)
    }
  } else{
    if(model=="BS"){
      if(fact_dep==TRUE){
        if(parameters==0){
          if(n_factors==2){
            parameters <- list(x0=sv_default$BSd$x0[1:2], delta=sv_default$BSd$delta[1:3], kappa=sv_default$BSd$kappa[1:2], sigma_dg=sv_default$BSd$sigma_dg[1:2], Sigma_cov=sv_default$BSd$Sigma_cov[1], r=c(sv_default$BSd$r1, sv_default$BSd$r2, sv_default$BSd$rc)) # - change in order to take into account the number of factors
            std_res <- residuals_std_BSd_2F(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc))
            return(std_res)
          } else{ # - the factors must be three
            parameters <- sv_default$BSd # - change in order to take into account the number of factors
            std_res <- residuals_std_BSd_3F(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc))
            return(std_res)

          }
        } else{
          if(n_factors==2){
            std_res <- residuals_std_BSd_2F(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc))
            return(std_res)
          } else{
            std_res <- residuals_std_BSd_3F(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc))
            return(std_res)
          }
        }
      } else{ # - it is the Blackburn-Sherris model with independent factors
        if(parameters==0){
          parameters <- sv_default$BSi
          std_res <- residuals_std_BSi(data, x0=parameters$x0[1:n_factors], delta=parameters$delta[1:n_factors], kappa=parameters$kappa[1:n_factors], sigma = parameters$sigma[1:n_factors], r=c(parameters$r1, parameters$r2, parameters$rc))
          return(std_res)
        } else{
          std_res <- residuals_std_BSi(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc))
          return(std_res)
        }
      }
    } else{
      if(model=="CIR"){
        if(parameters==0){
          parameters <- sv_default$CIR
        }
          std_res <- residuals_std_CIR(data, x0=parameters$x0[1:n_factors], delta=parameters$delta[1:n_factors], kappa=parameters$kappa[1:n_factors], sigma=parameters$sigma[1:n_factors], theta_P=parameters$theta_P[1:n_factors], r=c(parameters$r1, parameters$r2, parameters$rc))
          return(std_res)
      } else{
        if(model=="AFGNS"){
        if(fact_dep==TRUE){
          if(parameters==0){
            parameters <- sv_default$AFGNSd
          }
            std_res <- residuals_std_AFGNSd(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc))
            return(std_res)
        } else{
          if(parameters==0){
            parameters <- sv_default$AFGNSi
          }
            std_res <- residuals_std_AFGNSi(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc))
            return(std_res)
        }}else{
          if(model=="AFRNS"){
            if(fact_dep==TRUE){
              if(parameters==0){
                parameters <- sv_default$AFRNSd
              }
              std_res <- residuals_std_AFRNSd(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc))
              return(std_res)
            } else{
              if(parameters==0){
                parameters <- sv_default$AFRNSi
              }
              std_res <- residuals_std_AFRNSi(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc))
              return(std_res)
            }
          }else{
            if(model=="AFUNS"){
              if(fact_dep==TRUE){
                if(parameters==0){
                  parameters <- sv_default$AFUNSd
                }
                std_res <- residuals_std_AFUNSd(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc))
                return(std_res)
              } else{
                if(parameters==0){
                  parameters <- sv_default$AFUNSi
                }
                std_res <- residuals_std_AFUNSi(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc))
                return(std_res)
              }
            }else{
              if(model=="GMk"){
              if(fact_dep==TRUE){
                if(parameters==0){
                  parameters <- sv_default$GMkd
                }
                std_res <- residuals_std_GMkd(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc))
                return(std_res)
              } else{
                if(parameters==0){
                  parameters <- sv_default$GMki
                }
                std_res <- residuals_std_GMki(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc))
                return(std_res)
              }
            }}
          }
        }
      }
    }
  }
}

## - Probability of negative rates
#' @title prob_neg_mu
#'
#' @description Probability of negative rates when projected over n years, based on simulated values of the latent variables
#'
#' @param model Affine model. E.g. for the Blackburn-Sherris model we have model="BS", and so on
#' @param fact_dep Boolean parameter indicating whether estimate models with factor dependence (fact_dep=TRUE) or independence
#' @param n_factors Number of factors. For some models, these are set by default
#' @param parameters Starting value for the parameters. If not set, then default values will be used
#' @param data Table with the average mortality rates
#' @param years_proj Number of years of ahead-projected rates
#' @param n_simulations Number of simulations
#'
#' @return A data frame object that contains a summary of a sample that
#'     can later be converted to a TeX output using \code{overview_print}
#' @examples
#' data(toydata)
#' output_table <- overview_tab(dat = toydata, id = ccode, time = year)
#' @export
prob_neg_mu <- function(model="BS", fact_dep=FALSE, n_factors=3, parameters=0, data=data_default, years_proj=1, n_simulations=100000){
  if(model=="AFNS"){
    if(fact_dep==TRUE){
      if(parameters==0){ # - then use default starting value
        parameters <- sv_default$AFNSd
      }
      pr_neg <- pr_neg_rates_f_AFNSd(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
      return(pr_neg)
    } else{
      if(parameters==0){
        parameters <- sv_default$AFNSi
      }
        pr_neg <- pr_neg_rates_f_AFNSi(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
        return(pr_neg)
    }
  } else{
    if(model=="BS"){
      if(fact_dep==TRUE){
        if(parameters==0){
          if(n_factors==2){
            parameters <- list(x0=sv_default$BSd$x0[1:2], delta=sv_default$BSd$delta[1:3], kappa=sv_default$BSd$kappa[1:2], sigma_dg=sv_default$BSd$sigma_dg[1:2], Sigma_cov=sv_default$BSd$Sigma_cov[1], r=c(sv_default$BSd$r1, sv_default$BSd$r2, sv_default$BSd$rc)) # - change in order to take into account the number of factors
            pr_neg <- pr_neg_rates_f_BSd_2F(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
            return(pr_neg)
          } else{ # - the factors must be three

            parameters <- sv_default$BSd # - change in order to take into account the number of factors
            pr_neg <- pr_neg_rates_f_BSd_3F(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
            return(pr_neg)

          }
        } else{
          if(n_factors==2){
            pr_neg <- pr_neg_rates_f_BSd_2F(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
            return(pr_neg)
          } else{
            pr_neg <- pr_neg_rates_f_BSd_3F(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
            return(pr_neg)
          }
        }
      } else{ # - it is the Blackburn-Sherris model with independent factors
        if(parameters==0){
          parameters <- sv_default$BSi
        } else{
          pr_neg <- pr_neg_rates_f_BSi(n_sim=n_simulations, x0=parameters$x0[1:n_factors], delta=parameters$delta[1:n_factors], kappa=parameters$kappa[1:n_factors], sigma = parameters$sigma[1:n_factors], r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
          return(pr_neg)
        }
      }
    } else{
      if(model=="CIR"){
          pr_neg <- 0
          return(pr_neg)
      } else{
        if(model=="AFGNS"){
        if(fact_dep==TRUE){
          if(parameters==0){
            parameters <- sv_default$AFGNSd
          }
            pr_neg <- pr_neg_rates_f_AFGNSd(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
            return(pr_neg)
        } else{
          if(parameters==0){
            parameters <- sv_default$AFGNSi
          }
            pr_neg <- pr_neg_rates_f_AFGNSi(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
            return(pr_neg)
        }}else{
        if(model=="AFRNS"){
          if(fact_dep==TRUE){
            if(parameters==0){
              parameters <- sv_default$AFRNSd
            }
            pr_neg <- pr_neg_rates_f_AFRNSd(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
            return(pr_neg)
          } else{
            if(parameters==0){
              parameters <- sv_default$AFRNSi
            }
            pr_neg <- pr_neg_rates_f_AFRNSi(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
            return(pr_neg)
          }
        }else{
          if(model=="AFUNS"){
            if(fact_dep==TRUE){
              if(parameters==0){
                parameters <- sv_default$AFUNSd
              }
              pr_neg <- pr_neg_rates_f_AFUNSd(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
              return(pr_neg)
            } else{
              if(parameters==0){
                parameters <- sv_default$AFUNSi
              }
              pr_neg <- pr_neg_rates_f_AFUNSi(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
              return(pr_neg)
            }
          }else{
            if(model=="GMk"){
              if(fact_dep==TRUE){
                if(parameters==0){
                  parameters <- sv_default$GMkd
                }
                pr_neg <- pr_neg_rates_f_GMkd(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
                return(pr_neg)
              } else{
                if(parameters==0){
                  parameters <- sv_default$GMki
                }
                pr_neg <- pr_neg_rates_f_GMki(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
                return(pr_neg)
              }
            }
          }
        }
      }
    }
  }
}}


#======================== - Projection - ===================================

# Improve output presentation with a smarter use of print
#' @title affine_project
#'
#' @description Projected survival curves and average mortality rates over n years ahead projections
#'
#' @param model Affine model. E.g. for the Blackburn-Sherris model we have model="BS", and so on
#' @param fact_dep Boolean parameter indicating whether estimate models with factor dependence (fact_dep=TRUE) or independence
#' @param n_factors Number of factors. For some models, these are set by default
#' @param parameters Starting value for the parameters. If not set, then default values will be used
#' @param data Table with the average mortality rates
#' @param years_proj Number of years of ahead-projected rates
#'
#' @return A data frame object that contains a summary of a sample that
#'     can later be converted to a TeX output using \code{overview_print}
#' @examples
#' data(toydata)
#' output_table <- overview_tab(dat = toydata, id = ccode, time = year)
#' @export
affine_project <- function(model="BS", fact_dep=FALSE, n_factors=3, parameters=0, data=data_default, years_proj=1){
  if(model=="AFNS"){
    if(fact_dep==TRUE){
      if(parameters==0){ # - then use default starting value
        return(projection)
      }  # - otherwise use user input starting values
        projection <- S_t_AFNSd_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
        return(projection)
    } else{
      if(parameters==0){
        parameters <- sv_default$AFNSi
       }
        projection <- S_t_AFNSi_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
        return(projection)
    }
  } else{
    if(model=="BS"){
      if(fact_dep==TRUE){
        if(parameters==0){
          if(n_factors==2){
            parameters <- list(x0=sv_default$BSd$x0[1:2], delta=sv_default$BSd$delta[1:3], kappa=sv_default$BSd$kappa[1:2], sigma_dg=sv_default$BSd$sigma_dg[1:2], Sigma_cov=sv_default$BSd$Sigma_cov[1], r=c(sv_default$BSd$r1, sv_default$BSd$r2, sv_default$BSd$rc)) # - change in order to take into account the number of factors
            projection <- S_t_BSd_2F_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
            return(projection)
          } else{ # - the factors must be three
            parameters <- sv_default$BSd # - change in order to take into account the number of factors
            projection <- S_t_BSd_3F_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
            return(projection)
          }
        } else{
          if(n_factors==2){
            projection <- S_t_BSd_2F_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
            return(projection)
          } else{
            projection <- S_t_BSd_3F_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
            return(projection)
          }
        }
      } else{ # - it is the Blackburn-Sherris model with independent factors
        if(parameters==0){
          parameters <- sv_default$BSi
        }
          projection <- S_t_BSi_proj(x0=parameters$x0[1:n_factors], delta=parameters$delta[1:n_factors], kappa=parameters$kappa[1:n_factors], sigma = parameters$sigma[1:n_factors], r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
          return(projection)
      }
    } else{
      if(model=="CIR"){
        if(parameters==0){
          parameters <- sv_default$CIR
        }
          projection <- S_t_CIR_proj(x0=parameters$x0[1:n_factors], delta=parameters$delta[1:n_factors], kappa=parameters$kappa[1:n_factors], sigma=parameters$sigma[1:n_factors], theta_P=parameters$theta_P[1:n_factors], r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
          return(projection)
      } else{
        if(model=="AFGNS"){
        # - if none of the previous model was selected, then it is an AFGNS
        if(fact_dep==TRUE){
          if(parameters==0){
            parameters <- sv_default$AFGNSd
          }
            projection <- S_t_AFGNSd_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
            return(projection)
        } else{
          if(parameters==0){
            parameters <- sv_default$AFGNSi
          }
            projection <- S_t_AFGNSi_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
            return(projection)
        }
        }else{
        if(model=="AFRNS"){
          if(fact_dep==TRUE){
            if(parameters==0){
              parameters <- sv_default$AFRNSd
            }
            projection <- S_t_AFRNSd_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
            return(projection)
          } else{
            if(parameters==0){
              parameters <- sv_default$AFRNSi
            }
            projection <- S_t_AFRNSi_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
            return(projection)
          }
        }else{
          if(model=="AFUNS"){
            if(fact_dep==TRUE){
              if(parameters==0){
                parameters <- sv_default$AFUNSd
              }
              projection <- S_t_AFUNSd_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
              return(projection)
            } else{
              if(parameters==0){
                parameters <- sv_default$AFUNSi
              }
              projection <- S_t_AFUNSi_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
              return(projection)
            }
          }else{
            if(model=="GMk"){
              if(fact_dep==TRUE){
                if(parameters==0){
                  parameters <- sv_default$GMkd
                }
                projection <- S_t_GMkd_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
                return(projection)
              } else{
                if(parameters==0){
                  parameters <- sv_default$GMki
                }
                projection <- S_t_GMki_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
                return(projection)
              }
            }
          }
        }
      }
    }
  }
}}


#======================= - Graphics - ======================================

## - Heatmap of residuals
# Improve output presentation with a smarter use of print
#' @title heatmap_res
#'
#' @description Returns the heatmaps created by using heatmaply
#'
#' @param residuals Table of residuals obtainable using std.res
#' @param color TRUE if colored (default) or FALSE if black and white
#'
#' @return A heatmap of the residuals
#' @examples
#' data(toydata)
#' output_table <- overview_tab(dat = toydata, id = ccode, time = year)
#' @export
heatmap_res <- function(residuals, color=TRUE){
  if(color==FALSE){
    heatmaply::heatmaply(residuals[c(nrow(residuals):1),], dendrogram = "none", xlab = "Year", ylab="Age", dynamicTicks=FALSE, hide_colorbar=FALSE, margins= c(0,0,0,0), color = c("white", "black"))
  } else{
    heatmaply::heatmaply(residuals[c(nrow(residuals):1),], dendrogram = "none", xlab = "Year", ylab="Age", dynamicTicks=FALSE, hide_colorbar=FALSE, margins= c(0,0,0,0))
  }
}

#======================= - Parameter uncertainty - ======================================

## - Fix by adding a print function for partial results
#' @title par_cov
#'
#' @description Estimation of the uncertainty of the parameters by Bootstrap or Multiple imputation
#'
#' @param method MI if Multiple Imputation or Bootstrap if the Bootstrap is desired
#' @param model Affine model. E.g. for the Blackburn-Sherris model we have model="BS", and so on
#' @param fact_dep Boolean parameter indicating whether estimate models with factor dependence (fact_dep=TRUE) or independence
#' @param n_factors Number of factors. For some models, these are set by default
#' @param parameters Starting value for the parameters. If not set, then default values will be used
#' @param data Table with the average mortality rates
#' @param D_se Number of Iterations if Multiple Imputation is used (default set to 50)
#' @param BS_s Number of Boostrap sample if this method is used
#' @param t_excl Number of time-periods to be excluded for stability of the Bootstrap method
#' @param max_iter Maximum number of iterations for the Coordinate Ascent estimation algorithm in absence of convergence
#' @param tolerance Minimum value of improvement of the log-likelihood value before convergence
#' @param wd Working directory for saving the partial output of the estimation process
#'
#' @return A data frame object that contains a summary of a sample that
#'     can later be converted to a TeX output using \code{overview_print}
#' @examples
#' data(toydata)
#' output_table <- overview_tab(dat = toydata, id = ccode, time = year)
#' @export
par_cov <- function(method="MI", model="BS", fact_dep=FALSE, n_factors=3, parameters, data=data_default, D_se=50, BS_s=500, t_excl=4, max_iter=200, tolerance=0.1, wd=0){

  if(method=="MI"){
    if(model=="AFNS"){
      if(fact_dep==TRUE){
          covariance <- CovEst_MI_AFNSd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
          return(covariance)

      } else{
          covariance <- CovEst_MI_AFNSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
          return(covariance)
      }
    } else{
      if(model=="BS"){
        if(fact_dep==TRUE){
            if(n_factors==2){
              covariance <- CovEst_MI_BSd_2F(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
              return(covariance)
            } else{ # - the factors must be three
              covariance <- CovEst_MI_BSd_3F(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
              return(covariance)
           }

        } else{ # - it is the Blackburn-Sherris model with independent factors

            covariance <- CovEst_MI_BSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
            return(covariance)

        }
      } else{
        if(model=="CIR"){
            covariance <- CovEst_MI_CIR(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, theta_P = parameters$theta_P, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
            return(covariance)

        } else{
          if(model=="AFGNS"){   # - if none of the previous model was selected, then it is an AFGNS
          if(fact_dep==TRUE){
              covariance <- CovEst_MI_AFGNSd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
              return(covariance)

          } else{
              covariance <- CovEst_MI_AFGNSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
              return(covariance)

          }
          }else{
          if(model=="AFRNS"){
            if(fact_dep==TRUE){
              covariance <- CovEst_MI_AFRNSd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
              return(covariance)

            } else{
              covariance <- CovEst_MI_AFRNSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
              return(covariance)

            }
          }else{
            if(model=="AFUNS"){
              if(fact_dep==TRUE){
                covariance <- CovEst_MI_AFUNSd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
                return(covariance)

              } else{
                covariance <- CovEst_MI_AFUNSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
                return(covariance)

              }
            }else{
              if(model=="GMk"){
                if(fact_dep==TRUE){
                  covariance <- CovEst_MI_GMkd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
                  return(covariance)

                } else{
                  covariance <- CovEst_MI_GMki(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
                  return(covariance)

                }
              }
            }
          }
        }
      }
    }}
    }else{ #### - Perform Bootstrap
    if(model=="AFNS"){
      if(fact_dep==TRUE){
          covariance <- CovEst_BS_AFNSd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
          return(covariance)

      } else{
          covariance <- CovEst_BS_AFNSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
          return(covariance)

      }
    } else{
      if(model=="BS"){
        if(fact_dep==TRUE){
            if(n_factors==2){
              covariance <- CovEst_BS_BSd_2F(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
              return(covariance)
            } else{ # - the factors must be three
              covariance <- CovEst_BS_BSd_3F(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
              return(covariance)

            }

        } else{ # - it is the Blackburn-Sherris model with independent factors
            covariance <- CovEst_BS_BSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
            return(covariance)

        }
      } else{
        if(model=="CIR"){
            covariance <- CovEst_BS_CIR(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, theta_P = parameters$theta_P, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
            return(covariance)

        } else{
          if(model=="AFGNS"){
          if(fact_dep==TRUE){
              covariance <- CovEst_BS_AFGNSd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
              return(covariance)

          } else{
              covariance <- CovEst_BS_AFGNSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
              return(covariance)
          }
          }else{
          if(model=="AFRNS"){
            if(fact_dep==TRUE){
              covariance <- CovEst_BS_AFRNSd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
              return(covariance)

            } else{
              covariance <- CovEst_BS_AFRNSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
              return(covariance)
            }
          }else{
            if(model=="AFUNS"){
              if(fact_dep==TRUE){
                covariance <- CovEst_BS_AFUNSd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
                return(covariance)

              } else{
                covariance <- CovEst_BS_AFUNSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
                return(covariance)
              }
            }else{
              if(model=="GMk"){
                if(fact_dep==TRUE){
                  covariance <- CovEst_BS_GMkd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
                  return(covariance)

                } else{
                  covariance <- CovEst_BS_GMki(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
                  return(covariance)
                }
              }
            }
          }
        }
      }
    }
  }}
}















