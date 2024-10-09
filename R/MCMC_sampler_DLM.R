#' Fit DLM with MCMC
#' @description Fit DLM with MCMC
#' @export
#' @importFrom dlnm crossbasis
#' @param x matrix of exposure history (columns) for individuals (rows) (class "\code{matrix}")
#' @param y vector of response values (class "\code{numeric}")
#' @param z matrix of covariates, not including the modifiers (class "\code{matrix}")
#' @param df_l degrees of freedom for exposure time basis (class "\code{numeric}")
#' @param tau2 prior variance for cross-basis regression coefficients (class "\code{numeric}")
#' @param xi2 prior variance for covariate and modifier regression coefficients (class "\code{numeric}")
#' @param a prior shape for error variance for Gaussian model class "\code{numeric}")
#' @param b prior scale for error variance for Gaussian model class "\code{numeric}")
#' @param niter number of MCMC iterations class "\code{numeric}")
#' @param burnin number of MCMC iterations to remove as warm-up class "\code{numeric}")
#' @return This function returns a list that is an object of class "\code{dlimIM}" with the following components
#' \item{posterior}{posterior samples including burn-in (class "\code{matrix}")}


Gibbs_sampler <- function(x, y, z, df_l, niter, burnin, tau2, xi2, a, b){


  ###DLIM set up and cross-basis###

  #cross-basis
  cb <- crossbasis(x = x,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", df = df_l))

  #create design matrix
  U <- as.matrix(data.frame(intercept = rep(1,nrow(x)),
                            CB = cb,
                            z = z))

  ###MCMC initialization ###

  #create sigma inverse matrix
  sigma_inv <- diag(c(0,rep(1/tau2,ncol(cb)),rep(1/xi2,ncol(z))))

  # create space for MCMC samples
  coef_samples = matrix(NA, nrow=niter,ncol=ncol(U))
  colnames(coef_samples) <- colnames(U)
  sigma2_samples = rep(NA, niter)
  names(sigma2_samples) <- "sigma2"

  #create vector to help calculate WAIC
  # loglikmat <- matrix(NA, nrow=niter-burnin, ncol=nrow(x))
  loglikmat <- data.frame(Sum_lik = numeric(nrow(x)),
                          Sum_loglik = numeric(nrow(x)),
                          Sum_loglik2 = numeric(nrow(x)))

  # initialize parameters arbitrarily from prior
  coef_samples[1,] = sigma_inv %*% matrix(rnorm(ncol(sigma_inv)), ncol=1) #sampling from MVN(0,(sigma_inv)^}{-1})
  sigma2_samples[1] = 1/rgamma(1, shape = a, rate = b)


  ###MCMC iterating###

  for(s in 2:niter) {
    # regular Gibbs update
    coef_samples[s,] = update_coefs(y = y,
                                    U = U,
                                    sigma_inv = sigma_inv,
                                    sigma2 = sigma2_samples[s-1])
    sigma2_samples[s] = update_sigma2(a = a,
                                      b = b,
                                      y = y,
                                      U = U,
                                      Psi = coef_samples[s,])

    #log likelihood matrix (iterations by observations)
    if(s > burnin){
      # loglikmat[(s-burnin),] <- loglikelihood(y = y,
      #                                         U = U,
      #                                         Psi = coef_samples[s,],
      #                                         sigma2 = sigma2_samples[s])
      #
      loglikvec <- loglikelihood(y = y,
                                 U = U,
                                 Psi = coef_samples[s,],
                                 sigma2 = sigma2_samples[s])
      loglikmat$Sum_lik <- exp(loglikvec) + loglikmat$Sum_lik
      loglikmat$Sum_loglik <- loglikvec + loglikmat$Sum_loglik
      loglikmat$Sum_loglik2 <- loglikvec^2 + loglikmat$Sum_loglik2
    }
  }

  #compute WAIC
  # WAIC_calc <- compute_WAIC(loglikmat = loglikmat,
  #                           niter = niter)
  WAIC_calc <- compute_WAIC(loglikmat = loglikmat,
                            niter = niter,
                            burnin = burnin)
  # WAIC_calc <- waic(loglikmat)$estimates[3,1]
  # rel_n_eff <- relative_eff(exp(loglikmat), chain_id = rep(1, ncol(loglikmat)))
  # LOO_calc <- loo(loglikmat, r_eff = rel_n_eff)$estimates[3,1]

  #save posterior
  posterior <- cbind(coef_samples, sigma2_samples)

  #add WAIC to posterior
  attr(posterior, "WAIC") <- WAIC_calc
  # attr(posterior, "LOO") <- LOO_calc

  print("Done!")

  return(posterior)

}

