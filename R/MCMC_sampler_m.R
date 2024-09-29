#' Fit DLIM-IM
#' @description Create posterior samples from a DLIM-IM
#' @seealso Type \code{vignette('dlimOverview')} for a detailed description.
#' @seealso \link[dlimIM]{pred_m}
#' @seealso \link[dlimIM]{plot_bdlim}
#' @export
#' @importFrom splines ns
#' @param x matrix of exposure history (columns) for individuals (rows) (class "\code{matrix}")
#' @param y vector of response values (class "\code{numeric}")
#' @param trials number of trials for each observation when \code{family = "binomial"} (class "\code{numeric}")
#' @param M numeric matrix of modifying values (columns) for individuals with values between 0 and 1 (rows) (class "\code{matrix}")
#' @param z matrix of covariates, not including the modifiers (class "\code{matrix}")
#' @param df_m degrees of freedom for modifier basis (class "\code{numeric}")
#' @param df_l degrees of freedom for exposure time basis (class "\code{numeric}")
#' @param tau2 prior variance for cross-basis regression coefficients (class "\code{numeric}")
#' @param xi2 prior variance for covariate and modifier regression coefficients (class "\code{numeric}")
#' @param a prior shape for error variance for Gaussian model class "\code{numeric}")
#' @param b prior scale for error variance for Gaussian model class "\code{numeric}")
#' @param niter number of MCMC iterations class "\code{numeric}")
#' @param burnin number of MCMC iterations to remove as warm-up class "\code{numeric}")
#' @param model_type modifier spline basis construction ("linear" or "ns" for natural splines) class "\code{character}")
#' @param var_select TRUE to indicate performing selection on the modifiers class "\code{logical}")
#' @param weights_prior vector of Dirichlet shape parameters for each modifier for their prior distribution class "\code{numeric}")
#' @param WAIC TRUE to indicate computing WAIC class "\code{logical}")
#' @param prior_inclusions vector of prior inclusion probabilities for each weight class "\code{numeric}")
#' @param family family object specifying the distribution in fitting, either 'gaussian' or 'binomial' class "\code{character}")
#' @return This function returns a list that is an object of class "\code{dlimIM}" with the following components
#' \item{posterior}{posterior samples including burn-in (class "\code{matrix}")}


#M needs to be a matrix
MCMC_sampler_m <- function(x, y, trials = NULL, M, z, df_m, df_l, tau2, xi2,
                           a=NULL, b=NULL, niter, burnin, model_type = "ns",
                           var_select = FALSE, weights_prior = NULL, WAIC=FALSE,
                           prior_inclusions = NULL, family = "gaussian"){

  #stop if a and b not provided but family is Gaussian
  if(family == "gaussian" & (is.null(a) & is.null(b))){
    stop("Family is Gaussian but hyperparameters for error variance not provided. Please provide values for arguments a and b.")
  }

  #for binomial model, set y to be number of successes - trials/2 (kappa in Polson)
  if(family == "binomial"){
    y <- y - trials/2
  }

  #check to see if modifiers are between 0 and 1
  if(sum(M > 1) + sum(M < 0) != 0){
    stop("Modifier values must be between 0 and 1.")
  }

  x <- as.matrix(x)
  L <- ncol(x)
  n_chains <- 1

  #### MCMC initialization ###

   if(model_type=="linear"){
     df_m <- 2
   }
    z <- as.data.frame(z)
    z <- model.matrix(~ 0+.,model.frame(~ ., z, na.action=na.pass)) #handles factor covariates and missing values

    #number of modifiers in model
    n_m <- ncol(M) #make ncol(M) when all modifiers are covariates

    #set prior shapes for weights if NULL
    if(is.null(weights_prior)){
      shapes <- rep(1,n_m)
    }else{
      shapes <- weights_prior
    }

    #set prior inclusion probabilities for weights if NULL
    if(is.null(prior_inclusions)){
      prior_inclusions <- rep(0.5,n_m)
    }

    #set delta
    delta <- rep(0.2, n_m)

    #create sigma inverse matrix
    #intercept, cross-basis coefficients, all modifiers, covariates
    sigma_inv <- diag(c(0,rep(1/tau2,df_l*df_m),rep(1/xi2,n_m+ncol(z))))

    # create space for MCMC samples
    n_coefs <- 1 + df_l*df_m + n_m + ncol(z) #intercept & cross-basis, modifier, covariate coefs
    coef_samples <- matrix(NA, nrow=niter, ncol=n_coefs)
    colnames(coef_samples) <- c("Intercept",
                                paste0("CB", 1:(df_l*df_m)),
                                paste0("Modifier", 1:n_m),
                                paste0("Z", 1:ncol(z)))
    UN_weight_samples <- matrix(NA, nrow=niter, ncol=ncol(M)) #un-normalized weights (i.e. a_l)
    colnames(UN_weight_samples) <- paste0("UNweight", 1:ncol(M))
    weight_samples <- matrix(NA, nrow=niter, ncol=ncol(M)) #normalized weights (i.e. a_l)
    colnames(weight_samples) <- paste0("weight", 1:ncol(M))
    if(family == "gaussian"){
      sigma2_samples <- rep(NA, niter)
      names(sigma2_samples) <- "sigma2"
    }else if(family == "binomial"){
      omega_samples <- matrix(NA, nrow=niter, ncol=nrow(x))
      colnames(omega_samples) <- paste0("omega", 1:nrow(x))
    }else{
      stop("Family specified not supported.")
    }


    #create vector to help calculate WAIC
    if(WAIC == TRUE){
      # loglikmat <- matrix(NA, nrow=niter-burnin, ncol=nrow(x))
      loglikmat <- data.frame(Sum_lik = numeric(nrow(x)),
                              Sum_loglik = numeric(nrow(x)),
                              Sum_loglik2 = numeric(nrow(x)))
    }

    # initialize parameters arbitrarily from prior
    coef_samples[1,] <- sigma_inv%*%rnorm(ncol(sigma_inv)) #sampling from MVN(0,(sigma_inv)^}{-1})
    UN_weight_samples[1,] <- rMVgamma(shapes = shapes) #draw from gamma(1,1)
    weight_samples[1,] <- UN_weight_samples[1,]/sum(UN_weight_samples[1,])
    if(family == "gaussian"){
      sigma2_samples[1] <- 1/rgamma(1, shape = a, rate = b) #draw from IG(a,b)
    }else if(family == "binomial"){
      omega_samples[1,] <- rep(1, nrow(x))
    }else{
      stop("Family specified not supported.")
    }

    #create constant lag basis
    B_lag <- ns(0:(L-1), df=df_l, intercept = T) #L+1xdf_l

    #initial m_star to make basis
    if(model_type=="ns"){
      B_mod <- ns(rowMeans(M), df=df_m, intercept = T, Boundary.knots = c(0,1))
    }else if(model_type == "linear"){
      B_mod <- cbind(rep(1,length(y)), rowMeans(M))
    }

    #create design for the first time, update rest of times
    U <- weight_create_update(UN_weight_sample = UN_weight_samples[1,],
                              M = M,
                              x = x,
                              B_lag = B_lag,
                              B_mod = B_mod,
                              z = z,
                              model_type = model_type)

    #acceptance rate tracker
    tracker <- vector(mode = "list", length = ncol(M))
    delta_vec <- c()
    acc_rate_vec <- c()

    #profvis({ #profiling
    ###MCMC iterating###
    for(s in 2:niter) {
      #set.seed(100223 + s)
      print(s)

      # regular Gibbs update for reg coefs
      if(family == "gaussian"){
        coef_samples[s,] = update_coefs(y = y,
                                        U = U,
                                        sigma_inv = sigma_inv,
                                        sigma2 = sigma2_samples[s-1])
      }else if(family == "binomial"){
        coef_samples[s,] = update_coefs_logit(y = y,
                                              U = U,
                                              sigma_inv = sigma_inv,
                                              omegas = omega_samples[s-1,])
      }else{
        stop("Family not supported.")
      }


      if(family == "gaussian"){
        #Gibbs update for sigma2
        sigma2_samples[s] = update_sigma2(a = a,
                                          b = b,
                                          y = y,
                                          U = U,
                                          Psi = coef_samples[s,])
      }else if(family == "binomial"){
        #Gibbs update for omegas
        omega_samples[s,] = rpg(nrow(x), trials, U%*%coef_samples[s,])  #pgdraw(trials, U%*%coef_samples[s,]) #
      }else{
        stop("Family not supported.")
      }

      if(s > 500 & s <= burnin){
        acc_rate <- unlist(lapply(tracker, calc_acc_rate))
        tracker <- update_tracker(tracker)
      }else{
        acc_rate <- rep(NA, ncol(UN_weight_samples))
      }

      #MH step for un-normalized weights (doing last so that U does not need to be recalculated)
      updated <-  update_UN_weights_MH(UN_weights = UN_weight_samples[s-1,],
                                       delta = delta,
                                       y = y,
                                       sigma2 = sigma2_samples[s],
                                       omegas = omega_samples[s,],
                                       U_current = U,
                                       Psi = coef_samples[s,],
                                       x = x,
                                       M = M,
                                       z = z,
                                       B_lag = B_lag,
                                       B_mod = B_mod,
                                       acc_rate = acc_rate,
                                       model_type = model_type,
                                       var_select = var_select,
                                       tracker = tracker,
                                       shapes = shapes,
                                       prior_inclusions = prior_inclusions,
                                       family = family)
      UN_weight_samples[s,] <- updated$UN_weights_curr
      delta <- updated$delta
      if(sum(UN_weight_samples[s,]) != 0){
        weight_samples[s,] <- UN_weight_samples[s,]/sum(UN_weight_samples[s,])
      }else{
        weight_samples[s,] <- UN_weight_samples[s,]
      }
      tracker <- updated$tracker

      #remove
      delta_vec <- rbind(delta_vec, delta)
      acc_rate_vec <- rbind(acc_rate_vec, acc_rate)

      #weight modifiers, update cross-basis, and design matrix with intercept
      #I made this a function since I recreate U in the MH step, and I want changes to be consistent
      U <- weight_create_update(UN_weight_sample = UN_weight_samples[s-1,],
                                M = M,
                                x = x,
                                B_lag = B_lag,
                                B_mod = B_mod,
                                z = z,
                                U = U,
                                model_type = model_type)

      #log likelihood matrix (iterations by observations)
      if(WAIC & s > burnin){
        # loglikmat[(s-burnin),] <- loglikelihood(y = y,
        #                                         U = U,
        #                                         Psi = coef_samples[s,],
        #                                         sigma2 = sigma2_samples[s])
        loglikvec <- loglikelihood(y = y,
                                   U = U,
                                   Psi = coef_samples[s,],
                                   sigma2 = sigma2_samples[s])
        loglikmat$Sum_lik <- exp(loglikvec) + loglikmat$Sum_lik
        loglikmat$Sum_loglik <- loglikvec + loglikmat$Sum_loglik
        loglikmat$Sum_loglik2 <- loglikvec^2 + loglikmat$Sum_loglik2
      }

    }
    #})
    #print(delta)
    if(family == "gaussian"){
      posterior <- list(cbind(coef_samples, weight_samples, UN_weight_samples, sigma2_samples))
    }else if(family == "binomial"){
      posterior <- list(cbind(coef_samples, weight_samples, UN_weight_samples))
    }else{

    }
  #}#end foreach

  if(WAIC){
    WAIC_calc <- compute_WAIC(loglikmat = loglikmat,
                              niter = niter,
                              burnin = burnin)
    # WAIC_calc <- waic(loglikmat)$estimates[3,1]
    # rel_n_eff <- relative_eff(exp(loglikmat), chain_id = rep(1, ncol(loglikmat)))
    # LOO_calc <- loo(loglikmat, r_eff = rel_n_eff)$estimates[3,1]
  }

  names(posterior) <- paste0("chain",1:n_chains)
  attr(posterior, "model_type") <- model_type
  attr(posterior, "L") <- L
  attr(posterior, "df_l") <- df_l
  attr(posterior, "df_m") <- df_m
  attr(posterior, "niter") <- niter
  attr(posterior, "burnin") <- burnin
  attr(posterior, "delta") <- delta
  attr(posterior, "M") <- M
  attr(posterior, "var_select") <- var_select
  attr(posterior, "weights_prior") <- weights_prior
  attr(posterior, "family") <- family
  if(WAIC){
    attr(posterior, "WAIC") <- WAIC_calc
    # attr(posterior, "LOO") <- LOO_calc
  }

  message("Done!")

  class(posterior) <- "dlimIM"

  #stopCluster(cl)
  return(posterior)

}

