#' Prediction for DLIM-IM
#' @description Predicted values based on a \code{dlimIM} object.
#' @seealso \link[dlimIM]{plot_bdlim}
#' @seealso \link[dlimIM]{MCMC_sampler_m}
#' @export
#' @importFrom splines ns
#' @param posterior_list posterior samples  (class "\code{dlimIM}")
#' @param chain if \code{posterior_list} has multiple chains, indicates which chain (class "\code{numeric}")
#' @param burnin number of MCMC samples to remove as warm-up (class "\code{numeric}")
#' @param thin number post-burn-in MCMC samples to thin by (class "\code{numeric}")
#' @param new_mods matrix of new modifier values for prediction ((class "\code{matrix}"))
#' @param m_star vector of new modifier index values for prediction (class "\code{numeric}")
#' @param alpha significance level (class "\code{numeric}")
#' @return This function returns a list of predicted values and credible bounds \code{list}
#' \item{betas_cumul}{cumulative effect estimates for each modifier index value (class "\code{numeric}")}
#' \item{betas_LB}{lower bound for cumulative effect estimates for each modifier index value (class "\code{numeric}")}
#' \item{betas_UB}{upper bound for cumulative effect estimates for each modifier index value (class "\code{numeric}")}
#' \item{betas}{point-wise effect estimates for each modifier index value and exposure-time point (class "\code{numeric}")}
#' \item{betas}{lower bound for point-wise effect estimates for each modifier index value and exposure-time point (class "\code{numeric}")}
#' \item{betas}{upper bound for point-wise effect estimates for each modifier index value and exposure-time point (class "\code{numeric}")}


pred_m <- function(posterior_list,
                   chain=1,
                   burnin,
                   thin = 1,
                   new_mods = NULL,
                   m_star = NULL,
                   alpha=0.05){

  #set up
  niter <- nrow(posterior_list[[chain]])
  inf_idx <- seq(burnin, niter, thin)
  posterior <- posterior_list[[chain]][inf_idx,]
  L <- attr(posterior_list, "L")
  df_l <- attr(posterior_list, "df_l")
  df_m <- attr(posterior_list, "df_m")
  model_type <- attr(posterior_list, "model_type")
  M <- attr(posterior_list, "M")
  est_dlim <- list()

  #obtain indices of coefficients that are for the cross-basis
  idx <- grep("CB",colnames(posterior))#includes only cross-basis elements

  #indices for weights
  idx_w <- grep("^weight",colnames(posterior))#includes only cross-basis elements

  #reconstruct B_lag
  B_lag <- ns(0:(L-1), df = df_l, intercept = T)
  B <- ns(rowMeans(M), df=df_m, intercept = T, Boundary.knots = c(0,1))


  ### Estimate Cumulative effects ###

  #function for computing matrix multiplication across posterior samples
  construct_beta_post <- function(post, m_star=NULL){

    #weight modifier if grid of weighted modifiers not passed
    if(is.null(m_star)){
      weights <- post[idx_w]
      names(weights) <- paste0("weight",1:length(idx_w))
      m_star <- new_mods %*% matrix(weights,ncol=1)
    }

    #number of weighted modifier values
    m <- length(m_star)

    #reconstruct B_mod for given modifiers
    if(model_type == "linear"){
      B_mod <- cbind(rep(1,m),m_star)
    }else if(model_type == "ns"){
      orig_m_star <- M%*%matrix(post[idx_w],ncol=1)
      B_mod <- predict(B, m_star)
    }

    #construct cross-basis for inference on new modifiers
    oneC <- matrix(1, ncol=L, nrow=m)%*%B_lag
    m_expd <- B_mod %x% matrix(rep(1), ncol=df_l) #n x df_m*df_l
    oneC_rep <- matrix(rep(oneC, df_m), nrow=m) #n x df_m*df_l
    cb_est <- m_expd * oneC_rep #n x df_m*df_l

    #compute effect estimate
    cb_est %*% matrix(post[idx],ncol=1)
  }

  #obtain posterior samples for the cumulative beta for each modifier
  betas_cumul_post <- matrix(apply(posterior,1,construct_beta_post, m_star = m_star), nrow=length(m_star))#rows are cumul betas, columns are samples

  #average over samples to obtain estimates for each cumulative beta
  betas_cumul <- matrix(rowMeans(betas_cumul_post),ncol=1)

  #take the quantiles to obtain CIs
  cumul_UB <- apply(betas_cumul_post, 1, quantile, probs=1-alpha/2)
  cumul_LB <- apply(betas_cumul_post, 1, quantile, probs=alpha/2)

  #construct list of cumulative beta info to return
  est_dlim$betas_cumul <- betas_cumul
  est_dlim$cumul_LB <- cumul_LB
  est_dlim$cumul_UB <- cumul_UB

  ### Estimate Pointwise effects ###

  #function for computing matrix multiplication across posterior samples
  construct_beta_pw_post <- function(post, m_star=NULL){

    #weight modifier if grid of weighted modifiers not passed
    if(is.null(m_star)){
      weights <- post[idx_w]
      names(weights) <- paste0("weight",1:length(idx_w))
      m_star <- new_mods %*% matrix(weights,ncol=1)
    }

    #number of modifiers
    m <- length(m_star)

    #reconstruct B_mod for given modifiers
    if(model_type == "linear"){
      B_mod <- cbind(rep(1,m),m_star)
    }else if(model_type == "ns"){
      orig_m_star <- M%*%matrix(post[idx_w],ncol=1)
      B_mod <- predict(B, m_star)
    }

    #construct cross-basis for inference on new modifiers
    coefs <- matrix(post[idx],ncol=1)
    t(matrix((B_mod %x% B_lag) %*% coefs, ncol=m)) #individuals x time
  }

  #obtain posterior samples for the cumulative beta for each modifier
  betas_post <- simplify2array(apply(posterior,1,construct_beta_pw_post, m_star = m_star, simplify = F)) #individual x time x posterior sample

  #average over samples to obtain estimates for each beta_t(m*_i)
  #apply will return class matrix/array but str will be list, need to unlist and then put back into matrix, otherwise cannot do calculations
  betas <- matrix(as.numeric(apply(betas_post, 1:2, mean, simplify = F)), nrow=dim(betas_post)[1]) #individuals x time

  #take the quantiles to obtain CIs
  UB <- matrix(as.numeric(apply(betas_post, 1:2, quantile, probs=1-alpha/2)), nrow=dim(betas_post)[1])
  LB <- matrix(as.numeric(apply(betas_post, 1:2, quantile, probs=alpha/2)), nrow=dim(betas_post)[1])

  #construct list of beta info to return
  est_dlim$betas <- betas
  est_dlim$LB <- LB
  est_dlim$UB <- UB

  return(est_dlim)
}
