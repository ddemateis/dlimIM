### updating functions ### 
library(dlim)
library(Matrix)
library(dlnm)

#compute log likelihood
loglikelihood <- function(y, U, Psi, sigma2){

  Psi <- matrix(Psi, ncol=1)
  ll_mat <- dnorm(y - U%*%Psi, sd = sigma2, log = T)
  #ll_mat <-  -1/2*log(2*pi) - log(sqrt(sigma2)) - (y - U%*%Psi)^2/(2*sigma2)
  return(ll_mat)
}

#compute WAIC
compute_WAIC <- function(loglikmat, niter,burnin){
  # #computed log point-wise predictive density
  # # \text{computed lppd} = \sum_{i=1}^n \text{log}\left[ \frac{1}{S} \sum_{s=1}^S p\left(y_i | \boldsymbol{\phi}^{(s)} \right) \right] 
  # lppd <- sum(log(colMeans(exp(loglikmat)))) #should be sum not mean
  # 
  # #computed p_WAIC_2
  # #\text{computed } p_{\text{WAIC}_2} = \sum_{i=1}^n \frac{1}{S-1} \sum_{s=1}^S \left[ \text{log } p \left(y_i |  \boldsymbol{\phi}^{(s)} \right) - \frac{1}{S} \sum_{s=1}^S \text{log } p \left(y_i |  \boldsymbol{\phi}^{(s)} \right) \right] 
  # lp_mean <- colMeans(loglikmat) #mean across iterations of loglikelihood
  # lp_mean_mat <-matrix(rep(lp_mean,nrow(loglikmat)), nrow=nrow(loglikmat), byrow = T)
  # p_waic_2 <- sum(colSums((loglikmat - lp_mean_mat)^2) / (niter - 1))
  # 
  # #computed WAIC
  # WAIC_calc <- -2*lppd + 2*p_waic_2
  
  S <- niter - burnin #number of posterior samples post-burn-in
  lppd <- sum(log(loglikmat$Sum_lik/S)) #should be sum not mean
  
  #computed p_WAIC_2
  #\text{computed } p_{\text{WAIC}_2} = \sum_{i=1}^n \frac{1}{S-1} \sum_{s=1}^S \left[ \text{log } p \left(y_i |  \boldsymbol{\phi}^{(s)} \right) - \frac{1}{S} \sum_{s=1}^S \text{log } p \left(y_i |  \boldsymbol{\phi}^{(s)} \right) \right] 
  p_waic_2 <- sum((loglikmat$Sum_loglik2 - (loglikmat$Sum_loglik)^2/S)/(S - 1))
  
  #computed WAIC
  WAIC_calc <- -2*lppd + 2*p_waic_2
  
  return(WAIC_calc)
}

#Gibbs update of regression coefficients
update_coefs <- function(y, U, sigma_inv, sigma2){
  
  y <- matrix(y,ncol=1)
  
  V <- chol2inv(chol(sigma_inv + sigma2^(-1)*t(U)%*%U))
  m <- sigma2^(-1)*V%*%t(U)%*%y
  p <- ncol(V)
  
  return(m + t(chol(V))%*%rnorm(p)) #chol(V) is upper triangular
}

#Gibbs updtae of reg coeffs for logit
update_coefs_logit <- function(y, U, sigma_inv, omegas){
  
  y <- matrix(y,ncol=1)
  
  V <- chol2inv(chol(sigma_inv + t(U)%*%diag(omegas)%*%U))
  m <- V%*%t(U)%*%y
  p <- ncol(V)
  
  return(m + t(chol(V))%*%rnorm(p)) #chol(V) is upper triangular
}

#Gibbs update of model errors
update_sigma2 <- function(a, b, y, U, Psi){
  n <- length(y)
  y <- matrix(y,ncol=1)
  Psi <- matrix(Psi,ncol=1)
  
  sigma2_shape <- a + n/2
  sigma2_rate <- b + 0.5*sum((y - U%*%Psi)^2)
  
  return(rinvgamma(1, shape=sigma2_shape, rate=sigma2_rate))
}

#Gibbs update of tau2
update_tau2 <- function(a, b, theta, K_prec){
  
  theta <- matrix(theta, ncol=1)
  theta_p <- matrix(theta, nrow=1)
  
  shape <- a + 1/2
  rate <- b + 0.5*theta_p %*% K_prec %*% theta
  
  return(rinvgamma(1, shape=shape, rate=rate))
}

#MH update of un-normalized (UN) weights
update_UN_weights_MH <- function(UN_weights, delta, y, sigma2=NULL,
                                 omegas=NULL, U_current, Psi, x, M,
                                 z, B_lag, B_mod, acc_rate, model_type,
                                 var_select, tracker, shapes, prior_inclusions,
                                 family) {
  
  n <- length(y)
  UN_weights_curr <- UN_weights
  
  #make Psi a column vector
  Psi <- matrix(Psi, ncol=1)
  
  #current log likelihood
  if(family == "gaussian"){
    loglik <- -n * log(sqrt(sigma2)) - 
      sum((y - U_current%*%Psi)^2) / (2 * sigma2)
  }else if(family == "binomial"){
    Zs <- y / omegas
    loglik <- -1/2 * t(Zs - U_current%*%Psi)%*%diag(omegas)%*%(Zs - U_current%*%Psi)
  }else{
    stop("Family not supported.")
  }
  
  
  #update weights individually, not simultaneously 
  for(j in 1:length(UN_weights)){
    
    #adaptive step-size
    #increase if acceptance is high, decrease if low
    if(!is.na(acc_rate[[j]])){
      delta[j] <- delta[j]*(1 - (0.35 - acc_rate[j]))
    }
    
    #first, initialize proposed weights with weights that will be "new" if proposal is not accepted
    UN_weights_pro <- UN_weights_curr
    
    #propose a new value for the jth weight
    nu <- ifelse(UN_weights_curr[j] == 0, 1, 0) #0 or 1, do not reset after step 1
    
    #step 1
    if(var_select){ #if performing variable selection
      if(nu==1 & UN_weights_curr[j] == 0){ 
        #proposing non-zero, current is zero
        UN_weights_pro[j] <- rgamma(1,shape = shapes[j],1) #prior
      }else if(nu==0 & UN_weights_curr[j] != 0){ 
        #proposing zero, current is non-zero
        UN_weights_pro[j] <- 0
      }
      
      #accept/reject and update
      update_info <- acc_reject(UN_weights_pro = UN_weights_pro,
                                B_lag = B_lag,
                                U_current = U_current,
                                sigma2 = sigma2,
                                omegas = omegas,
                                Psi = Psi, 
                                UN_weights_curr = UN_weights_curr,
                                j = j,
                                var_select = var_select,
                                loglik = loglik,
                                x = x,
                                M = M,
                                z = z,
                                B_mod = B_mod,
                                model_type = model_type,
                                y = y,
                                shapes = shapes, 
                                prior_inclusion = prior_inclusions[j],
                                family = family)
      UN_weights_curr <- update_info$UN_weights_curr
      U_current <- update_info$U_current
      loglik <- update_info$loglik
    }
    
    
    #step 2
    if((nu==0 & UN_weights_curr[j]!=0) | !var_select){ #proposed a zero in step 1 and it was rejected OR not doing variable selection
      UN_weights_pro[j] <- abs(rnorm(1, UN_weights_curr[j], delta[j]))
      
      #accept/reject and update
      update_info <- acc_reject(UN_weights_pro = UN_weights_pro,
                                B_lag = B_lag,
                                U_current = U_current,
                                sigma2 = sigma2,
                                omegas = omegas,
                                Psi = Psi, 
                                UN_weights_curr = UN_weights_curr,
                                j = j,
                                var_select = var_select,
                                loglik = loglik,
                                x = x,
                                M = M,
                                z = z,
                                B_mod = B_mod,
                                model_type = model_type,
                                y = y,
                                shapes = shapes, 
                                prior_inclusion = prior_inclusions[j],
                                family = family)
      UN_weights_curr <- update_info$UN_weights_curr
      U_current <- update_info$U_current
      loglik <- update_info$loglik
      
      tracker[[j]] <- c(tracker[[j]],UN_weights_curr[j])
    }
    
    # #keep track of the movement from non-zero to non-zero weights for improving acceptance
    # if(UN_weights[j] != 0 & UN_weights_curr[j] != 0){
    #   tracker[[j]] <- c(tracker[[j]],UN_weights_curr[j])
    # }
    
  }
  
  result <- list(UN_weights_curr = UN_weights_curr,
                 delta = delta,
                 tracker = tracker)
  return(result)
}

acc_reject <- function(UN_weights_pro, B_lag, U_current, sigma2, omegas, Psi, UN_weights_curr,
                       j, var_select, loglik, x, M, z, B_mod, model_type, y, shapes,
                       prior_inclusion,
                       family = family){
  #make proposed design
  U_pro <- weight_create_update(UN_weight_sample = UN_weights_pro, 
                                x = x,
                                M = M,
                                z = z,
                                B_lag = B_lag,
                                B_mod = B_mod,
                                U = U_current,
                                model_type = model_type)
  
  # log likelihood with proposed weight
  if(family == "gaussian"){
    loglik_pro <- -n * log(sqrt(sigma2)) - sum((y - U_pro%*%Psi)^2) / (2 * sigma2)
  }else if(family == "binomial"){
    Zs <- y / omegas
    loglik_pro <- -1/2 * t(Zs - U_pro%*%Psi)%*%diag(omegas)%*%(Zs - U_pro%*%Psi)
  }else{
    stop("Family not supported.")
  }
  
  #ratio to determine acceptance (see model write up for derivations)
  if((UN_weights_curr[j] != 0 & UN_weights_pro[j] != 0) | !var_select){ #was non-zero, still non-zero or no variable selection
    r <- loglik_pro - loglik + 
      ((shapes[j]-1)*log(UN_weights_pro[j]) - UN_weights_pro[j]) - 
      ((shapes[j]-1)*log(UN_weights_curr[j]) - UN_weights_curr[j])
  }else if(UN_weights_pro[j] == 0){#proposing switch from non-zero to zero, and also doing variable selection
    r <- log(1 - prior_inclusion) - log(prior_inclusion) + loglik_pro - loglik
  }else{#proposing switch from zero to non-zero, and also doing variable selection
    r <- log(prior_inclusion) - log(1 - prior_inclusion) + loglik_pro - loglik
  }
  
  # accept with probability min(r, 1)
  u = runif(1, 0, 1)
  if(log(u) < r) {
    UN_weights_curr <- UN_weights_pro #update the weights
    U_current <- U_pro #update design
    loglik <- loglik_pro #update log likelihood
  }
  
  result <- list(UN_weights_curr = UN_weights_curr,
                 U_current = U_current,
                 loglik = loglik)
  
  return(result)
}

#not actually an updating function. Weights modifiers (creates m_star), creates cross-basis, and design matrix
#if you pass U, then it will just update the cross-basis
#if you do not pass U, it will create U
weight_create_update <- function(UN_weight_sample, M, x, B_lag, B_mod, z, U=NULL, model_type){
  
  #weight the modifiers into m*
  if(sum(UN_weight_sample)!=0){
    rho <- matrix(UN_weight_sample,ncol=1)/sum(UN_weight_sample) #check sums to 1
  }else{
    rho <- UN_weight_sample #all zeros
  }
  m_star <- M%*%rho
  df_l <- ncol(B_lag)
  
  #cross-basis
  XC <- as.matrix(x)%*%B_lag
  if(model_type == "linear"){
    XC_m <- sweep(XC, MARGIN=1, m_star, "*")
    cb <- cbind(XC,XC_m)
  }else if(model_type == "ns"){
    Bmod <- predict(B_mod, m_star)
    #B_mod <- ns(m_star, df=df_m, intercept = T, Boundary.knots = c(0,1))
    m_expd <- Bmod %x% matrix(rep(1), ncol=df_l) #n x df_m*df_l
    XC_rep <- matrix(rep(XC, df_m), nrow=n) #n x df_m*df_l
    cb <- m_expd * XC_rep #n x df_m*df_l
  }
  
  
  #create or update design
  if(is.null(U)){
    #create design matrix
    colnames(cb) <- paste0("CB", 1:ncol(cb))
    colnames(M) <- paste0("M", 1:ncol(M))
    colnames(z) <- paste0("Z", 1:ncol(z))
    intercept <- matrix(rep(1,length(m_star)),ncol=1)
    colnames(intercept) <- "Intercept"
    U <- cbind(intercept, cb, M, z)
  }else{
    cb_idx <- grep("CB",colnames(U)) #identify columns of U that are for the cross-basis
    U[,cb_idx] <- cb #replace cross-basis columns with new cross-basis
  }
  
  return(U)
}

calc_acc_rate <- function(x){
  q <- length(x)
  p <- ifelse(q < 500, 1, q-499)
  acc_rate <- ifelse(q == 500, mean(diff(x[p:q])!=0), NA)
  return(acc_rate)
}

update_tracker <- function(tracker){
  new_tracker <- vector(mode = "list", length = length(tracker))
  for(i in 1:length(tracker)){
    if(length(tracker[[i]])<500 & !is.null(tracker[[i]])){
      new_tracker[[i]] <- tracker[[i]]
    }
  }
  return(new_tracker)
}

#sample from gamma prior
rMVgamma <- function(shapes){
  samps <- c()
  for(shp in shapes){
    samps <- c(samps,rgamma(1, shape=shp))
  }
  return(samps)
}

#MH update of alpha (no longer can use, for penalized splines)
update_alpha <- function(alpha, K_current, theta, tau2, delta_a, acc_rate_a, Smod, Slag){
  
  alpha_curr <- alpha
  
  #make theta a column vector
  theta <- matrix(theta, ncol=1)
  theta_p <- matrix(theta, nrow=1)
  #current log likelihood
  loglik <- -0.5*determinant(K_current)$modulus - (theta_p%*%K_current%*%theta) / 2
  
  #adaptive step-size
  #decrease if acceptance is high, increase if low (this delta has inverse relationship)
  # if(!is.null(acc_rate_a)){
  #   delta_a <- delta_a*(1 + (0.35 - acc_rate_a))
  #   plot(seq(0,1,0.01), dbeta(seq(0,1,0.01), delta_a*alpha_curr, delta_a*(1-alpha_curr)), type="l")
  # }
  
  #propose alpha
  # p_curr <- delta_a*alpha_curr
  # l_curr <- delta_a*(1-alpha_curr)
  alpha_pro <- runif(1) #propose from prior #rbeta(1, p_curr, l_curr) #beta proposal
  
  #calculate proposed K
  K_pro <- (alpha_pro*Smod + (1-alpha_pro)*Slag)/tau2
  
  # proposed log likelihood
  loglik_pro <- -0.5*determinant(K_pro)$modulus - (theta_p%*%K_pro%*%theta) / 2
  
  #log proposal ratio (proposal distribution not symmetric)
  # p_pro <- delta_a*alpha_pro
  # l_pro <- delta_a*(1-alpha_pro)
  log_prop <- 0 #proposing from prior #dbeta(alpha_curr, p_pro, l_pro, log = T) - dbeta(alpha_pro, p_curr, l_curr, log = T) #proposing from beta
  
  #ratio to determine acceptance (see model write up for derivations)
  r <- loglik_pro - loglik + log_prop
  
  # accept with probability min(r, 1)
  u = runif(1, 0, 1)
  if(log(u) < r) {
    alpha_curr <- alpha_pro #update alpha
    K_current <- K_pro #update K
  } 
  
  result <- list(alpha = alpha_curr,
                 K = K_current,
                 delta_a = delta_a)
  
  return(result)
  
}
