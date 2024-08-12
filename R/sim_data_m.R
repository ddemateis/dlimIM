#' Simulate Data
#' @description Simulate data to use with the \pkg{dlimIM} package. There are different weight modification scenarios to choose for simulation.
#' @seealso \link[dlimIM]{sim_dlf_m}
#' @export
#' @param x a time series vector of length \code{n} or matrix of lagged exposures for \code{n} individuals (class "\code{numeric}", "\code{matrix}")
#' @param L a vector of length 1 containing the number of lag terms. This is required if \code{x} is vector, and is not used if \code{x} is a matrix (class "\code{numeric}")
#' @param M matrix containing modifying values for each individual (class "\code{matrix}")
#' @param w vector containing true weights for the modifier index (class "\code{numeric}")
#' @param noise a vector of length 1 containing the standard deviation for a normal distribution with mean 0 used to add noise to the simulated response values. Must proivde if \code{SNR} is not provided (class "\code{numeric}")
#' @param type a vector containing the number 1, 2, 3, or 4 for simulation modification type: none, linear, non-linear shift, non-linear shift with linear scale (class "\code{numeric}")
#' @param SNR The signal-to-noise ratio. If \code{SNR} is provided, but \code{noise} is not, \code{noise} is reset to be the standard deviation of the response, before adding noise.   (class "\code{numeric}")
#' @param ncovariates number of covariates to add to the model, numeric vector of length 1.
#' @param gamma True coefficient for the main effect of the modifier (class "\code{numeric}")
#' @param family Family object specifying the likelihood distribution for simulating. "gaussian" and "binomial" supported (class "\code{character}")
#' @param sf scale factor for logistic simulation (class "\code{numeric}")
#' @return This returns a list of 8 items:
#' \item{x}{a lagged exposure matrix. If \code{x} was a matrix, it is unchanged. (class "\code{matrix}")}
#' \item{L}{a numeric vector of length 1 containing the number of lag terms (class "\code{numeric}")}
#' \item{M}{the modifier matrix "\code{matrix}")}
#' \item{w}{the true modifier weights (class "\code{numeric}")}
#' \item{y}{a numeric vector of length \code{nrow(x)} containing the perturbed simulated response values. (class "\code{numeric}")}
#' \item{trials}{vector of number of trials for binomial model (class "\code{numeric}")}
#' \item{betas}{a matrix containing true coefficients for each lag/modifier combination, with each row representing a lag and each column a modifier (class "\code{matrix}")}
#' \item{betas_cumul}{a numeric vector of length \code{L+1} containing cumulative true coefficients for the lag terms, summed over modifiers (class "\code{numeric}")}
#' \item{Z}{covariates (class "\code{matrix}")}
#' \item{gammas}{true coefficients for the covariates (class "\code{numeric}")}
#' \item{noise2}{simulation variance, or square of noise (class "\code{numeric}")}
#' \item{family}{Family object specifying the likelihood distribution used for simulating. (class "\code{character}")}


sim_data_m <- function(x, L=NULL, M, w, noise=1, type=2, SNR,
                       ncovariates=0, gamma,
                       family="gaussian", sf=1){
  source("Functions/sim_dlf_m.R")
  #create lagged structure
  if(is.vector(x)){
    X <- Lag(x,0:L)[-c(1:L),]
    M <- M[-c(1:L),]
  }else{
    L <- ncol(x)-1
    X <- x
    M <- M
  }

  #create weighted modifier
  m_star <- M%*%w

  #Create Betas
  width=5
  if(family == "gaussian"){
    betas <- sim_dlf_m(L,m_star,type,width)
    betas_cumul <- colSums(betas)
    y_mean <- colSums(t(X)*betas)

    #if SNR is provided, but noise is not, reset noise to the SD based on data
    #if SNR is not provided, noise must be provided (default to noise of 1)
    if(!missing(SNR)){
      noise <- sd(y_mean)/SNR
    }
  }else if(family == "binomial"){
    betas <- sim_dlf_m(L,m_star,type,width)
    betas_cumul <- colSums(betas)
    y_mean <- colSums(t(X)*betas)
  }



  #Create gammas and covariates
  if(ncovariates!=0){
    Z <- matrix(rnorm(nrow(X)*(ncovariates)), ncol = ncovariates)
    mod_Z <- cbind(M,Z)
    if(family == "gaussian"){
      gammas <- c(gamma,matrix(rnorm(ncovariates),ncol=1))
      y <- y_mean + mod_Z%*%gammas + rnorm(length(y_mean),0,noise)
      trials <- NULL
    }else if(family == "binomial"){
      gammas <- c(gamma,matrix(rnorm(ncovariates),ncol=1))
      v <- y_mean + mod_Z%*%gammas
      v <- (v - mean(v))/sf #change scale to control logistic curve shape
      pr <- c(1/(1+exp(-v)))
      trials <- rep(1, nrow(X))#round(runif(nrow(x)),2)*100 + 1 #can't be 0
      y <- rbinom(nrow(X), as.integer(trials), pr)
      noise <- NULL
    }else{
      stop("Family not supported.")
    }
  }else{
    Z <- NULL
    gammas <- gamma
    if(family == "gaussian"){
      y <- y_mean + matrix(M,nrow=nrow(X))*gammas + rnorm(length(y_mean),0,noise)
      trials <- NULL
    }else if(family == "binomial"){
      gammas <- c(gamma,matrix(rnorm(ncovariates),ncol=1))
      v <- y_mean + matrix(M,nrow=nrow(X))*gammas
      v <- (v - mean(v))/sf #change scale to control logistic curve shape
      pr <- c(1/(1+exp(-v)))
      trials <- rep(1, nrow(x))#round(runif(nrow(x)),2)*100 + 1 #can't be 0
      y <- rbinom(nrow(x), as.integer(trials), pr)
      noise <- NULL
    }else{
      stop("Family not supported.")
    }
  }


  result <- list(x=X,
                 L=L,
                 M=M,
                 w = w,
                 m_star=m_star,
                 y=y,
                 trials=trials,
                 betas=betas,
                 betas_cumul=betas_cumul,
                 Z = Z,
                 gammas = gammas,
                 noise2 = noise^2,
                 family = family)

  return(result)
}
