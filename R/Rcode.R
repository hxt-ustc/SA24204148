#' @import latex2exp
#' @import ggplot2
#' @import bootstrap
#' @import GGally
#' @import boot
#' @import DAAG
#' @import coda
#' @import gridExtra
#' @import lpSolve
#' @import RcppArmadillo
NULL

vec_norm <- function(x) {sqrt(sum(x^2))}   

spec_norm_diff <- function(x, y, scale = T) {   
  if(length(x) != length(y)) stop('x and y are of different length')
  if(length(x) != 1 & scale) {
    x <- x/sqrt(sum(x^2))
    y <- y/sqrt(sum(y^2))
  }
  norm(x %*% t(x) - y %*% t(y), "2")
}

A_shift_unobserved_last <- function(A, weights) {
  
  n <- nrow(A)
  n_obs = sum(weights == 1)
  
  A_ori <- A
  unobs <- which(!weights) 
  obs <- which(weights == 1) 
  n_obs <- sum(weights == 1) 
  
  A[1:n_obs, 1:n_obs] <- A_ori[obs, obs]
  A[(1+n_obs):n, (1+n_obs):n] <- A_ori[unobs, unobs]
  A[1:n_obs, (1+n_obs):n] <- A_ori[obs, unobs]
  A[(1+n_obs):n, 1:n_obs] <- A_ori[unobs, obs]
  
  A
}

unobserved_shift_back <- function(xx, weights) {
  
  n <- length(xx)
  
  unobs <- which(!weights)
  obs <- which(weights == 1)
  n_obs <- sum(weights == 1)
  
  xx_tmp <- xx
  xx[obs] <- xx_tmp[1:n_obs]
  xx[unobs] <- xx_tmp[(1+n_obs):n]
  
  xx
}

adjust_sign <- function(ret) {
  
  if(is.null(ret$v)) {
    k <- length(ret$beta) - 1
    
    u_change_condition <- (sign(ret$u[which.max(abs(ret$u))]) < 0)
    
    if(u_change_condition) {
      ret$u <- -ret$u
      ret$beta[k+1] <- -ret$beta[k+1]
    }
  } else {
    k <- length(ret$beta) - 2
    
    u_condition <- (max(abs(ret$u)) > max(abs(ret$v)))
    u_change_condition <- (sign(ret$u[which.max(abs(ret$u))]) < 0)
    
    v_condition <- !u_condition
    v_change_condition <- (sign(ret$v[which.max(abs(ret$v))]) < 0)
    
    if( (u_condition & u_change_condition) | 
        (v_condition & v_change_condition)) {
      
      ret$u <- -ret$u; ret$v <- -ret$v;
      ret$beta[k+1] <- -ret$beta[k+1]
      ret$beta[k+2] <- -ret$beta[k+2]
      
    }
  }
  
  ret
}

fix_sign <- function(ret){
  if(is.null(ret$v)) {
    k <- length(ret$beta) - 1
    
    u_change_condition <- (sign(ret$beta[k+1]) < 0)
    
    if(u_change_condition) {
      ret$u <- -ret$u
      ret$beta[k+1] <- -ret$beta[k+1]
    }
  } else {
    k <- length(ret$beta) - 2
    
    u_condition <- ret$beta[k+1] < 0 
    
    if( u_condition ) {
      ret$u <- -ret$u; ret$v <- -ret$v;
      ret$beta[k+1] <- -ret$beta[k+1]
      ret$beta[k+2] <- -ret$beta[k+2]
      
    }
    
  }
  
  ret
}

epsa_hat <- function(ret) {
  n <- length(ret$u)
  if(!is.null(ret$v)) {
    epsa2 <- sum((ret$d*ret$u%*%t(ret$v) - ret$A)^2)/n^2
  } else {
    epsa2 <- sum((ret$d*ret$u%*%t(ret$u) - ret$A)^2)/n^2
  }
  
  sqrt(epsa2)
}


epsa_hat <- function(ret) {
  n <- length(ret$u)
  if(!is.null(ret$v)) {
    epsa2 <- sum((ret$d*ret$u%*%t(ret$v) - ret$A)^2)/n^2
  } else {
    epsa2 <- sum((ret$d*ret$u%*%t(ret$u) - ret$A)^2)/n^2
  }
  
  sqrt(epsa2)
}

epsy_hat <- function(ret) {
  n <- length(ret$residuals)
  p <- length(ret$beta)
  
  sqrt(sum(ret$residuals^2)/(n - p))
}

lopt_estimate <- function(A_ori, X, y, weights){
  
  n <- length(y)
  n_obs <- sum(weights == 1)
  
  ret_two_stage <- two_stage(A_ori, X, y, weights = weights)
  sigmayhat2 <- epsy_hat(ret_two_stage)^2
  sigmaahat2 <- epsa_hat(ret_two_stage)^2
  
  l <- n*sigmayhat2/sigmaahat2*n/n_obs
  
  l
}

#' @title Two-Stage Regression Analysis
#' @description This function implements a two-stage regression analysis. 
#' The first stage performs a truncated singular value decomposition (SVD) of the input matrix \(A\),
#' and the second stage fits an Ordinary Least Squares (OLS) regression model using the SVD components along with 
#' the optional design matrix \(X\).
#'
#' @param A A numeric matrix of size n and p, representing the input data matrix.
#' @param X An optional numeric matrix of size n and k, representing the design matrix. If `NULL`, no covariates are included.
#' @param y A numeric vector of length \(n\), representing the response variable.
#' @param r An integer specifying the rank for the truncated singular value decomposition (default is 1).
#' @param scaled A logical value indicating whether to scale the singular value decomposition components and the regression coefficients (default is `TRUE`).
#' @param weights A numeric vector of length \(n\), representing the weights for the observations. Default is a vector of 1's.
#' @param ... Additional arguments passed to other methods.
#' 
#' @return A list of class `"two_stage"` containing the following components:
#' \item{beta}{A numeric vector of estimated regression coefficients.}
#' \item{u}{A numeric vector of left singular vectors from the SVD of matrix \(A\).}
#' \item{v}{A numeric vector of right singular vectors from the SVD of matrix \(A\).}
#' \item{d}{The first singular value from the truncated SVD of matrix \(A\).}
#' \item{method}{A character string indicating the method used, i.e., "two_stage".}
#' \item{A}{The original input matrix \(A\).}
#' \item{X}{The optional design matrix \(X\).}
#' \item{y}{The response vector \(y\).}
#' \item{u_distance}{Currently set to `NA`. Placeholder for distance diagnostics.}
#' \item{iter}{Currently set to `NA`. Placeholder for the number of iterations.}
#' \item{l}{Currently set to `NA`. Placeholder for regularization parameters.}
#'
#' @importFrom irlba irlba
#' @importFrom stats lm.fit
#' @examples
#' \dontrun{
#' # Generate example data
#' A <- matrix(rnorm(100), nrow = 20, ncol = 5)
#' X <- matrix(rnorm(40), nrow = 20, ncol = 2)
#' y <- rnorm(20)
#' 
#' # Perform two-stage regression
#' result <- two_stage(A, X, y, r = 2)
#' print(result)
#' }
#' 
#' @export
two_stage <- function(A, X, y, r = 1, scaled = 1, weights = rep(1, length(y)), ...) {
  
  n = nrow(A)
  k = ifelse(is.null(X), 0, ncol(X)) 
  n_obs = sum(weights == 1)
  
  # move the unobserved to the end
  A_ori <- A
  
  if(n_obs < n) A <- A_shift_unobserved_last(A_ori, weights)
  
  # SVD
  svda <- irlba::irlba(A, nu = r, nv = r) 
  uu <- svda$u[,1]
  vv <- svda$v[,1]
  d <- svda$d[1]
  
  # OLS
  ret_two_stage <- lm.fit(cbind(X, uu[1:n_obs], vv[1:n_obs]), y)
  ret_two_stage$beta <- ret_two_stage$coefficients
  
  if(scaled) {
    uu <- uu*sqrt(n); vv <- vv*sqrt(n); d <- d/n
    ret_two_stage$beta[k+1] <- ret_two_stage$beta[k+1]/sqrt(n)
    ret_two_stage$beta[k+2] <- ret_two_stage$beta[k+2]/sqrt(n)
  }
  
  if(n_obs < n) {
    uu <- unobserved_shift_back(uu, weights)
    vv <- unobserved_shift_back(vv, weights)
  }
  
  ret_two_stage$u <- uu
  ret_two_stage$v <- vv
  ret_two_stage$d <- d
  
  # ret_two_stage <- adjust_sign(ret_two_stage) 
  ret_two_stage <- fix_sign(ret_two_stage) 
  
  ret_two_stage$u_distance <- NA
  ret_two_stage$method <- "two_stage"
  
  ret_two_stage$A <- A_ori
  ret_two_stage$X <- X
  ret_two_stage$y <- y
  
  #ret_two_stage$epsa <- epsa_hat(ret_two_stage)  
  #ret_two_stage$epsy <- epsy_hat(ret_two_stage)
  
  ret_two_stage$iter <- NA
  ret_two_stage$l <- NA
  
  class(ret_two_stage) <- "two_stage"
  
  ret_two_stage
}



#' @title QrNet Function
#' @description This function implements a method for Quantile Regression using ADMM (Alternating Direction Method of Multipliers).
#' @param A A matrix or data frame representing the design matrix.
#' @param X An optional matrix of covariates, where `ncol(X)` is the number of predictors.
#' @param y A vector of the response variable.
#' @param l Regularization parameter (lambda) for the penalty term. If `NULL`, it is estimated from the data.
#' @param sigma A parameter controlling the strength of the regularization (default is 1).
#' @param tau Quantile level for regression (default is 0.5 for median regression).
#' @param tol Tolerance for convergence in the ADMM optimization algorithm (default is 1e-4).
#' @param max_iter Maximum number of iterations for the ADMM algorithm (default is 200).
#' @param weights A vector of weights for the observations (default is equal weights).
#' @param verbose If non-zero, provides verbose output during function execution.
#' @param ... Additional arguments passed to the lower-level functions.
#' @return A list of results from the quantile regression, including estimated coefficients, residuals, and other diagnostics.
#' @examples
#' \dontrun{
#' # Example of usage with simulated data:
#' A <- matrix(rnorm(100), nrow = 20)
#' X <- matrix(rnorm(40), nrow = 20)
#' y <- rnorm(20)
#' result <- QrNet(A, X, y)
#' print(result) 
#' }
#' @export
QrNet <- function(A, X, y, l = NULL, sigma = 1, tau = 0.5, tol = 1e-4, max_iter = 200, 
                  weights = rep(1, length(y)), verbose = 0, ...) {
  
  n = nrow(A)
  k = ifelse(is.null(X), 0, ncol(X))
  n_obs = sum(weights == 1)
  
  # move the unobserved to the end
  A_ori <- A
  if(n_obs < n) {
    print("Semi")
    A <- A_shift_unobserved_last(A_ori, weights)
  }
  
  if(is.null(l)) l <- lopt_estimate(A_ori, X, y, weights)
  
  ret <- lr(A, X, y, tau, l, sigma, tol, max_iter, verbose, 1) 
  
  # Adjust sign
  # ret <- adjust_sign(ret)
  ret <- fix_sign(ret) 
  
  # shift the order back
  if(n_obs < n) {
    ret$u <- unobserved_shift_back(ret$u, weights)
    ret$v <- unobserved_shift_back(ret$v, weights)
  }
  
  ret$A <- A_ori
  ret$X <- X
  ret$y <- y
  
  ret$epsa <- epsa_hat(ret)
  ret$epsy <- epsy_hat(ret)
  
  class(ret) <- "QrNet"
  
  ret
}




















