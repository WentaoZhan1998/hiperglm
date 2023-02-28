#' log-likelihood under linear model
#'
#' Calculate the log-likelihood
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#' @param beta vector, the linear coefficients
#' @param noise_var number, the variance of the noise
#'
#' @return log-likelihood under linear model

log_likelihood <- function(design, outcome, beta, noise_var) {
  .5 * crossprod(outcome - crossprod(t(design), beta))[1, 1] / noise_var
}

#' Gradient for the likelihood under linear model
#'
#' Calculate the gradient
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#' @param beta vector, the linear coefficients
#' @param noise_var number, the variance of the noise
#'
#' @return Gradient for the likelihood under linear model

grad_linear <- function(design, outcome, beta, noise_var) {
  -(crossprod(design, outcome) - crossprod(crossprod(design), beta)) / noise_var
}

#' MLE for linear coefficients
#'
#' Calculate the MLE of linear coefficients using BFGS
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#'
#' @return Estimate of the linear coefficients
#'
mle_BFGS <- function(design, outcome) {
  log_lkl <- function(beta, noise_var = 1) {
    log_likelihood(design, outcome, beta, noise_var)
  }
  grad <- function(beta, noise_var = 1) {
    grad_linear(design, outcome, beta, noise_var)
  }
  beta <- optim(rep(0, ncol(design)), log_lkl, grad)$par
  return(beta)
}

#' MLE for linear coefficients
#'
#' Calculate the MLE of linear coefficients using pseudo-inverse
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#'
#' @return Estimate of the linear coefficients
#'
mle_pinv <- function(design, outcome) {
  solve(crossprod(design), crossprod(design, outcome))
}
