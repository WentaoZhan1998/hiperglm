#' log-likelihood under linear model
#'
#' Calculate the log-likelihood for the linear model
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#' @param beta vector, the linear coefficients
#' @param noise_var number, the variance of the noise
#'
#' @return log-likelihood under linear model

log_likelihood_linear <- function(design, outcome, beta, noise_var) {
  -.5 * crossprod(outcome - crossprod(t(design), beta))[1, 1] / noise_var
  #-sum((outcome - design %*% beta)^2) / 2 / noise_var
}

#' Gradient for the likelihood under linear model
#'
#' Calculate the gradient for the linear model
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
  (crossprod(design, outcome) - crossprod(crossprod(design), beta)) / noise_var
}

#' MLE for linear coefficients
#'
#' Calculate the MLE of linear coefficients using BFGS for the linear model
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#'
#' @return Estimate of the linear coefficients
#'
mle_BFGS_linear <- function(design, outcome) {
  return(
    mle_BFGS(design = design, outcome = outcome,
             fun = log_likelihood_linear, grad = grad_linear,
             noise_var = 1)
    )
}

#' MLE for linear coefficients
#'
#' Calculate the MLE of linear coefficients using pseudo-inverse for the linear model
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#'
#' @return Estimate of the linear coefficients
#'
mle_pinv <- function(design, outcome) {
  L <- chol(crossprod(design))
  beta <- backsolve(L, forwardsolve(t(L), crossprod(design, outcome)))
  return(beta)
}
