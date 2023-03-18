#' log-likelihood under linear model
#'
#' Calculate the log-likelihood for the logit model
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#' @param beta vector, the linear coefficients
#'
#' @return log-likelihood under linear model

log_likelihood_logit <- function(design, outcome, beta) {
  sum(outcome * crossprod(t(design), beta) - log(1 + exp(crossprod(t(design), beta))))
}

#' Gradient for the likelihood under linear model
#'
#' Calculate the gradient for the logit model
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#' @param beta vector, the linear coefficients
#'
#' @return Gradient for the likelihood under linear model
#'
grad_logit <- function(design, outcome, beta) {
  p <- outcome_est(design, beta)
  gradient <- as.vector(crossprod(outcome - p, design))
  return(gradient)
}

#' MLE for linear coefficients
#'
#' Calculate the MLE of linear coefficients using BFGS for the logit model
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#'
#' @return Estimate of the linear coefficients.
#'
mle_BFGS_logit <- function(design, outcome) {
  return(
    mle_BFGS(design = design, outcome = outcome,
             fun = log_likelihood_logit, grad = grad_logit)
  )
}

outcome_est <- function(design, beta) {
  1 - 1 / (1 + exp(crossprod(t(design), beta)))
}

update_step <- function(grad, Hessian) {
  step <- solve(Hessian, grad)
  return(step)
}

make_W <- function(design, beta){
  p <- as.vector(1 / (1 + exp(crossprod(t(design), beta))))
  W <- p*(1-p)
  return(W)
}

make_Hessian <- function(design, beta) {
  W = make_W(design, beta)
  Hessian <- -crossprod(design, W * design)
}

take_one_newton_step <-  function(design, outcome, beta, method = 'QR'){
  W = make_W(design, beta)
  p = outcome_est(design, beta)
  X_temp = sqrt(W) * design
  Y_temp = (outcome - p)/sqrt(W)
  if (method == 'QR'){
    step = QR_solve_rcpp(X_temp, Y_temp)
  }
  if (method == 'chol'){
    step = chol_solve(X_temp, Y_temp)
  }
  beta_new <- beta + step

  log_likelihood_new <- log_likelihood_logit(design, outcome, beta_new)
  return(list(
    beta_new = beta_new,
    log_likelihood_new = log_likelihood_new
  ))
}

#' MLE for linear coefficients
#'
#' Calculate the MLE of linear coefficients using pseudo-inverse for the logit model
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#'
#' @return Estimate of the linear coefficients
#'
mle_Newton_logit <- function(design, outcome, max_iter = 1000, method = 'QR') {
  beta_initial <- rep(0, ncol(design))
  beta <- beta_initial
  log_likelihood <- log_likelihood_logit(design, outcome, beta)

  for (i in 1:max_iter) {
    new <- take_one_newton_step(design, outcome, beta, method)
    beta_new <- new$beta_new
    log_likelihood_new <- new$log_likelihood
    if (are_all_close(log_likelihood_new, log_likelihood, abs_tol = 1e-6, rel_tol = 1e-6)) {
      print(paste0("Newton method converges at iteration ", as.character(i)))
      break
    }
    log_likelihood <- log_likelihood_new
    beta <- beta_new
    if (i == max_iter) {
      warning(paste0(
        "Potential non-converging behavior observed with ",
        as.character(max_iter),
        " times of iteration!"
      ))
    }
  }
  return(beta_new)
}
