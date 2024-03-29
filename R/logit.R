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

#' Response probability
#'
#' Calculate the response probability with the design and coefficient in the logit model.
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param beta vector, the linear coefficients.
#'
#' @return A probability vector for the response.
#'
outcome_est <- function(design, beta) {
  1 - 1 / (1 + exp(crossprod(t(design), beta)))
}

#' Update step for beta
#'
#' Calculate the update step for the linear coefficient beta using Newton's method.
#'
#' @details
#'
#' @param grad vector, the gradient of the logit likelihood.
#' @param Hessian matrix, the Hessian matrix of the likelihood
#'
#' @return A moving step for the linear coefficient beta.
#'
update_step <- function(grad, Hessian) {
  step <- solve(Hessian, grad)
  return(step)
}

#' Hessian matrix calculation
#'
#' Calculate the Hessian matrix for the current beta.
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param beta vector, the linear coefficients.
#'
#' @return The Hessian matrix.
#'
make_Hessian <- function(design, beta) {
  p <- as.vector(1 / (1 + exp(crossprod(t(design), beta))))
  Hessian <- -crossprod(design, p * design)
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
mle_Newton_logit <- function(design, outcome, max_iter = 1000) {
  beta_initial <- rep(0, ncol(design))
  beta <- beta_initial
  log_likelihood <- log_likelihood_logit(design, outcome, beta)

  for (i in 1:max_iter) {
    Hessian <- make_Hessian(design, beta)
    grad = grad_logit(design, outcome, beta)
    step <- update_step(grad, Hessian)
    beta_new <- beta - step
    log_likelihood_new <- log_likelihood_logit(design, outcome, beta_new)
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
