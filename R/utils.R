#' Test two vectors are close
#'
#' Test two vectors are close in the sense of both relative error and absolute error.
#'
#' @details
#'
#' @param v, vector.
#' @param w, vector of the same length as v.
#' @param abs_tol, number, tolerance for the absolute error.
#' @param rel_tol, number, tolerance for the relative error.
#'
#' @return T&F for the testing.
#'
are_all_close <- function(v, w, abs_tol = 1e-6, rel_tol = 1e-6) {
  abs_diff <- abs(v - w)
  are_all_within_atol <- all(abs_diff < abs_tol)
  are_all_within_rtol <- all(abs_diff < rel_tol * pmax(abs(v), abs(w)))
  return(are_all_within_atol && are_all_within_rtol)
}

#' MLE for linear coefficients
#'
#' Calculate the MLE of linear coefficients using BFGS for general
#' log-likelihood function and gradient.
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#' @param fun function, log-likelihood function taking design, outcome, ... as inputs.
#' @param grad function, gradient of the log-likelihood function.
#'
#' @return Estimate of the linear coefficients
#'
mle_BFGS <- function(design, outcome, fun, grad, ...) {
  beta <- optim(rep(0, ncol(design)), fun, grad,
                design = design, outcome = outcome, ...,
                method = "BFGS",
                control = list(fnscale = -1)
  )$par
  return(beta)
}
