#' Sequence optimize
#'
#' Do sequence optimize
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#' @param model the model to be used. The default is set as 'linear' to use the linear model.
#' @param option a list specifies the options for model solving, currently supporting 'mle_solver' only.
#'
#' @return A list including the information of the fitted model
#'
#' @export
hiper_glm <- function(design, outcome, model = "linear", option = NULL) {
  supported_model <- c("linear", "logit")
  if (!(model %in% supported_model)) {
    stop(sprintf("The model %s is not supported.", model))
  }
  if (model == "linear") {
    if (option$mle_solver == "pinv") {
      beta <- mle_pinv(design, outcome)
    }
    if (option$mle_solver == "BFGS") {
      beta <- mle_BFGS(design, outcome)
    }
    offsets <- crossprod(t(design), beta)
    res <- outcome - offsets
    beta <- as.vector(beta)
  }
  if (model == "logit") {
    if (option$mle_solver == "Newton") {
      beta <- mle_Newton_logit(design, outcome)
    }
    if (option$mle_solver == "BFGS") {
      beta <- mle_BFGS_logit(design, outcome)
    }
    offsets <- Outcome_est(design, beta)
    res <- outcome - offsets
    beta <- as.vector(beta)
  }
  hglm_out <- list()
  hglm_out$model <- model
  hglm_out$coefficients <- beta
  hglm_out$residuals <- as.vector(res)
  hglm_out$fitted.values <- as.vector(offsets)
  hglm_out$df_res <- nrow(design) - ncol(design)
  class(hglm_out) <- "hglm"
  return(hglm_out)
}
