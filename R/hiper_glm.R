#' Sequence optimize
#'
#' Do sequence optimize
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#' @param model the model to be used. The default is set as 'linear' to use the linear model.
#'
#' @return A list including the information of the fitted model
#'
#' @export
hiper_glm <- function(design, outcome, model = "linear", option = NULL) {
  supported_model <- c("linear")
  if (!(model %in% supported_model)) {
    stop(sprintf("The model %s is not supported.", model))
  }
  int <- F
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
  hglm_out <- list()
  hglm_out$model <- model
  hglm_out$coefficients <- beta
  hglm_out$residuals <- as.vector(res)
  hglm_out$fitted.values <- as.vector(offsets)
  hglm_out$df_res <- nrow(design) - ncol(design)
  class(hglm_out) <- "hglm"
  return(hglm_out)
}
