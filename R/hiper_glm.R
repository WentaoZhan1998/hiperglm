#' Gradient for the likelihood under linear model
#'
#' Calculate the gradient
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#' @param beta vector, the potent8ial
#'
#' @return Gradient for the likelihood under linear model

grad_linear = function(design, outcome, beta, noise_var = 1){
  -(crossprod(design, outcome) - crossprod(crossprod(design), beta))/noise_var
}

mle_BFGS = function(outcome, design){
  log_lkl = function(beta, noise_var = 1){
    .5*crossprod(outcome - crossprod(t(design), beta))[1,1]/noise_var
  }
  grad = function(beta, noise_var = 1){
    grad_linear(design, outcome, beta = beta, noise_var = noise_var)
  }
  beta = optim(rep(0, ncol(design)), log_lkl, grad)$par
  return(beta)
}

mle_pinv = function(outcome, design){
  solve(crossprod(design), crossprod(design, outcome))
}

#' Sequence optimize
#'
#' Do sequence optimize
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#' @param model the model to be used. The default is set as 'linear' to use the linear model;
#' otherwise use 'logit' for
#'
#' @return A list including the information of the fitted model
#'
#' @export
hiper_glm <- function(design, outcome, model = "linear", option = NULL) {
  supported_model <- c("linear")
  if (!(model %in% supported_model)) {
    stop(sprintf("The model %s is not supported.", model))
  }
  int = F
  if (model == "linear"){
    if(option$mle_solver == "pinv"){
      beta = mle_pinv(outcome, design)
    }
    else if(option$mle_solver == "BFGS"){
      beta = mle_BFGS(outcome, design)
    }
    offsets = crossprod(t(design), beta)
    res = outcome - offsets
    beta = as.vector(beta)
  }
  hglm_out <- list()
  hglm_out$model = model
  hglm_out$coefficients = beta
  hglm_out$residuals = as.vector(res)
  hglm_out$fitted.values = as.vector(offsets)
  hglm_out$df_res = nrow(design) - ncol(design)
  class(hglm_out) <- "hglm"
  return(hglm_out)
}
