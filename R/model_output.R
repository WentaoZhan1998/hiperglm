#'
#' @export
print.hglm <- function(hglm_out) {
  cat("`hiper_glm` output\n")
}

#' @export
coef.hglm <- function(hglm_out) {
  print(hglm_out$coefficients)
}

#' @export
vcov.hglm <- function(hglm_out) {
  names = names(hglm_out$coefficients)
  if (model == 'linear'){
    vcov = hglm_out$XtX_inv*sum(hglm_out$residuals^2)/hglm_out$df_res
  }

  rownames(vcov) = names
  colnames(vcov) = names
  sprintf('Model: ', hglm_out$model)
  return(vcov)
}
