#'
#' @export
print.hglm <- function(hglm_out) {
  cat("`hiper_glm` output\n")
}

#' @export
coef.hglm <- function(hglm_out) {
  hglm_out$coefficients
}

#' @export
vcov.hglm <- function(hglm_out) {
  warning("Yet to be implemented.")
}
