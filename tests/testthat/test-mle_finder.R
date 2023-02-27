n_obs <- 32; n_pred <- 4
data <- simulate_data(n_obs, n_pred, model = 'linear', seed = 1918)
design <- data$design; outcome <- data$outcome
set.seed(1)
beta = rnorm(ncol(design))

test_that("Gradients least-sq coincide", {
  beta_theory = grad_linear(design, outcome, beta, noise_var = 1)
  log_lkl = function(beta, noise_var = 1){
    .5*crossprod(outcome - crossprod(t(design), beta))[1,1]/noise_var
  }
  beta_numerical = approx_grad(log_lkl, beta)
  expect_equal(as.vector(beta_theory), as.vector(beta_numerical))
})

test_that("linalg and optim least-sq coincide", {
  via_linalg_out <- hiper_glm(design, outcome, model = 'linear')
  via_bfgs_out <- hiper_glm(
    design, outcome, model = 'linear', option = list(mle_solver = 'BFGS')
  )
  expect_true(are_all_close(
    coef(via_linalg_out), coef(via_bfgs_out), abs_tol = 1e-2, rel_tol = 1e-2
  ))
})
