test_that("Gradients least-sq coincide", {
  n_obs <- 32; n_pred <- 4
  data <- simulate_data(n_obs, n_pred, model = "linear", seed = 1918)
  design <- data$design; outcome <- data$outcome
  beta <- rep(0, n_pred)
  beta_theory <- grad_linear(design, outcome, beta, noise_var = 1)
  log_lkl <- function(beta, noise_var = 1) {
    log_likelihood_linear(design, outcome, beta, noise_var)
  }
  beta_numerical <- approx_grad(log_lkl, beta)
  expect_equal(as.vector(beta_theory), as.vector(beta_numerical))
})

test_that("Gradients logit coincide", {
  n_obs <- 32; n_pred <- 4
  data <- simulate_data(n_obs, n_pred, model = "logit", seed = 1918)
  design <- data$design; outcome <- data$outcome
  beta <- rep(0, n_pred)
  beta_theory <- grad_logit(design, outcome, beta)
  log_lkl <- function(beta) {
    log_likelihood_logit(design, outcome, beta)
  }
  beta_numerical <- approx_grad(log_lkl, beta)
  expect_equal(as.vector(beta_theory), as.vector(beta_numerical))
})

test_that("linalg and optim least-sq coincide", {
  n_obs <- 32; n_pred <- 4
  data <- simulate_data(n_obs, n_pred, model = "linear", seed = 1918)
  design <- data$design; outcome <- data$outcome
  beta <- rep(0, n_pred)
  via_linalg_out <- hiper_glm(
    design, outcome,
    model = "linear", option = list(mle_solver = "pinv")
  )
  via_bfgs_out <- hiper_glm(
    design, outcome,
    model = "linear", option = list(mle_solver = "BFGS")
  )
  expect_true(are_all_close(
    coef(via_linalg_out), coef(via_bfgs_out),
    abs_tol = 1e-6, rel_tol = 1e-6
  ))
})

test_that("Newton and bfgs outputs coincide on logit model", {
  n_obs <- 32; n_pred <- 4
  data <- simulate_data(n_obs, n_pred, model = 'logit', seed = 1918)
  design <- data$design; outcome <- data$outcome
  via_Newton_out <- hiper_glm(design, outcome, model = 'logit', option = list(mle_solver = 'Newton'))
  via_bfgs_out <- hiper_glm(
    design, outcome, model = 'logit', option = list(mle_solver = 'BFGS')
  )
  expect_true(are_all_close(
    coef(via_Newton_out), coef(via_bfgs_out),
    abs_tol = 1e-2, rel_tol = 1e-2
  ))
})
