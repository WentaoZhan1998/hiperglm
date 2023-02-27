test_that("Vectors are all close", {
  x = rep(1, 100)
  y = rep(1, 100)
  expect_true(are_all_close(
    x, y, abs_tol = 1e-2, rel_tol = 1e-2
  ))
})

test_that("Too large absolute error", {
  x = rep(1, 100)
  y = rep(1, 100) + 0.1
  expect_false(are_all_close(
    x, y, abs_tol = 1e-2, rel_tol = 1e16
  ))
})

test_that("Too large relative error", {
  x = rep(1e-2, 100)
  y = rep(1e-2, 100) + 1e-3
  expect_false(are_all_close(
    x, y, abs_tol = 1e16, rel_tol = 1e-2
  ))
})
