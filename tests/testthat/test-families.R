# Each family: simulate data where x1 and x2 matter and x3 is noise, then check
# that AQD-opt3 keeps x1 and x2, and that the output has the expected shape.

test_that("gaussian: all AQD criteria agree and match the true covariates", {
  set.seed(11)
  d <- make_X(300)
  d$y <- 1 + 0.6 * d$x1 + 0.4 * d$x2 + rnorm(300)
  fit <- fit_quiet(d, "y", c("x1", "x2", "x3"), family = "gaussian")

  expect_s3_class(fit, "aqd")
  expect_equal(nrow(fit$models), 8)
  expect_named(fit$nuisance, "sigma2")
  expect_equal(fit$selection[["AQD(y)"]], fit$selection[["AQD(mu)"]])
  expect_equal(fit$selection[["AQD(mu)"]], fit$selection[["AQD-opt3"]])
  expect_true(all(c("x1", "x2") %in% selected(fit, "AQD-opt3")))
  expect_equal(fit$best$Method, c("Gaussian AIC", "AQD(y)", "AQD(mu)", "AQD-opt3"))
  expect_equal(sum(fit$failures[["Failed fits"]]), 0)
})

test_that("tweedie: nuisance estimates are sensible and true covariates are kept", {
  set.seed(12)
  d <- make_X(500)
  mu <- exp(0.5 + 0.6 * d$x1 + 0.4 * d$x2)
  # compound Poisson-gamma with p = 1.5, phi = 2
  N <- rpois(500, mu^0.5 / (2 * 0.5))
  d$y <- ifelse(N == 0, 0, rgamma(500, shape = N * 1, scale = 2 * 0.5 * mu^0.5))
  fit <- fit_quiet(d, "y", c("x1", "x2", "x3"), family = "tweedie")

  expect_named(fit$nuisance, c("p", "phi"))
  expect_gt(fit$nuisance[["p"]], 1)
  expect_lt(fit$nuisance[["p"]], 2)
  expect_true(all(c("x1", "x2") %in% selected(fit, "AQD-opt3")))
  expect_true(all(c("P-val (LM)", "P-val (quasi-Poisson)") %in% names(fit$pvalues)))
})

test_that("poisson (quasi): overdispersion is estimated", {
  set.seed(13)
  d <- make_X(300)
  d$y <- rnbinom(300, mu = exp(1 + 0.5 * d$x1 + 0.3 * d$x2), size = 3)
  fit <- fit_quiet(d, "y", c("x1", "x2", "x3"), family = "poisson")

  expect_gt(fit$nuisance[["phi"]], 1)
  expect_true(all(c("x1", "x2") %in% selected(fit, "AQD-opt3")))
})

test_that("negbin: theta is estimated and the oracle runs", {
  set.seed(14)
  d <- make_X(300, p = 2)
  d$y <- rnbinom(300, mu = exp(1 + 0.5 * d$x1 + 0.3 * d$x2), size = 2)
  fit <- fit_quiet(d, "y", c("x1", "x2"), family = "negbin", oracle = TRUE)

  expect_named(fit$nuisance, "theta")
  expect_gt(fit$nuisance[["theta"]], 0.5)
  expect_true("NB Oracle" %in% names(fit$selection))
  expect_equal(selected(fit, "AQD-opt3"), c("x1", "x2"))
})

test_that("gamma and inverse gaussian: dispersion near the truth", {
  set.seed(15)
  d <- make_X(400)
  mu <- exp(0.5 + 0.5 * d$x1 + 0.3 * d$x2)

  d$y <- rgamma(400, shape = 2, scale = mu / 2)                      # phi = 0.5
  fit_g <- fit_quiet(d, "y", c("x1", "x2", "x3"), family = "gamma")
  expect_equal(fit_g$nuisance[["phi"]], 0.5, tolerance = 0.2)
  expect_true(all(c("x1", "x2") %in% selected(fit_g, "AQD-opt3")))

  # inverse Gaussian with phi = 0.2 (Michael-Schucany-Haas sampler)
  lambda <- 1 / 0.2
  v <- rnorm(400)^2
  x <- mu + mu^2 * v / (2 * lambda) - mu / (2 * lambda) * sqrt(4 * mu * lambda * v + mu^2 * v^2)
  d$y <- ifelse(runif(400) <= mu / (mu + x), x, mu^2 / x)
  fit_ig <- fit_quiet(d, "y", c("x1", "x2", "x3"), family = "inverse.gaussian")
  expect_equal(fit_ig$nuisance[["phi"]], 0.2, tolerance = 0.2)
  expect_true(all(c("x1", "x2") %in% selected(fit_ig, "AQD-opt3")))
})

test_that("binomial: no AQD(y) column, no nuisance, AQD-opt3 works", {
  set.seed(16)
  d <- make_X(400)
  d$y <- rbinom(400, 1, plogis(-0.3 + 0.9 * d$x1 + 0.6 * d$x2))
  fit <- fit_quiet(d, "y", c("x1", "x2", "x3"), family = "binomial")

  expect_false("AQD(y)" %in% names(fit$selection))
  expect_false("aqd_y" %in% names(fit$models))
  expect_length(fit$nuisance, 0)
  expect_true(all(c("x1", "x2") %in% selected(fit, "AQD-opt3")))

  # TRUE/FALSE outcomes are accepted
  d$y <- d$y == 1
  expect_s3_class(fit_quiet(d, "y", c("x1", "x2"), family = "binomial"), "aqd")
})

test_that("beta: precision is estimated", {
  set.seed(17)
  d <- make_X(400)
  mu <- plogis(-0.5 + 0.5 * d$x1 + 0.3 * d$x2)
  d$y <- rbeta(400, mu * 10, (1 - mu) * 10)
  fit <- fit_quiet(d, "y", c("x1", "x2", "x3"), family = "beta")

  expect_named(fit$nuisance, "precision")
  expect_equal(fit$nuisance[["precision"]], 10, tolerance = 0.3)
  expect_true(all(c("x1", "x2") %in% selected(fit, "AQD-opt3")))
})

test_that("cmp: nu is estimated (needs COMPoissonReg)", {
  skip_if_not_installed("COMPoissonReg")
  set.seed(18)
  d <- data.frame(x1 = rbinom(200, 1, 0.5), x2 = rbinom(200, 1, 0.5))
  d$y <- COMPoissonReg::rcmp(200, lambda = exp(1 + 0.5 * d$x1 + 0.3 * d$x2), nu = 0.5)
  fit <- fit_quiet(d, "y", c("x1", "x2"), family = "cmp")

  expect_named(fit$nuisance, "nu")
  expect_equal(fit$nuisance[["nu"]], 0.5, tolerance = 0.3)
  expect_equal(selected(fit, "AQD-opt3"), c("x1", "x2"))
})

test_that("optional oracles run when their packages are installed", {
  skip_on_cran()

  set.seed(19)
  d <- make_X(200, p = 2)
  mu <- exp(0.5 + 0.6 * d$x1)
  N <- rpois(200, mu^0.5)
  d$y_tw <- ifelse(N == 0, 0, rgamma(200, shape = N, scale = mu^0.5))
  d$y_beta <- rbeta(200, plogis(d$x1) * 10, (1 - plogis(d$x1)) * 10)

  if (requireNamespace("cplm", quietly = TRUE)) {
    fit <- fit_quiet(d, "y_tw", c("x1", "x2"), family = "tweedie", oracle = TRUE)
    expect_true("Tweedie Oracle" %in% names(fit$selection))
    expect_true("P-val (Tweedie)" %in% names(fit$pvalues))
  }
  if (requireNamespace("betareg", quietly = TRUE)) {
    fit <- fit_quiet(d, "y_beta", c("x1", "x2"), family = "beta", oracle = TRUE)
    expect_true("Beta Oracle" %in% names(fit$selection))
  }
  expect_true(TRUE)  # the test still counts if neither package is installed
})

test_that("oracle = TRUE is ignored with a message where no oracle exists", {
  set.seed(20)
  d <- make_X(100, p = 1)
  d$y <- rnorm(100)
  expect_message(dev_approximation(d, "y", "x1", family = "gaussian", oracle = TRUE), "No oracle")
})
