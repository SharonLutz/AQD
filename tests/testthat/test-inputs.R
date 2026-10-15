set.seed(1)
d <- make_X(100)
d$y <- rnorm(100)
d$count <- rpois(100, 3)

test_that("family is required and must be valid", {
  expect_error(dev_approximation(d, "y", "x1"), "Please choose a family")
  expect_error(dev_approximation(d, "y", "x1", family = "weibull"), "should be one of")
})

test_that("missing and misused columns are caught", {
  expect_error(dev_approximation(d, "y", c("x1", "nope"), family = "gaussian"), "not found")
  expect_error(dev_approximation(d, "y", c("x1", "y"), family = "gaussian"), "must not also be")
  expect_error(dev_approximation(d, "y", c("x1", "x1"), family = "gaussian"), "duplicate")
})

test_that("outcome values must fit the family", {
  expect_error(dev_approximation(transform(d, y = y - 10), "y", "x1", family = "tweedie"), "non-negative")
  expect_error(dev_approximation(transform(d, y = abs(y) + 0.5), "y", "x1", family = "poisson"), "whole numbers")
  expect_error(dev_approximation(transform(d, y = y - 10), "y", "x1", family = "gamma"), "positive")
  expect_error(dev_approximation(d, "count", "x1", family = "binomial"), "0/1")
  expect_error(dev_approximation(transform(d, y = abs(y) + 1), "y", "x1", family = "beta"), "between 0 and 1")
})

test_that("arguments are validated", {
  expect_error(dev_approximation(d, "y", c("x1", "x2", "x3"), family = "gaussian", max_covariates = 2), "8 models")
  expect_error(dev_approximation(d, "y", "x1", family = "gaussian", n_cores = 0), "n_cores")
  expect_error(dev_approximation(d, "y", "x1", family = "gaussian", penalty = -1), "penalty")
  expect_error(dev_approximation(d, "y", "x1", family = "gaussian", c_grid = 2), "c_grid")
})

test_that("constant and collinear covariates are caught", {
  expect_error(dev_approximation(transform(d, x4 = 1), "y", c("x1", "x4"), family = "gaussian"),
               "only one value")
  expect_error(dev_approximation(transform(d, x4 = x1 + x2), "y", c("x1", "x2", "x4"), family = "gaussian"),
               "collinear")
})

test_that("zero outcomes trigger the AQD(y) warning", {
  expect_warning(suppressMessages(dev_approximation(d, "count", "x1", family = "poisson")),
                 "outcomes are zero")
})

test_that("missing rows are dropped, character covariates become factors, odd names work", {
  d2 <- d
  d2$y[1:5] <- NA
  d2$sex <- ifelse(d2$x3 > 0, "M", "F")
  names(d2)[names(d2) == "x1"] <- "BMI (kg/m2)"

  expect_message(dev_approximation(d2, "y", c("BMI (kg/m2)", "sex"), family = "gaussian"),
                 "Dropping 5 row")
  fit <- fit_quiet(d2, "y", c("BMI (kg/m2)", "sex"), family = "gaussian")
  expect_equal(fit$n, 95)
  expect_equal(fit$selection$Covariate, c("BMI (kg/m2)", "sex"))
  expect_true("sexM" %in% fit$pvalues$Coefficient)
})
