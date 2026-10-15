# Helpers shared by the tests (testthat loads helper-*.R files automatically)

# n rows of standard normal covariates x1, x2, ..., xp
make_X <- function(n, p = 3) {
  X <- data.frame(matrix(rnorm(n * p), n, p))
  names(X) <- paste0("x", seq_len(p))
  X
}

# Run dev_approximation() without printing progress messages or warnings
fit_quiet <- function(...) suppressWarnings(suppressMessages(dev_approximation(...)))

# Covariates marked "X" by one method in the selection table
selected <- function(fit, method) fit$selection$Covariate[fit$selection[[method]] %in% "X"]
