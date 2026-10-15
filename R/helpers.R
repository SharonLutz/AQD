# Fits a full model once (section 2): stops with R's reason if the fit fails,
# and warns if it did not converge, since the nuisance estimates use its means
fit_full_model <- function(fit_expr, what) {
  fit <- tryCatch(fit_expr, error = function(e) stop(what, " failed: ", conditionMessage(e), call. = FALSE))
  if (isFALSE(fit$converged))
    warning(what, " did not converge; the nuisance estimates may be unreliable.", call. = FALSE)
  fit
}


# Fits one candidate model (section 4): returns the fit (NULL if it failed)
# and whether it raised any warnings, without printing them
quiet_fit <- function(fit_expr) {
  warned <- FALSE
  fit <- withCallingHandlers(
    tryCatch(fit_expr, error = function(e) NULL),
    warning = function(w) {
      warned <<- TRUE
      invokeRestart("muffleWarning")
    }
  )
  list(fit = fit, warned = warned)
}


# CMP variance lookup: mean and variance on a log-spaced grid of lambda,
# plus finite-difference V'(mu) and V''(mu) (adapted from Simulation 2)
cmp_lookup <- function(nu, mu_max, n_grid = 4000) {

  # Widen the lambda grid until the means cover the data (with a 20% margin)
  lambda_max <- max(50, mu_max * 10)
  for (iter in 1:20) {
    lambda_grid <- exp(seq(log(1e-10), log(lambda_max), length.out = n_grid))
    mu_grid  <- cummax(COMPoissonReg::ecmp(lambda = lambda_grid, nu = nu))  # remove tiny numerical wiggles
    var_grid <- COMPoissonReg::vcmp(lambda = lambda_grid, nu = nu)
    if (max(mu_grid) >= 1.2 * mu_max) break
    lambda_max <- lambda_max * 2
  }

  # First derivative on midpoints, second derivative on midpoints of midpoints
  dV_dmu <- diff(var_grid) / pmax(diff(mu_grid), 1e-14)
  mu1    <- (mu_grid[-1] + mu_grid[-n_grid]) / 2
  d2V    <- diff(dV_dmu) / pmax(diff(mu1), 1e-14)
  mu2    <- (mu1[-1] + mu1[-length(mu1)]) / 2

  list(mu = mu_grid, V = var_grid, mu1 = mu1, Vp = dV_dmu, mu2 = mu2, Vpp = d2V)
}


# Third-order AQD with the expansion point alpha = mu + c * (y - mu) chosen
# to minimize the fourth-order remainder (one c for all observations)
aqd_opt3 <- function(y, mu, wfun, wp, wpp, c_grid, eps) {
  Delta <- y - mu
  best_obj <- Inf
  best_c <- 0.0

  for (c in c_grid) {
    alpha <- pmax(mu + c * Delta, eps)

    # 1. The exact 4th-order remainder polynomial
    poly_c <- (c^2 / 2) - (c / 3) + (1 / 12)
    rem <- (Delta^4) * wpp(alpha) * poly_c

    obj <- sum(rem^2, na.rm = TRUE)
    if (!is.nan(obj) && obj < best_obj) {
      best_obj <- obj
      best_c <- c
    }
  }

  alpha_opt <- pmax(mu + best_c * Delta, eps)
  w_a  <- wfun(alpha_opt)
  wp_a <- wp(alpha_opt)

  d3 <- (Delta^2 * w_a) + (Delta^2 * wp_a * (mu - alpha_opt)) + ((Delta^3) * wp_a / 3)

  sum(d3, na.rm = TRUE)
}

