#' Model selection using the Approximate Quasi-Deviance (AQD)
#'
#' Evaluates every subset of the candidate covariates (exhaustive search) and
#' selects the best mean model under three AQD criteria: AQD(y), AQD(mu) and
#' AQD-opt3. Nuisance parameters (e.g. a dispersion phi, the Tweedie power p,
#' the CMP nu or the negative binomial theta) are estimated once from the full
#' model and held fixed across all candidate models ("anchoring").
#'
#' @param data A data.frame containing the outcome and covariates.
#' @param outcome Name of the outcome column (character).
#' @param covariates Names of the candidate covariate columns (character vector).
#' @param family The mean-variance model. One of
#'   "gaussian" (continuous, constant variance),
#'   "tweedie" (non-negative, right-skewed, zeros allowed),
#'   "cmp" (Conway-Maxwell-Poisson counts; needs the COMPoissonReg package),
#'   "poisson" (quasi-Poisson counts, V = phi * mu),
#'   "negbin" (negative binomial counts, V = mu + mu^2 / theta; uses MASS),
#'   "gamma" (positive, V = phi * mu^2),
#'   "inverse.gaussian" (positive, V = phi * mu^3),
#'   "binomial" (0/1 outcomes, V = mu * (1 - mu)) or
#'   "beta" (proportions strictly between 0 and 1, V = mu * (1 - mu) / (1 + phi)).
#' @param oracle If TRUE, also fit the exact-likelihood model to every subset
#'   for comparison: cplm::cpglm ("tweedie"), COMPoissonReg::glm.cmp ("cmp"),
#'   MASS::glm.nb ("negbin") or betareg::betareg ("beta"). This is slow. For
#'   the other families the baseline AIC is already the exact likelihood (or,
#'   for quasi-Poisson, no likelihood exists).
#' @param penalty Penalty per parameter. 2 gives AIC-type criteria.
#' @param c_grid Grid of expansion-point weights (between 0 and 1) searched by AQD-opt3.
#' @param eps Small positive bound used to keep means and variances away from zero.
#' @param max_covariates Maximum number of candidate covariates (2^p models are fit).
#' @param n_cores Number of cores for evaluating models in parallel (macOS/Linux only).
#'
#' @return An object of class "aqd", a list with components:
#' \describe{
#'   \item{selection}{One row per covariate and one column per method (the
#'     baseline AIC, the oracle if run, AQD(y), AQD(mu), AQD-opt3). "X" marks
#'     covariates in that method's selected model; NA means every fit failed.}
#'   \item{best}{One row per method: the selected model (Model), its number of
#'     mean parameters (k) and its penalized criterion value (Criterion).}
#'   \item{models}{Every candidate model (2^p rows) with columns m (model
#'     number), model (right-hand side), k, aic_base, aqd_y, aqd_mu, aqd_o3
#'     (AQD values without the penalty), aic_oracle (if run), the per-model
#'     times t_base, t_work, t_aqd_y, t_aqd_mu, t_aqd_o3, t_oracle (seconds),
#'     and warn_base, warn_work, warn_oracle (whether that fit gave a warning).
#'     Inf means the fit failed or did not converge.}
#'   \item{nuisance}{The anchored nuisance estimates from the full model.}
#'   \item{pvalues}{Full-model p-values, one row per coefficient.}
#'   \item{timing}{Total seconds per method, summed over all models.}
#'   \item{failures}{Per method: number of failed fits and fits with warnings.}
#'   \item{family, n, call}{The family, the sample size after dropping missing
#'     values, and the function call.}
#' }
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' sim <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
#' sim$y <- 1 + 0.5 * sim$x1 + 0.3 * sim$x2 + rnorm(n)
#' fit <- dev_approximation(sim, "y", c("x1", "x2", "x3"), family = "gaussian")
#' fit
#'
#' @export
dev_approximation <- function(data, outcome, covariates, family,
                              oracle = FALSE,
                              penalty = 2,
                              c_grid = seq(0, 1, by = 0.05),
                              eps = 1e-4,
                              max_covariates = 15,
                              n_cores = 1) {

  families <- c("gaussian", "tweedie", "cmp", "poisson", "negbin",
                "gamma", "inverse.gaussian", "binomial", "beta")
  if (missing(family)) stop("Please choose a family: ", paste0("\"", families, "\"", collapse = ", "))
  family <- match.arg(family, families)

  # 1. Check inputs (returns complete cases of outcome + covariates)
  dat <- check_inputs(data, outcome, covariates, family, oracle, penalty,
                      c_grid, eps, max_covariates, n_cores)

  y_raw <- dat[[outcome]]
  X_df  <- dat[, covariates, drop = FALSE]
  n     <- length(y_raw)

  # Naive baselines on a transformed outcome: log(y + 1) for Tweedie, logit(y) for beta
  y_trans <- y_raw
  if (family == "tweedie") y_trans <- log(y_raw + 1)
  if (family == "beta")    y_trans <- qlogis(y_raw)

  # Data frames for formula parsing
  df_trans <- cbind(y_trans = y_trans, X_df)
  df_raw   <- cbind(y_raw = y_raw, X_df)

  # Baseline model for each family: a GLM, or (gaussian, tweedie, beta) a linear model
  glm_family <- NULL
  if (family %in% c("cmp", "poisson", "negbin")) glm_family <- poisson(link = "log")
  if (family == "gamma")            glm_family <- Gamma(link = "log")
  if (family == "inverse.gaussian") glm_family <- inverse.gaussian(link = "log")
  if (family == "binomial")         glm_family <- binomial(link = "logit")

  # Labels used in the output tables
  base_label <- c(gaussian = "Gaussian AIC", tweedie = "LM AIC", cmp = "Poisson AIC",
                  poisson = "Poisson AIC", negbin = "Poisson AIC", gamma = "Gamma AIC",
                  inverse.gaussian = "Inv. Gaussian AIC", binomial = "Logistic AIC",
                  beta = "LM-logit AIC")[[family]]

  # Exact-likelihood oracle: only where the baseline is not already exact
  oracle_labels <- c(tweedie = "Tweedie Oracle", cmp = "CMP Oracle",
                     negbin = "NB Oracle", beta = "Beta Oracle")
  has_oracle <- oracle && family %in% names(oracle_labels)
  oracle_label <- if (has_oracle) oracle_labels[[family]] else NA
  if (oracle && !has_oracle) {
    reason <- if (family == "poisson") "quasi-Poisson has no likelihood" else "the baseline AIC is already the exact likelihood"
    message(sprintf("No oracle for family = \"%s\" (%s); oracle = TRUE is ignored.", family, reason))
  }

  # Keeps means inside (0, 1) for the binomial and beta variance functions
  clamp01 <- function(u) pmin(pmax(u, eps), 1 - eps)

  # 2. Global nuisance parameters, estimated once from the full model.
  #    Each family sets: nuisance, V_fun / Vp_fun / Vpp_fun, y_tilde (where AQD(y)
  #    evaluates V; NULL if AQD(y) is undefined) and the full fit used for p-values.
  if (family == "gaussian") {
    full_base <- fit_full_model(lm(y_trans ~ ., data = df_trans), "Full linear model")

    # Residual variance from the full model
    sigma2_hat <- sum(residuals(full_base)^2) / full_base$df.residual
    message(sprintf("Estimated Gaussian variance: sigma^2 = %.3f", sigma2_hat))
    nuisance <- c(sigma2 = sigma2_hat)

    # Constant variance, so all three AQD criteria equal RSS / sigma^2
    V_fun   <- function(u) rep(sigma2_hat, length(u))
    Vp_fun  <- function(u) rep(0, length(u))
    Vpp_fun <- function(u) rep(0, length(u))

    y_tilde <- y_raw
    p_fit <- full_base; p_disp <- NULL; p_label <- "P-val (LM)"

  } else if (family == "tweedie") {
    message("Estimating global p and phi via Method of Moments from full quasi-Poisson GLM...")
    full_base <- fit_full_model(lm(y_trans ~ ., data = df_trans), "Full linear model on log(y + 1)")

    # A) Fit the full quasi-Poisson log-link model
    full_glm <- fit_full_model(
      glm(y_raw ~ ., family = quasipoisson(link = "log"), data = df_raw,
          start = coef(full_base)),
      "Full quasi-Poisson GLM"
    )
    mu_full <- pmax(as.numeric(fitted(full_glm)), 1e-8)

    # B) Method of Moments (log-linear regression of squared residuals on fitted means)
    R2 <- (y_raw - mu_full)^2
    Z <- log(pmax(R2, 1e-8))
    U <- log(mu_full)

    S_U  <- sum(U)
    S_U2 <- sum(U^2)
    S_Z  <- sum(Z)
    S_UZ <- sum(U * Z)

    denom <- n * S_U2 - S_U^2
    if (denom < 1e-8) {
      p_hat <- 1.5
      A <- S_Z / n
    } else {
      p_hat <- (n * S_UZ - S_U * S_Z) / denom
      A <- (S_Z - p_hat * S_U) / n
    }

    # Bias correction for log of chi-sq error
    phi_hat <- exp(A + 1.27)

    # Constrain p_hat to strictly valid Tweedie bounds
    p_hat <- min(max(p_hat, 1.05), 1.95)

    message(sprintf("MoM Estimated Tweedie parameters: p = %.3f, phi = %.3f", p_hat, phi_hat))
    nuisance <- c(p = p_hat, phi = phi_hat)

    # Tweedie variance V(u) = phi * u^p and its derivatives
    V_fun   <- function(u) pmax(phi_hat * (pmax(u, eps)^p_hat), eps)
    Vp_fun  <- function(u) p_hat * phi_hat * (pmax(u, eps)^(p_hat - 1))
    Vpp_fun <- function(u) p_hat * (p_hat - 1) * phi_hat * (pmax(u, eps)^(p_hat - 2))

    # AQD(y) evaluates V at y, which must be positive
    y_tilde <- pmax(y_raw, eps)
    p_fit <- full_base; p_disp <- NULL; p_label <- "P-val (LM)"

  } else if (family == "cmp") {
    message("Estimating global nu from the Pearson dispersion of the full Poisson GLM...")
    full_base <- fit_full_model(glm(y_raw ~ ., family = glm_family, data = df_raw), "Full Poisson GLM")
    mu_full <- pmax(as.numeric(fitted(full_base)), eps)

    # CMP variance is roughly mu / nu, so nu ~ 1 / (Pearson dispersion)
    phi_pois <- sum((y_raw - mu_full)^2 / mu_full) / full_base$df.residual
    nu_hat <- 1 / pmax(phi_pois, 1e-12)

    message(sprintf("Estimated CMP dispersion: nu = %.3f", nu_hat))
    nuisance <- c(nu = nu_hat)

    # CMP has no closed-form variance, so V, V' and V'' come from a lookup table in mu
    lookup <- cmp_lookup(nu_hat, mu_max = max(c(y_raw, mu_full, 1)))
    V_fun   <- function(u) pmax(approx(lookup$mu, lookup$V, xout = pmax(u, eps), rule = 2, ties = "ordered")$y, eps)
    Vp_fun  <- function(u) approx(lookup$mu1, lookup$Vp, xout = pmax(u, eps), rule = 2, ties = "ordered")$y
    Vpp_fun <- function(u) approx(lookup$mu2, lookup$Vpp, xout = pmax(u, eps), rule = 2, ties = "ordered")$y

    # AQD(y) evaluates V at y, which must be positive
    y_tilde <- pmax(y_raw, eps)

    # Quasi-Poisson p-values, using the same anchored dispersion as AQD
    p_fit <- full_base; p_disp <- phi_pois; p_label <- "P-val (quasi-Poisson)"

  } else if (family %in% c("poisson", "gamma", "inverse.gaussian")) {
    # Power-variance families: V(u) = phi * u^q, with q = 1 (Poisson), 2 (Gamma),
    # 3 (inverse Gaussian). Tweedie is the same formula with the power estimated.
    q <- c(poisson = 1, gamma = 2, inverse.gaussian = 3)[[family]]

    message("Estimating global phi from the Pearson statistic of the full GLM...")
    full_base <- fit_full_model(glm(y_raw ~ ., family = glm_family, data = df_raw), "Full GLM")
    mu_full <- pmax(as.numeric(fitted(full_base)), eps)

    # Pearson dispersion: phi = sum((y - mu)^2 / mu^q) / (n - k)
    phi_hat <- sum((y_raw - mu_full)^2 / mu_full^q) / full_base$df.residual

    message(sprintf("Estimated dispersion: phi = %.3f", phi_hat))
    nuisance <- c(phi = phi_hat)

    V_fun   <- function(u) pmax(phi_hat * (pmax(u, eps)^q), eps)
    Vp_fun  <- function(u) q * phi_hat * (pmax(u, eps)^(q - 1))
    Vpp_fun <- function(u) q * (q - 1) * phi_hat * (pmax(u, eps)^(q - 2))

    # AQD(y) evaluates V at y, which must be positive (zero counts are clamped)
    y_tilde <- pmax(y_raw, eps)

    p_fit <- full_base; p_disp <- phi_hat
    p_label <- c(poisson = "P-val (quasi-Poisson)", gamma = "P-val (Gamma)",
                 inverse.gaussian = "P-val (Inv. Gaussian)")[[family]]

  } else if (family == "negbin") {
    message("Estimating global theta from the full negative binomial model (MASS::glm.nb)...")
    full_base <- fit_full_model(glm(y_raw ~ ., family = glm_family, data = df_raw), "Full Poisson GLM")
    full_nb   <- fit_full_model(MASS::glm.nb(y_raw ~ ., data = df_raw), "Full negative binomial model")
    theta_hat <- full_nb$theta

    message(sprintf("Estimated negative binomial theta = %.3f", theta_hat))
    nuisance <- c(theta = theta_hat)

    # Negative binomial (NB2) variance V(u) = u + u^2 / theta
    V_fun   <- function(u) pmax(pmax(u, eps) + pmax(u, eps)^2 / theta_hat, eps)
    Vp_fun  <- function(u) 1 + 2 * pmax(u, eps) / theta_hat
    Vpp_fun <- function(u) rep(2 / theta_hat, length(u))

    # AQD(y) evaluates V at y, which must be positive (zero counts are clamped)
    y_tilde <- pmax(y_raw, eps)
    p_fit <- full_nb; p_disp <- NULL; p_label <- "P-val (NB)"

  } else if (family == "binomial") {
    full_base <- fit_full_model(glm(y_raw ~ ., family = glm_family, data = df_raw), "Full logistic regression")

    # Binary outcomes: Var = mu (1 - mu) exactly, so there is no dispersion to estimate
    nuisance <- numeric(0)

    V_fun   <- function(u) clamp01(u) * (1 - clamp01(u))
    Vp_fun  <- function(u) 1 - 2 * clamp01(u)
    Vpp_fun <- function(u) rep(-2, length(u))

    # V(0) = V(1) = 0, so AQD(y) is undefined for 0/1 outcomes. AQD(mu) is the Pearson
    # statistic, which is close to n for every model with binary data, so it has
    # little power to separate models.
    message("AQD(y) is not computed for binary outcomes (V(y) = 0 at y = 0 and 1). ",
            "AQD(mu) is nearly constant for binary data; prefer AQD-opt3.")
    y_tilde <- NULL
    p_fit <- full_base; p_disp <- NULL; p_label <- "P-val (Logistic)"

  } else if (family == "beta") {
    message("Estimating global phi from the Pearson statistic of the full quasi-binomial GLM...")
    full_base <- fit_full_model(lm(y_trans ~ ., data = df_trans), "Full linear model on logit(y)")
    full_qb   <- fit_full_model(glm(y_raw ~ ., family = quasibinomial(link = "logit"), data = df_raw,
                                    start = coef(full_base)),
                                "Full quasi-binomial GLM")
    mu_full <- clamp01(as.numeric(fitted(full_qb)))

    # Beta variance is mu (1 - mu) / (1 + phi), so the quasi-binomial
    # Pearson dispersion s estimates 1 / (1 + phi)
    s_hat <- sum((y_raw - mu_full)^2 / (mu_full * (1 - mu_full))) / full_qb$df.residual
    if (s_hat >= 1) stop("The outcome is too dispersed for a beta model (Pearson dispersion >= 1).")
    phi_hat <- 1 / s_hat - 1

    # phi is a precision here (larger = less variable), unlike the dispersion phi elsewhere
    message(sprintf("Estimated beta precision: phi = %.3f", phi_hat))
    nuisance <- c(precision = phi_hat)

    V_fun   <- function(u) clamp01(u) * (1 - clamp01(u)) / (1 + phi_hat)
    Vp_fun  <- function(u) (1 - 2 * clamp01(u)) / (1 + phi_hat)
    Vpp_fun <- function(u) rep(-2 / (1 + phi_hat), length(u))

    y_tilde <- y_raw
    p_fit <- full_qb; p_disp <- s_hat; p_label <- "P-val (quasi-binomial)"
  }

  # AQD weight w(u) = 1 / V(u) and its derivatives
  wfun <- function(u) 1.0 / V_fun(u)
  wp   <- function(u) -Vp_fun(u) / (V_fun(u)^2)
  wpp  <- function(u) {
    Vv   <- V_fun(u)
    Vpv  <- Vp_fun(u)
    Vppv <- Vpp_fun(u)
    -Vppv / (Vv^2) + 2 * (Vpv^2) / (Vv^3)
  }

  # 3. Model space setup
  p_vars <- length(covariates)
  inc_mat <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), p_vars)))
  colnames(inc_mat) <- covariates
  M <- nrow(inc_mat)

  # 4. Evaluate all models
  evaluate_model <- function(m) {

    inc <- inc_mat[m, ]
    vars <- covariates[inc]

    # Backticks let column names contain spaces or symbols, e.g. `BMI (kg/m2)`
    rhs <- if (length(vars) == 0) "1" else paste0("`", vars, "`", collapse = " + ")
    model_label <- if (length(vars) == 0) "1" else paste(vars, collapse = " + ")
    fml_trans <- as.formula(paste("y_trans ~", rhs))
    fml_raw   <- as.formula(paste("y_raw ~", rhs))

    # Number of mean parameters (defined even if a fit below fails)
    k <- ncol(model.matrix(fml_trans, data = df_trans))

    # --------------------------------------------------
    # A) Baseline: LM (gaussian; tweedie on log(y + 1); beta on logit(y)) or GLM
    # --------------------------------------------------
    # quiet_fit() returns the fit (NULL on error) and whether it gave a warning
    t0_base <- proc.time()[["elapsed"]]
    if (is.null(glm_family)) {
      res_base <- quiet_fit(lm(fml_trans, data = df_trans))
    } else {
      res_base <- quiet_fit(glm(fml_raw, family = glm_family, data = df_raw))
    }
    fit_base <- res_base$fit
    warn_base <- res_base$warned

    # A failed or non-converged fit can never be selected
    base_ok <- !is.null(fit_base) && !isFALSE(fit_base$converged)
    aic_base <- if (base_ok) AIC(fit_base) + (penalty - 2) * k else Inf
    t_base <- proc.time()[["elapsed"]] - t0_base

    # --------------------------------------------------
    # B) AQD (working mean model + scoring)
    # --------------------------------------------------
    aqd_y <- if (is.null(y_tilde)) NA else Inf
    aqd_mu <- Inf; aqd_o3 <- Inf
    t_aqd_y <- 0; t_aqd_mu <- 0; t_aqd_o3 <- 0

    if (family == "tweedie") {
      # Quasi-Poisson log-link GLM, started from the log-scale LM
      t0_work <- proc.time()[["elapsed"]]
      res_work <- quiet_fit(
        glm(fml_raw, family = quasipoisson(link = "log"), data = df_raw, start = coef(fit_base))
      )
      fit_glm <- res_work$fit
      warn_work <- res_work$warned
      t_work <- proc.time()[["elapsed"]] - t0_work

      converged <- !is.null(fit_glm) && fit_glm$converged
      if (converged) muhat <- pmax(as.numeric(fitted(fit_glm)), eps)

    } else if (family == "beta") {
      # Quasi-binomial logit GLM, started from the logit-scale LM
      t0_work <- proc.time()[["elapsed"]]
      res_work <- quiet_fit(
        glm(fml_raw, family = quasibinomial(link = "logit"), data = df_raw, start = coef(fit_base))
      )
      fit_glm <- res_work$fit
      warn_work <- res_work$warned
      t_work <- proc.time()[["elapsed"]] - t0_work

      converged <- !is.null(fit_glm) && fit_glm$converged
      if (converged) muhat <- as.numeric(fitted(fit_glm))

    } else if (family == "gaussian") {
      # The baseline LM is the working model
      t_work <- t_base
      warn_work <- warn_base
      converged <- TRUE
      muhat <- as.numeric(fitted(fit_base))

    } else {
      # The baseline GLM is the working model
      t_work <- t_base
      warn_work <- warn_base
      converged <- base_ok
      if (converged) muhat <- pmax(as.numeric(fitted(fit_base)), eps)
    }

    if (converged) {
      if (!is.null(y_tilde)) {
        t0_y <- proc.time()[["elapsed"]]
        aqd_y <- sum((y_raw - muhat)^2 / V_fun(y_tilde))
        t_aqd_y <- proc.time()[["elapsed"]] - t0_y
      }

      t0_mu <- proc.time()[["elapsed"]]
      aqd_mu <- sum((y_raw - muhat)^2 / V_fun(muhat))
      t_aqd_mu <- proc.time()[["elapsed"]] - t0_mu

      t0_o3 <- proc.time()[["elapsed"]]
      aqd_o3 <- aqd_opt3(y_raw, muhat, wfun, wp, wpp, c_grid, eps)
      t_aqd_o3 <- proc.time()[["elapsed"]] - t0_o3
    }

    # --------------------------------------------------
    # C) Exact-likelihood oracle (optional)
    # --------------------------------------------------
    aic_oracle <- NA; t_oracle <- NA; warn_oracle <- NA
    if (has_oracle) {
      t0_oracle <- proc.time()[["elapsed"]]
      if (family == "tweedie") {
        res_oracle <- quiet_fit(cplm::cpglm(fml_raw, link = "log", data = df_raw))
      } else if (family == "cmp") {
        res_oracle <- quiet_fit(COMPoissonReg::glm.cmp(fml_raw, data = df_raw))
      } else if (family == "negbin") {
        res_oracle <- quiet_fit(MASS::glm.nb(fml_raw, data = df_raw))
      } else if (family == "beta") {
        res_oracle <- quiet_fit(betareg::betareg(fml_raw, data = df_raw))
      }
      fit_oracle <- res_oracle$fit
      warn_oracle <- res_oracle$warned

      if (is.null(fit_oracle)) {
        aic_oracle <- Inf
      } else if (family == "tweedie") {
        # cplm defines its AIC method for the stats4 generic
        aic_oracle <- stats4::AIC(fit_oracle) + (penalty - 2) * k
      } else {
        aic_oracle <- AIC(fit_oracle) + (penalty - 2) * k
      }
      t_oracle <- proc.time()[["elapsed"]] - t0_oracle
    }

    # Return row for this model (including tracked times)
    data.frame(m = m, model = model_label, k = k, aic_base = aic_base, aqd_y = aqd_y,
               aqd_mu = aqd_mu, aqd_o3 = aqd_o3, aic_oracle = aic_oracle,
               t_base = t_base, t_work = t_work, t_aqd_y = t_aqd_y,
               t_aqd_mu = t_aqd_mu, t_aqd_o3 = t_aqd_o3, t_oracle = t_oracle,
               warn_base = warn_base, warn_work = warn_work, warn_oracle = warn_oracle)
  }

  # Forking (mclapply) is not available on Windows
  if (n_cores > 1 && .Platform$OS.type == "windows") {
    message("Parallel evaluation is not supported on Windows; using 1 core.")
    n_cores <- 1
  }

  message(sprintf("Evaluating %d models across %d core(s)...", M, n_cores))
  res_list <- parallel::mclapply(seq_len(M), evaluate_model, mc.cores = n_cores)

  failed <- vapply(res_list, inherits, logical(1), what = "try-error")
  if (any(failed)) stop("Model evaluation failed: ", res_list[[which(failed)[1]]])
  res_df <- do.call(rbind, res_list)

  # 5. Penalize and find best models
  has_aqd_y <- !is.null(y_tilde)

  # Penalized criterion for each method (the AICs already include the penalty)
  criteria <- list()
  criteria[[base_label]] <- res_df$aic_base
  if (has_oracle) criteria[[oracle_label]] <- res_df$aic_oracle
  if (has_aqd_y)  criteria[["AQD(y)"]] <- res_df$aqd_y + penalty * res_df$k
  criteria[["AQD(mu)"]]  <- res_df$aqd_mu + penalty * res_df$k
  criteria[["AQD-opt3"]] <- res_df$aqd_o3 + penalty * res_df$k

  selection <- data.frame(Covariate = covariates, check.names = FALSE, row.names = NULL)
  best <- data.frame()
  for (method in names(criteria)) {
    crit <- criteria[[method]]
    if (all(!is.finite(crit))) {
      # which.min() would silently return the intercept-only model here
      message(sprintf("Every fit failed for %s, so no model is selected.", method))
      selection[[method]] <- NA
      best <- rbind(best, data.frame(Method = method, Model = NA, k = NA, Criterion = NA))
    } else {
      m_best <- which.min(crit)
      selection[[method]] <- ifelse(inc_mat[m_best, ], "X", "")
      best <- rbind(best, data.frame(Method = method, Model = res_df$model[m_best],
                                     k = res_df$k[m_best], Criterion = crit[m_best]))
    }
  }

  # Failed fits (criterion Inf) and fits that gave warnings, per method
  warn_by_method <- function(method) {
    if (method == base_label) return(sum(res_df$warn_base))
    if (has_oracle && method == oracle_label) return(sum(res_df$warn_oracle))
    sum(res_df$warn_work)
  }
  failures <- data.frame(
    Method = names(criteria),
    `Failed fits` = vapply(criteria, function(crit) sum(!is.finite(crit)), numeric(1)),
    `Fits with warnings` = vapply(names(criteria), warn_by_method, numeric(1)),
    check.names = FALSE, row.names = NULL
  )
  if (any(failures[, -1] > 0)) {
    message("Some model fits failed or gave warnings; see the 'failures' component of the result.")
  }

  # 6. P-values from the full models (one row per coefficient)
  coef_names <- names(coef(p_fit))[-1]
  pvalues <- data.frame(Coefficient = coef_names, check.names = FALSE, row.names = NULL)
  pvalues[[p_label]] <- summary(p_fit, dispersion = p_disp)$coefficients[-1, 4]

  # Tweedie: the LM p-values are on the log(y + 1) scale, so also report the
  # quasi-Poisson fit of log E[y], the mean model the AQD criteria select on
  if (family == "tweedie") {
    pvalues[["P-val (quasi-Poisson)"]] <- summary(full_glm)$coefficients[-1, 4]
  }

  if (has_oracle && family == "tweedie") {
    full_tweedie <- tryCatch(cplm::cpglm(y_raw ~ ., link = "log", data = df_raw),
                             error = function(e) NULL)
    if (is.null(full_tweedie)) {
      pvalues[["P-val (Tweedie)"]] <- NA
    } else {
      # cplm's summary method is S4 and prints itself, so look it up and silence it
      cpglm_summary <- methods::selectMethod("summary", "cpglm")
      invisible(capture.output(sum_tw <- cpglm_summary(full_tweedie)))
      pvalues[["P-val (Tweedie)"]] <- sum_tw$coefficients[-1, 4]
    }
  }

  if (has_oracle && family == "cmp") {
    full_cmp <- tryCatch(COMPoissonReg::glm.cmp(y_raw ~ ., data = df_raw),
                         error = function(e) NULL)
    if (is.null(full_cmp)) {
      pvalues[["P-val (CMP)"]] <- NA
    } else {
      # glm.cmp labels the mean-model coefficients "X:<name>" and stores p-values as text
      pvalues[["P-val (CMP)"]] <- as.numeric(summary(full_cmp)$DF[paste0("X:", coef_names), 4])
    }
  }

  if (has_oracle && family == "beta") {
    full_beta <- tryCatch(betareg::betareg(y_raw ~ ., data = df_raw),
                          error = function(e) NULL)
    pvalues[["P-val (Beta)"]] <- if (is.null(full_beta)) NA else
      summary(full_beta)$coefficients$mean[-1, 4]
  }

  # 7. Time summary and results
  time_summary <- data.frame(
    Method = c(base_label, "AQD(y)", "AQD(mu)", "AQD-opt3"),
    `Total Time (sec)` = c(
      sum(res_df$t_base),
      sum(res_df$t_work) + sum(res_df$t_aqd_y),
      sum(res_df$t_work) + sum(res_df$t_aqd_mu),
      sum(res_df$t_work) + sum(res_df$t_aqd_o3)
    ),
    check.names = FALSE
  )
  if (!has_aqd_y) {
    time_summary <- time_summary[time_summary$Method != "AQD(y)", ]
    res_df$aqd_y <- NULL
    res_df$t_aqd_y <- NULL
  }
  if (has_oracle) {
    time_summary <- rbind(time_summary,
                          data.frame(Method = oracle_label, `Total Time (sec)` = sum(res_df$t_oracle),
                                     check.names = FALSE))
  } else {
    res_df$aic_oracle <- NULL
    res_df$t_oracle <- NULL
    res_df$warn_oracle <- NULL
  }

  # Backticks from unusual column names are only needed inside formulas
  pvalues$Coefficient <- gsub("`", "", pvalues$Coefficient)

  structure(
    list(selection = selection, best = best, models = res_df, nuisance = nuisance,
         pvalues = pvalues, timing = time_summary, failures = failures,
         family = family, n = n, call = match.call()),
    class = "aqd"
  )
}

