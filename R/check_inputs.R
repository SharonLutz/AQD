# Checks user inputs and returns the complete cases of outcome + covariates
check_inputs <- function(data, outcome, covariates, family, oracle, penalty,
                         c_grid, eps, max_covariates, n_cores) {

  if (!is.data.frame(data)) stop("'data' must be a data.frame.")
  if (!is.character(outcome) || length(outcome) != 1)
    stop("'outcome' must be a single column name.")
  if (!is.character(covariates) || length(covariates) == 0)
    stop("'covariates' must be a character vector of column names.")

  missing_cols <- setdiff(c(outcome, covariates), names(data))
  if (length(missing_cols) > 0)
    stop("Column(s) not found in 'data': ", paste(missing_cols, collapse = ", "))
  if (outcome %in% covariates) stop("'outcome' must not also be in 'covariates'.")
  if (anyDuplicated(covariates)) stop("'covariates' contains duplicate names.")
  if (any(covariates %in% c("y_raw", "y_trans")))
    stop("Please rename covariates called 'y_raw' or 'y_trans' (these names are used internally).")

  if (length(covariates) > max_covariates)
    stop(sprintf("%d covariates means %s models. Use fewer covariates, or raise 'max_covariates' if you really want this.",
                 length(covariates), format(2^length(covariates), big.mark = ",")))

  # Drop rows with missing values
  dat <- data[, c(outcome, covariates), drop = FALSE]
  complete <- complete.cases(dat)
  if (sum(!complete) > 0) {
    message(sprintf("Dropping %d row(s) with missing values.", sum(!complete)))
    dat <- dat[complete, , drop = FALSE]
  }

  y <- dat[[outcome]]
  if (family == "binomial" && is.logical(y)) y <- as.numeric(y)  # TRUE/FALSE -> 1/0
  if (!is.numeric(y)) stop("The outcome must be numeric.")
  dat[[outcome]] <- y <- as.vector(y)  # e.g. a one-column matrix becomes a plain vector
  # Character covariates become factors; unused factor levels are dropped
  is_chr <- vapply(dat[covariates], is.character, logical(1))
  if (any(is_chr)) {
    message("Treating character covariate(s) as factors: ", paste(covariates[is_chr], collapse = ", "))
    dat[covariates[is_chr]] <- lapply(dat[covariates[is_chr]], factor)
  }
  dat[covariates] <- lapply(dat[covariates], function(x) if (is.factor(x)) droplevels(x) else x)

  ok_type <- vapply(dat[covariates], function(x) is.numeric(x) || is.factor(x) || is.logical(x), logical(1))
  if (!all(ok_type))
    stop("Covariates must be numeric, logical, factor or character: ",
         paste(covariates[!ok_type], collapse = ", "))

  # Constant or collinear covariates would give NA coefficients in the fits
  is_constant <- vapply(dat[covariates], function(x) length(unique(x)) < 2, logical(1))
  if (any(is_constant))
    stop("Covariate(s) with only one value (after dropping missing rows): ",
         paste(covariates[is_constant], collapse = ", "))

  X_full <- model.matrix(~ ., data = dat[covariates])
  n_params <- ncol(X_full)
  if (nrow(dat) <= n_params)
    stop(sprintf("Not enough observations (%d) for the full model (%d parameters).", nrow(dat), n_params))

  qr_full <- qr(X_full)
  if (qr_full$rank < n_params) {
    aliased <- colnames(X_full)[qr_full$pivot[(qr_full$rank + 1):n_params]]
    stop("Covariates are collinear (exact linear combinations of other covariates): ",
         paste(gsub("`", "", aliased), collapse = ", "))
  }

  # Outcome values allowed by each family
  if (family == "tweedie") {
    if (any(y < 0)) stop("Tweedie outcomes must be non-negative.")
    if (all(y == 0)) stop("The outcome is all zeros.")
  }
  if (family %in% c("cmp", "poisson", "negbin")) {
    if (any(y < 0 | y != round(y))) stop("Count outcomes must be non-negative whole numbers.")
    if (all(y == 0)) stop("The outcome is all zeros.")
  }
  if (family %in% c("gamma", "inverse.gaussian")) {
    if (any(y <= 0)) stop(sprintf("Outcomes for family = \"%s\" must be positive.", family))
  }
  if (family == "binomial") {
    if (!all(y %in% c(0, 1))) stop("Binomial outcomes must be 0/1 (or TRUE/FALSE).")
    if (length(unique(y)) < 2) stop("The binomial outcome has only one value.")
  }
  if (family == "beta") {
    if (any(y <= 0 | y >= 1)) stop("Beta outcomes must be strictly between 0 and 1.")
  }

  if (!is.numeric(penalty) || length(penalty) != 1 || penalty < 0)
    stop("'penalty' must be a single non-negative number.")
  if (!is.numeric(c_grid) || length(c_grid) == 0 || any(c_grid < 0 | c_grid > 1))
    stop("'c_grid' must be numbers between 0 and 1.")
  if (!is.numeric(eps) || length(eps) != 1 || eps <= 0)
    stop("'eps' must be a single positive number.")
  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores < 1 || n_cores != round(n_cores))
    stop("'n_cores' must be a positive whole number.")
  if (!is.logical(oracle) || length(oracle) != 1 || is.na(oracle))
    stop("'oracle' must be TRUE or FALSE.")

  # Packages needed by some families
  if (family == "cmp" && !requireNamespace("COMPoissonReg", quietly = TRUE))
    stop("family = \"cmp\" requires the 'COMPoissonReg' package. Install it with install.packages(\"COMPoissonReg\").")
  if (family == "negbin" && !requireNamespace("MASS", quietly = TRUE))
    stop("family = \"negbin\" requires the 'MASS' package. Install it with install.packages(\"MASS\").")
  if (oracle && family == "tweedie" && !requireNamespace("cplm", quietly = TRUE))
    stop("oracle = TRUE requires the 'cplm' package. Install it with install.packages(\"cplm\").")
  if (oracle && family == "beta" && !requireNamespace("betareg", quietly = TRUE))
    stop("oracle = TRUE requires the 'betareg' package. Install it with install.packages(\"betareg\").")

  if (family %in% c("tweedie", "cmp", "poisson", "negbin") && any(y == 0))
    warning(sprintf("%.0f%% of outcomes are zero. AQD(y) evaluates the variance at y, so it is unreliable here; prefer AQD(mu) or AQD-opt3.",
                    100 * mean(y == 0)), call. = FALSE)

  dat
}

