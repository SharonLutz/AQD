# AQD

Approximate quasi deviance function implementation: model selection with the Approximate Quasi-Deviance (AQD), from "A Class of Approximations of the Quasi-Deviance Function" (Lutz et al., 2026).

`dev_approximation()` fits every subset of a set of candidate covariates and selects the best mean model under three AQD criteria (AQD(y), AQD(mu) and AQD-opt3), alongside a baseline AIC and, optionally, the exact-likelihood AIC. AQD needs only the mean-variance relationship of the outcome.

*Development version (0.0.0.9000).*

## Installation

``` r
install.packages("devtools")   # only if devtools is not installed yet
devtools::install_github("SharonLutz/AQD")
```

Optional packages, needed only for some families or for `oracle = TRUE`:

``` r
install.packages("COMPoissonReg")   # family = "cmp"
install.packages("cplm")            # family = "tweedie", oracle = TRUE
install.packages("betareg")         # family = "beta", oracle = TRUE
```

## Example

``` r
library(AQD)

set.seed(1)
n <- 500
sim <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
sim$y <- 1 + 0.5 * sim$x1 + 0.3 * sim$x2 + rnorm(n)

fit <- dev_approximation(sim, outcome = "y", covariates = c("x1", "x2", "x3"),
                         family = "gaussian")
fit             # selection table, selected models, p-values, timing
fit$best        # selected formula, k and criterion for each method
fit$models      # every candidate model with its criteria
fit$failures    # failed fits / fits with warnings per method

# Top 5 models by AQD-opt3 (the AQD columns exclude the penalty, so add it)
m <- fit$models
head(m[order(m$aqd_o3 + 2 * m$k), c("model", "k", "aqd_o3")], 5)
```

See `?dev_approximation` for all arguments and the full description of the output.

## What it does

1.  Checks the inputs and drops rows with missing values.
2.  Fits the full model once and estimates the nuisance parameters, which are then held fixed for every candidate model ("anchoring").
3.  Fits every subset of the covariates (2\^p models, intercept-only included).
4.  For each model, computes AQD(y), AQD(mu) and AQD-opt3, plus a baseline AIC and, optionally, the exact-likelihood AIC ("oracle").
5.  Picks the model minimizing criterion + penalty × k for each method.

## Families

| family | outcome | link | variance V(μ) | nuisance (full model) | baseline | oracle (optional) |
|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
| `gaussian` | any real | identity | σ² | σ² = RSS/(n−k) | lm AIC (exact) | not needed |
| `tweedie` | ≥ 0, zeros allowed | log | φ μ\^p | p, φ (moments) | lm on log(y+1) | `cplm::cpglm` |
| `cmp` | counts | log | CMP, via lookup table | ν = 1/Pearson dispersion | Poisson AIC | `COMPoissonReg::glm.cmp` |
| `poisson` | counts | log | φ μ | φ = Pearson/(n−k) | Poisson AIC | none (no likelihood) |
| `negbin` | counts | log | μ + μ²/θ | θ (`MASS::glm.nb`, ML) | Poisson AIC | `MASS::glm.nb` |
| `gamma` | \> 0 | log | φ μ² | φ = Pearson/(n−k) | Gamma AIC (exact) | not needed |
| `inverse.gaussian` | \> 0 | log | φ μ³ | φ = Pearson/(n−k) | Inv. Gaussian AIC (exact) | not needed |
| `binomial` | 0/1 (or TRUE/FALSE) | logit | μ(1−μ) | none | logistic AIC (exact) | not needed |
| `beta` | strictly in (0, 1) | logit | μ(1−μ)/(1+φ) | precision φ (quasi-binomial) | lm on logit(y) | `betareg::betareg` |

The working model (which gives the fitted means used by AQD) is a quasi-Poisson GLM for `tweedie`, a quasi-binomial GLM for `beta`, and the baseline fit itself for every other family (a Poisson GLM for `cmp` and `negbin`).

## Notes

-   **AQD(y) and zeros.** AQD(y) evaluates the variance at the observed y, so zero outcomes get very large weights; the function warns. Prefer AQD(mu) or AQD-opt3.
-   **Binary outcomes.** AQD(y) is not computed (V(0) = V(1) = 0), and AQD(mu) is nearly constant across models, so it tends to select too few covariates. Use AQD-opt3.
-   **Input checks.** Rows with missing values are dropped; constant or collinear covariates stop with a message naming them; character covariates are treated as factors; column names may contain spaces or symbols.
-   **Failed fits.** Failed or non-converged fits get criterion Inf (never selected) and are counted in `fit$failures`; per-model warnings are counted, not printed.
-   **Run time.** Exhaustive search grows as 2\^p: 12 covariates = 4,096 models. Use `n_cores` (macOS/Linux) to run models in parallel.
