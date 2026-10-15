#' AQD: Model Selection with the Approximate Quasi-Deviance
#'
#' Exhaustive variable selection with approximate quasi-deviance (AQD)
#' criteria, which need only the mean-variance relationship of the outcome.
#' The main function is [dev_approximation()].
#'
#' Based on "A Class of Approximations of the Quasi-Deviance Function"
#' (Lutz et al., Biostatistics, 2026).
#'
#' @keywords internal
#' @importFrom stats AIC as.formula coef complete.cases fitted residuals approx qlogis
#' @importFrom stats lm glm model.matrix
#' @importFrom stats Gamma binomial inverse.gaussian poisson quasibinomial quasipoisson
#' @importFrom utils capture.output
"_PACKAGE"
