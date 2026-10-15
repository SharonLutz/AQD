#' @export
print.aqd <- function(x, ...) {
  cat(sprintf("AQD model selection (%s family, n = %d, %d models)\n\n",
              x$family, x$n, nrow(x$models)))

  cat("--- Nuisance parameters (from full model) ---\n")
  if (length(x$nuisance) == 0) cat("(none)\n") else print(round(x$nuisance, 3))

  cat("\n--- Variable Selection Results ---\n")
  print(x$selection, row.names = FALSE)

  # One method at a time, so long formulas don't break the table apart
  cat("\n--- Selected Model per Method ---\n")
  for (i in seq_len(nrow(x$best))) {
    b <- x$best[i, ]
    if (is.na(b$Model)) {
      cat(sprintf("%-16s no model (every fit failed)\n", b$Method))
    } else {
      cat(sprintf("%-16s k = %-3d criterion = %s\n", b$Method, as.integer(b$k), format(b$Criterion, digits = 6)))
      cat(strwrap(b$Model, width = 76, indent = 4, exdent = 4), sep = "\n")
    }
  }

  cat("\n--- Full-Model P-values ---\n")
  pv <- x$pvalues
  pv[-1] <- lapply(pv[-1], formatC, format = "e", digits = 2)
  print(pv, row.names = FALSE)

  cat(sprintf("\n--- Computational Cost Summary (sum over %d models) ---\n", nrow(x$models)))
  print(x$timing, row.names = FALSE)

  if (any(x$failures[, -1] > 0)) {
    cat("\n--- Failed Fits and Warnings ---\n")
    print(x$failures, row.names = FALSE)
  }

  invisible(x)
}
