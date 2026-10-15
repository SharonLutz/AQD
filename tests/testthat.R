# This file runs the package tests during R CMD check.
# Run them interactively with devtools::test().

library(testthat)
library(AQD)

test_check("AQD")
