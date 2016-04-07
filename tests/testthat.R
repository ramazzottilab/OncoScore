Sys.setenv("R_TESTS" = "")

library(testthat)
library(OncoScore)

test_check("OncoScore")
