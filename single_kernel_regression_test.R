single_kernel_regression_test <- function(K, state) {
  y <- K %*% state$alpha + state$b
  
  prediction <- list(y = y)
}