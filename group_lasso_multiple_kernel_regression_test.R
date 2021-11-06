group_lasso_multiple_kernel_regression_test <- function(Km, state) {
  Keta <- calculate_Keta(Km, state$eta)
  y <- Keta %*% state$alpha + state$b
  
  prediction <- list(y = y)
}