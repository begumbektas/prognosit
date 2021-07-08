single_kernel_regression_train <- function(K, y, parameters) {
  model <- solve_svr(K, y, parameters$C, parameters$tube, parameters$epsilon)
  
  state <- list(alpha = model$alpha, b = model$b, parameters = parameters)
}