group_lasso_multiple_kernel_regression_train <- function(Km, y, parameters) {
  P <- dim(Km)[3]
  eta <- rep(1 / P, P)
  Keta <- calculate_Keta(Km, eta)
  
  model <- solve_svr(Keta, y, parameters$C, parameters$tube, parameters$epsilon)
  objectives <- model$objective
  print(c(length(objectives), model$objective))
  while(1) {
    start_objective <- model$objective
    start_eta <- eta
    
    for (m in 1:P) {
      eta[m] <- eta[m] * sqrt(t(model$alpha) %*% Km[,,m] %*% model$alpha)
    }
    eta <- eta / sum(eta)
    eta[eta < parameters$epsilon] <- 0
    eta <- eta / sum(eta)
    Keta <- calculate_Keta(Km, eta)
    
    model <- solve_svr(Keta, y, parameters$C, parameters$tube, parameters$epsilon)
    
    objectives <- c(objectives, model$objective)
    print(c(length(objectives), model$objective))
    
    if (length(objectives) == parameters$iteration_count) {
      break
    }
  }
  
  state <- list(alpha = model$alpha, b = model$b, eta = eta, objectives = objectives, parameters = parameters)
}