source("helper.R")
source("solve_svr_cplex.R")
source("group_lasso_multiple_kernel_regression_train.R")
source("group_lasso_multiple_kernel_regression_test.R")
library(tidyr)

cohorts <- c("TCGA-THCA") 
data_path <- "./data"
result_path <- "./results"
pathway <- "hallmark" #replace with "pid" if you would like to use PID pathways

for(cohort in cohorts) {
  if (dir.exists(sprintf("%s/%s", result_path, cohort)) == FALSE) {
    dir.create(sprintf("%s/%s", result_path, cohort))
  }
  
  for(replication in 1:100) {
    if (file.exists(sprintf("%s/%s/glmkl_pathway_%s_measure_NRMSE_replication_%d_result.RData", result_path, cohort, pathway, replication)) == FALSE) {
      load(sprintf("%s/%s.RData", data_path, cohort))
      
      TCGA$clinical$tumour_volume <- TCGA$clinical$neoplasm_length * TCGA$clinical$neoplasm_width * TCGA$clinical$neoplasm_depth
      valid_patients <- intersect(rownames(TCGA$clinical)[which((TCGA$clinical$tumour_volume > 0) == TRUE)], rownames(TCGA$mrna))
      
      X <- log2(TCGA$mrna[valid_patients,] + 1)
      valid_features <- as.numeric(which(apply(X, 2, sd) != 0))
      X <- X[,valid_features]
      Y <- TCGA$clinical[valid_patients, c("tumour_volume")]
      Y_new <- Y^(1/3)
      
      C_set <- c(0.001, 0.01, 0.1, 1, 10, 100, 1000)
      epsilon <- 1e-5
      fold_count <- 4
      train_ratio <- 0.8
      iteration_count <- 200
      tube_set <- c(0, 1/4, 1/2, 1, 2)
      
      C_set_tube_set <- as.data.frame(crossing(C_set, tube_set))
      C_set_tube_set_list <- function(C_set, tube_set){
        C_set_tube_set_list <- list()
        i <- 0
        while(i < nrow(C_set_tube_set)){
          i <- i + 1
          C_set_tube_set_list[i] <- sprintf("%g_%g", C_set_tube_set[i, 1], C_set_tube_set[i, 2])
        }
        return(unlist(C_set_tube_set_list))
      }
      
      C_set_tube_set <- C_set_tube_set_list(C_set = C_set, tube_set = tube_set)
      
      pathways <- read_pathways(pathway)
      gene_names <- sort(unique(unlist(sapply(1:length(pathways), FUN = function(x) {pathways[[x]]$symbols}))))
      X <- X[, which(colnames(X) %in% gene_names)]
      
      set.seed(1606 * replication)
      train_indices <- sample(1:length(Y), ceiling(train_ratio * length(Y)))
      test_indices <- setdiff(1:length(Y), train_indices)
      error_matrix <- matrix(NA, nrow = fold_count, ncol = length(C_set_tube_set), dimnames = list(1:fold_count, sprintf("%s", C_set_tube_set)))
      allocation <- sample(rep(1:fold_count, ceiling(length(train_indices) / fold_count)), length(train_indices))
      
      for (fold in 1:fold_count) {
        
        train_indices_fold <- train_indices[which(allocation != fold)]
        test_indices_fold <-  setdiff(train_indices, train_indices_fold)
        
        X_train <- X[train_indices_fold,]
        X_test <- X[test_indices_fold,]
        X_train <- scale(X_train)
        X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
        
        N_train <- nrow(X_train)
        N_test <- nrow(X_test)
        N_pathway <- length(pathways)
        K_train <- array(0, dim = c(N_train, N_train, N_pathway))
        K_test <- array(0, dim = c(N_test, N_train, N_pathway))
        
        for (m in 1:N_pathway) {
          feature_indices <- which(colnames(X_train) %in% pathways[[m]]$symbols)
          D_train <- pdist(X_train[, feature_indices], X_train[, feature_indices])
          D_test <- pdist(X_test[, feature_indices], X_train[, feature_indices])
          sigma <- mean(D_train)
          K_train[,,m] <- exp(-D_train^2 / (2 * sigma^2))
          K_test[,,m] <- exp(-D_test^2 / (2 * sigma^2))
        }
        
        y_train <- Y_new[train_indices_fold]
        y_test <- Y[test_indices_fold]
        
        for (C in C_set) {
          for (tube in tube_set){
            print(sprintf("running fold = %d, C = %g, tube = %g", fold, C, tube * sd(y_train)))
            parameters <- list()
            parameters$C <- C
            parameters$epsilon <- epsilon
            parameters$tube <- tube * sd(y_train)
            parameters$iteration_count <- iteration_count
            
            state <- group_lasso_multiple_kernel_regression_train(K_train, y_train, parameters)
            prediction_test <- group_lasso_multiple_kernel_regression_test(K_test, state)
            
            error_matrix[fold, sprintf("%g_%g", C, tube)] <- NRMSE((prediction_test$y)^3, y_test)
          }
        }
      }
      
      C_star_tube_star_NRMSE <- C_set_tube_set[which.min(t(colMeans(error_matrix)))]    
      C_star_NRMSE <- as.numeric(sub("_.*", "", C_star_tube_star_NRMSE))
      tube_star_NRMSE_multiplier <- as.numeric(sub('.*\\_', '', C_star_tube_star_NRMSE)) 
      tube_star_NRMSE <- tube_star_NRMSE_multiplier * sd(Y_new[train_indices])
      
      X_train <- X[train_indices,]
      X_test <- X[test_indices,]
      X_train <- scale(X_train)
      X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
      
      N_train <- nrow(X_train)
      N_test <- nrow(X_test)
      N_pathway <- length(pathways)
      K_train <- array(0, dim = c(N_train, N_train, N_pathway))
      K_test <- array(0, dim = c(N_test, N_train, N_pathway))
      
      for (m in 1:N_pathway) {
        feature_indices <- which(colnames(X_train) %in% pathways[[m]]$symbols)
        D_train <- pdist(X_train[, feature_indices], X_train[, feature_indices])
        D_test <- pdist(X_test[, feature_indices], X_train[, feature_indices])
        sigma <- mean(D_train)
        K_train[,,m] <- exp(-D_train^2 / (2 * sigma^2))
        K_test[,,m] <- exp(-D_test^2 / (2 * sigma^2))
      }
      
      y_train <- Y_new[train_indices]
      y_test <- Y[test_indices]
      
      parameters <- list()
      parameters$C <- C_star_NRMSE
      parameters$epsilon <- epsilon
      parameters$tube <- tube_star_NRMSE
      parameters$iteration_count <- iteration_count
      
      state <- group_lasso_multiple_kernel_regression_train(K_train, y_train, parameters)
      prediction_test <- group_lasso_multiple_kernel_regression_test(K_test, state)
      prediction_train <- group_lasso_multiple_kernel_regression_test(K_train, state)
      
      result <- list()
      result$train_NRMSE <- NRMSE((prediction_train$y)^3, Y[train_indices])
      result$predicted <- (prediction_test$y)^3
      result$actual_y <- y_test
      result$test_NRMSE <- NRMSE((prediction_test$y)^3, y_test)
      result$tube_multiplier <- tube_star_NRMSE_multiplier
      result$C <- C_star_NRMSE
      
      save("state", file = sprintf("%s/%s/glmkl_pathway_%s_measure_NRMSE_replication_%d_state.RData", result_path, cohort, pathway, replication))
      save("result", file = sprintf("%s/%s/glmkl_pathway_%s_measure_NRMSE_replication_%d_result.RData", result_path, cohort, pathway, replication))
    }
  }
}
