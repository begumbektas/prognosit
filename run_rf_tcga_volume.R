library(randomForestSRC)
source("helper.R")

cohorts <- c("TCGA-THCA") 
data_path <- "./data"
result_path <- "./results"
pathway <- "none" #replace with "hallmark" or with "pid" if you would like to restrict RF to Hallmark genes or to PID genes only

for(cohort in cohorts) {
  if (dir.exists(sprintf("%s/%s", result_path, cohort)) == FALSE) {
    dir.create(sprintf("%s/%s", result_path, cohort)) 
  }
  
  for(replication in 1:100) {
    if (file.exists(sprintf("%s/%s/random_forest_pathway_%s_NRMSE_replication_%d_result.RData", result_path, cohort, pathway, replication)) == FALSE) {
      load(sprintf("%s/%s.RData", data_path, cohort))
      
      TCGA$clinical$tumour_volume <- TCGA$clinical$neoplasm_length * TCGA$clinical$neoplasm_width * TCGA$clinical$neoplasm_depth
      valid_patients <- intersect(rownames(TCGA$clinical)[which((TCGA$clinical$tumour_volume > 0) == TRUE)], rownames(TCGA$mrna))
      
      X <- log2(TCGA$mrna[valid_patients,] + 1)
      valid_features <- as.numeric(which(apply(X, 2, sd) != 0))
      X <- X[,valid_features]
      Y <- TCGA$clinical[valid_patients, c("tumour_volume")]
      Y_new <- Y^(1/3)
      
      train_ratio <- 0.8
      ntree_set <- c(1:5) * 500
      fold_count <- 4
      
      if (pathway != "none") {
        pathways <- read_pathways(pathway)
        gene_names <- sort(unique(unlist(sapply(1:length(pathways), FUN = function(x) {pathways[[x]]$symbols}))))
        X <- X[, which(colnames(X) %in% gene_names)]
      }
      
      set.seed(1606 * replication)
      train_indices <- sample(1:length(Y), ceiling(train_ratio * length(Y)))
      test_indices <- setdiff(1:length(Y), train_indices)
      error_matrix <- matrix(NA, nrow = fold_count, ncol = length(ntree_set), dimnames = list(1:fold_count, sprintf("%g", ntree_set)))
      allocation <- sample(rep(1:fold_count, ceiling(length(train_indices) / fold_count)), length(train_indices))
      
      for (fold in 1:fold_count) {
        train_indices_fold <- train_indices[which(allocation != fold)]
        test_indices_fold <-  setdiff(train_indices, train_indices_fold)
        
        X_train <- X[train_indices_fold,]
        X_test <- X[test_indices_fold,]
        X_train <- scale(X_train)
        X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
        X_train[is.na(X_train)] <- 0
        X_test[is.na(X_test)] <- 0
        
        y_train <- Y_new[train_indices_fold]
        y_test <- Y[test_indices_fold]
        
        train <- as.data.frame(cbind(X_train, y_train))
        test <- as.data.frame(cbind(X_test, y_test)) # Note: Labels of test data is not used. We calculate the error matrix as the following: 
        colnames(test)[length(test)] <- "y_train"
        
        for (ntree in ntree_set){
          rf <- rfsrc(y_train ~ ., data = train, ntree = ntree)
          prediction <- predict(rf, newdata = test)
          error_matrix[fold, sprintf("%g", ntree)] <- NRMSE((prediction$predicted)^3, y_test)
        }
      }
      
      ntree_star <- ntree_set[which.min(colMeans(error_matrix, na.rm = TRUE))]
      
      X_train <- X[train_indices,]
      X_test <- X[test_indices,]
      X_train <- scale(X_train)
      X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
      X_train[is.na(X_train)] <- 0
      X_test[is.na(X_test)] <- 0
      
      y_train <- Y_new[train_indices]
      y_test <- Y[test_indices]
      
      train <- as.data.frame(cbind(X_train, y_train))
      test <- as.data.frame(cbind(X_test, y_test))
      colnames(test)[length(test)] <- "y_train"
      
      rf <- rfsrc(y_train ~ ., data=train, ntree = ntree_star)
      prediction <- predict(rf, newdata = test)
      
      result <- list()
      result$NRMSE <- NRMSE((prediction$predicted)^3, y_test)
      result$predicted <- (prediction$predicted)^3
      result$actual_y <- y_test
      result$ntree_star <- ntree_star
      
      save("result", file = sprintf("%s/%s/random_forest_pathway_%s_NRMSE_replication_%d_result.RData", result_path, cohort, pathway, replication))
    }
  }
}
