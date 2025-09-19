library(caret)
library(kernlab)
library(pracma)


source("twostage_dataset.R")
source("minimax_core.R")

# Define the MinimaxRKHS class
MinimaxRKHS <- setRefClass(
  "MinimaxRKHS",
  fields = list(
    pidim = "numeric",
    hdim = "numeric",
    lambda_pi = "numeric",
    lambda_h = "numeric",
    lambda_score = "numeric",
    gamma = "numeric",
    gamma_score = "numeric",
    Kpi_ = "matrix",
    Kh_ = "matrix",
    alpha_ = "numeric",
    beta_ = "numeric",
    rpi_ = "matrix",
    rh_ = "matrix"
  ),
  
  methods = list(
    initialize = function(pidim, hdim, lambda_pi = 1, lambda_h = 1, lambda_score = 1, gamma = 0.1, gamma_score = 0.1) {
      .self$pidim <- pidim
      .self$hdim <- hdim
      .self$lambda_pi <- lambda_pi
      .self$lambda_h <- lambda_h
      .self$lambda_score <- lambda_score
      .self$gamma <- gamma
      .self$gamma_score <- gamma_score
    },
    
    .create_kernel_pi = function(X, Y = NULL) {
      kernelMatrix(kernel = rbfdot(sigma = .self$gamma), 
                   x = X, 
                   y = Y)
    },
    
    .create_kernel_h = function(X, Y = NULL) {
      kernelMatrix(kernel = rbfdot(sigma = .self$gamma), 
                   x = X, 
                   y = Y)
    },
    
    .unpack_data = function(data) {
      rpi <- data[, 1:.self$pidim]
      rh <- data[, (.self$pidim + 1):(.self$pidim + .self$hdim)]
      g1 <- data[, .self$pidim + .self$hdim + 1]
      g2 <- data[, .self$pidim + .self$hdim + 2]
      list(rpi = rpi, rh = rh, g1 = g1, g2 = g2)
    },
    
    fit = function(data) {
      unpacked <- .self$.unpack_data(data)
      rpi <- unpacked$rpi
      rh <- unpacked$rh
      g1 <- unpacked$g1
      g2 <- unpacked$g2
      
      .self$Kpi_ <- .self$.create_kernel_pi(rpi)
      .self$Kh_ <- .self$.create_kernel_h(rh)
      
      # Solve minimax problem
      solution <- minimax_solve(.self$Kpi_, .self$Kh_, g1, g2, .self$lambda_pi, .self$lambda_h)
      .self$alpha_ <- solution$alpha
      .self$beta_ <- solution$beta
      
      .self$rpi_ <- as.matrix(rpi)
      .self$rh_ <- as.matrix(rh)
      
      .self
    },
    
    predict = function(data) {
      if (is.null(.self$alpha_) | is.null(.self$beta_)) {
        stop("Model is not fitted.")
      }
      rpi <- .self$.unpack_data(data)$rpi
      .self$pi_(rpi)
    },
    
    pi_ = function(r){
      .self$alpha_ %*% .self$.create_kernel_pi(.self$rpi_, r)
    },
    
    h_ = function(r){
      .self$beta_ %*% .self$.create_kernel_h(.self$rh_, r)
    },
    
    score = function(data) {
      
      unpacked <- .self$.unpack_data(data)
      rpi <- unpacked$rpi
      rh <- unpacked$rh
      g1 <- unpacked$g1
      g2 <- unpacked$g2
      gammas <- .self$gamma_score
      
      pi_val <- .self$pi_(rpi)
      
      scrs <- rep(NA,length(gammas))
      for(i in 1:length(gammas)){
        Kh_score <- kernelMatrix(kernel = rbfdot(sigma = gammas[i]),
                                 x = rh,
                                 y = NULL)
        
        scrs[i] <- score_nuisance_function(pi_val, Kh_score, g1, g2, .self$lambda_score)
      }
      
      mean(scrs)
      
    }
  )
)


library(caret)

# Define the custom model

MinimaxRKHSCV <- function(data, pidim, hdim, lambdaspi = NULL, lambdash = NULL, gammas = NULL, split_ratio = 0.25) {
  if (is.null(lambdaspi)) {
    lambdaspi <- 10^seq(-5, 0)
  }
  if (is.null(lambdash)) {
    lambdash <- 10^seq(-5, 0)
  }
  if (is.null(gammas)) {
    gammas <- 1
  }
  
  tune_grid <- expand.grid(lambda_pi = lambdaspi, lambda_h = lambdash, gamma = gammas)
  
  minimax_models <- list()
  scores <- rep(NA,nrow(tune_grid))
  
  # initialize models
  for(i in 1:nrow(tune_grid)){
    
    train_model <- MinimaxRKHS$new(pidim = pidim, hdim = hdim, 
                                   lambda_pi = tune_grid$lambda_pi[i],
                                   lambda_h = tune_grid$lambda_h[i],
                                   lambda_score = min(lambdash),
                                   gamma = tune_grid$gamma[i],
                                   gamma_score = gammas)
    minimax_models[[i]] <- train_model
  }
  
  split_idx <- createDataPartition(1:nrow(data), times = 1, p = split_ratio, list = TRUE)
  scores <- rep(NA, nrow(tune_grid))
  if(split_ratio==1){
    for(i in 1:nrow(tune_grid)){
      minimax_models[[i]]$fit(data[split_idx[[1]],])
      # train_model$predict(data[flds[[k]],])
      scores[i] <- minimax_models[[i]]$score(data[split_idx[[1]],])
    }
  }
  else{
    for(i in 1:nrow(tune_grid)){
      minimax_models[[i]]$fit(data[-(split_idx[[1]]),])
      # train_model$predict(data[flds[[k]],])
      scores[i] <- minimax_models[[i]]$score(data[split_idx[[1]],])
    }
  }
  
  best_idx <- which.max(scores)
  
  best_model <- minimax_models[[best_idx]]
  best_params <- list(lambda_pi = tune_grid$lambda_pi[best_idx], lambda_h = tune_grid$lambda_h[best_idx], gamma = tune_grid$gamma[best_idx])
  best_score <- scores[best_idx]
  
  return(list(best_model_ = best_model, best_params_ = best_params, best_score_ = best_score))
}
