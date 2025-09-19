source("twostage_dataset.R")
source("twostage_minimax.R")

library(dplyr)
library(caret)
library(kernlab)
library(pracma)
library(evmix)
library(SuperLearner)
library(nnet)
library(randomForest)
library(gam)
library(ranger)

TwostageInference_real_dat <- setRefClass(
  "TwostageInference",
  fields = list(
    data = "TwostageData",
    crossfit_folds = "numeric",
    lambdaspi1 = "numeric",
    lambdash1 = "numeric",
    gammas1 = "numeric",
    lambdaspi2 = "numeric",
    lambdash2 = "numeric",
    gammas2 = "numeric",
    rf.min.node.size = "numeric",
    rf.num.trees = "numeric",
    split_ratio = "numeric",
    print_best_params = "logical",
    cf_inds = "list",
    pi1 = "list",
    pi2 = "list",
    pi1dim = "numeric",
    h1dim = "numeric",
    pi2dim = "numeric",
    h2dim = "numeric",
    mu1 = "list",
    mu2 = "list"
  ),
  methods = list(
    initialize = function(
    twostage_dataset,
    crossfit_folds = 1,
    lambdaspi1 = NULL,
    lambdash1 = NULL,
    gammas1 = NULL,
    lambdaspi2 = NULL,
    lambdash2 = NULL,
    gammas2 = NULL,
    rf.min.node.size = NULL,
    rf.num.trees = NULL,
    split_ratio = 0.25,
    print_best_params = FALSE
    ) {
      .self$data <- twostage_dataset
      .self$crossfit_folds <- crossfit_folds
      .self$lambdaspi1 <- lambdaspi1
      .self$lambdash1 <- lambdash1
      .self$gammas1 <- gammas1
      .self$lambdaspi2 <- lambdaspi2
      .self$lambdash2 <- lambdash2
      .self$gammas2 <- gammas2
      .self$rf.min.node.size <- rf.min.node.size
      .self$rf.num.trees <- rf.num.trees
      .self$split_ratio <- split_ratio
      .self$print_best_params <- print_best_params
      
      .self$pi1dim <- ncol(as.matrix(.self$data$r_pi1()))
      .self$h1dim <- ncol(as.matrix(.self$data$r_h1()))
      .self$pi2dim <- ncol(as.matrix(.self$data$r_pi2()))
      .self$h2dim <- ncol(as.matrix(.self$data$r_h2()))
      
      .self$cf_inds <- .self$data$create_crossfit_split(crossfit_folds)
      
      .self$pi1 <- lapply(1:crossfit_folds, function(i) {
        .self$estimate_pi1(fold = i)
      })
      
      .self$pi2 <- lapply(1:crossfit_folds, function(i) {
        .self$estimate_pi2(fold = i)
      })
      
      .self$mu1 <- lapply(1:crossfit_folds, function(i) {
        .self$estimate_mu1(fold = i)
      })
      
      .self$mu2 <- lapply(1:crossfit_folds, function(i) {
        .self$estimate_mu2(fold = i)
      })
      
    },
    
    estimate_pi1 = function(fold = 0) {
      g1 <- as.matrix(1-.self$data$A)
      g2 <- as.matrix(rep(1,length(.self$data$A)))
      dat <- cbind(.self$data$r_pi1(), .self$data$r_h1(), g1, g2)[.self$cf_inds[[fold]]$train, ]
      
      search <- MinimaxRKHSCV(
        dat,
        .self$pi1dim,
        .self$h1dim,
        lambdaspi = .self$lambdaspi1,
        lambdash = .self$lambdash1,
        gammas = .self$gammas1,
        split_ratio = .self$split_ratio
      )
      
      if (.self$print_best_params) {
        cat("pi1: gamma_1 = ", search$best_params_$gamma, ", lambda_h1 = ", search$best_params_$lambda_h, ", lambda_pi1 = ", search$best_params_$lambda_pi, ", best score: ", search$best_score_, "\n")
      }
      
      return(list(pi1_mod = search$best_model_$pi_,
                  h1_mod = search$best_model_$h_,
                  lambda_pi1 = search$best_params_$lambda_pi,
                  lambda_h1 = search$best_params_$lambda_h,
                  gamma_1 = search$best_params_$gamma,
                  score = search$best_score_))
    },
    
    estimate_pi2 = function(fold = 0) {
      
      r_pi1 <- .self$data$r_pi1()[.self$cf_inds[[fold]]$train, ]
      pi1.hat <- as.vector(.self$pi1[[fold]]$pi1_mod(r_pi1))
      # pi1.hat <- 1/(1-p.a.fun(.self$data$X[.self$cf_inds[[fold]]$train,1],.self$data$X[.self$cf_inds[[fold]]$train,2],.self$data$X[.self$cf_inds[[fold]]$train,3]))
      
      a <- as.matrix(.self$data$A[.self$cf_inds[[fold]]$train])
      g1 <- a
      g2 <- (1-a)*pi1.hat
      r_pi2 <- .self$data$r_pi2()[.self$cf_inds[[fold]]$train, ]
      r_h2 <- .self$data$r_h2()[.self$cf_inds[[fold]]$train, ]
      
      dat <- cbind(r_pi2, r_h2, g1, g2)
      
      search <- MinimaxRKHSCV(
        dat,
        .self$pi2dim,
        .self$h2dim,
        lambdaspi = .self$lambdaspi2,
        lambdash = .self$lambdash2,
        gammas = .self$gammas2,
        split_ratio = .self$split_ratio
      )
      
      if (.self$print_best_params) {
        cat("pi2: gamma_2 = ", search$best_params_$gamma, ", lambda_h2 = ", search$best_params_$lambda_h, ", lambda_pi2 = ", search$best_params_$lambda_pi, ", best score: ", search$best_score_, "\n")
      }
      
      return(list(pi2_mod = search$best_model_$pi_,
                  h2_mod = search$best_model_$h_,
                  lambda_pi2 = search$best_params_$lambda_pi,
                  lambda_h2 = search$best_params_$lambda_h,
                  gamma_2 = search$best_params_$gamma,
                  score = search$best_score_))
    },
    
    estimate_mu1 = function(fold = 0) {
      fold_train_idx <- .self$cf_inds[[fold]]$train
      x01 <- .self$data$X[fold_train_idx,1]
      x02 <- .self$data$X[fold_train_idx,2]
      x03 <- .self$data$X[fold_train_idx,3]
      x04 <- .self$data$X[fold_train_idx,4]
      x05 <- .self$data$X[fold_train_idx,5]
      x06 <- .self$data$X[fold_train_idx,6]
      x07 <- .self$data$X[fold_train_idx,7]
      x08 <- .self$data$X[fold_train_idx,8]
      x09 <- .self$data$X[fold_train_idx,9]
      x10 <- .self$data$X[fold_train_idx,10]
      x11 <- .self$data$X[fold_train_idx,11]
      m <- .self$data$M[fold_train_idx]
      a <- .self$data$A[fold_train_idx]
      y <- factor(.self$data$Y[fold_train_idx],levels = c(0,1))
      
      train_df1 <- data.frame(
        y1  = y[a == 1],
        m1 = m[a == 1],
        x011 = x01[a == 1],
        x021 = x02[a == 1],
        x031 = x03[a == 1],
        x041 = x04[a == 1],
        x051 = x05[a == 1],
        x061 = x06[a == 1],
        x071 = x07[a == 1],
        x081 = x08[a == 1],
        x091 = x09[a == 1],
        x101 = x10[a == 1],
        x111 = x11[a == 1]
      )
      
      mu1.fit <- ranger(
        y1 ~ .,
        data = train_df1,
        num.trees = .self$rf.num.trees,
        mtry = 3,
        min.node.size = .self$rf.min.node.size,
        probability = TRUE,
        splitrule = "gini",
        write.forest = TRUE
      )
      
      return(mu1.fit)
    },
    
    estimate_mu2 = function(fold = 0) {
      fold_train_idx <- .self$cf_inds[[fold]]$train
      x01 <- .self$data$X[fold_train_idx,1]
      x02 <- .self$data$X[fold_train_idx,2]
      x03 <- .self$data$X[fold_train_idx,3]
      x04 <- .self$data$X[fold_train_idx,4]
      x05 <- .self$data$X[fold_train_idx,5]
      x06 <- .self$data$X[fold_train_idx,6]
      x07 <- .self$data$X[fold_train_idx,7]
      x08 <- .self$data$X[fold_train_idx,8]
      x09 <- .self$data$X[fold_train_idx,9]
      x10 <- .self$data$X[fold_train_idx,10]
      x11 <- .self$data$X[fold_train_idx,11]
      m <- .self$data$M[fold_train_idx]
      a <- .self$data$A[fold_train_idx]
      # y <- .self$data$Y[fold_train_idx]
      
      x010 <- x01[a==0]
      x020 <- x02[a==0]
      x030 <- x03[a==0]
      x040 <- x04[a==0]
      x050 <- x05[a==0]
      x060 <- x06[a==0]
      x070 <- x07[a==0]
      x080 <- x08[a==0]
      x090 <- x09[a==0]
      x100 <- x10[a==0]
      x110 <- x11[a==0]
      m0 <- m[a==0]
      # y0 <- y[a==0]
      
      mu1.train <- predict(.self$mu1[[fold]],data = data.frame( m1  = m0,
                                                                x011 = x010,
                                                                x021 = x020,
                                                                x031 = x030,
                                                                x041 = x040,
                                                                x051 = x050,
                                                                x061 = x060,
                                                                x071 = x070,
                                                                x081 = x080,
                                                                x091 = x090,
                                                                x101 = x100,
                                                                x111 = x110))$predictions[, "1"]
      
      
      train_df0 <- data.frame(
        y2  = mu1.train,
        x010 = x010,
        x020 = x020,
        x030 = x030,
        x040 = x040,
        x050 = x050,
        x060 = x060,
        x070 = x070,
        x080 = x080,
        x090 = x090,
        x100 = x100,
        x110 = x110
      )
      
      mu2.fit <- ranger(
        y2 ~ .,
        data = train_df0,
        num.trees = .self$rf.num.trees,
        mtry = 3,
        min.node.size = .self$rf.min.node.size,
        splitrule = "variance",
        write.forest = TRUE
      )
      
      return(mu2.fit)
    },
    
    estimator = function() {
      estimates_naive <- c()
      estimates <- c()
      estimates.if <- list()
      for (fold in 1:(.self$crossfit_folds)) {
        fold_idx <- .self$cf_inds[[fold]]$eval
        
        r_pi1 <- .self$data$r_pi1()[fold_idx, ]
        pi1.hat <- .self$pi1[[fold]]$pi1_mod(r_pi1)
        
        r_pi2 <- .self$data$r_pi2()[fold_idx, ]
        pi2.hat <- .self$pi2[[fold]]$pi2_mod(r_pi2)
        
        mu1.hat <- predict(.self$mu1[[fold]], 
                           data = data.frame(m1 = .self$data$M[fold_idx],
                                             x011 = .self$data$X[fold_idx,1],
                                             x021 = .self$data$X[fold_idx,2],
                                             x031 = .self$data$X[fold_idx,3],
                                             x041 = .self$data$X[fold_idx,4],
                                             x051 = .self$data$X[fold_idx,5],
                                             x061 = .self$data$X[fold_idx,6],
                                             x071 = .self$data$X[fold_idx,7],
                                             x081 = .self$data$X[fold_idx,8],
                                             x091 = .self$data$X[fold_idx,9],
                                             x101 = .self$data$X[fold_idx,10],
                                             x111 = .self$data$X[fold_idx,11]))$predictions[, "1"]
        
        mu2.hat <- predict(.self$mu2[[fold]],
                           data = data.frame(x010 = .self$data$X[fold_idx,1],
                                             x020 = .self$data$X[fold_idx,2],
                                             x030 = .self$data$X[fold_idx,3],
                                             x040 = .self$data$X[fold_idx,4],
                                             x050 = .self$data$X[fold_idx,5],
                                             x060 = .self$data$X[fold_idx,6],
                                             x070 = .self$data$X[fold_idx,7],
                                             x080 = .self$data$X[fold_idx,8],
                                             x090 = .self$data$X[fold_idx,9],
                                             x100 = .self$data$X[fold_idx,10],
                                             x110 = .self$data$X[fold_idx,11]))$predictions
        
        ests <- as.numeric(mu2.hat) + (1-.self$data$A[fold_idx])*pi1.hat*(as.numeric(mu1.hat)-as.numeric(mu2.hat))+(.self$data$A[fold_idx]*pi2.hat)*(.self$data$Y[fold_idx]-as.numeric(mu1.hat))
        est <- mean(ests)
        estimates <- c(estimates, est)
        estimates.if[[fold]] <- ests
        
        ests_naive <- as.numeric(mu2.hat)
        est_naive <- mean(ests_naive)
        estimates_naive <- c(estimates_naive,est_naive)
      }
      
      estimates = mean(estimates)
      estimates_naive = mean(estimates_naive)
      
      estimates.sd <- c()
      for(fold in 1:(.self$crossfit_folds)) {
        estimates.sd <- c(estimates.sd, sd(estimates.if[[fold]]-estimates))
      }
      
      return(list(estimates = estimates,
                  estimates.sd = sqrt(mean(estimates.sd^2)),
                  estimates_naive = estimates_naive
                  
      ))
    }
    
  )
)