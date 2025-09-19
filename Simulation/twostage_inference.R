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

TwostageInference <- setRefClass(
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
    lambdaspi1_non_optimal = "numeric",
    lambdash1_non_optimal = "numeric",
    gammas1_non_optimal = "numeric",
    lambdaspi2_non_optimal = "numeric",
    lambdash2_non_optimal = "numeric",
    gammas2_non_optimal = "numeric",
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
    mu2 = "list",
    pi1_misspec = "list",
    pi2_misspec = "list",
    pi2_pi1misspec = "list",
    mu1_misspec = "list",
    mu2_misspec = "list",
    mu2_mu1misspec = "list"
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
    lambdaspi1_non_optimal = NULL,
    lambdash1_non_optimal = NULL,
    gammas1_non_optimal = NULL,
    lambdaspi2_non_optimal = NULL,
    lambdash2_non_optimal = NULL,
    gammas2_non_optimal = NULL,
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
      .self$lambdaspi1_non_optimal <- lambdaspi1_non_optimal
      .self$lambdash1_non_optimal <- lambdash1_non_optimal
      .self$gammas1_non_optimal <- gammas1_non_optimal
      .self$lambdaspi2_non_optimal <- lambdaspi2_non_optimal
      .self$lambdash2_non_optimal <- lambdash2_non_optimal
      .self$gammas2_non_optimal <- gammas2_non_optimal
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
      
      .self$pi1_misspec <- lapply(1:crossfit_folds, function(i) {
        .self$estimate_pi1_misspec(fold = i)
      })
      
      .self$pi2_misspec <- lapply(1:crossfit_folds, function(i) {
        .self$estimate_pi2_misspec(fold = i)
      })
      
      .self$pi2_pi1misspec <- lapply(1:crossfit_folds, function(i) {
        .self$estimate_pi2_pi1misspec(fold = i)
      })
      
      .self$mu1_misspec <- lapply(1:crossfit_folds, function(i) {
        .self$estimate_mu1_misspec(fold = i)
      })
      
      .self$mu2_misspec <- lapply(1:crossfit_folds, function(i) {
        .self$estimate_mu2_misspec(fold = i)
      })
      
      .self$mu2_mu1misspec <- lapply(1:crossfit_folds, function(i) {
        .self$estimate_mu2_mu1misspec(fold = i)
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
      x1 <- .self$data$X[fold_train_idx,1]
      x2 <- .self$data$X[fold_train_idx,2]
      x3 <- .self$data$X[fold_train_idx,3]
      m <- .self$data$M[fold_train_idx]
      a <- .self$data$A[fold_train_idx]
      y <- .self$data$Y[fold_train_idx]

      train_df1 <- data.frame(
        y1  = y[a == 1],
        m1 = m[a == 1],
        x11 = x1[a == 1],
        x21 = x2[a == 1],
        x31 = x3[a == 1]
      )
      
      mu1.fit <- ranger(
        y1 ~ .,
        data = train_df1,
        num.trees = .self$rf.num.trees,
        mtry = 2,
        min.node.size = .self$rf.min.node.size,
        splitrule = "variance",
        write.forest = TRUE
      )
      
      return(mu1.fit)
    },

    estimate_mu2 = function(fold = 0) {
      fold_train_idx <- .self$cf_inds[[fold]]$train
      x1 <- .self$data$X[fold_train_idx,1]
      x2 <- .self$data$X[fold_train_idx,2]
      x3 <- .self$data$X[fold_train_idx,3]
      m <- .self$data$M[fold_train_idx]
      a <- .self$data$A[fold_train_idx]
      y <- .self$data$Y[fold_train_idx]

      x10 <- x1[a==0]
      x20 <- x2[a==0]
      x30 <- x3[a==0]
      m0 <- m[a==0]
      y0 <- y[a==0]

      mu1.train <- predict(.self$mu1[[fold]],data = data.frame( m1  = m0,
                                                                x11 = x10,
                                                                x21 = x20,
                                                                x31 = x30))$predictions
         
      
      train_df0 <- data.frame(
        y2  = mu1.train,
        x10 = x10,
        x20 = x20,
        x30 = x30
      )
      
      mu2.fit <- ranger(
        y2 ~ .,
        data = train_df0,
        num.trees = .self$rf.num.trees,
        mtry = 2,
        min.node.size = .self$rf.min.node.size,
        splitrule = "variance",
        write.forest = TRUE
      )
      
      return(mu2.fit)
    },
    
    estimate_pi1_misspec = function(fold = 0) {
      g1 <- as.matrix(1-.self$data$A)
      g2 <- as.matrix(rep(1,length(.self$data$A)))
      dat <- cbind(.self$data$r_pi1(), .self$data$r_h1(), g1, g2)[.self$cf_inds[[fold]]$train, ]
      
      search <- MinimaxRKHSCV(
        dat,
        .self$pi1dim,
        .self$h1dim,
        lambdaspi = .self$lambdaspi1_non_optimal,
        lambdash = .self$lambdash1_non_optimal,
        gammas = .self$gammas1_non_optimal,
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
    
    estimate_pi2_misspec = function(fold = 0) {
      
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
        lambdaspi = .self$lambdaspi2_non_optimal,
        lambdash = .self$lambdash2_non_optimal,
        gammas = .self$gammas2_non_optimal,
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
    
    estimate_pi2_pi1misspec = function(fold = 0) {
      
      r_pi1 <- .self$data$r_pi1()[.self$cf_inds[[fold]]$train, ]
      pi1.hat <- as.vector(.self$pi1_misspec[[fold]]$pi1_mod(r_pi1))
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
        lambdaspi = .self$lambdaspi2_non_optimal,
        lambdash = .self$lambdash2_non_optimal,
        gammas = .self$gammas2_non_optimal,
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
    
    estimate_mu1_misspec = function(fold = 0) {
      fold_train_idx <- .self$cf_inds[[fold]]$train
      x1 <- .self$data$X[fold_train_idx,1]
      x2 <- .self$data$X[fold_train_idx,2]
      x3 <- .self$data$X[fold_train_idx,3]
      m <- .self$data$M[fold_train_idx]
      a <- .self$data$A[fold_train_idx]
      y <- .self$data$Y[fold_train_idx]
      
      train_df1 <- data.frame(
        y1  = y[a == 1],
        m1 = m[a == 1],
        x11 = x1[a == 1],
        x21 = x2[a == 1],
        x31 = x3[a == 1]
      )
      
      mu1.fit <- ranger(
        y1 ~ .,
        data = train_df1,
        num.trees = 50,
        mtry = 1,
        min.node.size = 2,
        splitrule = "variance",
        write.forest = TRUE
      )
      
      return(mu1.fit)
    },
    
    estimate_mu2_misspec = function(fold = 0) {
      fold_train_idx <- .self$cf_inds[[fold]]$train
      x1 <- .self$data$X[fold_train_idx,1]
      x2 <- .self$data$X[fold_train_idx,2]
      x3 <- .self$data$X[fold_train_idx,3]
      m <- .self$data$M[fold_train_idx]
      a <- .self$data$A[fold_train_idx]
      y <- .self$data$Y[fold_train_idx]
      
      x10 <- x1[a==0]
      x20 <- x2[a==0]
      x30 <- x3[a==0]
      m0 <- m[a==0]
      y0 <- y[a==0]
      
      mu1.train <- predict(.self$mu1[[fold]],data = data.frame(m1  = m0,
                                                               x11 = x10,
                                                               x21 = x20,
                                                               x31 = x30))$predictions
      
      train_df0 <- data.frame(
        y2  = mu1.train,
        x10 = x10,
        x20 = x20,
        x30 = x30
      )
      
      mu2.fit <- ranger(
        y2 ~ .,
        data = train_df0,
        num.trees = 50,
        mtry = 1,
        min.node.size = 2,
        splitrule = "variance",
        write.forest = TRUE
      )
      
      return(mu2.fit)
    },
    
    estimate_mu2_mu1misspec = function(fold = 0) {
      fold_train_idx <- .self$cf_inds[[fold]]$train
      x1 <- .self$data$X[fold_train_idx,1]
      x2 <- .self$data$X[fold_train_idx,2]
      x3 <- .self$data$X[fold_train_idx,3]
      m <- .self$data$M[fold_train_idx]
      a <- .self$data$A[fold_train_idx]
      y <- .self$data$Y[fold_train_idx]
      
      x10 <- x1[a==0]
      x20 <- x2[a==0]
      x30 <- x3[a==0]
      m0 <- m[a==0]
      y0 <- y[a==0]
      
      mu1.train <- predict(.self$mu1_misspec[[fold]],data = data.frame(m1  = m0,
                                                                       x11 = x10,
                                                                       x21 = x20,
                                                                       x31 = x30))$predictions
      
      train_df0 <- data.frame(
        y2  = mu1.train,
        x10 = x10,
        x20 = x20,
        x30 = x30
      )
      
      mu2.fit <- ranger(
        y2 ~ .,
        data = train_df0,
        num.trees = 50,
        mtry = 1,
        min.node.size = 2,
        splitrule = "variance",
        write.forest = TRUE
      )
      
      return(mu2.fit)
    },

    estimator = function() {
      estimates0naive <- c()
      estimates1naive <- c()
      estimates2naive <- c()
      estimates3naive <- c()
      estimates4naive <- c()
      estimates0 <- c()
      estimates0.if <- list()
      estimates1 <- c()
      estimates1.if <- list()
      estimates2 <- c()
      estimates2.if <- list()
      estimates3 <- c()
      estimates3.if <- list()
      estimates4 <- c()
      estimates4.if <- list()
      for (fold in 1:(.self$crossfit_folds)) {
        fold_idx <- .self$cf_inds[[fold]]$eval

        r_pi1 <- .self$data$r_pi1()[fold_idx, ]
        pi1.hat <- .self$pi1[[fold]]$pi1_mod(r_pi1)

        r_pi2 <- .self$data$r_pi2()[fold_idx, ]
        pi2.hat <- .self$pi2[[fold]]$pi2_mod(r_pi2)
        
        mu1.hat <- predict(.self$mu1[[fold]], 
                           data = data.frame(m1 = .self$data$M[fold_idx],
                                             x11 = .self$data$X[fold_idx,1],
                                             x21 = .self$data$X[fold_idx,2],
                                             x31 = .self$data$X[fold_idx,3]))$predictions
        
        mu2.hat <- predict(.self$mu2[[fold]],
                           data = data.frame(x10 = .self$data$X[fold_idx,1],
                                             x20 = .self$data$X[fold_idx,2],
                                             x30 = .self$data$X[fold_idx,3]))$predictions
                             

        pi1.hat.biased <- .self$pi1_misspec[[fold]]$pi1_mod(r_pi1)
        
        pi2.hat.biased <- .self$pi2_misspec[[fold]]$pi2_mod(r_pi2)
        
        pi2.hat.pi1biased <- .self$pi2_pi1misspec[[fold]]$pi2_mod(r_pi2)
        
        mu1.hat.biased <- predict(.self$mu1_misspec[[fold]],
                                  data = data.frame(x11 = .self$data$X[fold_idx,1],
                                                    x21 = .self$data$X[fold_idx,2],
                                                    x31 = .self$data$X[fold_idx,3],
                                                    m1 = .self$data$M[fold_idx]))$predictions
                                    
        
        mu2.hat.biased <- predict(.self$mu2_misspec[[fold]], 
                                  data = data.frame(x10 = .self$data$X[fold_idx,1],
                                                    x20 = .self$data$X[fold_idx,2],
                                                    x30 = .self$data$X[fold_idx,3]))$predictions
        
        mu2.hat.mu1biased <- predict(.self$mu2_mu1misspec[[fold]], 
                                  data = data.frame(x10 = .self$data$X[fold_idx,1],
                                                    x20 = .self$data$X[fold_idx,2],
                                                    x30 = .self$data$X[fold_idx,3]))$predictions
                                    
        
        # all true case
        ests0 <- as.numeric(mu2.hat) + (1-.self$data$A[fold_idx])*pi1.hat*(as.numeric(mu1.hat)-as.numeric(mu2.hat))+(.self$data$A[fold_idx]*pi2.hat)*(.self$data$Y[fold_idx]-as.numeric(mu1.hat))
        est0 <- mean(ests0)
        estimates0 <- c(estimates0, est0)
        estimates0.if[[fold]] <- ests0
        
        ests0naive <- as.numeric(mu2.hat)
        est0naive <- mean(ests0naive)
        estimates0naive <- c(estimates0naive,est0naive)
        
        # only pi1,pi2 true case
        ests1 <- as.numeric(mu2.hat.mu1biased) + (1-.self$data$A[fold_idx])*pi1.hat*(as.numeric(mu1.hat.biased)-as.numeric(mu2.hat.mu1biased))+(.self$data$A[fold_idx]*pi2.hat)*(.self$data$Y[fold_idx]-as.numeric(mu1.hat.biased))
        est1 <- mean(ests1)
        estimates1 <- c(estimates1, est1)
        estimates1.if[[fold]] <- ests1
        
        ests1naive <- as.numeric(mu2.hat.mu1biased)
        est1naive <- mean(ests1naive)
        estimates1naive <- c(estimates1naive,est1naive)
        
        # only mu1,mu2 true case
        ests2 <- as.numeric(mu2.hat) + (1-.self$data$A[fold_idx])*pi1.hat.biased*(as.numeric(mu1.hat)-as.numeric(mu2.hat))+(.self$data$A[fold_idx]*pi2.hat.pi1biased)*(.self$data$Y[fold_idx]-as.numeric(mu1.hat))
        est2 <- mean(ests2)
        estimates2 <- c(estimates2, est2)
        estimates2.if[[fold]] <- ests2
        
        ests2naive <- as.numeric(mu2.hat)
        est2naive <- mean(ests2naive)
        estimates2naive <- c(estimates2naive,est2naive)
        
        # only pi1,mu1 true case
        ests3 <- as.numeric(mu2.hat.biased) + (1-.self$data$A[fold_idx])*pi1.hat*(as.numeric(mu1.hat)-as.numeric(mu2.hat.biased))+(.self$data$A[fold_idx]*pi2.hat.biased)*(.self$data$Y[fold_idx]-as.numeric(mu1.hat))
        est3 <- mean(ests3)
        estimates3 <- c(estimates3, est3)
        estimates3.if[[fold]] <- ests3
        
        ests3naive <- as.numeric(mu2.hat.biased)
        est3naive <- mean(ests3naive)
        estimates3naive <- c(estimates3naive,est3naive)
        
        # all false case
        ests4 <- as.numeric(mu2.hat.mu1biased) + (1-.self$data$A[fold_idx])*pi1.hat.biased*(as.numeric(mu1.hat.biased)-as.numeric(mu2.hat.mu1biased))+(.self$data$A[fold_idx]*pi2.hat.pi1biased)*(.self$data$Y[fold_idx]-as.numeric(mu1.hat.biased))
        est4 <- mean(ests4)
        estimates4 <- c(estimates4, est4)
        estimates4.if[[fold]] <- ests4
        
        ests4naive <- as.numeric(mu2.hat.mu1biased)
        est4naive <- mean(ests4naive)
        estimates4naive <- c(estimates4naive,est4naive)
      }
      
      estimates0 = mean(estimates0)
      estimates1 = mean(estimates1)
      estimates2 = mean(estimates2)
      estimates3 = mean(estimates3)
      estimates4 = mean(estimates4)
      estimates0naive = mean(estimates0naive)
      estimates1naive = mean(estimates1naive)
      estimates2naive = mean(estimates2naive)
      estimates3naive = mean(estimates3naive)
      estimates4naive = mean(estimates4naive)
      
      estimates0.sd <- c()
      estimates1.sd <- c()
      estimates2.sd <- c()
      estimates3.sd <- c()
      estimates4.sd <- c()
      for(fold in 1:(.self$crossfit_folds)) {
        estimates0.sd <- c(estimates0.sd, sd(estimates0.if[[fold]]-estimates0))
        estimates1.sd <- c(estimates1.sd, sd(estimates1.if[[fold]]-estimates1))
        estimates2.sd <- c(estimates2.sd, sd(estimates2.if[[fold]]-estimates2))
        estimates3.sd <- c(estimates3.sd, sd(estimates3.if[[fold]]-estimates3))
        estimates4.sd <- c(estimates4.sd, sd(estimates4.if[[fold]]-estimates4))
      }
      
      return(list(estimates0 = estimates0,
                  estimates0.sd = sqrt(mean(estimates0.sd^2)),
                  estimates1 = estimates1,
                  estimates1.sd = sqrt(mean(estimates1.sd^2)),
                  estimates2 = estimates2,
                  estimates2.sd = sqrt(mean(estimates2.sd^2)),
                  estimates3 = estimates3,
                  estimates3.sd = sqrt(mean(estimates3.sd^2)),
                  estimates4 = estimates4,
                  estimates4.sd = sqrt(mean(estimates4.sd^2)),
                  estimates0naive = estimates0naive,
                  estimates1naive = estimates1naive,
                  estimates2naive = estimates2naive,
                  estimates3naive = estimates3naive,
                  estimates4naive = estimates4naive
                                      
      ))
    }
    
  )
)