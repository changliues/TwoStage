library(KernSmooth)
library(kernlab)
library(ranger)
source("minimax_core.R")
source("twostage_datagen.R")

kernel_cond_density_function <- function(x, y, bandwidth_x, bandwidth_y = 1) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  
  # Inverse of bandwidth matrix
  if (is.vector(bandwidth_x)) {
    if (length(bandwidth_x) != p) stop("bandwidth_x vector must match number of columns in x")
    H_inv <- diag(1 / (bandwidth_x^2), p)
  } else if (is.matrix(bandwidth_x)) {
    if (!all(dim(bandwidth_x) == c(p, p))) stop("bandwidth_x matrix must be square")
    H_inv <- solve(bandwidth_x)
  } else {
    stop("bandwidth_x must be vector or square matrix")
  }
  
  # Gaussian kernel for y
  gaussian_kernel <- function(u) {
    exp(-0.5 * (u / bandwidth_y)^2) / (sqrt(2 * pi) * bandwidth_y)
  }
  
  # Returned estimator
  function(y0, x0) {
    x0 <- as.numeric(x0)
    if (length(x0) != p) stop("x0 must have length ", p)
    
    # Compute Mahalanobis distances
    diff <- sweep(x, 2, x0)         # n x p
    qform <- rowSums((diff %*% H_inv) * diff)
    w_x <- exp(-0.5 * qform)
    
    # Compute weights over y0
    sapply(y0, function(yy) {
      w_y <- gaussian_kernel(y - yy)
      num <- sum(w_x * w_y)
      denom <- sum(w_x)
      num / denom
    })
  }
}

run_TTS <- function(dataset,cf_inds, alpha = 0.05, rf.min.node.size, rf.num.trees){
  
  n <- dataset$n
  
  y <- dataset$Y
  m <- dataset$M
  x1 <- dataset$X[,1]
  x2 <- dataset$X[,2]
  x3 <- dataset$X[,3]
  x1s <- dataset$Xs[,1]
  x2s <- dataset$Xs[,2]
  x3s <- dataset$Xs[,3]
  a <- dataset$A
  
  # Create dataframe
  df <- data.frame(y = y, m = m, x1 = x1, x2 = x2, x3 = x3, a = a,x1s = x1s, x2s = x2s, x3s = x3s)
  colnames(df) <- c("y", "m", "x1", "x2", "x3", "a", "x1s", "x2s", "x3s")
  
  # Estimation
  crossfit_folds <- length(cf_inds)
  # cf_inds <- dataset$create_crossfit_split(4)
  a.x.hat <- rep(NA,n)
  y.1mx.hat <- rep(NA,n)
  y.0mx.hat <- rep(NA,n)
  m.1x.hat <- rep(NA,n)
  m.0x.hat <- rep(NA,n)
  eta.10x.hat <- rep(NA,n)
  for(fold in 1:crossfit_folds){
    
    fold_train_idx <- cf_inds[[fold]]$train
    fold_idx <- cf_inds[[fold]]$eval
    
    df.train <- df[fold_train_idx,]
    df.train1 <- df.train[df.train$a==1,]
    df.train0 <- df.train[df.train$a==0,]
    df.train$a <- factor(df.train$a, levels = c(0, 1))
    
    ## P(A=1|X)
    a.x.fit <- ranger(
      a ~ x1 + x2 + x3,
      data = df.train,
      num.trees = rf.num.trees,
      mtry = 2,
      min.node.size = rf.min.node.size,
      probability = TRUE,
      splitrule = "gini",
      write.forest = TRUE
    )
    
    # probs on eval fold: P(A=1|X)
    a.x.hat[fold_idx] <- predict(a.x.fit, data = df[fold_idx, c("x1","x2","x3")])$predictions[, "1"]
    
    ## E(Y|A=1,M,X)
    y.1mx.fit <- ranger(
      y ~ x1 + x2 + x3 + m,
      data = df.train1,
      num.trees = rf.num.trees,
      mtry = 2,
      min.node.size = rf.min.node.size,
      splitrule = "variance",
      write.forest = TRUE
    )
    
    y.1mx.hat[fold_idx] <- as.numeric(predict(y.1mx.fit,data = df[fold_idx, c("x1", "x2", "x3", "m")])$predictions)
    
    ## E(Y|A=0,M,X)
    y.0mx.fit <- ranger(
      y ~ x1 + x2 + x3 + m,
      data = df.train0,
      num.trees = rf.num.trees,
      mtry = 2,
      min.node.size = rf.min.node.size,
      splitrule = "variance",
      write.forest = TRUE
    )
    
    y.0mx.hat[fold_idx] <- as.numeric(predict(y.0mx.fit,data = df[fold_idx, c("x1", "x2", "x3", "m")])$predictions)
    
    ## p(M|A=a,X)
    x_test <- cbind(x1,x2,x3)[fold_idx, , drop = FALSE]
    m_test <- m[fold_idx]
    a1_train <- cbind(x1,x2,x3,m)[fold_train_idx & (a[fold_train_idx]==1), , drop = FALSE]
    a0_train <- cbind(x1,x2,x3,m)[fold_train_idx & (a[fold_train_idx]==0), , drop = FALSE]
    d <- 3
    
    # Helper: safe covariance with ridge
    safe_cov <- function(X, eps = 1e-6) {
      S <- stats::cov(X)
      S + eps * diag(ncol(X))
    }
    
    # Standardize X within group (improves covariance bandwidths)
    zscore_cols <- function(M, cols = 1:3) {
      mu  <- colMeans(M[, cols, drop = FALSE])
      sdv <- apply(M[, cols, drop = FALSE], 2, sd)
      sdv[sdv == 0] <- 1   # avoid divide-by-zero
      M2 <- M
      M2[, cols] <- sweep(sweep(M[, cols, drop = FALSE], 2, mu, "-"), 2, sdv, "/")
      list(M = M2, center = as.numeric(mu), scale = as.numeric(sdv))
    }
    
    # Standardize X columns within each arm
    sx1 <- zscore_cols(a1_train, 1:3)
    sx0 <- zscore_cols(a0_train, 1:3)
    a1z <- sx1$M; a0z <- sx0$M
    
    # Scott's rule with anisotropic bandwidth (full matrix)
    scott_exp <- -2 / (d + 4)
    bandwidth_x_opt1 <- safe_cov(a1z[, 1:3, drop = FALSE]) * nrow(a1z)^scott_exp
    bandwidth_x_opt0 <- safe_cov(a0z[, 1:3, drop = FALSE]) * nrow(a0z)^scott_exp
    
    bandwidth_y_opt1 <- KernSmooth::dpik(a1z[, 4])
    bandwidth_y_opt0 <- KernSmooth::dpik(a0z[, 4])
    
    # --- fit KDEs on standardized X ---
    m.1x.fit <- kernel_cond_density_function(
      a1z[, 1:3, drop = FALSE], a1z[, 4],
      bandwidth_x = bandwidth_x_opt1,
      bandwidth_y = bandwidth_y_opt1
    )
    m.0x.fit <- kernel_cond_density_function(
      a0z[, 1:3, drop = FALSE], a0z[, 4],
      bandwidth_x = bandwidth_x_opt0,
      bandwidth_y = bandwidth_y_opt0
    )
    
    # --- build m grid + weights ---
    m_range <- seq(-8, 25, by = 0.1)
    dm <- c(diff(m_range), tail(diff(m_range), 1))
    
    # --- RF predictions on full (x, m) grid (raw scale; RF was fit on raw X) ---
    n_val <- length(fold_idx)
    n_m   <- length(m_range)
    x_rep <- x_test[seq_len(n_val), , drop = FALSE]
    newdata_mat <- data.frame(
      x1 = rep(x_rep[,1], each = n_m),
      x2 = rep(x_rep[,2], each = n_m),
      x3 = rep(x_rep[,3], each = n_m),
      m  = rep(m_range, times = n_val)
    )
    yhat_vec <- as.numeric(predict(y.1mx.fit, data = newdata_mat)$predictions)
    Yhat_mat <- matrix(yhat_vec, nrow = n_val, ncol = n_m, byrow = TRUE)
    
    # --- STANDARDIZE x_test for KDE prediction (use TRAINING-arm stats!) ---
    x_test_0z <- cbind(
      (x_test[,1] - sx0$center[1]) / sx0$scale[1],
      (x_test[,2] - sx0$center[2]) / sx0$scale[2],
      (x_test[,3] - sx0$center[3]) / sx0$scale[3]
    )
    x_test_1z <- cbind(
      (x_test[,1] - sx1$center[1]) / sx1$scale[1],
      (x_test[,2] - sx1$center[2]) / sx1$scale[2],
      (x_test[,3] - sx1$center[3]) / sx1$scale[3]
    )
    
    # --- density matrices over m_range ---
    D0 <- matrix(NA, nrow = n_val, ncol = n_m)
    for (i in seq_len(n_val)) D0[i, ] <- m.0x.fit(m_range, x_test_0z[i, ])
    
    # normalize rows so they integrate to 1 under dm
    row_sums <- as.numeric(D0 %*% dm)
    ok <- row_sums > .Machine$double.eps
    D0[ok, ] <- sweep(D0[ok, , drop = FALSE], 1, row_sums[ok], "/")
    
    # --- eta_10(x) = ∫ E[Y|M=m,X,A=1] p(m|X,A=0) dm ---
    eta_vals <- rowSums(D0 * Yhat_mat * rep(dm, each = n_val))
    
    # --- point densities at observed m (need standardized x for each arm) ---
    m_obs <- m[fold_idx]
    m0_point <- vapply(seq_len(n_val), function(i) m.0x.fit(m_obs[i], x_test_0z[i, ]), numeric(1))
    m1_point <- vapply(seq_len(n_val), function(i) m.1x.fit(m_obs[i], x_test_1z[i, ]), numeric(1))
    
    # write back
    eta.10x.hat[fold_idx] <- eta_vals
    m.0x.hat[fold_idx]    <- m0_point
    m.1x.hat[fold_idx]    <- m1_point
    
  }
  
  ## influence-function of psi
  if.psi0 <- ((a==1)*m.0x.hat)/(a.x.hat*m.1x.hat)*(y-y.1mx.hat)+((a==0)/(1-a.x.hat))*(y.1mx.hat-eta.10x.hat)+eta.10x.hat
  ## influence-function-based estimator of psi
  if.est.psi0  <- mean(if.psi0, na.rm = T)
  
  # 95% CI
  ### get standard derivation of each fold
  var0.vec <- rep(NA,crossfit_folds)
  
  for(fold in 1:crossfit_folds){
    fold_idx <- cf_inds[[fold]]$eval
    var0.vec[fold] <- var(if.psi0[fold_idx]-if.est.psi0, na.rm = T)
  }
  if.sd.psi0 <- sqrt(mean(var0.vec, na.rm = T))/sqrt(n)
  critical_t <- qnorm(1 - alpha/2)
  lower0 <- if.est.psi0  - critical_t * if.sd.psi0
  upper0 <- if.est.psi0  + critical_t * if.sd.psi0
  
  result <- list(if.est.psi0 = if.est.psi0,
                 if.sd.psi0 = if.sd.psi0,
                 if.lower.psi0 = lower0,
                 if.upper.psi0 = upper0)
  
  return(result)
}