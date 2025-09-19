# Load necessary library
library(MASS)

source("twostage_dataset.R")
source("twostage_dgp_fun.R")
twostage_true <- readRDS(file = "twostage_true.rds")

# Define the generate_data function
generate_data <- function(n) {
  
  # Generate covariates
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  X3 <- rbinom(n,1,0.5)
  X <- cbind(X1,X2,X3)
  colnames(X) <- c("X1","X2","X3")
  
  # generate fake covariates to test robustness
  X1s <- ifelse(X1 <= -0.2, -1, ifelse(X1>=0.3,1,0))
  X2s <- ifelse(X2 <= -0.3, 3, ifelse(X2>=0.2,1,2))
  X3s <- X3+rnorm(n,0,1)
  Xs <- cbind(X1s,X2s,X3s)
  colnames(Xs) <- c("X1s","X2s","X3s")
  
  # Generate treatment assignment
  p.A <- p.a.fun(X1,X2,X3)
  A <- rbinom(n,1,p.A)
  
  # Generate mediator
  M <- rM_bimodal(X1, X2, X3, A, sigma = sigma.m, delta = 4, p_left = 0.5)
  
  # Generate outcomes
  mu.y <- mu.y.fun(X1,X2,X3,A,M)
  Y <- rnorm(n,mu.y,sigma.y)
  
  dataset <- TwostageData$new(X = X, A = A, M = M, Y = Y,Xs = Xs)
  return(list(dataset = dataset, target = twostage_true))
}


