# Load necessary library
library(MASS)

source("twostage_dgp_fun.R")

set.seed(12345)
n <- 1e8
X1 <- rnorm(n)
X2 <- rnorm(n)
X3 <- rbinom(n,1,0.5)

# Generate treatment assignment
p.A <- p.a.fun(X1,X2,X3)
A <- rbinom(n,1,p.A)

# Generate mediator
M <- rM_bimodal(X1, X2, X3, 0, sigma = sigma.m, delta = 4, p_left = 0.5)

# Generate outcomes
mu.y <- mu.y.fun(X1,X2,X3,1,M)
Y <- rnorm(n,mu.y,sigma.y)

psi.true <- mean(Y)


twostage_true <- list(psi.true = psi.true)
saveRDS(twostage_true, file = "twostage_true.rds")
