expit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}

# p.a.fun <- function(x1,x2,x3){
#   expit(-1 + 0.1*x1)
# }

p.a.fun <- function(x1,x2,x3){
  expit(0.35 + 0.1*x1 - 0.2*x2 + 0.3*x3 - 0.15*x1^2 + 0.25*x2*x3 )
}

mu.m.fun <- function(x1,x2,x3,a){
  2.4 + 1.6*x1 - 1.2*x2 + 2.2*x3 + 0.8*a + 1.5*x1^2 - 1.8*x2*x3
}

sigma.m <- 1

## --- bimodal generator ---
rM_bimodal <- function(x1, x2, x3, a,
                       sigma = sigma.m,   # within‑mode SD  (keep = 1 to match original)
                       delta = 4,   # half‑distance between the two modes
                       p_left = 0.5 # probability of drawing from the “left” mode
){
  mu <- mu.m.fun(x1, x2, x3, a)
  
  # Draw ±delta with probability p_left / (1 - p_left)
  shift <- ifelse(runif(length(mu)) < p_left, -delta, +delta)
  
  rnorm(length(mu), mean = mu + shift, sd = sigma)
}

dM_bimodal <- function(m, x1, x2, x3, a,
                       sigma = sigma.m,
                       delta = 4,
                       p_left = 0.5) {
  mu <- mu.m.fun(x1, x2, x3, a)
  d_left <- dnorm(m, mean = mu - delta, sd = sigma)
  d_right <- dnorm(m, mean = mu + delta, sd = sigma)
  p_left * d_left + (1 - p_left) * d_right
}


mu.y.fun <- function(x1,x2,x3,a,m){
  6.4 + 2.7*x1 + 0.7*x2 - 3.6*x3 - 0.8*x2^2 + 2.5*x1*a + 3.2*m*a - 2.9*m + 4.5*a
}

sigma.y <- 1