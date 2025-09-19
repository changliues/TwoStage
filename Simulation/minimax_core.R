library(pracma)

minimax_solve <- function(gram1, gram2, g1, g2, lambda1, lambda2){
  n <- nrow(gram1)
  
  gamma_matrix <- (1/4) * gram2 %*% chol2inv(chol((1/(4*n)) * gram2 + lambda2 * diag(1, nrow = n)))
  minimax_1 <- gram1 %*% (diag(as.vector(g1))) %*% gamma_matrix %*% (diag(as.vector(g1))) %*% gram1 + n^2*lambda1*gram1
  
  alpha <- as.vector(ginv(minimax_1) %*% gram1 %*% (diag(as.vector(g1))) %*% gamma_matrix %*% g2)
  beta <- as.vector(1/2*chol2inv(chol((1/(4*n)) * gram2 + lambda2 * diag(1, nrow = n))) %*% t(1/n*(alpha %*% gram1 * g1 - g2)))
  
  return(list(alpha=alpha, beta=beta))
}


score_nuisance_function <- function(pi_values, gram2_score, g1, g2, lambda_score) {
  n <- nrow(gram2_score)
  
  beta <- as.vector(1/2*chol2inv(chol((1/(4*n)) * gram2_score + lambda_score * diag(1, nrow = n))) %*% t(1/n*(pi_values * g1 - g2)))
  h_values <- beta %*% gram2_score
  
  metric <- mean((g1 * pi_values - g2) * h_values - (1/4) *h_values^2)
  return(-metric)
}
