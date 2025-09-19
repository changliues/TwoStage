source("minimax_core.R")
source("twostage_dataset.R")
source("twostage_minimax.R")
source("twostage_inference.R")
source("twostage_datagen.R")
source("twostage_TTS.R")


# sample sizes
sample.size.for.sim = 1000
# # number of folds
split_ratio <- 0.25
n.flds <- 4
sample.size.for.tuning <- sample.size.for.sim/n.flds*(1-split_ratio)
# # CI parameters
alpha <- 0.05
critical_t <- qnorm(1 - alpha/2)
# # gammas and lambdas sequences
gammas1 <- c(0.08,0.1,0.12)
lambdash1 <- c(6,6.25,6.5)/sample.size.for.tuning^0.8
lambdaspi1 <- c(4600,4700,4800)/sample.size.for.tuning^3.6/lambdash1
gammas2 <- c(0.08,0.1,0.12)
lambdash2 <- c(6,6.25,6.5)/sample.size.for.tuning^0.8
lambdaspi2 <- c(4600,4700,4800)/sample.size.for.tuning^3.6/lambdash2
rf.min.node.size <- 10
rf.num.trees <- 800
# # non-optimal gammas and lambdas
gammas1_non_optimal <- 2
lambdash1_non_optimal <- 20/sample.size.for.tuning^0.8
lambdaspi1_non_optimal <- 10000/sample.size.for.tuning^3.6/lambdash1
gammas2_non_optimal <- 2
lambdash2_non_optimal <- 20/sample.size.for.tuning^0.8
lambdaspi2_non_optimal <- 10000/sample.size.for.tuning^3.6/lambdash2

params <- list()
for(j in 1:n.flds){
  params[[j]] <- list()
  params[[j]]$pi1_params <- as.data.frame(matrix(NA, nrow = 1, ncol = 4))
  names(params[[j]]$pi1_params) <- c("lambda_pi1", "lambda_h", "gamma", "score")
  params[[j]]$pi2_params <- as.data.frame(matrix(NA, nrow = 1, ncol = 4))
  names(params[[j]]$pi2_params) <- c("lambda_pi2", "lambda_h", "gamma", "score")
  params[[j]]$pi1_misspec_params <- as.data.frame(matrix(NA, nrow = 1, ncol = 4))
  names(params[[j]]$pi1_misspec_params) <- c("lambda_pi1", "lambda_h", "gamma", "score")
  params[[j]]$pi2_misspec_params <- as.data.frame(matrix(NA, nrow = 1, ncol = 4))
  names(params[[j]]$pi2_misspec_params) <- c("lambda_pi2", "lambda_h", "gamma", "score")
  params[[j]]$pi2_pi1misspec_params <- as.data.frame(matrix(NA, nrow = 1, ncol = 4))
  names(params[[j]]$pi2_pi1misspec_params) <- c("lambda_pi2", "lambda_h", "gamma", "score")
  
}

run_sim <- function(n.sim){
  
  data <- generate_data(sample.size.for.sim)
  dataset <- data$dataset
  target <- data$target

  estimator <- TwostageInference$new(
    dataset,
    crossfit_folds = n.flds,
    lambdash1 =lambdash1,
    lambdaspi1 = lambdaspi1,
    gammas1 = gammas1,
    lambdash2 =lambdash2,
    lambdaspi2 = lambdaspi2,
    gammas2 = gammas2,
    lambdash1_non_optimal =lambdash1_non_optimal,
    lambdaspi1_non_optimal = lambdaspi1_non_optimal,
    gammas1_non_optimal = gammas1_non_optimal,
    lambdash2_non_optimal =lambdash2_non_optimal,
    lambdaspi2_non_optimal = lambdaspi2_non_optimal,
    gammas2_non_optimal = gammas2_non_optimal,
    rf.min.node.size = rf.min.node.size,
    rf.num.trees = rf.num.trees,
    split_ratio = split_ratio,
    print_best_params = T)
  est <- estimator$estimator()

  critical_t <- qnorm(1 - alpha/2)
  est0.2stg <- est$estimates0
  est0.naive <- est$estimates0naive
  est0.2stg.sd <- est$estimates0.sd
  est0.2stg.lower <- est$estimates0-critical_t*est$estimates0.sd/sqrt(sample.size.for.sim)
  est0.2stg.upper <- est$estimates0+critical_t*est$estimates0.sd/sqrt(sample.size.for.sim)
  est1.2stg <- est$estimates1
  est1.naive <- est$estimates1naive
  est1.2stg.sd <- est$estimates1.sd
  est1.2stg.lower <- est$estimates1-critical_t*est$estimates1.sd/sqrt(sample.size.for.sim)
  est1.2stg.upper <- est$estimates1+critical_t*est$estimates1.sd/sqrt(sample.size.for.sim)
  est2.2stg <- est$estimates2
  est2.naive <- est$estimates2naive
  est2.2stg.sd <- est$estimates2.sd
  est2.2stg.lower <- est$estimates2-critical_t*est$estimates2.sd/sqrt(sample.size.for.sim)
  est2.2stg.upper <- est$estimates2+critical_t*est$estimates2.sd/sqrt(sample.size.for.sim)
  est3.2stg <- est$estimates3
  est3.naive <- est$estimates3naive
  est3.2stg.sd <- est$estimates3.sd
  est3.2stg.lower <- est$estimates3-critical_t*est$estimates3.sd/sqrt(sample.size.for.sim)
  est3.2stg.upper <- est$estimates3+critical_t*est$estimates3.sd/sqrt(sample.size.for.sim)
  est4.2stg <- est$estimates4
  est4.naive <- est$estimates4naive
  est4.2stg.sd <- est$estimates4.sd
  est4.2stg.lower <- est$estimates4-critical_t*est$estimates4.sd/sqrt(sample.size.for.sim)
  est4.2stg.upper <- est$estimates4+critical_t*est$estimates4.sd/sqrt(sample.size.for.sim)
  
  pi1 <- estimator$pi1
  pi2 <- estimator$pi2
  pi1_misspec <- estimator$pi1_misspec
  pi2_misspec <- estimator$pi2_misspec
  pi2_pi1misspec <- estimator$pi2_pi1misspec
  
  for(j in 1:n.flds){
    params[[j]]$pi1_params <- c(pi1[[j]]$lambda_pi1,pi1[[j]]$lambda_h1,pi1[[j]]$gamma_1, pi1[[j]]$score)
    params[[j]]$pi2_params <- c(pi2[[j]]$lambda_pi2,pi2[[j]]$lambda_h2,pi2[[j]]$gamma_2, pi2[[j]]$score)
    params[[j]]$pi1_misspec_params <- c(pi1_misspec[[j]]$lambda_pi1,pi1_misspec[[j]]$lambda_h1,pi1_misspec[[j]]$gamma_1, pi1_misspec[[j]]$score)
    params[[j]]$pi2_misspec_params <- c(pi2_misspec[[j]]$lambda_pi2,pi2_misspec[[j]]$lambda_h2,pi2_misspec[[j]]$gamma_2, pi2_misspec[[j]]$score)
    params[[j]]$pi2_pi1misspec_params <- c(pi2_pi1misspec[[j]]$lambda_pi2,pi2_pi1misspec[[j]]$lambda_h2,pi2_pi1misspec[[j]]$gamma_2, pi2_pi1misspec[[j]]$score)
  }
  
  result.2stg0 <- list(true.psi = target$psi.true,
                       est.2stg = est0.2stg,
                       est.naive = est0.naive,
                       est.2stg.sd = est0.2stg.sd,
                       est.2stg.lower = est0.2stg.lower,
                       est.2stg.upper = est0.2stg.upper,
                       ifb.cov = target$psi.true >= est0.2stg.lower & target$psi.true <= est0.2stg.upper,
                       sample.size = sample.size.for.sim
  )
  
  result.2stg1 <- list(true.psi = target$psi.true,
                       est.2stg = est1.2stg,
                       est.naive = est1.naive,
                       est.2stg.sd = est1.2stg.sd,
                       est.2stg.lower = est1.2stg.lower,
                       est.2stg.upper = est1.2stg.upper,
                       ifb.cov = target$psi.true >= est1.2stg.lower & target$psi.true <= est1.2stg.upper,
                       sample.size = sample.size.for.sim
  )
  
  result.2stg2 <- list(true.psi = target$psi.true,
                       est.2stg = est2.2stg,
                       est.naive = est2.naive,
                       est.2stg.sd = est2.2stg.sd,
                       est.2stg.lower = est2.2stg.lower,
                       est.2stg.upper = est2.2stg.upper,
                       ifb.cov = target$psi.true >= est2.2stg.lower & target$psi.true <= est2.2stg.upper,
                       sample.size = sample.size.for.sim
  )
  
  result.2stg3 <- list(true.psi = target$psi.true,
                       est.2stg = est3.2stg,
                       est.naive = est3.naive,
                       est.2stg.sd = est3.2stg.sd,
                       est.2stg.lower = est3.2stg.lower,
                       est.2stg.upper = est3.2stg.upper,
                       ifb.cov = target$psi.true >= est3.2stg.lower & target$psi.true <= est3.2stg.upper,
                       sample.size = sample.size.for.sim
  )
  
  result.2stg4 <- list(true.psi = target$psi.true,
                       est.2stg = est4.2stg,
                       est.naive = est4.naive,
                       est.2stg.sd = est4.2stg.sd,
                       est.2stg.lower = est4.2stg.lower,
                       est.2stg.upper = est4.2stg.upper,
                       ifb.cov = target$psi.true >= est4.2stg.lower & target$psi.true <= est4.2stg.upper,
                       sample.size = sample.size.for.sim
  )

  estimator_TTS <- run_TTS(dataset,estimator$cf_inds, alpha,rf.min.node.size, rf.num.trees)
  est0.TTS <- estimator_TTS$if.est.psi0
  est0.TTS.sd <- estimator_TTS$if.sd.psi0
  est0.TTS.lower <- estimator_TTS$if.lower.psi0
  est0.TTS.upper <- estimator_TTS$if.upper.psi0
  
  result.TTS0 <- list(true.psi = target$psi.true,
                      est.TTS = est0.TTS,
                      est.TTS.sd = est0.TTS.sd,
                      est.TTS.lower = est0.TTS.lower,
                      est.TTS.upper = est0.TTS.upper,
                      ifb.cov = target$psi.true >= est0.TTS.lower & target$psi.true <= est0.TTS.upper,
                      sample.size = sample.size.for.sim)
  
  
  return(list(result.2stg0 = result.2stg0,
              result.2stg1 = result.2stg1,
              result.2stg2 = result.2stg2,
              result.2stg3 = result.2stg3,
              result.2stg4 = result.2stg4,
              result.TTS0 = result.TTS0,
              params = params))
}


l <- as.integer( Sys.getenv("SGE_TASK_ID") )
if ( is.na(l) ) l <- 1

sim.result <- run_sim(l)

file_name <- paste("sim_result_size1000/job",l,".rds", sep = "")
saveRDS(sim.result, file_name)
