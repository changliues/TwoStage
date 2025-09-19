source("minimax_core.R")
source("twostage_dataset.R")
source("twostage_minimax.R")
source("twostage_inference_real_dat.R")
source("twostage_TTS_real_dat.R")

library(caret)
library(plyr)

alpha <- 0.05
critical_t <- qnorm(1 - alpha/2)

# functions to be used
expit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}

###############################################read data###############################################
analysis_dataset_complete <- readRDS(file = "RealData/analysis_dataset_complete.rds")
# analysis_dataset_complete <- analysis_dataset_complete[1:100,]
sample.size <- nrow(analysis_dataset_complete)
###############################################proxci###############################################
split_ratio <- 0.25
n.flds <- 4
sample.size.for.tuning <- sample.size/n.flds*(1-split_ratio)

# # gammas and lambdas sequences
gammas1 <- c(0.08,0.1,0.12)
lambdash1 <- c(6,6.25,6.5)/sample.size.for.tuning^0.8
lambdaspi1 <- c(4600,4700,4800)/sample.size.for.tuning^3.6/lambdash1
gammas2 <- c(0.08,0.1,0.12)
lambdash2 <- c(6,6.25,6.5)/sample.size.for.tuning^0.8
lambdaspi2 <- c(4600,4700,4800)/sample.size.for.tuning^3.6/lambdash2
rf.min.node.size <- 10
rf.num.trees <- 1000

run_real_dat <- function(n.sim){
  analysis_dataset_complete.new <- analysis_dataset_complete[sort(sample(sample.size,sample.size, replace = T)),]
  dataset <- TwostageData$new(X = as.matrix(analysis_dataset_complete.new[,1:11]), A = as.numeric(analysis_dataset_complete.new$obesity), M = analysis_dataset_complete.new$LBXGH, Y = analysis_dataset_complete.new$MCQ160C, Xs = as.matrix(analysis_dataset_complete.new[,1:11]))
  estimator_2stg <- TwostageInference_real_dat$new(
    dataset,
    crossfit_folds = n.flds,
    lambdash1 =lambdash1,
    lambdaspi1 = lambdaspi1,
    gammas1 = gammas1,
    lambdash2 =lambdash2,
    lambdaspi2 = lambdaspi2,
    gammas2 = gammas2,
    rf.min.node.size = rf.min.node.size,
    rf.num.trees = rf.num.trees,
    split_ratio = split_ratio,
    print_best_params = T)
  est_2stg <- estimator_2stg$estimator()
  
  
  result_2stg <- list(est = est_2stg$estimates,
                      est.sd = est_2stg$estimates.sd,
                      est_naive = est_2stg$estimates_naive,
                      sample.size = sample.size)
  
  est_TTS <- run_TTS_real_dat(dataset,estimator_2stg$cf_inds, alpha,rf.min.node.size, rf.num.trees)
  
  result_TTS <- list(est = est_TTS$if.est.psi10,
                     est.sd = est_TTS$if.sd.psi10,
                     est.psi11 = est_TTS$psi11,
                     est.psi00 = est_TTS$psi00)
  
  return(list(result_2stg = result_2stg,
              result_TTS = result_TTS))
  
}


l <- as.integer( Sys.getenv("SGE_TASK_ID") )
if ( is.na(l) ) l <- 1

real_dat_result <- run_real_dat(l)

file_name <- paste("real_dat_result_C/job",l,".rds", sep = "")
saveRDS(real_dat_result, file_name)
