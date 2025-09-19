library(caret)

# TwostageData Class
TwostageData <- setRefClass(
  "TwostageData",
  fields = list(
    X = "matrix",
    A = "numeric",
    M = "numeric",
    Y = "numeric",
    Xs = "matrix",
    n = "numeric",
    split_indices = "list"
  ),
  
  methods = list(
    initialize = function(..., X = matrix(), A = numeric(), M = numeric(), Y = numeric(), Xs = matrix()) {
      .self$X <<- X
      .self$A <<- A
      .self$M <<- M
      .self$Y <<- Y
      .self$Xs <<- Xs
      .self$n <<- nrow(X)
     },
    
    r_pi1 = function() {
      return(.self$X)
    },
    
    r_h1 = function() {
      return(.self$X)
    },
    
    r_pi2 = function() {
      return(cbind(.self$M, .self$X))
    },
    r_h2 = function() {
      return(cbind(.self$M, .self$X))
    },
    
    create_crossfit_split = function(n_splits) {
      if (n_splits <= 1) {
        idx <- seq_len(.self$n)
        split_indices <<- list(list(train = idx, eval = idx))
      } else {
        folds <- createFolds(seq_len(.self$n), n_splits)
        split_indices <<- list()
        for (i in seq_len(n_splits)) {
          train_idx <- unlist(folds[-i])
          eval_idx <- folds[[i]]
          split_indices[[i]] <<- list(train = train_idx, eval = eval_idx)
        }
      }
      
      return(split_indices=split_indices)
    }
  )
)
