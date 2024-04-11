#library(dplyr)
#library(leaps)
library(pcalg)
#library(glmnet)
##' rmvDAG2(N, randDAGobj):
##' a function that does the same thing as the pcalg::rmvDAG function
##' but the input DAG is not necessarily topologically ordered
##' @param N: number of samples to be drawn
##' @param randDAGobj: a graph object generated from the pcalg::randDAG function
##' @return a Gaussian dataset
rmvDAG2 <- function(N, randDAGobj) {
  AM <- as(randDAGobj, "matrix")
  sorted_ind <- ggm::topOrder((AM != 0))
  n <- nrow(AM)
  data <- matrix(nrow = N,ncol = n)
  for (j in sorted_ind) {
    parentnodes <- which(AM[,j] != 0)
    lp <- length(parentnodes)
    switch (as.character(lp),
            "0" = {data[,j] <- rnorm(N)},
            "1" = {data[,j] <- rnorm(N, mean = data[,parentnodes] * AM[parentnodes,j], sd = 1)},
            {data[,j] <- rnorm(N, mean = data[,parentnodes] %*% AM[parentnodes,j], sd = 1)}
    )
  }
  return(data)
}


##' cutfun(L, c):
##' a function that simulates the cell probabilities from a symmetric Dirichlet distribution
##' @param L: number of ordinal levels
##' @param c: Dirichlet concentration parameter
##' @return a list of probabilities of length L, summing up to 1
cutfun <- function(L,c) {
  p <- gtools::rdirichlet(1,rep(c,L))
  return(qnorm(cumsum(p)[1:(L-1)]))
}


##' mywFUN(m):
##' a function that samples the edge weights uniformly from the interval (-1,-0.4) U (0.4,1)
##' @param m: number of edges in the DAG
##' @return m edge weights
mywFUN <- function(m) {
  return(replicate(m,mywFUNhelper()))
}
mywFUNhelper <- function() {
  y <- runif(1, 0, 1.2)
  if( y < 0.6 ){
    x <- -1 + y
  }else{
    x <- 0.4 + y - 0.6
  }
  return(x)
}


##' convertToOrdinal(scaled_data, exp_levels, concent_param):
##' a function that converts standardized Gaussian data into ordinal data
##' @param scaled_data: Gaussian dataset with each dimension standardized
##' @param exp_levels: expected number of ordinal levels
##' @param concent_param: Dirichlet concentration parameter
##' @return an ordinal dataset
convertToOrdinal <- function(scaled_data, exp_levels = 5,concent_param = 2) {
  n <- ncol(scaled_data)
  if (exp_levels == 2) {
    ordinal_levels <- replicate(n,2)
  } else {
    #ordinal_levels <- replicate(n,sample(c(2:(2 * exp_levels - 2)),1))
    #ordinal_levels <- replicate(n,sample(c(3:(2 * exp_levels - 3)),1))
    ordinal_levels <- replicate(n,exp_levels)
  }
  ordinal_data <- scaled_data
  for (i in c(1:n)) {
    
    check_levels <- ordinal_levels[i] - 1
    while (check_levels != ordinal_levels[i]) {
      cuts <- c(-Inf,
                cutfun(ordinal_levels[i],concent_param),
                Inf)
      temp <- cut(scaled_data[,i], simplify2array(cuts), labels = FALSE) - 1
      check_levels <- length(unique(temp))
    }
    ordinal_data[,i] <- temp
    
  }
  colnames(ordinal_data) <- c(1:n)
  return(ordinal_data)
}


generateOrdinal <- function(N, n, trueDAG, exp_levels = 4, concent_param = 2) {
  
  hidden_data <- rmvDAG2(N, trueDAG)
  scaled_data <- t(t(hidden_data) - apply(hidden_data,2,mean))
  truecov <- trueCov(trueDAG)
  D <- diag(sqrt(diag(truecov)))
  D.inv <- chol2inv(chol(D))
  trueSigma <- D.inv %*% truecov %*% D.inv
  scaled_data <- t(D.inv %*% t(scaled_data))
  
  # Convert the Gaussian dataset into an ordinal dataset
  ordinal_data <- convertToOrdinal(scaled_data, exp_levels = exp_levels, concent_param = concent_param)
  
  return(ordinal_data)
  
}




#trueDAG <- randDAG(n = 30, d = 2, method = "er", wFUN = list(mywFUN))

# Convert the Gaussian dataset into an ordinal dataset
#ordinal_data <- generateOrdinal(500, 30, trueDAG, exp_levels = 5, concent_param = 2)
#ordinal_data_df <- as.data.frame(ordinal_data)
#ordinal_data_df[] <- lapply(ordinal_data_df[], as.ordered)
#ordinal_levels <- apply(ordinal_data, 2, function(x) length(unique(x)))