######### Taken from the source code of the pcalg package
## this function takes as parameter the adjacency matrix of a pdag (or cpdag)
## and returns the pattern of this pdag in the Meek sense, that is,
## it returns the adjacency matrix of the graph with the same skeleton where the only oriented
## edges are the v-structures (can be easily modified to work for MAGs/PAGs)
getPattern <- function(amat){
  
  ## makes the whole graph undirected
  tmp <- amat + t(amat)
  tmp[tmp == 2] <- 1
  
  ## find all v-structures i -> k <- j s.t. i not adj to k
  ## and make only those edges directed
  for (i in 1: (length(tmp[1,])-1)) {
    for (j in (i+1): length(tmp[1,])){
      if ((amat[j,i] ==0) & (amat[i,j] ==0) & (i!=j)){ ## if i no adjacent with j in G
        
        possible.k <- which(amat[i,]!= 0 & amat[,i]==0) ## finds all k such that i -> k is in G
        
        if (length(possible.k)!=0){    ## if there are any such k's then check whether j -> k for any of them
          for (k in 1: length(possible.k)){
            if ((amat[j,possible.k[k]] ==1) & (amat[possible.k[k],j]==0)) { ## if j -> k add the v-struc orientation to tmp
              tmp[possible.k[k],i] <- 0
              tmp[possible.k[k],j] <- 0
            }
          }
        }
      }
    }
  }
  tmp
}

convert_to_skeleton <- function(DAG) {
  DAG <- DAG + t(DAG)
  return(DAG != 0)
}


##' comparePatterns(estDAG,trueDAG):
##' a function that compares the patterns of two DAGs in the sense of Meek (1995)
##' @param estDAG: estimated DAG (need to be a matrix or a graphNEL object)
##' @param trueDAG: true DAG (need to be a matrix or a graphNEL object)
##' @param hardP2P: (default: FALSE) An edge in the estimated pattern is counted as 1 true positive,
##' if it has exactly the same direction (directed/undirected) as the corresponding edge in the true pattern.
##' Otherwise, it is counted as 1 false positive. When FALSE, an edge in the estimated pattern
##' is counted as 0.5 true positive and 0.5 false positive,
##' if exactly one of this edge and the corresponding edge in the true pattern is undirected.
##' @return an array of metrics
comparePatterns <- function(estDAG, trueDAG, hardP2P = FALSE, skeleton = FALSE, noPattern = FALSE) {
  
  if (skeleton) { # use skeletons
    trueSkel <- convert_to_skeleton(as(trueDAG, "matrix"))
    estSkel <- convert_to_skeleton(as(estDAG, "matrix"))
    
    temp1 <- estSkel[upper.tri(estSkel)]
    temp2 <- trueSkel[upper.tri(trueSkel)]
    
    pred_P <- sum(temp1 != 0)
    true_P <- sum(temp2 != 0)
    true_N <- sum(temp2 == 0)
    
    TP <- sum((temp1 != 0) * (temp2 != 0))
    FP <- pred_P - TP
    FN <- sum((temp1 == 0) * (temp2 != 0))
    TN <- sum((temp1 == 0) * (temp2 == 0))
    
  } else { # use patterns
    # Convert estimated DAG to CPDAG
    if (is.matrix(estDAG)) {
      estPDAG <- dag2cpdag(estDAG)
    } else if (class(estDAG) == "pcAlgo") {
      estPDAG <- t(as(estDAG,"matrix"))
    } else if (class(estDAG) == "graphNEL") {
      estPDAG <- dag2cpdag(as(estDAG,"matrix"))
    }
    
    # Convert true DAG to CPDAG
    if (is.matrix(trueDAG)) {
      truePDAG <- dag2cpdag(trueDAG)
    } else if (class(trueDAG) == "graphNEL") {
      truePDAG <- dag2cpdag(as(trueDAG,"matrix"))
    } else {
      stop("Please check if the true DAG is indeed a DAG...")
    }
    
    # Convert CPDAGs to patterns
    if (noPattern) {
      truePattern <- (as(trueDAG,"matrix") != 0)*1
      estPattern <- estDAG
    } else {
      truePattern <- getPattern(truePDAG)
      #truePattern <- trueDAG
      estPattern <- getPattern(estPDAG)
    }
    
    # 0: no edge; 1: -->; 2: <--; 3: <-->
    temp1 <- estPattern[upper.tri(estPattern)] + 2 * t(estPattern)[upper.tri(t(estPattern))]
    temp2 <- truePattern[upper.tri(truePattern)] + 2 * t(truePattern)[upper.tri(t(truePattern))]
    
    # Number of edges in the estimated pattern
    pred_P <- sum(temp1 != 0)
    
    # Number of edges in the true pattern
    true_P <- sum(temp2 != 0)
    
    # Number of non-edges in the true pattern
    true_N <- sum(temp2 == 0)
    
    # TP, FP, TN, FN, SHD
    if (hardP2P) {
      TP <- sum((temp1 != 0) * (temp1 == temp2))
    } else {
      TP <- sum((temp1 != 0) * (temp1 == temp2)) + 0.5 * sum((temp1 * temp2) == 3) + 0.5 * sum((temp1 * temp2) == 6)
    }
    FP <- pred_P - TP
    FN <- sum((temp1 == 0) * (temp2 != 0))
    TN <- sum((temp1 == 0) * (temp2 == 0))
  }
  
  # Structural hamming distance
  SHD <- FP + FN
  
  # Precision
  if ((TP + FP) == 0) {
    Precision <- 0
  } else {
    Precision <- TP / (TP + FP)
  }
  
  # TPR, FPR_P, FPR_N
  if (true_P == 0) { # true graph is empty
    if (FP >= 0) {
      TPR <- 0
      FPR_P <- 1
    } else {
      TPR <- 1
      FPR_P <- 0
    }
  } else { # true graph is non-empty
    TPR <- TP / true_P
    FPR_P <- FP / true_P
  }
  
  if (true_N == 0) { # true graph is full
    FPR_N <- 0
  } else { # true graph is not full
    FPR_N <- FP / true_N
  }
  
  compPattern <- c(SHD,TP,FP,TN,FN,Precision,TPR,FPR_N,FPR_P)
  names(compPattern) <- c("SHD","TP","FP","TN","FN","Precision","TPR","FPR_N","FPR_P")
  return(round(compPattern,2))
}

##' getSHD(estDAG, trueDAG):
##' a function that computes the structural difference between the PDAGs of two DAGs
##' @param estDAG: adjacency matrix of the estimated DAG
##' @param trueDAG: adjacency matrix of the true DAG
##' @return structural difference between two PDAGs
getSHD <- function(estDAG, trueDAG) {
  
  if (!is.matrix(estDAG)) {
    estDAG <- as(estDAG,"matrix") * 1
  }
  
  if (!is.matrix(trueDAG)) {
    trueDAG <- as(trueDAG,"matrix") * 1
  }
  
  DAGdiff <- dag2cpdag(estDAG) != dag2cpdag(trueDAG)
  return(sum(as.logical(DAGdiff + t(DAGdiff)))/2)
}