MBOCDmod <- function(dataset, alpha = 0.05, test = "zf", method = "probit") {
  if (!(test == "zf")) {
    dataset[] = lapply(dataset[], as.ordered)
  }
  mb = learn.mb(dataset, alpha = alpha, test = test)
  skel = skel_recov(dataset, mb, alpha = alpha, test = test)
  patternAndVstru = learn.vstruct(dataset, mb, test = test, alpha = alpha)
  if (test == "zf") {
    dataset[] = lapply(dataset[], as.ordered)
  }
  mbocd = oBNSearch(dataset, pattern = patternAndVstru$pattern, gam = patternAndVstru$vstruct, method = method)
  return(mbocd$gam)
}


learn.mb <- function(dataset,alpha = 0.05 , test = "zf", cluster = NULL,
                     max.sx = ncol(dataset), debug = FALSE){
  data.info = bnlearn:::check.data(dataset, allow.missing = TRUE)
  complete=data.info$complete.nodes
  
  `adjecencies` <- function(dataset, alpha = alpha){
    sepset <- lapply(colnames(dataset), function(.) {vec<-vector("list",ncol(dataset))
    names(vec)<-colnames(dataset)
    return(vec)})
    names(sepset)<-colnames(dataset)
    Adjs <-vector("list",ncol(dataset))
    names(Adjs)<-colnames(dataset)
    for (T in colnames(dataset)){
      #message("T = ", T)
      S = character(0) # empty set
      adjT = character(0) # empty set
      for (v in setdiff(colnames(dataset),T)) {
        pval <- bnlearn:::indep.test(v, T, sx = S, test = test, data = dataset,
                                     B = 0L, alpha = alpha, complete = complete)
        # #message("v = ",v)
        # #message("pval = ",pval)
        if (abs(pval[1]) >= alpha){
          sepset[[T]][[v]] <- S
        }else{
          adjT <- append(adjT, v)
        }
      }
      # Sort adj(T ) in increasing corr(Vi, T ) value
      # sorted_adj <- sort(corr[,T][adjT])
      # #print(sorted_adj)
      # #print(names(sorted_adj))
      nonPC = character(0)
      # sorted_adj <- names(sorted_adj)
      k = 1
      nNum = min(length(adjT),max.sx)
      if(nNum < 1){
        Adjs[[T]] <- character(0)
      }else{
        while(k <= nNum){
          # all subsets S \subset adj(T) needs to be tested whose
          # sizes do not exceed k
          for (v in adjT) {
            #message("v = ", v)
            condset=setdiff(adjT,v) 
            len <- length(condset)
            if(k > len){
              k <- nNum
              break
            }
            #for (j in 1:k) {
            a <- bnlearn:::allsubs.test(v, T, sx=condset, fixed = character(0), data=dataset, test=test, B = 0L,
                                        alpha = alpha, min = k, max = k, complete=complete, debug = FALSE)
            #print(a)
            if(a["p.value"] >= alpha){
              # adjT <- setdiff(adjT,v)
              nonPC = union(nonPC, v)
              sepset[[T]][[v]] <- attr(a,"dsep.set")
            }
            #} 
          }
          if (length(nonPC) > 0) {
            adjT = setdiff(adjT, nonPC)
            nonPC = character(0)
          }
          k <- k+1
          nNum <- length(adjT)
        }
        Adjs[[T]] <- adjT
      }
      
    }
    # message("before symmetry correction")
    # print(Adjs)
    ## symmetry correction for removing false positives from adjV(T)
    ## i.e., if X \in adjV(Y) but Y \not\in adjV(X) then
    ## adjV(Y) = adjV(Y)\{X} and sepset[[Y]][[X]] <- sepset[[X]][[Y]]
    for (var in colnames(dataset)) {
      for (u in Adjs[[var]]) {
        if(!(var %in% Adjs[[u]])){
          Adjs[[var]] <- setdiff(Adjs[[var]],u)
          sepset[[var]][[u]] <- sepset[[u]][[var]]
        }
      }
    }
    # message("after symmetry correction")
    # print(Adjs)
    return(list(adjs=Adjs,sepset = sepset))
  }#End adjecencies

  
  `findMb` <- function(T, dataset, adjs, sepset,alpha = alpha){
    adjT <- adjs[[T]]
    # bd(T) U ch(T) \subsetq Mb(T)
    mmbT <- adjT
    # chT <- character(0)
    # complex-spouses recovery phase
    # candidates of complex-spouses
    candids = character(0)
    for (v in mmbT) {
      candids = union(candids, adjs[[v]]) 
    }
    candids <- setdiff(candids,union(T,adjT))
    for (ch in adjT) {
      for (csp in candids) {
        pval1 <- bnlearn:::indep.test(csp, T, sx = sepset[[T]][[csp]], test = test, data = dataset,
                                      B = 0L, alpha = alpha, complete = complete)
        pval2 <- bnlearn:::indep.test(csp, T, sx = union(sepset[[T]][[csp]],ch), test = test, data = dataset,
                                      B = 0L, alpha = alpha, complete = complete)
        if(abs(pval1[1]) > alpha & abs(pval2[1]) <= alpha){
          mmbT <- union(mmbT,csp)
          #chT = union(chT, ch)
        }
      }
    }
    return(mmbT)
  }#End findMb
  
  # 1. [Compute Markov Blankets]
  result <- adjecencies(dataset = dataset,alpha = alpha)
  nam <- names(dataset)
  mb = bnlearn:::smartSapply(cluster, as.list(nam), findMb, dataset = dataset, 
                             adjs = result$adjs, sepset = result$sepset, alpha=alpha)
  names(mb) = nam
  p <- ncol(dataset)
  G <- matrix(0, p, p)
  colnames(G)<-rownames(G)<-nam

  # # check symmetry in the output of the algorithm
  for (var in colnames(dataset)) {
    for (u in mb[[var]]) {
      if(!(var %in% mb[[u]])){
        mb[[var]] <- setdiff(mb[[var]],u)
        #sepset[[var]][[u]] <- sepset[[u]][[var]]
      }
    }
  }#End check symmetry
  
  return(mb)
}


skel_recov <- function(dataset, mb, test = "zf", alpha = 0.05) {
  
  sepset <- lapply(colnames(dataset), function(.) {vec<-vector("list",ncol(dataset))
  names(vec)<-colnames(dataset)
  return(vec)})
  names(sepset)<-colnames(dataset)
  
  for (i in colnames(dataset)) {
    for (j in setdiff(colnames(dataset), i)) {
      if (!(i %in% mb[[j]]) & !(j %in% mb[[i]])) {
        if (length(mb[[i]]) < length(mb[[j]])) {
          sepset[[i]][[j]] = mb[[i]]
        } else {
          sepset[[i]][[j]] = mb[[j]]
        }
      }
    }
  }
  
  ## Skeleton Recovery
  skel = matrix(0, ncol(dataset), ncol(dataset))
  rownames(skel) = colnames(dataset)
  colnames(skel) = colnames(dataset)
  for (i in colnames(dataset)) {
    for (j in mb[[i]]) {
      skel[j,i] = 1
    }
  }
  
  for (i in 0:(ncol(dataset) - 2)) {
    for (u in colnames(dataset)) {
      for (v in setdiff(colnames(dataset), u)) {
        if (skel[u,v] == 1 & skel[v, u] == 1 & (sum(skel[, u] != 0) - 1) >= i) {
          
          adu = colnames(dataset)[which(skel[, u] == 1)]
          #names(adu) = NULL
          condset = setdiff(adu, v)
          a <- bnlearn:::allsubs.test(u, v, sx=condset, fixed = character(0), data=dataset, test=test, B = 0L,
                                      alpha = alpha, min = i, max = i, complete=TRUE, debug = FALSE)
          if(a["p.value"] >= alpha){
            # adjT <- setdiff(adjT,v)
            sepset[[u]][[v]] <- attr(a,"dsep.set")
            sepset[[v]][[u]] <- attr(a,"dsep.set")
            skel[u,v] = 0
            skel[v,u] = 0
          }
        }
      }
    }
  }
  
  return(list(skel=skel, sepset=sepset))
}



learn.vstruct <- function(dataset, mb, test = "zf",alpha = 0.05) {
  skel = skel_recov(dataset, mb, test = test, alpha = alpha)
  sepset = skel$sepset
  wmat = skel$skel
  vset = colnames(dataset)
  q = length(vset)
  for (i in 1:(q-1)) {
    for (j in (i+1):q) {
      for (l in 1:q) {
        u = vset[i]
        v = vset[j]
        w = vset[l]
        if (skel$skel[u,v] == 0 && wmat[u,w] == 1 && wmat[v,w] == 1 &&
            wmat[w,u] + wmat[w,v] != 0) {
          pval1 = bnlearn:::indep.test(u, v, sx = union(sepset[[u]][[v]],w), test = test, data = dataset,
                                       B = 0L, alpha = alpha, complete = TRUE)
          pval2 = bnlearn:::indep.test(u, v, sx = sepset[[u]][[v]], test = test, data = dataset,
                                       B = 0L, alpha = alpha, complete = TRUE)
          if (pval1 < alpha && pval2 > alpha) {
            wmat[w,u] = wmat[w,v] = 0
          }
        }
      }
    }
  }
  vstruct = t(skel$skel - wmat)
  return(list(pattern = wmat, vstruct = vstruct))
}


admissible = function(i, j, gam_short_old) {
  if (gam_short_old[i, j]) {
    #delete an edge is always admissible
    return(TRUE)
  } else{
    gam_short_old[i, j] = 1
    return(gRbase::is.DAG(igraph::graph_from_adjacency_matrix(gam_short_old)))
  }
}

admissible_rev = function(i, j, gam_short_old) {
  if (gam_short_old[i, j] == gam_short_old[j, i]) {
    return(FALSE)
  } else{
    tmp = gam_short_old[i, j]
    gam_short_old[j, i] = tmp
    gam_short_old[i, j] = !tmp
    return(gRbase::is.DAG(igraph::graph_from_adjacency_matrix(gam_short_old)))
  }
}

mypolr = function(formula, data, ic, method, nq_y, nq_x) {
  boolFalse = FALSE
  if (nq_y > 2) {
    if (ic == "bic") {
      tryCatch({
        #sometimes MASS::polr fails to initialize
        IC = stats::BIC(MASS::polr(formula, data = data, method = method))
        boolFalse <- TRUE
      }, error = function(e) {
        
      }, finally = {
        
      })
      while (!boolFalse) {
        tryCatch({
          IC = stats::BIC(MASS::polr(
            formula,
            data = data,
            start = sort(stats::rnorm(nq_y - 1 + sum(nq_x - 1))),
            method = method
          ))
          boolFalse <- TRUE
        }, error = function(e) {
          
        }, finally = {
          
        })
      }
    } else if (ic == "aic") {
      tryCatch({
        IC = stats::AIC(MASS::polr(formula, data = data, method = method))
        boolFalse <- TRUE
      }, error = function(e) {
        
      }, finally = {
        
      })
      while (!boolFalse) {
        tryCatch({
          IC = stats::AIC(MASS::polr(
            formula,
            data = data,
            start = sort(stats::rnorm(nq_y - 1 + sum(nq_x - 1))),
            method = method
          ))
          boolFalse <- TRUE
          
        }, error = function(e) {
          
        }, finally = {
          
        })
      }
    }
  } else{
    if (ic == "bic") {
      tryCatch({
        IC = stats::BIC(stats::glm(
          formula,
          data = data ,
          family = stats::binomial(link = method)
        ))
        boolFalse <- TRUE
      }, error = function(e) {
        
      }, finally = {
        
      })
      while (!boolFalse) {
        tryCatch({
          IC = stats::BIC(stats::glm(
            formula,
            data = data,
            start = sort(stats::rnorm(nq_y - 1 + sum(nq_x - 1))),
            family = stats::binomial(link = method)
          ))
          boolFalse <- TRUE
        }, error = function(e) {
          
        }, finally = {
          
        })
      }
    } else if (ic == "aic") {
      tryCatch({
        IC = stats::AIC(stats::glm(
          formula,
          data = data ,
          family = stats::binomial(link = method)
        ))
        boolFalse <- TRUE
      }, error = function(e) {
        
      }, finally = {
        
      })
      while (!boolFalse) {
        tryCatch({
          IC = stats::AIC(stats::glm(
            formula,
            data = data,
            start = sort(stats::rnorm(nq_y - 1 + sum(nq_x - 1))),
            family = stats::binomial(link = method)
          ))
          boolFalse <- TRUE
          
        }, error = function(e) {
          
        }, finally = {
          
        })
      }
    }
  }
  return(IC)
}


oBNSearch = function(y,pattern = NULL,gam = NULL,
                 ic = "bic",method = "probit",verbose = FALSE) {
  #hill-climbing
  
  vstru = gam
  n = nrow(y) # sample size
  q = ncol(y) # number of variables
  nq = rep(0, q)
  for (i in 1:q) {
    nq[i] = nlevels(y[, i]) # level of the ith variable
  }
  if (is.null(gam)) {
    gam = matrix(FALSE, q, q)
  } else{
    gam = (gam != 0) # boolean
  }
  
  # ind_q:4个变量形式
  # 1] 2 3 4
  # 2] 1 3 4
  # 3] 1 2 4
  # 4] 1 2 3
  ind_q = matrix(0, q, q - 1)
  for (i in 1:q) {
    if (i == 1) {
      ind_noi = 2:q
    } else if (i == q) {
      ind_noi = 1:(q - 1)
    } else{
      ind_noi = c(1:(i - 1), (i + 1):q)
    }
    ind_q[i, ] = ind_noi
  }
  
  iter = 0
  ic_improv = 1
  act_ind = c(NA, NA)
  state = "add" # or "del"
  if (ic == "bic") {
    ic_best = rep(0, q)
    for (i in 1:q) {
      if (sum(gam[,i]) > 0) { ########## sum(gam[i,]) > 0
        ic_best[i] = mypolr(
          y[, i] ~ .,
          data = y[, gam[,i] == 1], ############ gam[i,]
          ic = ic,
          method = method,
          nq_y = nq[i],
          nq_x = nq[gam[,i] == 1] # nq[gam[i,]]
        )
      } else{
        if (nq[i] > 2) {
          ic_best[i] = stats::BIC(MASS::polr(y[, i] ~ 1, method = method))
        } else{
          ic_best[i] = stats::BIC(stats::glm(y[, i] ~ 1, family = stats::binomial(link = method)))
        }
      }
    }
    ############################################################
    while (ic_improv > 0) {
      iter = iter + 1
      ic_improv = -Inf
      ic_improv_rev = rep(-Inf, 2)
      gam_new = gam
      ic_improv_new = -Inf
      ic_improv_rev_new = rep(-Inf, 2)
      ic_best_new = -Inf
      ic_rev_best_new = rep(-Inf, 2)
      for (i in 1:q) {
        for (j in 1:(q - 1)) {
          if (pattern[ind_q[i, j], i] + pattern[i, ind_q[i, j]] == 2 & admissible(ind_q[i, j], i, gam)) {# 有边或者加上边都没有影响
            if (gam[ind_q[i, j], i]) {
              #delete
              gam_new[ind_q[i, j], i] = FALSE # 如果有边就删掉边
              if (sum(gam_new[,i]) > 0) { # sum(gam_new[i,]) > 0
                ic_best_new = mypolr(
                  y[, i] ~ .,
                  data = y[, gam_new[,i] == 1],
                  ic = ic,
                  method = method,
                  nq_y = nq[i],
                  nq_x = nq[gam_new[,i] == 1] # nq[gam_new[i,]]
                )
              } else{
                if (nq[i] > 2) {
                  ic_best_new = stats::BIC(MASS::polr(y[, i] ~ 1, method = method))
                } else{
                  ic_best_new = stats::BIC(stats::glm(y[, i] ~ 1, family = stats::binomial(link = method)))
                }
              }
              ic_improv_new = ic_best[i] - ic_best_new
              if (ic_improv_new > ic_improv) { # 删边后ic值更小
                ic_improv = ic_improv_new
                act_ind = c(ind_q[i, j], i)
                state = "del"
              }
              gam_new[ind_q[i, j], i] = TRUE
            } else{
              #add
              gam_new[ind_q[i, j], i] = TRUE
              if (sum(gam_new[,i]) > 0) { # sum(gam_new[i,])
                ic_best_new = mypolr(
                  y[, i] ~ .,
                  data = y[, gam_new[,i] == 1],   ## y[, gam_new[i,]]
                  ic = ic,
                  method = method,
                  nq_y = nq[i],
                  nq_x = nq[gam_new[,i] == 1]     ### nq[gam_new[i,]]
                )
              } else{
                if (nq[i] > 2) {
                  ic_best_new = stats::BIC(MASS::polr(y[, i] ~ 1, method = method))
                } else{
                  ic_best_new = stats::BIC(stats::glm(y[, i] ~ 1, family = stats::binomial(link = method)))
                }
              }
              ic_improv_new = ic_best[i] - ic_best_new
              if (ic_improv_new > ic_improv) {
                ic_improv = ic_improv_new
                act_ind = c(ind_q[i, j], i)
                state = "add"
              }
              gam_new[ind_q[i, j], i] = FALSE
            }
          }
        }
      }
      #reverse edge
      for (i in 1:q) {
        for (j in 1:(q - 1)) {
          if (pattern[ind_q[i, j], i] + pattern[i, ind_q[i, j]] == 2  & admissible_rev(ind_q[i, j], i, gam)) {     # ind_q[i, j] -> i
            tmp = gam_new[ind_q[i, j], i] # temp = 1
            gam_new[i, ind_q[i, j]] = tmp # gam_new[i, ind_q[i, j]] = 1, i -> ind_q[i, j]
            gam_new[ind_q[i, j], i] = !tmp # gam_new[ind_q[i, j], i] = 0
            if (sum(gam_new[,i]) > 0) {
              ic_rev_best_new[1] = mypolr( # 存储指向i的
                y[, i] ~ .,
                data = y[, gam_new[,i] == 1],
                ic = ic,
                method = method,
                nq_y = nq[i],
                nq_x = nq[gam_new[,i] == 1]
              )
            } else{
              if (nq[i] > 2) {
                ic_rev_best_new[1] = stats::BIC(MASS::polr(y[, i] ~ 1, method = method))
              } else{
                ic_rev_best_new[1] = stats::BIC(stats::glm(y[, i] ~ 1, family = stats::binomial(link = method)))
              }
            }
            if (sum(gam_new[, ind_q[i, j]]) > 0) {
              ic_rev_best_new[2] = mypolr(
                y[, ind_q[i, j]] ~ .,
                data = y[, gam_new[,ind_q[i, j]] == 1],
                ic = ic,
                method = method,
                nq_y = nq[ind_q[i, j]],
                nq_x = nq[gam_new[,ind_q[i, j]] == 1]
              )
            } else{
              if (nq[ind_q[i, j]] > 2) { #######modify###########
                ic_rev_best_new[2] = stats::BIC(MASS::polr(y[, ind_q[i, j]] ~ 1, method = method))
              } else{
                ic_rev_best_new[2] = stats::BIC(stats::glm(y[, ind_q[i, j]] ~ 1, family = stats::binomial(link = method)))
              }
            }
            ic_improv_rev_new[1] = ic_best[i] - ic_rev_best_new[1]
            ic_improv_rev_new[2] = ic_best[ind_q[i, j]] - ic_rev_best_new[2]
            ic_improv_new = ic_improv_rev_new[1] + ic_improv_rev_new[2]
            if (ic_improv_new > ic_improv) {
              ic_improv = ic_improv_new
              ic_improv_rev = ic_improv_rev_new
              act_ind = c(ind_q[i, j], i)
              state = "rev"
            }
            gam_new[ind_q[i, j], i] = tmp
            gam_new[i, ind_q[i, j]] = !tmp
          }
        }
      }
      
      if (ic_improv > 0) {
        if (state == "add") {
          gam[act_ind[1], act_ind[2]] = TRUE
          ic_best[act_ind[2]] = ic_best[act_ind[2]] - ic_improv
        } else if (state == "del") {
          gam[act_ind[1], act_ind[2]] = FALSE
          ic_best[act_ind[2]] = ic_best[act_ind[2]] - ic_improv
        } else if (state == "rev") {
          tmp = gam[act_ind[1], act_ind[2]]
          gam[act_ind[2], act_ind[1]] = tmp
          gam[act_ind[1], act_ind[2]] = !tmp
          ic_best[act_ind[1]] = ic_best[act_ind[1]] - ic_improv_rev[2]
          ic_best[act_ind[2]] = ic_best[act_ind[2]] - ic_improv_rev[1]
        }
      }
      if (verbose&&iter%%1==0){
        print(paste(iter," iterations have completed",sep=""))
        print("The current DAG adjacency matrix is")
        print(gam)
        print(paste("with ",  ic, " = ",sum(ic_best),sep=""))
      }
    }
  } else if (ic == "aic") {
    
    ic_best = rep(0, q)
    for (i in 1:q) {
      if (sum(gam[,i]) > 0) { ########## sum(gam[i,]) > 0
        ic_best[i] = mypolr(
          y[, i] ~ .,
          data = y[, gam[,i] == 1], ############ gam[i,]
          ic = ic,
          method = method,
          nq_y = nq[i],
          nq_x = nq[gam[,i] == 1] # nq[gam[i,]]
        )
      } else{
        if (nq[i] > 2) {
          ic_best[i] = stats::AIC(MASS::polr(y[, i] ~ 1, method = method))
        } else{
          ic_best[i] = stats::AIC(stats::glm(y[, i] ~ 1, family = stats::binomial(link = method)))
        }
      }
    }
    ############################################################
    while (ic_improv > 0) {
      iter = iter + 1
      ic_improv = -Inf
      ic_improv_rev = rep(-Inf, 2)
      gam_new = gam
      ic_improv_new = -Inf
      ic_improv_rev_new = rep(-Inf, 2)
      ic_best_new = -Inf
      ic_rev_best_new = rep(-Inf, 2)
      for (i in 1:q) {
        for (j in 1:(q - 1)) {
          if (colnames(y)[ind_q[i, j]] %in% mb[[colnames(y)[i]]] & admissible(ind_q[i, j], i, gam)) {# 有边或者加上边都没有影响
            if (!gam[ind_q[i, j], i]) {
              #add
              gam_new[ind_q[i, j], i] = TRUE
              if (sum(gam_new[,i]) > 0) { # sum(gam_new[i,])
                ic_best_new = mypolr(
                  y[, i] ~ .,
                  data = y[, gam_new[,i] == 1],   ## y[, gam_new[i,]]
                  ic = ic,
                  method = method,
                  nq_y = nq[i],
                  nq_x = nq[gam_new[,i] == 1]     ### nq[gam_new[i,]]
                )
              } else{
                if (nq[i] > 2) {
                  ic_best_new = stats::AIC(MASS::polr(y[, i] ~ 1, method = method))
                } else{
                  ic_best_new = stats::AIC(stats::glm(y[, i] ~ 1, family = stats::binomial(link = method)))
                }
              }
              ic_improv_new = ic_best[i] - ic_best_new
              if (ic_improv_new > ic_improv) {
                ic_improv = ic_improv_new
                act_ind = c(ind_q[i, j], i)
                state = "add"
              }
              gam_new[ind_q[i, j], i] = FALSE
            }
          }
        }
      }
      #reverse edge
      for (i in 1:q) {
        for (j in 1:(q - 1)) {
          if (colnames(y)[ind_q[i, j]] %in% mb[[colnames(y)[i]]] & admissible_rev(ind_q[i, j], i, gam)) {     # ind_q[i, j] -> i
            tmp = gam_new[ind_q[i, j], i] # temp = 1
            gam_new[i, ind_q[i, j]] = tmp # gam_new[i, ind_q[i, j]] = 1, i -> ind_q[i, j]
            gam_new[ind_q[i, j], i] = !tmp # gam_new[ind_q[i, j], i] = 0
            if (sum(gam_new[,i]) > 0) {    
              ic_rev_best_new[1] = mypolr( # 存储指向i的
                y[, i] ~ .,
                data = y[, gam_new[,i] == 1],
                ic = ic,
                method = method,
                nq_y = nq[i],
                nq_x = nq[gam_new[,i] == 1]
              )
            } else{
              if (nq[i] > 2) {
                ic_rev_best_new[1] = stats::AIC(MASS::polr(y[, i] ~ 1, method = method))
              } else{
                ic_rev_best_new[1] = stats::AIC(stats::glm(y[, i] ~ 1, family = stats::binomial(link = method)))
              }
            }
            if (sum(gam_new[, ind_q[i, j]]) > 0) {
              ic_rev_best_new[2] = mypolr(
                y[, ind_q[i, j]] ~ .,
                data = y[, gam_new[,ind_q[i, j]] == 1],
                ic = ic,
                method = method,
                nq_y = nq[ind_q[i, j]],
                nq_x = nq[gam_new[,ind_q[i, j]] == 1]
              )
            } else{
              if (nq[ind_q[i, j]] > 2) { #######modify###########
                ic_rev_best_new[2] = stats::AIC(MASS::polr(y[, ind_q[i, j]] ~ 1, method = method))
              } else{
                ic_rev_best_new[2] = stats::AIC(stats::glm(y[, ind_q[i, j]] ~ 1, family = stats::binomial(link = method)))
              }
            }
            ic_improv_rev_new[1] = ic_best[i] - ic_rev_best_new[1]
            ic_improv_rev_new[2] = ic_best[ind_q[i, j]] - ic_rev_best_new[2]
            ic_improv_new = ic_improv_rev_new[1] + ic_improv_rev_new[2]
            if (ic_improv_new > ic_improv) {
              ic_improv = ic_improv_new
              ic_improv_rev = ic_improv_rev_new
              act_ind = c(ind_q[i, j], i)
              state = "rev"
            }
            gam_new[ind_q[i, j], i] = tmp
            gam_new[i, ind_q[i, j]] = !tmp
          }
        }
      }
      
      if (ic_improv > 0) {
        if (state == "add") {
          gam[act_ind[1], act_ind[2]] = TRUE
          ic_best[act_ind[2]] = ic_best[act_ind[2]] - ic_improv
        } else if (state == "del") {
          gam[act_ind[1], act_ind[2]] = FALSE
          ic_best[act_ind[2]] = ic_best[act_ind[2]] - ic_improv
        } else if (state == "rev") {
          tmp = gam[act_ind[1], act_ind[2]]
          gam[act_ind[2], act_ind[1]] = tmp
          gam[act_ind[1], act_ind[2]] = !tmp
          ic_best[act_ind[1]] = ic_best[act_ind[1]] - ic_improv_rev[2]
          ic_best[act_ind[2]] = ic_best[act_ind[2]] - ic_improv_rev[1]
        }
      }
      if (verbose&&iter%%1==0){
        print(paste(iter," iterations have completed",sep=""))
        print("The current DAG adjacency matrix is")
        print(gam)
        print(paste("with ",  ic, " = ",sum(ic_best),sep=""))
      }
    }
  }
  return(list(gam = gam * 1, ic_best = sum(ic_best)))
}