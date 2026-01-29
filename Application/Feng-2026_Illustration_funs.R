###############################################
## Causal inference in nonlinear factor models#
######## Supporting functions #################
##### for empirical illustration ##############
# Note: some functions are slightly different #
#     from those for simulation             ###
###############################################
# Searching KNN
dist.inf <-  function(v, A2) {
    tmp   <- abs(A2 - v[-1])
    diag(tmp)  <- -1        # implicitly delete diagonal elements
    tmp[v[1],] <- -1        # implicitly delete the v[1]th element
    tmp.d <- colMaxs(tmp, value = T)
    return(tmp.d)
}

dist.knn <- function(A) {
    p    <- ncol(A)
    A2   <- tcrossprod(A)/p   
    dis  <- apply(cbind(1:nrow(A), A2), 1, function(v) dist.inf(v=v, A2=A2))
    return(dis)      # n by n distance matrix
}

findknn <- function(v, A2, K) {
     tmp   <- abs(A2 - v[-1])
     diag(tmp)  <- -1        # implicitly delete diagonal elements
     tmp[v[1],] <- -1        # implicitly delete the v[1]th element
     tmp.d <- colMaxs(tmp, value = T)
     out   <- ifelse(tmp.d <= nth(tmp.d, K), TRUE, FALSE) 
     return(out)
}

knn.index <- function(A, K=NULL) {    # A: n by p; # K: tuning param.
    p    <- ncol(A)
    A2   <- tcrossprod(A)/p
    Kmat <- apply(cbind(1:nrow(A), A2), 1, function(v) findknn(v=v, A2=A2, K=K))
    return(Kmat)           # neighborhood saved in each column
}


# sel.r <- function(sval, n) {
#   d <- length(sval)
#   ratio <- (sval[-d] / sval[-1]) < log(log(n))
#   if (any(ratio)) {
#     r <- max(which.max(ratio)-1, 1)
#   } else {
#     r <- d - 1
#   }
#   return(r)
# }

sel.r <- function(sval, n, p, sd) {
  check <- (sval > log(log(n)) * (sqrt(n)+sqrt(p)) * sd)   # make it proportional to sd of data
  d <- length(sval)
  if (all(check)) {
    r <- d
  } else {
    r <- max(which.min(check)-1, 1)   # at least 1 factor
  }
  return(r)
}

# Local PCA (for each point)
# A: n by p, input matrix; index: n by 1, KNN index; nlam: MAXimum number of vectors to be extracted; i: unit of interest
lpca <- function(A, index, nlam, n, K) {      
  A <- cbind(1:n, A)   # add a row number
  A <- A[index,]
  svd <- irlba(A[,-1], nu=nlam, nv=nlam)
  # pos <- which(A[,1]==i)
  r <- sel.r(svd$d, K, ncol(A)-1, sd(A[,-1]))
  return(list(Lam=svd$u[,1:r,drop=F], no.neigh=A[,1], sv=svd$d, d.i=r))    # Lam: K by r matrix
}


pred.cons <- function(i, y, d, subset=NULL, kmat, shutdown.d=F) {   
    y.fitval <- d.fitval <- NA
    
    index  <- kmat[,i]   # of length n
    sub    <- index & subset   # of length n
    y.fitval <- mean(y[sub])
    if (!shutdown.d) {
      d.fitval <- mean(d[index])
    }
    return(c(y.fitval, d.fitval))   # two scalars
}

# prediction function
pred <- function(i, y, d, x, subset=NULL, kmat, nlam, n, p, ctrlvar=NULL, shutdown.d=F, K) {    
    
    y.fitval <- d.fitval <- NA
    
    index <- kmat[,i]   # length n
    PC    <- lpca(A=x[,(p/2+1):p], index=index, nlam=nlam, n=n, K=K)

    # prepare design
    sub    <- index & subset   # of length n
    y.sub  <- y[sub]
    design <- cbind(1, PC$Lam[sub[index],,drop=F])
    #design <- PC$Lam[sub[index],,drop=F]
    if (!is.null(ctrlvar)) {
      design <- cbind(design, ctrlvar[sub,])
    }
    
    # eval
    design.d <- cbind(1, PC$Lam)
    eval     <- design.d[which(PC$no.neigh==i),]
    #eval <- PC$Lam[which(PC$no.neigh==i),]
    if (!is.null(ctrlvar)) {
      eval <- c(eval, ctrlvar[i,])
    }
    
    # run a local regression of y at the ith obs.
    y.fitval <- sum(.lm.fit(design, y.sub)$coeff * eval)
    
    # local logit for treatment model
    if (!shutdown.d) {
        if (!is.null(ctrlvar)) {
          design.d <- cbind(design.d, ctrlvar[index,])
        }  
        model    <- glm.fit(x=design.d, y=d[index], family=binomial())
        d.fitval <- model$family$linkinv(sum(model$coeff * eval))      
        
        # run a local regression of d at ith obs.
        # if (is.null(ctrlvar)) {
        #   d.fitval <- mean(d[index])
        # } else {
        #   design.d <- cbind(PC$Lam, ctrlvar[index,])
        #   model    <- glm.fit(x=design.d, y=d[index], family=binomial())
        #   d.fitval <- model$family$linkinv(sum(model$coeff * eval))
        # }
       
        # trimming
        #if (d.fitval<0) d.fitval <- 0
        #if (d.fitval>1) d.fitval <- 1
    }
    
    return(c(y.fitval, d.fitval))   # two scalars
}


# Computation
compute <- function(range, y, d, x, subset=NULL, K, nlam, n, p, const=F, ctrlvar=NULL, shutdown.d=F) {
    if (const) {
      kmat <- knn.index(A=x, K=K)
      fit  <- sapply(range, function(i) pred.cons(i=i, y=y, d=d, subset=subset, kmat=kmat)) 
    } else {
      kmat <- knn.index(A=x[,1:(p/2)], K=K)
      fit  <- sapply(range, function(i) pred(i=i, y=y, d=d, x=x, subset=subset,
                                             kmat=kmat, nlam=nlam, n=n, p=p, ctrlvar=ctrlvar, shutdown.d=shutdown.d, K=K))
    }
    return(list(yfit=fit[1,], ps=fit[2,]))
}
