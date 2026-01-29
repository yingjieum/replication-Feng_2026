###########################################
### Define supporting functions ###########
####### for simulation ####################
####### Last updated: 01/29/2026 ##########
###########################################

# treatment model
mu0 <- function(alpha) {
    return(alpha + alpha^2)
}
mu1 <- function(alpha) {
    return(2*alpha + alpha^2 + 1)
}

mu0.m <- function(alpha1, alpha2) {
  return(2 + 10*sin(pi*alpha1*alpha2) + 2*(alpha2 - 0.5)^2 + 2*alpha1)
}
mu1.m <- function(alpha1, alpha2) {
  tau <- 1 + 0.6*(alpha1 - 0.5) + 0.4*sin(2*pi*alpha2) + 0.3*(alpha1 - 0.5)*(alpha2 - 0.5)
  return(mu0.m(alpha1, alpha2) + tau)
}

px <- function(alpha) {
  alpha <- alpha - 0.5
  return(exp(alpha+alpha^2)/(1 + exp(alpha+alpha^2)))
}

px.m <- function(alpha1, alpha2) {
  ps <- 1/(1+exp(-(-0.4 + 1.2*(alpha1 - 0.5) + 1.2*(alpha2 - 0.5) +
                   0.8*sin(2*pi*alpha1) - 0.8*(alpha2 - 0.5)^2)))
  return(ps)
}

# large-dim x
# 1D nonlinear factor model
hdfun0 <- function(u,v) {
    out <- outer(u, v, "*")
    return(out)
}

# rank-3, exactly
hdfun1 <- function(u,v) {
  s   <- outer(u, v, "*")
  out <- -4.5*s + 18.4*s^2 - 16.2*s^3
  return(out)
}

# approximately rank-3
hdfun2 <- function(u, v, sig=.1) {
  s <- outer(u, v, "-")
  out <- 2*exp(-s^2/sig)
  return(out)
}

# the following 2 is the same as hdfun2, but with different smoothness
hdfun7 <- function(u, v, sig=.5) {
  s <- outer(u, v, "-")
  out <- 2*exp(-s^2/sig)
  return(out)
}

hdfun8 <- function(u, v, sig=.02) {
  s <- outer(u, v, "-")
  out <- 2*exp(-s^2/sig)
  return(out)
}

# approximately high-rank >=5
hdfun3 <- function(u, v, sig=.1) {
  s <- outer(u, v, "-")
  out <- 2*exp(-abs(s)/sig)
  return(out)
}

# 2D nonlinear factor model
# rank-9, approximately rank-4
hdfun4 <- function(u1, u2, v1, v2) {
  s <- outer(2*(u1-0.5), 2*(v1-0.5), "*") + outer(2*(u2-0.5), 2*(v2-0.5), "*")
  out <- 2*s + (s*s) + 0.5*(s*s*s)
  return(out)
}

# approximately high-rank >=5
hdfun5 <- function(u1, u2, v1, v2, sig=.1) {
  d1 <- outer(u1, v1, "-")
  d2 <- outer(u2, v2, "-")
  out <- 2*exp(-(d1*d1+d2*d2)/sig)
  return(out)
}

# approximately high-rank >=3
hdfun6 <- function(u1, u2, v1, v2, sig=.5) {
  d1 <- outer(u1, v1, "-")
  d2 <- outer(u2, v2, "-")
  out <- 2*exp(-(abs(d1)+abs(d2))/sig)
  return(out)
}


funlist <- list(hdfun1=hdfun1, hdfun2=hdfun2, hdfun3=hdfun3, hdfun4=hdfun4, hdfun5=hdfun5, hdfun6=hdfun6, hdfun7=hdfun7, hdfun8=hdfun8)

# error generator:
gen_err <- function(n, p, sigma=1, rho) {
  U <- matrix(0, n, p)
  
  # stationary start u_{i,1} ~ N(0,1)
  U[, 1] <- rnorm(n)
  
  sd_eps <- sigma * sqrt(1 - rho^2)
  Eps <- matrix(rnorm(n * (p - 1), sd = sd_eps), n, p - 1)
  
  for (t in 2:p) {
    U[, t] <- rho * U[, t - 1] + Eps[, t - 1]
  }
  return(U)
}


# DGP generator
# sigma: s.d. for error in x
# err: heteroskedasticity or autocorrelated error in x? 0: hom; 1 hsk; 2 AR(1)
dgp <- function(n, p, hdmodel, addz=FALSE, sigma=1, err=0, rho=0.2, dim.alp=1) {
    # latent variable
    alpha1 <- runif(n, 0, 1)
    alpha2 <- NULL
    if (dim.alp==2) {
      alpha2 <- runif(n, 0, 1)
    }
    
    z <- NULL
    if (addz) {
        z <- runif(n, 0, 1)
    }
    
    # HD covariates
    eta1 <- runif(p, 0, 1)
    eta2 <- NULL
    
    # generate error
    if (err==0) {
      error <- matrix(rnorm(n*p, 0, sigma), n, p)
    } else if (err==1) {
      error <- exp(0.1*outer(alpha1-0.5, eta1-0.5, FUN="+")) * matrix(rnorm(n*p, 0, sigma), n, p)
    } else {
      error <- gen_err(n, p, sigma, rho)     
    }
    
    if (dim.alp==1) {
      x <- hdmodel(alpha1, eta1) + error   # n by p matrix
    } else {
      eta2 <- runif(p, 0, 1)
      x <- hdmodel(alpha1, alpha2, eta1, eta2) + error
    }
    
    # outcome and propensity
    if (dim.alp==1) {
      if (err!=2) {
        Y0 <- mu0(alpha1) + rnorm(n)
        Y1 <- mu1(alpha1) + rnorm(n)
      } else {
        Y0 <- mu0(alpha1) + (rho*error[,p] + rnorm(n, 0, sigma*sqrt(1-rho^2)))
        Y1 <- mu1(alpha1) + (rho*error[,p] + rnorm(n, 0, sigma*sqrt(1-rho^2)))
      }
      prob  <- px(alpha1)
    } else {
      if (err!=2) {
        Y0 <- mu0.m(alpha1, alpha2) + rnorm(n)
        Y1 <- mu1.m(alpha1, alpha2) + rnorm(n)
      } else {
        Y0 <- mu0.m(alpha1, alpha2) + (rho*error[,p] + rnorm(n, 0, sigma*sqrt(1-rho^2)))
        Y1 <- mu1.m(alpha1, alpha2) + (rho*error[,p] + rnorm(n, 0, sigma*sqrt(1-rho^2)))
      }
      prob <- px.m(alpha1, alpha2)
    }
    
    runis <- runif(n,0,1)
    D     <- ifelse(runis < prob,1,0)
    
    y <- D*Y1+(1-D)*Y0
    
    return(list(y=y, y0=Y0, y1=Y1, x=x, z=z, alpha1=alpha1, alpha2=alpha2, eta1=eta1, eta2=eta2, d=D))
}

################################################
# Note the following KNN can accept a vector of Ks and report a list of matrices
# Searching KNN
findknn <- function(v, A2, Kseq=NULL) {
     tmp   <- abs(A2 - v[-1])
     diag(tmp)  <- -1        # implicitly delete diagonal elements
     tmp[v[1],] <- -1        # implicitly delete the v[1]th element
     tmp.d <- colMaxs(tmp, value = T)
     out   <- sapply(Kseq, function(K) ifelse(tmp.d <= nth(tmp.d, K), TRUE, FALSE)) # n by length(Kseq) matrix
     return(out)
}

knn.index <- function(A, Kseq=NULL) {    # A: n by p; # Kseq: vector of tuning param.
    p    <- ncol(A)
    A2   <- tcrossprod(A)/p
    Kmat <- apply(cbind(1:nrow(A), A2), 1, function(v) findknn(v=v, A2=A2, Kseq=Kseq))
    return(Kmat)           # a matrix: neighborhood for each unit saved in each column, each n rows -> one choice of K in Kseq 
}

# select # of local factors, n should >=2
# sel.r <- function(sval, n) {
#   d <- length(sval)
#   ratio <- (sval[-d] / sval[-1]) < log(log(n))
#   if (any(ratio)) {
#     r <- max(which.max(ratio)-1, 1)   # at least 1 factor
#   } else {
#     r <- d - 1
#   }
#   return(r)
# }

sel.r <- function(sval, n, p) {
  check <- (sval > log(log(n)) * (sqrt(n)+sqrt(p)))
  d <- length(sval)
  if (all(check)) {
    r <- d
  } else {
    r <- max(which.min(check)-1, 1)   # at least 1 factor
  }
  return(r)
}

# Local PCA (for each point)
# A: n by p, input matrix; index: n by 1, KNN index; nlam: Maximum number of vectors to be extracted
lpca <- function(A, index, nlam, n, K) {      
    A <- cbind(1:n, A)   # add a row number
    A <- A[index,]
    svd <- irlba(A[,-1], nu=nlam, nv=nlam)
    r <- sel.r(svd$d, K, ncol(A)-1)
    return(list(Lam=svd$u[,1:r,drop=F], no.neigh=A[,1], sv=svd$d, d.i=r))    #Lam: K by r
}


###### Estimation ####################
# i: evaluate at ith obs. (order as in the whole sample); y: outcome of interest; 
# LAM: lambda matrix around i (NOT depend on d)!!; z: additional controls
# fullkmat: index matrix, relative to the whole sample, contains results for all K's
# subset: additional subsetting, relative to the whole sample
# i.Lam: PCs around i; i.no.neigh: no. of units around i; BOTH from lpca, of length K
linkinv.d <- binomial(link="logit")$linkinv

# local regression
locreg <- function(i, y, d, x, subset=NULL, kmat, nlam, n, p, shutdown.d=F, K) {
  y.fitval <- d.fitval <- NA
  index <- kmat[,i]   # length n
  PC    <- lpca(A=x[,(p/2+1):p], index=index, nlam=nlam, n=n, K=K)
  
  # prepare design
  if (is.null(subset)) sub <- index
  else                 sub <- index & subset   # of length n
  y.sub  <- y[sub]
  #design <- PC$Lam[sub[index],,drop=F]
  design <- cbind(1, PC$Lam[sub[index],,drop=F])
  
  # eval
  #design.d <- PC$Lam
  design.d <- cbind(1, PC$Lam)
  eval <- design.d[which(PC$no.neigh==i),]
  
  # run a local GLM of y at the i-th obs.
  y.fitval <- sum(.lm.fit(design, y.sub)$coeff * eval)
  
  if (!shutdown.d) {
    # run a local GLM of d at ith obs.
    d.fitval <- linkinv.d(sum(fastLR(design.d, d[index])$coeff * eval))
  }
  
  return(c(y.fitval, d.fitval, PC$d.i))   # three scalars
}

# local regression, only for local constant approx.
locreg.cons <- function(i, y, d, subset=NULL, kmat, shutdown.d=F) {   
  y.fitval <- d.fitval <- NA  
  index  <- kmat[,i]   # length n
  if (is.null(subset)) sub <- index
  else                 sub <- index & subset   # of length n
  y.fitval <- mean(y[sub])
  if (!shutdown.d) {
    d.fitval <- mean(d[index])
  }
  
  return(c(y.fitval, d.fitval, NA))   # three scalars
}

# Prediction, for a sequence of K's
pred <- function(i, y, d, x, subset=NULL, fullkmat, nlam, n, p, shutdown.d=F, const=F, Kseq) {   
    L <- nrow(fullkmat)/n
    
    if (const) {
       out <- sapply(c(1:L), function(l) locreg.cons(i, y, d, subset, kmat=fullkmat[((l-1)*n+1):(l*n),], shutdown.d))
    } else {
       out <- sapply(c(1:L), function(l) locreg(i, y, d, x, subset, kmat=fullkmat[((l-1)*n+1):(l*n),], nlam, n, p, shutdown.d, Kseq[l]))
    }
    
    return(out)   # 3 by Kseq matrix or a vector
}


# Computation, for all units in range
compute <- function(range, y, d, x, subset=NULL, Kseq, nlam, n, p, const=F, shutdown.d=F) {
    fullkmat <- knn.index(A=x[,1:(p/2)], Kseq=Kseq)
    fit  <- sapply(range, function(i) pred(i=i, y=y, d=d, x=x, subset=subset,
                                           fullkmat=fullkmat, nlam=nlam, n=n, p=p, 
                                           shutdown.d=shutdown.d, const=const, Kseq=Kseq))
    no <- rep(1:3, length(Kseq))
    
    return(list(yfit=fit[no==1,,drop=F], ps=fit[no==2,,drop=F], d.i=fit[no==3,,drop=F]))  # each is a length(Kseq) by length(range) matrix
}

#################################################

# calculate stats related to treatment effect
te.stats <- function(y, d, pr, theta0, result, l) {
    yfit <- result$yfit[l,]; ps <- result$ps[l,]
    
    # psi <- (d * yfit + (1-d)*(y-yfit)*ps/(1-ps))/pr       # NOT influence fun, mu01
    psi <- (d * (y-yfit) - (1-d)*(y-yfit)*ps/(1-ps))/pr     # ATT
    theta <- mean(psi)    # point estimate
    
    # varphi <- (d * (yfit-theta) + (1-d)*(y-yfit)*ps/(1-ps))/pr     # IF for mu01
    varphi <- (d * (y-yfit-theta) - (1-d)*(y-yfit)*ps/(1-ps))/pr     # IF for ATT
    se    <- sqrt(mean(varphi^2))/sqrt(n)
    
    rej <- (((theta-theta0)/se > qnorm(0.975)) | ((theta-theta0)/se < qnorm(0.025))) * 1
    ci  <- se*qnorm(0.975)*2
    
    return(c(theta, se, rej, ci))
}

# apply across K's, calculate necessary quantities
fitting <- function(range, y, d, x, subset=NULL, Kseq, nlam, n, p, pr, theta0, const=F) {
  result <- compute(range=range, y=y, d=d, x=x, subset=subset, 
                    Kseq=Kseq, nlam=nlam, n=n, p=p, const=const)
    
  # calculate counterfactual mean
  out <- sapply(1:length(Kseq), function(l) te.stats(y, d, pr, theta0, result, l=l))
  
  out <- rbind(out, rowmeans(result$d.i))
  
  return(out)   # 5 by Kseq matrix
}

# sim function
sim <- function(i, n, p, model, Kseq, nlam, theta0, const, err, rho, dim.alp) {
  data   <- dgp(n=n, p=p, hdmodel=model, err=err, rho=rho, dim.alp=dim.alp)
  subset <- (data$d==0)   # so compute mu.0 in the following
  pr     <- mean(data$d)
  range  <- 1:n
  
  # 20-fold CV choice
  Klist <- seq(Kseq[1], Kseq[length(Kseq)], 5)
  if (dim.alp==1) {
    K0  <- Klist[which.min(knn.cv(nfolds=20, y=data$y, x=data$x[,p,drop=F], k=Klist, type="R")$crit)]
  } else {
    K0  <- Klist[which.min(knn.cv(nfolds=20, y=data$y, x=data$x[,c(p-1,p),drop=F], k=Klist, type="R")$crit)]
  }
  
  # fitting
  Kseq <- c(Kseq, K0)   # add the DPI choice
  result <- fitting(range=range, y=data$y, d=data$d, x=data$x, subset=subset, 
                    Kseq=Kseq, nlam=nlam, n=n, p=p, pr=pr, theta0=theta0, const=const)
  
  result <- rbind(Kseq, result)
  
  # double lasso
  lasso.fit <- rlassoATET(x=data$x, d=data$d, y=data$y, bootstrap = "none", model = FALSE,
                         penalty   = list(homoscedastic = TRUE))
  lasso.theta <- lasso.fit$te
  lasso.se    <- lasso.fit$se
  lasso.rej   <- (((lasso.theta-theta0)/lasso.se > qnorm(0.975)) | ((lasso.theta-theta0)/lasso.se < qnorm(0.025))) * 1
  lasso.ci    <- lasso.se*qnorm(0.975)*2
  
  # select on most recent two x
  dr.fit   <- dr_att(y=data$y, d=data$d, x=data$x[,c(p-1,p)])
  dr.theta <- dr.fit$te
  dr.se    <- dr.fit$se
  dr.rej   <- (((dr.theta-theta0)/dr.se > qnorm(0.975)) | ((dr.theta-theta0)/dr.se < qnorm(0.025))) * 1
  dr.ci    <- dr.se*qnorm(0.975)*2
 
  # combine results
  result <- cbind(result, c(-1, dr.theta, dr.se, dr.rej, dr.ci, NA),
                          c(-2, lasso.theta, lasso.se, lasso.rej, lasso.ci, NA))
  
  return(result)  # (4+1+1) by length(Kseq)+1+1+1 matrix, -1: select on two x, -2: double lasso
}


# simple DR estimator for ATE
dr_att <- function(y, d, x, ps_tol = 1e-6) {
  # y: numeric vector, outcome
  # d: 0/1 treatment indicator
  # x: matrix or data.frame of controls (low-dimensional)
  # ps_tol: truncation for propensity scores
  
  n <- length(y)
  
  # Design matrix with intercept
  X <- cbind(1, x)
  
  ## 1. Propensity score: fast logistic regression
  e_hat <- as.vector(1 / (1 + exp(- (X %*% fastLR(X, d)$coeff))))
  e_hat <- pmin(pmax(e_hat, ps_tol), 1 - ps_tol)  # truncate to avoid extreme weights
  
  ## 2. Outcome regression for controls: E[Y | D=0, X]
  ctrl <- (d == 0)
  m0_coef <- .lm.fit(x = X[ctrl, , drop = FALSE], y = y[ctrl])$coefficients
  m0_hat  <- as.vector(X %*% m0_coef)
  
  ## 3. DR ATT estimator
  p_hat <- mean(d)  # estimate of P(D=1)
  
  # score_i = inside of expectation in the DR representation
  score <- d / p_hat * (y - m0_hat) -
          (1 - d) / p_hat * (e_hat / (1 - e_hat)) * (y - m0_hat)
  
  tau_hat <- mean(score)
  
  ## 4. IF-based standard error
  psi_hat <- score - tau_hat
  se_hat <- sqrt(mean(psi_hat^2) / n)
  
  return(list(te = tau_hat, se  = se_hat))
}

# small simulation function for DR ATT only
sim.selonx <- function(i, n, p, model, theta0, err, rho, dim.alp, nx=1) {
  data   <- dgp(n=n, p=p, hdmodel=model, err=err, rho=rho, dim.alp=dim.alp)
  
  # double lasso
  lasso.fit <- rlassoATET(x=data$x, d=data$d, y=data$y, bootstrap = "none", model = FALSE,
                          penalty   = list(homoscedastic = TRUE))
  lasso.theta <- lasso.fit$te
  lasso.se    <- lasso.fit$se
  lasso.rej   <- (((lasso.theta-theta0)/lasso.se > qnorm(0.975)) | ((lasso.theta-theta0)/lasso.se < qnorm(0.025))) * 1
  lasso.ci    <- lasso.se*qnorm(0.975)*2
  
  # select on last nx controls
  dr.fit   <- dr_att(y=data$y, d=data$d, x=data$x[, (p-nx+1):p])
  dr.theta <- dr.fit$te
  dr.se    <- dr.fit$se
  dr.rej   <- (((dr.theta-theta0)/dr.se > qnorm(0.975)) | ((dr.theta-theta0)/dr.se < qnorm(0.025))) * 1
  dr.ci    <- dr.se*qnorm(0.975)*2
  
  # combine results
  result <- cbind(c(-1, dr.theta, dr.se, dr.rej, dr.ci),
                  c(-2, lasso.theta, lasso.se, lasso.rej, lasso.ci))
  
  return(result)  # (4+1) by 1+1 matrix, -1: select on nx controls, -2: double lasso
}
