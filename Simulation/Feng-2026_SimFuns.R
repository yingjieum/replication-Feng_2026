# Supporting functions for the simulation study in Feng (2026).

# Potential-outcome functions and treatment probabilities.
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

# One-dimensional factor models.
# Model 1 is exactly rank 3.
hdfun1 <- function(u,v) {
  s   <- outer(u, v, "*")
  out <- -4.5*s + 18.4*s^2 - 16.2*s^3
  return(out)
}

# Model 2 is approximately rank 3.
hdfun2 <- function(u, v, sig=.1) {
  s <- outer(u, v, "-")
  out <- 2*exp(-s^2/sig)
  return(out)
}

# Main-paper Gaussian-kernel designs with smoother and rougher signals.
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

# Model 3 is approximately rank 5 or higher.
hdfun3 <- function(u, v, sig=.1) {
  s <- outer(u, v, "-")
  out <- 2*exp(-abs(s)/sig)
  return(out)
}

# Two-dimensional factor models.
# Model 4 is exactly rank 9 and approximately rank 4.
hdfun4 <- function(u1, u2, v1, v2) {
  s <- outer(2*(u1-0.5), 2*(v1-0.5), "*") + outer(2*(u2-0.5), 2*(v2-0.5), "*")
  out <- 2*s + (s*s) + 0.5*(s*s*s)
  return(out)
}

# Model 5 is approximately rank 5 or higher.
hdfun5 <- function(u1, u2, v1, v2, sig=.1) {
  d1 <- outer(u1, v1, "-")
  d2 <- outer(u2, v2, "-")
  out <- 2*exp(-(d1*d1+d2*d2)/sig)
  return(out)
}

# Model 6 is approximately rank 3 or higher.
hdfun6 <- function(u1, u2, v1, v2, sig=.5) {
  d1 <- outer(u1, v1, "-")
  d2 <- outer(u2, v2, "-")
  out <- 2*exp(-(abs(d1)+abs(d2))/sig)
  return(out)
}


funlist <- list(hdfun1=hdfun1, hdfun2=hdfun2, hdfun3=hdfun3, hdfun4=hdfun4, hdfun5=hdfun5, hdfun6=hdfun6, hdfun7=hdfun7, hdfun8=hdfun8)

# Generate stationary AR(1) measurement errors.
gen_err <- function(n, p, sigma=1, rho) {
  U <- matrix(0, n, p)
  
  # Start in the stationary distribution N(0, sigma^2).
  U[, 1] <- rnorm(n, sd = sigma)
  
  sd_eps <- sigma * sqrt(1 - rho^2)
  Eps <- matrix(rnorm(n * (p - 1), sd = sd_eps), n, p - 1)
  
  for (t in 2:p) {
    U[, t] <- rho * U[, t - 1] + Eps[, t - 1]
  }
  return(U)
}


# Generate one simulated dataset.
# err: 0 for i.i.d., 1 for heteroskedastic, and 2 for AR(1) measurement errors.
dgp <- function(n, p, hdmodel, sigma=1, err=0, rho=0.2, dim.alp=1) {
    # latent variable
    alpha1 <- runif(n, 0, 1)
    alpha2 <- NULL
    if (dim.alp==2) {
      alpha2 <- runif(n, 0, 1)
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
    
    return(list(y=y, y0=Y0, y1=Y1, x=x, alpha1=alpha1, alpha2=alpha2, eta1=eta1, eta2=eta2, d=D))
}

# Find nearest-neighbor sets for each value of K.
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

# Select the local factor dimension by singular-value thresholding.
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

# Estimate local factors for one evaluation unit.
# A is n by p; index selects its neighborhood; nlam is the maximum extracted dimension.
lpca <- function(A, index, nlam, n, K) {      
    A <- cbind(1:n, A)   # add a row number
    A <- A[index,]
    svd <- irlba(A[,-1], nu=nlam, nv=nlam)
    r <- sel.r(svd$d, K, ncol(A)-1)
    return(list(Lam=svd$u[,1:r,drop=F], no.neigh=A[,1], sv=svd$d, d.i=r))    #Lam: K by r
}


# Local estimation and ATT inference.
# subset optionally restricts the outcome-regression sample.
linkinv.d <- binomial(link="logit")$linkinv

# Local outcome and propensity-score regressions.
locreg <- function(i, y, d, x, subset=NULL, kmat, nlam, n, p, K) {
  index <- kmat[,i]   # length n
  PC    <- lpca(A=x[,(p/2+1):p], index=index, nlam=nlam, n=n, K=K)
  
  # Local-factor regression design.
  if (is.null(subset)) sub <- index
  else                 sub <- index & subset   # of length n
  y.sub  <- y[sub]
  design <- cbind(1, PC$Lam[sub[index],,drop=F])
  
  # Regressor values for the evaluation unit.
  design.d <- cbind(1, PC$Lam)
  eval <- design.d[which(PC$no.neigh==i),]
  
  # Outcome regression at the evaluation unit.
  y.fitval <- sum(.lm.fit(design, y.sub)$coeff * eval)
  
  # Propensity score at the evaluation unit.
  d.fitval <- linkinv.d(sum(fastLR(design.d, d[index])$coeff * eval))
  
  return(c(y.fitval, d.fitval, PC$d.i))   # three scalars
}

# Local-constant outcome and propensity-score estimates.
locreg.cons <- function(i, y, d, subset=NULL, kmat) {
  index  <- kmat[,i]   # length n
  if (is.null(subset)) sub <- index
  else                 sub <- index & subset   # of length n
  y.fitval <- mean(y[sub])
  d.fitval <- mean(d[index])
  
  return(c(y.fitval, d.fitval, NA))   # three scalars
}

# Predict for a sequence of neighborhood sizes.
pred <- function(i, y, d, x, subset=NULL, fullkmat, nlam, n, p, const=F, Kseq) {
    L <- nrow(fullkmat)/n
    
    if (const) {
       out <- sapply(c(1:L), function(l) locreg.cons(i, y, d, subset, kmat=fullkmat[((l-1)*n+1):(l*n),]))
    } else {
       out <- sapply(c(1:L), function(l) locreg(i, y, d, x, subset, kmat=fullkmat[((l-1)*n+1):(l*n),], nlam=nlam, n=n, p=p, K=Kseq[l]))
    }
    
    return(out)   # 3 by Kseq matrix or a vector
}


# Compute predictions for all requested units.
compute <- function(range, y, d, x, subset=NULL, Kseq, nlam, n, p, const=F) {
    fullkmat <- knn.index(A=x[,1:(p/2)], Kseq=Kseq)
    fit  <- sapply(range, function(i) pred(i=i, y=y, d=d, x=x, subset=subset,
                                           fullkmat=fullkmat, nlam=nlam, n=n, p=p, 
                                           const=const, Kseq=Kseq))
    no <- rep(1:3, length(Kseq))
    
    return(list(yfit=fit[no==1,,drop=F], ps=fit[no==2,,drop=F], d.i=fit[no==3,,drop=F]))  # each is a length(Kseq) by length(range) matrix
}

# Compute the ATT estimate and inference statistics.
te.stats <- function(y, d, pr, theta0, result, l) {
    yfit <- result$yfit[l,]; ps <- result$ps[l,]
    
    psi <- (d * (y-yfit) - (1-d)*(y-yfit)*ps/(1-ps))/pr     # ATT
    theta <- mean(psi)    # point estimate
    
    varphi <- (d * (y-yfit-theta) - (1-d)*(y-yfit)*ps/(1-ps))/pr     # IF for ATT
    se    <- sqrt(mean(varphi^2))/sqrt(n)
    
    rej <- (((theta-theta0)/se > qnorm(0.975)) | ((theta-theta0)/se < qnorm(0.025))) * 1
    ci  <- se*qnorm(0.975)*2
    
    return(c(theta, se, rej, ci))
}

# Apply ATT estimation across neighborhood sizes.
fitting <- function(range, y, d, x, subset=NULL, Kseq, nlam, n, p, pr, theta0, const=F) {
  result <- compute(range=range, y=y, d=d, x=x, subset=subset, 
                    Kseq=Kseq, nlam=nlam, n=n, p=p, const=const)
    
  # Estimate the untreated conditional mean.
  out <- sapply(1:length(Kseq), function(l) te.stats(y, d, pr, theta0, result, l=l))
  
  out <- rbind(out, rowmeans(result$d.i))
  
  return(out)   # 5 by Kseq matrix
}

# Run one replication for a table design.
sim <- function(n, p, model, Kseq, nlam, theta0, const, err, rho, dim.alp) {
  data   <- dgp(n=n, p=p, hdmodel=model, err=err, rho=rho, dim.alp=dim.alp)
  subset <- (data$d==0)   # so compute mu.0 in the following
  pr     <- mean(data$d)
  range  <- 1:n
  
  # Select K by 20-fold cross-validation.
  Klist <- seq(Kseq[1], Kseq[length(Kseq)], 5)
  if (dim.alp==1) {
    K0  <- Klist[which.min(knn.cv(nfolds=20, y=data$y, x=data$x[,p,drop=F], k=Klist, type="R")$crit)]
  } else {
    K0  <- Klist[which.min(knn.cv(nfolds=20, y=data$y, x=data$x[,c(p-1,p),drop=F], k=Klist, type="R")$crit)]
  }
  
  # Proposed local-PCA estimates.
  Kseq <- c(Kseq, K0)   # add the cross-validated choice
  result <- fitting(range=range, y=data$y, d=data$d, x=data$x, subset=subset, 
                    Kseq=Kseq, nlam=nlam, n=n, p=p, pr=pr, theta0=theta0, const=const)
  
  result <- rbind(Kseq, result)
  
  # Double-Lasso benchmark.
  lasso.fit <- rlassoATET(x=data$x, d=data$d, y=data$y, bootstrap = "none", model = FALSE,
                         penalty   = list(homoscedastic = TRUE))
  lasso.theta <- lasso.fit$te
  lasso.se    <- lasso.fit$se
  lasso.rej   <- (((lasso.theta-theta0)/lasso.se > qnorm(0.975)) | ((lasso.theta-theta0)/lasso.se < qnorm(0.025))) * 1
  lasso.ci    <- lasso.se*qnorm(0.975)*2
  
  # Doubly robust benchmark using the last two proxies.
  dr.fit   <- dr_att(y=data$y, d=data$d, x=data$x[,c(p-1,p)])
  dr.theta <- dr.fit$te
  dr.se    <- dr.fit$se
  dr.rej   <- (((dr.theta-theta0)/dr.se > qnorm(0.975)) | ((dr.theta-theta0)/dr.se < qnorm(0.025))) * 1
  dr.ci    <- dr.se*qnorm(0.975)*2
 
  # Assemble results in the order used by the output tables.
  result <- cbind(result, c(-1, dr.theta, dr.se, dr.rej, dr.ci, NA),
                          c(-2, lasso.theta, lasso.se, lasso.rej, lasso.ci, NA))
  
  return(result)  # (4+1+1) by length(Kseq)+1+1+1 matrix, -1: select on two x, -2: double lasso
}


# Doubly robust ATT estimator with low-dimensional controls.
dr_att <- function(y, d, x, ps_tol = 1e-6) {
  n <- length(y)
  
  # Add an intercept to both nuisance regressions.
  X <- cbind(1, x)
  
  # Fit and truncate propensity scores.
  e_hat <- as.vector(1 / (1 + exp(- (X %*% fastLR(X, d)$coeff))))
  e_hat <- pmin(pmax(e_hat, ps_tol), 1 - ps_tol)  # truncate to avoid extreme weights
  
  # Fit E[Y | D=0, X] on the control observations.
  ctrl <- (d == 0)
  m0_coef <- .lm.fit(x = X[ctrl, , drop = FALSE], y = y[ctrl])$coefficients
  m0_hat  <- as.vector(X %*% m0_coef)
  
  # Compute the doubly robust ATT score.
  p_hat <- mean(d)  # estimate of P(D=1)
  
  score <- d / p_hat * (y - m0_hat) -
          (1 - d) / p_hat * (e_hat / (1 - e_hat)) * (y - m0_hat)
  
  tau_hat <- mean(score)
  
  # Center the influence function for the estimated treated share.
  psi_hat <- score - d * tau_hat / p_hat
  se_hat <- sqrt(mean(psi_hat^2) / n)
  
  return(list(te = tau_hat, se  = se_hat))
}

# Run one replication for the three estimators in Figure 1.
sim.selonx <- function(n, p, model, theta0, err, rho, dim.alp, nx=1) {
  data   <- dgp(n=n, p=p, hdmodel=model, err=err, rho=rho, dim.alp=dim.alp)
  
  # Double-Lasso benchmark.
  lasso.fit <- rlassoATET(x=data$x, d=data$d, y=data$y, bootstrap = "none", model = FALSE,
                          penalty   = list(homoscedastic = TRUE))
  lasso.theta <- lasso.fit$te
  lasso.se    <- lasso.fit$se
  lasso.rej   <- (((lasso.theta-theta0)/lasso.se > qnorm(0.975)) | ((lasso.theta-theta0)/lasso.se < qnorm(0.025))) * 1
  lasso.ci    <- lasso.se*qnorm(0.975)*2
  
  # Doubly robust benchmark using the last nx proxies.
  dr.fit   <- dr_att(y=data$y, d=data$d, x=data$x[, (p-nx+1):p])
  dr.theta <- dr.fit$te
  dr.se    <- dr.fit$se
  dr.rej   <- (((dr.theta-theta0)/dr.se > qnorm(0.975)) | ((dr.theta-theta0)/dr.se < qnorm(0.025))) * 1
  dr.ci    <- dr.se*qnorm(0.975)*2
  
  # Proposed estimator with the same tuning rule as sim().
  Kbase <- n^(4/(4+dim.alp))
  Klist <- seq(ceiling(0.5*Kbase), ceiling(2*Kbase), 5)
  x.cv <- if (dim.alp==1) data$x[,p,drop=F] else data$x[,c(p-1,p),drop=F]
  K0 <- Klist[which.min(knn.cv(nfolds=20, y=data$y, x=x.cv, k=Klist, type="R")$crit)]
  lpca.fit <- compute(range=1:n, y=data$y, d=data$d, x=data$x, subset=(data$d==0),
                      Kseq=K0, nlam=4, n=n, p=p, const=F)
  pr <- mean(data$d)
  lpca.score <- (data$d*(data$y-lpca.fit$yfit[1,]) -
                  (1-data$d)*(data$y-lpca.fit$yfit[1,])*lpca.fit$ps[1,]/(1-lpca.fit$ps[1,]))/pr
  lpca.theta <- mean(lpca.score)
  lpca.varphi <- lpca.score - data$d*lpca.theta/pr   # same ATT influence function as te.stats()
  lpca.se <- sqrt(mean(lpca.varphi^2)/n)
  lpca.rej <- as.numeric(abs((lpca.theta-theta0)/lpca.se) > qnorm(0.975))
  lpca.ci <- lpca.se*qnorm(0.975)*2

  # Assemble results in the order used by the output tables.
  result <- cbind(c(-1, dr.theta, dr.se, dr.rej, dr.ci),
                  c(-2, lasso.theta, lasso.se, lasso.rej, lasso.ci),
                  c(K0, lpca.theta, lpca.se, lpca.rej, lpca.ci))
  
  return(result)  # 5 by 3 matrix; first row: -1 (nx controls), -2 (double lasso), selected K (proposed)
}
