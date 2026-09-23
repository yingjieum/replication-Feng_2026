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

############################################################
# Additional empirical diagnostics, separate from the benchmark helpers.

# Run a stochastic check without changing the main analysis's RNG state.
check_seed <- function(seed, expr) {
  had.seed <- exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
  if (had.seed) old.seed <- get(".Random.seed", envir=.GlobalEnv)
  on.exit(if (had.seed) assign(".Random.seed", old.seed, envir=.GlobalEnv) else
    if (exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE))
      rm(".Random.seed", envir=.GlobalEnv))
  set.seed(seed)
  force(expr)
}

# Record numerical warnings locally rather than suppressing their evidence.
check_capture <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(tryCatch(expr, error=function(e) e),
    warning=function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(value=value, warnings=unique(warnings))
}

# Save an auditable table using the check- prefix.
check_write <- function(state, name, tab) {
  state$tables[[name]] <- tab
  write.csv(tab, file.path(state$output, paste0("check-", name, ".csv")),
            row.names=FALSE, na="NA")
  invisible(tab)
}

# Define the exact paper samples, unused outcomes, and blocked proxy splits.
check_setup <- function(data, x, y, d, z, K,
                        output=Sys.getenv("FENG_CHECK_OUTPUT", "output/intermediate/application-diagnostics")) {
  state <- new.env(parent=emptyenv())
  state$output <- output
  dir.create(output, showWarnings=FALSE, recursive=TRUE)
  state$cache <- new.env(parent=emptyenv())
  state$tables <- list()
  state$K <- K
  p <- ncol(x)
  stopifnot(p==250L, all(is.finite(x)), all(is.finite(y)),
            all(is.finite(z)), all(d %in% 0:1))
  base <- as.vector(data$CorrCiti <= as.numeric(data$CorrCitiTr))
  stopifnot(!anyNA(base))
  day <- -30:-1
  pre <- t(data$Re[as.integer(data$GeiNomDat)-1+day,,drop=FALSE])
  colnames(pre) <- paste0("day", day)
  windows <- sapply(list(-30:-21, -20:-11, -10:-1, -30:-1), function(days)
    rowSums(pre[,match(days, day),drop=FALSE]))
  colnames(windows) <- c("CAR[-30,-21]", "CAR[-20,-11]", "CAR[-10,-1]", "CAR[-30,-1]")
  colnames(z) <- c("Log assets", "ROE", "Leverage")
  outcomes <- cbind(event=y, pre, windows)
  state$samples <- lapply(list(full=seq_along(d), base=which(base)), function(ids)
    list(ids=ids, x=x[ids,,drop=FALSE], y=outcomes[ids,,drop=FALSE],
         d=d[ids], z=z[ids,,drop=FALSE], balance=cbind(z, pre, windows)[ids,,drop=FALSE]))
  # Each group is contiguous and ordered; swapping groups never reverses dates.
  groups <- list(benchmark=1:125, reversed=126:250,
    boundary100=1:100, boundary150=1:150)
  state$splits <- lapply(groups, function(m) list(match=m, pca=setdiff(1:p,m)))
  stopifnot(all(vapply(state$splits,function(g)
    all(diff(g$match)==1) && all(diff(g$pca)==1) &&
    !length(intersect(g$match,g$pca)),logical(1))))
  state$folds <- lapply(state$samples, function(s)
    check_seed(20260917, sample(rep(1:5, length.out=length(s$d)))))
  check_write(state, "samples", do.call(rbind, lapply(names(state$samples), function(nm) {
    s <- state$samples[[nm]]
    data.frame(sample=nm, n=length(s$d), treated=sum(s$d), controls=sum(s$d==0),
      proxies=p, K=K, outcome="CAR[0,1]", proxy_first_day=-280, proxy_last_day=-31,
      nomination_date=as.character(as.Date(as.numeric(data$Statadate[as.integer(data$GeiNomDat)]), origin="1960-01-01")))
  })))
  check_write(state, "proxy-splits", do.call(rbind, lapply(names(groups), function(nm)
    data.frame(split=nm, proxy=1:p, event_day=1:p-281,
      group=ifelse(1:p %in% groups[[nm]], "matching", "PCA")))))
  state
}

# Compute center-inclusive neighborhoods under the three implemented distances.
check_neighbors <- function(state, s, tag, split="benchmark", K=state$K,
                            distance="pseudo-max") {
  key <- paste("neighbors", tag, split, K, distance, sep="|")
  if (exists(key, state$cache, inherits=FALSE)) return(get(key, state$cache))
  A <- s$x[,state$splits[[split]]$match,drop=FALSE]
  n <- nrow(A)
  if (distance=="pseudo-max") {
    dis <- dist.knn(A)
  } else if (distance=="squared-Euclidean") {
    dis <- pmax(outer(rowSums(A^2), rowSums(A^2), "+") - 2*tcrossprod(A), 0)/ncol(A)
  } else if (distance=="averages") {
    dis <- abs(outer(rowMeans(A), rowMeans(A), "-"))
  } else stop("Unknown diagnostic distance")
  diag(dis) <- 0
  kmat <- apply(dis, 2, function(v) v <= sort(v, partial=K)[K])
  stopifnot(all(diag(kmat)), K<n, all(colSums(kmat)>=K))
  tab <- do.call(rbind, lapply(seq_len(n), function(i) {
    nb <- which(kmat[,i]); other <- setdiff(nb,i)
    maximum <- max(dis[other,i]); scale <- sd(dis[-i,i])
    data.frame(sample=tag, split=split, distance=distance, K=K, unit=s$ids[i],
      treated=s$d[i], neighbors=length(nb), local_treated=sum(s$d[nb]),
      local_controls=sum(s$d[nb]==0), max_discrepancy=maximum,
      normalized_max=maximum/scale)
  }))
  out <- list(kmat=kmat, dis=dis, tab=tab)
  assign(key, out, state$cache)
  out
}

# Step 1.1: describe matches and split stability for every center in both samples.
check_matching <- function(state) {
  tabs <- stability <- list()
  for (nm in names(state$samples)) {
    s <- state$samples[[nm]]
    ref <- check_neighbors(state, s, nm)
    for (split in names(state$splits)) {
      cur <- check_neighbors(state, s, nm, split)
      tabs[[paste(nm,split)]] <- cur$tab
      if (split!="benchmark") stability[[paste(nm,split)]] <-
        do.call(rbind, lapply(seq_along(s$d), function(i) {
          a <- which(ref$kmat[,i]); a <- setdiff(a,i)
          b <- which(cur$kmat[,i]); b <- setdiff(b,i)
          data.frame(sample=nm, split=split, unit=s$ids[i], treated=s$d[i],
            overlap_fraction=length(intersect(a,b))/length(a),
            jaccard=length(intersect(a,b))/length(union(a,b)),
            chance_overlap=length(b)/(length(s$d)-1))
        }))
    }
  }
  check_write(state, "neighborhoods", do.call(rbind,tabs))
  check_write(state, "split-stability", do.call(rbind,stability))
  invisible(state)
}

# Cache exact local SVDs, including the full spectrum to diagnose the four-PC cap.
check_basis <- function(state, s, tag, split="benchmark", K=state$K,
                        distance="pseudo-max") {
  key <- paste("basis", tag, split, K, distance, sep="|")
  if (exists(key, state$cache, inherits=FALSE)) return(get(key, state$cache))
  neighbors <- check_neighbors(state,s,tag,split,K,distance)
  cols <- state$splits[[split]]$pca
  out <- lapply(seq_along(s$d), function(i) {
    ids <- which(neighbors$kmat[,i]); A <- s$x[ids,cols,drop=FALSE]
    fit <- svd(A, nu=min(4,nrow(A),ncol(A)), nv=0)
    list(ids=ids, u=fit$u, sv=fit$d,
      threshold=log(log(K))*(sqrt(K)+sqrt(ncol(A)))*sd(A))
  })
  assign(key,out,state$cache)
  out
}

# Apply the existing forced-at-least-one rule or an explicitly named sensitivity.
check_dimension <- function(pc, multiplier=1, cap=4) {
  min(max(sum(pc$sv > multiplier*pc$threshold),1),cap,ncol(pc$u))
}

# Step 1.2: report all spectra, fixed-neighborhood subspace stability, and serial dependence.
check_factors <- function(state) {
  signals <- dimensions <- spaces <- serial <- list()
  for (nm in names(state$samples)) {
    s <- state$samples[[nm]]; pcs <- check_basis(state,s,nm)
    signals[[nm]] <- do.call(rbind,lapply(seq_along(pcs),function(i) {
      pc <- pcs[[i]]
      data.frame(sample=nm,unit=s$ids[i],treated=s$d[i],component=seq_along(pc$sv),
        singular_value=pc$sv,threshold=pc$threshold,
        singular_threshold_ratio=pc$sv/pc$threshold,
        normalized_eigenvalue=pc$sv^2/(length(pc$ids)*125),
        eigen_threshold_ratio=(pc$sv/pc$threshold)^2)
    }))
    dimensions[[nm]] <- do.call(rbind,lapply(seq_along(pcs),function(i) {
      pc <- pcs[[i]]; above <- sum(pc$sv>pc$threshold)
      data.frame(sample=nm,unit=s$ids[i],treated=s$d[i],selected_dimension=check_dimension(pc),
        above_threshold=above,forced_one=above==0,cap_binds=above>4,
        leading_ratio=pc$sv[1]/pc$threshold)
    }))
    for (split in setdiff(names(state$splits),"benchmark")) {
      # Keep benchmark neighborhoods fixed; only the PCA columns change.
      spaces[[paste(nm,split)]] <- do.call(rbind,lapply(seq_along(pcs),function(i) {
        a <- pcs[[i]]; A <- s$x[a$ids,state$splits[[split]]$pca,drop=FALSE]
        fit <- svd(A,nu=4,nv=0)
        b <- list(u=fit$u,sv=fit$d,threshold=log(log(state$K))*(sqrt(state$K)+sqrt(ncol(A)))*sd(A))
        ra <- check_dimension(a); rb <- check_dimension(b)
        Ua <- a$u[,1:ra,drop=FALSE]; Ub <- b$u[,1:rb,drop=FALSE]
        common <- min(ra,rb)
        cosine <- svd(crossprod(Ua,Ub),nu=0,nv=0)$d
        data.frame(sample=nm,split=split,unit=s$ids[i],treated=s$d[i],rank_reference=ra,rank_alternative=rb,
          projection_distance=sqrt(max(0,ra+rb-2*sum(crossprod(Ua,Ub)^2))/(ra+rb)),
          max_common_angle_degrees=acos(min(1,max(0,min(cosine))))*180/pi)
      }))
    }
    for (lag in c(1,2,5,10,25)) {
      ac <- apply(s$x,1,function(v) cor(v[1:(length(v)-lag)],v[(lag+1):length(v)]))
      serial[[paste(nm,lag)]] <- data.frame(sample=nm,series="returns",lag=lag,
        median_correlation=median(ac),median_absolute_correlation=median(abs(ac)),
        q90_absolute_correlation=as.numeric(quantile(abs(ac),.9)))
    }
  }
  check_write(state,"signals",do.call(rbind,signals))
  check_write(state,"dimensions",do.call(rbind,dimensions))
  check_write(state,"factor-space",do.call(rbind,spaces))
  check_write(state,"serial-dependence",do.call(rbind,serial))
  invisible(state)
}

# Fit local outcomes and logits once per center; retain original untrimmed predictions.
check_fit <- function(state,s,tag,split="benchmark",K=state$K,distance="pseudo-max",
                      multiplier=1,cap=4,covariates=FALSE,leave_center_out=FALSE) {
  pcs <- check_basis(state,s,tag,split,K,distance)
  predictions <- matrix(NA_real_,length(s$d),ncol(s$y),dimnames=list(NULL,colnames(s$y)))
  ps <- numeric(length(s$d)); records <- vector("list",length(s$d))
  for (i in seq_along(s$d)) {
    pc <- pcs[[i]]; nb <- pc$ids; r <- check_dimension(pc,multiplier,cap)
    design <- cbind(Intercept=1,pc$u[,1:r,drop=FALSE])
    if (covariates) design <- cbind(design,s$z[nb,,drop=FALSE])
    pos <- match(i,nb); train <- rep(TRUE,length(nb))
    if (leave_center_out) train[pos] <- FALSE
    controls <- train & s$d[nb]==0
    out <- check_capture(.lm.fit(design[controls,,drop=FALSE],s$y[nb[controls],,drop=FALSE]))
    logit <- check_capture(glm.fit(x=design[train,,drop=FALSE],y=s$d[nb[train]],family=binomial()))
    if (!inherits(out$value,"error")) predictions[i,] <- drop(design[pos,,drop=FALSE] %*% out$value$coefficients)
    if (!inherits(logit$value,"error")) {
      ps[i] <- logit$value$family$linkinv(sum(logit$value$coefficients*design[pos,]))
    } else ps[i] <- NA_real_
    scaled <- sweep(design,2,sqrt(colSums(design^2)),"/")
    records[[i]] <- data.frame(sample=tag,unit=s$ids[i],treated=s$d[i],K=K,split=split,
      covariates=covariates,leave_center_out=leave_center_out,dimension=r,
      local_treated=sum(s$d[nb[train]]),local_controls=sum(controls),parameters=ncol(design),
      outcome_rank=if (inherits(out$value,"error")) NA_integer_ else out$value$rank,
      logit_rank=if (inherits(logit$value,"error")) NA_integer_ else logit$value$rank,
      converged=if (inherits(logit$value,"error")) FALSE else isTRUE(logit$value$converged),
      upper_fitted_probability=if (inherits(logit$value,"error")) NA else
        any(logit$value$fitted.values>1-1e-8),
      design_condition=kappa(scaled,exact=TRUE),ps=ps[i],
      outcome_warning=paste(out$warnings,collapse="; "),
      logit_warning=paste(logit$warnings,collapse="; "),
      error=paste(c(if(inherits(out$value,"error"))conditionMessage(out$value),
                    if(inherits(logit$value,"error"))conditionMessage(logit$value)),collapse="; "))
  }
  list(yfit=predictions,ps=ps,records=do.call(rbind,records))
}

# Compute the paper's ATT and influence-function SE without clipping or deleting failures.
check_att <- function(y,d,fit,ps) {
  if (is.null(dim(y))) y <- matrix(y,ncol=1)
  if (is.null(dim(fit))) fit <- matrix(fit,ncol=1)
  pr <- mean(d); weights <- ifelse(d==1,0,ps/(1-ps))
  score <- sweep(y-fit,1,d/pr-weights/pr,"*")
  tau <- colMeans(score)
  influence <- score - outer(d/pr,tau)
  se <- sqrt(colMeans(influence^2)/length(d))
  list(att=tau,se=se,influence=influence)
}

# Step 2: validate prediction on unseen units AND unseen proxy dates.
# Predictions use held-out unit folds; covariate-only comparators never use proxies.
check_prediction_design <- function(state,s,nm,design,matching,pca,target,controls_only=FALSE) {
  fold <- state$folds[[nm]]
  kmat <- knn.index(s$x[,matching,drop=FALSE],K=state$K)
  cov.linear <- cov.knn <- matrix(NA_real_,nrow(target),ncol(target))
  zscaled <- vector("list",5)
  for (f in 1:5) {
    train <- which(fold!=f & (!controls_only | s$d==0))
    test <- which(fold==f)
    center <- colMeans(s$z[train,,drop=FALSE])
    scale <- apply(s$z[train,,drop=FALSE],2,sd)
    stopifnot(all(scale>0),length(train)>=state$K)
    Z <- sweep(sweep(s$z,2,center,"-"),2,scale,"/")
    zscaled[[f]] <- Z
    model <- .lm.fit(cbind(1,Z[train,,drop=FALSE]),target[train,,drop=FALSE])
    cov.linear[test,] <- cbind(1,Z[test,,drop=FALSE]) %*% model$coefficients
    for (i in test) {
      distance <- rowSums(sweep(Z[train,,drop=FALSE],2,Z[i,],"-")^2)
      nb <- train[order(distance)[seq_len(state$K)]]
      cov.knn[i,] <- colMeans(target[nb,,drop=FALSE])
    }
  }
  do.call(rbind,lapply(seq_along(s$d),function(i) {
    train <- which(kmat[,i] & fold!=fold[i])
    stopifnot(!i %in% train,length(train)>8)
    A <- s$x[train,pca,drop=FALSE]
    fit <- svd(A,nu=4,nv=4)
    pc <- list(u=fit$u,sv=fit$d,threshold=log(log(length(train)))*
      (sqrt(length(train))+sqrt(ncol(A)))*sd(A))
    r <- check_dimension(pc)
    score <- drop(s$x[i,pca,drop=FALSE] %*% fit$v[,1:r,drop=FALSE])/fit$d[1:r]
    reg <- if (controls_only) s$d[train]==0 else rep(TRUE,length(train))
    Z <- zscaled[[fold[i]]]
    design.local <- cbind(1,fit$u[reg,1:r,drop=FALSE])
    ytrain <- target[train[reg],,drop=FALSE]
    model <- .lm.fit(design.local,ytrain)
    pred <- drop(c(1,score) %*% model$coefficients)
    model.z <- .lm.fit(cbind(design.local,Z[train[reg],,drop=FALSE]),ytrain)
    pred.z <- drop(c(1,score,Z[i,]) %*% model.z$coefficients)
    actual <- target[i,]
    data.frame(sample=nm,design=design,unit=s$ids[i],treated=s$d[i],fold=fold[i],
      dimension=r,training_units=length(train),regression_units=sum(reg),heldout_proxies=ncol(target),
      factor_mse=mean((actual-pred)^2),factor_controls_mse=mean((actual-pred.z)^2),
      mean_mse=mean((actual-colMeans(ytrain))^2),
      covariates_linear_mse=mean((actual-cov.linear[i,])^2),
      covariates_knn_mse=mean((actual-cov.knn[i,])^2))
  }))
}

check_heldout <- function(state) {
  # All proxy dates enter three disjoint, contiguous groups; the ATT estimator is unchanged.
  designs <- list(chronological_gapless=list(match=1:100,pca=101:200,hold=201:250),
    reversed_gapless=list(match=151:250,pca=51:150,hold=1:50))
  tabs <- definitions <- outcome <- list()
  for (nm in names(state$samples)) {
    s <- state$samples[[nm]]
    for (name in names(designs)) {
      cols <- designs[[name]]
      stopifnot(all(vapply(cols,function(v) all(diff(v)==1),logical(1))),
        !length(intersect(cols$match,cols$pca)),!length(intersect(cols$match,cols$hold)),
        !length(intersect(cols$pca,cols$hold)),
        length(cols$match)==100L,length(cols$pca)==100L,length(cols$hold)==50L,
        identical(sort(unlist(cols,use.names=FALSE)),seq_len(ncol(s$x))))
      definitions[[paste(nm,name)]] <- do.call(rbind,lapply(names(cols),function(g)
        data.frame(sample=nm,design=name,group=g,proxy=cols[[g]],event_day=cols[[g]]-281)))
      message("Held-out prediction: ",nm," / ",name)
      tabs[[paste(nm,name)]] <- check_prediction_design(state,s,nm,name,cols$match,cols$pca,
        s$x[,cols$hold,drop=FALSE])
    }
    message("Untreated-outcome validation: ",nm," / K=",state$K)
    v <- check_prediction_design(state,s,nm,"untreated_event",state$splits$benchmark$match,
      state$splits$benchmark$pca,s$y[,1,drop=FALSE],controls_only=TRUE)
    outcome[[nm]] <- v[v$treated==0,]
  }
  check_write(state,"heldout-designs",do.call(rbind,definitions))
  check_write(state,"heldout-prediction",do.call(rbind,tabs))
  check_write(state,"outcome-prediction",do.call(rbind,outcome))
  invisible(state)
}

# ATT balance uses the treated sample SD as a fixed denominator for every scheme.
check_balance <- function(s,weights,scheme,sample) {
  treated <- s$d==1; control <- !treated
  reference <- apply(s$balance[treated,,drop=FALSE],2,sd)
  mt <- colMeans(s$balance[treated,,drop=FALSE])
  mc <- colSums(sweep(s$balance[control,,drop=FALSE],1,weights[control],"*"))/sum(weights[control])
  data.frame(sample=sample,scheme=scheme,variable=colnames(s$balance),
    treated_mean=mt,control_mean=mc,treated_sd=reference,smd=(mt-mc)/reference)
}

# Step 2: nuisance-fit warnings, overlap, ATT weights, and balance in unused returns/controls.
check_nuisance <- function(state) {
  records <- balance <- overlap <- list()
  state$baseline <- list()
  for (nm in names(state$samples)) {
    s <- state$samples[[nm]]
    zfit <- check_capture(glm.fit(cbind(1,s$z),s$d,family=binomial()))
    schemes <- list(Unadjusted=rep(1,length(s$d)),
      Controls.only=ifelse(s$d==1,1,zfit$value$fitted.values/(1-zfit$value$fitted.values)))
    for (cov in c(FALSE,TRUE)) {
      key <- paste(nm,cov,sep="|")
      fit <- check_fit(state,s,nm,covariates=cov)
      state$baseline[[key]] <- fit
      records[[key]] <- fit$records
      scheme <- if(cov) "Local.PCA.with.controls" else "Local.PCA"
      w <- ifelse(s$d==1,1,fit$ps/(1-fit$ps))
      schemes[[scheme]] <- w
      controls <- s$d==0; treated <- !controls
      lo <- min(fit$ps[controls]); hi <- max(fit$ps[controls])
      overlap[[key]] <- data.frame(sample=nm,covariates=cov,
        treated_ps_min=min(fit$ps[treated]),treated_ps_median=median(fit$ps[treated]),
        treated_ps_max=max(fit$ps[treated]),control_ps_min=lo,control_ps_median=median(fit$ps[controls]),control_ps_max=hi,
        treated_outside_control_range=sum(fit$ps[treated]<lo | fit$ps[treated]>hi),
        control_ESS=sum(w[controls])^2/sum(w[controls]^2),max_control_weight=max(w[controls]),
        control_weight_sum=sum(w[controls]),treated_count=sum(treated),
        top5_control_weight_share=sum(sort(w[controls],decreasing=TRUE)[1:5])/sum(w[controls]),
        near_one_control_score=sum(fit$ps[controls]>1-1e-8),
        control_ps_above_09=sum(fit$ps[controls]>.9),
        nonconverged=sum(!fit$records$converged),
        nonconverged_treated=sum(!fit$records$converged[treated]),
        logit_warning_count=sum(nzchar(fit$records$logit_warning)),nonfinite_predictions=sum(!is.finite(fit$yfit)),
        nonfinite_scores=sum(!is.finite(fit$ps)))
    }
    for (scheme in names(schemes)) balance[[paste(nm,scheme)]] <- check_balance(s,schemes[[scheme]],scheme,nm)
  }
  check_write(state,"local-fits",do.call(rbind,records))
  check_write(state,"overlap",do.call(rbind,overlap))
  check_write(state,"balance",do.call(rbind,balance))
  invisible(state)
}

# Step 3: pre-treatment ATT placebos with simultaneous multiplier bands over all 34 outcomes.
check_placebos <- function(state) {
  tabs <- list()
  for (nm in names(state$samples)) for (cov in c(FALSE,TRUE)) {
    s <- state$samples[[nm]]; fit <- state$baseline[[paste(nm,cov,sep="|")]]
    stat <- check_att(s$y[,-1,drop=FALSE],s$d,fit$yfit[,-1,drop=FALSE],fit$ps)
    scale <- stat$se*sqrt(length(s$d))
    multiplier <- check_seed(20260916,matrix(rnorm(length(s$d)*1999),length(s$d),1999))
    boot <- sweep(crossprod(stat$influence,multiplier)/sqrt(length(s$d)),1,scale,"/")
    maximum <- apply(abs(boot),2,max)
    critical <- as.numeric(quantile(maximum,.95))
    tstat <- stat$att/stat$se
    tabs[[paste(nm,cov)]] <- data.frame(sample=nm,covariates=cov,outcome=colnames(s$y)[-1],
      att=stat$att,se=stat$se,pointwise_lower=stat$att-1.96*stat$se,
      pointwise_upper=stat$att+1.96*stat$se,p_normal=2*pnorm(-abs(tstat)),
      simultaneous_lower=stat$att-critical*stat$se,simultaneous_upper=stat$att+critical*stat$se,
      p_family_adjusted=sapply(abs(tstat),function(t) (1+sum(maximum>=t))/2000),
      family_critical=critical,multiplier_draws=1999,seed=20260916)
  }
  check_write(state,"placebos",do.call(rbind,tabs))
  invisible(state)
}

# Step 3: broaden sensitivities, distinguishing control and treated-score trimming.
check_sensitivity <- function(state) {
  # Sensitivity choices at fixed K: dimension rules, distances, ordered splits, and trimming.
  specs <- list(benchmark=list(),threshold075=list(multiplier=.75),threshold125=list(multiplier=1.25),
    squaredEuclidean=list(distance="squared-Euclidean"),averages=list(distance="averages"),
    reversed=list(split="reversed"),boundary100=list(split="boundary100"),boundary150=list(split="boundary150"),
    trim90=list(trim=.9),trim95=list(trim=.95))
  rows <- list()
  summarize_fit <- function(nm, specification, covariates, s, fit,
                            dropped_controls=0L, dropped_treated=0L) {
    stat <- check_att(s$y[,1],s$d,fit$yfit[,1],fit$ps)
    ctrl <- s$d==0; treated <- !ctrl
    w <- fit$ps[ctrl]/(1-fit$ps[ctrl])
    lo <- min(fit$ps[ctrl]); hi <- max(fit$ps[ctrl])
    data.frame(sample=nm,specification=specification,covariates=covariates,
      n=length(s$d),treated=sum(s$d),dropped_controls=dropped_controls,
      dropped_treated=dropped_treated,K=state$K,
      att=stat$att,se=stat$se,lower=stat$att-1.96*stat$se,
      upper=stat$att+1.96*stat$se,
      treated_ps_mean=mean(fit$ps[treated]),treated_ps_min=min(fit$ps[treated]),
      treated_ps_max=max(fit$ps[treated]),control_ps_mean=mean(fit$ps[ctrl]),
      control_ps_min=lo,control_ps_max=hi,
      treated_outside_control_range=sum(fit$ps[treated]<lo | fit$ps[treated]>hi),
      control_ESS=sum(w)^2/sum(w^2),max_control_weight=max(w),
      nonconverged=sum(!fit$records$converged),
      nonconverged_treated=sum(!fit$records$converged[treated]),
      rank_deficient=sum(fit$records$outcome_rank<fit$records$parameters |
                         fit$records$logit_rank<fit$records$parameters),
      nonfinite_predictions=sum(!is.finite(fit$yfit[,1])),
      nonfinite_scores=sum(!is.finite(fit$ps)),
      near_one_control_score=sum(fit$ps[ctrl]>1-1e-8),
      logit_warning_count=sum(nzchar(fit$records$logit_warning)),
      mean_dimension=mean(fit$records$dimension))
  }
  for (nm in names(state$samples)) {
    original <- state$samples[[nm]]
    for (spec in names(specs)) {
      message("Sensitivity: ",nm," / ",spec)
      options <- specs[[spec]]; s <- original; tag <- nm; dropped <- 0L
      if (!is.null(options$trim)) {
        discrepancy <- check_neighbors(state,s,nm)$tab$normalized_max
        cutoff <- quantile(discrepancy[s$d==0],options$trim)
        keep <- s$d==1 | discrepancy<cutoff
        dropped <- sum(!keep); tag <- paste(nm,spec,sep="_")
        s <- lapply(s,function(v) if(is.null(dim(v))) v[keep] else v[keep,,drop=FALSE])
        stopifnot(sum(s$d)==sum(original$d)); options$trim <- NULL
      }
      for (cov in c(FALSE,TRUE)) {
        fit <- if(spec=="benchmark") state$baseline[[paste(nm,cov,sep="|")]] else
          do.call(check_fit,c(list(state=state,s=s,tag=tag,covariates=cov),options))
        rows[[paste(nm,spec,cov)]] <- summarize_fit(nm,spec,cov,s,fit,
          dropped_controls=dropped)
      }
    }
    # Treated-score sensitivity: remove benchmark-treated firms with a
    # covariate-adjusted propensity above 0.9, then refit everything.
    baseline <- state$baseline[[paste(nm,TRUE,sep="|")]]
    drop <- original$d==1 & baseline$ps>.9
    stopifnot(sum(drop)==1L)
    keep <- !drop
    trimmed <- lapply(original,function(v)
      if(is.null(dim(v))) v[keep] else v[keep,,drop=FALSE])
    fit <- check_fit(state,trimmed,paste0(nm,"_propensity_trim09"),covariates=TRUE)
    rows[[paste(nm,"propensity_trim09",TRUE)]] <- summarize_fit(
      nm,"propensity_trim09",TRUE,trimmed,fit,dropped_treated=sum(drop))
  }
  check_write(state,"sensitivity",do.call(rbind,rows))
  invisible(state)
}

# Diagnostic figures show every active specification with explicit failure labels.
check_figures <- function(state,only=NULL) {
  draw <- function(name,expr,layout=c(2,2)) {
    if(!is.null(only) && !name %in% only) return(invisible(NULL))
    png(file.path(state$output,paste0("check-",name,".png")),width=1600,height=1100,res=160)
    on.exit(dev.off())
    par(mfrow=layout,mar=c(4,4,3,1),oma=c(0,0,2,0),cex.axis=.8,cex.main=.9)
    force(expr)
  }
  draw("matching",{
    for(nm in names(state$samples)) {
      v <- subset(state$tables$neighborhoods,sample==nm & split=="benchmark")
      boxplot(normalized_max~treated,data=v,names=c("Control","Treated"),xlab="",ylab="Normalized maximum distance",main=paste(nm,"matching"))
      boxplot(local_treated~treated,data=v,names=c("Control center","Treated center"),xlab="",ylab="Treated firms among K=99",main=paste(nm,"local treatment cells"))
    }
    mtext("Matching and treatment-cell counts",outer=TRUE,font=2)
  })
  draw("signals",{
    for(nm in names(state$samples)) {
      v <- subset(state$tables$dimensions,sample==nm)
      rmax <- max(v$selected_dimension)
      counts <- rbind(Control=tabulate(v$selected_dimension[v$treated==0],rmax),Treated=tabulate(v$selected_dimension[v$treated==1],rmax))
      barplot(t(counts/rowSums(counts)),beside=TRUE,col=gray.colors(rmax),legend.text=paste("r =",seq_len(rmax)),ylab="Fraction of centers",main=paste(nm,"selected dimensions"))
      v <- subset(state$tables$signals,sample==nm & component<=8)
      boxplot(eigen_threshold_ratio~component,data=v,outline=FALSE,xlab="Component",ylab="Eigenvalue / threshold squared",main=paste(nm,"local spectrum")); abline(h=1,lty=2)
    }
    mtext("Complete spectra in CSV; selection threshold is a scale heuristic",outer=TRUE,font=2)
  })
  draw("splits",{
    for(nm in names(state$samples)) {
      v <- subset(state$tables[["split-stability"]],sample==nm)
      v$split <- factor(v$split,levels=c("reversed","boundary100","boundary150"))
      boxplot(overlap_fraction~split,data=v,names=c("Swap halves","100 / 150","150 / 100"),xlab="",ylim=c(0,1),ylab="Neighbor overlap",main=paste(nm,"ordered proxy splits")); abline(h=median(v$chance_overlap),lty=2)
      v <- subset(state$tables[["factor-space"]],sample==nm)
      v$split <- factor(v$split,levels=c("reversed","boundary100","boundary150"))
      boxplot(projection_distance~split,data=v,names=c("Swap halves","100 / 150","150 / 100"),xlab="",ylim=c(0,1),ylab="Normalized projection distance",main=paste(nm,"fixed-neighborhood spaces"))
    }
    mtext("Dates remain chronological within each contiguous block",outer=TRUE,font=2)
  })
  draw("prediction",{
    for(nm in names(state$samples)) {
      v <- subset(state$tables[["prediction-gains"]],sample==nm & group=="treated")
      ratios <- cbind(Linear=v$factor_controls_mse/v$covariates_linear_mse,KNN=v$factor_controls_mse/v$covariates_knn_mse)
      matplot(1:2,ratios,pch=c(16,17),type="b",col=c("black","gray50"),xaxt="n",xlim=c(.7,2.3),ylim=range(c(1,ratios)),ylab="MSE ratio: PCA + controls / covariates",xlab="",main=paste(nm,"treated-firm proxy prediction"))
      axis(1,1:2,c("Later: [-80,-31]","Earlier: [-280,-231]"));abline(h=1,lty=2);legend("center",c("Covariate linear","Covariate KNN"),pch=c(16,17),col=c("black","gray50"),bty="n",cex=.75)
    }
    for(nm in names(state$samples)) {
      v <- subset(state$tables[["outcome-gains"]],sample==nm)
      barplot(sqrt(c(v$covariates_linear_mse,v$covariates_knn_mse,v$factor_controls_mse)),names.arg=c("Cov linear","Cov KNN","PCA + controls"),col=c("gray75","gray50","gray25"),ylab="Control-outcome RMSE",main=paste(nm,"event-outcome validation"))
    }
    mtext("Five-fold prediction: treated pre-returns above; untreated event outcomes below",outer=TRUE,font=2)
  })
  draw("overlap",{
    for(nm in names(state$samples)) for(cov in c(FALSE,TRUE)) {
      v <- subset(state$tables[["local-fits"]],sample==nm & covariates==cov)
      ov <- subset(state$tables$overlap,sample==nm & covariates==cov)
      boxplot(ps~treated,data=v,names=c("Control","Treated"),xlab="",ylim=c(0,1),ylab="Estimated propensity",main=paste(nm,if(cov)"with controls" else "without controls",sprintf("(control ESS %.1f)",ov$control_ESS)))
    }
    mtext("Near-zero control scores are not an ATT failure",outer=TRUE,font=2)
  })
  draw("preoutcomes",{
    for(nm in names(state$samples)) for(cov in c(FALSE,TRUE)) {
      v <- subset(state$tables$placebos,sample==nm & covariates==cov & grepl("^day",outcome))
      day <- as.numeric(sub("day","",v$outcome))
      plot(day,v$att,type="b",pch=16,ylim=range(v$simultaneous_lower,v$simultaneous_upper),xlab="Pre-event day",ylab="Adjusted daily-return difference",main=paste(nm,if(cov)"with controls" else "without controls"))
      segments(day,v$simultaneous_lower,day,v$simultaneous_upper,col="gray50");abline(h=0,lty=2)
    }
    mtext("These 30 days were excluded from proxies; bands cover all 34 pre-outcomes",outer=TRUE,font=2)
  })
  draw("sensitivity",{
    par(mar=c(4,8,3,1))
    for(nm in names(state$samples)) {
      v <- subset(state$tables$sensitivity,sample==nm)
      specs <- unique(v$specification); idx <- match(v$specification,specs); shift <- ifelse(v$covariates,.13,-.13)
      good <- v$near_one_control_score==0 & is.finite(v$att) & is.finite(v$se)
      lim <- range(c(0,v$lower[good],v$upper[good]),finite=TRUE)
      plot(v$att[good],(idx+shift)[good],xlim=lim,ylim=c(.5,length(specs)+.5),pch=ifelse(v$covariates[good],17,16),col=ifelse(v$covariates[good],"gray45","black"),yaxt="n",xlab="ATT and nominal 95% interval",ylab="",main=paste(nm,"sample"))
      axis(2,seq_along(specs),specs,las=2,cex.axis=.8)
      segments(v$lower[good],(idx+shift)[good],v$upper[good],(idx+shift)[good],col=ifelse(v$covariates[good],"gray45","black"));abline(v=0,lty=2)
      if(any(!good)) text(mean(lim),(idx+shift)[!good],"Numerical issue: see full CSV",cex=.7)
    }
    mtext("Sensitivity analyses at K=99; black: without controls, gray: with controls",outer=TRUE,font=2)
  },layout=c(1,2))
  invisible(state)
}

check_prediction_summary <- function(tab,outcome=FALSE) {
  groups <- list()
  for(nm in unique(tab$sample)) for(des in unique(tab$design)) {
    v <- subset(tab,sample==nm & design==des)
    for(group in if(outcome) "controls" else c("all","treated")) {
      x <- if(group=="treated") v[v$treated==1,] else v
      cols <- c("factor_mse","factor_controls_mse","mean_mse","covariates_linear_mse","covariates_knn_mse")
      stopifnot(nrow(x)>0,all(is.finite(as.matrix(x[,cols]))))
      m <- colMeans(x[,cols])
      groups[[paste(nm,des,group)]] <- data.frame(sample=nm,design=des,group=group,units=nrow(x),
        factor_mse=m[1],factor_controls_mse=m[2],mean_mse=m[3],covariates_linear_mse=m[4],covariates_knn_mse=m[5],
        gain_vs_covariates_linear=100*(1-m[2]/m[4]),gain_vs_covariates_knn=100*(1-m[2]/m[5]),
        secondary_factor_gain_vs_neighbor_mean=100*(1-m[1]/m[3]),row.names=NULL)
    }
  }
  do.call(rbind,groups)
}

# Concise result report; all active results remain in the accompanying CSVs.
check_report <- function(state) {
  check_write(state,"prediction-gains",check_prediction_summary(state$tables[["heldout-prediction"]]))
  check_write(state,"outcome-gains",check_prediction_summary(state$tables[["outcome-prediction"]],TRUE))
  # Primary reporting focus: treated firms in the later block of the gapless design.
  focus <- subset(state$tables[["prediction-gains"]],design=="chronological_gapless" & group=="treated")
  focus$prediction_first_day <- -80L
  focus$prediction_last_day <- -31L
  focus$prediction_days <- 50L
  stopifnot(nrow(focus)==length(state$samples),all(focus$prediction_days==50L))
  check_write(state,"treated-proxy-focus",focus)
  tabs <- state$tables
  summary <- do.call(rbind,lapply(names(state$samples),function(nm) {
    dims <- subset(tabs$dimensions,sample==nm); event <- subset(tabs$sensitivity,sample==nm & specification=="benchmark")
    ov <- subset(tabs$overlap,sample==nm); pre <- subset(tabs$placebos,sample==nm & grepl("^day",outcome))
    data.frame(sample=nm,n=nrow(dims),treated=sum(dims$treated),dimension1=sum(dims$selected_dimension==1),dimension2=sum(dims$selected_dimension==2),
      forced_one=sum(dims$forced_one),cap_binds=sum(dims$cap_binds),
      ATT_no_controls=event$att[!event$covariates],SE_no_controls=event$se[!event$covariates],
      ATT_with_controls=event$att[event$covariates],SE_with_controls=event$se[event$covariates],
      ESS_no_controls=ov$control_ESS[!ov$covariates],ESS_with_controls=ov$control_ESS[ov$covariates],
      near_one_control_scores_with_controls=ov$near_one_control_score[ov$covariates],
      significant_daily_preoutcomes_no_controls=sum(!pre$covariates & pre$p_family_adjusted<.05),
      significant_daily_preoutcomes_with_controls=sum(pre$covariates & pre$p_family_adjusted<.05))
  }))
  check_write(state,"summary",summary)
  lines <- c("# Empirical diagnostic audit","",
    "Run from Application/: `FENG_CHECKS_ONLY=1 Rscript --vanilla Feng-2026_Illustration.R`. These intermediate files provide the numerical provenance for the application results reported in the paper and supplement.","",
    "The original full/base samples, CAR[0,1], K=99, controls, dimension threshold, intercepts and proxy window [-280,-31] are preserved. No new K grid or forced fixed-factor experiment is included.","",
    "## Benchmark and local signals","",
    "| Sample | Firms / treated | Selected dimensions 1 / 2 | ATT without controls (SE) | ATT with controls (SE) |","| --- | ---: | ---: | ---: | ---: |")
  for(i in seq_len(nrow(summary))) {
    v <- summary[i,]
    lines <- c(lines,sprintf("| %s | %d / %d | %d / %d | %.4f (%.4f) | %.4f (%.4f) |",v$sample,v$n,v$treated,v$dimension1,v$dimension2,v$ATT_no_controls,v$SE_no_controls,v$ATT_with_controls,v$SE_with_controls))
  }
  lines <- c(lines,"","All treated centers select one factor. The four-component cap does not bind; thirteen control neighborhoods in each sample retain one component despite no signal above the rule's threshold. Exact singular values and eigenvalues normalized by neighborhood size times PCA-column count are available for every center. The selection threshold is a local-data-scale heuristic, not an independently estimated noise variance.","",
    "## Ordered splits and sensitivity analyses","",
    "Every matching/PCA group is one contiguous chronological block. The benchmark uses 125/125 columns; alternatives swap halves or move the boundary to 100/150 and 150/100. Dates are never shuffled, interleaved, or reversed within a group. Fixed-neighborhood projection distances compare factor spaces without depending on rotations.","")
  for(nm in names(state$samples)) {
    v <- subset(tabs[["split-stability"]],sample==nm & treated==1)
    med <- tapply(v$overlap_fraction,v$split,median)
    lines <- c(lines,sprintf("- %s: treated-center median neighbor overlaps range %.1f-%.1f%%, versus chance overlap %.1f%%.",nm,100*min(med),100*max(med),100*median(v$chance_overlap)))
    for(cov in c(FALSE,TRUE)) {
      v <- subset(tabs$sensitivity,sample==nm & covariates==cov)
      good <- v$near_one_control_score==0 & is.finite(v$att) & is.finite(v$se)
      lines <- c(lines,sprintf("- %s, %s controls: ATT across the %d specifications ranges %.4f-%.4f; %d nominal 95%% intervals exclude zero. Specifications with nonfinite ATT/SE or control scores above 1-1e-8: %d (full results retained in CSV).",nm,if(cov)"with" else "without",nrow(v),min(v$att[good]),max(v$att[good]),sum(v$lower[good]>0 | v$upper[good]<0),sum(!good)))
    }
  }
  lines <- c(lines,"","Threshold multipliers 0.75/1.25, the two alternative implemented distances, and 90th/95th-percentile poor-match control trimming supplement split sensitivity. Poor-match trimming retains every treated firm. The separate propensity-score check drops the one treated firm above 0.9 in each benchmark sample. Each trim recomputes neighborhoods and nuisances. These are illustrative sensitivity comparisons, not evidence that every alternative distance meets the identifying assumptions.","",
    "## Prediction beyond observed covariates","",
    "Five-fold held-out prediction compares local PCA plus the same controls with (i) covariate-only linear regression and (ii) K=99 nearest-neighbor means using standardized observed covariates alone. Covariate-only models never use proxy-based neighborhoods. PCA and regressions exclude the entire test fold. The three proxy groups are disjoint contiguous chronological blocks and use all 250 dates.","",
    "A positive MSE gain means improvement over that covariate-only comparator; a negative gain means higher prediction error. These are descriptive gains, without significance claims.","",
    "The primary reported comparison focuses on treated firms in the later prediction block, days [-80,-31]. This is relevant to the ATT target population and uses returns closer to treatment. Matching uses [-280,-181], PCA uses [-180,-81], and the last 30 pre-event days remain excluded. Each prediction date has its own regression; MSE averages squared errors over the 50 dates and treated firms, rather than summing returns into a CAR.","",
    "| Sample | Treated firms | Prediction days | MSE gain vs covariate linear | MSE gain vs covariate KNN |","| --- | ---: | --- | ---: | ---: |")
  for(i in seq_len(nrow(focus))) {
    v <- focus[i,]
    lines <- c(lines,sprintf("| %s | %d | [-80,-31] | %.2f%% | %.2f%% |",v$sample,v$units,v$gain_vs_covariates_linear,v$gain_vs_covariates_knn))
  }
  lines <- c(lines,"","The later-block focus and the gapless 100/100/50 allocation were fixed before these gains were computed. Folds, K and model settings were preserved; no alternative windows were searched. [Treated-focus CSV](check-treated-proxy-focus.csv) records absolute MSEs and sample/date definitions. The reversed-block treated and all-firm results remain in the prediction CSV and figure regardless of sign.","",
    "All-firm averages remain supporting comparisons:","",
    "| Sample | Held-out design | MSE gain vs covariate linear | MSE gain vs covariate KNN |","| --- | --- | ---: | ---: |")
  for(i in which(tabs[["prediction-gains"]]$group=="all")) {
    v <- tabs[["prediction-gains"]][i,]
    lines <- c(lines,sprintf("| %s | %s | %.2f%% | %.2f%% |",v$sample,v$design,v$gain_vs_covariates_linear,v$gain_vs_covariates_knn))
  }
  lines <- c(lines,"","[Prediction-gain CSV](check-prediction-gains.csv) also reports treated-firm results, absolute MSEs, and the secondary factor-versus-neighbor-mean decomposition. That secondary comparison measures PCA's increment after matching, rather than the full gain beyond observed controls.","",
    "At the single benchmark K, validation of untreated event outcomes uses control firms' observed outcomes; treated counterfactual errors are not observable:","",
    "| Sample | RMSE: covariate linear | RMSE: covariate KNN | RMSE: PCA + controls |","| --- | ---: | ---: | ---: |")
  for(i in seq_len(nrow(tabs[["outcome-gains"]]))) {
    v <- tabs[["outcome-gains"]][i,]
    lines <- c(lines,sprintf("| %s | %.4f | %.4f | %.4f |",v$sample,sqrt(v$covariates_linear_mse),sqrt(v$covariates_knn_mse),sqrt(v$factor_controls_mse)))
  }
  lines <- c(lines,"","## ATT overlap, balance and pre-event differences","",
    "Near-zero control propensity is not counted as an ATT problem. The summaries focus on treated comparability and upper-tail control weights. ESS describes weight concentration in the residual correction, not a literal sample size for the full doubly robust estimator.","",
    "| Sample, with controls | Treated outside control-score range | Max control score | Max control weight | Control ESS | Top-five weight share |","| --- | ---: | ---: | ---: | ---: | ---: |")
  for(i in which(tabs$overlap$covariates)) {
    v <- tabs$overlap[i,]
    lines <- c(lines,sprintf("| %s | %d / %d | %.3f | %.2f | %.1f | %.1f%% |",v$sample,v$treated_outside_control_range,v$treated_count,v$control_ps_max,v$max_control_weight,v$control_ESS,100*v$top5_control_weight_share))
  }
  lines <- c(lines,"","These concentrations and sample-range comparisons should be described briefly rather than called a definitive overlap failure. Full convergence records, matching discrepancies, treatment-cell counts, and balance in controls and unused returns remain in CSVs.","",
    "With controls, 10/38 benchmark local logits report nonconvergence in the full/base samples, all at control centers; center predictions remain finite. This audit flag is separate from treating near-zero control scores as an ATT problem.","")
  for(nm in names(state$samples)) {
    v <- subset(tabs$balance,sample==nm & scheme=="Local.PCA.with.controls" & variable %in% c("Log assets","ROE","Leverage"))
    lines <- c(lines,sprintf("- %s: adjusted absolute standardized differences in log assets / ROE / leverage are %.3f / %.3f / %.3f; balance is not uniformly close to zero.",nm,abs(v$smd[1]),abs(v$smd[2]),abs(v$smd[3])))
  }
  lines <- c(lines,"",
    "The 30 daily and four cumulative pre-outcomes lie in [-30,-1], the period excluded from the proxy sample. Significance can reflect anticipation, earlier news, or other differences; these comparisons are not an automatic rejection of the event-outcome specification. We retain the original window and adjust jointly across all 34 outcomes within each sample/control specification.","")
  for(i in seq_len(nrow(summary))) {
    v <- summary[i,]
    lines <- c(lines,sprintf("- %s: %d significant daily pre-outcomes without controls and %d with controls after family adjustment.",v$sample,v$significant_daily_preoutcomes_no_controls,v$significant_daily_preoutcomes_with_controls))
  }
  significant <- subset(tabs$placebos,covariates & p_family_adjusted<.05)
  for(i in seq_len(nrow(significant))) {
    v <- significant[i,]
    lines <- c(lines,sprintf("- %s with controls, %s: difference %.4f (SE %.4f), family-adjusted p=%.4f.",v$sample,v$outcome,v$att,v$se,v$p_family_adjusted))
  }
  lines <- c(lines,"","No cumulative window rejects at 5% in the current results. Nonrejection is not proof of no pre-event differences.","",
    "## Outputs and interpretation","",
    "The examples assess observed stability and fit; they do not directly verify latent assumptions. The paper-ready tables and figures are written to the parent `output/` directory.","",
    "Entry inputs and IDs: check-samples.csv / check-results.rds. Main tables: check-summary.csv, check-sensitivity.csv, check-treated-proxy-focus.csv, check-prediction-gains.csv, check-outcome-gains.csv, check-overlap.csv, check-balance.csv, check-placebos.csv. Geometry: check-neighborhoods.csv, check-signals.csv, check-dimensions.csv, check-factor-space.csv, check-split-stability.csv. Seven PNGs visualize these results. check-local-fits.csv records numerical warnings. check-sessionInfo.txt records package/platform dependence.","",
    "Seeds retained: 20260917 for five unit folds and 20260916 for 1,999 independent-unit Gaussian multipliers. Proxy splits are deterministic. Multipliers use the existing influence-function formula without refitting nuisances; the inference does not establish cross-unit independence. No package installation or estimator regularization was performed.")
  writeLines(lines,file.path(state$output,"check-report.md"))
  invisible(summary)
}

# Preserve audit inputs, seeds, and predictions alongside the reports.
check_finish <- function(state) {
  check_report(state)
  check_figures(state)
  write_application_outputs(state)
  saveRDS(list(tables=state$tables,baseline=state$baseline,splits=state$splits,
    folds=state$folds,ids=lapply(state$samples,`[[`,"ids"),K=state$K,
    seeds=c(placebos=20260916,unit_folds=20260917)),file.path(state$output,"check-results.rds"))
  writeLines(capture.output(sessionInfo()),file.path(state$output,"check-sessionInfo.txt"))
  message("Application diagnostics saved under ",state$output)
  invisible(state)
}

# Write a LaTeX tabular fragment that can be copied directly into the paper input folder.
application_write_table <- function(output, name, lines) {
  path <- file.path(output,name)
  writeLines(lines,path)
  path
}

# Recreate the main-paper distribution of the second local singular-value ratio.
application_figure_singular_ratio <- function(signals, output) {
  values <- subset(signals,sample=="full" & component==2L)
  stopifnot(nrow(values)==583L,!anyDuplicated(values$unit),
    all(is.finite(values$singular_threshold_ratio)),
    sum(values$singular_threshold_ratio>1)==23L)
  figure <- ggplot(data.frame(ratio=values$singular_threshold_ratio),aes(x=ratio))+
    geom_histogram(aes(y=after_stat(count/sum(count))),binwidth=.05,boundary=1,
      fill="grey80",color="grey35",linewidth=.3)+
    labs(x="Second singular value / selection threshold",
      y="Fraction of neighborhoods")+
    theme_bw()+theme(panel.grid.minor=element_blank())
  path <- file.path(output,"singular_ratio_full.pdf")
  ggsave(path,figure,width=5,height=3,device="pdf",useDingbats=FALSE)
  path
}

# Recreate the supplemental matching-discrepancy distributions for both samples.
application_figure_matching <- function(neighborhoods, output) {
  neighborhoods <- subset(neighborhoods,sample %in% c("full","base") &
    split=="benchmark" & distance=="pseudo-max" & K==99)
  get_values <- function(sample,treated_only) {
    x <- neighborhoods[neighborhoods$sample==sample,]
    expected <- if(sample=="full") 583L else 526L
    stopifnot(nrow(x)==expected,!anyDuplicated(x$unit),
      all(x$treated %in% 0:1),all(is.finite(x$normalized_max)),
      all(x$normalized_max>=0))
    if(treated_only) {
      x <- x[x$treated==1,]
      stopifnot(nrow(x)==if(sample=="full") 22L else 12L)
    }
    x$normalized_max
  }
  samples <- c("full","base"); columns <- c("All firms","Treated firms")
  values <- lapply(c(FALSE,TRUE),function(treated_only)
    lapply(samples,get_values,treated_only=treated_only))
  breaks <- list(pretty(range(unlist(values[[1]])),n=20),
    pretty(range(unlist(values[[2]])),n=8))
  heights <- lapply(seq_along(breaks),function(column)
    max(unlist(lapply(values[[column]],function(x)
      hist(x,breaks=breaks[[column]],plot=FALSE)$density))))
  path <- file.path(output,"matching_discrepancy_dist.pdf")
  pdf(path,width=6.6,height=4.2,useDingbats=FALSE)
  tryCatch({
    par(mfrow=c(2,2),mar=c(3,3.5,1.6,.6),mgp=c(2,.6,0),
      cex.axis=.8,cex.lab=.85,cex.main=.9)
    for(row in seq_along(samples)) for(column in seq_along(columns)) {
      x <- values[[column]][[row]]
      hist(x,breaks=breaks[[column]],freq=FALSE,xlim=range(breaks[[column]]),
        ylim=c(0,1.06*heights[[column]]),col="gray75",border="black",
        xlab="Normalized maximum discrepancy",ylab="Density",
        main=sprintf("%s: %s (n=%d)",if(row==1L) "Full" else "Base",
          columns[[column]],length(x)))
    }
  },finally=dev.off())
  path
}

# Recreate the supplemental distributions of the first two local eigenvalue ratios.
application_figure_eigen_ratios <- function(signals, output) {
  path <- file.path(output,"eigen_ratio_dist.pdf")
  pdf(path,width=6.6,height=2.7,useDingbats=FALSE)
  tryCatch({
    par(mfrow=c(1,2),mar=c(3.1,3.5,1.7,.6),mgp=c(2,.6,0),
      cex.axis=.85,cex.lab=.9,cex.main=.9)
    for(component in 1:2) {
      values <- lapply(c("full","base"),function(sample) {
        x <- signals[signals$sample==sample & signals$component==component,]
        expected <- if(sample=="full") 583L else 526L
        stopifnot(nrow(x)==expected,all(is.finite(x$eigen_threshold_ratio)),
          !anyDuplicated(x$unit))
        x$eigen_threshold_ratio
      })
      breaks <- pretty(range(c(unlist(values),1)),n=if(component==1L) 22 else 18)
      full_hist <- hist(values[[1]],breaks=breaks,plot=FALSE)
      base_hist <- hist(values[[2]],breaks=breaks,plot=FALSE)
      hist(values[[1]],breaks=breaks,freq=FALSE,xlim=range(breaks),
        ylim=c(0,1.08*max(full_hist$density,base_hist$density)),
        col="gray80",border="gray55",xlab="Eigenvalue ratio",ylab="Density",
        main=if(component==1L) "First eigenvalue ratio" else "Second eigenvalue ratio")
      hist(values[[2]],breaks=breaks,freq=FALSE,add=TRUE,
        col=NA,border="black",lwd=1.2)
      legend("topright",c("Full","Base"),fill=c("gray80",NA),
        border=c("gray55","black"),bty="n",cex=.75)
    }
  },finally=dev.off())
  path
}

# Write the paper-ready ATT sensitivity table from freshly computed checks.
application_write_sensitivity <- function(sensitivity, output="output") {
  pick_sensitivity <- function(sample_name,spec,with_covariates)
    sensitivity[sensitivity$sample==sample_name & sensitivity$specification==spec &
      sensitivity$covariates==with_covariates,]
  specs <- c("benchmark","threshold075","threshold125","reversed","boundary100",
    "boundary150","trim90","trim95")
  labels <- c("Benchmark 125/125","$0.75\\times$ dimension threshold",
    "$1.25\\times$ dimension threshold","Swapped halves","Matching 100, PCA 150",
    "Matching 150, PCA 100","Trim top 10\\%","Trim top 5\\%")
  rows <- vapply(seq_along(specs),function(i) {
    v <- lapply(c("full","base"),function(sample)
      lapply(c(FALSE,TRUE),function(covariates) {
        x <- pick_sensitivity(sample,specs[i],covariates)
        stopifnot(nrow(x)==1L)
        sprintf("%.3f (%.3f)",x$att,x$se)
      }))
    sprintf("%s & %s & %s & %s & %s \\\\",labels[i],v[[1]][[1]],v[[2]][[1]],
      v[[1]][[2]],v[[2]][[2]])
  },character(1))
  ps <- lapply(c("full","base"),function(sample) {
    x <- pick_sensitivity(sample,"propensity_trim09",TRUE)
    stopifnot(nrow(x)==1L,x$dropped_treated==1L)
    sprintf("%.3f (%.3f)",x$att,x$se)
  })
  rows <- c(rows,sprintf("Propensity trim $>0.9$ & -- & -- & %s & %s \\\\",ps[[1]],ps[[2]]))
  application_write_table(output,"SA_Table_Application_Sensitivity.txt",c(
    "\\begin{tabular}{lcccc}","\\hline\\hline",
    "& \\multicolumn{2}{c}{No covariates} & \\multicolumn{2}{c}{With covariates} \\\\",
    "\\cline{2-3}\\cline{4-5}","Specification & Full & Base & Full & Base \\\\",
    "\\hline",rows,"\\hline","\\end{tabular}"))
}

# Export every empirical table and figure currently reported in the paper or supplement.
write_application_outputs <- function(state, output="output") {
  dir.create(output,showWarnings=FALSE,recursive=TRUE)
  tabs <- state$tables
  required <- c("sensitivity","split-stability","factor-space","prediction-gains",
    "outcome-gains","placebos","dimensions","signals","balance","neighborhoods")
  stopifnot(all(required %in% names(tabs)))
  files <- c(
    application_figure_singular_ratio(tabs$signals,output),
    application_figure_matching(tabs$neighborhoods,output),
    application_figure_eigen_ratios(tabs$signals,output))

  files <- c(files,application_write_sensitivity(tabs$sensitivity,output))

  sensitivity <- tabs$sensitivity
  pick_sensitivity <- function(sample_name,spec,with_covariates)
    sensitivity[sensitivity$sample==sample_name & sensitivity$specification==spec &
      sensitivity$covariates==with_covariates,]

  # Neighbor and loading-space stability across ordered proxy splits.
  overlap <- tabs[["split-stability"]]; spaces <- tabs[["factor-space"]]
  stability_value <- function(tab,sample,split,treated,value) {
    keep <- tab$sample==sample & tab$split==split
    if(treated) keep <- keep & tab$treated==1
    median(tab[[value]][keep])
  }
  split_names <- c(reversed="Swapped halves",boundary100="Matching 100, PCA 150",
    boundary150="Matching 150, PCA 100")
  rows <- unlist(lapply(c("full","base"),function(sample)
    vapply(names(split_names),function(split) sprintf(
      "%s & %s & %.3f & %.3f & %.3f & %.3f \\\\",
      if(sample=="full") "Full" else "Base",split_names[[split]],
      stability_value(overlap,sample,split,TRUE,"overlap_fraction"),
      stability_value(overlap,sample,split,FALSE,"overlap_fraction"),
      stability_value(spaces,sample,split,TRUE,"projection_distance"),
      stability_value(spaces,sample,split,FALSE,"projection_distance")),character(1))))
  files <- c(files,application_write_table(output,"SA_Table_Application_Split_Stability.txt",c(
    "\\begin{tabular}{llcccc}","\\hline\\hline",
    "& & \\multicolumn{2}{c}{Neighbor overlap} & \\multicolumn{2}{c}{Loading-space distance} \\\\",
    "\\cline{3-4}\\cline{5-6}",
    "Sample & Split & Treated & All & Treated & All \\\\","\\hline",rows,
    "\\hline","\\end{tabular}")))

  # Held-out proxy and untreated-event prediction gains.
  proxy <- subset(tabs[["prediction-gains"]],design=="chronological_gapless" & group=="treated")
  event <- tabs[["outcome-gains"]]
  stopifnot(nrow(proxy)==2L,nrow(event)==2L)
  rows <- c(vapply(c("full","base"),function(sample) {
      x <- proxy[proxy$sample==sample,]
      sprintf("%s & Treated daily proxies & %d & %.2f & %.2f \\\\",
        if(sample=="full") "Full" else "Base",x$units,
        x$gain_vs_covariates_linear,x$gain_vs_covariates_knn)
    },character(1)),vapply(c("full","base"),function(sample) {
      x <- event[event$sample==sample,]
      sprintf("%s & Untreated CAR[0,1] & %d & %.2f & %.2f \\\\",
        if(sample=="full") "Full" else "Base",x$units,
        x$gain_vs_covariates_linear,x$gain_vs_covariates_knn)
    },character(1)))
  files <- c(files,application_write_table(output,"SA_Table_Application_Prediction.txt",c(
    "\\begin{tabular}{llrcc}","\\hline\\hline",
    "Sample & Prediction target & Firms & Linear gain (\\%) & $K$-NN gain (\\%) \\\\",
    "\\hline",rows,"\\hline","\\end{tabular}")))

  # Simultaneously adjusted pre-event outcome checks.
  rows <- unlist(lapply(c("full","base"),function(sample)
    vapply(c(FALSE,TRUE),function(covariates) {
      x <- tabs$placebos[tabs$placebos$sample==sample &
        tabs$placebos$covariates==covariates,]
      sprintf("%s & %s & %d & %d \\\\",if(sample=="full") "Full" else "Base",
        if(covariates) "Yes" else "No",
        sum(grepl("^day",x$outcome) & x$p_family_adjusted<.05),
        sum(!grepl("^day",x$outcome) & x$p_family_adjusted<.05))
    },character(1))))
  files <- c(files,application_write_table(output,"SA_Table_Application_Placebos.txt",c(
    "\\begin{tabular}{lccc}","\\hline\\hline",
    "Sample & Covariates & Daily rejections (30) & Cumulative rejections (4) \\\\",
    "\\hline",rows,"\\hline","\\end{tabular}")))

  # Selected dimensions and quartiles of the first two signal ratios.
  rows <- vapply(c("full","base"),function(sample) {
    dims <- tabs$dimensions[tabs$dimensions$sample==sample,]
    ratios <- vapply(1:2,function(component) {
      x <- tabs$signals$eigen_threshold_ratio[tabs$signals$sample==sample &
        tabs$signals$component==component]
      paste(sprintf("%.2f",quantile(x,c(.25,.5,.75),names=FALSE)),collapse=" / ")
    },character(1))
    sprintf("%s & %d & %d & %s & %s \\\\",if(sample=="full") "Full" else "Base",
      sum(dims$selected_dimension==1),sum(dims$selected_dimension==2),ratios[1],ratios[2])
  },character(1))
  files <- c(files,application_write_table(output,"SA_Table_Application_Signals.txt",c(
    "\\begin{tabular}{lcccc}","\\hline\\hline",
    "Sample & $\\ttd_i=1$ & $\\ttd_i=2$ & First eigenvalue ratio & Second eigenvalue ratio \\\\",
    "\\hline",rows,"\\hline","\\end{tabular}")))

  # Absolute standardized differences for the three observed covariates.
  variables <- c("Log assets","ROE","Leverage")
  rows <- unlist(lapply(c("full","base"),function(sample)
    vapply(c("Unadjusted","Local.PCA.with.controls"),function(scheme) {
      x <- tabs$balance[tabs$balance$sample==sample & tabs$balance$scheme==scheme &
        tabs$balance$variable %in% variables,]
      x <- x[match(variables,x$variable),]
      stopifnot(nrow(x)==3L)
      sprintf("%s & %s & %.3f & %.3f & %.3f \\\\",
        if(sample=="full") "Full" else "Base",
        if(scheme=="Unadjusted") "Unweighted" else "Reweighted",abs(x$smd[1]),
        abs(x$smd[2]),abs(x$smd[3]))
    },character(1))))
  files <- c(files,application_write_table(output,"SA_Table_Application_Balance.txt",c(
    "\\begin{tabular}{llccc}","\\hline\\hline",
    "Sample & Control weights & Log assets & ROE & Leverage \\\\","\\hline",rows,
    "\\hline","\\end{tabular}")))

  # ATT-relevant propensity summaries before and after 5-percent poor-match trimming.
  ps_range <- function(mean,min,max,control=FALSE) {
    lo <- if(control && min<.001) "$<0.001$" else sprintf("%.3f",min)
    hi <- if(max>.999) "$>0.999$" else sprintf("%.3f",max)
    sprintf("%.3f (%s, %s)",mean,lo,hi)
  }
  rows <- unlist(lapply(c("benchmark","trim95"),function(spec)
    vapply(c("full","base"),function(sample) {
      x <- pick_sensitivity(sample,spec,TRUE)
      stopifnot(nrow(x)==1L)
      label <- if(spec=="benchmark") {
        if(sample=="full") "Full" else "Base"
      } else sprintf("Trim 5\\%%, %s",sample)
      sprintf("%s & %s & %s & %d/%d & %.1f \\\\",label,
        ps_range(x$treated_ps_mean,x$treated_ps_min,x$treated_ps_max),
        ps_range(x$control_ps_mean,x$control_ps_min,x$control_ps_max,TRUE),
        x$treated_outside_control_range,x$treated,x$control_ESS)
    },character(1))))
  files <- c(files,application_write_table(output,"SA_Table_Application_Overlap.txt",c(
    "\\begin{tabular}{lcccc}","\\hline\\hline",
    "& \\multicolumn{2}{c}{Propensity: mean (min, max)} & \\multicolumn{2}{c}{ATT diagnostics} \\\\",
    "\\cline{2-3}\\cline{4-5}",
    "Sample & Treated & Controls & Beyond range & Control ESS \\\\","\\hline",rows,
    "\\hline","\\end{tabular}")))
  invisible(files)
}
