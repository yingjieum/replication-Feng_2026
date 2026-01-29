#####################################################
##### Causal Inference in Nonlinear Factor Models ###
############ Empirical Application ##################
########### Last updated: 01/21/2026 ################
#####################################################
# This file replicates the results in Feng (2026)
# Data source: Acemoglu et al (JFE, 2016) 
#####################################################
rm(list=ls())
library(R.matlab)
library(Rfast)
library(irlba)
library(RcppNumerical)
library(ggplot2)
library(reshape2)
library(Hmisc)
library(cobalt)
library(hdm)

source("Feng-2026_Illustration_funs.R")

# Load data
data <- readMat("Data.mat")
ret <- data$Re
connect <- data$num[,3] #1: Shared Board 2: NY Connection 3: Geithner Schedule 4: Geithner Schedule 2007 
d <- 1*(connect>0)
date1 <- data$GeiNomDat
x <- t(ret[(date1-281):(date1-32),])   # n=583 by p=250, pre-treatment history
y.t0 <- ret[date1-1,]                  # CAR0
y.t1 <- y.t0+ret[date1,]               # CAR[0,1]
p <- dim(x)[2]
ctrlvar <- data$num[,8:10]             # 3 control vars

# Estimation
# define pars
nlam   <- 4             # maximum number of local factors to be extracted
y.full <- y.t1
d.full <- d
x.full <- x
n      <- length(d.full)
range  <- 1:n
subset <- (d.full==0)
pr     <- mean(d.full)
# A data-driven benchmark choice
Klist  <- 80:200
K0     <- Klist[which.min(knn.cv(nfolds=n, y=y.full, x=x.full[,p,drop=F], k=Klist, type="R")$crit)]
# several choices of K based on K0
Kseq   <- ceiling(K0 * c(.5, 1, 2))

##########################################################
# Only for Step 1 and 2 illustration, choose one K=Ksea[3]

# Step 1: latent variables extraction
# Step 1.1: illustrate KNN matching
kmat      <- knn.index(A=x.full[,1:(p/2)], K=K0)
distmat   <- dist.knn(A=x.full[,1:(p/2)])
i         <- which(d.full==1)[1]
pos       <- kmat; diag(pos) <- F    # remove itself in summary
dist.norm <- sapply(1:n, function(i) max(distmat[pos[,i],i])/sd(distmat[-i,i])) # normalized distance

# summary of distance
treated   <- summary(dist.norm[d.full==1]); ctrl <- summary(dist.norm[d.full==0])
sum.table <- round(rbind(treated, ctrl), 3)

# Generate Table 1
colheads <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
rowname  <- c("Treated", "Control")
latex(sum.table, file=paste("output/MainPaper_Table_SumKNN", ".txt", sep = ""), 
      append=FALSE, table.env=FALSE, center="none", title="",
      n.cgroup=NULL, cgroup=NULL, colheads=colheads,
      n.rgroup=NULL, rgroup=NULL, rowname=rowname
      )

# label outliers
outlier <- which(dist.norm>=quantile(dist.norm[d.full==0], .9))


# Step 1.2: illustrate local PCA
index <- kmat[,i]   # length n
PC    <- lpca(A=x.full[,-(1:(p/2))], index=index, nlam=5, n=n, K=K0)
ggplot()+geom_line(data=data.frame(x=1:5, y=PC$sv), aes(x=x,y=y))+
         xlab("Component Number")+ylab("Singular value")+
         theme_bw() + theme(panel.grid.minor = element_blank())
ggsave("output/scree.pdf", width = 5, height = 3)


##################################################
# Step 2: local least squares

# fitted sequence from day -30 to 1
N1 <- sum(d.full)
mat.y <- mat.fit <- matrix(NA, 32, N1)
for (j in 1:32) {
   y.t    <- ret[date1-32+j,]
   result <- compute(range=range, y=y.t, d=d.full, x=x.full, subset=subset, 
                     K=K0, nlam=nlam, n=n, p=p, const=F, ctrlvar=ctrlvar)
   mat.fit[j,] <- result$yfit[d.full==1]
   mat.y[j,]   <- y.t[d.full==1]
}

##################################################
# Balance check of pre-treatment returns for the full sample with covariates
# returned PS above
w_att.me <- ifelse(d.full == 1, 1, result$ps/(1-result$ps))

# PS based on z
ps.z     <- glm.fit(x=cbind(1, ctrlvar), y=d.full, family=binomial())$fitted.values
w_att.z  <- ifelse(d.full == 1, 1, ps.z/(1-ps.z))

# Balance table
x.pre    <- t(ret[(date1-31):(date1-27),])
dat <- data.frame(x1=x.pre[,1], x2=x.pre[,2], x3=x.pre[,3], x4=x.pre[,4], x5=x.pre[,5])

bt <- bal.tab(d.full ~ x1 + x2 + x3 + x4 + x5, data=dat, 
              weights = data.frame(PS1 = w_att.me, PS2 = w_att.z), estimand = "ATT", un = TRUE)
print(bt)
baplot <- love.plot(bt, abs = TRUE, shapes = c(16, 17, 15),
                    var.names = c(x1 = "x[-30]", x2 = "x[-29]", x3 ="x[-28]", x4 ="x[-27]", x5 ="x[-26]"))  
baplot + scale_shape_discrete(name = "Weighting scheme", labels = c("Unadjusted",
                                                                    expression(ps(alpha, z)), expression(ps(z)))) + labs(x = "Absolute standardized mean difference") +
  theme(
    legend.position = c(0.98, 0.98),   # inside plot (x,y in [0,1])
    legend.justification = c(1, 1),    # anchor at top-right of legend box
    legend.background = element_rect(fill = "white", color = NA)
  ) + ggtitle(NULL)
ggsave(paste("output/love_plot", ".pdf", sep=""), width = 5, height = 4, dpi = 300)


########## plot fitted values for y
plot <- ggplot()
for (j in 1:N1) {
  plot.fit <- data.frame(y=mat.fit[-(1:10),j], t=-20:1, name="Fit") 
  plot.y   <- data.frame(y=mat.y[-(1:10),j],   t=-20:1, name="Y")
  if (j==1) {
    plot <- plot + geom_line(data=plot.fit, aes(x=t, y=y, colour=name))
    plot <- plot + geom_line(data=plot.y,   aes(x=t, y=y, colour=name))
    plot <- plot + scale_colour_manual(name="", values=c("black","grey"))
  } else {
    plot <- plot + geom_line(data=plot.fit, aes(x=t, y=y), colour="black")
    plot <- plot + geom_line(data=plot.y,   aes(x=t, y=y), colour="grey")
  }
}
plot <- plot + geom_vline(xintercept = -.5, linetype="dashed", color="black")+ ylim(-.3, .4) +
               xlab("date") + ylab("return") + 
               theme_bw() + 
               theme(legend.position = c(.07,0.15),
                     legend.background = element_rect(fill="transparent"), 
                     panel.grid.minor = element_blank(),panel.grid.major = element_blank())
ggsave("output/fit.pdf", width = 5, height = 4)

##################################################################
# Plot the result for one firm (NOT longer reported in the paper)
plot.single <- ggplot()
plot.fit <- data.frame(y=mat.fit[-(1:10),1], t=-20:1, name="Fit") 
plot.y   <- data.frame(y=mat.y[-(1:10),1],   t=-20:1, name="Y")
plot.single <- plot.single + geom_line(data=plot.fit, aes(x=t, y=y, colour=name))
plot.single <- plot.single + geom_line(data=plot.y,   aes(x=t, y=y, colour=name))
plot.single <- plot.single + scale_colour_manual(name="", values=c("black","grey"))

plot.single <- plot.single + geom_vline(xintercept = -.5, linetype="dashed", color="black")+ ylim(-.3, .4) +
               xlab("date") + ylab("return") + 
               theme_bw() + 
               theme(legend.position = c(.07,0.15),
                     legend.background = element_rect(fill="transparent"),
                     panel.grid.minor = element_blank(),panel.grid.major = element_blank())
ggsave("output/fit_single.pdf", width = 5, height = 4)


################################################
# Step 3: calculate ATT

# full sample
n <- length(y.full); subset <- (d.full==0); pr <- mean(d.full)
theta1 <- mean(y.full[d.full==1]); range <- 1:n


reg.table <- matrix(NA, 12, 4)
reg.table[c(1,3,5),1] <- Kseq
for  (j in 1:3) {
   K <- Kseq[j]
   result <- compute(range=range, y=y.full, d=d.full, x=x.full, subset=subset, 
                     K=K, nlam=nlam, n=n, p=p, const=F)
   psi <- (d.full * result$yfit + 
          (1-d.full)*(y.full-result$yfit)*result$ps/(1-result$ps))/pr # NOT influence fun
   theta0 <- mean(psi)   # point estimate
   tau    <- theta1 - theta0
   reg.table[2*j-1,2] <- format(round(tau, 3), nsmall=3)
   
   # influence fun.    
   varphi <- (d.full*(y.full-theta1)/pr) - 
             ((d.full*(result$yfit-theta0)+
              (1-d.full)*(y.full-result$yfit)*result$ps/(1-result$ps))/pr)
   se     <- sqrt(mean(varphi^2))/sqrt(n)
   reg.table[2*j,2] <- paste("(", round(se, 3),  ")", sep="")
   #reg.table[2*j,2] <- paste("[", -round(se*1.96,3), ", ", round(se*1.96,3), "]", sep="")
}
reg.table[c(2,4,6),1] <- ""

lasso.fit <- rlassoATET(x=x.full, d=d.full, y=y.full, bootstrap = "none", model = FALSE)
reg.table[7,2] <- format(round(lasso.fit$te, 3), nsmall=3)
reg.table[8,2] <- paste("(", round(lasso.fit$se, 3),  ")", sep="")

# subsample: drop control units with large matching discrepancy
y.sub <- y.full[-outlier]; d.sub <- d.full[-outlier]; x.sub <- x.full[-outlier,]
n <- length(y.sub); subset <- (d.sub==0); pr <- mean(d.sub)
theta1 <- mean(y.sub[d.sub==1]); range <- 1:n

for  (j in 1:3) {
   K <- Kseq[j]
   result <- compute(range=range, y=y.sub, d=d.sub, x=x.sub, subset=subset, 
                     K=K, nlam=nlam, n=n, p=p, const=F)
   psi <- (d.sub * result$yfit + 
          (1-d.sub)*(y.sub-result$yfit)*result$ps/(1-result$ps))/pr #NOT influence fun
   theta0 <- mean(psi)   # point estimate
   tau    <- theta1 - theta0
   reg.table[2*j-1,3] <- format(round(tau, 3), nsmall=3)
   
   # influence fn.    
   varphi <- (d.sub*(y.sub-theta1)/pr) - 
             ((d.sub*(result$yfit-theta0)+
              (1-d.sub)*(y.sub-result$yfit)*result$ps/(1-result$ps))/pr)
   se     <- sqrt(mean(varphi^2))/sqrt(n)
   reg.table[2*j,3] <- paste("(", round(se, 3),  ")", sep="")
   #reg.table[2*j,3] <- paste("[", -round(se*1.96,3), ", ", round(se*1.96,3), "]", sep="")
}
lasso.fit <- rlassoATET(x=x.sub, d=d.sub, y=y.sub, bootstrap = "none", model = FALSE)
reg.table[7,3] <- format(round(lasso.fit$te, 3), nsmall=3)
reg.table[8,3] <- paste("(", round(lasso.fit$se, 3),  ")", sep="")


# base sample: drop Citi-related
ind <- (data$CorrCiti<=as.numeric(data$CorrCitiTr))
y.sub <- y.full[ind]; d.sub <- d.full[ind]; x.sub <- x.full[ind,]
n <- length(y.sub); subset <- (d.sub==0); pr <- mean(d.sub)
theta1 <- mean(y.sub[d.sub==1]); range <- 1:n

for  (j in 1:3) {
   K <- Kseq[j]
   result <- compute(range=range, y=y.sub, d=d.sub, x=x.sub, subset=subset, 
                     K=K, nlam=nlam, n=n, p=p, const=F)
   psi <- (d.sub * result$yfit + 
          (1-d.sub)*(y.sub-result$yfit)*result$ps/(1-result$ps))/pr #NOT influence fun
   theta0 <- mean(psi)   # point estimate
   tau    <- theta1 - theta0
   reg.table[2*j-1,4] <- format(round(tau, 3), nsmall=3)
   
   # influence fn.    
   varphi <- (d.sub*(y.sub-theta1)/pr) - 
             ((d.sub*(result$yfit-theta0)+
              (1-d.sub)*(y.sub-result$yfit)*result$ps/(1-result$ps))/pr)
   se     <- sqrt(mean(varphi^2))/sqrt(n)
   reg.table[2*j,4] <- paste("(", round(se, 3),  ")", sep="")
   #reg.table[2*j,4] <- paste("[", -round(se*1.96,3), ", ", round(se*1.96,3), "]", sep="")
}
lasso.fit <- rlassoATET(x=x.sub, d=d.sub, y=y.sub, bootstrap = "none", model = FALSE)
reg.table[7,4] <- format(round(lasso.fit$te, 3), nsmall=3)
reg.table[8,4] <- paste("(", round(lasso.fit$se, 3),  ")", sep="")

reg.table[9,]  <- c("", "0.005", "", "0.060")
reg.table[10,] <- c("","[-0.029, 0.014]","", "[-0.068, 0.036]")
reg.table[11,] <- c("", "-", "", "0.061")
reg.table[12,] <- c("", "-", "", "[-0.050, 0.061]")


################################################
############ Add control vars ##################

# full sample
n <- length(y.full); subset <- (d.full==0); pr <- mean(d.full)
theta1 <- mean(y.full[d.full==1]); range <- 1:n

reg.table.ctr <- matrix(NA, 8, 4)
reg.table.ctr[c(1,3,5),1] <- Kseq; reg.table.ctr[c(2,4,6),1] <- ""


for  (j in 1:3) {
   K <- Kseq[j]
   result <- compute(range=range, y=y.full, d=d.full, x=x.full, subset=subset, 
                     K=K, nlam=nlam, n=n, p=p, const=F, ctrlvar=ctrlvar)
   psi <- (d.full * result$yfit + 
          (1-d.full)*(y.full-result$yfit)*result$ps/(1-result$ps))/pr #NOT influence fun
   theta0 <- mean(psi)   # point estimate
   tau    <- theta1 - theta0
   reg.table.ctr[2*j-1,2] <- format(round(tau, 3), nsmall=3)
   
# influence fun.    
   varphi <- (d.full*(y.full-theta1)/pr) - 
             ((d.full*(result$yfit-theta0)+
              (1-d.full)*(y.full-result$yfit)*result$ps/(1-result$ps))/pr)
   se     <- sqrt(mean(varphi^2))/sqrt(n)
   reg.table.ctr[2*j,2] <- paste("(", round(se, 3),  ")", sep="")
   #reg.table.ctr[2*j,2] <- paste("[", -round(se*1.96,3), ", ", round(se*1.96,3), "]", sep="")
}

lasso.fit <- rlassoATET(x=cbind(x.full, ctrlvar), d=d.full, y=y.full, bootstrap = "none", model = FALSE)
reg.table.ctr[7,2] <- format(round(lasso.fit$te, 3), nsmall=3)
reg.table.ctr[8,2] <- paste("(", round(lasso.fit$se, 3),  ")", sep="")

###################################
# subsample: drop control units with large matching discrepancy
y.sub <- y.full[-outlier]; d.sub <- d.full[-outlier]; x.sub <- x.full[-outlier,]
n <- length(y.sub); subset <- (d.sub==0); pr <- mean(d.sub); ctrl.sub <- ctrlvar[-outlier,]
theta1 <- mean(y.sub[d.sub==1]); range <- 1:n


for  (j in 1:3) {
   K <- Kseq[j]
   result <- compute(range=range, y=y.sub, d=d.sub, x=x.sub, subset=subset, 
                     K=K, nlam=nlam, n=n, p=p, const=F, ctrlvar=ctrl.sub)
   psi <- (d.sub * result$yfit + 
          (1-d.sub)*(y.sub-result$yfit)*result$ps/(1-result$ps))/pr #NOT influence fun
   theta0 <- mean(psi)   # point estimate
   tau    <- theta1 - theta0
   reg.table.ctr[2*j-1,3] <- format(round(tau, 3), nsmall=3)
   
   # influence fun.    
   varphi <- (d.sub*(y.sub-theta1)/pr) - 
             ((d.sub*(result$yfit-theta0)+
              (1-d.sub)*(y.sub-result$yfit)*result$ps/(1-result$ps))/pr)
   se     <- sqrt(mean(varphi^2))/sqrt(n)
   reg.table.ctr[2*j,3] <- paste("(", round(se, 3),  ")", sep="")
   #reg.table.ctr[2*j,3] <- paste("[", -round(se*1.96,3), ", ", round(se*1.96,3), "]", sep="")
}

lasso.fit <- rlassoATET(x=cbind(x.sub, ctrl.sub), d=d.sub, y=y.sub, bootstrap = "none", model = FALSE)
reg.table.ctr[7,3] <- format(round(lasso.fit$te, 3), nsmall=3)
reg.table.ctr[8,3] <- paste("(", round(lasso.fit$se, 3),  ")", sep="")


# base sample: drop Citi-related
ind <- (data$CorrCiti<=as.numeric(data$CorrCitiTr))
y.sub <- y.full[ind]; d.sub <- d.full[ind]; x.sub <- x.full[ind,]; ctrl.sub <- ctrlvar[ind,]
n <- length(y.sub); subset <- (d.sub==0); pr <- mean(d.sub)
theta1 <- mean(y.sub[d.sub==1]); range <- 1:n

for  (j in 1:3) {
   K <- Kseq[j]
   result <- compute(range=range, y=y.sub, d=d.sub, x=x.sub, subset=subset, 
                     K=K, nlam=nlam, n=n, p=p, const=F, ctrlvar=ctrl.sub)
   psi <- (d.sub * result$yfit + 
          (1-d.sub)*(y.sub-result$yfit)*result$ps/(1-result$ps))/pr #NOT influence fun
   theta0 <- mean(psi)   # point estimate
   tau    <- theta1 - theta0
   reg.table.ctr[2*j-1,4] <- format(round(tau, 3), nsmall=3)
   
   # influence fun.    
   varphi <- (d.sub*(y.sub-theta1)/pr) - 
             ((d.sub*(result$yfit-theta0)+
              (1-d.sub)*(y.sub-result$yfit)*result$ps/(1-result$ps))/pr)
   se     <- sqrt(mean(varphi^2))/sqrt(n)
   reg.table.ctr[2*j,4] <- paste("(", round(se, 3),  ")", sep="")
   #reg.table.ctr[2*j,4] <- paste("[", -round(se*1.96,3), ", ", round(se*1.96,3), "]", sep="")
}
lasso.fit <- rlassoATET(x=cbind(x.sub, ctrl.sub), d=d.sub, y=y.sub, bootstrap = "none", model = FALSE)
reg.table.ctr[7,4] <- format(round(lasso.fit$te, 3), nsmall=3)
reg.table.ctr[8,4] <- paste("(", round(lasso.fit$se, 3),  ")", sep="")


#######################################
# generate Table 2 in the main paper ##
#######################################
output <- cbind(reg.table[1:8,2], 
                reg.table[1:8,4],
                reg.table.ctr[,2], 
                reg.table.ctr[,4])

n.cgroup <- c(2, 2)
colheads <- c("Full Sample", "Base Sample", "Full Sample", "Base Sample")
cgroup   <- c("No Covariates", "Add Covariates")

n.rgroup <- c(6, 2)
rgroup   <- c("Local PCA, $K=$", "Double lasso")
rowname  <- c(Kseq[1], "", Kseq[2], "", Kseq[3], "",  "", "")
latex(output, file=paste("output/MainPaper_Table_ATT", ".txt", sep = ""), 
      append=FALSE, table.env=FALSE, center="none", title="", col.just=rep("c",4),
      n.cgroup=n.cgroup, cgroup=cgroup, colheads=colheads,
      n.rgroup=n.rgroup, rgroup=rgroup, rowname=rowname
      )

############################################################
# Uniform Inference
# Test stochastic dominance F(1|1) over F(0|1)
n <- length(y.full); subset <- (d.full==0); pr <- mean(d.full); 
range <- 1:n
eval <- sort(y.full); L <- length(eval)   # evaluation point
F.11 <- matrix(NA, L, 1); F.01 <- matrix(NA, L, 1); Fx <- matrix(NA, L, n)
score.11 <- matrix(NA, L, n); score.01 <- matrix(NA, L, n); score <- matrix(NA, L, n)

for (j in 1:L) {
   z <- eval[j]
   y.dist <- 1*(y.full <= z)
   result <- compute(range=range, y=y.dist, d=d.full, x=x.full, subset=subset, 
                     K=K0, nlam=nlam, n=n, p=p, const=F)
   # trim
   result$yfit[result$yfit > 1] <- 1
   result$yfit[result$yfit < 0] <- 0
   Fx[j,] <- result$yfit
   if (j>1) {
      index       <- (Fx[j,] < Fx[j-1,]) 
      Fx[j,index] <- Fx[j-1, index]
   }
   
   psi <- (d.full * Fx[j,] + 
          (1-d.full)*(y.dist-Fx[j,])*result$ps/(1-result$ps))/pr #NOT influence fun
   theta0 <- mean(psi)
   # trim
   if (theta0>1) {
      theta0 <- 1
   } else if (theta0<0) {
      theta0 <- 0
   }
   if (j>1) {
      if (theta0 < F.01[j-1,1]) {
          theta0 <- F.01[j-1,1]
      }
   } 
   F.01[j,1] <- theta0   # point estimate
   F.11[j,1] <- theta1 <- mean(y.dist[d.full==1])
   
   score.11[j,] <- (d.full*(y.dist-theta1)/pr)
   score.01[j,] <- ((d.full*(Fx[j,]-theta0) +
                    (1-d.full)*(y.dist-Fx[j,])*result$ps/(1-result$ps))/pr)
   score[j,]    <- score.11[j,] - score.01[j,]
}

## Test SD: Simulate the process
set.seed(12345)
sim <- 500; supre <- matrix(NA, sim, 1)
for (l in 1:sim) {
   eta <- rnorm(n)
   score.mul <- rowSums(sweep(score, 2, eta, "*"))/sqrt(n)
   supre[l]  <- max(score.mul)
}
cval <- quantile(supre, .95)
stat <- sqrt(n)*max(F.11-F.01)

# confidence band
sim <- 500; supre.11 <- supre.01 <- matrix(NA, sim, 1)

t1.11 <- min(which(F.11>=.1)); t9.11 <- max(which(F.11<=.9))
t1.01 <- min(which(F.01>=.1)); t9.01 <- max(which(F.01<=.9))

sd.11 <- sqrt(rowMeans(score.11[t1.11:t9.11,]^2))
sd.01 <- sqrt(rowMeans(score.01[t1.01:t9.01,]^2))
stu.score.11 <- sweep(score.11[t1.11:t9.11,], 1, sd.11, "/")
stu.score.01 <- sweep(score.01[t1.01:t9.01,], 1, sd.01, "/")

set.seed(12345)
for (l in 1:sim) {
   eta <- rnorm(n)
   score.mul.11 <- abs(rowSums(sweep(stu.score.11, 2, eta, "*"))/sqrt(n))
   score.mul.01 <- abs(rowSums(sweep(stu.score.01, 2, eta, "*"))/sqrt(n))
   supre.11[l]  <- max(score.mul.11)
   supre.01[l]  <- max(score.mul.01) 
}
cval.11 <- quantile(supre.11, .95)
cval.01 <- quantile(supre.01, .95)

F11.lb <- F11.ub <- F01.lb <- F01.ub <- matrix(NA, L, 1)
F11.lb[t1.11:t9.11,1] <- F.11[t1.11:t9.11,1] - cval.11*sd.11/sqrt(n)
F11.ub[t1.11:t9.11,1] <- F.11[t1.11:t9.11,1] + cval.11*sd.11/sqrt(n)
F01.lb[t1.01:t9.01,1] <- F.01[t1.01:t9.01,1] - cval.01*sd.01/sqrt(n)
F01.ub[t1.01:t9.01,1] <- F.01[t1.01:t9.01,1] + cval.01*sd.01/sqrt(n)
plot.dist <- data.frame(x=eval, F11=F.11, F01=F.01, 
                        F11.lb=F11.lb, F11.ub=F11.ub,
                        F01.lb=F01.lb, F01.ub=F01.ub)
plot.dist <- plot.dist[order(plot.dist$x),]
plot <- ggplot(data=plot.dist)
plot <- plot + geom_line(aes(x=x, y=F11, colour="Y(1)|D=1", linetype="Y(1)|D=1")) + 
               geom_line(aes(x=x, y=F01, colour="Y(0)|D=1", linetype="Y(0)|D=1")) +
               geom_ribbon(aes(x=x, ymin=F11.lb, ymax=F11.ub), fill="grey70", alpha=.2) +
               geom_ribbon(aes(x=x, ymin=F01.lb, ymax=F01.ub), fill="blue", alpha=.2) +
               scale_color_manual(name = 'CDF', 
                                  values = c('Y(1)|D=1' = 'black', 'Y(0)|D=1' = 'blue')) +
               scale_linetype_manual(name = 'CDF', 
                                     values = c('Y(1)|D=1' = 'solid', 'Y(0)|D=1' = 'dashed')) +
               xlab("") + ylab("") +theme_bw() + xlim(-0.2, 0.6) +
               theme(legend.position = c(.85, 0.15),
                     legend.background = element_rect(fill="transparent"),
                     panel.grid.major=element_blank(), panel.grid.minor = element_blank())
ggsave("output/SD.pdf", width = 6, height=3)


