#########################################################
##### This file produces Tables in the main paper#######
########################################################
library(Hmisc)
library(Rfast)
source("Feng-2026_SimFuns.R")
par <- read.csv("Feng-2026_SimModels.csv", header = T, colClasses=c(rep("numeric", 5), "logical", rep("numeric", 3)))

####################################
# Table used in the main paper
####################################
mat <- matrix(NA, 9, 7)
#sum.K <- NULL

for (j in c(2, 31, 32)) {
n       <- par$n[j]
p       <- par$p[j]
model   <- funlist[[par$hdmodel[j]]]
nlam    <- par$nlam[j]
dim.alp <- par$dim.alp[j]
K       <- n^(4/(4+dim.alp))
Kseq    <- ceiling(K * c(0.5, 1, 1.5, 2))
###############################################
# calculate true value
if (dim.alp==1) {
  p1 <- integrate(px, 0, 1)$value
  integrand <- function(alpha) {
    return((mu1(alpha) - mu0(alpha)) * px(alpha))     # True ATT
  }
  theta0 <- (integrate(integrand, 0, 1)$value)/p1
} else {
  p1 <- integral2(px.m, 0, 1, 0, 1)$Q
  integrand <- function(alpha1, alpha2) {
    return((1 + 0.6*(alpha1 - 0.5) + 0.4*sin(2*pi*alpha2) + 0.3*(alpha1 - 0.5)*(alpha2 - 0.5))*px.m(alpha1, alpha2))     # True ATT
  }
  theta0 <- (integral2(integrand, 0, 1, 0, 1)$Q)/p1
}
#####################################################

output <- as.matrix(read.table(paste("output/rawoutput_par", j ,"txt", sep="."), sep = ","))

rep <- 2000
rownum <- rep(1:6, rep)

est   <- output[rownum==2,] 
se    <- output[rownum==3,]
rej   <- output[rownum==4,]
ci    <- output[rownum==5,]
k.cv  <- output[rownum==1,5]


bias  <- colMeans(est) - theta0
sd    <- sqrt(colVars(est))
rmse  <- sqrt(colMeans((est - theta0)^2))
CR    <- 1-colMeans(rej)
AL    <- colMeans(ci)
#table <- round(rbind(bias, sd, rmse, CR, AL), 3)
table <- round(rbind(rmse, CR, AL), 3)

sum.K <- rbind(sum.K, c(summary(k.cv), sd(k.cv)))

if (j==31) {         # sig=0.5
  mat[1:3,] <- table
} else if (j==2) {   # sig=0.1
  mat[4:6,] <- table
} else if (j==32) {  # sig=0.02
  mat[7:9,] <- table
}



# if (j==1) {
#    mat[6:10, 1:4] <- table  
# } else if (j==2) {
#    mat[6:10, 5:8] <- table
# } else if (j==3) {
#    mat[11:15, 1:4] <- table
# } else if (j==4) {
#    mat[11:15, 5:8] <- table
# } else if (j==5) {
#    mat[1:5, 1:4] <- table
# } else if (j==6) {
#    mat[1:5, 5:8] <- table
# }

}

# Report main simulation results
n.rgroup <- c(3, 3, 3)
rowname  <- rep(c("RMSE", "CR", "AL"), 3)
rgroup   <- c("low-rank", "mid-rank", "high-rank")

n.cgroup <- c(5, 1, 1)
cgroup   <- c("Local PCA, $K=$", "DR", "DR-DL")
colheads <- c(paste(Kseq), "$\\widehat{K}_{\\mathtt{CV}}$", "", "")
latex(mat, file=paste("Table_Pointwise_MainPaper", ".txt", sep = ""),
      col.just = rep("c", 7),
      append=FALSE, table.env=FALSE, center="none", title="",
      n.cgroup=n.cgroup, cgroup=cgroup, colheads=colheads,
      n.rgroup=n.rgroup, rgroup=rgroup, rowname=rowname
      )




########################################################
# Report tables in the Supplemental Appendix ###########
########################################################
m1   <- c(13, 1, 7, 33)
hom  <- 19:24
ar1  <- 25:30
lcon <- 39:44 
m.ls <- list(m1, m1+1, m1+2, m1+3, m1+4, m1+5, hom, ar1, lcon)

for (i in 1:9) {
  if (i<=6) {
    mat <- matrix(NA, 12, 7)
  } else {
    mat <- NULL    # finally, it should be matrix(, 18, 7)
  }
  
for (j in m.ls[[i]]) {
  n       <- par$n[j]
  p       <- par$p[j]
  model   <- funlist[[par$hdmodel[j]]]
  nlam    <- par$nlam[j]
  dim.alp <- par$dim.alp[j]
  K       <- n^(4/(4+dim.alp))
  Kseq    <- ceiling(K * c(0.5, 1, 1.5, 2))
  ###############################################
  # calculate true value
  if (dim.alp==1) {
    p1 <- integrate(px, 0, 1)$value
    integrand <- function(alpha) {
      return((mu1(alpha) - mu0(alpha)) * px(alpha))     # True ATT
    }
    theta0 <- (integrate(integrand, 0, 1)$value)/p1
  } else {
    p1 <- integral2(px.m, 0, 1, 0, 1)$Q
    integrand <- function(alpha1, alpha2) {
      return((1 + 0.6*(alpha1 - 0.5) + 0.4*sin(2*pi*alpha2) + 0.3*(alpha1 - 0.5)*(alpha2 - 0.5))*px.m(alpha1, alpha2))     # True ATT
    }
    theta0 <- (integral2(integrand, 0, 1, 0, 1)$Q)/p1
  }
  #####################################################
  
  output <- as.matrix(read.table(paste("output/rawoutput_par", j ,"txt", sep="."), sep = ","))
  
  rep <- 2000
  rownum <- rep(1:6, rep)
  
  est   <- output[rownum==2,] 
  se    <- output[rownum==3,]
  rej   <- output[rownum==4,]
  ci    <- output[rownum==5,]
  k.cv  <- output[rownum==1,5]
  
  
  bias  <- colMeans(est) - theta0
  sd    <- sqrt(colVars(est))
  rmse  <- sqrt(colMeans((est - theta0)^2))
  CR    <- 1-colMeans(rej)
  AL    <- colMeans(ci)
  #table <- round(rbind(bias, sd, rmse, CR, AL), 3)
  table <- round(rbind(rmse, CR, AL), 3)
  
  if (i<=6) {
    if (p==125) {         
      mat[1:3,]   <- table
    } else if (p==250) {   
      mat[4:6,]   <- table
    } else if (p==500) {  
      mat[7:9,]   <- table
    } else if (p==750) {
      mat[10:12,] <- table
    }
  } else {
    mat <- rbind(mat, table)
  }
  
}

# Report main simulation results
if (i<=6) {
  n.rgroup <- c(3, 3, 3, 3)
  rowname <- rep(c("RMSE", "CR", "AL"), 4)
  rgroup   <- c("$(n,p)=(500, 125)$", "$(n,p)=(500, 250)$", "$(n,p)=(500, 500)$", "$(n,p)=(500, 750)$")

  n.cgroup <- c(5, 1, 1)
  cgroup   <- c("Local PCA, $K=$", "DR", "DR-DL")
  colheads <- c(paste(Kseq), "$\\widehat{K}_{\\mathtt{CV}}$", "", "")
  latex(mat, file=paste("Table_Pointwise_SA_", i, ".txt", sep = ""), 
        col.just = rep("c", 7),
        append=FALSE, table.env=FALSE, center="none", title="",
        n.cgroup=n.cgroup, cgroup=cgroup, colheads=colheads,
        n.rgroup=n.rgroup, rgroup=rgroup, rowname=rowname
  )
} else {
  n.rgroup <- rep(3, 6)
  rowname <- rep(c("RMSE", "CR", "AL"), 6)
  rgroup   <- c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6")
  
  n.cgroup <- c(5, 1, 1)
  if (i==9) {
    cgroup   <- c("Local constant, $K=$", "DR", "DR-DL")
  } else {
    cgroup   <- c("Local PCA, $K=$", "DR", "DR-DL")
  }
  colheads <- c(paste(Kseq), "$\\widehat{K}_{\\mathtt{CV}}$", "", "")
  latex(mat, file=paste("Table_Pointwise_SA_", i, ".txt", sep = ""), 
        col.just = rep("c", 7),
        append=FALSE, table.env=FALSE, center="none", title="",
        n.cgroup=n.cgroup, cgroup=cgroup, colheads=colheads,
        n.rgroup=n.rgroup, rgroup=rgroup, rowname=rowname
  )
}
  
}


########## Generate figure ############
rep <- 2000
rownum <- rep(1:5, rep)
theta0 <- 1.539319
CR <- c()
for (j in 1:30) {
  output <- as.matrix(read.table(paste("output/selectX/selonx_output_nx", j ,"txt", sep="."), sep = ","))
  rej   <- output[rownum==4,]
  CR    <- c(CR, 1-colMeans(rej))
}
######

