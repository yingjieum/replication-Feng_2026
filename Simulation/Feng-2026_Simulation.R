###################################################
#### This file replicates the simulation results ##
############## in Feng (2026)  ####################
########### Author: Yingjie Feng ##################
########## Last updated: 01/29/2026 ################

rm(list = ls())
library(Rfast)
library(foreach)
library(doParallel)
library(doRNG)
library(irlba)
library(RcppNumerical)
library(pracma)
library(hdm)
source("Feng-2026_SimFuns.R")

# param
rep <- 2000
par <- read.csv("Feng-2026_SimModels.csv", header = T, colClasses=c(rep("numeric", 5), "logical", rep("numeric", 3)))

for (j in 1:44) {
n       <- par$n[j]
p       <- par$p[j]
model   <- funlist[[par$hdmodel[j]]]
nlam    <- par$nlam[j]
const   <- par$const[j]
err     <- par$err[j]
rho     <- par$rho[j]
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

#ptm <- proc.time()
## simulation
cl <- makeCluster(23)
registerDoParallel(cl)

output <- foreach (i = 1:rep, .options.RNG=1234, .packages=c('Rfast','irlba','RcppNumerical','hdm'),
                   .combine=rbind) %dorng% {
                     output <- sim(i, n, p, model, Kseq, nlam, theta0, const, err, rho, dim.alp)
                     output   # (6*rep) by (length(Kseq)+1+1+1) matrix
                   }

stopCluster(cl)

#proc.time() - ptm
###################
write.table(output, paste("output/rawoutput_par", j, "txt", sep = "."), sep = ",", 
            row.names = F, col.names = F)

}



###### Only used to generate Figure 1 in the paper #############################
theta0 <- 1.539319 # dim 1
#theta0 <- 0.9459911 # dim 2
rep <- 2000

model.f1 <- function(u, v) {
  return(outer(u, rep(1, length(v)), "*"))
}


cl <- makeCluster(23)
registerDoParallel(cl)

output <- foreach (i = 1:rep, .options.RNG=1234, .packages=c('Rfast','RcppNumerical','hdm'),
                   .combine=rbind) %dorng% {
                    output <- sim.selonx(i, n=500, p=250, model=model.f1, theta0=theta0, err=0, rho=0, dim.alp=1, nx=10)
                    output   # (5*rep) by (1+1) matrix
                   }

stopCluster(cl)

rownum <- rep(1:5, rep)
est.dr   <- output[rownum==2,1] 
se.dr    <- output[rownum==3,1]
t.dr     <- (est.dr-theta0)/se.dr
est.la   <- output[rownum==2,2] 
se.la    <- output[rownum==3,2]
t.la     <- (est.la-theta0)/se.la

df <- rbind(
  data.frame(value = t.dr, series = "10 proxies"),
  data.frame(value = t.la, series = "double lasso")
)

hist <- ggplot(df, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30, fill = "grey80", color = "white") +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ series, nrow = 1, strip.position = "bottom") +
  theme_classic() +
  theme(
    strip.placement = "outside",
    strip.background = element_blank()
  ) +
  labs(x = NULL, y = "Density")

ggsave("hist_dr_la.pdf", plot = hist, width = 8, height = 3, dpi = 300)
