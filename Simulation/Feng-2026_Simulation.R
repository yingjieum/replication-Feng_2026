# Replicate all simulation results reported in Feng (2026).
# Run this file from the Simulation directory.

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

dir.create("raw-output", showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)

n.rep <- 2000
available.cores <- parallel::detectCores()
n.workers <- if (is.na(available.cores)) 1L else max(1L, available.cores - 1L)
designs <- read.csv(
  "Feng-2026_SimModels.csv",
  header = TRUE,
  colClasses = c(rep("numeric", 5), "logical", rep("numeric", 3))
)

# Population ATT for designs with one or two latent confounders.
true.att <- function(dim.alp) {
  if (dim.alp == 1) {
    p1 <- integrate(px, 0, 1)$value
    integrand <- function(alpha) {
      (mu1(alpha) - mu0(alpha)) * px(alpha)
    }
    theta0 <- integrate(integrand, 0, 1)$value / p1
  } else {
    p1 <- integral2(px.m, 0, 1, 0, 1)$Q
    integrand <- function(alpha1, alpha2) {
      (1 + 0.6 * (alpha1 - 0.5) + 0.4 * sin(2 * pi * alpha2) +
         0.3 * (alpha1 - 0.5) * (alpha2 - 0.5)) * px.m(alpha1, alpha2)
    }
    theta0 <- integral2(integrand, 0, 1, 0, 1)$Q / p1
  }
  theta0
}

# Generate the raw results for the 44 designs in Feng-2026_SimModels.csv.
for (j in 1:44) {
  n       <- designs$n[j]
  p       <- designs$p[j]
  model   <- funlist[[designs$hdmodel[j]]]
  nlam    <- designs$nlam[j]
  const   <- designs$const[j]
  err     <- designs$err[j]
  rho     <- designs$rho[j]
  dim.alp <- designs$dim.alp[j]

  K <- n^(4 / (4 + dim.alp))
  Kseq <- ceiling(K * c(0.5, 1, 1.5, 2))
  theta0 <- true.att(dim.alp)

  cl <- makeCluster(n.workers)
  registerDoParallel(cl)

  sim.output <- foreach(
    i = 1:n.rep,
    .options.RNG = 1234,
    .packages = c("Rfast", "irlba", "RcppNumerical", "hdm"),
    .combine = rbind
  ) %dorng% {
    sim(n, p, model, Kseq, nlam, theta0, const, err, rho, dim.alp)
  }

  stopCluster(cl)

  write.table(
    sim.output,
    paste("raw-output/rawoutput_par", j, "txt", sep = "."),
    sep = ",",
    row.names = FALSE,
    col.names = FALSE
  )
}

# Figure 1: t statistics for the three ATT estimators.
library(ggplot2)
theta0 <- 1.539319

model.f1 <- function(u, v) {
  outer(u, rep(1, length(v)), "*")
}

cl <- makeCluster(n.workers)
registerDoParallel(cl)

figure.output <- foreach(
  i = 1:n.rep,
  .options.RNG = 1234,
  .packages = c("Rfast", "irlba", "RcppNumerical", "hdm"),
  .combine = rbind
) %dorng% {
  sim.selonx(
    n = 500, p = 250, model = model.f1, theta0 = theta0,
    err = 0, rho = 0, dim.alp = 1, nx = 10
  )
}

stopCluster(cl)

rownum <- rep(1:5, n.rep)
est.dr <- figure.output[rownum == 2, 1]
se.dr <- figure.output[rownum == 3, 1]
t.dr <- (est.dr - theta0) / se.dr
est.la <- figure.output[rownum == 2, 2]
se.la <- figure.output[rownum == 3, 2]
t.la <- (est.la - theta0) / se.la
est.lpca <- figure.output[rownum == 2, 3]
se.lpca <- figure.output[rownum == 3, 3]
t.lpca <- (est.lpca - theta0) / se.lpca

plot.data <- rbind(
  data.frame(value = t.dr, series = "10 proxies"),
  data.frame(value = t.la, series = "double lasso"),
  data.frame(value = t.lpca, series = "Proposed approach")
)
plot.data$series <- factor(
  plot.data$series,
  levels = c("10 proxies", "double lasso", "Proposed approach")
)

figure.plot <- ggplot(plot.data, aes(x = value)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30, fill = "grey80", color = "white"
  ) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1), linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ series, nrow = 1, strip.position = "bottom") +
  theme_classic(base_size = 12) +
  theme(
    strip.placement = "outside",
    strip.background = element_blank()
  ) +
  labs(x = NULL, y = "Density")

ggsave(
  "output/hist_dr_la_lpca.pdf",
  plot = figure.plot,
  width = 9,
  height = 3,
  dpi = 300
)

# Summarize a raw simulation file as RMSE, coverage rate, and interval length.
summarize.output <- function(j, theta0) {
  sim.output <- as.matrix(read.table(
    paste("raw-output/rawoutput_par", j, "txt", sep = "."),
    sep = ","
  ))
  rownum <- rep(1:6, n.rep)
  est <- sim.output[rownum == 2, ]
  rej <- sim.output[rownum == 4, ]
  ci <- sim.output[rownum == 5, ]

  rmse <- sqrt(colMeans((est - theta0)^2))
  CR <- 1 - colMeans(rej)
  AL <- colMeans(ci)
  round(rbind(rmse, CR, AL), 3)
}

library(Hmisc)

# Table 2 in the main paper.
mat <- matrix(NA, 9, 7)
for (j in c(2, 31, 32)) {
  n <- designs$n[j]
  dim.alp <- designs$dim.alp[j]
  K <- n^(4 / (4 + dim.alp))
  Kseq <- ceiling(K * c(0.5, 1, 1.5, 2))
  table.values <- summarize.output(j, true.att(dim.alp))

  if (j == 31) {
    mat[1:3, ] <- table.values
  } else if (j == 2) {
    mat[4:6, ] <- table.values
  } else if (j == 32) {
    mat[7:9, ] <- table.values
  }
}

n.rgroup <- c(3, 3, 3)
rowname <- rep(c("RMSE", "CR", "AL"), 3)
rgroup <- c("low-rank", "mid-rank", "high-rank")
n.cgroup <- c(5, 1, 1)
cgroup <- c("Local PCA, $K=$", "DR", "DR-DL")
colheads <- c(paste(Kseq), "$\\widehat{K}_{\\mathtt{CV}}$", "", "")
latex(
  mat,
  file = "output/Table_Pointwise_MainPaper.txt",
  col.just = rep("c", 7),
  append = FALSE,
  table.env = FALSE,
  center = "none",
  title = "",
  n.cgroup = n.cgroup,
  cgroup = cgroup,
  colheads = colheads,
  n.rgroup = n.rgroup,
  rgroup = rgroup,
  rowname = rowname
)

# Tables SA-9 through SA-17 in the supplement.
model1.rows <- c(13, 1, 7, 33)
iid.rows <- 19:24
ar1.rows <- 25:30
local.constant.rows <- 39:44
supplement.rows <- list(
  model1.rows,
  model1.rows + 1,
  model1.rows + 2,
  model1.rows + 3,
  model1.rows + 4,
  model1.rows + 5,
  iid.rows,
  ar1.rows,
  local.constant.rows
)

for (i in 1:9) {
  if (i <= 6) {
    mat <- matrix(NA, 12, 7)
  } else {
    mat <- NULL
    Kseq.by.model <- list()
  }

  for (j in supplement.rows[[i]]) {
    n <- designs$n[j]
    p <- designs$p[j]
    dim.alp <- designs$dim.alp[j]
    K <- n^(4 / (4 + dim.alp))
    Kseq <- ceiling(K * c(0.5, 1, 1.5, 2))
    if (i > 6) Kseq.by.model[[length(Kseq.by.model) + 1]] <- Kseq

    table.values <- summarize.output(j, true.att(dim.alp))

    if (i <= 6) {
      if (p == 125) {
        mat[1:3, ] <- table.values
      } else if (p == 250) {
        mat[4:6, ] <- table.values
      } else if (p == 500) {
        mat[7:9, ] <- table.values
      } else if (p == 750) {
        mat[10:12, ] <- table.values
      }
    } else {
      mat <- rbind(mat, table.values)
    }
  }

  if (i <= 6) {
    n.rgroup <- c(3, 3, 3, 3)
    rowname <- rep(c("RMSE", "CR", "AL"), 4)
    rgroup <- c(
      "$(n,p)=(500, 125)$", "$(n,p)=(500, 250)$",
      "$(n,p)=(500, 500)$", "$(n,p)=(500, 750)$"
    )
    n.cgroup <- c(5, 1, 1)
    cgroup <- c("Local PCA, $K=$", "DR", "DR-DL")
    colheads <- c(paste(Kseq), "$\\widehat{K}_{\\mathtt{CV}}$", "", "")
  } else {
    n.rgroup <- rep(3, 6)
    rowname <- rep(c("RMSE", "CR", "AL"), 6)
    rgroup <- paste("Model", 1:6)

    stopifnot(
      length(Kseq.by.model) == 6,
      all(vapply(Kseq.by.model[1:3], identical, logical(1), Kseq.by.model[[1]])),
      all(vapply(Kseq.by.model[4:6], identical, logical(1), Kseq.by.model[[4]])),
      identical(Kseq.by.model[[1]], c(73, 145, 217, 289)),
      identical(Kseq.by.model[[4]], c(32, 63, 95, 126))
    )
    n.cgroup <- c(5, 1, 1)
    method <- if (i == 9) "Local constant" else "Local PCA"
    cgroup <- c(paste0(method, ", $K=$"), "DR", "DR-DL")
    colheads <- c(
      "$K_1$", "$K_2$", "$K_3$", "$K_4$",
      "$\\widehat{K}_{\\mathtt{CV}}$", "", ""
    )
  }

  latex(
    mat,
    file = paste("output/Table_Pointwise_SA_", i, ".txt", sep = ""),
    col.just = rep("c", 7),
    append = FALSE,
    table.env = FALSE,
    center = "none",
    title = "",
    n.cgroup = n.cgroup,
    cgroup = cgroup,
    colheads = colheads,
    n.rgroup = n.rgroup,
    rgroup = rgroup,
    rowname = rowname
  )
}
