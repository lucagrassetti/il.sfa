## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ilsfa)

## -----------------------------------------------------------------------------
N <- 100
T <- 5
para <- c(1,-1.5,1,0.25)
dat <- DGP(N, T, list(beta=para[1], sigmav=para[4], gamma=para[2:3]), model = "Exp",
           seed.val = 22, nsim = 10)

## -----------------------------------------------------------------------------
mle <- il_sfa(X = dat$x, y = dat$daty[1,], distr = "Exp", het = TRUE, z = dat$datz,
              group = dat$g, useHess = TRUE)
mle

## -----------------------------------------------------------------------------
mle <- il_sfa(X = dat$x, y = dat$daty[1,], distr = "HN", het = TRUE, z = dat$datz,
              group = dat$g, useHess = TRUE)

mle

## -----------------------------------------------------------------------------
N <- 100
T <- 5
para <- c(1,-1.5,1,0.25)
nsimul <- 10
dat <- DGP(N, T, list(beta=para[1], sigmav=para[4], gamma=para[2:3]), model = "Exp",
           seed.val = 22, nsim = nsimul)

## -----------------------------------------------------------------------------
res <- matrix(0, nrow = nsimul, ncol=4)
for(i in 1:nrow(res)){
#  print(i)
  mle <- il_sfa(X = dat$x, y = dat$daty[i,], distr = "Exp", het = TRUE, z = dat$datz,
                group = dat$g, useHess = TRUE)
  res[i,] <- mle$par
  }

## -----------------------------------------------------------------------------
for(i in 1:4) print(MSE(res[,i], para[i]))

## -----------------------------------------------------------------------------
data(res_N100_Stata)

for(i in 1:4) print(MSE(res_N100_Stata[1:nsimul,i], para[i]))

