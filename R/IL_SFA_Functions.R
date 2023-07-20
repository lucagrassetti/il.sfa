# Da fare:
# metti in ordine anche le altre funzioni oltre a exp per il calcolo di tutto
# sistema gli essiani creando un paio di funzioni farlocche peril solo calcolo degli essiani.

##### Functions for IL estimation

MSE <- function(x, para) mean( (x-para)^2)

# DGP <- function(n, m, param, model = "HN", seed.val = 0, nsim = 1000, zvar = FALSE)
# {
#     g <- sort(rep(1:n, m))
#     daty <- datu <- matrix(NA, nrow = nsim, ncol = n * m)
#     set.seed(seed.val)
#     alpha <- rnorm(n)
#     w <- rnorm(n * m)
#     beta <- param$beta
#     s2v <- param$sigmav^2
#     x <- 0.5 * alpha[g] + sqrt(0.5^2) * w
#     z <- if (zvar) rnorm(n * m, sqrt(0.0625)) else rnorm(n, sqrt(0.0625))[g]
#     for(i in 1:nsim) {
#         if(model == "HN")  {
#           s2u <- param$sigmau^2
#           sigma <- sqrt(s2u + s2v)
#           lambda <- param$sigmau / param$sigmav
#           u <- abs(rnorm(n * m, s = sqrt(s2u)))
#           }
#         if(model == "Exp")  {
#            su <- exp(param$gamma[1] +  param$gamma[2] * z)
#            u <- rexp(n *m, rate = 1 / su)
#         }
#         daty[i,] <-  alpha[g] + x * beta + rnorm(n * m, s = sqrt(s2v)) - u
#         datu[i,] <- u
#        }
#     return(list(daty = daty, datu = datu, x = as.matrix(x, ncol = 1), g = g, datz = z))
# }
#

##########################################################################
#    HN case
##########################################################################
likHNall <- function(para, X, group, list_y, alphainit, niter, ws, nodes){
  p <- ncol(X)
  eta <- as.vector(as.matrix(X) %*% para[1:p])
  list_eta <- split(eta, group)
  ll <- likHN(para[p + 1], para[p + 2],  list_eta, list_y, niter, alphainit, ws, nodes)
  return(-ll)
}

likHNall_orig <- function(para, X, group, list_y, alphainit, niter, ws, nodes){
  p <- ncol(X)
  eta <- as.vector(as.matrix(X) %*% para[1:p])
  list_eta <- split(eta, group)
  log_sig <- log(sqrt(para[p + 1]^2+para[p+2]^2))
  log_lam <- log(para[p+1]/para[p+2])
  ll <- likHN(log_lam, log_sig,  list_eta, list_y, niter, alphainit, ws, nodes)
  return(-ll)
}

## ucminf
estim_HN <- function(X, y, group, ols, nq = 25, niter = 10,
                         init = NULL, useHess = TRUE,
                        trace = 0, grtol = 10^-6,
                         grad = "central"){
  list_y <- split(y, group)
  p <- ncol(X)
  obj.gh <- statmod::gauss.quad(nq, "hermite")
  ws <- obj.gh$weights * exp(obj.gh$nodes^2)
  if(is.null(init)) init <- c(ols$beta, 1, 0) # changed initial values
  alpha <- ols[1:(length(ols)-p)]
  H <- if(useHess) {numDeriv::hessian(likHNall, init,  X = X, list_y = list_y,  group = group,
                                     alphainit = ols$alpha, niter = niter,
                                     ws = ws, nodes = obj.gh$nodes)} else diag(length(init))
  L1 <- solve(H)[lower.tri(H, diag=TRUE)]
  ogg <- ucminf::ucminf(init, likHNall, X = X, list_y = list_y,  group = group,
                        alphainit = ols$alpha, niter = niter,
                        ws = ws, nodes = obj.gh$nodes,
                        control = list(grad = grad, trace = trace, invhessian.lt = L1, grtol = grtol))
  out  <- list(par = getsig(ogg$par))
  names(out$par) <- c(names(X), "sigmau", "sigmav")
  hes <- numDeriv::hessian(likHNall_orig, out$par,  X = X, list_y = list_y,  group = group,
                           alphainit = ols$alpha, niter = niter,
                           ws = ws, nodes = obj.gh$nodes)
  colnames(hes) <- rownames(hes) <- c(names(X), "sigmau", "sigmav")
  out$invhes <- solve(hes)
  out$se <- sqrt(diag(out$invhes))
  names(out$se) <- c(names(X), "sigmau", "sigmav")
  out$alphai <- alphaiHN(para = ogg$par, X = X, group = group, list_y = list_y,
                           niter=niter, uinit = rep(0, length(y)))
  out$uit <- extractUit_HN(para = out$par, X = X, y = y, group = group, alphai = out$alphai)
  return(out)
}

## conditional expected value of ui a la Jondrow et al. (1982)

extractUit_HN <- function(para, X, y, group, alphai)
{
  p <- ncol(X)
  sigmau <- para[p+1]
  sigmav <- para[p+2]
  sig2 <- sigmau^2+sigmau^2
  res <- y - (alphai[group] + as.vector(as.matrix(X) %*% para[1:p]))
  mu_star <- -(sigmau^2*res)/sig2
  sig2_star <- ((sigmau*sigmav)^2/sig2)

  u_hat <- mu_star + sqrt(sig2_star) * (dnorm(-mu_star/sqrt(sig2_star))/(1 - pnorm(-mu_star/sqrt(sig2_star))))
}


##beta lambda sigma
getsig <- function(para)
{
  p <- length(para) - 2
  out <- para
  lambda <- exp(para[p + 1])
  sigma <- exp(para[p + 2])
  sigmav <- sigma / sqrt(1 + lambda^2)
  sigmau <- lambda * sigmav
  out[p + 1] <- sigmau
  out[p + 2] <- sigmav
  return(out)
}

alphaiHN <- function(para, X, group, list_y, niter, uinit){
  p <- ncol(X)
  lnlambda <- para[p+1]
  lnsigma <- para[p+2]
  eta <- as.vector(as.matrix(X) %*% para[1:p])
  list_eta <- split(eta, group)
  alphai <- getEffectsHN(para[p+1], para[p+2], list_eta, list_y,
                         niter, uinit)
  return(alphai)
}

#############################################################
#    EXP case
#############################################################
likEHETall <- function(para, X, group, list_y, list_z, alphainit, niter, ws, nodes){
  p <- ncol(X)
  eta <- as.vector(as.matrix(X) %*% para[1:p])
  list_eta <- split(eta, group)
  ll <- likEHET(para[p + (1:2)], para[p + 3],  list_eta, list_y, list_z, niter, alphainit, ws, nodes)
  return(-ll)
}

likEHETall_orig <- function(para, X, group, list_y, list_z, alphainit, niter, ws, nodes){
  p <- ncol(X)
  eta <- as.vector(as.matrix(X) %*% para[1:p])
  list_eta <- split(eta, group)
  ll <- likEHET(para[p + (1:2)], log(para[p + 3]),  list_eta, list_y, list_z, niter, alphainit, ws, nodes)
  return(-ll)
}


likEXPall <- function(para, X, group, list_y, alphainit, niter, ws, nodes){
  p <- ncol(X)
  eta <- as.vector(as.matrix(X) %*% para[1:p])
  list_eta <- split(eta, group)
  ll <- likEHET(c(para[p + 1],0), para[p + 2],  list_eta, list_y, list_y, niter, alphainit, ws, nodes)
  return(-ll)
}

likEXPall_orig <- function(para, X, group, list_y, alphainit, niter, ws, nodes){
  p <- ncol(X)
  eta <- as.vector(as.matrix(X) %*% para[1:p])
  list_eta <- split(eta, group)
  ll <- likEHET(c(para[p + 1],0), log(para[p + 2]),  list_eta, list_y, list_y, niter, alphainit, ws, nodes)
  return(-ll)
}

alphaiEHET <- function(para, X, group, list_y, list_z, niter, uinit){
  p <- ncol(X)
  gamma <- para[p+1:2]
  lnsigmav <- para[p+3]
  eta <- as.vector(as.matrix(X) %*% para[1:p])
  list_eta <- split(eta, group)
  alphai <- getEffectsEHET(para[p+1:2], para[p+3], list_eta, list_y,
                       list_z, niter, uinit)
  return(alphai)
}

## ucminf
estim_EHET <- function(X, y, z, group, ols, nq = 25, niter = 10,
                      trace = 0, grtol = 10^-6,
                      init = NULL, useHess = TRUE,
                      grad = "central"){
  list_y <- split(y, group)
  list_z <- split(z, group)
  p <- ncol(X)
  obj.gh <- statmod::gauss.quad(nq, "hermite")
  ws <- obj.gh$weights * exp(obj.gh$nodes^2)
  if(is.null(init)) init <- c(ols$beta, 0, 0, 0)
  H <- if(useHess) numDeriv::hessian(likEHETall, init,  X = X, list_y = list_y, list_z = list_z, group = group,
                                     alphainit = ols$alpha, niter = niter,
                                     ws = ws, nodes = obj.gh$nodes)
       else diag(length(init))
  L1 <- solve(H)[lower.tri(H, diag=TRUE)]
  ogg <- ucminf::ucminf(init, likEHETall, X = X, list_y = list_y,  list_z = list_z, group = group,
                        alphainit = ols$alpha, niter = niter,
                        ws = ws, nodes = obj.gh$nodes,
                        control = list(grad = grad, trace = trace, invhessian.lt = L1, grtol = grtol))
  out  <- list(par = ogg$par)
  names(out$par) <- c(names(X), "gamma1", "gamma2", "sigma_v")
  out$par[p + 3] <- exp(out$par[p + 3])
  hes <- numDeriv::hessian(likEHETall_orig, out$par,  X = X, list_y = list_y, list_z = list_z, group = group,
                           alphainit = ols$alpha, niter = niter,
                           ws = ws, nodes = obj.gh$nodes)
  colnames(hes) <- rownames(hes) <- c(names(X), "gamma1", "gamma2", "sigma_v")
  out$invhes <- solve(hes)
  out$se <- sqrt(diag(out$invhes))
#  out$se[p + 3] <- sqrt(exp(out$par[p + 3]*2) * out$se[p + 3]^2)
  names(out$se) <- c(names(X), "gamma1", "gamma2", "sigma_v")
  out$alphai <- alphaiEHET(para = ogg$par, X = X, group = group, list_y = list_y,
                    list_z = list_z, niter=niter, uinit = rep(0, length(y)))
  out$uit <- extractUit_EHET(para = out$par, X = X, y = y, z = z, group = group, alphai = out$alphai)
  return(out)
}

## conditional expected value of ui a la Jondrow et al. (1982)

extractUit_EHET <- function(para, X, y, z, group, alphai)
{
  p <- ncol(X)
  gamma <- para[p+(1:2)]
  sigmav <- para[p+3]

  lambda_it <- gamma[1] + gamma[2] * z

  res <- y - (alphai[group] + as.vector(as.matrix(X) %*% para[1:p]))

  A <- res/sigmav + sigmav*lambda_it

  u_hat <- sigmav * (dnorm(A)/(1 - pnorm(A)) - A)
}

## ucminf
estim_EXP <- function(X, y, group, ols, nq = 25, niter = 10,
                     init = NULL, useHess = TRUE,
                     trace = 0, grtol = 10^-6,
                     grad = "central"){
  list_y <- split(y, group)
  list_z <- split(rep(1,length(y)), group)
  p <- ncol(X)
  obj.gh <- statmod::gauss.quad(nq, "hermite")
  ws <- obj.gh$weights * exp(obj.gh$nodes^2)
  if(is.null(init)) init <- c(ols$beta, 0, 0)
  H <- if(useHess) numDeriv::hessian(likEXPall, init,  X = X, list_y = list_y,  group = group,
                                     alphainit = ols$alpha, niter = niter,
                                     ws = ws, nodes = obj.gh$nodes)
  else diag(length(init))
  L1 <- solve(H)[lower.tri(H, diag=TRUE)]
  ogg <- ucminf::ucminf(init, likEXPall, X = X, list_y = list_y,  group = group,
                        alphainit = ols$alpha, niter = niter,
                        ws = ws, nodes = obj.gh$nodes,
                        control = list(grad = grad, trace = trace, invhessian.lt = L1, grtol = grtol))
  out  <- list(par = ogg$par)
  names(out$par) <- c(names(X), "gamma", "sigmav")
  out$par[p + 2] <- exp(out$par[p + 2])
  hes <- numDeriv::hessian(likEXPall_orig, out$par,  X = X, list_y = list_y,  group = group,
                           alphainit = ols$alpha, niter = niter,
                           ws = ws, nodes = obj.gh$nodes)
  colnames(hes) <- rownames(hes) <- c(names(X), "gamma", "sigmav")
  out$invhes = solve(hes)
  out$se <- sqrt(diag(out$invhes))
#  out$se[p + 2] <- sqrt(exp(out$par[p + 2]*2) * out$se[p + 2]^2)
  names(out$se) <- c(names(X), "gamma", "sigma_v")
  out$alphai <- alphaiEHET(para = c(ogg$par[1:p],ogg$par[p+1], 0,ogg$par[p+2]), X = X, group = group, list_y = list_y,
                           list_z = list_z, niter=niter, uinit = rep(0, length(y)))
  out$uit <- extractUit_EHET(para = c(out$par[1:p],out$par[p+1],0,out$par[p+2]), X = X, y = y, z =z, group = group, alphai = out$alphai)

  return(out)
}


