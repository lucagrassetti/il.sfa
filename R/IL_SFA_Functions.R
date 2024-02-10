MSE <- function(x, para) mean( (x-para)^2)

##### Functions for IL estimation

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
  names(out$par) <- c(colnames(X), "sigmau", "sigmav")
  hes <- numDeriv::hessian(likHNall_orig, out$par,  X = X, list_y = list_y,  group = group,
                           alphainit = ols$alpha, niter = niter,
                           ws = ws, nodes = obj.gh$nodes)
  colnames(hes) <- rownames(hes) <- c(colnames(X), "sigmau", "sigmav")
  out$invhes <- solve(hes)
  out$se <- sqrt(diag(out$invhes))
  names(out$se) <- c(colnames(X), "sigmau", "sigmav")
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
  names(out$par) <- c(colnames(X), "gamma1", "gamma2", "sigma_v")
  out$par[p + 3] <- exp(out$par[p + 3])
  hes <- numDeriv::hessian(likEHETall_orig, out$par,  X = X, list_y = list_y, list_z = list_z, group = group,
                           alphainit = ols$alpha, niter = niter,
                           ws = ws, nodes = obj.gh$nodes)
  colnames(hes) <- rownames(hes) <- c(colnames(X), "gamma1", "gamma2", "sigma_v")
  out$invhes <- solve(hes)
  out$se <- sqrt(diag(out$invhes))
#  out$se[p + 3] <- sqrt(exp(out$par[p + 3]*2) * out$se[p + 3]^2)
  names(out$se) <- c(colnames(X), "gamma1", "gamma2", "sigma_v")
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
  names(out$par) <- c(colnames(X), "gamma", "sigmav")
  out$par[p + 2] <- exp(out$par[p + 2])
  hes <- numDeriv::hessian(likEXPall_orig, out$par,  X = X, list_y = list_y,  group = group,
                           alphainit = ols$alpha, niter = niter,
                           ws = ws, nodes = obj.gh$nodes)
  colnames(hes) <- rownames(hes) <- c(colnames(X), "gamma", "sigmav")
  out$invhes = solve(hes)
  out$se <- sqrt(diag(out$invhes))
#  out$se[p + 2] <- sqrt(exp(out$par[p + 2]*2) * out$se[p + 2]^2)
  names(out$se) <- c(colnames(X), "gamma", "sigma_v")
  out$alphai <- alphaiEHET(para = c(ogg$par[1:p],ogg$par[p+1], 0,ogg$par[p+2]), X = X, group = group, list_y = list_y,
                           list_z = list_z, niter=niter, uinit = rep(0, length(y)))
  out$uit <- extractUit_EHET(para = c(out$par[1:p],out$par[p+1],0,out$par[p+2]), X = X, y = y, z =z, group = group, alphai = out$alphai)

  return(out)
}

#############################################################
#    Gamma case
#############################################################
###Gong-Samaniego: para[1] is log(sigmaA), para[2] is log(muA), para[3] is log(lambda)
# psGong <- function(para, list_eta, list_y, alphastar){
#   out <- likGong(para[1], para[2], para[3], list_eta, list_y, alphastar)
#   alpha <- exp(para[2] - para[3])
#   sigma <- sqrt(exp(para[1] * 2) - alpha * exp(para[3] * 2))
#   cat(sigma, alpha, exp(para[3]), "-log=", out, "\n")
#  return(out)
# }

#minus log-likelihood for fixed alpha
# likGamma_alpha <- function(para, alpha, X, group, list_y, ymeans,
#                            ws, nodes,  Kinit = 5 , eps = 0.0001,
#                            niter = 10, umeth = "GS",
#                            trace = FALSE){
#   p <- ncol(X)
#   eta <- as.vector(X %*% para[1:p])
#   list_eta <- split(eta, group)
#   m0 <- if(umeth == "GS") 0 else 1
#   ll <- likG(para[p + 1], log(alpha), para[p + 2], list_eta, list_y,
#               ymeans - Kinit, ymeans + Kinit, eps, ws, nodes, niter, ymeans, m0)
#   if(trace) cat(c(para[1:p], exp(para[p+1]), alpha, exp(para)[p+2]), "-log=", ll, "\n")
#   return(-ll)
# }
#

## minus log-likelihood: para = beta, sigmaA, muA, lambda
likGamma <- function(para, X, group, list_y, alphacenter, ws, nodes,
                     Kinit = 5 , eps = 0.0001, trace = FALSE,
                     niter = 10, umeth = "GS"){
  p <- ncol(X)
  eta <- as.vector(X %*% para[1:p])
  list_eta <- split(eta, group)
  m0 <- if(umeth == "GS") 0 else 1
  alpha <- exp(para[2 + p] - para[3 + p])
  sigma2 <- exp(para[p + 1] * 2) - alpha * exp(para[p + 3] * 2)
  sigma <- if(sigma2 > 0) sqrt(sigma2) else eps * 100
  ll <- if(sigma == eps * 100) -10^6
  else likG(log(sigma), log(alpha), para[p + 3], list_eta, list_y, alphacenter - Kinit,
            alphacenter + Kinit, eps, ws, nodes, niter, alphacenter, m0)
  if(trace) cat(c(para[1:p], sigma, alpha, exp(para[p+3])), "-log=", ll, "\n")
  return(-ll)
}

## only muA and lambda
likGamma.cons <- function(para, b, sA, X, group, list_y, alphacenter, ws, nodes,
                          Kinit = 5 , eps = 0.0001, trace = FALSE,
                          niter = 10, umeth = "GS"){
  parcomp <- c(b, log(sA), para)
  out <- likGamma(parcomp, X, group, list_y, alphacenter, ws, nodes,  Kinit, eps,
                  trace, niter, umeth)
  return(out)
}


estim_G <- function(X, y, group, ols, nq = 25, eps = 10^-4, Kinit = 5,
                    niter = 10, umeth = "GS", initdelta = 0.5,
                    init = NULL,  trace = TRUE, int.trace = FALSE,
                    initNM = TRUE, useHess = TRUE){
  p <- ncol(X)
  list_y <- split(y, group)
  sA <- sqrt(sum(ols2$residuals^2) / (nrow(X) - max(group) - 1))
  ols$sA <- sA
  list_eta <- split(X * ols$beta, group)
  obj.gh <- statmod::gauss.quad(nq, "hermite")
  ws <- obj.gh$weights * exp(obj.gh$nodes^2)
  initA <- if(is.null(init)) c(ols$beta, log(ols$sA), log(ols$sA^2 - initdelta), 0)
  else init
  if(initNM){
    if(trace) cat("Searching for initial points...", "\n")
    mleinit <-  optim(c(log(sA^2 - initdelta), 0), b = ols$beta, sA = ols$sA,
                      likGamma.cons,  X = X, list_y = list_y,
                      alphacenter = ols$alpha + (sA^2 - initdelta),
                      group = mydat$g, ws = ws, nodes = obj.gh$nodes,
                      Kinit = Kinit, eps = eps, trace = int.trace,
                      control = list(trace = trace, reltol = eps * 10))
    initA <- c(ols$beta, log(ols$sA), mleinit$par)
    cat("done", "\n")
  }
  H <- diag(p + 3)
  if(useHess) {
    if(trace) cat("Computing Hessian at starting values...", "\t")
    H <- pracma::hessian(likGamma, initA, X = X, list_y = list_y,
                         alphacenter = ols$alpha + exp(initA[3]),
                         group = mydat$g, ws = ws, nodes = obj.gh$nodes,
                         Kinit = Kinit, eps = eps, trace = int.trace)
    cat("done", "\n")
  }
  E <- eigen(H)
  if(min(E$values)<=0) H <- diag(p+3)
  L1 <- solve(H)[lower.tri(H, diag = TRUE)]
  if(trace) cat("Starting final optimization...", "\n")
  mle <- ucminf::ucminf(initA,
                        likGamma,  X = X, list_y = list_y,
                        alphacenter = ols$alpha + exp(mleinit$par[1]),
                        group = mydat$g, ws = ws, nodes = obj.gh$nodes,
                        Kinit = Kinit, eps = eps, trace = int.trace,
                        control = list(trace = trace, grtol = eps * 10^3, invhessian.lt = L1))
  cat("done", "\n")
  out <- mle
  alpha <- exp(mle$par[2 + p] - mle$par[3 + p])
  sigma <- sqrt(exp(mle$par[p + 1] * 2) - alpha * exp(mle$par[p + 3] * 2))
  out$paraorig <- c(mle$par[1], sigma, alpha, exp(mle$par[p + 3]))
  attr(out$paraorig, "names") <- c("beta", "sigma", "gamma", "lambda")
  attr(out$par, "names") <- c("beta", "lnsigmaA", "lnmuA", "lnlambda")
  return(out)
}



