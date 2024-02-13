#' True Fixed Effects Stochastic Frontier Model estimation via Integrated Likelihood.
#'
#' This function adopts the IL approach to the estimation of some kind of SF models
#' (see Bellio and Grassetti, 2023) such
#' as those based on Normal-Half Normal, N-Exponential, and N-Gamma error terms convolutions.
#' The function adopts a consistent and computationally efficient approach to the
#' model estimation.
#' The model specification is as
#' \deqn{y_{ij} = \alpha_i + x_{it}^T\beta - u_{it} + v{it},}
#' where the efficiency-related error term can be specified as
#' \deqn{u_{it} \sim HN(0, \sigma_u^2),}
#' \deqn{u_{it} \sim Exp(\sigma_u^2),}
#' \deqn{u_{it} \sim Gamma(\gamma_1, \gamma_2)\mbox{, or}}
#' \deqn{u_{it} \sim Exp(\sigma_{uit}^2) \mbox{ with } \sigma_{uit} = e^{\gamma_0+\gamma_1z_it}.}
#' The idiosyncratic error term is considered Gaussian and homoscedastic.
#'
#' @param X the matrix of independent variables.
#' @param y the vector of the response variable.
#' @param distr the distributional form of the efficiency-related error terms (default
#' is \code{Exp}).
#' @param het a logical indicating the case of heteroscedastic efficiency related error
#' terms (default is \code{TRUE}, but it works for the Exponential case only).
#' @param z the vector of values of the variable influencing the efficiency (for the Exponential case only).
#' @param group the vector of the grouping variable defining the panel.
#' @param nq the number of quadrature points used in the numerical integration routine
#' (default is 25).
#' @param niter the number of iterations used in the numerical integration routine for
#' locating the mode (default is 10).
#' @param init the optional vector of initial values for the maximization of
#' integrated likelihood (default is \code{NULL}).
#' @param useHess a logical: if \code{TRUE}, the likelihood maximization is initialized
#' with the numerical Hessian; otherwise, the Hessian is initialized as a unit diagonal.
#' @param trace trace parameter used in \code{ucminf} function (if \code{trace} is positive, then
#' detailed tracing information is printed for each iteration).
#' @param grtol grtol parameter used in ucminf function (the algorithm stops when
#' the largest absolute value of the gradient is less than \code{grtol}).
#' @param grad grad parameter used in \code{ucminf} function (controls the type
#' of finite difference approximation to be used for the gradient).
#' @param eps tolerance for the gamma estimation process, the default is 10^-4.
#' @param Kinit radius of the interval for the golden-section search for locating
#' the mode with respect to \eqn{\alpha_i}. Default is 5.
#' @param umeth method for locating the mode with respect to \eqn{\alpha_i}.
#' Default is \code{GS} for golden-section search, any other value for Newton-Raphson.
#' @param initdelta initial guess for \eqn{sigma^2}, only used to determine the
#' initial parameter \eqn{\mu_G} when \code{initNM = FALSE}. Default is 0.5.
#' @param int.trace a logical for internal tracing. Default is \code{FALSE}.
#' @param initNM logical for activating the initial parameter search using first
#' the exponential model fit result and then the Nelder-Mead method for the
#' gamma model (strongly suggested). Detault is \code{TRUE}.
#'
#' @return A list including
#' @return -- the parameters point estimates (\code{par}),
#' @return -- the inverse of the numerical hessian (\code{invhes}),
#' @return -- the standard errors of the estimated parameters (\code{se})
#' @return -- the firm specific effects (\code{alphai}, still not available for the gamma model), and
#' @return -- the estimated efficiency related terms (\code{uit}, still not available for the gamma model).
#'
#' @references Belotti, F., & Ilardi, G. (2018). Consistent inference in
#' fixed-effects stochastic frontier models. Journal of Econometrics, 202(2), 161--177.
#' @references Bellio, R., & Grassetti, L. (2023). Efficient estimation of true
#' fixed-effects stochastic frontier models. Submitted
#'
#' @examples
#' # Generate one single dataset with the setting used in Table 1, (a)
#' # Belotti and Ilardi (2018) "Consistent inference in fixed-effects
#' # stochastic frontier models" - N=100, T=5.
#'
#' N <- 100
#' T <- 5
#' para <- c(1,-1.5,1,0.25)
#' dat <- DGP(N, T, list(beta=para[1], sigmav=para[4], gamma=para[2:3]), model = "Exp",
#'            seed.val = 22, nsim = 1)
#' mydat <- list(y = dat$daty[1,], x = data.frame(dat$x), g = dat$g, z = dat$datz)
#' mle <- il_sfa(X = mydat$x, y = mydat$y, distr = "Exp", het = TRUE, z = mydat$z,
#'               group = mydat$g, useHess = TRUE)
#' mle
#'
#' @export
il_sfa <- function(X, y, distr = "Exp", het = TRUE, z = NULL, group,nq = 25, niter = 10, init = NULL, Kinit = 5, initdelta = 0.5, useHess = TRUE, trace = 0, int.trace = FALSE, initNM = TRUE, grtol = 10^-6, eps = 10^-4, grad = "central", umeth = "GS")
{
  X_plus <- model.matrix(~factor(group)-1+X)
  ols2 <- RcppEigen::fastLmPure(X_plus, y)
  ols <- list(alpha = ols2$coef[1:max(group)], beta = ols2$coef[max(group)+1], residuals = ols2$residuals)
  ols$sA <- sqrt(sum(ols$residuals^2) / (nrow(X) - max(group) - 1))
  if(distr == "Exp")
  {
    if(het == TRUE)
    {
      mle <- estim_EHET(X = X, y = y, z = z, group = group,
                        ols = ols, nq = nq, niter = niter,
                        init = init, useHess = useHess,
                        trace = trace, grtol = grtol,
                        grad = grad)
      p <- ncol(X)
      if(mle$par[1:p] == ols$beta)
      {
        print("Hessian computation failed: using useHess = FALSE")
        mle <- estim_EHET(X = X, y = y, z = z, group = group,
                          ols = ols, nq = nq, niter = niter,
                          init = init, useHess = FALSE,
                          trace = trace, grtol = grtol,
                          grad = grad)
      }
    } else
    {
      mle <- estim_EXP(X = X, y = y, group = group,
                       ols = ols, nq = nq, niter = niter,
                       init = init, useHess = useHess,
                       trace = trace, grtol = grtol,
                       grad = grad)
      p <- ncol(X)
      if(mle$par[1:p] == ols$beta)
      {
        print("Hessian computation failed: using useHess = FALSE")
        mle <- estim_EXP(X = X, y = y, group = group,
                         ols = ols, nq = nq, niter = niter,
                         init = init, useHess = FALSE,
                         trace = trace, grtol = grtol,
                         grad = grad)
      }
    }
  }
  if(distr == "HN")
  {
    if(het == TRUE)
    {
      print("Heteroschedastic errors are still not supported for Half-Normal distribution.")
      print("The estimation results are obtained under homoscedasticity assumption.")
      het <- FALSE
      mle <- estim_HN(X = X, y = y, group = group,
                      ols = ols, nq = nq, niter = niter,
                      init = init, useHess = useHess,
                      trace = trace, grtol = grtol,
                      grad = grad)
      p <- ncol(X)
      if(mle$par[1:p] == ols$beta)
      {
        print("Hessian computation failed: using useHess = FALSE")
        mle <- estim_HN(X = X, y = y, group = group,
                        ols = ols, nq = nq, niter = niter,
                        init = init, useHess = FALSE,
                        trace = trace, grtol = grtol,
                        grad = grad)
      }
    } else
    {
      mle <- estim_HN(X = X, y = y, group = group,
                      ols = ols, nq = nq, niter = niter,
                      init = init, useHess = useHess,
                      trace = trace, grtol = grtol,
                      grad = grad)
      p <- ncol(X)
      if(mle$par[1:p] == ols$beta)
      {
        print("Hessian computation failed: using useHess = FALSE")
        mle <- estim_HN(X = X, y = y, group = group,
                        ols = ols, nq = nq, niter = niter,
                        init = init, useHess = FALSE,
                        trace = trace, grtol = grtol,
                        grad = grad)
      }
    }
  }
  if(distr == "Gamma")
  {
    if(het == TRUE)
    {
      print("Heteroschedastic errors are still not supported for Gamma distribution.")
      print("The estimation results are obtained under homoscedasticity assumption.")
      het <- FALSE
      mle <- estim_G(X = X, y = y, group = group,
                     ols = ols, nq = nq, eps = eps, Kinit = Kinit,
                     niter = niter, umeth = umeth, initdelta = initdelta,
                     int.trace = int.trace, init = init,  trace = trace,
                     initNM = initNM,
                     useHess = useHess)
      p <- ncol(X)
      if(mle$par[1:p] == ols$beta)
      {
        print("Hessian computation failed: using useHess = FALSE")
        mle <- estim_G(X = X, y = y, group = group,
                       ols = ols, nq = nq, eps = eps, Kinit = Kinit,
                       niter = niter, umeth = umeth, initdelta = initdelta,
                       int.trace = int.trace, init = init,  trace = trace,
                       initNM = initNM,
                       useHess = FALSE)
      }
    } else
    {
      mle <- estim_G(X = X, y = y, group = group,
                     ols = ols, nq = nq, eps = eps, Kinit = Kinit,
                     niter = niter, umeth = umeth, initdelta = initdelta,
                     int.trace = int.trace, init = init,  trace = trace,
                     initNM = initNM,
                     useHess = useHess)
      p <- ncol(X)
      if(mle$par[1:p] == ols$beta)
      {
        print("Hessian computation failed: using useHess = FALSE")
        mle <- estim_G(X = X, y = y, group = group,
                       ols = ols, nq = nq, eps = eps, Kinit = Kinit,
                       niter = niter, umeth = umeth, initdelta = initdelta,
                       int.trace = int.trace, init = init,  trace = trace,
                       initNM = initNM,
                       useHess = FALSE)
      }
    }
  }
  if(sum(distr %in% c("HN", "Exp", "Gamma"))==0)
  {
    print("Error: distr must be HN (Half-Normal), Gamma (Gamma), or Exponential (Exp).")
    mle <- NULL
  }
  mle
}
