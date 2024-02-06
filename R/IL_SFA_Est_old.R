#' True Fixed Effects Stochastic Frontier Model estimation via Integrated Likelihood.
#'
#' This function adopts the IL approach to the estimation of some kind of SF models
#' (see Bellio and Grassetti, 2023) such
#' as those based on N-HN and N-Exponential error terms convolutions.
#' The function adopts a consistent and computationally efficient approach to the
#' model estimation.
#' The model specification is as
#' \deqn{y_{ij} = \alpha_i + x_{it}^T\beta - u_{it} + v{it},}
#' where the efficiency-related error term can be estimated as \eqn{u_{it} \sim HN(0, \sigma_u^2)},
#' \eqn{u_{it} \sim Exp(\sigma_u^2)}, or \eqn{u_{it} \sim Exp(\sigma_{uit}^2)} and
#' \eqn{\sigma_{uit} = e^{\gamma_0+\gamma_1z_it}}.
#' The idiosyncratic error term is considered Gaussian and homoscedastic.
#'
#' @param X the matrix of independent variables.
#' @param y the vector of the response variable.
#' @param distr the distributional form of the efficiency-related error terms (default is Exp).
#' @param het a logical indicating the case of heteroscedastic efficiency related error
#' terms (default is TRUE, but it works for the Exponential case only).
#' @param z the vector of values of the variable influencing the efficiency (for the Exponential case only).
#' @param group the vector of the grouping variable defining the panel.
#' @param nq the number of quadrature points used in the numerical integration routine (default is 25).
#' @param niter the number of iterations used in the numerical integration routine for locating the mode (default is 10).
#' @param init the optional vector of initial values for the maximization of integrated likelihood (default is NULL).
#' @param useHess a logical: if TRUE, the likelihood maximization is initialized with the numerical
#' Hessian; otherwise, the Hessian is initialized as a unit diagonal.
#' @param trace trace parameter used in ucminf function (if trace is positive, then
#' detailed tracing information is printed for each iteration).
#' @param grtol grtol parameter used in ucminf function (the algorithm stops when
#' the largest absolute value of the gradient is less than grtol).
#' @param grad grad parameter used in ucminf function (controls the type
#' of finite difference approximation to be used for the gradient).
#'
#' @return A list including
#' @return -- the parameters point estimates - \eqn{par},
#' @return -- the inverse of the numerical hessian - \eqn{invhes},
#' @return -- the standard errors of the estimated parameters - \eqn{se}
#' @return -- the firm specific effects - \eqn{alphai}, and
#' @return -- the estimated efficiency related terms -  - \eqn{uit}.
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
il_sfa <- function(X, y, distr = "Exp", het = TRUE, z, group, nq = 25, niter = 10,
                   init = NULL, useHess = TRUE,
                   trace = 0, grtol = 10^-6,
                   grad = "central")
{
  formula <- paste(c("y", paste(c(names(X), "(1|group)"), collapse = " + ")), collapse = " ~ ")
  ols1 <- lme4::lmer(formula, data.frame(y = y, X, group = group, z = z))

  ols <- list(alpha = tapply(y, group, mean), beta = lme4::fixef(ols1)[1+(1:ncol(X))])
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
  if(sum(distr %in% c("HN", "Exp"))==0)
  {
    print("Error: distr must be HN - Half-Normal or Exponential - Exp.")
    mle <- NULL
  }
  mle
}
