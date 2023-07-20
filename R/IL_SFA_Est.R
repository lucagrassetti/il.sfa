#' Integrated Likelihood estimation function
#'
#' Adopts the IL approach to the estimation of Stochastic Frontier Model via marginal likelihood.
#'
#' @param X the matrix of independent variables
#' @param y the response variable vector
#' @param distr the distributional form of the efficiency related error terms (default is Exp)
#' @param het a logical indicating the
#' @param z the matrix of variables influencing the efficiency
#' @param group the vector of grouping variable
#' @param nq the number of quadrature points used in the numerical integration routine
#' @param niter the number of iterations used in the numerical integration routine
#' @param init the vector of initial values for likelihood parameters
#' @param useHess a logical: if TRUE the likelihood maximization is initialized with the numerical Hessian otherwise the Hessian is initialized as a unit diagonal
#' @param trace trace parameter used in ucminf maximization algorithm
#' @param grtol grtol parameter used in ucminf maximization algorithm (gradient tollerance)
#' @param grad grad parameter used in ucminf maximization algorithm
#'
#' @return a list including the point estimates, the numerical hessian, the standard errors of the estimates and the firm specific effects (alpha_i) and the estimated efficiency related terms (u_it).
#'
#' @examples
#' ### Reproduce one single dataset with the setting used in Table 1, (a) lambda = 1
#' ### Belotti TFE - N=100, T=5
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
