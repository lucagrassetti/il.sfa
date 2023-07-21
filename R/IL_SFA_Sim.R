#' The Data Generating Process used in the development of simulation studies.
#'
#' The DGP is defined as
#' \deqn{y_{it} = \alpha_i + x_{it}^T \beta - u_{it} + v_{it}}
#' where \eqn{u_{it}} can be Half-Normal or Exponential and in this last
#' case can also consider the mean of the Exponential distirbution as
#' \deqn{\lambda_{it} = \gamma_0 + \gamma_1 z_{it}.}
#' This last definition also brings to time-constant cases when
#' \eqn{z_{it}} is considered constant over time.
#'
#' @param n the cross-sectional dimension.
#' @param m the time dimension.
#' @param param the vector of model parameters (the frontier regression parameters - \eqn{\beta},
#' the idiosyncratic error variance - \eqn{\sigma_v},
#' the vector of the parameters describing the efficiency behaviour - \eqn{\gamma}).
#' @param model the kind of distribution considered to generate the efficiency related terms.
#' @param seed.val the seed value used to make the simulations reproducible.
#' @param nsim an integer indicating the number of simulated dataset.
#' @param zvar a logical idenfying the case in which the efficiency is time variant.
#'
#' @return a list object including:
#' @return -- the simulated y vectors (a matrix daty of dimension \eqn{nsim \times n\cdot m}),
#' @return -- the simulated of efficiency related terms (a matrix datu of dimension \eqn{nsim \times n\cdot m}),
#' @return -- the simulated of independent variables (a matrix \eqn{x} of dimension \eqn{n\cdot m \times p}),
#' @return -- the simulated of grouping identifiers (a vector \eqn{g} of length \eqn{n\cdot m}), and
#' @return -- the simulated values of the variable influencing the efficiency level (a vector \eqn{datz} of length \eqn{n\cdot m}).
#'
#' @references Belotti, F., & Ilardi, G. (2018). Consistent inference in
#' fixed-effects stochastic frontier models. Journal of Econometrics, 202(2), 161--177.
#' @references Bellio, R., & Grassetti, L. (2023). Efficient estimation of true
#' fixed-effects stochastic frontier models. WP
#'
#' @examples
#' # This code reproduce the results in Table 1 panel (a)
#' # Belotti and Ilardi (2018) "Consistent inference in fixed-effects
#' # stochastic frontier models" - N=100, T=5.
#' # The example is limited to 100 replications to make the code
#' # fairly executable. Changing the nsim parameter, one can consider the
#' # original simulation with 1000 replications.
#'
#' N <- 100
#' T <- 5
#' para <- c(1,-1.5,1,0.25)
#' nsimul <- 100
#' dat <- DGP(N, T, list(beta=para[1], sigmav=para[4], gamma=para[2:3]), model = "Exp",
#'            seed.val = 22, nsim = nsimul)
#' res100 <- matrix(0, nrow = nsimul, ncol=4)
#' for(i in 1:nrow(res100)){
#'   print(i)
#'   mydat <- list(y = dat$daty[i,], x = data.frame(dat$x), g = dat$g, z = dat$datz)
#'   mle <- il_sfa(X = mydat$x, y = mydat$y, distr = "Exp", het = TRUE, z = mydat$z,
#'                 group = mydat$g, useHess = TRUE)
#'   res100[i,] <- mle$par
#'   }
#'   # Compute the MSE of the estimates
#'   for(i in 1:4) print(MSE(res100[,i], para[i]))
#'   # Compare them with the estimated values included in the
#'   # package data (res_N100_Stata.rda)
#'   data(res_N100_Stata)
#'   for(i in 1:4) print(MSE(res_N100_Stata[1:nsimul,i], para[i]))
#'
#' @export
DGP <- function(n, m, param, model = "HN", seed.val = 0, nsim = 1000, zvar = FALSE)
{
  g <- sort(rep(1:n, m))
  daty <- datu <- matrix(NA, nrow = nsim, ncol = n * m)
  set.seed(seed.val)
  alpha <- rnorm(n)
  w <- rnorm(n * m)
  beta <- param$beta
  s2v <- param$sigmav^2
  x <- 0.5 * alpha[g] + sqrt(0.5^2) * w
  z <- if (zvar) rnorm(n * m, sqrt(0.0625)) else rnorm(n, sqrt(0.0625))[g]
  for(i in 1:nsim) {
    if(model == "HN")  {
      s2u <- param$sigmau^2
      sigma <- sqrt(s2u + s2v)
      lambda <- param$sigmau / param$sigmav
      u <- abs(rnorm(n * m, sd = sqrt(s2u)))
    }
    if(model == "Exp")  {
      su <- exp(param$gamma[1] +  param$gamma[2] * z)
      u <- rexp(n *m, rate = 1 / su)
    }
    daty[i,] <-  alpha[g] + x * beta + rnorm(n * m, sd = sqrt(s2v)) - u
    datu[i,] <- u
  }
  return(list(daty = daty, datu = datu, x = as.matrix(x, ncol = 1), g = g, datz = z))
}

