#' The Data Generating Process used in the simulation studies.
#'
#' The DGP is defined as
#' \deqn{y_{it} = \alpha_i + x_{it}^T \beta - u_{it} + v_{it}}
#' where \eqn{u_{it}} can be Half-Normal, Gamma, or Exponential. In this last
#' case, it can also consider the mean of the Exponential distribution as
#' \deqn{\sigma_{uit} = \exp\left(\gamma_0 + \gamma_1 z_{it}\right).}
#' This last definition also brings to time-constant cases when
#' \eqn{z_{it}} is considered constant over time.
#'
#' @param n the cross-sectional dimension.
#' @param m the time dimension.
#' @param param the vector of model parameters (the frontier regression parameters - \eqn{\beta},
#' the idiosyncratic error variance - \eqn{\sigma_v},
#' the vector of the parameters describing the efficiency behavior - \eqn{\gamma}).
#' \eqn{\gamma} depends on \eqn{u_{it}} distribution choice.
#' In the half-normal case, it corresponds to \eqn{\sigma_u}. In the exponential case, it collects
#' \eqn{\gamma_0} and \eqn{\gamma_1} defined in the formula above.
#' In the gamma case, it corresponds to the shape and scale parameters of the gamma distribution.
#' @param model the kind of distribution considered to generate the efficiency-related terms.
#' \code{HN}, \code{Exp}, or \code{Gamma} options are available. Default is \code{HN}.
#' @param seed.val the seed value used to make the simulations reproducible. Default value is 0.
#' @param nsim an integer indicating the number of simulated datasets. Default value is 1000.
#' @param zvar a logical identifying the case in which the efficiency is time variant. Default is \code{FALSE}.
#'
#' @return a list object including:
#' @return -- the simulated y vectors (a matrix daty of dimension \eqn{nsim \times n\cdot m}),
#' @return -- the simulated of efficiency-related terms (a matrix datu of dimension \eqn{nsim \times n\cdot m}),
#' @return -- the simulated of independent variables (a matrix \eqn{x} of dimension \eqn{n\cdot m \times p}),
#' @return -- the simulated of grouping identifiers (a vector \eqn{g} of length \eqn{n\cdot m}), and
#' @return -- the simulated values of the variable influencing the efficiency level (a vector \eqn{datz} of length \eqn{n\cdot m}).
#'
#' @description This code allows to reproduce the results in Table 1 panel (a)
#' Belotti and Ilardi (2018) "Consistent inference in fixed-effects
#' stochastic frontier models". The sample size is fixed to \code{n=100} (individuals)
#' and \code{m=5} (time points).
#' The example is limited to 10 replications to make the code
#' fairly executable. Changing the nsim parameter, one can consider the
#' original simulation with 1000 replications. \code{set.seed = 22} is used to generate the
#' data analysed with Stata and R and compared in Bellio and Grassetti (2023).
#' In the cited papers, the simulations use a number of replications equal to 1000.
#' The available distributions for the efficiency term are homoscedastic half normal (\code{HN}),
#' homoscedastic and heteroscedastic Exponential (\code{Exp}),
#' and homoscedastic Gamma (\code{Gamma}).
#'
#' @references Bellio, R. & Grassetti, L. (2023). Efficient estimation of true
#' fixed-effects stochastic frontier models. Submitted
#' @references Belotti, F. & Ilardi, G. (2018). Consistent inference in
#' fixed-effects stochastic frontier models. Journal of Econometrics, 202(2), 161--177.
#'
#' @examples
#'
#' N <- 100
#' T <- 5
#' para <- c(1,-1.5,1,0.25)
#' nsimul <- 10
#' dat <- DGP(N, T, list(beta = para[1], sigmav = para[4], gamma = para[2:3]), model = "Exp",
#'            seed.val = 22, nsim = nsimul)
#' res <- matrix(0, nrow = nsimul, ncol=4)
#' for(i in 1:nrow(res)){
#'   print(i)
#'   mle <- il_sfa(X = dat$x, y = dat$daty[i,], distr = "Exp", het = TRUE, z = dat$datz,
#'                 group = dat$g, useHess = TRUE)
#'   res[i,] <- mle$par
#'   }
#' # Compute the MSE of the estimates
#' for(i in 1:4) print(MSE(res[,i], para[i]))
#' # Compare them with the estimated values included in the
#' # package data (res_N100_Stata.rda)
#' data(res_N100_Stata)
#' for(i in 1:4) print(MSE(res_N100_Stata[1:nsimul,i], para[i]))
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
    if(model == "Gamma")  {
      u <- rgamma(n * m, shape = param$gamma[1], scale = param$gamma[2])
    }
    daty[i,] <-  alpha[g] + x * beta + rnorm(n * m, sd = sqrt(s2v)) - u
    datu[i,] <- u
  }
  x <- as.matrix(x, ncol = 1)
  colnames(x) <- paste0("x", 1:ncol(x))
  return(list(daty = daty, datu = datu, x = x, g = g, datz = z))
}

