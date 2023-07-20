#' Results from the Stata estimation function using PDE method
#'
#' Contains the estimation results of the 1000 replications of the simulation experiment also developed in this package
#'
#' @format A data frame with 1000 rows and 4 variables (the model parameters)
#'  \describe{
#'      \item{b1}{The slope parameter in the frontier function}
#'      \item{b2}{The heteroschedasticity parameter gamma_0}
#'      \item{b3}{The heteroschedasticity parameter gamma_1}
#'      \item{b4}{The idiosyncratic error variance}
#'      }
#'
#' @source {Stata estimation results obtained considering the DGP described in the paper}
#'
#' @examples
#' data(res_N100_Stata)
"res_N100_Stata"
