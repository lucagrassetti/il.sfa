#' Estimation results from Stata (a comparison set of estimates)
#'
#' Collects the estimation results obtained in Stata using sftfe function
#' (Belotti et al., 2013) and the PDE method. These results are used in the
#' reference paper Bellio and Grassetti (2023) to compare ILE and PDE results.
#' The reported results refer to the heteroscedastic exponential SFM case. Data
#' are obtained considering the default arguments in the DGP() function included
#' in this library.
#'
#' @format A data frame with 1000 rows and 4 variables (the model parameters)
#'  \describe{
#'      \item{b1}{The slope parameter in the frontier function}
#'      \item{b2}{The heteroschedasticity parameter gamma_0}
#'      \item{b3}{The heteroschedasticity parameter gamma_1}
#'      \item{b4}{The idiosyncratic error variance}
#'      }
#'
#' @source {Stata estimation results obtained considering the DGP described in the
#' reference papers.
#' N = 100, T = 5, beta = 1, gamma_0 = 1.5, gamma_1 = 1 and sigma_v = 0.25}
#'
#' @references Bellio, R. & Grassetti, L. (2023). Efficient estimation of true
#' fixed-effects stochastic frontier models. Submitted.
#' @references Belotti, F., Daidone, S., Ilardi, G. & Atella, V. (2013). Stochastic frontier
#' analysis using Stata. The Stata Journal 13, 719-758
#'
#' @examples
#' data(res_N100_Stata)
"res_N100_Stata"
