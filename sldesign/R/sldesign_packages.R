
#' The 'sldesign' package for skip-Levels design for phase I oncology study
#'
#' @description 
#' A tool to calculate the sufficient range of the response rates 
#' that the skip-levels  (SL) design performs better than the uniformly-assign 
#' (UA) design, as well as conduct simulation studies for the phase 1 
#' dose-escalation and dose-expansion using the SL design or the UA design 
#' and find the recommended phase 2 dose admissible region using Bayesian 
#' 2-parameter logistic regression.
#'
#' @docType package
#' @name sldesign_packages
#' @aliases sldesign
#' @useDynLib sldesign, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom BOIN get.boundary select.mtd
#' @import dplyr
#' @import tidyr
#' @import Iso
#' @import coda
#' @import DoseFinding
#' 
#' @references
#' Stan Development Team (NA). RStan: the R interface to Stan. 
#' R package version 2.32.6. https://mc-stan.org
#' 
NULL
