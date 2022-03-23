#' @title \code{scMET}: Bayesian modelling of DNA methylation at single-cell
#'   resolution.
#'
#' @description Package for analysing single-cell DNA methylation datasets.
#'   scMET performs feature selection, by identifying highly variable features,
#'   and also differential testing, based on mean but also more importantly on
#'   variability between two groups of cells.
#'
#' @return scMET main package documentation.
#'
#' @seealso \code{\link{scmet}}, \code{\link{scmet_differential}},
#'   \code{\link{scmet_hvf}}
#'
#' @author C.A.Kapourani \email{kapouranis.andreas@@gmail.com}
#'
#' @docType package
#' @name scMET-package
#' @useDynLib scMET, .registration = TRUE
#' @import methods
#' @import Rcpp BiocStyle rstantools
#' @importFrom rstan sampling vb
#' @importFrom cowplot plot_grid
#' @rawNamespace importFrom(magrittr,"%>%")
#' @rawNamespace importFrom(data.table,":=")
#'
#' @references Stan Development Team (2020). RStan: the R interface to Stan. R
#'   package version 2.19.3. https://mc-stan.org
#'
.datatable.aware <- TRUE
NULL

.onLoad <- function(libname = find.package("scMET"), pkgname = "scMET"){
  # CRAN Note avoidance
  if (getRversion() >= "2.15.1")
    utils::globalVariables(
      # sample file names from taxstats
      c(# we use the magrittr pipe
        "."
      )
    )
  invisible()
}


#' @title Synthetic methylation data from a single population
#'
#' @description Small synthetic data for quick analysis, mostly useful for
#'   performing feature selection and capturing mean-variance relationship with
#'   scMET.
#'
#' @return A list object with simulated data.
#'
"scmet_dt"



#' @title Synthetic methylation data from two groups of cells
#'
#' @description Small synthetic data for quick analysis, mostly useful for
#'   showing the differential analysis one can perform using scMET.
#'
#' @return A list object with simulated data.
#'
"scmet_diff_dt"
