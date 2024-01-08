#' @name ZoopFit-package
#' @aliases ZoopFit
#' @docType package
#' 
#' @title \strong{ZoopFit}: Statistical Modeling for Zooplankton Availability
#' 
#' @description \strong{ZoopFit} is the companion package for the paper 
#'   Castillo-Mateo et al. (2023). The package includes functions to fit and
#'   predict with the double data fusion and calibration models from the paper.
#'   
#' @details The back-end code of this package is built under \code{C++} 
#'   language. 
#'   
#'   Main functions used:
#' 
#'   > \code{\link{GibbsZoop}}
#'   
#'   > \code{\link{predictZoop}}
#'   
#'   Some other important functions:
#'   
#'   > \code{\link{meanGP.GibbsZoop}}
#'   
#'   > \code{\link{predictZoopMap}}
#'   
#'   > \code{\link{recoverZoop}}
#'   
#'   > \code{\link{validation}}
#'   
#' @note \strong{Warning:} Please note that the code provided is for the paper only 
#'   and a generalization is under development. I do not guarantee that it 
#'   will work outside of the data from the paper and it will likely lead to 
#'   errors if used with other datasets.
#' 
#' @references 
#'   Castillo-Mateo J, Gelfand AE, Hudak CA, Mayo CA, Schick RS (2023).
#'   “Space-time multi-level modeling for zooplankton abundance employing double data fusion and calibration.” 
#'   \emph{Environmental and Ecological Statistics}, \strong{30}(4), 769-795.
#'   \doi{10.1007/s10651-023-00583-6}.
#'
#' @author Jorge Castillo-Mateo <jorgecastillomateo@gmail.com>
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib ZoopFit
NULL  