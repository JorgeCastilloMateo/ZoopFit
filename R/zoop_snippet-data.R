#' @title Snippet of zooplankton data
#' 
#' @description Snippet of 200 zooplankton data.
#' 
#' @note Note that the data were collected by and are property of the Center 
#'   for Coastal Studies. We share a small 200 row subset of these data here to 
#'   show how the model fitting works. To fully reproduce this analysis, you 
#'   will need to contact CCS to get access via a data sharing agreement.
#' 
#' @docType data
#' @keywords datasets
#' @name zoop_snippet
#' @usage data("zoop_snippet")
#' @format Two matrices with coordinates, \code{coords}, and whale and regime 
#'   covariates, \code{X}. Nine vectors with the measurement \code{date}, 
#'   indicators of the associated measurements being 1 (available) or 0 (NA),
#'   \code{oneT1}, \code{oneT2}, \code{oneY1}, and \code{oneY2}; and 
#'   temperature and zooplankton measurements \code{T1} (CTD), \code{T2} 
#'   (Therm), \code{Y1} (log surface), and \code{Y2} (log oblique).
"coords"
"X"
"date"
"oneT1"
"oneT2"
"oneY1"
"oneY2"
"T1"
"T2"
"Y1"
"Y2"
