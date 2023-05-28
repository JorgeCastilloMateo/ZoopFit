#' Validation RMSE MAE CRPS CVG 
#' 
#' @importFrom stats quantile
#' 
#' @description This function validates the model
#' 
#' @param Y actual data
#' @param Yhat replicated or simulated data, rows are replicates 
#'   columns are data
#' @param alpha CVG level
#' @return matrix with RMSE, MAE, CRPS and CVG
#' 
#' @author Jorge Castillo-Mateo
#' @export 
validation <- function(Y, Yhat, alpha = 0.1) {
  
  x <- data.frame("RMSE" = .RMSE(Y, Yhat),
                  "MAE"  = .MAE(Y, Yhat),
                  "CRPS" = .CRPS(Y, Yhat),
                  "CVG"  = .CVG(Y, Yhat, alpha = alpha))
  
  return(x)
}

.RMSE <- function(Y, Yhat) sqrt(mean((Y - colMeans(Yhat))^2, na.rm = TRUE))
.MAE  <- function(Y, Yhat) mean(abs(Y - colMeans(Yhat)), na.rm = TRUE)
.CRPS <- function(Y, Yhat) {
  n <- length(Y)
  crps <- 0
  for (i in 1:n) {
    if (!is.na(Y[i])) {
      crps <- crps + mean(abs(Yhat[,i] - Y[i])) - mean(abs(outer(Yhat[,i], Yhat[,i], FUN = `-`))) / 2
    }
  }
  crps <- crps / sum(!is.na(Y))
  return(crps)
}
.CVG <- function(Y, Yhat, alpha = 0.1) {
  n <- length(Y)
  q <- c(alpha / 2, 1 - alpha / 2)
  cvg <- 0
  for (i in 1:n) {
    if (!is.na(Y[i])) {
      LU  <- stats::quantile(Yhat[,i], probs = q)
      cvg <- cvg + as.numeric(LU[1] <= Y[i] && Y[i] <= LU[2])
    }
  }
  cvg <- cvg / sum(!is.na(Y))
  return(cvg)
}
