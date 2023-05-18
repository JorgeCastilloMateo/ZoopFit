#' Recover ZOOP
#' 
#' @description This function recovers ZOOP from a fitted model
#' 
#' @param model model fitted with the function \code{\link{GibbsZoop}}
#' @param X covariates (whale, r2, r3?): \eqn{N \times p}
#' @param date \eqn{N \times 1}
#' @param coords \eqn{N \times 2}
#' @param do.mean return the mean? TRUE Or samples? FALSE
#' @param log.scale return log-scale? TRUE Or original scale? FALSE
#' @return matrix: rows are simuations columns are (log) ZOOP
#' @examples
#' #ZOOP <- recoverZoop(model, X, date, coords)
#' #plot(x=ZOOP, y=auxY1)
#' #plot(x=ZOOP, y=auxY2)
#' # 
#' #beta <- matrix(nrow = 17, ncol = 4)
#' #for (i in 1:17) {
#' #   beta[i,] <- coefficients(lm(auxY2 ~ ZOOP, subset = year(date) == i + 2002))
#' #}
#' #plot(beta[,1])
#' #plot(beta[,2])
#' #plot(beta[,3])
#' #plot(beta[,4])
#' 
#' @export 
recoverZoop <- function(model, X, date, coords, do.mean = FALSE, log.scale = FALSE) {
  
  N <- nrow(X)
  p <- ncol(X)
  B <- nrow(model)
  
  # recover Xeta
  year <- lubridate::year(date)
  T <- length(unique(year))
  firstYear <- min(year) - 1
  nt <- c() 
  for (tInd in 1:T) 
    nt[tInd] <- nrow(unique(coords[year == firstYear + tInd,]))
  Ent <- sum(nt)
  my_fun <- function(x) match(data.frame(t(x)), data.frame(t(unique(x))))
  ti <- c()
  for (tInd in 1:T) 
    ti <- c(ti, my_fun(coords[year == firstYear + tInd,]) + sum(nt[0:(tInd-1)]))
  Xeta <- matrix(nrow = N, ncol = Ent)
  for (nInd in 1:Ent)
    Xeta[,nInd] <- ti == nInd
  Xeta <- 1 * Xeta

  # ZOOP
  ZOOP <- matrix(nrow = B, ncol = N)
  
  for (i in 1:N) {
    ZOOP[,i] <- model[,paste0("beta", 1:p)] %*% X[i,] + 
      model[, "betaTemp"] * model[, paste0("temp", i)] + 
      model[, which(Xeta[i,] == 1)]
  }
  
  if (!log.scale) {
    ZOOP <- exp(ZOOP)
  }
  
  if (do.mean) {
    ZOOP <- colMeans(ZOOP)
  }
  
  return(ZOOP)
}
