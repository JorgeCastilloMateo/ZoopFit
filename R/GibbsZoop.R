#' Double Fusion-Calibration Zoop Models
#' 
#' @importFrom lubridate year
#' @importFrom stats dist
#' 
#' @description This function fits the double fusion and calibration zoop model
#' 
#' @param Y1 log Surface zoop: \eqn{N_{Y1} \times 1}
#' @param Y2 log Oblique zoop: \eqn{N_{Y2} \times 1}
#' @param T1 CTD temperature: \eqn{N_{T1} \times 1}
#' @param T2 therm temperature: \eqn{N_{T2} \times 1}
#' @param X covariates (whale, r2, r3): \eqn{N \times p}
#' @param date measurements date: \eqn{N \times 1}
#' @param coords measurements coordinate: \eqn{N \times 2}
#' @param oneY1 binary wether \code{Y1} exists or not: \eqn{N \times 1}
#' @param oneY2 binary wether \code{Y2} exists or not: \eqn{N \times 1}
#' @param oneT1 binary wether \code{T1} exists or not: \eqn{N \times 1}
#' @param oneT2 binary wether \code{T2} exists or not: \eqn{N \times 1}
#' @param calibration one of 
#'   \code{c("lm", "GP", "GPs", "cor", "corTime")}
#'   indicating a linear model (lm) in the calibration, a lm with
#'   a spatially varying intercept, a lm with spatially varying 
#'   intercept and slope, the same with coregionalization, and
#'   the same with time varying coefs. in ZOOP calibration,
#'   respectively
#' @param na,nb mean precision prior normal's
#' @param ga,gb shape rate prior gamma's
#' @param rangeTimesDmax the decay for the eta's is this value over dmax.
#'   Suggestion, keep between 3 and 30.
#' @param sd tuned sd in rwMetropolis
#' @param n.burnin,max.iter,n.thin,n.report features MCMC
#' @return Matrix with posterior samples from model parameters: 
#'   rows are iterations and columns are parameters
#'   
#' @author Jorge Castillo-Mateo
#' @export 
GibbsZoop <- function(Y1, Y2, T1, T2, X, date, coords, oneY1, oneY2, oneT1, oneT2,
                      calibration = c("lm", "GP", "GPs", "cor", "corTime"),
                      na = 0, nb = 1e-4, ga = 2, gb = 1, rangeTimesDmax = 3, sd = 0.0001,
                      n.burnin = 10000, max.iter = 10000, n.thin = 10, n.report = 1000) {
  
  calibration <- match.arg(calibration)
  if (length(rangeTimesDmax) == 1) {
    rangeTimesDmax <- rep(rangeTimesDmax, 2)
  } 
  
  ## PARAMETERS ##
  NY1 <- length(Y1)
  NY2 <- length(Y2)
  NT1 <- length(T1)
  NT2 <- length(T2)

  N <- nrow(X)
  n <- nrow(unique(coords))
  p <- ncol(X)
  year <- lubridate::year(date)
  T <- length(unique(year))
  
  firstYear <- min(year) - 1
  
  nt <- c() # number of sites per year
  for (tInd in 1:T) 
    nt[tInd] <- nrow(unique(coords[year == firstYear + tInd,]))
  
  Ent <- sum(nt)

  indY1 <- which(oneY1 == 1)
  indY2 <- which(oneY2 == 1)
  indT1 <- which(oneT1 == 1)
  indT2 <- which(oneT2 == 1)
  
  auxY1 <- rep(0, N)
  auxY2 <- rep(0, N)
  auxT1 <- rep(0, N)
  auxT2 <- rep(0, N)
  
  auxY1[indY1] <- Y1
  auxY2[indY2] <- Y2
  auxT1[indT1] <- T1
  auxT2[indT2] <- T2
  
  #NY1t <- tapply(Y1, year[indY1], length)
  #NY2t <- tapply(Y2, year[indY2], length)
  #NT1t <- tapply(T1, year[indT1], length)
  #NT2t <- tapply(T2, year[indT2], length)
  
  ## SPATIAL MATRIX ##
  allDist <- stats::dist(coords)
  dMax <- max(allDist)
  phi <- rangeTimesDmax[1] / dMax 
  R.phi             <- list()
  R.phi.inverse     <- array(dim = c(max(nt), max(nt), T))
  R.phi.inverse.sum <- rep(NA, T)
  for (tInd in 1:T) {
    distances <- stats::dist(unique(coords[year == firstYear + tInd,]))
    d <- matrix(0, nrow = nt[tInd], ncol = nt[tInd])
    d[lower.tri(d)] <- distances
    d <- d + t(d)
    R.phi[[tInd]] <- exp(- phi * d)
    R.phi.inverse[1:nt[tInd],1:nt[tInd],tInd] <- solve(R.phi[[tInd]])
    R.phi.inverse.sum[tInd] <- sum(R.phi.inverse[1:nt[tInd],1:nt[tInd],tInd])
  }
  
  ## MATRIX ETA ##
  my_fun <- function(x) match(data.frame(t(x)), data.frame(t(unique(x))))
  
  ti <- c()
  for (tInd in 1:T) 
    ti <- c(ti, my_fun(coords[year == firstYear + tInd,]) + sum(nt[0:(tInd-1)]))
  
  Xeta <- matrix(nrow = N, ncol = Ent)
  for (nInd in 1:Ent)
    Xeta[,nInd] <- ti == nInd
  
  Xeta <- 1 * Xeta
  
  ## ALPHA ##
  nT2 <- nrow(unique(coords[indT2,]))
  
  distT2 <- stats::dist(unique(coords[indT2,]))
  dT2 <- matrix(0, nrow = nT2, ncol = nT2)
  dT2[lower.tri(dT2)] <- distT2
  dT2 <- dT2 + t(dT2)
  
  sT2 <- my_fun(coords[indT2,])
  
  IT2 <- matrix(nrow = NT2, ncol = nT2)
  for (nInd in 1:nT2) IT2[,nInd] <- sT2 == nInd
  IT2 <- 1 * IT2
  
  ## LAMBDA ##
  nY2 <- nrow(unique(coords[indY2,]))
  
  distY2 <- stats::dist(unique(coords[indY2,]))
  dY2 <- matrix(0, nrow = nY2, ncol = nY2)
  dY2[lower.tri(dY2)] <- distY2
  dY2 <- dY2 + t(dY2)
  
  sY2 <- my_fun(coords[indY2,])
  
  IY2 <- matrix(nrow = NY2, ncol = nY2)
  for (nInd in 1:nY2) IY2[,nInd] <- sY2 == nInd
  IY2 <- 1 * IY2
  
  ## TEMP ##
  XTEMP <- cbind(sin(2*pi*lubridate::yday(date) / 365), cos(2*pi*lubridate::yday(date) / 365))

  distances <- stats::dist(unique(coords))
  d <- matrix(0, nrow = n, ncol = n)
  d[lower.tri(d)] <- distances
  d <- d + t(d)
  phi <- rangeTimesDmax[2] / dMax 
  RinvTEMP <- solve(exp(- phi * d))
  RsumTEMP <- sum(RinvTEMP)
  
  sTEMP <- my_fun(coords)
  
  ITEMP <- matrix(nrow = N, ncol = n)
  for (nInd in 1:n) ITEMP[,nInd] <- sTEMP == nInd
  ITEMP <- 1 * ITEMP
  
  Xt <- matrix(nrow = N, ncol = T)
  for (tInd in 1:T) Xt[,tInd] <- lubridate::year(date) == firstYear + tInd
  Xt <- 1 * Xt
  
  ## GIBBS SAMPLER ##
  if (calibration == "lm") {
    keep <- GibbsZoop0Cpp(
      auxY1, auxY2, auxT1, auxT2,
      X, Xeta, XTEMP, 
      R.phi.inverse.sum, RsumTEMP, R.phi.inverse, RinvTEMP,
      oneY1, oneY2, oneT1, oneT2, 
      indY1 - 1, indY2 - 1, indT1 - 1, indT2 - 1,
      ITEMP, sTEMP - 1, Xt,
      N, NY1, NY2, NT1, NT2, T, nt, Ent, n, p,
      na, nb, ga, gb, n.burnin, max.iter, n.thin, n.report)
    
    colnames(keep) <- 
      c(paste0("etat", rep(1:T, times = nt), "s", unlist(lapply(nt, seq_len))),
        paste0("etat", 1:T), paste0("beta", 1:p), "betaTemp", "beta0", "precY1", "precT1", "precEtas", "precEtat",
        "a0", "a1", "precT2", "l0", "l1", "precY2",
        paste0("temp", 1:N), paste0("phis", 1:n), "psi", "precphi", "precpsi", "sine", "cosine", "precTemp", paste0("psit", 1:T), "rhopsi")
  } else if (calibration == "GP") {
    keep <- GibbsZoop1Cpp(
      auxY1, auxY2, auxT1, auxT2,
      X, Xeta, XTEMP, R.phi.inverse.sum, RsumTEMP, R.phi.inverse, RinvTEMP, dMax,
      dY2, dT2, sd,
      oneY1, oneY2, oneT1, oneT2, 
      indY1 - 1, indY2 - 1, indT1 - 1, indT2 - 1,
      IY2, IT2, ITEMP, sY2 - 1, sT2 - 1, sTEMP - 1, Xt,
      N, NY1, NY2, NT1, NT2, T, nt, Ent, nT2, nY2, n, p,
      na, nb, ga, gb, n.burnin, max.iter, n.thin, n.report)
    
    colnames(keep) <- 
      c(paste0("etat", rep(1:T, times = nt), "s", unlist(lapply(nt, seq_len))),
        paste0("etat", 1:T), paste0("beta", 1:p), "betaTemp", "beta0", "precY1", "precT1", "precEtas", "precEtat",
        paste0("a0s", 1:nT2), "a1", "a0", "preca0", "decaya0", "precT2",
        paste0("l0s", 1:nY2), "l1", "l0", "precl0", "decayl0", "precY2",
        paste0("temp", 1:N), paste0("phis", 1:n), "psi", "precphi", "precpsi", "sine", "cosine", "precTemp", paste0("psit", 1:T), "rhopsi")
  } else if (calibration == "GPs") {
    stop("This model is not implemented. Not updated since Version 0.0.1")
    keep <- GibbsZoop2Cpp(
      auxY1, auxY2, auxT1, auxT2,
      X, Xeta, XTEMP, R.phi.inverse.sum, RsumTEMP, R.phi.inverse, RinvTEMP, dMax,
      dY2, dT2, sd,
      oneY1, oneY2, oneT1, oneT2, 
      indY1 - 1, indY2 - 1, indT1 - 1, indT2 - 1,
      IY2, IT2, ITEMP, sY2 - 1, sT2 - 1, sTEMP - 1, Xt,
      N, NY1, NY2, NT1, NT2, T, nt, Ent, nT2, nY2, n, p,
      na, nb, ga, gb, n.burnin, max.iter, n.thin, n.report)
    
    colnames(keep) <- 
      c(paste0("etat", rep(1:T, times = nt), "s", unlist(lapply(nt, seq_len))),
        paste0("etat", 1:T), paste0("beta", 1:p), "betaTemp", "beta0", "precY1", "precT1", "precEtas", "precEtat",
        paste0("a0s", 1:nT2), paste0("a1s", 1:nT2), "a0", "a1", "preca0", "preca1", "decaya0", "decaya1", "precT2", 
        paste0("l0s", 1:nY2), paste0("l1s", 1:nY2), "l0", "l1", "precl0", "precl1", "decayl0", "decayl1", "precY2", 
        paste0("temp", 1:N), paste0("phis", 1:n), "psi", "precphi", "precpsi", "sine", "cosine", "precTemp", paste0("psit", 1:T), "rhopsi")
  } else if (calibration == "cor")  {
    keep <- GibbsZoop3Cpp(
      auxY1, auxY2, auxT1, auxT2,
      X, Xeta, XTEMP, R.phi.inverse.sum, RsumTEMP, R.phi.inverse, RinvTEMP, dMax,
      dY2, dT2, sd,
      oneY1, oneY2, oneT1, oneT2, 
      indY1 - 1, indY2 - 1, indT1 - 1, indT2 - 1,
      IY2, IT2, ITEMP, sY2 - 1, sT2 - 1, sTEMP - 1, Xt,
      N, NY1, NY2, NT1, NT2, T, nt, Ent, nT2, nY2, n, p,
      na, nb, ga, gb, n.burnin, max.iter, n.thin, n.report)
    
    colnames(keep) <- 
      c(paste0("etat", rep(1:T, times = nt), "s", unlist(lapply(nt, seq_len))),
        paste0("etat", 1:T), paste0("beta", 1:p), "betaTemp", "beta0", "precY1", "precT1", "precEtas", "precEtat",
        paste0("a0s", 1:nT2), paste0("a1s", 1:nT2), "a0", "a1", "a11a", "a22a", "decaya0", "decaya1", "precT2", "a21a",
        paste0("l0s", 1:nY2), paste0("l1s", 1:nY2), "l0", "l1", "a11l", "a22l", "decayl0", "decayl1", "precY2", "a21l",
        paste0("temp", 1:N), paste0("phis", 1:n), "psi", "precphi", "precpsi", "sine", "cosine", "precTemp", paste0("psit", 1:T), "rhopsi")
  } else { #if (calibration == "corTime") 
    keep <- GibbsZoop4Cpp(
      auxY1, auxY2, auxT1, auxT2,
      X, Xeta, XTEMP, R.phi.inverse.sum, RsumTEMP, R.phi.inverse, RinvTEMP, dMax,
      dY2, dT2, sd,
      oneY1, oneY2, oneT1, oneT2, 
      indY1 - 1, indY2 - 1, indT1 - 1, indT2 - 1,
      IY2, IT2, ITEMP, sY2 - 1, sT2 - 1, sTEMP - 1, Xt,
      N, NY1, NY2, NT1, NT2, T, nt, Ent, nT2, nY2, n, p,
      na, nb, ga, gb, n.burnin, max.iter, n.thin, n.report)
    
    colnames(keep) <- 
      c(paste0("etat", rep(1:T, times = nt), "s", unlist(lapply(nt, seq_len))),
        paste0("etat", 1:T), paste0("beta", 1:p), "betaTemp", "beta0", "precY1", "precT1", "precEtas", "precEtat",
        paste0("a0s", 1:nT2), paste0("a1s", 1:nT2), "a0", "a1", "a11a", "a22a", "decaya0", "decaya1", "precT2", "a21a",
        paste0("l0s", 1:nY2), paste0("l1s", 1:nY2), "l0", "l1", "a11l", "a22l", "decayl0", "decayl1", "precY2", "a21l",
        paste0("l0t", 1:T), paste0("l1t", 1:T), "rhol0", "rhol1", "v11", "v22", "v21",
        paste0("temp", 1:N), paste0("phis", 1:n), "psi", "precphi", "precpsi", "sine", "cosine", "precTemp", paste0("psit", 1:T), "rhopsi")
  }
  
  return(keep)
}