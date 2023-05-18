#' Zoop prediction
#' 
#' @importFrom lubridate yday
#' @importFrom lubridate year
#' @importFrom stats dist
#' 
#' @description This function predicts Y1 Y2 T1 T2 temp and log ZOOP
#' 
#' @param model model fitted with the function \code{\link{GibbsZoop}}
#' @param X covariates (whale, r2, r3): \eqn{1 \times p}
#' @param dateNew new date
#' @param coordsNew new coords
#' @param dateData dates used to fit the model
#' @param coordsData coords used to fit the model
#' @param indY2 positions with Y2 not NA: \eqn{NY2 \times 1}
#' @param indT2 positions with T2 not NA: \eqn{NT2 \times 1}
#' @param calibration one of 
#'   \code{c("lm", "GP", "GPs", "cor", "corTime")}
#'   indicating a linear model (lm) in the calibration, a lm with
#'   a spatially varying intercept, a lm with spatially varying 
#'   intercept and slope, the same with coregionalization, and
#'   the same with time varying coefs. in ZOOP calibration,
#'   respectively
#' @param rangeTimesDmax the same as used to fit the model
#' @examples
#' #a <- predict(model, 
#' #  X = c(1,0,0), 
#' #  dateNew = as.Date("2005-01-31"), coordsNew = c(-30000, 0),
#' #  dateData = date, coordsData = coords)
#' 
#' #b <- predict(model, 
#' #  X = c(1,0,0), 
#' #  dateNew = as.Date("2005-01-31"), coordsNew = c(-30000, 0),
#' #  dateData = date, coordsData = coords,
#' #  indY2, indT2, "GP")
#' 
#' @return matrix rows are simuations columns are:
#'   Y1 Y2 T1 T2 temp and log ZOOP
#' @export 
predictZoop <- function(model, X, 
                    dateNew, coordsNew, dateData, coordsData, 
                    indY2 = NULL, indT2 = NULL, 
                    calibration = c("lm", "GP", "GPs", "cor", "corTime"),
                    rangeTimesDmax = 3) {
  
  calibration <- match.arg(calibration)
  if (length(rangeTimesDmax) == 1) {
    rangeTimesDmax <- rep(rangeTimesDmax, 2)
  } 
  
  B <- nrow(model)
  p <- length(X)
  ydayY <- lubridate::yday(dateNew)
  yearY <- lubridate::year(dateNew)
  yearData <- lubridate::year(dateData)
  firstYear <- min(yearData) - 1
  t <- yearY - firstYear
  
  Ynew <- matrix(nrow = B, ncol = 6) ## Y1 Y2 T1 T2 temp ZOOP
  colnames(Ynew) <- c("Y1", "Y2", "T1", "T2", "temp", "ZOOP")
  
  # SIMULACION CAMPOS 
  nt <- nrow(unique(coordsData[yearData == yearY,]))
  
  n  <- nrow(unique(coordsData))

  ## CAMPO ETA ##
  allDist <- stats::dist(coordsData)
  dMax <- max(allDist)
  phi <- rangeTimesDmax[1] / dMax 
  
  distances <- stats::dist(unique(coordsData[yearData == yearY,]))
  d <- matrix(0, nrow = nt, ncol = nt)
  d[lower.tri(d)] <- distances
  d <- d + t(d)
  Rinv <- solve(exp(- phi * d))

  # Sigma, Sigma00, Sigmai0
  distances <- dist(rbind(coordsNew, unique(coordsData[yearData == yearY,])))
  d <- matrix(0, nrow = nt + 1, ncol = nt + 1)
  d[lower.tri(d)] <- distances
  d <- d + t(d)
  
  if (any(d[1, 2:(nt + 1)] == 0)) {
    SIGMA_pre <- matrix(0, nrow = 1, ncol = nt)
    SIGMA_pre[1,d[1, 2:(nt + 1)] == 0] <- 1
    SIGMA_cond <- 0
  } else {
    SIGMA    <- exp(- phi * d)
    Sigma00  <- SIGMA[1, 1]
    Sigmai0T <- SIGMA[1, 2:(nt + 1)]
    
    SIGMA_pre  <- Sigmai0T %*% Rinv
    SIGMA_cond <- Sigma00 - SIGMA_pre %*% Sigmai0T
  }
  
  ## CAMPO TEMP ##
  distances <- stats::dist(unique(coordsData))
  d <- matrix(0, nrow = n, ncol = n)
  d[lower.tri(d)] <- distances
  d <- d + t(d)
  phi <- rangeTimesDmax[2] / dMax 
  RinvTEMP <- solve(exp(- phi * d))
  
  # Sigma, Sigma00, Sigmai0
  distances <- dist(rbind(coordsNew, unique(coordsData)))
  d <- matrix(0, nrow = n + 1, ncol = n + 1)
  d[lower.tri(d)] <- distances
  d <- d + t(d)
  
  if (any(d[1, 2:(n + 1)] == 0)) {
    SIGMA_preTEMP <- matrix(0, nrow = 1, ncol = n)
    SIGMA_preTEMP[1,d[1, 2:(n + 1)] == 0] <- 1
    SIGMA_condTEMP <- 0
  } else {
    SIGMA    <- exp(- phi * d)
    Sigma00  <- SIGMA[1, 1]
    Sigmai0T <- SIGMA[1, 2:(n + 1)]
    
    SIGMA_preTEMP  <- Sigmai0T %*% RinvTEMP
    SIGMA_condTEMP <- Sigma00 - SIGMA_preTEMP %*% Sigmai0T
  }
  
  if (calibration == "lm") {
    
    ## ITERATIONS
    for (b in 1:B) {
      print(b)
      campoEspacialTEMP <- stats::rnorm(1,
                                 SIGMA_preTEMP %*% matrix(model[b, paste0("phis", 1:n)], ncol = 1),
                                 sqrt(SIGMA_condTEMP / model[b,"precphi"]))
      
      campoEspacialEtats <- stats::rnorm(1,
                                  model[b, paste0("etat",t)] + SIGMA_pre %*% matrix(model[b, paste0("etat", t, "s", 1:nt)] - model[b, paste0("etat",t)], ncol = 1),
                                  sqrt(SIGMA_cond / model[b,"precEtas"]))
      
      Ynew[b,"temp"] <- stats::rnorm(1, 
                              campoEspacialTEMP + model[b,"sine"] * sin(2*pi*ydayY/365) + model[b,"cosine"] * cos(2*pi*ydayY/365) + model[b,paste0("psit",t)], 
                              1 / sqrt(model[b,"precTemp"]))
      
      Ynew[b,"ZOOP"] <- campoEspacialEtats + X %*% model[b, c(paste0("beta",1:p))] + model[b,"betaTemp"] * Ynew[b,"temp"]
      
      Ynew[b,"Y1"] <- stats::rnorm(1, Ynew[b,"ZOOP"], 1 / sqrt(model[b,"precY1"]))
      Ynew[b,"Y2"] <- stats::rnorm(1, model[b,"l0"] + model[b,"l1"] * Ynew[b,"ZOOP"], 1 / sqrt(model[b,"precY2"]))
      
      Ynew[b,"T1"] <- stats::rnorm(1, Ynew[b,"temp"], 1 / sqrt(model[b,"precT1"]))
      Ynew[b,"T2"] <- stats::rnorm(1, model[b,"a0"] + model[b,"a1"] * Ynew[b,"temp"], 1 / sqrt(model[b,"precT2"]))
    }
    
  } else if (calibration == "GP") {
    
    ## CAMPO a0 ##
    nT2 <- nrow(unique(coordsData[indT2,]))
    
    distT2 <- stats::dist(unique(coordsData[indT2,]))
    dT2 <- matrix(0, nrow = nT2, ncol = nT2)
    dT2[lower.tri(dT2)] <- distT2
    dT2 <- dT2 + t(dT2)
    
    distA <- dist(rbind(coordsNew, unique(coordsData[indT2,])))
    dA <- matrix(0, nrow = nT2 + 1, ncol = nT2 + 1)
    dA[lower.tri(dA)] <- distA
    dA <- dA + t(dA)
    
    ## CAMPO l0 ##
    nY2 <- nrow(unique(coordsData[indY2,]))
    
    distY2 <- stats::dist(unique(coordsData[indY2,]))
    dY2 <- matrix(0, nrow = nY2, ncol = nY2)
    dY2[lower.tri(dY2)] <- distY2
    dY2 <- dY2 + t(dY2)
    
    distL <- dist(rbind(coordsNew, unique(coordsData[indY2,])))
    dL <- matrix(0, nrow = nY2 + 1, ncol = nY2 + 1)
    dL[lower.tri(dL)] <- distL
    dL <- dL + t(dL)
    
    ## ITERATIONS
    for (b in 1:B) {
      print(b)
      campoEspacialTEMP <- stats::rnorm(1,
        SIGMA_preTEMP %*% matrix(model[b, paste0("phis", 1:n)], ncol = 1),
        sqrt(SIGMA_condTEMP / model[b,"precphi"]))
      
      campoEspacialEtats <- stats::rnorm(1,
        model[b, paste0("etat",t)] + SIGMA_pre %*% matrix(model[b, paste0("etat", t, "s", 1:nt)] - model[b, paste0("etat",t)], ncol = 1),
        sqrt(SIGMA_cond / model[b,"precEtas"]))
      
      ### a0
      if (any(dA[1, 2:(nT2 + 1)] == 0)) {
        SIGMA_prea0 <- matrix(0, nrow = 1, ncol = nT2)
        SIGMA_prea0[1,dA[1, 2:(nT2 + 1)] == 0] <- 1
        SIGMA_conda0 <- 0
      } else {
        Rinva0 <- solve(exp(- model[b,"decaya0"] * dT2))
        
        SIGMA    <- exp(- model[b,"decaya0"] * dA)
        Sigma00  <- SIGMA[1, 1]
        Sigmai0T <- SIGMA[1, 2:(nT2 + 1)]
        
        SIGMA_prea0  <- Sigmai0T %*% Rinva0
        SIGMA_conda0 <- Sigma00 - SIGMA_prea0 %*% Sigmai0T
      }
      
      campoEspaciala0 <- stats::rnorm(1,
        model[b,"a0"] + SIGMA_prea0 %*% matrix(model[b, paste0("a0s", 1:nT2)] - model[b,"a0"], ncol = 1),
        sqrt(SIGMA_conda0 / model[b,"preca0"]))
      
      ### l0
      if (any(dL[1, 2:(nY2 + 1)] == 0)) {
        SIGMA_prel0 <- matrix(0, nrow = 1, ncol = nY2)
        SIGMA_prel0[1,dL[1, 2:(nY2 + 1)] == 0] <- 1
        SIGMA_condl0 <- 0
      } else {
        Rinvl0 <- solve(exp(- model[b,"decayl0"] * dY2))
        
        SIGMA    <- exp(- model[b,"decayl0"] * dL)
        Sigma00  <- SIGMA[1, 1]
        Sigmai0T <- SIGMA[1, 2:(nY2 + 1)]
        
        SIGMA_prel0  <- Sigmai0T %*% Rinvl0
        SIGMA_condl0 <- Sigma00 - SIGMA_prel0 %*% Sigmai0T
      }
      
      campoEspaciall0 <- stats::rnorm(1,
        model[b,"l0"] + SIGMA_prel0 %*% matrix(model[b, paste0("l0s", 1:nY2)] - model[b,"l0"], ncol = 1),
        sqrt(SIGMA_condl0 / model[b,"precl0"]))
      
      ## continue
      Ynew[b,"temp"] <- stats::rnorm(1, 
                              campoEspacialTEMP + model[b,"sine"] * sin(2*pi*ydayY/365) + model[b,"cosine"] * cos(2*pi*ydayY/365) + model[b,paste0("psit",t)], 
                              1 / sqrt(model[b,"precTemp"]))
      
      Ynew[b,"ZOOP"] <- campoEspacialEtats + X %*% model[b, c(paste0("beta",1:p))] + model[b,"betaTemp"] * Ynew[b,"temp"]
      
      Ynew[b,"Y1"] <- stats::rnorm(1, Ynew[b,"ZOOP"], 1 / sqrt(model[b,"precY1"]))
      Ynew[b,"Y2"] <- stats::rnorm(1, campoEspaciall0 + model[b,"l1"] * Ynew[b,"ZOOP"], 1 / sqrt(model[b,"precY2"]))
      
      Ynew[b,"T1"] <- stats::rnorm(1, Ynew[b,"temp"], 1 / sqrt(model[b,"precT1"]))
      Ynew[b,"T2"] <- stats::rnorm(1, campoEspaciala0 + model[b,"a1"] * Ynew[b,"temp"], 1 / sqrt(model[b,"precT2"]))
    }
    
  } else if (calibration == "GPs") {
    
    ## CAMPO a's ##
    nT2 <- nrow(unique(coordsData[indT2,]))
    
    distT2 <- stats::dist(unique(coordsData[indT2,]))
    dT2 <- matrix(0, nrow = nT2, ncol = nT2)
    dT2[lower.tri(dT2)] <- distT2
    dT2 <- dT2 + t(dT2)
    
    distA <- dist(rbind(coordsNew, unique(coordsData[indT2,])))
    dA <- matrix(0, nrow = nT2 + 1, ncol = nT2 + 1)
    dA[lower.tri(dA)] <- distA
    dA <- dA + t(dA)
    
    ## CAMPO l's ##
    nY2 <- nrow(unique(coordsData[indY2,]))
    
    distY2 <- stats::dist(unique(coordsData[indY2,]))
    dY2 <- matrix(0, nrow = nY2, ncol = nY2)
    dY2[lower.tri(dY2)] <- distY2
    dY2 <- dY2 + t(dY2)
    
    distL <- dist(rbind(coordsNew, unique(coordsData[indY2,])))
    dL <- matrix(0, nrow = nY2 + 1, ncol = nY2 + 1)
    dL[lower.tri(dL)] <- distL
    dL <- dL + t(dL)
    
    ## ITERATIONS
    for (b in 1:B) {
      print(b)
      campoEspacialTEMP <- stats::rnorm(1,
                                        SIGMA_preTEMP %*% matrix(model[b, paste0("phis", 1:n)], ncol = 1),
                                        sqrt(SIGMA_condTEMP / model[b,"precphi"]))
      
      campoEspacialEtats <- stats::rnorm(1,
                                         model[b, paste0("etat",t)] + SIGMA_pre %*% matrix(model[b, paste0("etat", t, "s", 1:nt)] - model[b, paste0("etat",t)], ncol = 1),
                                         sqrt(SIGMA_cond / model[b,"precEtas"]))
      
      ### a's
      if (any(dA[1, 2:(nT2 + 1)] == 0)) {
        SIGMA_prea0 <- SIGMA_prea1 <- matrix(0, nrow = 1, ncol = nT2)
        SIGMA_prea0[1,dA[1, 2:(nT2 + 1)] == 0] <- SIGMA_prea1[1,dA[1, 2:(nT2 + 1)] == 0] <- 1
        SIGMA_conda0 <- SIGMA_conda1 <- 0
      } else {
        #a0
        Rinva0 <- solve(exp(- model[b,"decaya0"] * dT2))
        
        SIGMA    <- exp(- model[b,"decaya0"] * dA)
        Sigma00  <- SIGMA[1, 1]
        Sigmai0T <- SIGMA[1, 2:(nT2 + 1)]
        
        SIGMA_prea0  <- Sigmai0T %*% Rinva0
        SIGMA_conda0 <- Sigma00 - SIGMA_prea0 %*% Sigmai0T
        #a1
        Rinva1 <- solve(exp(- model[b,"decaya1"] * dT2))
        
        SIGMA    <- exp(- model[b,"decaya1"] * dA)
        Sigma00  <- SIGMA[1, 1]
        Sigmai0T <- SIGMA[1, 2:(nT2 + 1)]
        
        SIGMA_prea1  <- Sigmai0T %*% Rinva1
        SIGMA_conda1 <- Sigma00 - SIGMA_prea1 %*% Sigmai0T
      }
      
      campoEspaciala0 <- stats::rnorm(1,
        model[b,"a0"] + SIGMA_prea0 %*% matrix(model[b, paste0("a0s", 1:nT2)] - model[b,"a0"], ncol = 1),
        sqrt(SIGMA_conda0 / model[b,"preca0"]))
      
      campoEspaciala1 <- stats::rnorm(1,
        model[b,"a1"] + SIGMA_prea1 %*% matrix(model[b, paste0("a1s", 1:nT2)] - model[b,"a1"], ncol = 1),
        sqrt(SIGMA_conda1 / model[b,"preca1"]))
      
      ### l's
      if (any(dL[1, 2:(nY2 + 1)] == 0)) {
        SIGMA_prel0 <- SIGMA_prel1 <- matrix(0, nrow = 1, ncol = nY2)
        SIGMA_prel0[1,dL[1, 2:(nY2 + 1)] == 0] <- SIGMA_prel1[1,dL[1, 2:(nY2 + 1)] == 0] <- 1
        SIGMA_condl0 <- SIGMA_condl1 <- 0
      } else {
        #l0
        Rinvl0 <- solve(exp(- model[b,"decayl0"] * dY2))
        
        SIGMA    <- exp(- model[b,"decayl0"] * dL)
        Sigma00  <- SIGMA[1, 1]
        Sigmai0T <- SIGMA[1, 2:(nY2 + 1)]
        
        SIGMA_prel0  <- Sigmai0T %*% Rinvl0
        SIGMA_condl0 <- Sigma00 - SIGMA_prel0 %*% Sigmai0T
        #l1
        Rinvl1 <- solve(exp(- model[b,"decayl1"] * dY2))
        
        SIGMA    <- exp(- model[b,"decayl1"] * dL)
        Sigma00  <- SIGMA[1, 1]
        Sigmai0T <- SIGMA[1, 2:(nY2 + 1)]
        
        SIGMA_prel1  <- Sigmai0T %*% Rinvl1
        SIGMA_condl1 <- Sigma00 - SIGMA_prel1 %*% Sigmai0T
      }
      
      campoEspaciall0 <- stats::rnorm(1,
        model[b,"l0"] + SIGMA_prel0 %*% matrix(model[b, paste0("l0s", 1:nY2)] - model[b,"l0"], ncol = 1),
        sqrt(SIGMA_condl0 / model[b,"precl0"]))
      
      campoEspaciall1 <- stats::rnorm(1,
        model[b,"l1"] + SIGMA_prel1 %*% matrix(model[b, paste0("l1s", 1:nY2)] - model[b,"l1"], ncol = 1),
        sqrt(SIGMA_condl1 / model[b,"precl1"]))
      
      ## continue
      Ynew[b,"temp"] <- stats::rnorm(1, 
                                     campoEspacialTEMP + model[b,"sine"] * sin(2*pi*ydayY/365) + model[b,"cosine"] * cos(2*pi*ydayY/365) + model[b,paste0("psit",t)], 
                                     1 / sqrt(model[b,"precTemp"]))
      
      Ynew[b,"ZOOP"] <- campoEspacialEtats + X %*% model[b, c(paste0("beta",1:p))] + model[b,"betaTemp"] * Ynew[b,"temp"]
      
      Ynew[b,"Y1"] <- stats::rnorm(1, Ynew[b,"ZOOP"], 1 / sqrt(model[b,"precY1"]))
      Ynew[b,"Y2"] <- stats::rnorm(1, campoEspaciall0 + campoEspaciall1 * Ynew[b,"ZOOP"], 1 / sqrt(model[b,"precY2"]))
      
      Ynew[b,"T1"] <- stats::rnorm(1, Ynew[b,"temp"], 1 / sqrt(model[b,"precT1"]))
      Ynew[b,"T2"] <- stats::rnorm(1, campoEspaciala0 + campoEspaciala1 * Ynew[b,"temp"], 1 / sqrt(model[b,"precT2"]))
    }
    
  } else if (calibration == "cor") {
    
    ## CAMPO a's ##
    nT2 <- nrow(unique(coordsData[indT2,]))
    
    distT2 <- stats::dist(unique(coordsData[indT2,]))
    dT2 <- matrix(0, nrow = nT2, ncol = nT2)
    dT2[lower.tri(dT2)] <- distT2
    dT2 <- dT2 + t(dT2)
    
    distA <- dist(rbind(coordsNew, unique(coordsData[indT2,])))
    dA <- matrix(0, nrow = nT2 + 1, ncol = nT2 + 1)
    dA[lower.tri(dA)] <- distA
    dA <- dA + t(dA)
    
    ## CAMPO l's ##
    nY2 <- nrow(unique(coordsData[indY2,]))
    
    distY2 <- stats::dist(unique(coordsData[indY2,]))
    dY2 <- matrix(0, nrow = nY2, ncol = nY2)
    dY2[lower.tri(dY2)] <- distY2
    dY2 <- dY2 + t(dY2)
    
    distL <- dist(rbind(coordsNew, unique(coordsData[indY2,])))
    dL <- matrix(0, nrow = nY2 + 1, ncol = nY2 + 1)
    dL[lower.tri(dL)] <- distL
    dL <- dL + t(dL)
    
    ## ITERATIONS
    for (b in 1:B) {
      print(b)
      campoEspacialTEMP <- stats::rnorm(1,
        SIGMA_preTEMP %*% matrix(model[b, paste0("phis", 1:n)], ncol = 1),
        sqrt(SIGMA_condTEMP / model[b,"precphi"]))
      
      campoEspacialEtats <- stats::rnorm(1,
        model[b, paste0("etat",t)] + SIGMA_pre %*% matrix(model[b, paste0("etat", t, "s", 1:nt)] - model[b, paste0("etat",t)], ncol = 1),
        sqrt(SIGMA_cond / model[b,"precEtas"]))
      
      ### a's
      if (any(dA[1, 2:(nT2 + 1)] == 0)) {
        SIGMA_prea0 <- SIGMA_prea1 <- matrix(0, nrow = 1, ncol = nT2)
        SIGMA_prea0[1,dA[1, 2:(nT2 + 1)] == 0] <- SIGMA_prea1[1,dA[1, 2:(nT2 + 1)] == 0] <- 1
        SIGMA_conda0 <- SIGMA_conda1 <- 0
      } else {
        #a0
        Rinva0 <- solve(exp(- model[b,"decaya0"] * dT2))
        
        SIGMA    <- exp(- model[b,"decaya0"] * dA)
        Sigma00  <- SIGMA[1, 1]
        Sigmai0T <- SIGMA[1, 2:(nT2 + 1)]
        
        SIGMA_prea0  <- Sigmai0T %*% Rinva0
        SIGMA_conda0 <- Sigma00 - SIGMA_prea0 %*% Sigmai0T
        #a1
        Rinva1 <- solve(exp(- model[b,"decaya1"] * dT2))
        
        SIGMA    <- exp(- model[b,"decaya1"] * dA)
        Sigma00  <- SIGMA[1, 1]
        Sigmai0T <- SIGMA[1, 2:(nT2 + 1)]
        
        SIGMA_prea1  <- Sigmai0T %*% Rinva1
        SIGMA_conda1 <- Sigma00 - SIGMA_prea1 %*% Sigmai0T
      }
      
      campoEspacialv0 <- stats::rnorm(1,
        SIGMA_prea0 %*% matrix(model[b, paste0("a0s", 1:nT2)], ncol = 1),
        sqrt(SIGMA_conda0))
      
      campoEspacialv1 <- stats::rnorm(1,
        SIGMA_prea1 %*% matrix(model[b, paste0("a1s", 1:nT2)], ncol = 1),
        sqrt(SIGMA_conda1))
      
      campoEspaciala0 <- model[b,"a0"] + 
        model[b,"a11a"] * campoEspacialv0
      
      campoEspaciala1 <- model[b,"a1"] + 
        model[b,"a21a"] * campoEspacialv0 +
        model[b,"a22a"] * campoEspacialv1
      
      ### l's
      if (any(dL[1, 2:(nY2 + 1)] == 0)) {
        SIGMA_prel0 <- SIGMA_prel1 <- matrix(0, nrow = 1, ncol = nY2)
        SIGMA_prel0[1,dL[1, 2:(nY2 + 1)] == 0] <- SIGMA_prel1[1,dL[1, 2:(nY2 + 1)] == 0] <- 1
        SIGMA_condl0 <- SIGMA_condl1 <- 0
      } else {
        #l0
        Rinvl0 <- solve(exp(- model[b,"decayl0"] * dY2))
        
        SIGMA    <- exp(- model[b,"decayl0"] * dL)
        Sigma00  <- SIGMA[1, 1]
        Sigmai0T <- SIGMA[1, 2:(nY2 + 1)]
        
        SIGMA_prel0  <- Sigmai0T %*% Rinvl0
        SIGMA_condl0 <- Sigma00 - SIGMA_prel0 %*% Sigmai0T
        #l1
        Rinvl1 <- solve(exp(- model[b,"decayl1"] * dY2))
        
        SIGMA    <- exp(- model[b,"decayl1"] * dL)
        Sigma00  <- SIGMA[1, 1]
        Sigmai0T <- SIGMA[1, 2:(nY2 + 1)]
        
        SIGMA_prel1  <- Sigmai0T %*% Rinvl1
        SIGMA_condl1 <- Sigma00 - SIGMA_prel1 %*% Sigmai0T
      }
      
      campoEspacialw0 <- stats::rnorm(1,
        SIGMA_prel0 %*% matrix(model[b, paste0("l0s", 1:nY2)], ncol = 1),
        sqrt(SIGMA_condl0))
      
      campoEspacialw1 <- stats::rnorm(1,
        SIGMA_prel1 %*% matrix(model[b, paste0("l1s", 1:nY2)], ncol = 1),
        sqrt(SIGMA_condl1))
      
      campoEspaciall0 <- model[b,"l0"] + 
        model[b,"a11l"] * campoEspacialw0
      
      campoEspaciall1 <- model[b,"l1"] + 
        model[b,"a21l"] * campoEspacialw0 +
        model[b,"a22l"] * campoEspacialw1
      
      ## continue
      Ynew[b,"temp"] <- stats::rnorm(1, 
        campoEspacialTEMP + model[b,"sine"] * sin(2*pi*ydayY/365) + model[b,"cosine"] * cos(2*pi*ydayY/365) + model[b,paste0("psit",t)], 
        1 / sqrt(model[b,"precTemp"]))
      
      Ynew[b,"ZOOP"] <- campoEspacialEtats + X %*% model[b, c(paste0("beta",1:p))] + model[b,"betaTemp"] * Ynew[b,"temp"]
      
      Ynew[b,"Y1"] <- stats::rnorm(1, Ynew[b,"ZOOP"], 1 / sqrt(model[b,"precY1"]))
      Ynew[b,"Y2"] <- stats::rnorm(1, campoEspaciall0 + campoEspaciall1 * Ynew[b,"ZOOP"], 1 / sqrt(model[b,"precY2"]))
      
      Ynew[b,"T1"] <- stats::rnorm(1, Ynew[b,"temp"], 1 / sqrt(model[b,"precT1"]))
      Ynew[b,"T2"] <- stats::rnorm(1, campoEspaciala0 + campoEspaciala1 * Ynew[b,"temp"], 1 / sqrt(model[b,"precT2"]))
    }
    
  } else {# if (calibration == "corTime")
    
    ## CAMPO a's ##
    nT2 <- nrow(unique(coordsData[indT2,]))
    
    distT2 <- stats::dist(unique(coordsData[indT2,]))
    dT2 <- matrix(0, nrow = nT2, ncol = nT2)
    dT2[lower.tri(dT2)] <- distT2
    dT2 <- dT2 + t(dT2)
    
    distA <- dist(rbind(coordsNew, unique(coordsData[indT2,])))
    dA <- matrix(0, nrow = nT2 + 1, ncol = nT2 + 1)
    dA[lower.tri(dA)] <- distA
    dA <- dA + t(dA)
    
    ## CAMPO l's ##
    nY2 <- nrow(unique(coordsData[indY2,]))
    
    distY2 <- stats::dist(unique(coordsData[indY2,]))
    dY2 <- matrix(0, nrow = nY2, ncol = nY2)
    dY2[lower.tri(dY2)] <- distY2
    dY2 <- dY2 + t(dY2)
    
    distL <- dist(rbind(coordsNew, unique(coordsData[indY2,])))
    dL <- matrix(0, nrow = nY2 + 1, ncol = nY2 + 1)
    dL[lower.tri(dL)] <- distL
    dL <- dL + t(dL)
    
    ## ITERATIONS
    for (b in 1:B) {
      print(b)
      campoEspacialTEMP <- stats::rnorm(1,
                                        SIGMA_preTEMP %*% matrix(model[b, paste0("phis", 1:n)], ncol = 1),
                                        sqrt(SIGMA_condTEMP / model[b,"precphi"]))
      
      campoEspacialEtats <- stats::rnorm(1,
                                         model[b, paste0("etat",t)] + SIGMA_pre %*% matrix(model[b, paste0("etat", t, "s", 1:nt)] - model[b, paste0("etat",t)], ncol = 1),
                                         sqrt(SIGMA_cond / model[b,"precEtas"]))
      
      ### a's
      if (any(dA[1, 2:(nT2 + 1)] == 0)) {
        SIGMA_prea0 <- SIGMA_prea1 <- matrix(0, nrow = 1, ncol = nT2)
        SIGMA_prea0[1,dA[1, 2:(nT2 + 1)] == 0] <- SIGMA_prea1[1,dA[1, 2:(nT2 + 1)] == 0] <- 1
        SIGMA_conda0 <- SIGMA_conda1 <- 0
      } else {
        #a0
        Rinva0 <- solve(exp(- model[b,"decaya0"] * dT2))
        
        SIGMA    <- exp(- model[b,"decaya0"] * dA)
        Sigma00  <- SIGMA[1, 1]
        Sigmai0T <- SIGMA[1, 2:(nT2 + 1)]
        
        SIGMA_prea0  <- Sigmai0T %*% Rinva0
        SIGMA_conda0 <- Sigma00 - SIGMA_prea0 %*% Sigmai0T
        #a1
        Rinva1 <- solve(exp(- model[b,"decaya1"] * dT2))
        
        SIGMA    <- exp(- model[b,"decaya1"] * dA)
        Sigma00  <- SIGMA[1, 1]
        Sigmai0T <- SIGMA[1, 2:(nT2 + 1)]
        
        SIGMA_prea1  <- Sigmai0T %*% Rinva1
        SIGMA_conda1 <- Sigma00 - SIGMA_prea1 %*% Sigmai0T
      }
      
      campoEspacialv0 <- stats::rnorm(1,
                                      SIGMA_prea0 %*% matrix(model[b, paste0("a0s", 1:nT2)], ncol = 1),
                                      sqrt(SIGMA_conda0))
      
      campoEspacialv1 <- stats::rnorm(1,
                                      SIGMA_prea1 %*% matrix(model[b, paste0("a1s", 1:nT2)], ncol = 1),
                                      sqrt(SIGMA_conda1))
      
      campoEspaciala0 <- model[b,"a0"] + 
        model[b,"a11a"] * campoEspacialv0
      
      campoEspaciala1 <- model[b,"a1"] + 
        model[b,"a21a"] * campoEspacialv0 +
        model[b,"a22a"] * campoEspacialv1
      
      ### l's
      if (any(dL[1, 2:(nY2 + 1)] == 0)) {
        SIGMA_prel0 <- SIGMA_prel1 <- matrix(0, nrow = 1, ncol = nY2)
        SIGMA_prel0[1,dL[1, 2:(nY2 + 1)] == 0] <- SIGMA_prel1[1,dL[1, 2:(nY2 + 1)] == 0] <- 1
        SIGMA_condl0 <- SIGMA_condl1 <- 0
      } else {
        #l0
        Rinvl0 <- solve(exp(- model[b,"decayl0"] * dY2))
        
        SIGMA    <- exp(- model[b,"decayl0"] * dL)
        Sigma00  <- SIGMA[1, 1]
        Sigmai0T <- SIGMA[1, 2:(nY2 + 1)]
        
        SIGMA_prel0  <- Sigmai0T %*% Rinvl0
        SIGMA_condl0 <- Sigma00 - SIGMA_prel0 %*% Sigmai0T
        #l1
        Rinvl1 <- solve(exp(- model[b,"decayl1"] * dY2))
        
        SIGMA    <- exp(- model[b,"decayl1"] * dL)
        Sigma00  <- SIGMA[1, 1]
        Sigmai0T <- SIGMA[1, 2:(nY2 + 1)]
        
        SIGMA_prel1  <- Sigmai0T %*% Rinvl1
        SIGMA_condl1 <- Sigma00 - SIGMA_prel1 %*% Sigmai0T
      }
      
      campoEspacialw0 <- stats::rnorm(1,
                                      SIGMA_prel0 %*% matrix(model[b, paste0("l0s", 1:nY2)], ncol = 1),
                                      sqrt(SIGMA_condl0))
      
      campoEspacialw1 <- stats::rnorm(1,
                                      SIGMA_prel1 %*% matrix(model[b, paste0("l1s", 1:nY2)], ncol = 1),
                                      sqrt(SIGMA_condl1))
      
      campoEspaciall0 <- model[b,paste0("l0t",t)] + 
        model[b,"a11l"] * campoEspacialw0
      
      campoEspaciall1 <- model[b,paste0("l1t",t)] + 
        model[b,"a21l"] * campoEspacialw0 +
        model[b,"a22l"] * campoEspacialw1
      
      ## continue
      Ynew[b,"temp"] <- stats::rnorm(1, 
                                     campoEspacialTEMP + model[b,"sine"] * sin(2*pi*ydayY/365) + model[b,"cosine"] * cos(2*pi*ydayY/365) + model[b,paste0("psit",t)], 
                                     1 / sqrt(model[b,"precTemp"]))
      
      Ynew[b,"ZOOP"] <- campoEspacialEtats + X %*% model[b, c(paste0("beta",1:p))] + model[b,"betaTemp"] * Ynew[b,"temp"]
      
      Ynew[b,"Y1"] <- stats::rnorm(1, Ynew[b,"ZOOP"], 1 / sqrt(model[b,"precY1"]))
      Ynew[b,"Y2"] <- stats::rnorm(1, campoEspaciall0 + campoEspaciall1 * Ynew[b,"ZOOP"], 1 / sqrt(model[b,"precY2"]))
      
      Ynew[b,"T1"] <- stats::rnorm(1, Ynew[b,"temp"], 1 / sqrt(model[b,"precT1"]))
      Ynew[b,"T2"] <- stats::rnorm(1, campoEspaciala0 + campoEspaciala1 * Ynew[b,"temp"], 1 / sqrt(model[b,"precT2"]))
    }
  }
  
  return(Ynew)
}