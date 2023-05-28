#' Posterior Mean Surfaces of the GP's
#' 
#' @importFrom lubridate year
#' @importFrom scales muted
#' @importFrom stats dist
#' 
#' @description This function saves the map of the posterior mean of the
#'   GP's in \code{\link{GibbsZoop}} when \code{calibration = "corTime"}
#' 
#' @param model parameters of \code{\link{GibbsZoop}}
#' @param coords \eqn{N \times 2}
#' @param date \eqn{N \times 1}
#' @param oneY2 binary wether \code{Y2} exists or not: \eqn{N \times 1}
#' @param oneT2 binary wether \code{T2} exists or not: \eqn{N \times 1}
#' @param process one of \code{"eta"} (zoop GP's), \code{"alpha"} (temp 
#'   calibration GP's), \code{"lambda"} (zoop calibration GP's), or 
#'   \code{"phi"} (temp GP)
#' @param pred.coords coords where to predict
#' @param background background of CCB 
#' @param xlim,ylim,zlim x y z limits in the map
#' @return Saved png
#' 
#' @author Jorge Castillo-Mateo
#' @export 
meanGP.GibbsZoop <- function(model, coords, date = NULL, 
  oneY2 = NULL, oneT2 = NULL, 
  process = c("eta", "alpha", "lambda", "phi"), 
  pred.coords, background, 
  xlim = c(-61771.464, 1228.536), ylim = c(-34392.144, 13607.86), zlim = NULL) {
  
  process <- match.arg(process)
  
  B <- nrow(model)
  
  ## SPATIAL GRID ##
  xy <- pred.coords
  colnames(xy) <- c("X", "Y")
  
  if (process == "eta") {
    year <- lubridate::year(date)
    T <- length(unique(year))
    
    firstYear <- min(year) - 1
    
    nt <- c()
    for (tInd in 1:T) 
      nt[tInd] <- nrow(unique(coords[year == firstYear + tInd,]))
    
    ## SPATIAL DISTANCE ##
    allDist <- stats::dist(coords)
    dMax <- max(allDist)
    phi  <- 3 / dMax 
    R.phi         <- list()
    R.phi.inverse <- array(dim = c(max(nt), max(nt), T))
    for (tInd in 1:T) {
      distances <- stats::dist(unique(coords[year == firstYear + tInd,]))
      d <- matrix(0, nrow = nt[tInd], ncol = nt[tInd])
      d[lower.tri(d)] <- distances
      d <- d + t(d)
      R.phi[[tInd]] <- exp(- phi * d)
      R.phi.inverse[1:nt[tInd],1:nt[tInd],tInd] <- solve(R.phi[[tInd]])
    }
    
    for (tInd in 1:T) {
      
      print(paste0("Year ", tInd, " of ", T, " in process..."))
      
      distances <- stats::dist(rbind(xy, unique(coords[year == firstYear + tInd,])))
      d <- matrix(0, nrow = nrow(xy) + nt[tInd], ncol = nrow(xy) + nt[tInd])
      d[lower.tri(d)] <- distances
      d <- d + t(d)
      
      z <- matrix(nrow = 1000, ncol = nrow(xy))
      index <- 1:1000 * as.integer(B / 1000) + (B %% 1000)
      k <- 0
      
      for (nInd in index) {
        k <- k + 1
        Sigma12 <- exp(- phi * d[1:(nrow(xy)), (nrow(xy) + 1):(nrow(xy) + nt[tInd])])
        z[k,]   <- #model[nInd, paste0("etat",tInd)] + 
          Sigma12 %*% R.phi.inverse[1:nt[tInd],1:nt[tInd],tInd] %*% (model[nInd, paste0("etat",tInd,"s",1:nt[tInd])] - model[nInd, paste0("etat",tInd)])
      }
      
      #df <- akima::interp(xy[,1], xy[,2], 
      #  colMeans(z), extrap = TRUE, linear = FALSE,
      #  xo = seq(min(pred.coords[,1]), max(pred.coords[,1]), length = 1000), 
      #  yo = seq(min(pred.coords[,2]), max(pred.coords[,2]), length = 1000))
      
      df <- data.frame("lon" = xy[,1], 
                       "lat" = xy[,2], 
                       "z" = colMeans(z))
      
      #inout <- !is.na(sp::over(
      #  sp::SpatialPoints(df[,1:2], proj4string=sp::CRS(raster::projection(poly))),
      #  methods::as(poly,"SpatialPolygons")
      #))
      
      map.ccb <- ggplot2::ggplot(data = background) + 
        ggplot2::geom_sf(fill= "antiquewhite") + 
        ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) + 
        ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude") + 
        ggplot2::theme(panel.background = ggplot2::element_rect(fill = "aliceblue"),
              axis.text.x=ggplot2::element_text(size=6),
              axis.text.y=ggplot2::element_text(size=6,angle=90),
              axis.title=ggplot2::element_text(size=10,face="bold"))
      
      map.ccb <- map.ccb + 
        ggplot2::geom_tile(data=df, ggplot2::aes(x=df$lon, y=df$lat, fill=z)) + 
        ggplot2::scale_fill_gradient2(midpoint=0, low=scales::muted("blue"), mid="white",
                             high=scales::muted("red"), space ="Lab", limits = zlim) +
        ggplot2::labs(title = bquote(eta[.(tInd+firstYear)](s) - eta[.(tInd+firstYear)]), fill="")
      
      ggplot2::ggsave(paste0("MAPetat",tInd,"s.png"), map.ccb, width = 8.27 / 2, height = 11.69 / 4)
    }
  } else if (process == "alpha") {
    
    indT2 <- which(oneT2 == 1)
    
    nT2 <- nrow(unique(coords[indT2,]))
    
    distT2 <- stats::dist(unique(coords[indT2,]))
    dT2 <- matrix(0, nrow = nT2, ncol = nT2)
    dT2[lower.tri(dT2)] <- distT2
    dT2 <- dT2 + t(dT2)
    
    distances <- stats::dist(rbind(xy, unique(coords[indT2,])))
    d <- matrix(0, nrow = nrow(xy) + nT2, ncol = nrow(xy) + nT2)
    d[lower.tri(d)] <- distances
    d <- d + t(d)
    
    v0 <- matrix(nrow = 1000, ncol = nrow(xy))
    v1 <- matrix(nrow = 1000, ncol = nrow(xy))
    index <- 1:1000 * as.integer(B / 1000) + (B %% 1000)
    k <- 1
    
    for (nInd in index) {
      if (nInd %% 50 == 0) print(paste0("Iteration ", k, " out of 1000..."))
      R.phi.inverse <- solve(exp(- model[nInd,"decaya0"] * dT2))
      Sigma12 <- exp(- model[nInd,"decaya0"] * d[1:(nrow(xy)), (nrow(xy) + 1):(nrow(xy) + nT2)])
      v0[k,]  <- Sigma12 %*% R.phi.inverse %*% model[nInd, paste0("a0s",1:nT2)]
      R.phi.inverse <- solve(exp(- model[nInd,"decaya1"] * dT2))
      Sigma12 <- exp(- model[nInd,"decaya1"] * d[1:(nrow(xy)), (nrow(xy) + 1):(nrow(xy) + nT2)])
      v1[k,]  <- Sigma12 %*% R.phi.inverse %*% model[nInd, paste0("a1s",1:nT2)]
      k <- k + 1
    }
    
    a0 <- model[index, "a11a"] * v0
    a1 <- model[index, "a21a"] * v0 + model[index, "a22a"] * v1
    
    # a0
    df <- data.frame("lon" = xy[,1], 
                     "lat" = xy[,2], 
                     "z" = colMeans(a0))
    
    map.ccb <- ggplot2::ggplot(data = background) + 
      ggplot2::geom_sf(fill= "antiquewhite") + 
      ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) + 
      ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude") + 
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "aliceblue"),
            axis.text.x=ggplot2::element_text(size=6),
            axis.text.y=ggplot2::element_text(size=6,angle=90),
            axis.title=ggplot2::element_text(size=10,face="bold"))
    
    map.ccb <- map.ccb + 
      ggplot2::geom_tile(data=df, ggplot2::aes(x=df$lon, y=df$lat, fill=z)) + 
      ggplot2::scale_fill_gradient2(midpoint=0, low=scales::muted("blue"), mid="white",
                           high=scales::muted("red"), space ="Lab", limits = zlim[1,]) +
      ggplot2::labs(title = bquote(alpha[0](s)), fill="")
    
    ggplot2::ggsave(paste0("MAPa0s.png"), map.ccb, width = 8.27 / 2, height = 11.69 / 4)
    
    # a1
    df <- data.frame("lon" = xy[,1], 
                     "lat" = xy[,2], 
                     "z" = colMeans(a1))
    
    map.ccb <- ggplot2::ggplot(data = background) + 
      ggplot2::geom_sf(fill= "antiquewhite") + 
      ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) + 
      ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude") + 
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "aliceblue"),
            axis.text.x=ggplot2::element_text(size=6),
            axis.text.y=ggplot2::element_text(size=6,angle=90),
            axis.title=ggplot2::element_text(size=10,face="bold"))
    
    map.ccb <- map.ccb + 
      ggplot2::geom_tile(data=df, ggplot2::aes(x=df$lon, y=df$lat, fill=z)) + 
      ggplot2::scale_fill_gradient2(midpoint=0, low=scales::muted("blue"), mid="white",
                           high=scales::muted("red"), space ="Lab", limits = zlim[2,]) +
      ggplot2::labs(title = bquote(alpha[1](s)), fill="")
    
    ggplot2::ggsave(paste0("MAPa1s.png"), map.ccb, width = 8.27 / 2, height = 11.69 / 4)
    
  } else if (process == "lambda") {
    
    indY2 <- which(oneY2 == 1)
    
    nY2 <- nrow(unique(coords[indY2,]))
    
    distY2 <- stats::dist(unique(coords[indY2,]))
    dY2 <- matrix(0, nrow = nY2, ncol = nY2)
    dY2[lower.tri(dY2)] <- distY2
    dY2 <- dY2 + t(dY2)
    
    distances <- stats::dist(rbind(xy, unique(coords[indY2,])))
    d <- matrix(0, nrow = nrow(xy) + nY2, ncol = nrow(xy) + nY2)
    d[lower.tri(d)] <- distances
    d <- d + t(d)
    
    v0 <- matrix(nrow = 1000, ncol = nrow(xy))
    v1 <- matrix(nrow = 1000, ncol = nrow(xy))
    index <- 1:1000 * as.integer(B / 1000) + (B %% 1000)
    k <- 1
    
    for (nInd in index) {
      if (nInd %% 50 == 0) print(paste0("Iteration ", k, " out of 1000..."))
      R.phi.inverse <- solve(exp(- model[nInd,"decayl0"] * dY2))
      Sigma12 <- exp(- model[nInd,"decayl0"] * d[1:(nrow(xy)), (nrow(xy) + 1):(nrow(xy) + nY2)])
      v0[k,]  <- Sigma12 %*% R.phi.inverse %*% model[nInd, paste0("l0s",1:nY2)]
      R.phi.inverse <- solve(exp(- model[nInd,"decayl1"] * dY2))
      Sigma12 <- exp(- model[nInd,"decayl1"] * d[1:(nrow(xy)), (nrow(xy) + 1):(nrow(xy) + nY2)])
      v1[k,]  <- Sigma12 %*% R.phi.inverse %*% model[nInd, paste0("l1s",1:nY2)]
      k <- k + 1
    }
    
    l0 <- model[index, "a11l"] * v0
    l1 <- model[index, "a21l"] * v0 + model[index, "a22l"] * v1
    
    # l0
    df <- data.frame("lon" = xy[,1], 
                     "lat" = xy[,2], 
                     "z" = colMeans(l0))
    
    map.ccb <- ggplot2::ggplot(data = background) + 
      ggplot2::geom_sf(fill= "antiquewhite") + 
      ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) + 
      ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude") + 
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "aliceblue"),
            axis.text.x=ggplot2::element_text(size=6),
            axis.text.y=ggplot2::element_text(size=6,angle=90),
            axis.title=ggplot2::element_text(size=10,face="bold"))
    
    map.ccb <- map.ccb + 
      ggplot2::geom_tile(data=df, ggplot2::aes(x=df$lon, y=df$lat, fill=z)) + 
      ggplot2::scale_fill_gradient2(midpoint=0, low=scales::muted("blue"), mid="white",
                           high=scales::muted("red"), space ="Lab", limits = zlim[1,]) +
      ggplot2::labs(title = bquote(lambda[0](s)), fill="")
    
    ggplot2::ggsave(paste0("MAPl0s.png"), map.ccb, width = 8.27 / 2, height = 11.69 / 4)
    
    # l1
    df <- data.frame("lon" = xy[,1], 
                     "lat" = xy[,2], 
                     "z" = colMeans(l1))
    
    map.ccb <- ggplot2::ggplot(data = background) + 
      ggplot2::geom_sf(fill= "antiquewhite") + 
      ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) + 
      ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude") + 
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "aliceblue"),
            axis.text.x=ggplot2::element_text(size=6),
            axis.text.y=ggplot2::element_text(size=6,angle=90),
            axis.title=ggplot2::element_text(size=10,face="bold"))
    
    map.ccb <- map.ccb + 
      ggplot2::geom_tile(data=df, ggplot2::aes(x=df$lon, y=df$lat, fill=z)) + 
      ggplot2::scale_fill_gradient2(midpoint=0, low=scales::muted("blue"), mid="white",
                           high=scales::muted("red"), space ="Lab", limits = zlim[2,]) +
      ggplot2::labs(title = bquote(lambda[1](s)), fill="")
    
    ggplot2::ggsave(paste0("MAPl1s.png"), map.ccb, width = 8.27 / 2, height = 11.69 / 4)
  } else { # process == "phi"
  
    n <- nrow(unique(coords))
    
    distances <- stats::dist(unique(coords))
    dMax <- max(distances)
    phi  <- 3 / dMax 
    d <- matrix(0, nrow = n, ncol = n)
    d[lower.tri(d)] <- distances
    d <- d + t(d)
    Rinv <- solve(exp(- phi * d))

    distances <- stats::dist(rbind(xy, unique(coords)))
    d <- matrix(0, nrow = nrow(xy) + n, ncol = nrow(xy) + n)
    d[lower.tri(d)] <- distances
    d <- d + t(d)
    
    z <- matrix(nrow = 1000, ncol = nrow(xy))
    index <- 1:1000 * as.integer(B / 1000) + (B %% 1000)
    k <- 0
    
    for (nInd in index) {
      k <- k + 1
      if (k %% 10 == 0) print(paste("Iteration", k, "out of 1000 ..."))
      Sigma12 <- exp(- phi * d[1:(nrow(xy)), (nrow(xy) + 1):(nrow(xy) + n)])
      z[k,]   <- Sigma12 %*% Rinv %*% model[nInd, paste0("phis",1:n)]
    }
    
    df <- data.frame("lon" = xy[,1], 
                     "lat" = xy[,2], 
                     "z" = colMeans(z))
    
    map.ccb <- ggplot2::ggplot(data = background) + 
      ggplot2::geom_sf(fill= "antiquewhite") + 
      ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) + 
      ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude") + 
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "aliceblue"),
            axis.text.x=ggplot2::element_text(size=6),
            axis.text.y=ggplot2::element_text(size=6,angle=90),
            axis.title=ggplot2::element_text(size=10,face="bold"))
    
    map.ccb <- map.ccb + 
      ggplot2::geom_tile(data=df, ggplot2::aes(x=df$lon, y=df$lat, fill=z)) + 
      ggplot2::scale_fill_gradient2(midpoint=0, low=scales::muted("blue"), mid="white",
                           high=scales::muted("red"), space ="Lab", limits = zlim) +
      ggplot2::labs(title = bquote(gamma[0](s)), fill="")
  
    ggplot2::ggsave(paste0("MAPphis.png"), map.ccb, width = 8.27 / 2, height = 11.69 / 4)
  }
}
