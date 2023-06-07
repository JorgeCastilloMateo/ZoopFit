## This is a working script to reproduce some of the results
## It's not clean nor updated or organized, 
## Some results can be the result of changing the code below slightly
## But this can be useful for reproducibility

## LOAD PACKAGE AND DATA ##
library("ZoopFit")
data("zoop_snippet.rda")

# FIT TWO CHAINS OF THE MODEL
set.seed(12345)
MODEL <- list()
MODEL[[1]] <- GibbsZoop(
  Y1, Y2, T1, T2, 
  X[,1:3], date, coords, 
  oneY1, oneY2, oneT1, oneT2,
  calibration = "corTime",
  n.burnin = 250000, max.iter = 250000, n.thin = 250, n.report = 1000
)
MODEL[[2]] <- GibbsZoop(
  Y1, Y2, T1, T2, 
  X[,1:3], date, coords, 
  oneY1, oneY2, oneT1, oneT2,
  calibration = "corTime",
  n.burnin = 250000, max.iter = 250000, n.thin = 250, n.report = 1000
)

# CHECK POTENTIAL SCALE REDUCTION FACTOR
library(coda)
model <- MODEL
R <- rep(NA, ncol(model[[1]]))
for (i in 1:ncol(model[[1]])) {
  print(i)
  R[i] <- try(gelman.diag(mcmc.list(mcmc(model[[1]][,i]), mcmc(model[[2]][,i])))$psrf[1])
}
max(R[!is.nan(R)])

saveRDS(MODEL, file = "MODEL.rds")

MODEL <- rbind(MODEL[[1]], MODEL[[2]])

# SIMULATE GP's (fixed decay)
# @param GP values of the GP at the observed locations
# @param mean,prec posterior mean and precision of the GP
# @param coords coordinates where GP was fitted
# @param grid new grid of locations to krige the GP
GP <- function(GP, mean, prec, coords, grid) {
  
  ### PREDICTION ###
  n <- nrow(coords)
  B <- length(prec)
  
  # Sigma
  distances       <- dist(coords)
  d               <- matrix(0, nrow = n, ncol = n)
  d[lower.tri(d)] <- distances
  d               <- d + t(d)
  
  dmax            <- max(distances)
  phi             <- 3 / dmax
  
  R.phi           <- exp(-phi * d)
  R.phi.inverse   <- solve(R.phi)
  
  # Sigma, Sigma00, Sigmai0
  xy <- grid
  colnames(xy) <- c("X", "Y")
  
  coords <- rbind(xy, coords)  
  distances <- dist(coords)
  d <- matrix(0, nrow = nrow(coords), ncol = nrow(coords))
  d[lower.tri(d)] <- distances
  d <- d + t(d)
  
  # Sample field
  m <- nrow(xy)
  SIGMA    <- exp(- phi * d)
  Sigma00  <- SIGMA[1:m, 1:m]
  Sigmai0T <- SIGMA[1:m, (m + 1):nrow(coords)]
  
  SIGMA_pre  <- Sigmai0T %*% R.phi.inverse
  SIGMA_cond <- Sigma00 - SIGMA_pre %*% t(Sigmai0T)
  
  MVN <- ZoopFit:::mvrnormArma(B, Sigma = SIGMA_cond)
  
  campos <- matrix(MVN, nrow = B, ncol = m)
  
  for (n.sim in 1:B) {
    print(n.sim)
    media <- mean[n.sim] + SIGMA_pre %*% matrix(GP[n.sim,] - mean[n.sim], ncol = 1)
    campos[n.sim,] <- media + 1 / sqrt(prec[n.sim]) * campos[n.sim,]
  }
  
  return(campos)
}

# BAYESIAN KRIGING FOR THE GP's
year <- lubridate::year(date)
T <- length(unique(year))
firstYear <- min(year) - 1
nt <- c()
for (tInd in 1:T) 
  nt[tInd] <- nrow(unique(coords[year == firstYear + tInd, ]))

set.seed(12345)
gamma0s <- GP(MODEL[,paste0("phis",1:435)], rep(0, 2000), MODEL[,"precphi"], unique(coords), grid)
etats <- list()
for (t in 1:17) {
  print(t)
  etats[[t]] <- GP(MODEL[,paste0("etat", t, "s", 1:nt[t])], MODEL[,paste0("etat",t)], MODEL[,"precEtas"], unique(coords[year == firstYear + t,]), grid)
}

# 4/29/2011 (regime 3)
yday("2011-04-29")
# 3/21/2017 (regime 2)
yday("2017-03-21")
# 4/17/2019 (regime 3)
yday("2019-04-17")

# MAPS OF POSTERIOR MEAN OF ZOOP and OF
# PROBABILITY OF EXCEEDING 1000 AND 1500
# @param model matrix of posterior values of the model parameters
# @param date date to predict
# @param predictWhales same as in predictZoopMap()
# @param gamma0s,etats kriged GP's
# @param grid,background grid of points within CCB and the background of eeuu
#   the packages rnaturalearthdata and rnaturalearth
# @param xlim,ylim,zlim coodinates and zoop values limits
predictMAPZOOP <- function(model, date, predictWhales, gamma0s, etats, 
                           grid, background = eeuu, 
                           xlim = c(-61771.464, 1228.536), ylim = c(-34392.144, 13607.86), zlim = NULL) {
  
  j <- lubridate::yday(date)
  t <- lubridate::year(date) - 2002
  
  ZOOP <- predictZoopMap(model, date, predictWhales, gamma0s, etats[[t]])
  
  # plot mean
  df <- data.frame("lon" = grid[,1], 
                   "lat" = grid[,2], 
                   "z" = colMeans(ZOOP))
  
  map.ccb <- ggplot2::ggplot(data = background) + 
    ggplot2::geom_sf(fill= "antiquewhite") + 
    ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) + 
    ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude") + 
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "aliceblue"),
                   axis.text.x=ggplot2::element_text(size=6),
                   axis.text.y=ggplot2::element_text(size=6,angle=90),
                   axis.title=ggplot2::element_text(size=10,face="bold"))
  
  map.ccb <- map.ccb + 
    ggplot2::geom_tile(data=df, ggplot2::aes(x=lon, y=lat, fill=z)) + 
    ggplot2::scale_fill_gradient2(midpoint=mean(df$z), low=scales::muted("blue"), mid="white",
                                  high=scales::muted("red"), space ="Lab", limits = zlim[1,]) +
    ggplot2::labs(title = bquote("Mean " * ZOOP[.(t+2002) * "," * .(j)]^{"true"}*(s)), subtitle = date, fill="")
  
  ggplot2::ggsave(paste0("MAPmeanZOOP",date,".png"), map.ccb, width = 8.27 / 2, height = 11.69 / 4)
  
  print(mean(df$z))
  
  # plot sd
  df <- data.frame("lon" = grid[,1], 
                   "lat" = grid[,2], 
                   "z" = apply(ZOOP, 2, sd))
  
  map.ccb <- ggplot2::ggplot(data = background) + 
    ggplot2::geom_sf(fill= "antiquewhite") + 
    ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) + 
    ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude") + 
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "aliceblue"),
                   axis.text.x=ggplot2::element_text(size=6),
                   axis.text.y=ggplot2::element_text(size=6,angle=90),
                   axis.title=ggplot2::element_text(size=10,face="bold"))
  
  map.ccb <- map.ccb + 
    ggplot2::geom_tile(data=df, ggplot2::aes(x=lon, y=lat, fill=z)) + 
    ggplot2::scale_fill_gradient2(midpoint=mean(df$z), low=scales::muted("blue"), mid="white",
                                  high=scales::muted("red"), space ="Lab", limits = zlim[1,]) +
    ggplot2::labs(title = bquote("Sd " * ZOOP[.(t+2002) * "," * .(j)]^{"true"}*(s)), subtitle = date, fill="")
  
  ggplot2::ggsave(paste0("MAPsdZOOP",date,".png"), map.ccb, width = 8.27 / 2, height = 11.69 / 4)
  
  # plot probability ZOOP > 1000
  df <- data.frame("lon" = grid[,1], 
                   "lat" = grid[,2], 
                   "z" = colMeans(ZOOP > 1000))
  
  map.ccb <- ggplot2::ggplot(data = background) + 
    ggplot2::geom_sf(fill= "antiquewhite") + 
    ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) + 
    ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude") + 
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "aliceblue"),
                   axis.text.x=ggplot2::element_text(size=6),
                   axis.text.y=ggplot2::element_text(size=6,angle=90),
                   axis.title=ggplot2::element_text(size=10,face="bold"))
  
  map.ccb <- map.ccb + 
    ggplot2::geom_tile(data=df, ggplot2::aes(x=lon, y=lat, fill=z)) + 
    ggplot2::scale_fill_gradient2(midpoint=mean(df$z), low=scales::muted("blue"), mid="white",
                                  high=scales::muted("red"), space ="Lab", limits = zlim[2,]) +
    ggplot2::labs(title = bquote(P(ZOOP[.(t+2002) * "," * .(j)]^{"true"}*(s) > 1000*" | "*"data")), subtitle = date, fill="")
  
  ggplot2::ggsave(paste0("MAPprobZOOP1000",date,".png"), map.ccb, width = 8.27 / 2, height = 11.69 / 4)
  
  # plot probability ZOOP > 1500
  df <- data.frame("lon" = grid[,1], 
                   "lat" = grid[,2], 
                   "z" = colMeans(ZOOP > 1500))
  
  map.ccb <- ggplot2::ggplot(data = background) + 
    ggplot2::geom_sf(fill= "antiquewhite") + 
    ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) + 
    ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude") + 
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "aliceblue"),
                   axis.text.x=ggplot2::element_text(size=6),
                   axis.text.y=ggplot2::element_text(size=6,angle=90),
                   axis.title=ggplot2::element_text(size=10,face="bold"))
  
  map.ccb <- map.ccb + 
    ggplot2::geom_tile(data=df, ggplot2::aes(x=lon, y=lat, fill=z)) + 
    ggplot2::scale_fill_gradient2(midpoint=mean(df$z), low=scales::muted("blue"), mid="white",
                                  high=scales::muted("red"), space ="Lab", limits = zlim[2,]) +
    ggplot2::labs(title = bquote(P(ZOOP[.(t+2002) * "," * .(j)]^{"true"}*(s) > 1500*" | "*"data")), subtitle = date, fill="")
  
  ggplot2::ggsave(paste0("MAPprobZOOP1500",date,".png"), map.ccb, width = 8.27 / 2, height = 11.69 / 4)
}

# 4/29/2011 (regime 3)
predictMAPZOOP(MODEL, "2011-04-29", predictWhales, 
               gamma0s, etats, grid, eeuu, zlim = rbind(c(250, 2650), c(0,1)))
# 3/21/2017 (regime 2)
predictMAPZOOP(MODEL, "2017-03-21", predictWhales, 
               gamma0s, etats, grid, eeuu, zlim = rbind(c(150, 700), c(0,.2)))
# 4/17/2019 (regime 3)
predictMAPZOOP(MODEL, "2019-04-17", predictWhales, 
               gamma0s, etats, grid, eeuu, zlim = rbind(c(150, 700), c(0,.2)))


# MAPS OF AVERAGE ZOOP BY YEAR
predictMAPZOOPyear <- function(model, year, dday = 1:181, predictWhales, gamma0s, etats, 
                               grid, background = eeuu, 
                               xlim = c(-61771.464, 1228.536), ylim = c(-34392.144, 13607.86), zlim = NULL) {
  
  t <- year - 2002
  ZOOP <- matrix(0, nrow = nrow(gamma0s), ncol = ncol(gamma0s))
  df <- data.frame(dayyear = dday, year = year)
  df$origin <- as.Date(paste0(df$year, "-01-01"),tz = "UTC") - days(1)
  my.date <- as.Date(df$dayyear, origin = df$origin, tz = "UTC")
  for (j in (dday - dday[1] + 1)) {
    ZOOP <- ZOOP + predictZoopMap(model, my.date[j], predictWhales, gamma0s, etats[[t]])
  }
  ZOOP <- ZOOP / length(dday)
  # plot mean
  df <- data.frame("lon" = grid[,1], 
                   "lat" = grid[,2], 
                   "z" = colMeans(ZOOP))
  
  map.ccb <- ggplot2::ggplot(data = background) + 
    ggplot2::geom_sf(fill= "antiquewhite") + 
    ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) + 
    ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude") + 
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "aliceblue"),
                   axis.text.x=ggplot2::element_text(size=6),
                   axis.text.y=ggplot2::element_text(size=6,angle=90),
                   axis.title=ggplot2::element_text(size=10,face="bold"))
  
  map.ccb <- map.ccb + 
    ggplot2::geom_tile(data=df, ggplot2::aes(x=lon, y=lat, fill=z)) + 
    ggplot2::scale_fill_gradient2(midpoint=mean(df$z), low=scales::muted("blue"), mid="white",
                                  high=scales::muted("red"), space ="Lab", limits = zlim[1,]) +
    ggplot2::labs(title = bquote("Mean " * ZOOP[.(t+2002) * "," * .(min(dday)) * ":" * .(max(dday))]^{"true"}*(s)), fill="")
  
  ggplot2::ggsave(paste0("MAPmeanZOOPyear",year,"days",paste(range(dday), collapse = "-"),".png"), map.ccb, width = 8.27 / 2, height = 11.69 / 4)
  
  print(mean(df$z))
}

for (i in 2003:2019) {
  print(i)
  predictMAPZOOPyear(MODEL, i, 1:181, predictWhales, 
                     gamma0s, etats, grid, eeuu)
}

predictMAPZOOPyear(MODEL, 2017, 43:56, predictWhales, 
                   gamma0s, etats, grid, eeuu)
predictMAPZOOPyear(MODEL, 2017, 71:84, predictWhales, 
                   gamma0s, etats, grid, eeuu)
predictMAPZOOPyear(MODEL, 2017, 127:140, predictWhales, 
                   gamma0s, etats, grid, eeuu)

# DAILY TIME SERIES OF ZOOP
tsZOOP <- function() {
  
  B <- nrow(gamma0s)
  m <- ncol(gamma0s)
  
  predictTemp <- array(dim = c(B, 181, 17))
  predictZOOP <- array(dim = c(B, 181, 17))
  IndZOOPabove1000 <- array(dim = c(B, 181, 17))
  for (j in 1:181) {
    print(j)
    for (t in 1:17) {
      
      pT <- 
        MODEL[,"sine"]   * sin(2 * pi * j / 365) + 
        MODEL[,"cosine"] * cos(2 * pi * j / 365) +
        MODEL[,paste0("psit",t)] +
        gamma0s
      
      predictTemp[,j,t] <- rowMeans(pT)
      
      pZ0 <- 
        MODEL[,"betaTemp"] * pT +
        MODEL[,"beta2"] * ifelse(32 < j & j < 90, 1, 0) + MODEL[,"beta3"] * ifelse(90 <= j, 1, 0) +
        etats[[t]]
      
      pZ1 <- pZ0 + MODEL[,"beta1"]
      
      p <- predictWhales[,j + 181 * (t - 1)]
      
      aux <- p * exp(pZ1) + (1 - p) * exp(pZ0)
      predictZOOP[,j,t] <- rowMeans(aux)
      IndZOOPabove1000[,j,t] <- rowMeans(aux > 1000)
    }
  }
  
  return(predictZOOP)
}

Ynew <- probWhales
Ynew <- predictTemp
Ynew <- predictZOOP
Ynew <- IndZOOPabove1000

# A4 landscape HALF HIGH
Ynew <- predictZOOP
m <- 181*17
YnewM   <- c(apply(Ynew, 2:3, mean))
Ynewq05 <- c(apply(Ynew, 2:3, quantile, prob = 0.05))
Ynewq95 <- c(apply(Ynew, 2:3, quantile, prob = 0.95))
plot(x = 1:m, YnewM, type = "n", #ylim = c(0, 10000),
     ylab = expression(paste("Averaged ZOOP abundance / ", m^3)), 
     xlab = "Day within year", xaxt='n', cex.main=1.25, cex.lab=1.25) 

for (i in 1:17) {
  lines(x = 1:181 + (i-1) * 181, y = YnewM[1:181 + (i-1) * 181])
}
axis(side = 1, at = c(1, 1:17 * 181), labels = c(1, rep(181, 17)), cex=1.25)
mtext("Year", side=3, line=1, cex=1.25)
title(main=expression("Mean of "*widehat(ZOOP)["t,j"]^{"true"}*"(CCB)"), line=3, cex=1.5)
mtext(2003:2019, side = 3, line = 0, at = 91 + 181 * 0:16, cex=1.25)
abline(v = 0:17 * 181, col = "black", lty = 2)
abline(v = 0:17 * 181 + 32, col = "red", lty = 1)
abline(v = 0:17 * 181 + 90, col = "red", lty = 1)
abline(v = 0:17 * 181 + 140, col = "blue", lty = 1)
par(mar = c(5,5,5,2) + 0.1) 

# A4 landscape HALF HIGH
Ynew <- IndZOOPabove1000
m <- 181*17
YnewM   <- c(apply(Ynew, 2:3, mean))
plot(x = 1:m, YnewM, type = "n", ylim = c(0, 1),
     ylab = "Extent", 
     xlab = "Day within year", xaxt='n', cex.main=1.25, cex.lab=1.25) 
for (i in 1:17) {
  lines(x = 1:181 + (i-1) * 181, y = YnewM[1:181 + (i-1) * 181])
}
axis(side = 1, at = c(1, 1:17 * 181), labels = c(1, rep(181, 17)), cex=1.25)
mtext("Year", side=3, line=1, cex=1.25)
title(main=expression("Mean of the extent of 1("*ZOOP["t,j"]^{"true"}*"(s) > 1000)"), line=3, cex=1.5)
mtext(2003:2019, side = 3, line = 0, at = 91 + 181 * 0:16, cex=1.25)
abline(v = 0:17 * 181, col = "black", lty = 2)
par(mar = c(5,5,5,2) + 0.1) ## default is c(5,4,4,2) + 0.1


pdf("blockaverageZOOPyear.pdf", 5, 5)
m <- 17
Ynew <- predictZOOP
par(mar = c(5,4.5,4,2) + 0.1)
plot(x = 1:m, c(apply(Ynew, 3, mean)), type = "l",
     main = "Annual block average",
     ylim = c(200, 2000),
     ylab = expression(paste("Averaged ZOOP abundance / ", m^3)), xlab = "Year", xaxt = "n") 
axis(1, at=1:m, labels=2003:2019, las = 2)

YnewM <- apply(Ynew[,91:181,], 3, mean)
lines(x = 1:m, YnewM, lty = 2)
legend("topleft", 
       legend=rev(c(expression(widehat(ZOOP)["t."]^{"true"}*"(CCB)"),
                    expression(widehat(ZOOP)[t*","*r[j]*"="*3]^{"true"}*"(CCB)"))),
       lty=rev(1:2))
dev.off()

pdf("blockaverageZOOPyearProb.pdf", 5, 5)
m <- 17
Ynew <- IndZOOPabove1000
par(mar = c(5,4.5,4,2) + 0.1)
plot(x = 1:m, c(apply(Ynew, 3, mean)), type = "l",
     main = "Annual extent",
     ylim = c(0, 1),
     ylab = "Extent", xlab = "Year", xaxt = "n") 
axis(1, at=1:m, labels=2003:2019, las = 2)

YnewM <- apply(Ynew[,91:181,], 3, mean)
lines(x = 1:m, YnewM, lty = 2)
dev.off()

pdf("blockaverageZOOPmonth.pdf", 5, 5)
m <- 6 * 17
YnewM <- matrix(nrow = 6, ncol = 17)
for (i in 1:17) {
  YnewM[1,i] <- mean(Ynew[,1:31,i])
  YnewM[2,i] <- mean(Ynew[,32:59,i])
  YnewM[3,i] <- mean(Ynew[,60:90,i])
  YnewM[4,i] <- mean(Ynew[,91:120,i])
  YnewM[5,i] <- mean(Ynew[,121:151,i])
  YnewM[6,i] <- mean(Ynew[,152:181,i])
}
plot(x = 1:m, c(YnewM), type = "l",
     main = expression("Monthly average "*widehat(ZOOP)["t,j"]^{"true"}*"(CCB)"),
     ylim = c(200, 2250),
     ylab = expression(paste("Averaged ZOOP abundance / ", m^3)), xlab = "Year", xaxt = "n") 
axis(1, at=1:m, labels=2003:2019, las = 2)
polygon(c(1:m, m:1), c(apply(apply(Ynew,c(1,3),mean), 2, quantile, prob = 0.05), rev(apply(apply(Ynew,c(1,3),mean), 2, quantile, prob = 0.95))), col = "gray", border = "gray")
lines(x = 1:m, apply(Ynew, 3, mean))
dev.off()

pdf("blockaverageZOOPregime3.pdf", 5, 5)
m <- 17
YnewM <- apply(Ynew[,91:181,], 3, mean)
plot(x = 1:m, c(YnewM), type = "l",
     main = expression(widehat(ZOOP)[t*","*r[j]*"="*3]^{"true"}*"(CCB)"),
     ylim = c(200, 2250),
     ylab = expression(paste("Averaged ZOOP abundance / ", m^3)), xlab = "Year", xaxt = "n") 
axis(1, at=1:m, labels=2003:2019, las = 2)
polygon(c(1:m, m:1), c(apply(apply(Ynew[,91:181,],c(1,3),mean), 2, quantile, prob = 0.05), rev(apply(apply(Ynew[,91:181,],c(1,3),mean), 2, quantile, prob = 0.95))), col = "gray", border = "gray")
lines(x = 1:m, c(YnewM))
dev.off()

