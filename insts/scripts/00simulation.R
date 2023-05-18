# 

# GIBBS SAMPLER
model1 <- function(Y, X, coords, newX, newcoords, 
                   n.burnin = 10, max.iter = 10, n.thin = 1) {
  
  n    <- nrow(coords)
  nnew <- nrow(newcoords)
  
  d <- dist(rbind(newcoords, coords))
  D <- matrix(0, nrow = nnew + n, ncol = nnew + n)
  D[lower.tri(D)] <- d
  D <- D + t(D)
  R <- exp(- 3 * 0.7 / max(d) * D)
  Rinv <- solve(R[nnew + 1:n, nnew + 1:n])
  R11 <- R[1:nnew, 1:nnew]
  R12 <- R[1:nnew, nnew + 1:n]
  R21 <- t(R12)
  R12Rinv <- R12 %*% Rinv
  Rnew <- R11 - R12Rinv %*% R21
  
  keep   <- matrix(nrow = max.iter / n.thin, 
                   ncol = 4 + n + nnew)
  params <- rep(0, 4 + n + nnew); params[3:4] <- 1
  
  I2 <- diag(2)
  In <- diag(n)
  
  for (b in (-n.burnin + 1):max.iter) {
    if (b %% 100 == 0)
      print(paste("Iteration number:", b))
    
    # beta
    Xt_Rinv <- t(X) %*% Rinv
    delta <- solve(Xt_Rinv %*% X * params[3] + I2 * 1e-04)
    chi   <- Xt_Rinv %*% params[4 + 1:n] * params[3]
    params[1:2] <- rmvnorm(1, delta %*% chi, delta, 
                           checkSymmetry = FALSE)
    Xb <- X %*% params[1:2]
    
    # prec true
    vect <- params[4 + 1:n] - Xb
    params[3] <- rgamma(1, n / 2 + 2, t(vect) %*% Rinv %*% vect / 2 + 1)
    
    # prec 1
    vect <- Y - params[4 + 1:n]
    params[4] <- rgamma(1, n / 2 + 2, t(vect) %*% vect / 2 + 1)
    
    # Y true
    delta <- solve(In * params[4] + Rinv * params[3])
    chi   <- Y * params[4] + Rinv %*% Xb * params[3]
    params[4 + 1:n] <- rmvnorm(1, delta %*% chi, delta,
                               checkSymmetry = FALSE)
    
    # Y true interpolation
    mu <- newX %*% params[1:2] + R12Rinv %*% (params[4 + 1:n] - Xb)
    params[4 + n + 1:nnew] <- rmvnorm(1, mu, Rnew / params[3],
                               checkSymmetry = FALSE)
    
    if ((b > 0) && (b %% n.thin == 0))
      keep[b / n.thin, ] <- params
  }
  
  colnames(keep) <- c("beta0", "beta1", "precTrue", "prec1",
                      paste0("Ytrue", 1:n), 
                      paste0("YtruePred", 1:nnew))
  
  return(keep)
}

model2 <- function(Y1, Y2, ind, X, coords, newX, newcoords, 
                   n.burnin = 10, max.iter = 10, n.thin = 1) {
  
  ind1 <- which(ind == 1)
  ind2 <- which(ind == 2)
  
  n    <- nrow(coords)
  n1   <- length(Y1)
  n2   <- length(Y2)
  nnew <- nrow(newcoords)
  
  d <- dist(rbind(newcoords, coords))
  D <- matrix(0, nrow = nnew + n, ncol = nnew + n)
  D[lower.tri(D)] <- d
  D <- D + t(D)
  R <- exp(- 3 * 0.7 / max(d) * D)
  Rinv <- solve(R[nnew + 1:n, nnew + 1:n])
  R11 <- R[1:nnew, 1:nnew]
  R12 <- R[1:nnew, nnew + 1:n]
  R21 <- t(R12)
  R12Rinv <- R12 %*% Rinv
  Rnew <- R11 - R12Rinv %*% R21
  
  keep   <- matrix(nrow = max.iter / n.thin, 
                   ncol = 7 + n + nnew)
  params <- rep(0, 7 + n + nnew); params[c(3:4, 7)] <- 1
  
  I2 <- diag(2)
  In <- diag(n)
  M  <- diag(n)
  Z  <- rep(NA, n)
  Z[ind1] <- Y1
  
  for (b in (-n.burnin + 1):max.iter) {
    if (b %% 100 == 0)
      print(paste("Iteration number:", b))
    
    # beta
    #Xt_Rinv <- t(X) %*% Rinv
    #delta <- solve(Xt_Rinv %*% X * params[3] + I2 * 1e-04)
    #chi   <- Xt_Rinv %*% params[7 + 1:n] * params[3]
    #params[1:2] <- rmvnorm(1, delta %*% chi, delta, 
    #                       checkSymmetry = FALSE)
    params[1:2] <- c(10, -5)
    Xb <- X %*% params[1:2]
    
    # prec true
    vect <- params[7 + 1:n] - Xb
    params[3] <- rgamma(1, n / 2 + 2, t(vect) %*% Rinv %*% vect / 2 + 1)
    #params[3] <- 1 / 16
    
    # prec 1
    vect <- Y1 - params[7 + ind1]
    params[4] <- rgamma(1, n1 / 2 + 2, t(vect) %*% vect / 2 + 1)
    #params[4] <- 1 / 2
    
    # lambda
    XX    <- cbind(1, params[7 + ind2])
    #delta <- solve(t(XX) %*% XX * params[7] + I2 * 1e-04)
    #chi   <- t(XX) %*% Y2 * params[7] + c(0, 1) * 1e-04
    #params[5:6] <- rmvnorm(1, delta %*% chi, delta, 
    #                       checkSymmetry = FALSE)
    params[5:6] <- c(1, .5)
    
    # prec 2
    vect <- Y2 - XX %*% params[5:6]
    params[7] <- rgamma(1, n2 / 2 + 2, t(vect) %*% vect / 2 + 1)
    #params[7] <- 1 / 3
    
    # Y true
    Z[ind2] <- (Y2 - params[5]) / params[6]
    diag(M)[ind1] <- params[4]
    diag(M)[ind2] <- params[7] * params[6]^2
    delta <- solve(M + Rinv * params[3])
    chi   <- M %*% Z + Rinv %*% Xb * params[3]
    params[7 + 1:n] <- rmvnorm(1, delta %*% chi, delta,
                               checkSymmetry = FALSE)
    
    # Y true interpolation
    mu <- newX %*% params[1:2] + R12Rinv %*% (params[7 + 1:n] - Xb)
    params[7 + n + 1:nnew] <- rmvnorm(1, mu, Rnew / params[3],
                                      checkSymmetry = FALSE)
    
    if ((b > 0) && (b %% n.thin == 0))
      keep[b / n.thin, ] <- params
  }
  
  colnames(keep) <- c("beta0", "beta1", "precTrue", "prec1",
                      "lambda0", "lambda1", "prec2",
                      paste0("Ytrue", 1:n), 
                      paste0("YtruePred", 1:nnew))
  
  return(keep)
}

# GRID
set.seed(23)
grid <- expand.grid(x = 1:10, y = 1:10)
ind  <- 1:100
out  <- sort(sample(ind, size = 20))
ind1 <- sort(sample(ind[-out], size = 40))
ind2 <- sort(ind[-c(out, ind1)])

d <- dist(grid)
D <- matrix(0, nrow = 100, ncol = 100)
D[lower.tri(D)] <- d
D <- D + t(D)
R <- exp(- 3 * 0.7 / max(d) * D)

# PARAMETERS and DATA
library("sf")
library("coda")
library("ZoopFit")
library("ggplot2")
library("mvtnorm")
set.seed(23)
X     <- cbind(1, rnorm(100, 5, 5))
beta  <- c(10, -5)
Ytrue <- c(rmvnorm(1, X %*% beta, 16 * R))
Y1    <- rnorm(100, Ytrue, 2)
Y2    <- rnorm(100, 1 + .5 * Ytrue, 3)
my_sf <- st_as_sf(data.frame(grid, Ytrue), coords = c("x", "y"))
ggplot(my_sf) + 
  geom_tile(data = data.frame(grid, Z = Ytrue - X %*% beta), aes(x = x, y = y, fill = Z))

# FIT MODELS
m1 <- list()
m1[[1]] <- model1(Y1[ind1], X[ind1,], grid[ind1,],
                  X[out,], grid[out,],
                  10000, 10000, 10)
m1[[2]] <- model1(Y1[ind1], X[ind1,], grid[ind1,],
                  X[out,], grid[out,],
                  10000, 10000, 10)

m2 <- list()
index <- rep(NA, 80)
index[match(ind1, sort(c(ind1, ind2)))] <- 1
index[match(ind2, sort(c(ind1, ind2)))] <- 2
m2[[1]] <- model2(Y1[ind1], Y2[ind2], index, X[-out,], 
                  grid[-out,], X[out,], grid[out,], 
                  n.burnin = 10000, max.iter = 10000, n.thin = 10) 
m2[[2]] <- model2(Y1[ind1], Y2[ind2], index, X[-out,], 
                  grid[-out,], X[out,], grid[out,], 
                  n.burnin = 10000, max.iter = 10000, n.thin = 10) 

colMeans(m2[[1]])[1:8]
colMeans(m2[[2]])[1:8]

min(effectiveSize(m1[[1]]) + effectiveSize(m1[[2]]))
gelman.diag(mcmc.list(mcmc(m1[[1]]), mcmc(m1[[2]])))

min(effectiveSize(m2[[1]]) + effectiveSize(m2[[2]]))
gelman.diag(mcmc.list(mcmc(m2[[1]]), mcmc(m2[[2]])))

ZoopFit::validation(Ytrue[out], m1[[1]][,paste0("YtruePred", 1:20)])
ZoopFit::validation(Ytrue[out], m1[[2]][,paste0("YtruePred", 1:20)])
ZoopFit::validation(Ytrue[out], m2[[1]][,paste0("YtruePred", 1:20)])
ZoopFit::validation(Ytrue[out], m2[[2]][,paste0("YtruePred", 1:20)])

plot(Ytrue[out], ylim = c(-70, 50))
lines(colMeans(m1[[1]][,paste0("YtruePred", 1:20)]), col="red")
lines(apply(m1[[1]][,paste0("YtruePred", 1:20)], 2, quantile, prob = 0.05), col="red")
lines(apply(m1[[1]][,paste0("YtruePred", 1:20)], 2, quantile, prob = 0.95), col="red")

plot(Ytrue[out], ylim = c(-70, 50))
lines(colMeans(m2[[2]][,paste0("YtruePred", 1:20)]))
lines(apply(m2[[2]][,paste0("YtruePred", 1:20)], 2, quantile, prob = 0.05))
lines(apply(m2[[2]][,paste0("YtruePred", 1:20)], 2, quantile, prob = 0.95))

plot(Ytrue[-out], ylim = c(-70, 50))
lines(colMeans(m2[[1]][,paste0("Ytrue", 1:80)]))
lines(apply(m2[[1]][,paste0("Ytrue", 1:80)], 2, quantile, prob = 0.05))
lines(apply(m2[[1]][,paste0("Ytrue", 1:80)], 2, quantile, prob = 0.95))

plot(Ytrue[ind1], ylim = c(-70, 50))
lines(Y1[ind1], col = "red")
lines(colMeans(m2[[1]][,paste0("Ytrue", which(index==1))]))
lines(apply(m2[[1]][,paste0("Ytrue", which(index==1))], 2, quantile, prob = 0.05))
lines(apply(m2[[1]][,paste0("Ytrue", which(index==1))], 2, quantile, prob = 0.95))

