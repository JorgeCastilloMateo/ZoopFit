### PAPER EES -- SECTION 3 ###
###   SIMULATION EXAMPLE   ###

# PACKAGES
library("coda")
library("ZoopFit")
library("mvtnorm")

# GIBBS SAMPLER model 1: Spatial LM
# Y = Ytrue + e1, 
# Ytrue = X %*% beta + W
model1 <- function(Y, X, coords, newX, newcoords, 
                   inits,
                   n.burnin = 1000, max.iter = 1000, n.thin = 1) {
  
  k    <- ncol(X)
  n    <- nrow(coords)
  nnew <- nrow(newcoords)
  
  d <- dist(rbind(newcoords, coords))
  D <- matrix(0, nrow = nnew + n, ncol = nnew + n)
  D[lower.tri(D)] <- d
  D <- D + t(D)
  R <- exp(- 3 * 0.6 / max(d) * D)
  Rinv <- solve(R[nnew + 1:n, nnew + 1:n])
  R11 <- R[1:nnew, 1:nnew]
  R12 <- R[1:nnew, nnew + 1:n]
  R21 <- t(R12)
  R12Rinv <- R12 %*% Rinv
  Rnew <- R11 - R12Rinv %*% R21
  
  XtRinv  <- t(X) %*% Rinv
  XtRinvX <- XtRinv %*% X
  
  keep   <- matrix(nrow = max.iter / n.thin, 
                   ncol = k + 2 + n + nnew)
  
  if (missing(inits) | any(is.na(inits))) {
    params <- rep(0, k + 2 + n + nnew)
    params[k + 1:2] <- 1
  } else {
    params <- inits
  }
  
  Ik <- diag(k)
  In <- diag(n)
  
  for (b in (-n.burnin + 1):max.iter) {
    if (b %% 1000 == 0)
      print(paste("Iteration number:", b))
    
    # beta
    delta <- solve(XtRinvX * params[k + 1] + Ik * 1e-04)
    chi   <- XtRinv %*% params[k + 2 + 1:n] * params[k + 1]
    params[1:k] <- mvtnorm::rmvnorm(1, delta %*% chi, delta, 
      checkSymmetry = FALSE)
    Xb <- X %*% params[1:k]
    
    # prec true
    vect <- params[k + 2 + 1:n] - Xb
    params[k + 1] <- rgamma(1, n / 2 + 2, t(vect) %*% Rinv %*% vect / 2 + 1)
    
    # prec 1
    vect <- Y - params[k + 2 + 1:n]
    params[k + 2] <- rgamma(1, n / 2 + 2, t(vect) %*% vect / 2 + 1)
    
    # Y true
    delta <- solve(In * params[k + 2] + Rinv * params[k + 1])
    chi   <- Y * params[k + 2] + Rinv %*% Xb * params[k + 1]
    params[k + 2 + 1:n] <- mvtnorm::rmvnorm(1, delta %*% chi, delta,
      checkSymmetry = FALSE)
    
    if ((b > 0) && (b %% n.thin == 0)) {
      # Y true interpolation
      mu <- newX %*% params[1:k] + R12Rinv %*% (params[k + 2 + 1:n] - Xb)
      params[k + 2 + n + 1:nnew] <- mvtnorm::rmvnorm(1, mu, Rnew / params[k + 1],
        checkSymmetry = FALSE)
      
      keep[b / n.thin, ] <- params
    }
  }
  
  colnames(keep) <- c(paste0("beta", 1:k - 1), 
                      "precTrue", "prec1",
                      paste0("Ytrue", 1:n), 
                      paste0("YtruePred", 1:nnew))
  
  return(keep)
}

# GIBBS SAMPLER model 2: Calibration model
# Y1 = Ytrue + e1, 
# Y2 = lambda0 + lambda1 * Ytrue + e2,
# Ytrue = X %*% beta + W
model2 <- function(Y1, Y2, ind, X, coords, newX, newcoords,
                   inits,
                   n.burnin = 1000, max.iter = 1000, n.thin = 1) {
  
  ind1 <- which(ind == 1)
  ind2 <- which(ind == 2)
  
  k    <- ncol(X)
  n    <- nrow(coords)
  n1   <- length(Y1)
  n2   <- length(Y2)
  nnew <- nrow(newcoords)
  
  d <- dist(rbind(newcoords, coords))
  D <- matrix(0, nrow = nnew + n, ncol = nnew + n)
  D[lower.tri(D)] <- d
  D <- D + t(D)
  R <- exp(- 3 * 0.6 / max(d) * D)
  Rinv <- solve(R[nnew + 1:n, nnew + 1:n])
  R11 <- R[1:nnew, 1:nnew]
  R12 <- R[1:nnew, nnew + 1:n]
  R21 <- t(R12)
  R12Rinv <- R12 %*% Rinv
  Rnew <- R11 - R12Rinv %*% R21
  
  XtRinv  <- t(X) %*% Rinv
  XtRinvX <- XtRinv %*% X
  
  keep   <- matrix(nrow = max.iter / n.thin, 
                   ncol = k + 5 + n + nnew)
  
  if (missing(inits) | any(is.na(inits))) {
    params <- rep(0, k + 5 + n + nnew)
    params[k + c(1:2, 5)] <- 1
  } else {
    params <- inits
  }
  
  I2 <- diag(2)
  Ik <- diag(k)
  In <- diag(n)
  M  <- rep(NA, n)
  Z  <- rep(NA, n)
  Z[ind1] <- Y1
  
  for (b in (-n.burnin + 1):max.iter) {
    if (b %% 1000 == 0)
      print(paste("Iteration number:", b))
    
    # beta
    delta <- solve(XtRinvX * params[k + 1] + Ik * 1e-04)
    chi   <- XtRinv %*% params[k + 5 + 1:n] * params[k + 1]
    params[1:k] <- mvtnorm::rmvnorm(1, delta %*% chi, delta, 
      checkSymmetry = FALSE)
    #params[1] <- 6
    Xb <- X %*% params[1:k]
    
    # prec true
    vect <- params[k + 5 + 1:n] - Xb
    params[k + 1] <- rgamma(1, n / 2 + 2, t(vect) %*% Rinv %*% vect / 2 + 1)
    #params[k + 1] <- 1/9
      
    # prec 1
    vect <- Y1 - params[k + 5 + ind1]
    params[k + 2] <- rgamma(1, n1 / 2 + 2, t(vect) %*% vect / 2 + 1)
    #params[k + 2] <- 1/2^2
      
    # lambda
    XX    <- cbind(1, params[k + 5 + ind2])
    delta <- solve(t(XX) %*% XX * params[k + 5] + I2 * 1e-04)
    chi   <- t(XX) %*% Y2 * params[k + 5] + c(0, 1) * 1e-04
    params[k + 3:4] <- mvtnorm::rmvnorm(1, delta %*% chi, delta, 
      checkSymmetry = FALSE)
    #params[k + 3:4] <- c(5, 2)
      
    # prec 2
    vect <- Y2 - XX %*% params[k + 3:4]
    params[k + 5] <- rgamma(1, n2 / 2 + 2, t(vect) %*% vect / 2 + 1)
    #params[k + 5] <- 1
    
    # Y true
    Z[ind2] <- (Y2 - params[k + 3]) / params[k + 4]
    M[ind1] <- params[k + 2]
    M[ind2] <- params[k + 5] * params[k + 4]^2
    delta <- Rinv * params[k + 1]
    diag(delta) <- diag(delta) + M
    delta <- solve(delta)
    chi   <- M * Z + Rinv %*% Xb * params[k + 1]
    params[k + 5 + 1:n] <- mvtnorm::rmvnorm(1, delta %*% chi, delta,
      checkSymmetry = FALSE)

    if ((b > 0) && (b %% n.thin == 0)) {
      # Y true interpolation
      mu <- newX %*% params[1:k] + R12Rinv %*% (params[k + 5 + 1:n] - Xb)
      params[k + 5 + n + 1:nnew] <- mvtnorm::rmvnorm(1, mu, Rnew / params[k + 1],
        checkSymmetry = FALSE)
    
      keep[b / n.thin, ] <- params
    }
  }
  
  colnames(keep) <- c(paste0("beta", 1:k - 1), 
                      "precTrue", "prec1",
                      "lambda0", "lambda1", "prec2",
                      paste0("Ytrue", 1:n), 
                      paste0("YtruePred", 1:nnew))
  
  return(keep)
}

# GRID
nx <- ny <- 15
nxy <- nx * ny
grid <- expand.grid(x = 1:nx, y = 1:ny)

d <- dist(grid)
D <- matrix(0, nrow = nxy, ncol = nxy)
D[lower.tri(D)] <- d
D <- D + t(D)
R <- exp(- 3 * 0.6 / max(d) * D)

# FUNCTION
sims <- function(iter) {
  set.seed(3 * iter + 123)
  
  # TRAIN AND TEST DATA
  ind  <- 1:nxy
  ind1 <- sort(sample(ind, size = nxy * .4))
  ind2 <- sort(sample(ind[-ind1], size = nxy * .4))
  out  <- sort(ind[-c(ind1, ind2)])
  
  index <- rep(NA, nxy * .8)
  index[match(ind1, sort(c(ind1, ind2)))] <- 1
  index[match(ind2, sort(c(ind1, ind2)))] <- 2
  
  # PARAMETERS and DATA
  X     <- cbind(rep(1, nxy))
  beta  <- c(6)
  Ytrue <- c(mvtnorm::rmvnorm(1, X %*% beta, 9 * R))
  Y1    <- rnorm(nxy, Ytrue, 2)
  Y2    <- rnorm(nxy, 5 + 2 * Ytrue, 1)
  
  # FIT MODELS
  m1 <- list()
  k <- 0
  repeat {
    m1[[1]] <- model1(Y1[ind1], X[ind1,,drop=FALSE], grid[ind1,],
                      X[out,,drop=FALSE], grid[out,],
                      if(k < 1) NA else tail(m1[[1]], 1),
                      ifelse(k < 1, 25000, 0), 25000, 25)
    m1[[2]] <- model1(Y1[ind1], X[ind1,,drop=FALSE], grid[ind1,],
                      X[out,,drop=FALSE], grid[out,],
                      if(k < 1) NA else tail(m1[[2]], 1),
                      ifelse(k < 1, 25000, 0), 25000, 25)
    
    if(
      (min(coda::effectiveSize(m1[[1]]) + coda::effectiveSize(m1[[2]])) > 200) &
      (coda::gelman.diag(coda::mcmc.list(coda::mcmc(m1[[1]]), coda::mcmc(m1[[2]])))$mpsrf < 1.2)
    ){
      break
    }
    
    k <- k + 1
  }
  
  m2 <- list()
  k <- 0
  repeat {
    m2[[1]] <- model2(Y1[ind1], Y2[ind2], index, X[-out,,drop=FALSE], 
                      grid[-out,], X[out,,drop=FALSE], grid[out,], 
                      if(k < 1) NA else tail(m2[[1]], 1),
                      ifelse(k < 1, 25000, 0), 25000, 25) 
    m2[[2]] <- model2(Y1[ind1], Y2[ind2], index, X[-out,,drop=FALSE], 
                      grid[-out,], X[out,,drop=FALSE], grid[out,],
                      if(k < 1) NA else tail(m2[[2]], 1),
                      ifelse(k < 1, 25000, 0), 25000, 25) 

    if(
      (min(coda::effectiveSize(m2[[1]]) + coda::effectiveSize(m2[[2]])) > 200) &
      (coda::gelman.diag(coda::mcmc.list(coda::mcmc(m2[[1]]), coda::mcmc(m2[[2]])))$mpsrf < 1.2)
    ){
      break
    }
    
    k <- k + 1
  }
  
  # COMPUTE METRICS
  metrics <- rep(NA, 8)
  metrics[1:4] <- ZoopFit::validation(Ytrue[out], rbind(
    m1[[1]][,paste0("YtruePred", 1:(nxy * .2))],
    m1[[2]][,paste0("YtruePred", 1:(nxy * .2))]))
  metrics[5:8] <- ZoopFit::validation(Ytrue[out], rbind(
    m2[[1]][,paste0("YtruePred", 1:(nxy * .2))],
    m2[[2]][,paste0("YtruePred", 1:(nxy * .2))]))
                      
  return(metrics)
}

# PARALLEL COMPUTING
time <- Sys.time()
n.cores <- max(parallel::detectCores() - 2, 2)
cl <- parallel::makeCluster(n.cores)
parallel::clusterExport(cl, c("model1", "model2", "nxy", "grid", "R"))
metrics <- parallel::parSapply(cl = cl, X = 1:100, FUN = sims)
parallel::stopCluster(cl = cl)
Sys.time() - time
