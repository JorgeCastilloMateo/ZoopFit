### PAPER EES -- SECTION 3 ###
###   SIMULATION EXAMPLE   ###

# PACKAGES
library("coda")
library("ZoopFit")
library("mvtnorm")

# GIBBS SAMPLER model 1: Spatial LM
# Y = Ytrue + e1, 
# Ytrue = X %*% beta + W
model1 <- function(Y, X, coords, Xnew, newcoords, 
                   inits,
                   n.burnin = 10, max.iter = 10, n.thin = 1) {
  
  k    <- ncol(X)
  n    <- nrow(coords)
  nnew <- nrow(newcoords)
  
  d <- dist(rbind(newcoords, coords))
  D <- matrix(0, nrow = nnew + n, ncol = nnew + n)
  D[lower.tri(D)] <- d
  D <- D + t(D)
  R <- exp(- 3 / 0.5 * D)
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
    params <- rnorm(k + 2 + n + nnew, mean = 0, sd = 10)
    params[k + 1:2] <- 1 / rgamma(2, shape = 0.1, rate = 0.1)
  } else {
    params <- inits
  }
  
  Ik <- diag(k)
  In <- diag(n)
  
  for (b in (-n.burnin + 1):max.iter) {
    if (b %% 1000 == 0)
      print(paste("Iteration number:", b))
    
    # beta
    delta <- solve(XtRinvX * params[k + 1] + Ik * 1e-02)
    chi   <- XtRinv %*% params[k + 2 + 1:n] * params[k + 1]
    params[1:k] <- mvtnorm::rmvnorm(1, delta %*% chi, delta, 
                                    checkSymmetry = FALSE)
    Xb <- X %*% params[1:k]
    
    # prec true
    vect <- params[k + 2 + 1:n] - Xb
    params[k + 1] <- rgamma(1, n / 2 + 0.1, t(vect) %*% Rinv %*% vect / 2 + 0.1)
    
    # prec 1
    vect <- Y - params[k + 2 + 1:n]
    params[k + 2] <- rgamma(1, n / 2 + 0.1, t(vect) %*% vect / 2 + 0.1)
    
    # Y true
    delta <- solve(In * params[k + 2] + Rinv * params[k + 1])
    chi   <- Y * params[k + 2] + Rinv %*% Xb * params[k + 1]
    params[k + 2 + 1:n] <- mvtnorm::rmvnorm(1, delta %*% chi, delta,
                                            checkSymmetry = FALSE)
    
    if ((b > 0) && (b %% n.thin == 0)) {
      # Y true interpolation
      mu <- Xnew %*% params[1:k] + R12Rinv %*% (params[k + 2 + 1:n] - Xb)
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

# GIBBS SAMPLER model 2: Spatial Calibration
# Y1 = Ytrue + e1, 
# Y2 = lambda0 + lambda1 * Ytrue + e2,
# Ytrue = X %*% beta + W
Rcpp::cppFunction('
  arma::mat model2Rcpp(
    arma::vec Y1, arma::vec Y2,
    arma::mat X, arma::mat Xnew, arma::mat keep,
    int k, int n, int n1, int n2, int nnew,
    arma::uvec ind1, arma::uvec ind2,
    arma::mat Rinv, arma::mat Rnew,
    arma::mat XtRinvX, arma::mat XtRinv, arma::mat R12Rinv,
    arma::vec params,
    int nBurnin, int maxIter, int nThin
  ) {
  
    arma::mat I2 = arma::mat(2, 2, arma::fill::eye);
    arma::mat Ik = arma::mat(k, k, arma::fill::eye);
    arma::mat XX(n2, 2, arma::fill::ones);
    arma::vec Z(n);
    Z(ind1) = Y1;
    arma::vec Xb(n);
    arma::vec vectn(n);
    arma::vec vectn1(n1);
    arma::vec vectn2(n2);
    arma::vec l0 = { 0, 1 };
    arma::vec M(n);
  
    arma::mat delta2(2, 2);
    arma::vec chi2(2);
    arma::mat deltak(k, k);
    arma::vec chik(k);
    arma::mat deltan(n, n);
    arma::vec chin(n);
    
    arma::mat cholRnew = arma::chol(Rnew, "lower");
  
    for (int b = 1 - nBurnin; b <= maxIter; ++b) {
      
      if (b % 1000 == 0) {
        Rcpp::Rcout << "The value of b : " << b << "\\n";
      }
      
      // beta
      deltak = arma::inv_sympd(XtRinvX * params(k) + Ik * 1e-02);
      chik   = XtRinv * params(arma::span(k + 5, k + 4 + n)) * params(k);
      params.head(k) = deltak * chik + arma::chol(deltak, "lower") * arma::randn(k);
      Xb = X * params.head(k);
    
      // prec true
      vectn     = params(arma::span(k + 5, k + 4 + n)) - Xb;
      params(k) = arma::randg(arma::distr_param(
        n / 2.0 + 0.1, 
        1 / arma::as_scalar(vectn.t() * Rinv * vectn / 2 + 0.1)
      ));

      // prec 1
      vectn1        = Y1 - params(k + 5 + ind1);
      params(k + 1) = arma::randg(arma::distr_param(
        n1 / 2.0 + 0.1, 
        1 / arma::as_scalar(vectn1.t() * vectn1 / 2 + 0.1)
      ));

      // lambda
      XX.col(1) = params(k + 5 + ind2);
      delta2 = arma::inv_sympd(XX.t() * XX * params(k + 4) + I2 * 1e-02);
      chi2   = XX.t() * Y2 * params(k + 4) + l0 * 1e-02;
      params(arma::span(k + 2, k + 3)) = 
        delta2 * chi2 + arma::chol(delta2, "lower") * arma::randn(2);

      // prec 2
      vectn2        = Y2 - XX * params(arma::span(k + 2, k + 3));
      params(k + 4) = arma::randg(arma::distr_param(
        n2 / 2.0 + 0.1, 
        1 / arma::as_scalar(vectn2.t() * vectn2 / 2 + 0.1)
      ));
    
      // Y true
      Z(ind2) = (Y2 - params(k + 2)) / params(k + 3);
      M(ind1).fill(params(k + 1));
      M(ind2).fill(params(k + 4) * params(k + 3) * params(k + 3));
      deltan = Rinv * params(k);
      deltan.diag() += M;
      deltan = arma::inv_sympd(deltan);
      chin   = M % Z + Rinv * Xb * params(k);
      params(arma::span(k + 5, k + 4 + n)) = 
        deltan * chin + arma::chol(deltan, "lower") * arma::randn(n);
    
      if ((b > 0) && (b % nThin == 0)) {
        // Y true interpolation
        params(arma::span(k + 5 + n, k + 4 + n + nnew)) = 
          Xnew * params.head(k) + R12Rinv * (params(arma::span(k + 5, k + 4 + n)) - Xb) +
          cholRnew * arma::randn(nnew) / sqrt(params(k));
      
        keep.row(b / nThin - 1) = params.t();
      }
    }
    
    return keep;
  
}', depends = "RcppArmadillo")

model2 <- function(Y1, Y2, ind, X, coords, Xnew, newcoords,
                   inits,
                   n.burnin = 10, max.iter = 10, n.thin = 1) {
  
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
  R <- exp(- 3 / 0.5 * D)
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
    params <- rnorm(k + 5 + n + nnew, mean = 0, sd = 10)
    params[k + c(1:2, 5)] <- 1 / rgamma(3, shape = 0.1, rate = 0.1)
  } else {
    params <- inits
  }
  
  keep <- model2Rcpp(
    Y1, Y2, X, Xnew, keep,
    k, n, n1, n2, nnew,
    ind1 - 1, ind2 - 1,
    Rinv, Rnew,
    XtRinvX, XtRinv, R12Rinv,
    params,
    n.burnin, max.iter, n.thin)
  
  colnames(keep) <- c(paste0("beta", 1:k - 1), 
                      "precTrue", "prec1",
                      "lambda0", "lambda1", "prec2",
                      paste0("Ytrue", 1:n), 
                      paste0("YtruePred", 1:nnew))
  
  return(keep)
}

sims <- function(iter) {
  set.seed(3* iter + 123)
  
  # GRID
  nx <- ny <- 10
  nxy <- nx * ny
  grid <- cbind(runif(nxy), runif(nxy))
  
  d <- dist(grid)
  D <- matrix(0, nrow = nxy, ncol = nxy)
  D[lower.tri(D)] <- d
  D <- D + t(D)
  R <- exp(- 3 / 0.5 * D)
  
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
  beta  <- 4
  Ytrue <- c(mvtnorm::rmvnorm(1, X %*% beta, 4 * R))
  Y1    <- rnorm(nxy, Ytrue, 2)
  Y2    <- rnorm(nxy, 2 + 2 * Ytrue, 1)
  
  # FIT MODELS
  m1 <- list()
  k <- 0
  repeat {
    m1[[1]] <- model1(Y1[ind1], X[ind1,,drop=FALSE], grid[ind1,],
                      X[out,,drop=FALSE], grid[out,],
                      if(k < 1) NA else tail(m1[[1]], 1),
                      ifelse(k < 1, 25000, 0), (k + 1) * 25000, (k + 1) * 25)
    m1[[2]] <- model1(Y1[ind1], X[ind1,,drop=FALSE], grid[ind1,],
                      X[out,,drop=FALSE], grid[out,],
                      if(k < 1) NA else tail(m1[[2]], 1),
                      ifelse(k < 1, 25000, 0), (k + 1) * 25000, (k + 1) * 25)
    
    print(max(coda::gelman.diag(coda::mcmc.list(coda::mcmc(m1[[1]]), coda::mcmc(m1[[2]])))$psrf[,2]))
    
    if(
      (min(coda::effectiveSize(m1[[1]]) + coda::effectiveSize(m1[[2]])) > 200) &
      (max(coda::gelman.diag(coda::mcmc.list(coda::mcmc(m1[[1]]), coda::mcmc(m1[[2]])))$psrf[,2]) < 1.1)
    ){
      break
    }
    
    k <- min(k + 1, 3)
  }
  
  m2 <- list()
  k <- 0
  repeat {
    m2[[1]] <- model2(Y1[ind1], Y2[ind2], index, X[-out,,drop=FALSE], 
                      grid[-out,], X[out,,drop=FALSE], grid[out,], 
                      if(k < 1) NA else tail(m2[[1]], 1),
                      ifelse(k < 1, 25000, 0), (k + 1) * 25000, (k + 1) * 25) 
    m2[[2]] <- model2(Y1[ind1], Y2[ind2], index, X[-out,,drop=FALSE], 
                      grid[-out,], X[out,,drop=FALSE], grid[out,],
                      if(k < 1) NA else tail(m2[[2]], 1),
                      ifelse(k < 1, 25000, 0), (k + 1) * 25000, (k + 1) * 25) 
    
    print(max(coda::gelman.diag(coda::mcmc.list(coda::mcmc(m2[[1]]), coda::mcmc(m2[[2]])))$psrf[,2]))
    
    if(
      (min(coda::effectiveSize(m2[[1]]) + coda::effectiveSize(m2[[2]])) > 200) &
      (max(coda::gelman.diag(coda::mcmc.list(coda::mcmc(m2[[1]]), coda::mcmc(m2[[2]])))$psrf[,2]) < 1.1)
    ){
      break
    }
    
    k <- min(k + 1, 3)
  }
  
  metrics <- rep(NA, 8)
  metrics[1:4] <- unlist(ZoopFit::validation(Ytrue[out], rbind(
    m1[[1]][,paste0("YtruePred", 1:(nxy * .2))],
    m1[[2]][,paste0("YtruePred", 1:(nxy * .2))])))
  metrics[5:8] <- unlist(ZoopFit::validation(Ytrue[out], rbind(
    m2[[1]][,paste0("YtruePred", 1:(nxy * .2))],
    m2[[2]][,paste0("YtruePred", 1:(nxy * .2))])))
  
  return(metrics)
}

metrics <- matrix(nrow = 100, ncol = 8)
for (i in 1:100) {
  print(paste0("Sim: ", i))
  time <- Sys.time()
  metrics[i,] <- sims(i)
  print(Sys.time() - time)
}

## parallel requires to define model2Rcpp within the package
#time <- Sys.time()
#n.cores <- parallel::detectCores() - 2
#cl <- parallel::makeCluster(n.cores)
#parallel::clusterExport(cl, c("model1", "model2", "nxy", "grid", "R"))
#metrics <- parallel::parSapply(cl = cl, X = 1:100, FUN = sims)
#parallel::stopCluster(cl = cl)
#Sys.time() - time

saveRDS(metrics, file = "simsMetrics.rds")
