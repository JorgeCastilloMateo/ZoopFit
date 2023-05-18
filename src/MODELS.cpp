#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

///////////////////
//TRUNCATE NORMAL//
///////////////////

// [[Rcpp::export]]
double rtnorm(double mu, double sigma, double a, double b) {
  double alpha = (a - mu) / sigma;
  double beta  = (b - mu) / sigma;
  double pAlpha = R::pnorm(alpha, 0, 1, 1, 0);
  double pBeta  = R::pnorm(beta, 0, 1, 1, 0);
  return R::qnorm(pAlpha + R::runif(0, 1) * (pBeta - pAlpha), 0, 1, 1, 0) * sigma + mu;
}

// [[Rcpp::export]]
double dtnorm(double x, double mu, double sigma, double a, double b) { // log
  return R::dnorm(x, mu, sigma, 1) - log(R::pnorm(b, mu, sigma, 1, 0) - R::pnorm(a, mu, sigma, 1, 0));
}

// [[Rcpp::export]]
double rfnorm(double mu, double sigma) {
  double alpha = - mu / sigma;
  double pAlpha = R::pnorm(alpha, 0, 1, 1, 0);
  return R::qnorm(pAlpha + R::runif(0, 1) * (1 - pAlpha), 0, 1, 1, 0) * sigma + mu;
}

// [[Rcpp::export]]
arma::mat rWishart(arma::mat V, int p, double n) {
  
  // int p = arma::size(V)(0);
  arma::mat A(p, p, arma::fill::zeros);
  arma::mat U(p, p);
  arma::mat X(p, p);
  
  int j;
  for (int i = 0; i < p; ++i) {
    j = 0;
    A(i,i) = sqrt(R::rchisq(n - i));
    while (j < i) {
      A(i,j++) = arma::randn();
    }
  }
  
  U = arma::chol(V); // upper by default
  
  X = U.t() * A * A.t() * U;
  
  return X;
}

//////////////////
//    NORMAL    //
//////////////////

//' COMPUTE PRODUCT OF NORMAL DENSITIES
//' @param mu vector mean
//' @param prec vector precision
//' @return vector(0) mean, vector(1) variance
// [[Rcpp::export]]
arma::vec productNormals(arma::vec mu, arma::vec prec) {
  arma::vec x(2);
  x(1) = 1 / arma::accu(prec); // variance
  x(0) = arma::accu(mu % prec) * x(1); // mean
  return x;
}

//////////////////
//      GP      //
//////////////////

// [[Rcpp::export]]
double meanGP(arma::vec vect, double sigma2inv, arma::mat Rinv, double a, double b2inv) {
  
  int n = vect.n_elem;
  arma::vec one(n, arma::fill::ones);
  
  double sd2 = 1 / (arma::accu(Rinv) * sigma2inv + b2inv);
  double mu = ((one.t() * Rinv * vect).eval()(0,0) * sigma2inv + a * b2inv) * sd2;
  
  return R::rnorm(mu, sqrt(sd2));
}

// [[Rcpp::export]]
double precGP(arma::vec vect, double mean, arma::mat Rinv, double a, double b) {
  
  vect -= mean;
  
  double ga = vect.n_elem / 2 + a;
  double gb = (vect.t() * Rinv * vect).eval()(0,0) / 2 + b;
  
  return R::rgamma(ga, 1 / gb);
}

//////////////////
// FINAL MODEL  //
//////////////////

// SIMPLE LINEAR REGRESSION (x is the covariate)
// [[Rcpp::export]]
arma::mat lmS(
    arma::vec Y, arma::vec X, arma::vec inits,
    double na, double nb, double ga, double gb,
    int N, int nBurnin, int maxIter, int nThin, int nReport) {
  
  arma::mat keep(maxIter / nThin, 3);
  arma::vec parameters = inits; // beta0 beta1 prec
  
  arma::vec Z = Y - (parameters(0) + parameters(1) *  X);
  
  double delta;
  double chi;
  double na1 = 1;
  
  // ITERATIONS //
  for (int b = 1 - nBurnin; b <= maxIter; ++b) {
    
    // beta0
    Z += parameters(0);
    delta = 1 / (N * parameters(2) + nb);
    chi   = arma::accu(Z) * parameters(2) + na * nb;
    parameters(0) = R::rnorm(delta * chi, sqrt(delta));
    
    // beta1
    Z += parameters(1) * X - parameters(0);
    delta = 1 / (arma::accu(pow(X, 2)) * parameters(2) + nb);
    chi   = (Z.t() * X).eval()(0) * parameters(2) + na1 * nb;
    parameters(1) = R::rnorm(delta * chi, sqrt(delta));
    
    // prec
    Z -= parameters(1) * X;
    parameters(2) = R::rgamma(N / 2 + ga, 1 / ((Z.t() * Z).eval()(0) / 2 + gb));
    
    //
    if ((b != 0) && (b % nReport == 0)) {
      Rcpp::Rcout << "The value of b : " << b << "\n";
    }
    
    if ((b > 0) && (b % nThin == 0)) {
      keep.row(b / nThin - 1) = parameters.t();
    }
  }
  
  return keep;
}

// MULTIPLE LINEAR REGRESSION (x is the design matrix, p number of covariates)
// [[Rcpp::export]]
arma::mat lmM(
    arma::vec Y, arma::mat X, arma::vec inits,
    double na, double nb, double ga, double gb,
    int p, int N, 
    int nBurnin, int maxIter, int nThin, int nReport) {
  
  arma::mat keep(maxIter / nThin, p + 2);
  arma::vec parameters = inits; // beta0 beta1 ... betap prec
  
  arma::vec Z = Y - (X * parameters(arma::span(0, p)));
  
  double delta;
  double chi;
  
  // ITERATIONS //
  for (int b = 1 - nBurnin; b <= maxIter; ++b) {
    
    // betaj
    for (int j = 0; j <= p; ++j) {
      Z += parameters(j) * X.col(j);
      delta = 1 / (arma::accu(pow(X.col(j), 2)) * parameters(p + 1) + nb);
      chi   = (Z.t() * X.col(j)).eval()(0) * parameters(p + 1) + na * nb;
      parameters(j) = R::rnorm(delta * chi, sqrt(delta));
      Z -= parameters(j) * X.col(j);
    }
    
    // prec
    parameters(p + 1) = R::rgamma(N / 2 + ga, 1 / ((Z.t() * Z).eval()(0) / 2 + gb));
    
    //
    if ((b != 0) && (b % nReport == 0)) {
      Rcpp::Rcout << "The value of b : " << b << "\n";
    }
    
    if ((b > 0) && (b % nThin == 0)) {
      keep.row(b / nThin - 1) = parameters.t();
    }
  }
  
  return keep;
}

// SPATIAL MULTIPLE LINEAR REGRESSION 
// (x is the design matrix without intercept, p number of covariates)
// (decay parameters known)
// [[Rcpp::export]]
arma::mat splmM(
    arma::vec Y, arma::mat X, 
    arma::mat I, arma::uvec si, arma::mat Rb0, double Rb0sum,
    arma::vec inits, // beta0s's beta0 precbeta0 beta1 ... betap prec
    double na, double nb, double ga, double gb,
    int p, int n, int N, 
    int nBurnin, int maxIter, int nThin, int nReport) {
  
  arma::mat keep(maxIter / nThin, n + p + 3);
  arma::vec beta0 = inits(arma::span(0, n - 1));
  arma::vec hiper = inits(arma::span(n, n + 1));// beta0 precbeta0
  arma::vec parameters = inits(arma::span(n + 2, n + p + 2)); // beta1 ... betap prec
  
  arma::vec Z = Y - (beta0.elem(si) + X * parameters(arma::span(0, p - 1)));
  
  // HELP //
  arma::vec onen(n, arma::fill::ones);

  double delta;
  double chi;
  arma::vec process;
  arma::mat r;
  
  double sumita;
  
  // ITERATIONS //
  for (int b = 1 - nBurnin; b <= maxIter; ++b) {
    
    // beta0 si i=1,...,n
    Z += beta0.elem(si);
    for (int iInd = 0; iInd < n; ++iInd){
      process = beta0;
      process.shed_row(iInd);
      r = Rb0.row(iInd);
      r.shed_col(iInd);
      delta = 1 / (arma::accu(I.col(iInd)) * parameters(p) + Rb0(iInd, iInd) * hiper(1));
      chi   = (Z.t() * I.col(iInd)).eval()(0) * parameters(p) +
        (hiper(0) * Rb0(iInd, iInd) + r * (hiper(0) - process)).eval()(0,0) * hiper(1);
      beta0(iInd) = R::rnorm(delta * chi, sqrt(delta));
    }
    Z -= beta0.elem(si);
    
    // beta0
    delta = 1 / (Rb0sum * hiper(1) + nb); 
    chi   = (onen.t() * Rb0 * beta0).eval()(0,0) * hiper(1) + na * nb;
    hiper(0) = R::rnorm(delta * chi, sqrt(delta));
    
    // precbeta0
    process = beta0 - hiper(0);
    sumita = (process.t() * Rb0 * process).eval()(0,0); 
    hiper(1) = R::rgamma(n / 2 + ga, 1 / (sumita / 2 + gb));
    
    // betaj
    for (int j = 0; j < p; ++j) {
      Z += parameters(j) * X.col(j);
      delta = 1 / (arma::accu(pow(X.col(j), 2)) * parameters(p) + nb);
      chi   = (Z.t() * X.col(j)).eval()(0) * parameters(p) + na * nb;
      parameters(j) = R::rnorm(delta * chi, sqrt(delta));
      Z -= parameters(j) * X.col(j);
    }
    
    // prec
    parameters(p) = R::rgamma(N / 2 + ga, 1 / ((Z.t() * Z).eval()(0) / 2 + gb));
    
    //
    if ((b != 0) && (b % nReport == 0)) {
      Rcpp::Rcout << "The value of b : " << b << "\n";
    }
    
    if ((b > 0) && (b % nThin == 0)) {
      keep(b / nThin - 1, arma::span(0, n - 1)) = beta0.t();
      keep(b / nThin - 1, arma::span(n, n + 1)) = hiper.t();
      keep(b / nThin - 1, arma::span(n + 2, n + p + 2)) = parameters.t();
    }
  }
  
  return keep;
}

// SPATIAL SIMPLE LINEAR REGRESSION 
// (x is the covariate)
// (decay parameters UNknown)
// [[Rcpp::export]]
arma::mat splmS(
    arma::vec Y, arma::vec X, 
    arma::mat I, arma::uvec si, 
    arma::mat dist,
    arma::vec inits, // beta0s's beta1 beta0 precbeta0 decaybeta0 prec
    double sd,
    double na, double nb, double ga, double gb, double Umin, double Umax,
    int n, int N, 
    int nBurnin, int maxIter, int nThin, int nReport) {
  
  // DATA //
  arma::mat keep(maxIter / nThin, n + 5);
  arma::vec beta0 = inits(arma::span(0, n - 1));
  arma::vec hiper = inits(arma::span(n + 1, n + 3));// beta0 precbeta0 decaybeta0
  arma::vec parameters(2); // beta1 prec
  parameters(0) = inits(n);
  parameters(1) = inits(n + 4);
  
  arma::vec Z = Y - (beta0.elem(si) + parameters(0) * X);
  
  // HELP //
  arma::vec onen(n, arma::fill::ones);
  arma::vec oneN(N, arma::fill::ones);
  
  arma::mat Delta;
  arma::mat Chi;
  double delta;
  double chi;
  arma::vec process;
  arma::mat r;
  double A;
  
  arma::vec aux(N);
  
  double phiAux;
  arma::mat RAux(n, n);
  double RlogdetAux;
  arma::mat Rb0    = arma::inv_sympd(exp(- hiper(2) * dist));
  double Rb0logdet = arma::log_det_sympd(Rb0);

  double sumita;
  double sumitaAux;
  
  double na1 = 1;
  
  // ITERATIONS //
  for (int b = 1 - nBurnin; b <= maxIter; ++b) {
    
    // beta0 si i=1,...,n
    Z += beta0.elem(si);
    for (int iInd = 0; iInd < n; ++iInd){
      process = beta0;
      process.shed_row(iInd);
      r = Rb0.row(iInd);
      r.shed_col(iInd);
      delta = 1 / (arma::accu(I.col(iInd)) * parameters(1) + Rb0(iInd, iInd) * hiper(1));
      chi   = (Z.t() * I.col(iInd)).eval()(0) * parameters(1) +
        (hiper(0) * Rb0(iInd, iInd) + r * (hiper(0) - process)).eval()(0,0) * hiper(1);
      beta0(iInd) = R::rnorm(delta * chi, sqrt(delta));
    }
    Z -= beta0.elem(si);
    
    // beta0
    delta = 1 / (arma::accu(Rb0) * hiper(1) + nb); 
    chi   = (onen.t() * Rb0 * beta0).eval()(0,0) * hiper(1) + na * nb;
    hiper(0) = R::rnorm(delta * chi, sqrt(delta));
    
    // precbeta0
    process = beta0 - hiper(0);
    sumita = (process.t() * Rb0 * process).eval()(0,0); 
    hiper(1) = R::rgamma(n / 2 + ga, 1 / (sumita / 2 + gb));
    
    // phi v0
    process = beta0;
    phiAux     = rtnorm(hiper(2), sd, Umin, Umax);
    RAux       = arma::inv_sympd(exp(- phiAux * dist));  

    RlogdetAux = arma::log_det_sympd(RAux);
    sumita    = (process.t() * Rb0 * process).eval()(0,0);  
    sumitaAux = (process.t() * RAux * process).eval()(0,0);  
    A = (RlogdetAux - sumitaAux) / 2 - 
      (Rb0logdet - sumita) / 2 + 
      dtnorm(hiper(2), phiAux, sd, Umin, Umax) - 
      dtnorm(phiAux, hiper(2), sd, Umin, Umax);
    if (log(R::runif(0, 1)) <= A) {
      hiper(2) = phiAux;
      Rb0 = RAux;
      Rb0logdet = RlogdetAux;
      sumita = sumitaAux;
    }
    
    // beta1
    Z += parameters(0) * X;
    delta = 1 / (arma::accu(pow(X, 2)) * parameters(1) + nb);
    chi   = (Z.t() * X).eval()(0) * parameters(1) + na1 * nb;
    parameters(0) = R::rnorm(delta * chi, sqrt(delta));
    Z -= parameters(0) * X;
    
    // prec
    parameters(1) = R::rgamma(N / 2 + ga, 1 / ((Z.t() * Z).eval()(0) / 2 + gb));
    
    //
    if ((b != 0) && (b % nReport == 0)) {
      Rcpp::Rcout << "The value of b : " << b << "\n";
    }
    
    if ((b > 0) && (b % nThin == 0)) {
      keep(b / nThin - 1, arma::span(0, n - 1)) = beta0.t();
      keep(b / nThin - 1, n) = parameters(0);
      keep(b / nThin - 1, arma::span(n + 1, n + 3)) = hiper.t();
      keep(b / nThin - 1, n + 4) = parameters(1);
    }
  }
  
  return keep;
}

// SPATIAL MULTIPLE LINEAR REGRESSION WITH RANDOM EFFECTS
// (x is the design matrix without intercept, p number of covariates)
// (decay parameters known)
// [[Rcpp::export]]
arma::mat splmMRnoAR(
    arma::vec Y, arma::mat X, 
    arma::mat I, arma::uvec si, arma::mat Rb0, double Rb0sum,
    arma::vec inits, // beta0s's beta0 precbeta0 precyear beta1 ... betap prec year1 ... yearT
    double na, double nb, double ga, double gb,
    int p, int n, int N, int T, arma::mat Xt, // T and t new & output T+1 elems new
    int nBurnin, int maxIter, int nThin, int nReport) {
  
  arma::mat keep(maxIter / nThin, n + p + T + 4);
  arma::vec beta0 = inits(arma::span(0, n - 1));
  arma::vec hiper = inits(arma::span(n, n + 2));// beta0 precbeta0 precyear
  arma::vec parameters = inits(arma::span(n + 3, n + p + 3)); // beta1 ... betap prec
  arma::vec year = inits(arma::span(n + p + 4, n + p + T + 3)); // random effects for each year
  
  arma::vec Z = Y - (beta0.elem(si) + X * parameters(arma::span(0, p - 1)) + Xt * year(arma::span(0, T - 1)));
  
  // HELP //
  arma::vec onen(n, arma::fill::ones);
  
  double delta;
  double chi;
  arma::vec process;
  arma::mat r;
  
  double sumita;
  
  // ITERATIONS //
  for (int b = 1 - nBurnin; b <= maxIter; ++b) {
    
    // beta0 si i=1,...,n
    Z += beta0.elem(si);
    for (int iInd = 0; iInd < n; ++iInd){
      process = beta0;
      process.shed_row(iInd);
      r = Rb0.row(iInd);
      r.shed_col(iInd);
      delta = 1 / (arma::accu(I.col(iInd)) * parameters(p) + Rb0(iInd, iInd) * hiper(1));
      chi   = (Z.t() * I.col(iInd)).eval()(0) * parameters(p) +
        (- r * process).eval()(0,0) * hiper(1);
      beta0(iInd) = R::rnorm(delta * chi, sqrt(delta));
    }
    Z -= beta0.elem(si);
    
    // beta0
    delta = 1 / (T * hiper(2) + nb); 
    chi   = arma::accu(year) * hiper(2) + na * nb;
    hiper(0) = R::rnorm(delta * chi, sqrt(delta));
    
    // precbeta0
    process = beta0;
    sumita = (process.t() * Rb0 * process).eval()(0,0); 
    hiper(1) = R::rgamma(n / 2 + ga, 1 / (sumita / 2 + gb));
    
    // betaj
    for (int j = 0; j < p; ++j) {
      Z += parameters(j) * X.col(j);
      delta = 1 / (arma::accu(pow(X.col(j), 2)) * parameters(p) + nb);
      chi   = (Z.t() * X.col(j)).eval()(0) * parameters(p) + na * nb;
      parameters(j) = R::rnorm(delta * chi, sqrt(delta));
      Z -= parameters(j) * X.col(j);
    }
    
    // prec
    parameters(p) = R::rgamma(N / 2 + ga, 1 / ((Z.t() * Z).eval()(0) / 2 + gb));
    
    // yeart
    for (int tInd = 0; tInd < T; ++tInd) {
      Z += Xt * year;
      delta = 1 / (arma::accu(Xt.col(tInd)) * parameters(p) + hiper(2));
      chi   = (Z.t() * Xt.col(tInd)).eval()(0) * parameters(p) + hiper(0) * hiper(2);
      year(tInd) = R::rnorm(delta * chi, sqrt(delta));
      Z -= Xt * year;
    }
    
    // precyear
    hiper(2) = R::rgamma(T / 2 + ga, 1 / (arma::accu(pow(year - hiper(0), 2)) / 2 + gb));
    
    //
    if ((b != 0) && (b % nReport == 0)) {
      Rcpp::Rcout << "The value of b : " << b << "\n";
    }
    
    if ((b > 0) && (b % nThin == 0)) {
      keep(b / nThin - 1, arma::span(0, n - 1)) = beta0.t();
      keep(b / nThin - 1, arma::span(n, n + 2)) = hiper.t();
      keep(b / nThin - 1, arma::span(n + 3, n + p + 3)) = parameters.t();
      keep(b / nThin - 1, arma::span(n + p + 4, n + p + T + 3)) = year.t();
    }
  }
  
  return keep;
}

// SPATIAL MULTIPLE LINEAR REGRESSION WITH AR RANDOM EFFECTS
// (x is the design matrix without intercept, p number of covariates)
// (decay parameters known)
// [[Rcpp::export]]
arma::mat splmMR(
    arma::vec Y, arma::mat X, 
    arma::mat I, arma::uvec si, arma::mat Rb0, double Rb0sum,
    arma::vec inits, // beta0s's beta0 precbeta0 precyear beta1 ... betap prec year1 ... yearT rhoyear
    double na, double nb, double ga, double gb,
    int p, int n, int N, int T, arma::mat Xt, // T and t new & output T+1 elems new
    int nBurnin, int maxIter, int nThin, int nReport) {
  
  arma::mat keep(maxIter / nThin, n + p + T + 5);
  arma::vec beta0 = inits(arma::span(0, n - 1));
  arma::vec hiper = inits(arma::span(n, n + 2));// beta0 precbeta0 precyear
  arma::vec parameters = inits(arma::span(n + 3, n + p + 3)); // beta1 ... betap prec 
  arma::vec year = inits(arma::span(n + p + 4, n + p + T + 3)); // random effects for each year
  double rho = inits(n + p + T + 4);

  arma::vec Z = Y - (beta0.elem(si) + X * parameters(arma::span(0, p - 1)) + Xt * year(arma::span(0, T - 1)));
  
  // HELP //
  arma::vec onen(n, arma::fill::ones);
  
  double delta;
  double chi;
  arma::vec process;
  arma::mat r;
  double rho2;
  double sumita;
  arma::vec innovation(T);
  
  // ITERATIONS //
  for (int b = 1 - nBurnin; b <= maxIter; ++b) {
    
    // beta0 si i=1,...,n
    Z += beta0.elem(si);
    for (int iInd = 0; iInd < n; ++iInd){
      process = beta0;
      process.shed_row(iInd);
      r = Rb0.row(iInd);
      r.shed_col(iInd);
      delta = 1 / (arma::accu(I.col(iInd)) * parameters(p) + Rb0(iInd, iInd) * hiper(1));
      chi   = (Z.t() * I.col(iInd)).eval()(0) * parameters(p) +
        (- r * process).eval()(0,0) * hiper(1);
      beta0(iInd) = R::rnorm(delta * chi, sqrt(delta));
    }
    Z -= beta0.elem(si);
    
    // precbeta0
    process = beta0;
    sumita = (process.t() * Rb0 * process).eval()(0,0); 
    hiper(1) = R::rgamma(n / 2 + ga, 1 / (sumita / 2 + gb));
    
    // betaj
    for (int j = 0; j < p; ++j) {
      Z += parameters(j) * X.col(j);
      delta = 1 / (arma::accu(pow(X.col(j), 2)) * parameters(p) + nb);
      chi   = (Z.t() * X.col(j)).eval()(0) * parameters(p) + na * nb;
      parameters(j) = R::rnorm(delta * chi, sqrt(delta));
      Z -= parameters(j) * X.col(j);
    }
    
    //// prec
    parameters(p) = R::rgamma(N / 2 + ga, 1 / ((Z.t() * Z).eval()(0) / 2 + gb));
    
    // beta0
    delta  = 1 / (((T - 1) * pow(1 - rho, 2) + 1) * hiper(2) + nb); 
    sumita = year(T - 1) + (1 - rho) * arma::accu(year(arma::span(1, T - 2))) - rho * year(0);
    chi    = (sumita * (1 - rho) + year(0)) * hiper(2) + na * nb;
    hiper(0) = R::rnorm(delta * chi, sqrt(delta));
    
    // precyear
    innovation = year - hiper(0);
    hiper(2) = R::rgamma(
      T / 2 + ga, 
      1 / ((arma::accu(pow(innovation(arma::span(1, T - 1)) - rho * innovation(arma::span(0, T - 2)), 2)) + pow(innovation(0), 2)) / 2 + gb));
    
    // rhoyear
    sumita = arma::accu(pow(innovation(arma::span(0, T - 2)), 2));
    rho = rtnorm(
      arma::accu(innovation(arma::span(0, T - 2)) % innovation(arma::span(1, T - 1))) / sumita, 
      1 / sqrt(hiper(2) * sumita), 
      -1, 1);
    
    // yeart
    rho2 = pow(rho, 2);
    // tInd = 0
    Z += Xt * year;
    delta = 1 / (arma::accu(Xt.col(0)) * parameters(p) + (1 + rho2) * hiper(2));
    chi   = (Z.t() * Xt.col(0)).eval()(0) * parameters(p) + (rho * (year(1) - hiper(0)) + (1 + rho2) * hiper(0)) * hiper(2);
    year(0) = R::rnorm(delta * chi, sqrt(delta));
    // tInd = 1, ..., T-2
    for (int tInd = 1; tInd < T - 1; ++tInd) {
      delta = 1 / (arma::accu(Xt.col(tInd)) * parameters(p) + (1 + rho2) * hiper(2));
      chi   = (Z.t() * Xt.col(tInd)).eval()(0) * parameters(p) + (rho * (year(tInd + 1) + year(tInd - 1)) + pow(1 - rho, 2) * hiper(0)) * hiper(2);
      year(tInd) = R::rnorm(delta * chi, sqrt(delta));
    }
    // tInd = T-1
    delta = 1 / (arma::accu(Xt.col(T - 1)) * parameters(p) + hiper(2));
    chi   = (Z.t() * Xt.col(T - 1)).eval()(0) * parameters(p) + (rho * year(T - 2) + (1 - rho) * hiper(0)) * hiper(2);
    year(T - 1) = R::rnorm(delta * chi, sqrt(delta));
    Z -= Xt * year;
    
    //
    if ((b != 0) && (b % nReport == 0)) {
      Rcpp::Rcout << "The value of b : " << b << "\n";
    }
    
    if ((b > 0) && (b % nThin == 0)) {
      keep(b / nThin - 1, arma::span(0, n - 1)) = beta0.t();
      keep(b / nThin - 1, arma::span(n, n + 2)) = hiper.t();
      keep(b / nThin - 1, arma::span(n + 3, n + p + 3)) = parameters.t();
      keep(b / nThin - 1, arma::span(n + p + 4, n + p + T + 3)) = year.t();
      keep(b / nThin - 1, n + p + T + 4) = rho;
    }
  }
  
  return keep;
}

// como la de arriba pero centrando en espacio
// [[Rcpp::export]]
arma::mat splmMR2(
    arma::vec Y, arma::mat X, 
    arma::mat I, arma::uvec si, arma::mat Rb0, double Rb0sum,
    arma::vec inits, // beta0s's beta0 precbeta0 precyear beta1 ... betap prec year1 ... yearT rhoyear
    double na, double nb, double ga, double gb,
    int p, int n, int N, int T, arma::mat Xt, // T and t new & output T+1 elems new
    int nBurnin, int maxIter, int nThin, int nReport) {
  
  arma::mat keep(maxIter / nThin, n + p + T + 5);
  arma::vec beta0 = inits(arma::span(0, n - 1));
  arma::vec hiper = inits(arma::span(n, n + 2));// beta0 precbeta0 precyear
  arma::vec parameters = inits(arma::span(n + 3, n + p + 3)); // beta1 ... betap prec 
  arma::vec year = inits(arma::span(n + p + 4, n + p + T + 3)); // random effects for each year
  double rho = inits(n + p + T + 4);
  
  arma::vec Z = Y - (beta0.elem(si) + X * parameters(arma::span(0, p - 1)) + Xt * year(arma::span(0, T - 1)));
  
  // HELP //
  arma::vec onen(n, arma::fill::ones);
  
  double delta;
  double chi;
  arma::vec process;
  arma::mat r;
  double rho2;
  double sumita;
  arma::vec innovation(T);
  
  // ITERATIONS //
  // lo comentado es si centro el espacio en lugar de los AR year en el intercepto
  for (int b = 1 - nBurnin; b <= maxIter; ++b) {
    
    // beta0 si i=1,...,n
    Z += beta0.elem(si);
    for (int iInd = 0; iInd < n; ++iInd){
      process = beta0;
      process.shed_row(iInd);
      r = Rb0.row(iInd);
      r.shed_col(iInd);
      delta = 1 / (arma::accu(I.col(iInd)) * parameters(p) + Rb0(iInd, iInd) * hiper(1));
      chi   = (Z.t() * I.col(iInd)).eval()(0) * parameters(p) +
        (hiper(0) * Rb0(iInd, iInd) + r * (hiper(0) - process)).eval()(0,0) * hiper(1);
      beta0(iInd) = R::rnorm(delta * chi, sqrt(delta));
    }
    Z -= beta0.elem(si);
    
    // beta0
    delta = 1 / (arma::accu(Rb0) * hiper(1) + nb); 
    chi   = (onen.t() * Rb0 * beta0).eval()(0,0) * hiper(1) + na * nb;
    hiper(0) = R::rnorm(delta * chi, sqrt(delta));
    
    // precbeta0
    process = beta0 - hiper(0);
    sumita = (process.t() * Rb0 * process).eval()(0,0); 
    hiper(1) = R::rgamma(n / 2 + ga, 1 / (sumita / 2 + gb));
    
    // betaj
    for (int j = 0; j < p; ++j) {
      Z += parameters(j) * X.col(j);
      delta = 1 / (arma::accu(pow(X.col(j), 2)) * parameters(p) + nb);
      chi   = (Z.t() * X.col(j)).eval()(0) * parameters(p) + na * nb;
      parameters(j) = R::rnorm(delta * chi, sqrt(delta));
      Z -= parameters(j) * X.col(j);
    }
    
    // prec
    parameters(p) = R::rgamma(N / 2 + ga, 1 / ((Z.t() * Z).eval()(0) / 2 + gb));
    
    // precyear
    hiper(2) = R::rgamma(
      T / 2 + ga, 
      1 / ((arma::accu(pow(year(arma::span(1, T - 1)) - rho * year(arma::span(0, T - 2)), 2)) + pow(year(0), 2)) / 2 + gb));
    
    // rhoyear
    sumita = arma::accu(pow(year(arma::span(0, T - 2)), 2));
    rho = rtnorm(
      (year(arma::span(1, T - 1)).t() * year(arma::span(0, T - 2))).eval()(0) / sumita, 
      1 / sqrt(hiper(2) * sumita), 
      -1, 1);
    
    // yeart
    rho2 = pow(rho, 2);
    // tInd = 0
    Z += Xt * year;
    delta = 1 / (arma::accu(Xt.col(0)) * parameters(p) + (1 + rho2) * hiper(2));
    chi   = (Z.t() * Xt.col(0)).eval()(0) * parameters(p) + rho * year(1) * hiper(2);
    year(0) = R::rnorm(delta * chi, sqrt(delta));
    Z -= Xt * year;
    // tInd = 1, ..., T-2
    for (int tInd = 1; tInd < T - 1; ++tInd) {
      Z += Xt * year;
      delta = 1 / (arma::accu(Xt.col(tInd)) * parameters(p) + (1 + rho2) * hiper(2));
      chi   = (Z.t() * Xt.col(tInd)).eval()(0) * parameters(p) + rho * (year(tInd + 1) + year(tInd - 1)) * hiper(2);
      year(tInd) = R::rnorm(delta * chi, sqrt(delta));
      Z -= Xt * year;
    }
    // tInd = T-1
    Z += Xt * year;
    delta = 1 / (arma::accu(Xt.col(T - 1)) * parameters(p) + hiper(2));
    chi   = (Z.t() * Xt.col(T - 1)).eval()(0) * parameters(p) + rho * year(T - 2) * hiper(2);
    year(T - 1) = R::rnorm(delta * chi, sqrt(delta));
    Z -= Xt * year;
    
    //
    if ((b != 0) && (b % nReport == 0)) {
      Rcpp::Rcout << "The value of b : " << b << "\n";
    }
    
    if ((b > 0) && (b % nThin == 0)) {
      keep(b / nThin - 1, arma::span(0, n - 1)) = beta0.t();
      keep(b / nThin - 1, arma::span(n, n + 2)) = hiper.t();
      keep(b / nThin - 1, arma::span(n + 3, n + p + 3)) = parameters.t();
      keep(b / nThin - 1, arma::span(n + p + 4, n + p + T + 3)) = year.t();
      keep(b / nThin - 1, n + p + T + 4) = rho;
    }
  }
  
  return keep;
}

// Calibration two GPs (slope and intercept)
// [[Rcpp::export]]
arma::mat splmSS(
  arma::vec Y, arma::vec X, arma::mat I, arma::uvec si,
  arma::mat dist, arma::vec init, double sd,
  double na, double nb, double ga, double gb, double Umin, double Umax,
  int n, int N, int nBurnin, int maxIter, int nThin, int nReport) {
  
  // PARAMETERS // // a0's a1's mean0 mean1 prec0 prec1 decay0 decay1 precT2 
  arma::mat keep(maxIter / nThin, 7 + 2 * n);
  arma::vec parameters(2 * n + 7);
  
  parameters(0) = init(2 * n);
  parameters(1) = init(2 * n + 2);
  parameters(2) = init(2 * n + 4);
  parameters(3) = init(2 * n + 1);
  parameters(4) = init(2 * n + 3);
  parameters(5) = init(2 * n + 5);
  parameters(6) = init(2 * n + 6);
  parameters(arma::span(7, 6 + n)) = init(arma::span(0, n - 1));
  parameters(arma::span(7 + n, 6 + 2 * n)) = init(arma::span(n, 2 * n - 1));
  
  arma::vec beta0 = parameters(arma::span(7, 6 + n));
  arma::vec beta1 = parameters(arma::span(7 + n, 6 + 2 * n));
  arma::vec Z = Y - (beta0.elem(si) + X % beta1.elem(si));
  
  // HELP //
  arma::vec onen(n, arma::fill::ones);
  arma::vec oneN(N, arma::fill::ones);
  
  arma::mat Delta;
  arma::mat Chi;
  double delta;
  double chi;
  arma::vec process;
  arma::mat r;
  double A;
  
  double phiAux;
  arma::mat RAux(n, n);
  double RlogdetAux;
  
  arma::mat Rb0    = arma::inv_sympd(exp(- parameters(2) * dist));
  double Rb0logdet = arma::log_det_sympd(Rb0);
  double Rb0sum    = arma::accu(Rb0);
  
  arma::mat Rb1    = arma::inv_sympd(exp(- parameters(5) * dist));
  double Rb1logdet = arma::log_det_sympd(Rb1);
  double Rb1sum    = arma::accu(Rb1);
  
  double sumita;
  double sumitaAux;
  
  double na1 = 1;
  
  // ITERATIONS //
  for (int b = 1 - nBurnin; b <= maxIter; ++b) {
    
    // beta0
    Delta = 1 / (Rb0sum * parameters(1) + nb); 
    Chi   = onen.t() * Rb0 * parameters(arma::span(7, 6 + n)) * parameters(1) + na * nb;
    parameters(0) = R::rnorm((Delta * Chi).eval()(0,0), sqrt(Delta(0,0)));
    
    // phi beta0
    process = parameters(arma::span(7, 6 + n)) - parameters(0);
    phiAux     = rtnorm(parameters(2), sd, Umin, Umax);
    RAux       = arma::inv_sympd(exp(- phiAux * dist));
    RlogdetAux = arma::log_det_sympd(RAux);
    sumita    = (process.t() * Rb0 * process).eval()(0,0);  
    sumitaAux = (process.t() * RAux * process).eval()(0,0);  
    A = (RlogdetAux - sumitaAux * parameters(1)) / 2 - 
      (Rb0logdet  - sumita    * parameters(1)) / 2 + 
      dtnorm(parameters(2), phiAux, sd, Umin, Umax) - 
      dtnorm(phiAux, parameters(2), sd, Umin, Umax);
    if (log(R::runif(0, 1)) <= A) {
      parameters(2) = phiAux;
      Rb0 = RAux;
      Rb0logdet = RlogdetAux;
      Rb0sum = arma::accu(Rb0);
      sumita = sumitaAux;
    }
    
    // sigma beta0
    parameters(1) = R::rgamma(n / 2 + ga, 1 / (sumita / 2 + gb));
    
    // beta1
    Delta = 1 / (Rb1sum * parameters(4) + nb); 
    Chi   = onen.t() * Rb1 * parameters(arma::span(7 + n, 6 + 2 * n)) * parameters(4) + na1 * nb;
    parameters(3) = R::rnorm((Delta * Chi).eval()(0,0), sqrt(Delta(0,0)));
    
    // phi beta1
    process = parameters(arma::span(7 + n, 6 + 2 * n)) - parameters(3);
    phiAux     = rtnorm(parameters(5), sd, Umin, Umax);
    RAux       = arma::inv_sympd(exp(- phiAux * dist));
    RlogdetAux = arma::log_det_sympd(RAux);
    sumita    = (process.t() * Rb1 * process).eval()(0,0);  
    sumitaAux = (process.t() * RAux * process).eval()(0,0);  
    A = (RlogdetAux - sumitaAux * parameters(4)) / 2 - 
      (Rb1logdet  - sumita    * parameters(4)) / 2 + 
      dtnorm(parameters(5), phiAux, sd, Umin, Umax) - 
      dtnorm(phiAux, parameters(5), sd, Umin, Umax);
    if (log(R::runif(0, 1)) <= A) {
      parameters(5) = phiAux;
      Rb1 = RAux;
      Rb1logdet = RlogdetAux;
      Rb1sum = arma::accu(Rb1);
      sumita = sumitaAux;
    }
    
    // sigma beta1
    parameters(4) = R::rgamma(n / 2 + ga, 1 / (sumita / 2 + gb));
    
    // sigma
    parameters(6) = R::rgamma(N / 2 + ga, 
               1 / (arma::accu(arma::pow(Z, 2)) / 2 + gb));
    
    // beta0 si i=1,...,n
    Z += beta0.elem(si);
    for (int iInd = 0; iInd < n; ++iInd){
      process = beta0;
      process.shed_row(iInd);
      r = Rb0.row(iInd);
      r.shed_col(iInd);
      delta = 1 / (arma::accu(I.col(iInd)) * parameters(6) + Rb0(iInd, iInd) * parameters(1));
      chi   = (Z.t() * I.col(iInd)).eval()(0) * parameters(6) +
        (parameters(0) * Rb0(iInd, iInd) + r * (parameters(0) - process)).eval()(0,0) * parameters(1);
      beta0(iInd) = R::rnorm(delta * chi, sqrt(delta));
    }
    parameters(arma::span(7, 6 + n)) = beta0;
    
    // beta1 si i=1,...,n
    Z += X % beta1.elem(si) - beta0.elem(si);
    for (int iInd = 0; iInd < n; ++iInd){
      process = beta1;
      process.shed_row(iInd);
      r = Rb1.row(iInd);
      r.shed_col(iInd);
      delta = 1 / (arma::accu(X % I.col(iInd)) * parameters(6) + Rb1(iInd, iInd) * parameters(4));
      chi   = (Z.t() * (X % I.col(iInd))).eval()(0) * parameters(6) +
        (parameters(3) * Rb1(iInd, iInd) + r * (parameters(3) - process)).eval()(0,0) * parameters(4);
      beta1(iInd) = R::rnorm(delta * chi, sqrt(delta));
    }
    parameters(arma::span(7 + n, 6 + 2 * n)) = beta1;  
    Z -= X % beta1.elem(si);
    
    //
    if ((b != 0) && (b % nReport == 0)) {
      Rcpp::Rcout << "The value of b : " << b << "\n";
    }
    
    if ((b > 0) && (b % nThin == 0)) {
      init(2 * n) = parameters(0);
      init(2 * n + 2) = parameters(1);
      init(2 * n + 4) = parameters(2);
      init(2 * n + 1) = parameters(3);
      init(2 * n + 3) = parameters(4);
      init(2 * n + 5) = parameters(5);
      init(2 * n + 6) = parameters(6);
      init(arma::span(0, n - 1)) = parameters(arma::span(7, 6 + n));
      init(arma::span(n, 2 * n - 1)) = parameters(arma::span(7 + n, 6 + 2 * n));
      keep.row(b / nThin - 1) = init.t();
    }
  }
  
  return keep;
}

// Calibration Coregionalization (slope and intercept)
// [[Rcpp::export]]
arma::mat calCor(
    arma::vec Y, arma::vec X, arma::mat I, arma::uvec si,
    arma::mat dist, arma::vec init, double sd,
    double na, double nb, double ga, double gb,
    int n, int N, int nBurnin, int maxIter, int nThin, int nReport) {
  
  // PRIORS //
  double Umin = 3 / dist.max();
  double Umax = 3 / (0.1 * dist.max());
  double na1 = 1;
  double nb1 = 0.01;
  
  // PARAMETERS //
  arma::mat keep(maxIter / nThin, 8 + 2 * n);
  arma::vec parameters(8 + 2 * n); 
  // para: beta0 a11 phiv0 beta1 a21 a22 phiv1 sigma v0s v1s
  // init: a0's a1's mean0 mean1 prec0 prec1 decay0 decay1 precT2 a21 
  parameters(0) = init(2 * n);
  parameters(1) = init(2 * n + 2);
  parameters(2) = init(2 * n + 4);
  parameters(3) = init(2 * n + 1);
  parameters(4) = init(2 * n + 7);
  parameters(5) = init(2 * n + 3);
  parameters(6) = init(2 * n + 5);
  parameters(7) = init(2 * n + 6);
  parameters(arma::span(8, 7 + n)) = init(arma::span(0, n - 1));
  parameters(arma::span(8 + n, 7 + 2 * n)) = init(arma::span(n, 2 * n - 1));
  
  arma::vec v0 = parameters(arma::span(8, 7 + n));
  arma::vec v1 = parameters(arma::span(8 + n, 7 + 2 * n));
  arma::vec Z = Y - (parameters(0) + parameters(1) * v0.elem(si) + (parameters(3) + parameters(4) * v0.elem(si) + parameters(5) * v1.elem(si)) % X);
  
  // HELP //
  arma::vec onen(n, arma::fill::ones);
  arma::vec oneN(N, arma::fill::ones);
  
  arma::mat Delta;
  arma::mat Chi;
  double delta;
  double chi;
  arma::vec process;
  arma::mat r;
  double A;
  
  arma::vec aux(N);
  
  double phiAux;
  arma::mat RAux(n, n);
  double RlogdetAux;
  
  arma::mat Rb0    = arma::inv_sympd(exp(- parameters(2) * dist));
  double Rb0logdet = arma::log_det_sympd(Rb0);
  
  arma::mat Rb1    = arma::inv_sympd(exp(- parameters(6) * dist));
  double Rb1logdet = arma::log_det_sympd(Rb1);
  
  double sumita;
  double sumitaAux;
  
  // ITERATIONS //
  for (int b = 1 - nBurnin; b <= maxIter; ++b) {
    
    // beta0
    Z += parameters(0);
    delta = 1 / (N * parameters(7) + nb);
    chi   = arma::accu(Z) * parameters(7) + na * nb;
    parameters(0) = R::rnorm(delta * chi, sqrt(delta));
    
    // a11
    aux = v0.elem(si);
    Z += parameters(1) * aux - parameters(0);
    delta = 1 / (arma::accu(pow(aux, 2)) * parameters(7) + nb1);
    chi   = (Z.t() * aux).eval()(0) * parameters(7) + na * nb1;
    parameters(1) = rfnorm(delta * chi, sqrt(delta));
    
    // phi v0
    process = parameters(arma::span(8, 7 + n));
    phiAux     = rtnorm(parameters(2), sd, Umin, Umax);
    RAux       = arma::inv_sympd(exp(- phiAux * dist));
    RlogdetAux = arma::log_det_sympd(RAux);
    sumita    = (process.t() * Rb0 * process).eval()(0,0);  
    sumitaAux = (process.t() * RAux * process).eval()(0,0);  
    A = (RlogdetAux - sumitaAux) / 2 - 
      (Rb0logdet - sumita) / 2 + 
      dtnorm(parameters(2), phiAux, sd, Umin, Umax) - 
      dtnorm(phiAux, parameters(2), sd, Umin, Umax);
    if (log(R::runif(0, 1)) <= A) {
      parameters(2) = phiAux;
      Rb0 = RAux;
      Rb0logdet = RlogdetAux;
      sumita = sumitaAux;
    }

    // beta1
    Z += parameters(3) * X - parameters(1) * aux;
    delta = 1 / (arma::accu(pow(X, 2)) * parameters(7) + nb);
    chi   = (Z.t() * X).eval()(0) * parameters(7) + na1 * nb;
    parameters(3) = R::rnorm(delta * chi, sqrt(delta));

    // a21
    aux = v0.elem(si) % X;
    Z += parameters(4) * aux - parameters(3) * X;
    delta = 1 / (arma::accu(pow(aux, 2)) * parameters(7) + nb);
    chi   = (Z.t() * aux).eval()(0) * parameters(7) + na * nb;
    parameters(4) = R::rnorm(delta * chi, sqrt(delta));
    Z -= parameters(4) * aux;
    
    // a22
    aux = v1.elem(si) % X;
    Z += parameters(5) * aux; 
    delta = 1 / (arma::accu(pow(aux, 2)) * parameters(7) + nb1);
    chi   = (Z.t() * aux).eval()(0) * parameters(7) + na * nb1;
    parameters(5) = rfnorm(delta * chi, sqrt(delta));
    Z -= parameters(5) * aux;
    
    // phi v1
    process = parameters(arma::span(8 + n, 7 + 2 * n));
    phiAux     = rtnorm(parameters(6), sd, Umin, Umax);
    RAux       = arma::inv_sympd(exp(- phiAux * dist));
    RlogdetAux = arma::log_det_sympd(RAux);
    sumita    = (process.t() * Rb1 * process).eval()(0,0);  
    sumitaAux = (process.t() * RAux * process).eval()(0,0);  
    A = (RlogdetAux - sumitaAux) / 2 - 
      (Rb1logdet - sumita) / 2 + 
      dtnorm(parameters(6), phiAux, sd, Umin, Umax) - 
      dtnorm(phiAux, parameters(6), sd, Umin, Umax);
    if (log(R::runif(0, 1)) <= A) {
      parameters(6) = phiAux;
      Rb1 = RAux;
      Rb1logdet = RlogdetAux;
      sumita = sumitaAux;
    }
    
    // sigma
    parameters(7) = R::rgamma(N / 2 + ga, 
               1 / (arma::accu(arma::pow(Z, 2)) / 2 + gb));
    
    // v0 si i=1,...,n
    aux = parameters(1) + parameters(4) * X;
    Z += aux % v0.elem(si);
    for (int iInd = 0; iInd < n; ++iInd){
      process = v0;
      process.shed_row(iInd);
      r = Rb0.row(iInd);
      r.shed_col(iInd);
      delta = 1 / (arma::accu(pow(aux % I.col(iInd), 2)) * parameters(7) + Rb0(iInd, iInd));
      chi   = (Z.t() * (aux % I.col(iInd))).eval()(0) * parameters(7) +
        (- r * process).eval()(0,0);
      v0(iInd) = R::rnorm(delta * chi, sqrt(delta));
    }
    parameters(arma::span(8, 7 + n)) = v0;
    Z -= aux % v0.elem(si);
    
    // v1 si i=1,...,n
    aux = parameters(5) * X;
    Z += aux % v1.elem(si);
    for (int iInd = 0; iInd < n; ++iInd){
      process = v1;
      process.shed_row(iInd);
      r = Rb1.row(iInd);
      r.shed_col(iInd);
      delta = 1 / (arma::accu(pow(aux % I.col(iInd), 2)) * parameters(7) + Rb1(iInd, iInd));
      chi   = (Z.t() * (aux % I.col(iInd))).eval()(0) * parameters(7) +
        (- r * process).eval()(0,0);
      v1(iInd) = R::rnorm(delta * chi, sqrt(delta));
    }
    parameters(arma::span(8 + n, 7 + 2 * n)) = v1;  
    Z -= aux % v1.elem(si);
    
    //
    if ((b != 0) && (b % nReport == 0)) {
      Rcpp::Rcout << "The value of b : " << b << "\n";
    }
    
    if ((b > 0) && (b % nThin == 0)) {
      init(2 * n) = parameters(0);
      init(2 * n + 2) = parameters(1);
      init(2 * n + 4) = parameters(2);
      init(2 * n + 1) = parameters(3);
      init(2 * n + 7) = parameters(4);
      init(2 * n + 3) = parameters(5);
      init(2 * n + 5) = parameters(6);
      init(2 * n + 6) = parameters(7);
      init(arma::span(0, n - 1)) = parameters(arma::span(8, 7 + n));
      init(arma::span(n, 2 * n - 1)) = parameters(arma::span(8 + n, 7 + 2 * n));
      keep.row(b / nThin - 1) = init.t();
    }
  }
  
  return keep;
}

// AR(1) (solo permito 1 iter)
// [[Rcpp::export]]
arma::vec AR1(
    arma::vec Y, arma::vec init,
    double na, double nb, double ga, double gb, 
    int T) {
  
  // PARAMETERS //
  arma::vec keep(3); // beta rho prec
  double beta = init(0);
  double rho  = init(1);
  double prec = init(2);
  
  arma::vec Z(T);
  Z(0) = Y(0) - beta;
  for (int t = 1; t < T; ++t) {
    Z(t) = Y(t) - beta - rho * (Y(t - 1) - beta);
  }
  
  double sumita;
  double delta;
  double chi;
  arma::vec innovation(T - 1);
  
  // ITERATIONS //
  // prec
  prec = R::rgamma(T / 2 + ga, 
    1 / ((Z.t() * Z).eval()(0) / 2 + gb));
  
  
  // beta
  Z(0) += beta;
  Z(arma::span(1, T - 1)) += (1 - rho) * beta;
  delta = 1 / ((pow(1 - rho, 2) * (T - 1) + 1) * prec + nb);
  chi   = ((1 - rho) * arma::accu(Z(arma::span(1, T - 1))) + Z(0)) * prec +  
    na * nb;
  beta = R::rnorm(delta * chi, sqrt(delta));
  Z(0) -= beta;
  Z(arma::span(1, T - 1)) -= (1 - rho) * beta;
  
  // rho
  innovation = Y(arma::span(0, T - 2)) - beta;
  sumita = arma::accu(pow(innovation, 2));
  Z(arma::span(1, T - 1)) += rho * innovation;
  delta = 1 / (sumita * prec);
  chi   = arma::accu(innovation % Z(arma::span(1, T - 1))) * prec;
  rho = rtnorm(delta * chi, sqrt(delta), -1, 1);
  //Z(arma::span(1, T - 1)) -= rho * innovation;
  
  // SAVE //
  keep(0) = beta;
  keep(1) = rho;
  keep(2) = prec;

  return keep;
}

// bivariate AR(1) (solo permito 1 iter) (validada, funciona bien)
// [[Rcpp::export]]
arma::vec biAR1(
    arma::vec Y0, arma::vec Y1, 
    arma::vec init,
    double na1, double na, double nb, int T) {
  
  // PARAMETERS //
  arma::vec keep(7); // beta0 beta1 rho0 rho1 v11 v22 v12
  arma::vec beta = init(arma::span(0, 1));
  arma::vec rho  = init(arma::span(2, 3));
  arma::mat V(2,2, arma::fill::zeros);
  V(0,0) = init(4);
  V(1,1) = init(5);
  V(0,1) = init(6);
  V(1,0) = init(6);
  
  arma::mat Z(T, 2);
  Z(0, 0) = Y0(0) - beta(0);
  Z(0, 1) = Y1(0) - beta(1);
  for (int t = 1; t < T; ++t) {
    Z(t, 0) = Y0(t) - beta(0) - rho(0) * (Y0(t - 1) - beta(0));
    Z(t, 1) = Y1(t) - beta(1) - rho(1) * (Y1(t - 1) - beta(1));
  }
  
  double sumita;
  double delta;
  double chi;
  arma::vec innovation(T - 1);
  
  // ITERATIONS //
  // V
  arma::mat W = arma::eye(2, 2);
  for (int t = 0; t < T; ++t) {
    W += Z.row(t).t() * Z.row(t);
  }
  V = rWishart(arma::inv_sympd(W), 2, T + 2);  
  
  // beta0
  Z(0, 0) += beta(0);
  Z(arma::span(1, T - 1), 0) += (1 - rho(0)) * beta(0);
  delta = 1 / ((pow(1 - rho(0), 2) * (T - 1) + 1) * V(0,0) + nb);
  chi   = ((1 - rho(0)) * arma::accu(Z(arma::span(1, T - 1), 0)) + Z(0, 0)) * V(0,0) + 
    ((1 - rho(0)) * arma::accu(Z(arma::span(1, T - 1), 1)) + Z(0, 1)) * V(0,1) +  
    na * nb;
  beta(0) = R::rnorm(delta * chi, sqrt(delta));
  Z(0, 0) -= beta(0);
  Z(arma::span(1, T - 1), 0) -= (1 - rho(0)) * beta(0);
  
  // beta1
  Z(0, 1) += beta(1);
  Z(arma::span(1, T - 1), 1) += (1 - rho(1)) * beta(1);
  delta = 1 / ((pow(1 - rho(1), 2) * (T - 1) + 1) * V(1,1) + nb);
  chi   = ((1 - rho(1)) * arma::accu(Z(arma::span(1, T - 1), 1)) + Z(0, 1)) * V(1,1) + 
    ((1 - rho(1)) * arma::accu(Z(arma::span(1, T - 1), 0)) + Z(0, 0)) * V(0,1) +  
    na1 * nb;
  beta(1) = R::rnorm(delta * chi, sqrt(delta));
  Z(0, 1) -= beta(1);
  Z(arma::span(1, T - 1), 1) -= (1 - rho(1)) * beta(1);
  
  // rho0
  innovation = Y0(arma::span(0, T - 2)) - beta(0);
  sumita = arma::accu(pow(innovation, 2));
  Z(arma::span(1, T - 1), 0) += rho(0) * innovation;
  delta = 1 / (sumita * V(0,0));
  chi   = arma::accu(innovation % (Z(arma::span(1, T - 1), 0) * V(0,0) + Z(arma::span(1, T - 1), 1) * V(0,1)));
  rho(0) = rtnorm(delta * chi, sqrt(delta), -1, 1);
  Z(arma::span(1, T - 1), 0) -= rho(0) * innovation;
  
  // rho1
  innovation = Y1(arma::span(0, T - 2)) - beta(1);
  sumita = arma::accu(pow(innovation, 2));
  Z(arma::span(1, T - 1), 1) += rho(1) * innovation;
  delta = 1 / (sumita * V(1,1));
  chi   = arma::accu(innovation % (Z(arma::span(1, T - 1), 1) * V(1,1) + Z(arma::span(1, T - 1), 0) * V(0,1)));
  rho(1) = rtnorm(delta * chi, sqrt(delta), -1, 1);
  //Z(arma::span(1, T - 1), 1) -= rho(1) * innovation;
  
  // SAVE //
  keep(arma::span(0, 1)) = beta;
  keep(arma::span(2, 3)) = rho;
  keep(4) = V(0,0);
  keep(5) = V(1,1);
  keep(6) = V(0,1);
  
  return keep;
}

// Calibration Coregionalization (slope and intercept)
// Time varying coefficients
// [[Rcpp::export]]
arma::mat calCorTime(
    arma::vec Y, arma::vec X, 
    arma::mat Xt, 
    arma::mat I, arma::uvec si,
    arma::mat dist, arma::vec init, double sd,
    double na, double nb, double ga, double gb,
    int n, int T, int N, 
    int nBurnin, int maxIter, int nThin, int nReport) {
  
  // PRIORS //
  double Umin = 3 / dist.max();
  double Umax = 3 / (0.1 * dist.max());
  double na1 = 1;
  double nb1 = 0.01;
  
  // PARAMETERS //
  arma::mat keep(maxIter / nThin, 13 + 2 * (n + T));
  arma::vec parameters(8 + 2 * n); 
  arma::vec dynamic = init(arma::span(8 + 2 * n, 12 + 2 * (n + T))); // a0t' a1t' rho0 rho1 v00 v11 v01
  // para: beta0 a11 phiv0 beta1 a21 a22 phiv1 sigma v0s v1s
  // init: a0's a1's mean0 mean1 prec0 prec1 decay0 decay1 precT2 a21 (a0t' a1t' rho0 rho1 v00 v11 v01)
  parameters(0) = init(2 * n);
  parameters(1) = init(2 * n + 2);
  parameters(2) = init(2 * n + 4);
  parameters(3) = init(2 * n + 1);
  parameters(4) = init(2 * n + 7);
  parameters(5) = init(2 * n + 3);
  parameters(6) = init(2 * n + 5);
  parameters(7) = init(2 * n + 6);
  parameters(arma::span(8, 7 + n)) = init(arma::span(0, n - 1));
  parameters(arma::span(8 + n, 7 + 2 * n)) = init(arma::span(n, 2 * n - 1));
  arma::vec v0 = parameters(arma::span(8, 7 + n));
  arma::vec v1 = parameters(arma::span(8 + n, 7 + 2 * n));
  arma::vec Z = Y - (Xt * dynamic(arma::span(0, T - 1)) + parameters(1) * v0.elem(si) + (Xt * dynamic(arma::span(T, 2 * T - 1)) + parameters(4) * v0.elem(si) + parameters(5) * v1.elem(si)) % X);
  
  // HELP //
  arma::vec onen(n, arma::fill::ones);
  arma::vec oneN(N, arma::fill::ones);
  
  arma::mat Delta;
  arma::mat Chi;
  double delta;
  double chi;
  arma::vec process;
  arma::mat r;
  double A;
  
  arma::vec aux(N);
  
  double phiAux;
  arma::mat RAux(n, n);
  double RlogdetAux;
  
  arma::mat Rb0    = arma::inv_sympd(exp(- parameters(2) * dist));
  double Rb0logdet = arma::log_det_sympd(Rb0);
  
  arma::mat Rb1    = arma::inv_sympd(exp(- parameters(6) * dist));
  double Rb1logdet = arma::log_det_sympd(Rb1);
  
  double sumita;
  double sumitaAux;
  
  double precb0t;
  double precb1t;
  arma::vec tildeb0t(T);
  arma::vec tildeb1t(T);
  // beta0 beta1 rho0 rho1 v00 v11 v01
  arma::vec auxSave(7);
  auxSave(0) = parameters(0);
  auxSave(1) = parameters(3);
  auxSave(arma::span(2, 6)) = dynamic(arma::span(2 * T, 2 * T + 4));
  
  // ITERATIONS //
  for (int b = 1 - nBurnin; b <= maxIter; ++b) {
    
    tildeb0t = dynamic(arma::span(0, T - 1)) - parameters(0);
    tildeb1t = dynamic(arma::span(T, 2 * T - 1)) - parameters(3);
    
    // beta0t
    precb0t = (1 + pow(dynamic(2 * T), 2) * dynamic(2 * T + 2));
    Z += Xt * dynamic(arma::span(0, T - 1));
    delta = 1 / (arma::accu(Xt.col(0)) * parameters(7) + precb0t);
    chi   = arma::accu(Xt.col(0) % Z) * parameters(7) + 
      parameters(0) * precb0t + 
      dynamic(2 * T + 2) * dynamic(2 * T) * tildeb0t(1) +
      dynamic(2 * T + 4) * (dynamic(2 * T) * tildeb1t(1) - (1 + dynamic(2 * T) * dynamic(2 * T + 1)) * tildeb1t(0));
    dynamic(0) = R::rnorm(delta * chi, sqrt(delta));
    tildeb0t(0) = dynamic(0) - parameters(0);
    for (int tInd = 1; tInd < T - 1; ++tInd) {
      delta = 1 / (arma::accu(Xt.col(tInd)) * parameters(7) + precb0t);
      chi   = arma::accu(Xt.col(tInd) % Z) * parameters(7) + 
        parameters(0) * precb0t + 
        dynamic(2 * T + 2) * dynamic(2 * T) * (tildeb0t(tInd - 1) + tildeb0t(tInd + 1)) +
        dynamic(2 * T + 4) * (dynamic(2 * T + 1) * tildeb1t(tInd - 1) - (1 + dynamic(2 * T) * dynamic(2 * T + 1)) * tildeb1t(tInd) + dynamic(2 * T) * tildeb1t(tInd + 1));
      dynamic(tInd) = R::rnorm(delta * chi, sqrt(delta));
      tildeb0t(tInd) = dynamic(tInd) - parameters(0);
    }
    delta = 1 / (arma::accu(Xt.col(T - 1)) * parameters(7) + dynamic(2 * T + 2));
    chi   = arma::accu(Xt.col(T - 1) % Z) * parameters(7) + 
      parameters(0) * dynamic(2 * T + 2) + 
      dynamic(2 * T + 2) * dynamic(2 * T) * tildeb0t(T - 2) +
      dynamic(2 * T + 4) * (dynamic(2 * T + 1) * tildeb1t(T - 2) - tildeb1t(T - 1));
    dynamic(T - 1) = R::rnorm(delta * chi, sqrt(delta));
    tildeb0t(T - 1) = dynamic(T - 1) - parameters(0);
    
    // a11
    aux = v0.elem(si);
    Z += parameters(1) * aux - Xt * dynamic(arma::span(0, T - 1));
    delta = 1 / (arma::accu(pow(aux, 2)) * parameters(7) + nb1);
    chi   = (Z.t() * aux).eval()(0) * parameters(7) + na * nb1;
    parameters(1) = rfnorm(delta * chi, sqrt(delta));
    
    // phi v0
    process = parameters(arma::span(8, 7 + n));
    phiAux     = rtnorm(parameters(2), sd, Umin, Umax);
    RAux       = arma::inv_sympd(exp(- phiAux * dist));
    RlogdetAux = arma::log_det_sympd(RAux);
    sumita    = (process.t() * Rb0 * process).eval()(0,0);  
    sumitaAux = (process.t() * RAux * process).eval()(0,0);  
    A = (RlogdetAux - sumitaAux) / 2 - 
      (Rb0logdet - sumita) / 2 + 
      dtnorm(parameters(2), phiAux, sd, Umin, Umax) - 
      dtnorm(phiAux, parameters(2), sd, Umin, Umax);
    if (log(R::runif(0, 1)) <= A) {
      parameters(2) = phiAux;
      Rb0 = RAux;
      Rb0logdet = RlogdetAux;
      sumita = sumitaAux;
    }
    
    // beta1t // dynamic: a0t' a1t' rho0 rho1 v00 v11 v01
    precb1t = (1 + pow(dynamic(2 * T + 1), 2) * dynamic(2 * T + 3));
    Z += (Xt * dynamic(arma::span(T, 2 * T - 1))) % X - parameters(1) * aux;
    delta = 1 / (arma::accu(pow(X % Xt.col(0), 2)) * parameters(7) + precb1t);
    chi   = arma::accu(Xt.col(0) % Z % X) * parameters(7) + 
      parameters(3) * precb1t + 
      dynamic(2 * T + 3) * dynamic(2 * T + 1) * tildeb1t(1) +
      dynamic(2 * T + 4) * (dynamic(2 * T + 1) * tildeb0t(1) - (1 + dynamic(2 * T) * dynamic(2 * T + 1)) * tildeb0t(0));
    dynamic(T) = R::rnorm(delta * chi, sqrt(delta));
    tildeb1t(0) = dynamic(T) - parameters(3);
    for (int tInd = 1; tInd < T - 1; ++tInd) {
      delta = 1 / (arma::accu(pow(X % Xt.col(tInd), 2)) * parameters(7) + precb1t);
      chi   = arma::accu(Xt.col(tInd) % Z % X) * parameters(7) + 
        parameters(3) * precb1t + 
        dynamic(2 * T + 3) * dynamic(2 * T + 1) * (tildeb1t(tInd - 1) + tildeb1t(tInd + 1)) +
        dynamic(2 * T + 4) * (dynamic(2 * T) * tildeb0t(tInd - 1) - (1 + dynamic(2 * T) * dynamic(2 * T + 1)) * tildeb0t(tInd) + dynamic(2 * T + 1) * tildeb0t(tInd + 1));
      dynamic(tInd + T) = R::rnorm(delta * chi, sqrt(delta));
      tildeb1t(tInd) = dynamic(tInd + T) - parameters(3);
    }
    delta = 1 / (arma::accu(pow(X % Xt.col(T - 1), 2)) * parameters(7) + dynamic(2 * T + 3));
    chi   = arma::accu(Xt.col(T - 1) % Z % X) * parameters(7) + 
      parameters(3) * dynamic(2 * T + 3) + 
      dynamic(2 * T + 3) * dynamic(2 * T + 1) * tildeb1t(T - 2) +
      dynamic(2 * T + 4) * (dynamic(2 * T) * tildeb0t(T - 2) - tildeb0t(T - 1));
    dynamic(2 * T - 1) = R::rnorm(delta * chi, sqrt(delta));
    tildeb1t(T - 1) = dynamic(2 * T - 1) - parameters(3);
    
    // a21
    aux = v0.elem(si) % X;
    Z += parameters(4) * aux - (Xt * dynamic(arma::span(T, 2 * T - 1))) % X;
    delta = 1 / (arma::accu(pow(aux, 2)) * parameters(7) + nb);
    chi   = (Z.t() * aux).eval()(0) * parameters(7) + na * nb;
    parameters(4) = R::rnorm(delta * chi, sqrt(delta));
    Z -= parameters(4) * aux;
    
    // a22
    aux = v1.elem(si) % X;
    Z += parameters(5) * aux; 
    delta = 1 / (arma::accu(pow(aux, 2)) * parameters(7) + nb1);
    chi   = (Z.t() * aux).eval()(0) * parameters(7) + na * nb1;
    parameters(5) = rfnorm(delta * chi, sqrt(delta));
    Z -= parameters(5) * aux;
    
    // phi v1
    process = parameters(arma::span(8 + n, 7 + 2 * n));
    phiAux     = rtnorm(parameters(6), sd, Umin, Umax);
    RAux       = arma::inv_sympd(exp(- phiAux * dist));
    RlogdetAux = arma::log_det_sympd(RAux);
    sumita    = (process.t() * Rb1 * process).eval()(0,0);  
    sumitaAux = (process.t() * RAux * process).eval()(0,0);  
    A = (RlogdetAux - sumitaAux) / 2 - 
      (Rb1logdet - sumita) / 2 + 
      dtnorm(parameters(6), phiAux, sd, Umin, Umax) - 
      dtnorm(phiAux, parameters(6), sd, Umin, Umax);
    if (log(R::runif(0, 1)) <= A) {
      parameters(6) = phiAux;
      Rb1 = RAux;
      Rb1logdet = RlogdetAux;
      sumita = sumitaAux;
    }
    
    // sigma
    parameters(7) = R::rgamma(N / 2 + ga, 
               1 / (arma::accu(arma::pow(Z, 2)) / 2 + gb));
    
    // v0 si i=1,...,n
    aux = parameters(1) + parameters(4) * X;
    Z += aux % v0.elem(si);
    for (int iInd = 0; iInd < n; ++iInd){
      process = v0;
      process.shed_row(iInd);
      r = Rb0.row(iInd);
      r.shed_col(iInd);
      delta = 1 / (arma::accu(pow(aux % I.col(iInd), 2)) * parameters(7) + Rb0(iInd, iInd));
      chi   = (Z.t() * (aux % I.col(iInd))).eval()(0) * parameters(7) +
        (- r * process).eval()(0,0);
      v0(iInd) = R::rnorm(delta * chi, sqrt(delta));
    }
    parameters(arma::span(8, 7 + n)) = v0;
    Z -= aux % v0.elem(si);
    
    // v1 si i=1,...,n
    aux = parameters(5) * X;
    Z += aux % v1.elem(si);
    for (int iInd = 0; iInd < n; ++iInd){
      process = v1;
      process.shed_row(iInd);
      r = Rb1.row(iInd);
      r.shed_col(iInd);
      delta = 1 / (arma::accu(pow(aux % I.col(iInd), 2)) * parameters(7) + Rb1(iInd, iInd));
      chi   = (Z.t() * (aux % I.col(iInd))).eval()(0) * parameters(7) +
        (- r * process).eval()(0,0);
      v1(iInd) = R::rnorm(delta * chi, sqrt(delta));
    }
    parameters(arma::span(8 + n, 7 + 2 * n)) = v1;  
    Z -= aux % v1.elem(si);
    
    // beta0 beta1 rho0 rho1 v00 v11 v01
    auxSave = biAR1(
      dynamic(arma::span(0, T - 1)), 
      dynamic(arma::span(T, 2 * T - 1)), 
      auxSave,
      na1, na, nb, T);
    parameters(0) = auxSave(0);
    parameters(3) = auxSave(1);
    dynamic(arma::span(2 * T, 2 * T + 4)) = auxSave(arma::span(2, 6));
    
    //
    if ((b != 0) && (b % nReport == 0)) {
      Rcpp::Rcout << "The value of b : " << b << "\n";
    }
    
    if ((b > 0) && (b % nThin == 0)) {
      init(2 * n) = parameters(0);
      init(2 * n + 2) = parameters(1);
      init(2 * n + 4) = parameters(2);
      init(2 * n + 1) = parameters(3);
      init(2 * n + 7) = parameters(4);
      init(2 * n + 3) = parameters(5);
      init(2 * n + 5) = parameters(6);
      init(2 * n + 6) = parameters(7);
      init(arma::span(0, n - 1)) = parameters(arma::span(8, 7 + n));
      init(arma::span(n, 2 * n - 1)) = parameters(arma::span(8 + n, 7 + 2 * n));
      init(arma::span(2 * n + 8, 2 * (n + T) + 12)) = dynamic;
      // keep: a0's a1's mean0 mean1 prec0 prec1 decay0 decay1 precT2 a21 (a0t' a1t' rho0 rho1 v00 v11 v01)
      keep.row(b / nThin - 1) = init.t();
    }
  }
  
  return keep;
}


// Both calibrations with coregionalization
// Time varying coefs for ZOOP calibration
// [[Rcpp::export]]
arma::mat GibbsZoop4Cpp(
    arma::vec Y1, arma::vec Y2, arma::vec T1, arma::vec T2,
    arma::mat X, arma::mat Xeta, arma::mat XTEMP, 
    arma::vec Rsum, double RsumTEMP, arma::cube R, arma::mat RTEMP, double dMax,
    arma::mat dY2, arma::mat dT2, double sd,
    arma::vec oneY1, arma::vec oneY2, arma::vec oneT1, arma::vec oneT2,
    arma::uvec indY1, arma::uvec indY2, arma::uvec indT1, arma::uvec indT2,
    arma::mat IY2, arma::mat IT2, arma::mat ITEMP, arma::uvec sY2, arma::uvec sT2, arma::uvec sTEMP, arma::mat Xt,
    int N, int NY1, int NY2, int NT1, int NT2, int T, arma::vec nt, double Ent, int nT2, int nY2, int nTEMP, int p,
    double na, double nb, double ga, double gb,
    int nBurnin, int maxIter, int nThin, int nReport) {
  
  // PARAMETERS //
  arma::vec temp(N + nTEMP + T + 7, arma::fill::ones); // temp'tls, mean's, mean, precmean, precyear, sin, cos, prec, year1 ... yearT, rhoyear
  arma::vec alpha(2 * nT2 + 8, arma::fill::ones); // a0's a1's mean0 mean1 prec0 prec1 decay0 decay1 precT2 a21
  arma::vec lambda(2 * nY2 + 2 * T + 13, arma::fill::ones); // l0's l1's mean0 mean1 prec0 prec1 decay0 decay1 precY2 a21 l0't l1't rho0 rho1 v00 v11 v01
  arma::vec parameters(Ent + T + p + 6, arma::fill::ones); // eta'ts eta't X' beta1 beta0 precY1 precT1 precEtas precEtat
  arma::mat keep(maxIter / nThin, (Ent + T + p + 6) + (2 * nT2 + 8) + (2 * nY2 + 2 * T + 13) + (N + nTEMP + T + 7));
  
  alpha(2 * nT2 + 4) = 10 / dMax; // decay0
  alpha(2 * nT2 + 5) = 10 / dMax; // decay1
  lambda(2 * nY2 + 4) = 10 / dMax; // decay0
  lambda(2 * nY2 + 5) = 10 / dMax; // decay1
  
  //parameters(Ent + T + p + 3) = 10000; // precT1
  temp(N + nTEMP + 5) = INFINITY; // precTemp
  temp(N + nTEMP + T + 6) = 0; // rhoyear
  lambda(2 * nY2 + 2 * T + 8) = 0; //rho0
  lambda(2 * nY2 + 2 * T + 9) = 0; //rho1
  
  arma::vec TRUEtemp = ITEMP * temp(arma::span(N, N + nTEMP - 1)) + XTEMP * temp(arma::span(N + nTEMP + 3, N + nTEMP + 4)) + Xt * temp(arma::span(N + nTEMP + 6, N + nTEMP + T + 5));
  arma::vec ZOOP = Xeta * parameters(arma::span(0, Ent - 1)) + X * parameters(arma::span(Ent + T, Ent + T + p - 1)) + parameters(Ent + T + p) * TRUEtemp;
  
  // HELP //
  arma::mat Delta;
  arma::mat Chi;
  double delta;
  double chi;
  
  arma::vec alphas(nT2);
  arma::vec lambdas(nY2);
  
  arma::vec auxAlpha0(N, arma::fill::zeros); 
  alphas = alpha(2 * nT2) + alpha(2 * nT2 + 2) * alpha(arma::span(0, nT2 - 1));
  auxAlpha0.elem(indT2) = alphas.elem(sT2);
  arma::vec auxAlpha1(N, arma::fill::zeros); 
  alphas = alpha(2 * nT2 + 1) + alpha(2 * nT2 + 7) * alpha(arma::span(0, nT2 - 1)) + alpha(2 * nT2 + 3) * alpha(arma::span(nT2, 2 * nT2 - 1));
  auxAlpha1.elem(indT2) = alphas.elem(sT2);
  
  arma::vec auxLambda0(N, arma::fill::zeros); 
  lambdas = lambda(2 * nY2 + 2) * lambda(arma::span(0, nY2 - 1));
  auxLambda0.elem(indY2) = Xt.rows(indY2) * lambda(arma::span(2 * nY2 + 8, 2 * nY2 + T + 7)) + lambdas.elem(sY2);
  arma::vec auxLambda1(N, arma::fill::zeros); 
  lambdas = lambda(2 * nY2 + 7) * lambda(arma::span(0, nY2 - 1)) + lambda(2 * nY2 + 3) * lambda(arma::span(nY2, 2 * nY2 - 1));
  auxLambda1.elem(indY2) = Xt.rows(indY2) * lambda(arma::span(2 * nY2 + T + 8, 2 * nY2 + 2 * T + 7)) + lambdas.elem(sY2);
  
  arma::vec oneN(N, arma::fill::ones);
  arma::vec auxN(N);
  arma::vec auxT(T);
  arma::vec auxNY1(NY1);
  arma::vec auxNY2(NY2);
  arma::vec auxNT1(NT1);
  arma::vec auxNT2(NT2);
  
  arma::mat r;
  arma::vec process;
  double sumita;
  double contador;
  double rho2;
  
  arma::uvec col0 = { 0 };
  arma::uvec col1 = { 1 };
  
  double sumY1Sin2 = arma::accu(pow(XTEMP(indY1, col0), 2));
  double sumY1Cos2 = arma::accu(pow(XTEMP(indY1, col1), 2));
  double sumT1Sin2 = arma::accu(pow(XTEMP(indT1, col0), 2));
  double sumT1Cos2 = arma::accu(pow(XTEMP(indT1, col1), 2));
  
  arma::mat XtauxY1 = Xt.rows(indY1);
  arma::mat XtauxY2 = Xt.rows(indY2);
  arma::mat XtauxT1 = Xt.rows(indT1);
  arma::mat XtauxT2 = Xt.rows(indT2);
  
  arma::vec NY1t = arma::sum(XtauxY1).t();
  arma::vec NT1t = arma::sum(XtauxT1).t();
  
  arma::mat XsauxY1 = ITEMP.rows(indY1);
  arma::mat XsauxY2 = ITEMP.rows(indY2);
  arma::mat XsauxT1 = ITEMP.rows(indT1);
  arma::mat XsauxT2 = ITEMP.rows(indT2);
  
  arma::vec NY1s = arma::sum(XsauxY1).t();
  arma::vec NY2s = arma::sum(XsauxY2).t();
  arma::vec NT1s = arma::sum(XsauxT1).t();
  arma::vec NT2s = arma::sum(XsauxT2).t();
  
  arma::vec initAR1(3);
  initAR1(0) = temp(N + nTEMP);
  initAR1(1) = temp(N + nTEMP + T + 6);
  initAR1(2) = temp(N + nTEMP + 2);
  
  arma::vec CUMSUMnt(T + 1);
  CUMSUMnt(0) = 0;
  CUMSUMnt(arma::span(1, T)) = arma::cumsum(nt);
  
  // ITERATIONS //
  for (int b = 1 - nBurnin; b <= maxIter; ++b) {
    
    //// temp t l s OK
    //ZOOP -= parameters(Ent + T + p) * temp(arma::span(0, N - 1));
    //for (int j = 0; j < N; ++j) {
    //  delta = 1 / (oneT1(j) * parameters(Ent + T + p + 3) +
    //    oneT2(j) * pow(auxAlpha1(j), 2) * alpha(2 * nT2 + 6) + 
    //    oneY1(j) * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
    //    oneY2(j) * pow(auxLambda1(j) * parameters(Ent + T + p), 2) * lambda(2 * nY2 + 6) + 
    //    temp(N + nTEMP + 5));
    //  chi   = T1(j) * oneT1(j) * parameters(Ent + T + p + 3) +
    //    (T2(j) - auxAlpha0(j)) * oneT2(j) * auxAlpha1(j) * alpha(2 * nT2 + 6) + 
    //    (Y1(j) - ZOOP(j)) * oneY1(j) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
    //    (Y2(j) - (auxLambda0(j) + auxLambda1(j) * ZOOP(j))) * oneY2(j) * auxLambda1(j) * parameters(Ent + T + p) * lambda(2 * nY2 + 6) +
    //    ((ITEMP.row(j) * temp(arma::span(N, N + nTEMP - 1))).eval()(0) + temp(N + nTEMP + 3) * XTEMP(j,0) + temp(N + nTEMP + 4) * XTEMP(j,1) + (Xt.row(j) * temp(arma::span(N + nTEMP + 6, N + nTEMP + T + 5))).eval()(0)) * temp(N + nTEMP + 5);
    //  temp(j) = R::rnorm(chi * delta, sqrt(delta));
    //}
    //ZOOP += parameters(Ent + T + p) * temp(arma::span(0, N - 1));
    //
    //temp(arma::span(N, N + nTEMP + T + 6)) = splmMR(
    //  temp(arma::span(0, N - 1)), XTEMP,
    //  ITEMP, sTEMP, RTEMP, RsumTEMP,
    //  temp(arma::span(N, N + nTEMP + T + 6)), 
    //  na, nb, ga, gb,
    //  2, nTEMP, N, T, Xt,
    //  0, 1, 1, 10
    //).t();
    
    // sin (temp)
    TRUEtemp -= temp(N + nTEMP + 3) * XTEMP.col(0);
    ZOOP     -= parameters(Ent + T + p) * temp(N + nTEMP + 3) * XTEMP.col(0);
    delta = 1 / (
      sumT1Sin2 * parameters(Ent + T + p + 3) +
        arma::accu(pow(XTEMP(indT2, col0) % auxAlpha1(indT2), 2)) * alpha(2 * nT2 + 6) + 
        sumY1Sin2 * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
        arma::accu(pow(XTEMP(indY2, col0) % auxLambda1(indY2) * parameters(Ent + T + p), 2)) * lambda(2 * nY2 + 6) + 
        nb);
    chi   = 
      arma::accu((T1(indT1) - TRUEtemp(indT1)) % XTEMP(indT1, col0)) * parameters(Ent + T + p + 3) +
      arma::accu((T2(indT2) - (auxAlpha0(indT2) + auxAlpha1(indT2) % TRUEtemp(indT2))) % XTEMP(indT2, col0) % auxAlpha1(indT2)) * alpha(2 * nT2 + 6) + 
      arma::accu((Y1(indY1) - ZOOP(indY1)) % XTEMP(indY1, col0)) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
      arma::accu((Y2(indY2) - (auxLambda0(indY2) + auxLambda1(indY2) % ZOOP(indY2))) % XTEMP(indY2, col0) % auxLambda1(indY2)) * parameters(Ent + T + p) * lambda(2 * nY2 + 6) +
      na * nb;
    temp(N + nTEMP + 3) = R::rnorm(chi * delta, sqrt(delta));
    TRUEtemp += temp(N + nTEMP + 3) * XTEMP.col(0);
    ZOOP     += parameters(Ent + T + p) * temp(N + nTEMP + 3) * XTEMP.col(0);
    
    // cos (temp)
    TRUEtemp -= temp(N + nTEMP + 4) * XTEMP.col(1);
    ZOOP     -= parameters(Ent + T + p) * temp(N + nTEMP + 4) * XTEMP.col(1);
    delta = 1 / (
      sumT1Cos2 * parameters(Ent + T + p + 3) +
        arma::accu(pow(XTEMP(indT2, col1) % auxAlpha1(indT2), 2)) * alpha(2 * nT2 + 6) + 
        sumY1Cos2 * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
        arma::accu(pow(XTEMP(indY2, col1) % auxLambda1(indY2) * parameters(Ent + T + p), 2)) * lambda(2 * nY2 + 6) + 
        nb);
    chi   = 
      arma::accu((T1(indT1) - TRUEtemp(indT1)) % XTEMP(indT1, col1)) * parameters(Ent + T + p + 3) +
      arma::accu((T2(indT2) - (auxAlpha0(indT2) + auxAlpha1(indT2) % TRUEtemp(indT2))) % XTEMP(indT2, col1) % auxAlpha1(indT2)) * alpha(2 * nT2 + 6) + 
      arma::accu((Y1(indY1) - ZOOP(indY1)) % XTEMP(indY1, col1)) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
      arma::accu((Y2(indY2) - (auxLambda0(indY2) + auxLambda1(indY2) % ZOOP(indY2))) % XTEMP(indY2, col1) % auxLambda1(indY2)) * parameters(Ent + T + p) * lambda(2 * nY2 + 6) +
      na * nb;
    temp(N + nTEMP + 4) = R::rnorm(chi * delta, sqrt(delta));
    TRUEtemp += temp(N + nTEMP + 4) * XTEMP.col(1);
    ZOOP     += parameters(Ent + T + p) * temp(N + nTEMP + 4) * XTEMP.col(1);
    
    // year t (temp)
    rho2 = pow(temp(N + nTEMP + T + 6), 2);
    // tInd = 0
    TRUEtemp -= Xt.col(0) * temp(N + nTEMP + 6);
    ZOOP     -= Xt.col(0) * temp(N + nTEMP + 6) * parameters(Ent + T + p);
    delta = 1 / (
      NT1t(0) * parameters(Ent + T + p + 3) +
        arma::accu(pow(XtauxT2.col(0) % auxAlpha1(indT2), 2)) * alpha(2 * nT2 + 6) + 
        NY1t(0) * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
        arma::accu(pow(XtauxY2.col(0) % auxLambda1(indY2) * parameters(Ent + T + p), 2)) * lambda(2 * nY2 + 6) + 
        (1 + rho2) * temp(N + nTEMP + T + 2));
    chi   = 
      arma::accu((T1(indT1) - TRUEtemp(indT1)) % XtauxT1.col(0)) * parameters(Ent + T + p + 3) +
      arma::accu((T2(indT2) - (auxAlpha0(indT2) + auxAlpha1(indT2) % TRUEtemp(indT2))) % XtauxT2.col(0) % auxAlpha1(indT2)) * alpha(2 * nT2 + 6) + 
      arma::accu((Y1(indY1) - ZOOP(indY1)) % XtauxY1.col(0)) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
      arma::accu((Y2(indY2) - (auxLambda0(indY2) + auxLambda1(indY2) % ZOOP(indY2))) % XtauxY2.col(0) % auxLambda1(indY2)) * parameters(Ent + T + p) * lambda(2 * nY2 + 6) +
      (temp(N + nTEMP + T + 6) * (temp(N + nTEMP + 7) - temp(N + nTEMP)) + (1 + rho2) * temp(N + nTEMP)) * temp(N + nTEMP + T + 2);
    temp(N + nTEMP + 6) = R::rnorm(chi * delta, sqrt(delta));
    TRUEtemp += Xt.col(0) * temp(N + nTEMP + 6);
    ZOOP     += Xt.col(0) * temp(N + nTEMP + 6) * parameters(Ent + T + p);
    // tInd = 1,...,T-2
    for (int tInd = 1; tInd < T - 1; ++tInd) {
      TRUEtemp -= Xt.col(tInd) * temp(N + nTEMP + 6 + tInd);
      ZOOP     -= Xt.col(tInd) * temp(N + nTEMP + 6 + tInd) * parameters(Ent + T + p);
      delta = 1 / (
        NT1t(tInd) * parameters(Ent + T + p + 3) +
          arma::accu(pow(XtauxT2.col(tInd) % auxAlpha1(indT2), 2)) * alpha(2 * nT2 + 6) + 
          NY1t(tInd) * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
          arma::accu(pow(XtauxY2.col(tInd) % auxLambda1(indY2) * parameters(Ent + T + p), 2)) * lambda(2 * nY2 + 6) + 
          (1 + rho2) * temp(N + nTEMP + T + 2));
      chi   = 
        arma::accu((T1(indT1) - TRUEtemp(indT1)) % XtauxT1.col(tInd)) * parameters(Ent + T + p + 3) +
        arma::accu((T2(indT2) - (auxAlpha0(indT2) + auxAlpha1(indT2) % TRUEtemp(indT2))) % XtauxT2.col(tInd) % auxAlpha1(indT2)) * alpha(2 * nT2 + 6) + 
        arma::accu((Y1(indY1) - ZOOP(indY1)) % XtauxY1.col(tInd)) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
        arma::accu((Y2(indY2) - (auxLambda0(indY2) + auxLambda1(indY2) % ZOOP(indY2))) % XtauxY2.col(tInd) % auxLambda1(indY2)) * parameters(Ent + T + p) * lambda(2 * nY2 + 6) +
        (temp(N + nTEMP + T + 6) * (temp(N + nTEMP + 7 + tInd) + temp(N + nTEMP + 5 + tInd)) + pow(1 - temp(N + nTEMP + T + 6), 2) * temp(N + nTEMP)) * temp(N + nTEMP + T + 2);
      temp(N + nTEMP + 6 + tInd) = R::rnorm(chi * delta, sqrt(delta));
      TRUEtemp += Xt.col(tInd) * temp(N + nTEMP + 6 + tInd);
      ZOOP     += Xt.col(tInd) * temp(N + nTEMP + 6 + tInd) * parameters(Ent + T + p);
    }
    // tInd = T-1
    TRUEtemp -= Xt.col(T - 1) * temp(N + nTEMP + 5 + T);
    ZOOP     -= Xt.col(T - 1) * temp(N + nTEMP + 5 + T) * parameters(Ent + T + p);
    delta = 1 / (
      NT1t(T - 1) * parameters(Ent + T + p + 3) +
        arma::accu(pow(XtauxT2.col(T - 1) % auxAlpha1(indT2), 2)) * alpha(2 * nT2 + 6) + 
        NY1t(T - 1) * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
        arma::accu(pow(XtauxY2.col(T - 1) % auxLambda1(indY2) * parameters(Ent + T + p), 2)) * lambda(2 * nY2 + 6) + 
        temp(N + nTEMP + T + 2));
    chi   = 
      arma::accu((T1(indT1) - TRUEtemp(indT1)) % XtauxT1.col(T - 1)) * parameters(Ent + T + p + 3) +
      arma::accu((T2(indT2) - (auxAlpha0(indT2) + auxAlpha1(indT2) % TRUEtemp(indT2))) % XtauxT2.col(T - 1) % auxAlpha1(indT2)) * alpha(2 * nT2 + 6) + 
      arma::accu((Y1(indY1) - ZOOP(indY1)) % XtauxY1.col(T - 1)) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
      arma::accu((Y2(indY2) - (auxLambda0(indY2) + auxLambda1(indY2) % ZOOP(indY2))) % XtauxY2.col(T - 1) % auxLambda1(indY2)) * parameters(Ent + T + p) * lambda(2 * nY2 + 6) +
      (temp(N + nTEMP + T + 6) * temp(N + nTEMP + 4 + T) + (1 - temp(N + nTEMP + T + 6)) * temp(N + nTEMP)) * temp(N + nTEMP + T + 2);
    temp(N + nTEMP + 5 + T) = R::rnorm(chi * delta, sqrt(delta));
    TRUEtemp += Xt.col(T - 1) * temp(N + nTEMP + 5 + T);
    ZOOP     += Xt.col(T - 1) * temp(N + nTEMP + 5 + T) * parameters(Ent + T + p);
    
    // beta rho prec (year t in temp)
    initAR1 = AR1(temp(arma::span(N + nTEMP + 6, N + nTEMP + 5 + T)),
                  initAR1, na, nb, ga, gb, T);
    
    temp(N + nTEMP) = initAR1(0);
    temp(N + nTEMP + T + 6) = initAR1(1);
    temp(N + nTEMP + 2) = initAR1(2);
    
    // space si (temp)
    for (int iInd = 0; iInd < nTEMP; ++iInd) {
      TRUEtemp -= ITEMP.col(iInd) * temp(N + iInd);
      ZOOP     -= ITEMP.col(iInd) * temp(N + iInd) * parameters(Ent + T + p);
      process = temp(arma::span(N, N + nTEMP - 1));
      process.shed_row(iInd);
      r = RTEMP.row(iInd);
      r.shed_col(iInd);
      delta = 1 / (
        NT1s(iInd) * parameters(Ent + T + p + 3) +
          arma::accu(pow(XsauxT2.col(iInd) % auxAlpha1(indT2), 2)) * alpha(2 * nT2 + 6) + 
          NY1s(iInd) * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
          arma::accu(pow(XsauxY2.col(iInd) % auxLambda1(indY2) * parameters(Ent + T + p), 2)) * lambda(2 * nY2 + 6) + 
          RTEMP(iInd, iInd) * temp(N + nTEMP + 1));
      chi   = 
        arma::accu((T1(indT1) - TRUEtemp(indT1)) % XsauxT1.col(iInd)) * parameters(Ent + T + p + 3) +
        arma::accu((T2(indT2) - (auxAlpha0(indT2) + auxAlpha1(indT2) % TRUEtemp(indT2))) % XsauxT2.col(iInd) % auxAlpha1(indT2)) * alpha(2 * nT2 + 6) + 
        arma::accu((Y1(indY1) - ZOOP(indY1)) % XsauxY1.col(iInd)) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
        arma::accu((Y2(indY2) - (auxLambda0(indY2) + auxLambda1(indY2) % ZOOP(indY2))) % XsauxY2.col(iInd) % auxLambda1(indY2)) * parameters(Ent + T + p) * lambda(2 * nY2 + 6) +
        (- r * process).eval()(0,0) * temp(N + nTEMP + 1);
      temp(N + iInd) = R::rnorm(chi * delta, sqrt(delta));
      TRUEtemp += ITEMP.col(iInd) * temp(N + iInd);
      ZOOP     += ITEMP.col(iInd) * temp(N + iInd) * parameters(Ent + T + p);
    }
    
    // prec (space s in temp)
    process = temp(arma::span(N, N + nTEMP - 1));
    sumita = (process.t() * RTEMP * process).eval()(0,0); 
    temp(N + nTEMP + 1) = R::rgamma(nTEMP / 2 + ga, 1 / (sumita / 2 + gb));
    
    // temp
    temp(arma::span(0, N - 1)) = TRUEtemp;
    
    // alpha's and hiper, precT2
    alpha = calCor(
      T2.elem(indT2), temp.elem(indT2), IT2, sT2,
      dT2, alpha, sd,
      na, nb, ga, gb,
      nT2, NT2, 0, 1, 1, 10).t();
    
    alphas = alpha(2 * nT2) + alpha(2 * nT2 + 2) * alpha(arma::span(0, nT2 - 1));
    auxAlpha0.elem(indT2) = alphas.elem(sT2);
    alphas = alpha(2 * nT2 + 1) + alpha(2 * nT2 + 7) * alpha(arma::span(0, nT2 - 1)) + alpha(2 * nT2 + 3) * alpha(arma::span(nT2, 2 * nT2 - 1));
    auxAlpha1.elem(indT2) = alphas.elem(sT2);
    
    // precT1 
    auxN = T1 - temp(arma::span(0, N - 1)); 
    auxNT1 = auxN.elem(indT1); 
    parameters(Ent + T + p + 3) = R::rgamma(NT1 / 2 + ga, 1 / ((auxNT1.t() * auxNT1).eval()(0) / 2 + gb));
    
    // lambda's and hiper, precY2
    // a0's a1's mean0 mean1 prec0 prec1 decay0 decay1 precT2 a21 (a0t' a1t' rho0 rho1 v00 v11 v01)
    lambda = calCorTime(
      Y2.elem(indY2), ZOOP.elem(indY2), 
      Xt.rows(indY2), IY2, sY2,
      dY2, lambda, sd,
      na, nb, ga, gb,
      nY2, T, NY2, 
      0, 1, 1, 10).t();
    
    lambdas = lambda(2 * nY2 + 2) * lambda(arma::span(0, nY2 - 1));
    auxLambda0.elem(indY2) = Xt.rows(indY2) * lambda(arma::span(2 * nY2 + 8, 2 * nY2 + T + 7)) + lambdas.elem(sY2);
    lambdas = lambda(2 * nY2 + 7) * lambda(arma::span(0, nY2 - 1)) + lambda(2 * nY2 + 3) * lambda(arma::span(nY2, 2 * nY2 - 1));
    auxLambda1.elem(indY2) = Xt.rows(indY2) * lambda(arma::span(2 * nY2 + T + 8, 2 * nY2 + 2 * T + 7)) + lambdas.elem(sY2);
    
    // precY1 OK
    auxN = Y1 - ZOOP; 
    auxNY1 = auxN.elem(indY1); 
    parameters(Ent + T + p + 2) = R::rgamma(NY1 / 2 + ga, 1 / ((auxNY1.t() * auxNY1).eval()(0) / 2 + gb));
    
    // beta's OK
    for (int j = 0; j < p; ++j) {
      ZOOP -= parameters(Ent + T + j) * X.col(j);
      
      auxN   = pow(X.col(j), 2);
      auxNY1 = auxN.elem(indY1); 
      auxN   = pow(auxLambda1 % X.col(j), 2);
      auxNY2 = auxN.elem(indY2); 
      
      delta = 1 / (arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(2 * nY2 + 6) + nb);
      
      auxN   = X.col(j) % (Y1 - ZOOP); 
      auxNY1 = auxN.elem(indY1); 
      auxN   = auxLambda1 % X.col(j) % (Y2 - (auxLambda0 + auxLambda1 % ZOOP)); 
      auxNY2 = auxN.elem(indY2); 
      
      chi   = arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(2 * nY2 + 6) + na * nb;
      
      parameters(Ent + T + j) = R::rnorm(chi * delta, sqrt(delta));
      
      ZOOP += parameters(Ent + T + j) * X.col(j);
    }
    
    // beta1 temp OK
    ZOOP -= parameters(Ent + T + p) * temp(arma::span(0, N - 1));
    
    auxN   = pow(temp(arma::span(0, N - 1)), 2);
    auxNY1 = auxN.elem(indY1); 
    auxN   = pow(auxLambda1 % temp(arma::span(0, N - 1)), 2);
    auxNY2 = auxN.elem(indY2); 
    
    delta = 1 / (arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(2 * nY2 + 6) + nb);
    
    auxN   = temp(arma::span(0, N - 1)) % (Y1 - ZOOP); 
    auxNY1 = auxN.elem(indY1); 
    auxN   = auxLambda1 % temp(arma::span(0, N - 1)) % (Y2 - (auxLambda0 + auxLambda1 % ZOOP)); 
    auxNY2 = auxN.elem(indY2);
    
    chi   = arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(2 * nY2 + 6) + na * nb;
    
    parameters(Ent + T + p) = R::rnorm(chi * delta, sqrt(delta));
    
    ZOOP += parameters(Ent + T + p) * temp(arma::span(0, N - 1));
    
    // etats
    contador = 0;
    for (int tInd = 0; tInd < T; ++tInd){
      for (int iInd = 0; iInd < nt(tInd); ++iInd){
        
        ZOOP -= Xeta.col(contador) * parameters(contador);
        
        process = parameters(arma::span(CUMSUMnt(tInd), CUMSUMnt(tInd + 1) - 1));
        process.shed_row(iInd);
        r = R(arma::span(iInd, iInd), arma::span(0, nt(tInd) - 1), arma::span(tInd, tInd));
        r.shed_col(iInd);
        
        auxN = Xeta.col(contador);
        auxNY1 = auxN.elem(indY1);
        auxN = pow(auxLambda1 % auxN, 2);
        auxNY2 = auxN.elem(indY2);
        
        delta = 1 / (arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(2 * nY2 + 6) + R.slice(tInd)(iInd, iInd) * parameters(Ent + T + p + 4));
        
        auxN   = Xeta.col(contador) % (Y1 - ZOOP); 
        auxNY1 = auxN.elem(indY1); 
        auxN   = auxLambda1 % Xeta.col(contador) % (Y2 - (auxLambda0 + auxLambda1 % ZOOP)); 
        auxNY2 = auxN.elem(indY2); 
        
        chi   = arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(2 * nY2 + 6) + 
          (parameters(Ent + tInd) * R.slice(tInd)(iInd, iInd) + r * (parameters(Ent + tInd) - process)).eval()(0,0) * parameters(Ent + T + p + 4);
        
        parameters(contador) = R::rnorm(delta * chi, sqrt(delta));
        
        ZOOP += Xeta.col(contador) * parameters(contador);
        
        ++contador;
      }
    }
    
    //// etat OK
    for (int tInd = 0; tInd < T; ++tInd) {
      delta = 1 / (Rsum(tInd) * parameters(Ent + T + p + 4) + parameters(Ent + T + p + 5)); 
      chi   = (oneN(arma::span(0, nt(tInd) - 1)).t() * R.slice(tInd)(arma::span(0, nt(tInd) - 1), arma::span(0, nt(tInd) - 1)) * parameters(arma::span(CUMSUMnt(tInd), CUMSUMnt(tInd + 1) - 1))).eval()(0,0) * parameters(Ent + T + p + 4) + parameters(Ent + T + p + 1) * parameters(Ent + T + p + 5);
      parameters(Ent + tInd) = R::rnorm(chi * delta, sqrt(delta));
    }
    
    ////// beta0 OK
    delta = 1 / (T * parameters(Ent + T + p + 5) + nb);
    chi   = arma::accu(parameters(arma::span(Ent, Ent + T - 1))) * parameters(Ent + T + p + 5) + na * nb;
    parameters(Ent + T + p + 1) = R::rnorm(chi * delta, sqrt(delta));
    
    ////// precEtat OK
    auxT = parameters(arma::span(Ent, Ent + T - 1)) - parameters(Ent + T + p + 1);
    parameters(Ent + T + p + 5) = R::rgamma(T / 2 + ga, 1 / ((auxT.t() * auxT).eval()(0) / 2 + gb));
    
    //// precEtas OK
    process = parameters(arma::span(0, nt(0) - 1)) - parameters(Ent);
    sumita  = (process.t() * R.slice(0)(arma::span(0, nt(0) - 1),arma::span(0, nt(0) - 1)) * process).eval()(0,0);  
    for (int tInd = 1; tInd < T; ++tInd) {
      process = parameters(arma::span(CUMSUMnt(tInd), CUMSUMnt(tInd + 1) - 1)) - parameters(Ent + tInd);
      sumita += (process.t() * R.slice(tInd)(arma::span(0, nt(tInd) - 1),arma::span(0, nt(tInd) - 1)) * process).eval()(0,0);  
    }
    parameters(Ent + T + p + 4) = R::rgamma(Ent / 2 + ga, 1 / (sumita / 2 + gb));
    
    if ((b != 0) && (b % nReport == 0)) {
      Rcpp::Rcout << "The value of b : " << b << "\n";
    }
    
    if ((b > 0) && (b % nThin == 0)) {
      keep(b / nThin - 1, arma::span(0, Ent + T + p + 5)) = parameters.t();
      keep(b / nThin - 1, arma::span(Ent + T + p + 6, (Ent + T + p + 6) + (2 * nT2 + 7))) = alpha.t();
      keep(b / nThin - 1, arma::span((Ent + T + p + 6) + (2 * nT2 + 8), (Ent + T + p + 6) + (2 * nT2 + 8) + (2 * nY2 + 2 * T + 12))) = lambda.t();
      keep(b / nThin - 1, arma::span((Ent + T + p + 6) + (2 * nT2 + 8) + (2 * nY2 + 2 * T + 13), (Ent + T + p + 6) + (2 * nT2 + 8) + (2 * nY2 + 2 * T + 13) + (N + nTEMP + T + 6))) = temp.t();
    }
  }
  
  return keep;
}

// Both calibrations with coregionalization
// [[Rcpp::export]]
arma::mat GibbsZoop3Cpp(
    arma::vec Y1, arma::vec Y2, arma::vec T1, arma::vec T2,
    arma::mat X, arma::mat Xeta, arma::mat XTEMP, 
    arma::vec Rsum, double RsumTEMP, arma::cube R, arma::mat RTEMP, double dMax,
    arma::mat dY2, arma::mat dT2, double sd,
    arma::vec oneY1, arma::vec oneY2, arma::vec oneT1, arma::vec oneT2,
    arma::uvec indY1, arma::uvec indY2, arma::uvec indT1, arma::uvec indT2,
    arma::mat IY2, arma::mat IT2, arma::mat ITEMP, arma::uvec sY2, arma::uvec sT2, arma::uvec sTEMP, arma::mat Xt,
    int N, int NY1, int NY2, int NT1, int NT2, int T, arma::vec nt, double Ent, int nT2, int nY2, int nTEMP, int p,
    double na, double nb, double ga, double gb,
    int nBurnin, int maxIter, int nThin, int nReport) {
  
  // PARAMETERS //
  arma::vec temp(N + nTEMP + T + 7, arma::fill::ones); // temp'tls, mean's, mean, precmean, precyear, sin, cos, prec, year1 ... yearT, rhoyear
  arma::vec alpha(2 * nT2 + 8, arma::fill::ones); // a0's a1's mean0 mean1 prec0 prec1 decay0 decay1 precT2 a21
  arma::vec lambda(2 * nY2 + 8, arma::fill::ones); // l0's l1's mean0 mean1 prec0 prec1 decay0 decay1 precY2 a21
  arma::vec parameters(Ent + T + p + 6, arma::fill::ones); // eta'ts eta't X' beta1 beta0 precY1 precT1 precEtas precEtat
  arma::mat keep(maxIter / nThin, (Ent + T + p + 6) + (2 * nT2 + 8) + (2 * nY2 + 8) + (N + nTEMP + T + 7));
  
  alpha(2 * nT2 + 4) = 10 / dMax; // decay0
  alpha(2 * nT2 + 5) = 10 / dMax; // decay1
  lambda(2 * nY2 + 4) = 10 / dMax; // decay0
  lambda(2 * nY2 + 5) = 10 / dMax; // decay1
  
  //parameters(Ent + T + p + 3) = 10000; // precT1
  temp(N + nTEMP + 5) = INFINITY; // precTemp
  temp(N + nTEMP + T + 6) = 0; // rhoyear
  
  arma::vec TRUEtemp = ITEMP * temp(arma::span(N, N + nTEMP - 1)) + XTEMP * temp(arma::span(N + nTEMP + 3, N + nTEMP + 4)) + Xt * temp(arma::span(N + nTEMP + 6, N + nTEMP + T + 5));
  arma::vec ZOOP = Xeta * parameters(arma::span(0, Ent - 1)) + X * parameters(arma::span(Ent + T, Ent + T + p - 1)) + parameters(Ent + T + p) * TRUEtemp;
  
  // HELP //
  arma::mat Delta;
  arma::mat Chi;
  double delta;
  double chi;
  
  arma::vec alphas(nT2);
  arma::vec lambdas(nY2);
  
  arma::vec auxAlpha0(N, arma::fill::zeros); 
  alphas = alpha(2 * nT2) + alpha(2 * nT2 + 2) * alpha(arma::span(0, nT2 - 1));
  auxAlpha0.elem(indT2) = alphas.elem(sT2);
  arma::vec auxAlpha1(N, arma::fill::zeros); 
  alphas = alpha(2 * nT2 + 1) + alpha(2 * nT2 + 7) * alpha(arma::span(0, nT2 - 1)) + alpha(2 * nT2 + 3) * alpha(arma::span(nT2, 2 * nT2 - 1));
  auxAlpha1.elem(indT2) = alphas.elem(sT2);
  
  arma::vec auxLambda0(N, arma::fill::zeros); 
  lambdas = lambda(2 * nY2) + lambda(2 * nY2 + 2) * lambda(arma::span(0, nY2 - 1));
  auxLambda0.elem(indY2) = lambdas.elem(sY2);
  arma::vec auxLambda1(N, arma::fill::zeros); 
  lambdas = lambda(2 * nY2 + 1) + lambda(2 * nY2 + 7) * lambda(arma::span(0, nY2 - 1)) + lambda(2 * nY2 + 3) * lambda(arma::span(nY2, 2 * nY2 - 1));
  auxLambda1.elem(indY2) = lambdas.elem(sY2);
  
  arma::vec oneN(N, arma::fill::ones);
  arma::vec auxN(N);
  arma::vec auxT(T);
  arma::vec auxNY1(NY1);
  arma::vec auxNY2(NY2);
  arma::vec auxNT1(NT1);
  arma::vec auxNT2(NT2);
  
  arma::mat r;
  arma::vec process;
  double sumita;
  double contador;
  double rho2;
  
  arma::uvec col0 = { 0 };
  arma::uvec col1 = { 1 };
  
  double sumY1Sin2 = arma::accu(pow(XTEMP(indY1, col0), 2));
  double sumY1Cos2 = arma::accu(pow(XTEMP(indY1, col1), 2));
  double sumT1Sin2 = arma::accu(pow(XTEMP(indT1, col0), 2));
  double sumT1Cos2 = arma::accu(pow(XTEMP(indT1, col1), 2));
  
  arma::mat XtauxY1 = Xt.rows(indY1);
  arma::mat XtauxY2 = Xt.rows(indY2);
  arma::mat XtauxT1 = Xt.rows(indT1);
  arma::mat XtauxT2 = Xt.rows(indT2);
  
  arma::vec NY1t = arma::sum(XtauxY1).t();
  arma::vec NT1t = arma::sum(XtauxT1).t();
  
  arma::mat XsauxY1 = ITEMP.rows(indY1);
  arma::mat XsauxY2 = ITEMP.rows(indY2);
  arma::mat XsauxT1 = ITEMP.rows(indT1);
  arma::mat XsauxT2 = ITEMP.rows(indT2);
  
  arma::vec NY1s = arma::sum(XsauxY1).t();
  arma::vec NY2s = arma::sum(XsauxY2).t();
  arma::vec NT1s = arma::sum(XsauxT1).t();
  arma::vec NT2s = arma::sum(XsauxT2).t();
  
  arma::vec initAR1(3);
  initAR1(0) = temp(N + nTEMP);
  initAR1(1) = temp(N + nTEMP + T + 6);
  initAR1(2) = temp(N + nTEMP + 2);
  
  arma::vec CUMSUMnt(T + 1);
  CUMSUMnt(0) = 0;
  CUMSUMnt(arma::span(1, T)) = arma::cumsum(nt);
  
  // ITERATIONS //
  for (int b = 1 - nBurnin; b <= maxIter; ++b) {
    
    //// temp t l s OK
    //ZOOP -= parameters(Ent + T + p) * temp(arma::span(0, N - 1));
    //for (int j = 0; j < N; ++j) {
    //  delta = 1 / (oneT1(j) * parameters(Ent + T + p + 3) +
    //    oneT2(j) * pow(auxAlpha1(j), 2) * alpha(2 * nT2 + 6) + 
    //    oneY1(j) * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
    //    oneY2(j) * pow(auxLambda1(j) * parameters(Ent + T + p), 2) * lambda(2 * nY2 + 6) + 
    //    temp(N + nTEMP + 5));
    //  chi   = T1(j) * oneT1(j) * parameters(Ent + T + p + 3) +
    //    (T2(j) - auxAlpha0(j)) * oneT2(j) * auxAlpha1(j) * alpha(2 * nT2 + 6) + 
    //    (Y1(j) - ZOOP(j)) * oneY1(j) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
    //    (Y2(j) - (auxLambda0(j) + auxLambda1(j) * ZOOP(j))) * oneY2(j) * auxLambda1(j) * parameters(Ent + T + p) * lambda(2 * nY2 + 6) +
    //    ((ITEMP.row(j) * temp(arma::span(N, N + nTEMP - 1))).eval()(0) + temp(N + nTEMP + 3) * XTEMP(j,0) + temp(N + nTEMP + 4) * XTEMP(j,1) + (Xt.row(j) * temp(arma::span(N + nTEMP + 6, N + nTEMP + T + 5))).eval()(0)) * temp(N + nTEMP + 5);
    //  temp(j) = R::rnorm(chi * delta, sqrt(delta));
    //}
    //ZOOP += parameters(Ent + T + p) * temp(arma::span(0, N - 1));
    
    //temp(arma::span(N, N + nTEMP + T + 6)) = splmMR(
    //  temp(arma::span(0, N - 1)), XTEMP,
    //  ITEMP, sTEMP, RTEMP, RsumTEMP,
    //  temp(arma::span(N, N + nTEMP + T + 6)), 
    //  na, nb, ga, gb,
    //  2, nTEMP, N, T, Xt,
    //  0, 1, 1, 10
    //).t();
    
    // sin (temp)
    TRUEtemp -= temp(N + nTEMP + 3) * XTEMP.col(0);
    ZOOP     -= parameters(Ent + T + p) * temp(N + nTEMP + 3) * XTEMP.col(0);
    delta = 1 / (
      sumT1Sin2 * parameters(Ent + T + p + 3) +
      arma::accu(pow(XTEMP(indT2, col0) % auxAlpha1(indT2), 2)) * alpha(2 * nT2 + 6) + 
      sumY1Sin2 * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
      arma::accu(pow(XTEMP(indY2, col0) % auxLambda1(indY2) * parameters(Ent + T + p), 2)) * lambda(2 * nY2 + 6) + 
      nb);
    chi   = 
      arma::accu((T1(indT1) - TRUEtemp(indT1)) % XTEMP(indT1, col0)) * parameters(Ent + T + p + 3) +
      arma::accu((T2(indT2) - (auxAlpha0(indT2) + auxAlpha1(indT2) % TRUEtemp(indT2))) % XTEMP(indT2, col0) % auxAlpha1(indT2)) * alpha(2 * nT2 + 6) + 
      arma::accu((Y1(indY1) - ZOOP(indY1)) % XTEMP(indY1, col0)) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
      arma::accu((Y2(indY2) - (auxLambda0(indY2) + auxLambda1(indY2) % ZOOP(indY2))) % XTEMP(indY2, col0) % auxLambda1(indY2)) * parameters(Ent + T + p) * lambda(2 * nY2 + 6) +
      na * nb;
    temp(N + nTEMP + 3) = R::rnorm(chi * delta, sqrt(delta));
    TRUEtemp += temp(N + nTEMP + 3) * XTEMP.col(0);
    ZOOP     += parameters(Ent + T + p) * temp(N + nTEMP + 3) * XTEMP.col(0);
    
    // cos (temp)
    TRUEtemp -= temp(N + nTEMP + 4) * XTEMP.col(1);
    ZOOP     -= parameters(Ent + T + p) * temp(N + nTEMP + 4) * XTEMP.col(1);
    delta = 1 / (
      sumT1Cos2 * parameters(Ent + T + p + 3) +
      arma::accu(pow(XTEMP(indT2, col1) % auxAlpha1(indT2), 2)) * alpha(2 * nT2 + 6) + 
      sumY1Cos2 * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
      arma::accu(pow(XTEMP(indY2, col1) % auxLambda1(indY2) * parameters(Ent + T + p), 2)) * lambda(2 * nY2 + 6) + 
      nb);
    chi   = 
      arma::accu((T1(indT1) - TRUEtemp(indT1)) % XTEMP(indT1, col1)) * parameters(Ent + T + p + 3) +
      arma::accu((T2(indT2) - (auxAlpha0(indT2) + auxAlpha1(indT2) % TRUEtemp(indT2))) % XTEMP(indT2, col1) % auxAlpha1(indT2)) * alpha(2 * nT2 + 6) + 
      arma::accu((Y1(indY1) - ZOOP(indY1)) % XTEMP(indY1, col1)) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
      arma::accu((Y2(indY2) - (auxLambda0(indY2) + auxLambda1(indY2) % ZOOP(indY2))) % XTEMP(indY2, col1) % auxLambda1(indY2)) * parameters(Ent + T + p) * lambda(2 * nY2 + 6) +
      na * nb;
    temp(N + nTEMP + 4) = R::rnorm(chi * delta, sqrt(delta));
    TRUEtemp += temp(N + nTEMP + 4) * XTEMP.col(1);
    ZOOP     += parameters(Ent + T + p) * temp(N + nTEMP + 4) * XTEMP.col(1);
    
    // year t (temp)
    rho2 = pow(temp(N + nTEMP + T + 6), 2);
    // tInd = 0
    TRUEtemp -= Xt.col(0) * temp(N + nTEMP + 6);
    ZOOP     -= Xt.col(0) * temp(N + nTEMP + 6) * parameters(Ent + T + p);
    delta = 1 / (
      NT1t(0) * parameters(Ent + T + p + 3) +
      arma::accu(pow(XtauxT2.col(0) % auxAlpha1(indT2), 2)) * alpha(2 * nT2 + 6) + 
      NY1t(0) * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
      arma::accu(pow(XtauxY2.col(0) % auxLambda1(indY2) * parameters(Ent + T + p), 2)) * lambda(2 * nY2 + 6) + 
      (1 + rho2) * temp(N + nTEMP + T + 2));
    chi   = 
      arma::accu((T1(indT1) - TRUEtemp(indT1)) % XtauxT1.col(0)) * parameters(Ent + T + p + 3) +
      arma::accu((T2(indT2) - (auxAlpha0(indT2) + auxAlpha1(indT2) % TRUEtemp(indT2))) % XtauxT2.col(0) % auxAlpha1(indT2)) * alpha(2 * nT2 + 6) + 
      arma::accu((Y1(indY1) - ZOOP(indY1)) % XtauxY1.col(0)) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
      arma::accu((Y2(indY2) - (auxLambda0(indY2) + auxLambda1(indY2) % ZOOP(indY2))) % XtauxY2.col(0) % auxLambda1(indY2)) * parameters(Ent + T + p) * lambda(2 * nY2 + 6) +
      (temp(N + nTEMP + T + 6) * (temp(N + nTEMP + 7) - temp(N + nTEMP)) + (1 + rho2) * temp(N + nTEMP)) * temp(N + nTEMP + T + 2);
    temp(N + nTEMP + 6) = R::rnorm(chi * delta, sqrt(delta));
    TRUEtemp += Xt.col(0) * temp(N + nTEMP + 6);
    ZOOP     += Xt.col(0) * temp(N + nTEMP + 6) * parameters(Ent + T + p);
    // tInd = 1,...,T-2
    for (int tInd = 1; tInd < T - 1; ++tInd) {
      TRUEtemp -= Xt.col(tInd) * temp(N + nTEMP + 6 + tInd);
      ZOOP     -= Xt.col(tInd) * temp(N + nTEMP + 6 + tInd) * parameters(Ent + T + p);
      delta = 1 / (
        NT1t(tInd) * parameters(Ent + T + p + 3) +
        arma::accu(pow(XtauxT2.col(tInd) % auxAlpha1(indT2), 2)) * alpha(2 * nT2 + 6) + 
        NY1t(tInd) * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
        arma::accu(pow(XtauxY2.col(tInd) % auxLambda1(indY2) * parameters(Ent + T + p), 2)) * lambda(2 * nY2 + 6) + 
        (1 + rho2) * temp(N + nTEMP + T + 2));
      chi   = 
        arma::accu((T1(indT1) - TRUEtemp(indT1)) % XtauxT1.col(tInd)) * parameters(Ent + T + p + 3) +
        arma::accu((T2(indT2) - (auxAlpha0(indT2) + auxAlpha1(indT2) % TRUEtemp(indT2))) % XtauxT2.col(tInd) % auxAlpha1(indT2)) * alpha(2 * nT2 + 6) + 
        arma::accu((Y1(indY1) - ZOOP(indY1)) % XtauxY1.col(tInd)) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
        arma::accu((Y2(indY2) - (auxLambda0(indY2) + auxLambda1(indY2) % ZOOP(indY2))) % XtauxY2.col(tInd) % auxLambda1(indY2)) * parameters(Ent + T + p) * lambda(2 * nY2 + 6) +
        (temp(N + nTEMP + T + 6) * (temp(N + nTEMP + 7 + tInd) + temp(N + nTEMP + 5 + tInd)) + pow(1 - temp(N + nTEMP + T + 6), 2) * temp(N + nTEMP)) * temp(N + nTEMP + T + 2);
      temp(N + nTEMP + 6 + tInd) = R::rnorm(chi * delta, sqrt(delta));
      TRUEtemp += Xt.col(tInd) * temp(N + nTEMP + 6 + tInd);
      ZOOP     += Xt.col(tInd) * temp(N + nTEMP + 6 + tInd) * parameters(Ent + T + p);
    }
    // tInd = T-1
    TRUEtemp -= Xt.col(T - 1) * temp(N + nTEMP + 5 + T);
    ZOOP     -= Xt.col(T - 1) * temp(N + nTEMP + 5 + T) * parameters(Ent + T + p);
    delta = 1 / (
      NT1t(T - 1) * parameters(Ent + T + p + 3) +
      arma::accu(pow(XtauxT2.col(T - 1) % auxAlpha1(indT2), 2)) * alpha(2 * nT2 + 6) + 
      NY1t(T - 1) * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
      arma::accu(pow(XtauxY2.col(T - 1) % auxLambda1(indY2) * parameters(Ent + T + p), 2)) * lambda(2 * nY2 + 6) + 
      temp(N + nTEMP + T + 2));
    chi   = 
      arma::accu((T1(indT1) - TRUEtemp(indT1)) % XtauxT1.col(T - 1)) * parameters(Ent + T + p + 3) +
      arma::accu((T2(indT2) - (auxAlpha0(indT2) + auxAlpha1(indT2) % TRUEtemp(indT2))) % XtauxT2.col(T - 1) % auxAlpha1(indT2)) * alpha(2 * nT2 + 6) + 
      arma::accu((Y1(indY1) - ZOOP(indY1)) % XtauxY1.col(T - 1)) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
      arma::accu((Y2(indY2) - (auxLambda0(indY2) + auxLambda1(indY2) % ZOOP(indY2))) % XtauxY2.col(T - 1) % auxLambda1(indY2)) * parameters(Ent + T + p) * lambda(2 * nY2 + 6) +
      (temp(N + nTEMP + T + 6) * temp(N + nTEMP + 4 + T) + (1 - temp(N + nTEMP + T + 6)) * temp(N + nTEMP)) * temp(N + nTEMP + T + 2);
    temp(N + nTEMP + 5 + T) = R::rnorm(chi * delta, sqrt(delta));
    TRUEtemp += Xt.col(T - 1) * temp(N + nTEMP + 5 + T);
    ZOOP     += Xt.col(T - 1) * temp(N + nTEMP + 5 + T) * parameters(Ent + T + p);
    
    // beta rho prec (year t in temp)
    initAR1 = AR1(temp(arma::span(N + nTEMP + 6, N + nTEMP + 5 + T)),
                  initAR1, na, nb, ga, gb, T);
    
    temp(N + nTEMP) = initAR1(0);
    temp(N + nTEMP + T + 6) = initAR1(1);
    temp(N + nTEMP + 2) = initAR1(2);
    
    // space si (temp)
    for (int iInd = 0; iInd < nTEMP; ++iInd) {
      TRUEtemp -= ITEMP.col(iInd) * temp(N + iInd);
      ZOOP     -= ITEMP.col(iInd) * temp(N + iInd) * parameters(Ent + T + p);
      process = temp(arma::span(N, N + nTEMP - 1));
      process.shed_row(iInd);
      r = RTEMP.row(iInd);
      r.shed_col(iInd);
      delta = 1 / (
        NT1s(iInd) * parameters(Ent + T + p + 3) +
        arma::accu(pow(XsauxT2.col(iInd) % auxAlpha1(indT2), 2)) * alpha(2 * nT2 + 6) + 
        NY1s(iInd) * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
        arma::accu(pow(XsauxY2.col(iInd) % auxLambda1(indY2) * parameters(Ent + T + p), 2)) * lambda(2 * nY2 + 6) + 
        RTEMP(iInd, iInd) * temp(N + nTEMP + 1));
      chi   = 
        arma::accu((T1(indT1) - TRUEtemp(indT1)) % XsauxT1.col(iInd)) * parameters(Ent + T + p + 3) +
        arma::accu((T2(indT2) - (auxAlpha0(indT2) + auxAlpha1(indT2) % TRUEtemp(indT2))) % XsauxT2.col(iInd) % auxAlpha1(indT2)) * alpha(2 * nT2 + 6) + 
        arma::accu((Y1(indY1) - ZOOP(indY1)) % XsauxY1.col(iInd)) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
        arma::accu((Y2(indY2) - (auxLambda0(indY2) + auxLambda1(indY2) % ZOOP(indY2))) % XsauxY2.col(iInd) % auxLambda1(indY2)) * parameters(Ent + T + p) * lambda(2 * nY2 + 6) +
        (- r * process).eval()(0,0) * temp(N + nTEMP + 1);
      temp(N + iInd) = R::rnorm(chi * delta, sqrt(delta));
      TRUEtemp += ITEMP.col(iInd) * temp(N + iInd);
      ZOOP     += ITEMP.col(iInd) * temp(N + iInd) * parameters(Ent + T + p);
    }
    
    // prec (space s in temp)
    process = temp(arma::span(N, N + nTEMP - 1));
    sumita = (process.t() * RTEMP * process).eval()(0,0); 
    temp(N + nTEMP + 1) = R::rgamma(nTEMP / 2 + ga, 1 / (sumita / 2 + gb));
    
    // temp
    temp(arma::span(0, N - 1)) = TRUEtemp;
    
    // alpha's and hiper, precT2
    alpha = calCor(
      T2.elem(indT2), temp.elem(indT2), IT2, sT2,
      dT2, alpha, sd,
      na, nb, ga, gb,
      nT2, NT2, 0, 1, 1, 10).t();
    
    alphas = alpha(2 * nT2) + alpha(2 * nT2 + 2) * alpha(arma::span(0, nT2 - 1));
    auxAlpha0.elem(indT2) = alphas.elem(sT2);
    alphas = alpha(2 * nT2 + 1) + alpha(2 * nT2 + 7) * alpha(arma::span(0, nT2 - 1)) + alpha(2 * nT2 + 3) * alpha(arma::span(nT2, 2 * nT2 - 1));
    auxAlpha1.elem(indT2) = alphas.elem(sT2);
    
    // precT1 
    auxN = T1 - temp(arma::span(0, N - 1)); 
    auxNT1 = auxN.elem(indT1); 
    parameters(Ent + T + p + 3) = R::rgamma(NT1 / 2 + ga, 1 / ((auxNT1.t() * auxNT1).eval()(0) / 2 + gb));
    
    // lambda's and hiper, precY2
    lambda = calCor(
      Y2.elem(indY2), ZOOP.elem(indY2), IY2, sY2,
      dY2, lambda, sd,
      na, nb, ga, gb,
      nY2, NY2, 0, 1, 1, 10).t();
    
    lambdas = lambda(2 * nY2) + lambda(2 * nY2 + 2) * lambda(arma::span(0, nY2 - 1));
    auxLambda0.elem(indY2) = lambdas.elem(sY2);
    lambdas = lambda(2 * nY2 + 1) + lambda(2 * nY2 + 7) * lambda(arma::span(0, nY2 - 1)) + lambda(2 * nY2 + 3) * lambda(arma::span(nY2, 2 * nY2 - 1));
    auxLambda1.elem(indY2) = lambdas.elem(sY2);
    
    // precY1 OK
    auxN = Y1 - ZOOP; 
    auxNY1 = auxN.elem(indY1); 
    parameters(Ent + T + p + 2) = R::rgamma(NY1 / 2 + ga, 1 / ((auxNY1.t() * auxNY1).eval()(0) / 2 + gb));
    
    // beta's OK
    for (int j = 0; j < p; ++j) {
      ZOOP -= parameters(Ent + T + j) * X.col(j);
      
      auxN   = pow(X.col(j), 2);
      auxNY1 = auxN.elem(indY1); 
      auxN   = pow(auxLambda1 % X.col(j), 2);
      auxNY2 = auxN.elem(indY2); 
      
      delta = 1 / (arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(2 * nY2 + 6) + nb);
      
      auxN   = X.col(j) % (Y1 - ZOOP); 
      auxNY1 = auxN.elem(indY1); 
      auxN   = auxLambda1 % X.col(j) % (Y2 - (auxLambda0 + auxLambda1 % ZOOP)); 
      auxNY2 = auxN.elem(indY2); 
      
      chi   = arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(2 * nY2 + 6) + na * nb;
      
      parameters(Ent + T + j) = R::rnorm(chi * delta, sqrt(delta));
      
      ZOOP += parameters(Ent + T + j) * X.col(j);
    }
    
    // beta1 temp OK
    ZOOP -= parameters(Ent + T + p) * temp(arma::span(0, N - 1));
    
    auxN   = pow(temp(arma::span(0, N - 1)), 2);
    auxNY1 = auxN.elem(indY1); 
    auxN   = pow(auxLambda1 % temp(arma::span(0, N - 1)), 2);
    auxNY2 = auxN.elem(indY2); 
    
    delta = 1 / (arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(2 * nY2 + 6) + nb);
    
    auxN   = temp(arma::span(0, N - 1)) % (Y1 - ZOOP); 
    auxNY1 = auxN.elem(indY1); 
    auxN   = auxLambda1 % temp(arma::span(0, N - 1)) % (Y2 - (auxLambda0 + auxLambda1 % ZOOP)); 
    auxNY2 = auxN.elem(indY2);
    
    chi   = arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(2 * nY2 + 6) + na * nb;
    
    parameters(Ent + T + p) = R::rnorm(chi * delta, sqrt(delta));
    
    ZOOP += parameters(Ent + T + p) * temp(arma::span(0, N - 1));
    
    // etats
    contador = 0;
    for (int tInd = 0; tInd < T; ++tInd){
      for (int iInd = 0; iInd < nt(tInd); ++iInd){
        
        ZOOP -= Xeta.col(contador) * parameters(contador);
        
        process = parameters(arma::span(CUMSUMnt(tInd), CUMSUMnt(tInd + 1) - 1));
        process.shed_row(iInd);
        r = R(arma::span(iInd, iInd), arma::span(0, nt(tInd) - 1), arma::span(tInd, tInd));
        r.shed_col(iInd);
        
        auxN = Xeta.col(contador);
        auxNY1 = auxN.elem(indY1);
        auxN = pow(auxLambda1 % auxN, 2);
        auxNY2 = auxN.elem(indY2);
        
        delta = 1 / (arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(2 * nY2 + 6) + R.slice(tInd)(iInd, iInd) * parameters(Ent + T + p + 4));
        
        auxN   = Xeta.col(contador) % (Y1 - ZOOP); 
        auxNY1 = auxN.elem(indY1); 
        auxN   = auxLambda1 % Xeta.col(contador) % (Y2 - (auxLambda0 + auxLambda1 % ZOOP)); 
        auxNY2 = auxN.elem(indY2); 
        
        chi   = arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(2 * nY2 + 6) + 
          (parameters(Ent + tInd) * R.slice(tInd)(iInd, iInd) + r * (parameters(Ent + tInd) - process)).eval()(0,0) * parameters(Ent + T + p + 4);
        
        parameters(contador) = R::rnorm(delta * chi, sqrt(delta));
        
        ZOOP += Xeta.col(contador) * parameters(contador);
        
        ++contador;
      }
    }
    
    //// etat OK
    for (int tInd = 0; tInd < T; ++tInd) {
      delta = 1 / (Rsum(tInd) * parameters(Ent + T + p + 4) + parameters(Ent + T + p + 5)); 
      chi   = (oneN(arma::span(0, nt(tInd) - 1)).t() * R.slice(tInd)(arma::span(0, nt(tInd) - 1), arma::span(0, nt(tInd) - 1)) * parameters(arma::span(CUMSUMnt(tInd), CUMSUMnt(tInd + 1) - 1))).eval()(0,0) * parameters(Ent + T + p + 4) + parameters(Ent + T + p + 1) * parameters(Ent + T + p + 5);
      parameters(Ent + tInd) = R::rnorm(chi * delta, sqrt(delta));
    }
    
    ////// beta0 OK
    delta = 1 / (T * parameters(Ent + T + p + 5) + nb);
    chi   = arma::accu(parameters(arma::span(Ent, Ent + T - 1))) * parameters(Ent + T + p + 5) + na * nb;
    parameters(Ent + T + p + 1) = R::rnorm(chi * delta, sqrt(delta));
    
    ////// precEtat OK
    auxT = parameters(arma::span(Ent, Ent + T - 1)) - parameters(Ent + T + p + 1);
    parameters(Ent + T + p + 5) = R::rgamma(T / 2 + ga, 1 / ((auxT.t() * auxT).eval()(0) / 2 + gb));
    
    //// precEtas OK
    process = parameters(arma::span(0, nt(0) - 1)) - parameters(Ent);
    sumita  = (process.t() * R.slice(0)(arma::span(0, nt(0) - 1),arma::span(0, nt(0) - 1)) * process).eval()(0,0);  
    for (int tInd = 1; tInd < T; ++tInd) {
      process = parameters(arma::span(CUMSUMnt(tInd), CUMSUMnt(tInd + 1) - 1)) - parameters(Ent + tInd);
      sumita += (process.t() * R.slice(tInd)(arma::span(0, nt(tInd) - 1),arma::span(0, nt(tInd) - 1)) * process).eval()(0,0);  
    }
    parameters(Ent + T + p + 4) = R::rgamma(Ent / 2 + ga, 1 / (sumita / 2 + gb));
    
    if ((b != 0) && (b % nReport == 0)) {
      Rcpp::Rcout << "The value of b : " << b << "\n";
    }
    
    if ((b > 0) && (b % nThin == 0)) {
      keep(b / nThin - 1, arma::span(0, Ent + T + p + 5)) = parameters.t();
      keep(b / nThin - 1, arma::span(Ent + T + p + 6, (Ent + T + p + 6) + (2 * nT2 + 7))) = alpha.t();
      keep(b / nThin - 1, arma::span((Ent + T + p + 6) + (2 * nT2 + 8), (Ent + T + p + 6) + (2 * nT2 + 8) + (2 * nY2 + 7))) = lambda.t();
      keep(b / nThin - 1, arma::span((Ent + T + p + 6) + (2 * nT2 + 8) + (2 * nY2 + 8), (Ent + T + p + 6) + (2 * nT2 + 8) + (2 * nY2 + 8) + (N + nTEMP + T + 6))) = temp.t();
    }
  }
  
  return keep;
}

// Both calibrations with 2 GPs (intercept and slope) 
// [[Rcpp::export]]
arma::mat GibbsZoop2Cpp(
    arma::vec Y1, arma::vec Y2, arma::vec T1, arma::vec T2,
    arma::mat X, arma::mat Xeta, arma::mat XTEMP, 
    arma::vec Rsum, double RsumTEMP, arma::cube R, arma::mat RTEMP, double dMax,
    arma::mat dY2, arma::mat dT2, double sd,
    arma::vec oneY1, arma::vec oneY2, arma::vec oneT1, arma::vec oneT2,
    arma::uvec indY1, arma::uvec indY2, arma::uvec indT1, arma::uvec indT2,
    arma::mat IY2, arma::mat IT2, arma::mat ITEMP, arma::uvec sY2, arma::uvec sT2, arma::uvec sTEMP, arma::mat Xt,
    int N, int NY1, int NY2, int NT1, int NT2, int T, arma::vec nt, double Ent, int nT2, int nY2, int nTEMP, int p,
    double na, double nb, double ga, double gb,
    int nBurnin, int maxIter, int nThin, int nReport) {
  
  // PARAMETERS //
  arma::vec temp(N + nTEMP + T + 7, arma::fill::ones); // temp'tls, mean's, mean, precmean, precyear, sin, cos, prec, year1 ... yearT, rhoyear
  arma::vec alpha(2 * nT2 + 7, arma::fill::ones); // a0's a1's mean0 mean1 prec0 prec1 decay0 decay1 precT2 
  arma::vec lambda(2 * nY2 + 7, arma::fill::ones); // l0's l1's mean0 mean1 prec0 prec1 decay0 decay1 precY2
  arma::vec parameters(Ent + T + p + 6, arma::fill::ones); // eta'ts eta't X' beta1 beta0 precY1 precT1 precEtas precEtat
  arma::mat keep(maxIter / nThin, (Ent + T + p + 6) + (2 * nT2 + 7) + (2 * nY2 + 7) + (N + nTEMP + T + 7));
  
  alpha(2 * nT2 + 4) = 10 / dMax; // decay0
  alpha(2 * nT2 + 5) = 10 / dMax; // decay1
  lambda(2 * nY2 + 4) = 10 / dMax; // decay0
  lambda(2 * nY2 + 5) = 10 / dMax; // decay1
  double Umin = 3 / dMax;
  double Umax = 3 / (0.1 * dMax);
  
  //parameters(Ent + T + p + 3) = 10000; // precT1
  temp(N + nTEMP + 5) = INFINITY; // precTemp
  temp(N + nTEMP + T + 6) = 0; // rhoyear
  
  arma::vec ZOOP = Xeta * parameters(arma::span(0, Ent - 1)) + X * parameters(arma::span(Ent + T, Ent + T + p - 1)) + parameters(Ent + T + p) * temp(arma::span(0, N - 1));
  
  // HELP //
  arma::mat Delta;
  arma::mat Chi;
  double delta;
  double chi;
  
  arma::vec alphas(nT2);
  arma::vec lambdas(nY2);
  
  arma::vec auxAlpha0(N, arma::fill::zeros); 
  alphas = alpha(arma::span(0, nT2 - 1));
  auxAlpha0.elem(indT2) = alphas.elem(sT2);
  arma::vec auxAlpha1(N, arma::fill::zeros); 
  alphas = alpha(arma::span(nT2, 2 * nT2 - 1));
  auxAlpha1.elem(indT2) = alphas.elem(sT2);
  
  arma::vec auxLambda0(N, arma::fill::zeros); 
  lambdas = lambda(arma::span(0, nY2 - 1));
  auxLambda0.elem(indY2) = lambdas.elem(sY2);
  arma::vec auxLambda1(N, arma::fill::zeros); 
  lambdas = lambda(arma::span(nY2, 2 * nY2 - 1));
  auxLambda1.elem(indY2) = lambdas.elem(sY2);
  
  arma::vec oneN(N, arma::fill::ones);
  arma::vec auxN(N);
  arma::vec auxT(T);
  arma::vec auxNY1(NY1);
  arma::vec auxNY2(NY2);
  arma::vec auxNT1(NT1);
  arma::vec auxNT2(NT2);
  
  arma::mat r;
  arma::vec process;
  double sumita;
  double contador;
  
  arma::vec CUMSUMnt(T + 1);
  CUMSUMnt(0) = 0;
  CUMSUMnt(arma::span(1, T)) = arma::cumsum(nt);
  
  // ITERATIONS //
  for (int b = 1 - nBurnin; b <= maxIter; ++b) {
    
    // temp t l s OK
    ZOOP -= parameters(Ent + T + p) * temp(arma::span(0, N - 1));
    for (int j = 0; j < N; ++j) {
      delta = 1 / (oneT1(j) * parameters(Ent + T + p + 3) +
        oneT2(j) * pow(auxAlpha1(j), 2) * alpha(2 * nT2 + 6) + 
        oneY1(j) * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
        oneY2(j) * pow(auxLambda1(j) * parameters(Ent + T + p), 2) * lambda(2 * nY2 + 6) + 
        temp(N + nTEMP + 5));
      chi   = T1(j) * oneT1(j) * parameters(Ent + T + p + 3) +
        (T2(j) - auxAlpha0(j)) * oneT2(j) * auxAlpha1(j) * alpha(2 * nT2 + 6) + 
        (Y1(j) - ZOOP(j)) * oneY1(j) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
        (Y2(j) - (auxLambda0(j) + auxLambda1(j) * ZOOP(j))) * oneY2(j) * auxLambda1(j) * parameters(Ent + T + p) * lambda(2 * nY2 + 6) +
        ((ITEMP.row(j) * temp(arma::span(N, N + nTEMP - 1))).eval()(0) + temp(N + nTEMP + 3) * XTEMP(j,0) + temp(N + nTEMP + 4) * XTEMP(j,1) + (Xt.row(j) * temp(arma::span(N + nTEMP + 6, N + nTEMP + T + 5))).eval()(0)) * temp(N + nTEMP + 5);
      temp(j) = R::rnorm(chi * delta, sqrt(delta));
    }
    ZOOP += parameters(Ent + T + p) * temp(arma::span(0, N - 1));
    
    // mean temp OK
    //delta = 1 / (N * temp(N + 1) + nb);
    //chi   = arma::accu(temp(arma::span(0, N - 1))) * temp(N + 1) + na * nb;
    //temp(N) = R::rnorm(chi * delta, sqrt(delta));
    
    // prec temp OK
    //auxN = temp(arma::span(0, N - 1)) - temp(N);
    //temp(N + 1) = R::rgamma(N / 2 + ga, 1 / ((auxN.t() * auxN).eval()(0) / 2 + gb));
    
    temp(arma::span(N, N + nTEMP + T + 6)) = splmMR(
      temp(arma::span(0, N - 1)), XTEMP,
      ITEMP, sTEMP, RTEMP, RsumTEMP,
      temp(arma::span(N, N + nTEMP + T + 6)), 
      na, nb, ga, gb,
      2, nTEMP, N, T, Xt,
      0, 1, 1, 10
    ).t();
    
    // alpha's and hiper, precT2
    alpha = splmSS(
      T2.elem(indT2), temp.elem(indT2), IT2, sT2,
      dT2, alpha, sd,
      na, nb, ga, gb, Umin, Umax,
      nT2, NT2, 0, 1, 1, 10).t();
    
    alphas = alpha(arma::span(0, nT2 - 1));
    auxAlpha0.elem(indT2) = alphas.elem(sT2);
    alphas = alpha(arma::span(nT2, 2 * nT2 - 1));
    auxAlpha1.elem(indT2) = alphas.elem(sT2);
    
    // precT1 
    auxN = T1 - temp(arma::span(0, N - 1)); 
    auxNT1 = auxN.elem(indT1); 
    parameters(Ent + T + p + 3) = R::rgamma(NT1 / 2 + ga, 1 / ((auxNT1.t() * auxNT1).eval()(0) / 2 + gb));
    
    // lambda's and hiper, precY2
    lambda = splmSS(
      Y2.elem(indY2), ZOOP.elem(indY2), IY2, sY2,
      dY2, lambda, sd,
      na, nb, ga, gb, Umin, Umax,
      nY2, NY2, 0, 1, 1, 10).t();
    
    lambdas = lambda(arma::span(0, nY2 - 1));
    auxLambda0.elem(indY2) = lambdas.elem(sY2);
    lambdas = lambda(arma::span(nY2, 2 * nY2 - 1));
    auxLambda1.elem(indY2) = lambdas.elem(sY2);
    
    // precY1 OK
    auxN = Y1 - ZOOP; 
    auxNY1 = auxN.elem(indY1); 
    parameters(Ent + T + p + 2) = R::rgamma(NY1 / 2 + ga, 1 / ((auxNY1.t() * auxNY1).eval()(0) / 2 + gb));
    
    // beta's OK
    for (int j = 0; j < p; ++j) {
      ZOOP -= parameters(Ent + T + j) * X.col(j);
      
      auxN   = pow(X.col(j), 2);
      auxNY1 = auxN.elem(indY1); 
      auxN   = pow(auxLambda1 % X.col(j), 2);
      auxNY2 = auxN.elem(indY2); 
      
      delta = 1 / (arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(2 * nY2 + 6) + nb);
      
      auxN   = X.col(j) % (Y1 - ZOOP); 
      auxNY1 = auxN.elem(indY1); 
      auxN   = auxLambda1 % X.col(j) % (Y2 - (auxLambda0 + auxLambda1 % ZOOP)); 
      auxNY2 = auxN.elem(indY2); 
      
      chi   = arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(2 * nY2 + 6) + na * nb;
      
      parameters(Ent + T + j) = R::rnorm(chi * delta, sqrt(delta));
      
      ZOOP += parameters(Ent + T + j) * X.col(j);
    }
    
    // beta1 temp OK
    ZOOP -= parameters(Ent + T + p) * temp(arma::span(0, N - 1));
    
    auxN   = pow(temp(arma::span(0, N - 1)), 2);
    auxNY1 = auxN.elem(indY1); 
    auxN   = pow(auxLambda1 % temp(arma::span(0, N - 1)), 2);
    auxNY2 = auxN.elem(indY2); 
    
    delta = 1 / (arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(2 * nY2 + 6) + nb);
    
    auxN   = temp(arma::span(0, N - 1)) % (Y1 - ZOOP); 
    auxNY1 = auxN.elem(indY1); 
    auxN   = auxLambda1 % temp(arma::span(0, N - 1)) % (Y2 - (auxLambda0 + auxLambda1 % ZOOP)); 
    auxNY2 = auxN.elem(indY2);
    
    chi   = arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(2 * nY2 + 6) + na * nb;
    
    parameters(Ent + T + p) = R::rnorm(chi * delta, sqrt(delta));
    
    ZOOP += parameters(Ent + T + p) * temp(arma::span(0, N - 1));
    
    // etats
    contador = 0;
    for (int tInd = 0; tInd < T; ++tInd){
      for (int iInd = 0; iInd < nt(tInd); ++iInd){
        
        ZOOP -= Xeta.col(contador) * parameters(contador);
        
        process = parameters(arma::span(CUMSUMnt(tInd), CUMSUMnt(tInd + 1) - 1));
        process.shed_row(iInd);
        r = R(arma::span(iInd, iInd), arma::span(0, nt(tInd) - 1), arma::span(tInd, tInd));
        r.shed_col(iInd);
        
        auxN = Xeta.col(contador);
        auxNY1 = auxN.elem(indY1);
        auxN = pow(auxLambda1 % auxN, 2);
        auxNY2 = auxN.elem(indY2);
        
        delta = 1 / (arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(2 * nY2 + 6) + R.slice(tInd)(iInd, iInd) * parameters(Ent + T + p + 4));
        
        auxN   = Xeta.col(contador) % (Y1 - ZOOP); 
        auxNY1 = auxN.elem(indY1); 
        auxN   = auxLambda1 % Xeta.col(contador) % (Y2 - (auxLambda0 + auxLambda1 % ZOOP)); 
        auxNY2 = auxN.elem(indY2); 
        
        chi   = arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(2 * nY2 + 6) + 
          (parameters(Ent + tInd) * R.slice(tInd)(iInd, iInd) + r * (parameters(Ent + tInd) - process)).eval()(0,0) * parameters(Ent + T + p + 4);
        
        parameters(contador) = R::rnorm(delta * chi, sqrt(delta));
        
        ZOOP += Xeta.col(contador) * parameters(contador);
        
        ++contador;
      }
    }
    
    //// etat OK
    for (int tInd = 0; tInd < T; ++tInd) {
      delta = 1 / (Rsum(tInd) * parameters(Ent + T + p + 4) + parameters(Ent + T + p + 5)); 
      chi   = (oneN(arma::span(0, nt(tInd) - 1)).t() * R.slice(tInd)(arma::span(0, nt(tInd) - 1), arma::span(0, nt(tInd) - 1)) * parameters(arma::span(CUMSUMnt(tInd), CUMSUMnt(tInd + 1) - 1))).eval()(0,0) * parameters(Ent + T + p + 4) + parameters(Ent + T + p + 1) * parameters(Ent + T + p + 5);
      parameters(Ent + tInd) = R::rnorm(chi * delta, sqrt(delta));
    }
    
    ////// beta0 OK
    delta = 1 / (T * parameters(Ent + T + p + 5) + nb);
    chi   = arma::accu(parameters(arma::span(Ent, Ent + T - 1))) * parameters(Ent + T + p + 5) + na * nb;
    parameters(Ent + T + p + 1) = R::rnorm(chi * delta, sqrt(delta));
    
    ////// precEtat OK
    auxT = parameters(arma::span(Ent, Ent + T - 1)) - parameters(Ent + T + p + 1);
    parameters(Ent + T + p + 5) = R::rgamma(T / 2 + ga, 1 / ((auxT.t() * auxT).eval()(0) / 2 + gb));
    
    //// precEtas OK
    process = parameters(arma::span(0, nt(0) - 1)) - parameters(Ent);
    sumita  = (process.t() * R.slice(0)(arma::span(0, nt(0) - 1),arma::span(0, nt(0) - 1)) * process).eval()(0,0);  
    for (int tInd = 1; tInd < T; ++tInd) {
      process = parameters(arma::span(CUMSUMnt(tInd), CUMSUMnt(tInd + 1) - 1)) - parameters(Ent + tInd);
      sumita += (process.t() * R.slice(tInd)(arma::span(0, nt(tInd) - 1),arma::span(0, nt(tInd) - 1)) * process).eval()(0,0);  
    }
    parameters(Ent + T + p + 4) = R::rgamma(Ent / 2 + ga, 1 / (sumita / 2 + gb));
    
    if ((b != 0) && (b % nReport == 0)) {
      Rcpp::Rcout << "The value of b : " << b << "\n";
    }
    
    if ((b > 0) && (b % nThin == 0)) {
      keep(b / nThin - 1, arma::span(0, Ent + T + p + 5)) = parameters.t();
      keep(b / nThin - 1, arma::span(Ent + T + p + 6, (Ent + T + p + 6) + (2 * nT2 + 6))) = alpha.t();
      keep(b / nThin - 1, arma::span((Ent + T + p + 6) + (2 * nT2 + 7), (Ent + T + p + 6) + (2 * nT2 + 7) + (2 * nY2 + 6))) = lambda.t();
      keep(b / nThin - 1, arma::span((Ent + T + p + 6) + (2 * nT2 + 7) + (2 * nY2 + 7), (Ent + T + p + 6) + (2 * nT2 + 7) + (2 * nY2 + 7) + (N + nTEMP + T + 6))) = temp.t();
    }
  }
  
  return keep;
}

// Both calibrations with 1 GP in the intercept
// [[Rcpp::export]]
arma::mat GibbsZoop1Cpp(
    arma::vec Y1, arma::vec Y2, arma::vec T1, arma::vec T2,
    arma::mat X, arma::mat Xeta, arma::mat XTEMP, 
    arma::vec Rsum, double RsumTEMP, arma::cube R, arma::mat RTEMP, double dMax,
    arma::mat dY2, arma::mat dT2, double sd,
    arma::vec oneY1, arma::vec oneY2, arma::vec oneT1, arma::vec oneT2,
    arma::uvec indY1, arma::uvec indY2, arma::uvec indT1, arma::uvec indT2,
    arma::mat IY2, arma::mat IT2, arma::mat ITEMP, arma::uvec sY2, arma::uvec sT2, arma::uvec sTEMP, arma::mat Xt,
    int N, int NY1, int NY2, int NT1, int NT2, int T, arma::vec nt, double Ent, int nT2, int nY2, int nTEMP, int p,
    double na, double nb, double ga, double gb,
    int nBurnin, int maxIter, int nThin, int nReport) {
  
  // PARAMETERS //
  arma::vec temp(N + nTEMP + T + 7, arma::fill::ones); // temp'tls, mean's, mean, precmean, precyear, sin, cos, prec, year1 ... yearT, rhoyear
  arma::vec alpha(nT2 + 5, arma::fill::ones); // a0's a1 mean0 prec0 decay0 precT2
  arma::vec lambda(nY2 + 5, arma::fill::ones); // l0's l1 mean0 prec0 decay0 precY2
  arma::vec parameters(Ent + T + p + 6, arma::fill::ones); // eta'ts eta't X' beta1 beta0 precY1 precT1 precEtas precEtat
  arma::mat keep(maxIter / nThin, (Ent + T + p + 6) + (nT2 + 5) + (nY2 + 5) + (N + nTEMP + T + 7));
  
  alpha(nT2 + 3) = 10 / dMax; // decay0
  lambda(nY2 + 3) = 10 / dMax; // decay0
  double Umin = 3 / dMax;
  double Umax = 3 / (0.1 * dMax);

  //parameters(Ent + T + p + 3) = 10000; // precT1
  temp(N + nTEMP + 5) = INFINITY; // precTemp
  temp(N + nTEMP + T + 6) = 0; // rho
  
  arma::vec TRUEtemp = ITEMP * temp(arma::span(N, N + nTEMP - 1)) + XTEMP * temp(arma::span(N + nTEMP + 3, N + nTEMP + 4)) + Xt * temp(arma::span(N + nTEMP + 6, N + nTEMP + T + 5));
  arma::vec ZOOP = Xeta * parameters(arma::span(0, Ent - 1)) + X * parameters(arma::span(Ent + T, Ent + T + p - 1)) + parameters(Ent + T + p) * TRUEtemp;
  
  // HELP //
  arma::mat Delta;
  arma::mat Chi;
  double delta;
  double chi;
  
  arma::vec alphas(nT2);
  arma::vec lambdas(nY2);
  
  arma::vec auxAlpha0(N, arma::fill::zeros); 
  alphas = alpha(arma::span(0, nT2 - 1));
  auxAlpha0.elem(indT2) = alphas.elem(sT2);
  
  arma::vec auxLambda0(N, arma::fill::zeros); 
  lambdas = lambda(arma::span(0, nY2 - 1));
  auxLambda0.elem(indY2) = lambdas.elem(sY2);
  
  arma::vec oneN(N, arma::fill::ones);
  arma::vec auxN(N);
  arma::vec auxT(T);
  arma::vec auxNY1(NY1);
  arma::vec auxNY2(NY2);
  arma::vec auxNT1(NT1);
  arma::vec auxNT2(NT2);
  
  arma::mat r;
  arma::vec process;
  double sumita;
  double contador;
  double rho2;
  
  arma::uvec col0 = { 0 };
  arma::uvec col1 = { 1 };
  
  double sumY1Sin2 = arma::accu(pow(XTEMP(indY1, col0), 2));
  double sumY1Cos2 = arma::accu(pow(XTEMP(indY1, col1), 2));
  double sumY2Sin2 = arma::accu(pow(XTEMP(indY2, col0), 2));
  double sumY2Cos2 = arma::accu(pow(XTEMP(indY2, col1), 2));
  double sumT1Sin2 = arma::accu(pow(XTEMP(indT1, col0), 2));
  double sumT1Cos2 = arma::accu(pow(XTEMP(indT1, col1), 2));
  double sumT2Sin2 = arma::accu(pow(XTEMP(indT2, col0), 2));
  double sumT2Cos2 = arma::accu(pow(XTEMP(indT2, col1), 2));
  
  arma::mat XtauxY1 = Xt.rows(indY1);
  arma::mat XtauxY2 = Xt.rows(indY2);
  arma::mat XtauxT1 = Xt.rows(indT1);
  arma::mat XtauxT2 = Xt.rows(indT2);
  
  arma::vec NY1t = arma::sum(XtauxY1).t();
  arma::vec NY2t = arma::sum(XtauxY2).t();
  arma::vec NT1t = arma::sum(XtauxT1).t();
  arma::vec NT2t = arma::sum(XtauxT2).t();
  
  arma::mat XsauxY1 = ITEMP.rows(indY1);
  arma::mat XsauxY2 = ITEMP.rows(indY2);
  arma::mat XsauxT1 = ITEMP.rows(indT1);
  arma::mat XsauxT2 = ITEMP.rows(indT2);
  
  arma::vec NY1s = arma::sum(XsauxY1).t();
  arma::vec NY2s = arma::sum(XsauxY2).t();
  arma::vec NT1s = arma::sum(XsauxT1).t();
  arma::vec NT2s = arma::sum(XsauxT2).t();
  
  arma::vec initAR1(3);
  initAR1(0) = temp(N + nTEMP);
  initAR1(1) = temp(N + nTEMP + T + 6);
  initAR1(2) = temp(N + nTEMP + 2);
  
  arma::vec CUMSUMnt(T + 1);
  CUMSUMnt(0) = 0;
  CUMSUMnt(arma::span(1, T)) = arma::cumsum(nt);
  
  // ITERATIONS //
  for (int b = 1 - nBurnin; b <= maxIter; ++b) {
    
    //// temp t l s OK
    //ZOOP -= parameters(Ent + T + p) * temp(arma::span(0, N - 1));
    //for (int j = 0; j < N; ++j) {
    //  delta = 1 / (oneT1(j) * parameters(Ent + T + p + 3) +
    //    oneT2(j) * pow(alpha(nT2), 2) * alpha(nT2 + 4) + 
    //    oneY1(j) * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
    //    oneY2(j) * pow(lambda(nY2) * parameters(Ent + T + p), 2) * lambda(nY2 + 4) + 
    //    temp(N + nTEMP + 5));
    //  chi   = T1(j) * oneT1(j) * parameters(Ent + T + p + 3) +
    //    (T2(j) - auxAlpha0(j)) * oneT2(j) * alpha(nT2) * alpha(nT2 + 4) + 
    //    (Y1(j) - ZOOP(j)) * oneY1(j) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
    //    (Y2(j) - (auxLambda0(j) + lambda(nY2) * ZOOP(j))) * oneY2(j) * lambda(nY2) * parameters(Ent + T + p) * lambda(nY2 + 4) +
    //    ((ITEMP.row(j) * temp(arma::span(N, N + nTEMP - 1))).eval()(0) + temp(N + nTEMP + 3) * XTEMP(j,0) + temp(N + nTEMP + 4) * XTEMP(j,1) + (Xt.row(j) * temp(arma::span(N + nTEMP + 6, N + nTEMP + T + 5))).eval()(0)) * temp(N + nTEMP + 5);
    //  temp(j) = R::rnorm(chi * delta, sqrt(delta));
    //}
    //ZOOP += parameters(Ent + T + p) * temp(arma::span(0, N - 1));
    
    //temp(arma::span(N, N + nTEMP + T + 6)) = splmMR(
    //  temp(arma::span(0, N - 1)), XTEMP,
    //  ITEMP, sTEMP, RTEMP, RsumTEMP,
    //  temp(arma::span(N, N + nTEMP + T + 6)), 
    //  na, nb, ga, gb,
    //  2, nTEMP, N, T, Xt,
    //  0, 1, 1, 10
    //).t();
    
    // sin (temp)
    TRUEtemp -= temp(N + nTEMP + 3) * XTEMP.col(0);
    ZOOP     -= parameters(Ent + T + p) * temp(N + nTEMP + 3) * XTEMP.col(0);
    delta = 1 / (
      sumT1Sin2 * parameters(Ent + T + p + 3) +
      sumT2Sin2 * pow(alpha(nT2), 2) * alpha(nT2 + 4) + 
      sumY1Sin2 * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
      sumY2Sin2 * pow(lambda(nY2) * parameters(Ent + T + p), 2) * lambda(nY2 + 4) + 
      nb);
    chi   = 
      arma::accu((T1(indT1) - TRUEtemp(indT1)) % XTEMP(indT1, col0)) * parameters(Ent + T + p + 3) +
      arma::accu((T2(indT2) - (auxAlpha0(indT2) + alpha(nT2) * TRUEtemp(indT2))) % XTEMP(indT2, col0)) * alpha(nT2) * alpha(nT2 + 4) + 
      arma::accu((Y1(indY1) - ZOOP(indY1)) % XTEMP(indY1, col0)) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
      arma::accu((Y2(indY2) - (auxLambda0(indY2) + lambda(nY2) * ZOOP(indY2))) % XTEMP(indY2, col0)) * lambda(nY2) * parameters(Ent + T + p) * lambda(nY2 + 4) +
      na * nb;
    temp(N + nTEMP + 3) = R::rnorm(chi * delta, sqrt(delta));
    TRUEtemp += temp(N + nTEMP + 3) * XTEMP.col(0);
    ZOOP     += parameters(Ent + T + p) * temp(N + nTEMP + 3) * XTEMP.col(0);
    
    // cos (temp)
    TRUEtemp -= temp(N + nTEMP + 4) * XTEMP.col(1);
    ZOOP     -= parameters(Ent + T + p) * temp(N + nTEMP + 4) * XTEMP.col(1);
    delta = 1 / (
      sumT1Cos2 * parameters(Ent + T + p + 3) +
      sumT2Cos2 * pow(alpha(nT2), 2) * alpha(nT2 + 4) + 
      sumY1Cos2 * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
      sumY2Cos2 * pow(lambda(nY2) * parameters(Ent + T + p), 2) * lambda(nY2 + 4) + 
      nb);
    chi   = 
      arma::accu((T1(indT1) - TRUEtemp(indT1)) % XTEMP(indT1, col1)) * parameters(Ent + T + p + 3) +
      arma::accu((T2(indT2) - (auxAlpha0(indT2) + alpha(nT2) * TRUEtemp(indT2))) % XTEMP(indT2, col1)) * alpha(nT2) * alpha(nT2 + 4) + 
      arma::accu((Y1(indY1) - ZOOP(indY1)) % XTEMP(indY1, col1)) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
      arma::accu((Y2(indY2) - (auxLambda0(indY2) + lambda(nY2) * ZOOP(indY2))) % XTEMP(indY2, col1)) * lambda(nY2) * parameters(Ent + T + p) * lambda(nY2 + 4) +
      na * nb;
    temp(N + nTEMP + 4) = R::rnorm(chi * delta, sqrt(delta));
    TRUEtemp += temp(N + nTEMP + 4) * XTEMP.col(1);
    ZOOP     += parameters(Ent + T + p) * temp(N + nTEMP + 4) * XTEMP.col(1);
    
    // year t (temp)
    rho2 = pow(temp(N + nTEMP + T + 6), 2);
    // tInd = 0
    TRUEtemp -= Xt.col(0) * temp(N + nTEMP + 6);
    ZOOP     -= Xt.col(0) * temp(N + nTEMP + 6) * parameters(Ent + T + p);
    delta = 1 / (
      NT1t(0) * parameters(Ent + T + p + 3) +
      NT2t(0) * pow(alpha(nT2), 2) * alpha(nT2 + 4) + 
      NY1t(0) * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
      NY2t(0) * pow(lambda(nY2) * parameters(Ent + T + p), 2) * lambda(nY2 + 4) + 
      (1 + rho2) * temp(N + nTEMP + T + 2));
    chi   = 
      arma::accu((T1(indT1) - TRUEtemp(indT1)) % XtauxT1.col(0)) * parameters(Ent + T + p + 3) +
      arma::accu((T2(indT2) - (auxAlpha0(indT2) + alpha(nT2) * TRUEtemp(indT2))) % XtauxT2.col(0)) * alpha(nT2) * alpha(nT2 + 4) + 
      arma::accu((Y1(indY1) - ZOOP(indY1)) % XtauxY1.col(0)) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
      arma::accu((Y2(indY2) - (auxLambda0(indY2) + lambda(nY2) * ZOOP(indY2))) % XtauxY2.col(0)) * lambda(nY2) * parameters(Ent + T + p) * lambda(nY2 + 4) +
      (temp(N + nTEMP + T + 6) * (temp(N + nTEMP + 7) - temp(N + nTEMP)) + (1 + rho2) * temp(N + nTEMP)) * temp(N + nTEMP + T + 2);
    temp(N + nTEMP + 6) = R::rnorm(chi * delta, sqrt(delta));
    TRUEtemp += Xt.col(0) * temp(N + nTEMP + 6);
    ZOOP     += Xt.col(0) * temp(N + nTEMP + 6) * parameters(Ent + T + p);
    // tInd = 1,...,T-2
    for (int tInd = 1; tInd < T - 1; ++tInd) {
      TRUEtemp -= Xt.col(tInd) * temp(N + nTEMP + 6 + tInd);
      ZOOP     -= Xt.col(tInd) * temp(N + nTEMP + 6 + tInd) * parameters(Ent + T + p);
      delta = 1 / (
        NT1t(tInd) * parameters(Ent + T + p + 3) +
        NT2t(tInd) * pow(alpha(nT2), 2) * alpha(nT2 + 4) + 
        NY1t(tInd) * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
        NY2t(tInd) * pow(lambda(nY2) * parameters(Ent + T + p), 2) * lambda(nY2 + 4) + 
        (1 + rho2) * temp(N + nTEMP + T + 2));
      chi   = 
        arma::accu((T1(indT1) - TRUEtemp(indT1)) % XtauxT1.col(tInd)) * parameters(Ent + T + p + 3) +
        arma::accu((T2(indT2) - (auxAlpha0(indT2) + alpha(nT2) * TRUEtemp(indT2))) % XtauxT2.col(tInd)) * alpha(nT2) * alpha(nT2 + 4) + 
        arma::accu((Y1(indY1) - ZOOP(indY1)) % XtauxY1.col(tInd)) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
        arma::accu((Y2(indY2) - (auxLambda0(indY2) + lambda(nY2) * ZOOP(indY2))) % XtauxY2.col(tInd)) * lambda(nY2) * parameters(Ent + T + p) * lambda(nY2 + 4) +
        (temp(N + nTEMP + T + 6) * (temp(N + nTEMP + 7 + tInd) + temp(N + nTEMP + 5 + tInd)) + pow(1 - temp(N + nTEMP + T + 6), 2) * temp(N + nTEMP)) * temp(N + nTEMP + T + 2);
      temp(N + nTEMP + 6 + tInd) = R::rnorm(chi * delta, sqrt(delta));
      TRUEtemp += Xt.col(tInd) * temp(N + nTEMP + 6 + tInd);
      ZOOP     += Xt.col(tInd) * temp(N + nTEMP + 6 + tInd) * parameters(Ent + T + p);
    }
    // tInd = T-1
    TRUEtemp -= Xt.col(T - 1) * temp(N + nTEMP + 5 + T);
    ZOOP     -= Xt.col(T - 1) * temp(N + nTEMP + 5 + T) * parameters(Ent + T + p);
    delta = 1 / (
      NT1t(T - 1) * parameters(Ent + T + p + 3) +
      NT2t(T - 1) * pow(alpha(nT2), 2) * alpha(nT2 + 4) + 
      NY1t(T - 1) * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
      NY2t(T - 1) * pow(lambda(nY2) * parameters(Ent + T + p), 2) * lambda(nY2 + 4) + 
      temp(N + nTEMP + T + 2));
    chi   = 
      arma::accu((T1(indT1) - TRUEtemp(indT1)) % XtauxT1.col(T - 1)) * parameters(Ent + T + p + 3) +
      arma::accu((T2(indT2) - (auxAlpha0(indT2) + alpha(nT2) * TRUEtemp(indT2))) % XtauxT2.col(T - 1)) * alpha(nT2) * alpha(nT2 + 4) + 
      arma::accu((Y1(indY1) - ZOOP(indY1)) % XtauxY1.col(T - 1)) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
      arma::accu((Y2(indY2) - (auxLambda0(indY2) + lambda(nY2) * ZOOP(indY2))) % XtauxY2.col(T - 1)) * lambda(nY2) * parameters(Ent + T + p) * lambda(nY2 + 4) +
      (temp(N + nTEMP + T + 6) * temp(N + nTEMP + 4 + T) + (1 - temp(N + nTEMP + T + 6)) * temp(N + nTEMP)) * temp(N + nTEMP + T + 2);
    temp(N + nTEMP + 5 + T) = R::rnorm(chi * delta, sqrt(delta));
    TRUEtemp += Xt.col(T - 1) * temp(N + nTEMP + 5 + T);
    ZOOP     += Xt.col(T - 1) * temp(N + nTEMP + 5 + T) * parameters(Ent + T + p);
    
    // beta rho prec (year t in temp)
    initAR1 = AR1(temp(arma::span(N + nTEMP + 6, N + nTEMP + 5 + T)),
                  initAR1, na, nb, ga, gb, T);
    
    temp(N + nTEMP) = initAR1(0);
    temp(N + nTEMP + T + 6) = initAR1(1);
    temp(N + nTEMP + 2) = initAR1(2);
    
    // space si (temp)
    for (int iInd = 0; iInd < nTEMP; ++iInd) {
      TRUEtemp -= ITEMP.col(iInd) * temp(N + iInd);
      ZOOP     -= ITEMP.col(iInd) * temp(N + iInd) * parameters(Ent + T + p);
      process = temp(arma::span(N, N + nTEMP - 1));
      process.shed_row(iInd);
      r = RTEMP.row(iInd);
      r.shed_col(iInd);
      delta = 1 / (
        NT1s(iInd) * parameters(Ent + T + p + 3) +
        NT2s(iInd) * pow(alpha(nT2), 2) * alpha(nT2 + 4) + 
        NY1s(iInd) * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
        NY2s(iInd) * pow(lambda(nY2) * parameters(Ent + T + p), 2) * lambda(nY2 + 4) + 
        RTEMP(iInd, iInd) * temp(N + nTEMP + 1));
      chi   = 
        arma::accu((T1(indT1) - TRUEtemp(indT1)) % XsauxT1.col(iInd)) * parameters(Ent + T + p + 3) +
        arma::accu((T2(indT2) - (auxAlpha0(indT2) + alpha(nT2) * TRUEtemp(indT2))) % XsauxT2.col(iInd)) * alpha(nT2) * alpha(nT2 + 4) + 
        arma::accu((Y1(indY1) - ZOOP(indY1)) % XsauxY1.col(iInd)) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
        arma::accu((Y2(indY2) - (auxLambda0(indY2) + lambda(nY2) * ZOOP(indY2))) % XsauxY2.col(iInd)) * lambda(nY2) * parameters(Ent + T + p) * lambda(nY2 + 4) +
        (- r * process).eval()(0,0) * temp(N + nTEMP + 1);
      temp(N + iInd) = R::rnorm(chi * delta, sqrt(delta));
      TRUEtemp += ITEMP.col(iInd) * temp(N + iInd);
      ZOOP     += ITEMP.col(iInd) * temp(N + iInd) * parameters(Ent + T + p);
    }
    
    // prec (space s in temp)
    process = temp(arma::span(N, N + nTEMP - 1));
    sumita = (process.t() * RTEMP * process).eval()(0,0); 
    temp(N + nTEMP + 1) = R::rgamma(nTEMP / 2 + ga, 1 / (sumita / 2 + gb));
    
    // temp
    temp(arma::span(0, N - 1)) = TRUEtemp;
    
    // alpha's and hiper, precT2
    auxNT2 = temp.elem(indT2);
    
    alpha = splmS(
      T2.elem(indT2), auxNT2, IT2, sT2,
      dT2, alpha, sd,
      na, nb, ga, gb, Umin, Umax,
      nT2, NT2, 0, 1, 1, 10).t();
    
    alphas = alpha(arma::span(0, nT2 - 1));
    auxAlpha0.elem(indT2) = alphas.elem(sT2);
    
    // precT1 
    auxN = T1 - temp(arma::span(0, N - 1)); 
    auxNT1 = auxN.elem(indT1); 
    parameters(Ent + T + p + 3) = R::rgamma(NT1 / 2 + ga, 1 / ((auxNT1.t() * auxNT1).eval()(0) / 2 + gb));
    
    // lambda's and hiper, precY2
    auxNY2 = ZOOP.elem(indY2);
    
    lambda = splmS(
      Y2.elem(indY2), auxNY2, IY2, sY2,
      dY2, lambda, sd,
      na, nb, ga, gb, Umin, Umax,
      nY2, NY2, 0, 1, 1, 10).t();
    
    lambdas = lambda(arma::span(0, nY2 - 1));
    auxLambda0.elem(indY2) = lambdas.elem(sY2);
    
    // precY1 OK
    auxN = Y1 - ZOOP; 
    auxNY1 = auxN.elem(indY1); 
    parameters(Ent + T + p + 2) = R::rgamma(NY1 / 2 + ga, 1 / ((auxNY1.t() * auxNY1).eval()(0) / 2 + gb));
    
    // beta's OK
    for (int j = 0; j < p; ++j) {
      ZOOP -= parameters(Ent + T + j) * X.col(j);
      
      auxN   = pow(X.col(j), 2);
      auxNY1 = auxN.elem(indY1); 
      auxN   = pow(lambda(nY2) * X.col(j), 2);
      auxNY2 = auxN.elem(indY2); 
      
      delta = 1 / (arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(nY2 + 4) + nb);
      
      auxN   = X.col(j) % (Y1 - ZOOP); 
      auxNY1 = auxN.elem(indY1); 
      auxN   = lambda(nY2) * X.col(j) % (Y2 - (auxLambda0 + lambda(nY2) * ZOOP)); 
      auxNY2 = auxN.elem(indY2); 
      
      chi   = arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(nY2 + 4) + na * nb;
      
      parameters(Ent + T + j) = R::rnorm(chi * delta, sqrt(delta));
      
      ZOOP += parameters(Ent + T + j) * X.col(j);
    }
    
    // beta1 temp OK
    ZOOP -= parameters(Ent + T + p) * temp(arma::span(0, N - 1));
    
    auxN   = pow(temp(arma::span(0, N - 1)), 2);
    auxNY1 = auxN.elem(indY1); 
    auxN   = pow(lambda(nY2) * temp(arma::span(0, N - 1)), 2);
    auxNY2 = auxN.elem(indY2); 
    
    delta = 1 / (arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(nY2 + 4) + nb);
    
    auxN   = temp(arma::span(0, N - 1)) % (Y1 - ZOOP); 
    auxNY1 = auxN.elem(indY1); 
    auxN   = lambda(nY2) * temp(arma::span(0, N - 1)) % (Y2 - (auxLambda0 + lambda(nY2) * ZOOP)); 
    auxNY2 = auxN.elem(indY2);
    
    chi   = arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(nY2 + 4) + na * nb;
    
    parameters(Ent + T + p) = R::rnorm(chi * delta, sqrt(delta));
    
    ZOOP += parameters(Ent + T + p) * temp(arma::span(0, N - 1));
    
    // etats
    contador = 0;
    for (int tInd = 0; tInd < T; ++tInd){
      for (int iInd = 0; iInd < nt(tInd); ++iInd){
        
        ZOOP -= Xeta.col(contador) * parameters(contador);
        
        process = parameters(arma::span(CUMSUMnt(tInd), CUMSUMnt(tInd + 1) - 1));
        process.shed_row(iInd);
        r = R(arma::span(iInd, iInd), arma::span(0, nt(tInd) - 1), arma::span(tInd, tInd));
        r.shed_col(iInd);
        
        auxN = Xeta.col(contador);
        auxNY1 = auxN.elem(indY1);
        auxN = pow(lambda(nY2) * auxN, 2);
        auxNY2 = auxN.elem(indY2);
        
        delta = 1 / (arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(nY2 + 4) + R.slice(tInd)(iInd, iInd) * parameters(Ent + T + p + 4));
        
        auxN   = Xeta.col(contador) % (Y1 - ZOOP); 
        auxNY1 = auxN.elem(indY1); 
        auxN   = lambda(nY2) * Xeta.col(contador) % (Y2 - (auxLambda0 + lambda(nY2) * ZOOP)); 
        auxNY2 = auxN.elem(indY2); 
        
        chi   = arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(nY2 + 4) + 
          (parameters(Ent + tInd) * R.slice(tInd)(iInd, iInd) + r * (parameters(Ent + tInd) - process)).eval()(0,0) * parameters(Ent + T + p + 4);
        
        parameters(contador) = R::rnorm(delta * chi, sqrt(delta));
        
        ZOOP += Xeta.col(contador) * parameters(contador);
        
        ++contador;
      }
    }
    
    //// etat OK
    for (int tInd = 0; tInd < T; ++tInd) {
      delta = 1 / (Rsum(tInd) * parameters(Ent + T + p + 4) + parameters(Ent + T + p + 5)); 
      chi   = (oneN(arma::span(0, nt(tInd) - 1)).t() * R.slice(tInd)(arma::span(0, nt(tInd) - 1), arma::span(0, nt(tInd) - 1)) * parameters(arma::span(CUMSUMnt(tInd), CUMSUMnt(tInd + 1) - 1))).eval()(0,0) * parameters(Ent + T + p + 4) + parameters(Ent + T + p + 1) * parameters(Ent + T + p + 5);
      parameters(Ent + tInd) = R::rnorm(chi * delta, sqrt(delta));
    }
    
    ////// beta0 OK
    delta = 1 / (T * parameters(Ent + T + p + 5) + nb);
    chi   = arma::accu(parameters(arma::span(Ent, Ent + T - 1))) * parameters(Ent + T + p + 5) + na * nb;
    parameters(Ent + T + p + 1) = R::rnorm(chi * delta, sqrt(delta));
    
    ////// precEtat OK
    auxT = parameters(arma::span(Ent, Ent + T - 1)) - parameters(Ent + T + p + 1);
    parameters(Ent + T + p + 5) = R::rgamma(T / 2 + ga, 1 / ((auxT.t() * auxT).eval()(0) / 2 + gb));
    
    //// precEtas OK
    process = parameters(arma::span(0, nt(0) - 1)) - parameters(Ent);
    sumita  = (process.t() * R.slice(0)(arma::span(0, nt(0) - 1),arma::span(0, nt(0) - 1)) * process).eval()(0,0);  
    for (int tInd = 1; tInd < T; ++tInd) {
      process = parameters(arma::span(CUMSUMnt(tInd), CUMSUMnt(tInd + 1) - 1)) - parameters(Ent + tInd);
      sumita += (process.t() * R.slice(tInd)(arma::span(0, nt(tInd) - 1),arma::span(0, nt(tInd) - 1)) * process).eval()(0,0);  
    }
    parameters(Ent + T + p + 4) = R::rgamma(Ent / 2 + ga, 1 / (sumita / 2 + gb));
    
    if ((b != 0) && (b % nReport == 0)) {
      Rcpp::Rcout << "The value of b : " << b << "\n";
    }
    
    if ((b > 0) && (b % nThin == 0)) {
      keep(b / nThin - 1, arma::span(0, Ent + T + p + 5)) = parameters.t();
      keep(b / nThin - 1, arma::span(Ent + T + p + 6, (Ent + T + p + 6) + (nT2 + 4))) = alpha.t();
      keep(b / nThin - 1, arma::span((Ent + T + p + 6) + (nT2 + 5), (Ent + T + p + 6) + (nT2 + 5) + (nY2 + 4))) = lambda.t();
      keep(b / nThin - 1, arma::span((Ent + T + p + 6) + (nT2 + 5) + (nY2 + 5), (Ent + T + p + 6) + (nT2 + 5) + (nY2 + 5) + (N + nTEMP + T + 6))) = temp.t();
    }
  }
  
  return keep;
}

// Both calibrations without GPs
// [[Rcpp::export]]
arma::mat GibbsZoop0Cpp(
    arma::vec Y1, arma::vec Y2, arma::vec T1, arma::vec T2,
    arma::mat X, arma::mat Xeta, arma::mat XTEMP, 
    arma::vec Rsum, double RsumTEMP, arma::cube R, arma::mat RTEMP,
    arma::vec oneY1, arma::vec oneY2, arma::vec oneT1, arma::vec oneT2,
    arma::uvec indY1, arma::uvec indY2, arma::uvec indT1, arma::uvec indT2,
    arma::mat ITEMP, arma::uvec sTEMP, arma::mat Xt,
    int N, int NY1, int NY2, int NT1, int NT2, int T, 
    arma::vec nt, double Ent, int nTEMP, int p,
    double na, double nb, double ga, double gb,
    int nBurnin, int maxIter, int nThin, int nReport) {
  
  // PARAMETERS //
  arma::vec temp(N + nTEMP + T + 7, arma::fill::ones); // temp'tls, mean's, mean, precmean, precyear, sin, cos, prec, year1 ... yearT, rhoyear
  arma::vec alpha(3, arma::fill::ones); // a0 a1 precT2
  arma::vec lambda(3, arma::fill::ones); // l0 l1 precY2
  arma::vec parameters(Ent + T + p + 6, arma::fill::ones); // eta'ts eta't X' beta1 beta0 precY1 precT1 precEtas precEtat
  arma::mat keep(maxIter / nThin, (Ent + T + p + 6) + (3) + (3) + (N + nTEMP + T + 7));
  
  //parameters(Ent + T + p + 3) = 10000; // precT1
  temp(N + nTEMP + 5) = INFINITY; // precTemp
  temp(N + nTEMP + T + 6) = 0; // rhoyear

  arma::vec TRUEtemp = ITEMP * temp(arma::span(N, N + nTEMP - 1)) + XTEMP * temp(arma::span(N + nTEMP + 3, N + nTEMP + 4)) + Xt * temp(arma::span(N + nTEMP + 6, N + nTEMP + T + 5));
  arma::vec ZOOP = Xeta * parameters(arma::span(0, Ent - 1)) + X * parameters(arma::span(Ent + T, Ent + T + p - 1)) + parameters(Ent + T + p) * TRUEtemp;

  // HELP //
  arma::mat Delta;
  arma::mat Chi;
  double delta;
  double chi;
  
  arma::vec oneN(N, arma::fill::ones);
  arma::vec auxN(N);
  arma::vec auxT(T);
  arma::vec auxNY1(NY1);
  arma::vec auxNY2(NY2);
  arma::vec auxNT1(NT1);
    
  arma::mat r;
  arma::vec process;
  double sumita;
  double contador;
  double rho2;
  
  arma::uvec col0 = { 0 };
  arma::uvec col1 = { 1 };

  double sumY1Sin2 = arma::accu(pow(XTEMP(indY1, col0), 2));
  double sumY1Cos2 = arma::accu(pow(XTEMP(indY1, col1), 2));
  double sumY2Sin2 = arma::accu(pow(XTEMP(indY2, col0), 2));
  double sumY2Cos2 = arma::accu(pow(XTEMP(indY2, col1), 2));
  double sumT1Sin2 = arma::accu(pow(XTEMP(indT1, col0), 2));
  double sumT1Cos2 = arma::accu(pow(XTEMP(indT1, col1), 2));
  double sumT2Sin2 = arma::accu(pow(XTEMP(indT2, col0), 2));
  double sumT2Cos2 = arma::accu(pow(XTEMP(indT2, col1), 2));

  arma::mat XtauxY1 = Xt.rows(indY1);
  arma::mat XtauxY2 = Xt.rows(indY2);
  arma::mat XtauxT1 = Xt.rows(indT1);
  arma::mat XtauxT2 = Xt.rows(indT2);
  
  arma::vec NY1t = arma::sum(XtauxY1).t();
  arma::vec NY2t = arma::sum(XtauxY2).t();
  arma::vec NT1t = arma::sum(XtauxT1).t();
  arma::vec NT2t = arma::sum(XtauxT2).t();

  arma::mat XsauxY1 = ITEMP.rows(indY1);
  arma::mat XsauxY2 = ITEMP.rows(indY2);
  arma::mat XsauxT1 = ITEMP.rows(indT1);
  arma::mat XsauxT2 = ITEMP.rows(indT2);
  
  arma::vec NY1s = arma::sum(XsauxY1).t();
  arma::vec NY2s = arma::sum(XsauxY2).t();
  arma::vec NT1s = arma::sum(XsauxT1).t();
  arma::vec NT2s = arma::sum(XsauxT2).t();
  
  arma::vec initAR1(3);
  initAR1(0) = temp(N + nTEMP);
  initAR1(1) = temp(N + nTEMP + T + 6);
  initAR1(2) = temp(N + nTEMP + 2);
  
  arma::vec CUMSUMnt(T + 1);
  CUMSUMnt(0) = 0;
  CUMSUMnt(arma::span(1, T)) = arma::cumsum(nt);

  // ITERATIONS //
  for (int b = 1 - nBurnin; b <= maxIter; ++b) {
    
    // temp t l s OK
    //ZOOP -= parameters(Ent + T + p) * temp(arma::span(0, N - 1));
    //for (int j = 0; j < N; ++j) {
    //  delta = 1 / (oneT1(j) * parameters(Ent + T + p + 3) +
    //    oneT2(j) * pow(alpha(1), 2) * alpha(2) + 
    //    oneY1(j) * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
    //    oneY2(j) * pow(lambda(1) * parameters(Ent + T + p), 2) * lambda(2) + 
    //    temp(N + nTEMP + 5));
    //  chi   = T1(j) * oneT1(j) * parameters(Ent + T + p + 3) +
    //    (T2(j) - alpha(0)) * oneT2(j) * alpha(1) * alpha(2) + 
    //    (Y1(j) - ZOOP(j)) * oneY1(j) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
    //    (Y2(j) - (lambda(0) + lambda(1) * ZOOP(j))) * oneY2(j) * lambda(1) * parameters(Ent + T + p) * lambda(2) +
    //    ((ITEMP.row(j) * temp(arma::span(N, N + nTEMP - 1))).eval()(0) + temp(N + nTEMP + 3) * XTEMP(j,0) + temp(N + nTEMP + 4) * XTEMP(j,1) + (Xt.row(j) * temp(arma::span(N + nTEMP + 6, N + nTEMP + T + 5))).eval()(0)) * temp(N + nTEMP + 5);
    //  temp(j) = R::rnorm(chi * delta, sqrt(delta));
    //}
    //ZOOP += parameters(Ent + T + p) * temp(arma::span(0, N - 1));
    
    //temp(arma::span(N, N + nTEMP + T + 6)) = splmMR(
    //  temp(arma::span(0, N - 1)), XTEMP,
    //  ITEMP, sTEMP, RTEMP, RsumTEMP,
    //  temp(arma::span(N, N + nTEMP + T + 6)), 
    //  na, nb, ga, gb,
    //  2, nTEMP, N, T, Xt,
    //  0, 1, 1, 10
    //).t();
    
    // sin (temp)
    TRUEtemp -= temp(N + nTEMP + 3) * XTEMP.col(0);
    ZOOP     -= parameters(Ent + T + p) * temp(N + nTEMP + 3) * XTEMP.col(0);
    delta = 1 / (
      sumT1Sin2 * parameters(Ent + T + p + 3) +
      sumT2Sin2 * pow(alpha(1), 2) * alpha(2) + 
      sumY1Sin2 * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
      sumY2Sin2 * pow(lambda(1) * parameters(Ent + T + p), 2) * lambda(2) + 
      nb);
    chi   = 
      arma::accu((T1(indT1) - TRUEtemp(indT1)) % XTEMP(indT1, col0)) * parameters(Ent + T + p + 3) +
      arma::accu((T2(indT2) - (alpha(0) + alpha(1) * TRUEtemp(indT2))) % XTEMP(indT2, col0)) * alpha(1) * alpha(2) + 
      arma::accu((Y1(indY1) - ZOOP(indY1)) % XTEMP(indY1, col0)) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
      arma::accu((Y2(indY2) - (lambda(0) + lambda(1) * ZOOP(indY2))) % XTEMP(indY2, col0)) * lambda(1) * parameters(Ent + T + p) * lambda(2) +
      na * nb;
    temp(N + nTEMP + 3) = R::rnorm(chi * delta, sqrt(delta));
    TRUEtemp += temp(N + nTEMP + 3) * XTEMP.col(0);
    ZOOP     += parameters(Ent + T + p) * temp(N + nTEMP + 3) * XTEMP.col(0);
    
    // cos (temp)
    TRUEtemp -= temp(N + nTEMP + 4) * XTEMP.col(1);
    ZOOP     -= parameters(Ent + T + p) * temp(N + nTEMP + 4) * XTEMP.col(1);
    delta = 1 / (
      sumT1Cos2 * parameters(Ent + T + p + 3) +
      sumT2Cos2 * pow(alpha(1), 2) * alpha(2) + 
      sumY1Cos2 * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
      sumY2Cos2 * pow(lambda(1) * parameters(Ent + T + p), 2) * lambda(2) + 
      nb);
    chi   = 
      arma::accu((T1(indT1) - TRUEtemp(indT1)) % XTEMP(indT1, col1)) * parameters(Ent + T + p + 3) +
      arma::accu((T2(indT2) - (alpha(0) + alpha(1) * TRUEtemp(indT2))) % XTEMP(indT2, col1)) * alpha(1) * alpha(2) + 
      arma::accu((Y1(indY1) - ZOOP(indY1)) % XTEMP(indY1, col1)) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
      arma::accu((Y2(indY2) - (lambda(0) + lambda(1) * ZOOP(indY2))) % XTEMP(indY2, col1)) * lambda(1) * parameters(Ent + T + p) * lambda(2) +
      na * nb;
    temp(N + nTEMP + 4) = R::rnorm(chi * delta, sqrt(delta));
    TRUEtemp += temp(N + nTEMP + 4) * XTEMP.col(1);
    ZOOP     += parameters(Ent + T + p) * temp(N + nTEMP + 4) * XTEMP.col(1);

    // year t (temp)
    rho2 = pow(temp(N + nTEMP + T + 6), 2);
    // tInd = 0
    TRUEtemp -= Xt.col(0) * temp(N + nTEMP + 6);
    ZOOP     -= Xt.col(0) * temp(N + nTEMP + 6) * parameters(Ent + T + p);
    delta = 1 / (
      NT1t(0) * parameters(Ent + T + p + 3) +
      NT2t(0) * pow(alpha(1), 2) * alpha(2) + 
      NY1t(0) * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
      NY2t(0) * pow(lambda(1) * parameters(Ent + T + p), 2) * lambda(2) + 
      (1 + rho2) * temp(N + nTEMP + T + 2));
    chi   = 
      arma::accu((T1(indT1) - TRUEtemp(indT1)) % XtauxT1.col(0)) * parameters(Ent + T + p + 3) +
      arma::accu((T2(indT2) - (alpha(0) + alpha(1) * TRUEtemp(indT2))) % XtauxT2.col(0)) * alpha(1) * alpha(2) + 
      arma::accu((Y1(indY1) - ZOOP(indY1)) % XtauxY1.col(0)) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
      arma::accu((Y2(indY2) - (lambda(0) + lambda(1) * ZOOP(indY2))) % XtauxY2.col(0)) * lambda(1) * parameters(Ent + T + p) * lambda(2) +
      (temp(N + nTEMP + T + 6) * (temp(N + nTEMP + 7) - temp(N + nTEMP)) + (1 + rho2) * temp(N + nTEMP)) * temp(N + nTEMP + T + 2);
    temp(N + nTEMP + 6) = R::rnorm(chi * delta, sqrt(delta));
    TRUEtemp += Xt.col(0) * temp(N + nTEMP + 6);
    ZOOP     += Xt.col(0) * temp(N + nTEMP + 6) * parameters(Ent + T + p);
    // tInd = 1,...,T-2
    for (int tInd = 1; tInd < T - 1; ++tInd) {
      TRUEtemp -= Xt.col(tInd) * temp(N + nTEMP + 6 + tInd);
      ZOOP     -= Xt.col(tInd) * temp(N + nTEMP + 6 + tInd) * parameters(Ent + T + p);
      delta = 1 / (
        NT1t(tInd) * parameters(Ent + T + p + 3) +
        NT2t(tInd) * pow(alpha(1), 2) * alpha(2) + 
        NY1t(tInd) * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
        NY2t(tInd) * pow(lambda(1) * parameters(Ent + T + p), 2) * lambda(2) + 
        (1 + rho2) * temp(N + nTEMP + T + 2));
      chi   = 
        arma::accu((T1(indT1) - TRUEtemp(indT1)) % XtauxT1.col(tInd)) * parameters(Ent + T + p + 3) +
        arma::accu((T2(indT2) - (alpha(0) + alpha(1) * TRUEtemp(indT2))) % XtauxT2.col(tInd)) * alpha(1) * alpha(2) + 
        arma::accu((Y1(indY1) - ZOOP(indY1)) % XtauxY1.col(tInd)) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
        arma::accu((Y2(indY2) - (lambda(0) + lambda(1) * ZOOP(indY2))) % XtauxY2.col(tInd)) * lambda(1) * parameters(Ent + T + p) * lambda(2) +
        (temp(N + nTEMP + T + 6) * (temp(N + nTEMP + 7 + tInd) + temp(N + nTEMP + 5 + tInd)) + pow(1 - temp(N + nTEMP + T + 6), 2) * temp(N + nTEMP)) * temp(N + nTEMP + T + 2);
      temp(N + nTEMP + 6 + tInd) = R::rnorm(chi * delta, sqrt(delta));
      TRUEtemp += Xt.col(tInd) * temp(N + nTEMP + 6 + tInd);
      ZOOP     += Xt.col(tInd) * temp(N + nTEMP + 6 + tInd) * parameters(Ent + T + p);
    }
    // tInd = T-1
    TRUEtemp -= Xt.col(T - 1) * temp(N + nTEMP + 5 + T);
    ZOOP     -= Xt.col(T - 1) * temp(N + nTEMP + 5 + T) * parameters(Ent + T + p);
    delta = 1 / (
      NT1t(T - 1) * parameters(Ent + T + p + 3) +
      NT2t(T - 1) * pow(alpha(1), 2) * alpha(2) + 
      NY1t(T - 1) * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
      NY2t(T - 1) * pow(lambda(1) * parameters(Ent + T + p), 2) * lambda(2) + 
      temp(N + nTEMP + T + 2));
    chi   = 
      arma::accu((T1(indT1) - TRUEtemp(indT1)) % XtauxT1.col(T - 1)) * parameters(Ent + T + p + 3) +
      arma::accu((T2(indT2) - (alpha(0) + alpha(1) * TRUEtemp(indT2))) % XtauxT2.col(T - 1)) * alpha(1) * alpha(2) + 
      arma::accu((Y1(indY1) - ZOOP(indY1)) % XtauxY1.col(T - 1)) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
      arma::accu((Y2(indY2) - (lambda(0) + lambda(1) * ZOOP(indY2))) % XtauxY2.col(T - 1)) * lambda(1) * parameters(Ent + T + p) * lambda(2) +
      (temp(N + nTEMP + T + 6) * temp(N + nTEMP + 4 + T) + (1 - temp(N + nTEMP + T + 6)) * temp(N + nTEMP)) * temp(N + nTEMP + T + 2);
    temp(N + nTEMP + 5 + T) = R::rnorm(chi * delta, sqrt(delta));
    TRUEtemp += Xt.col(T - 1) * temp(N + nTEMP + 5 + T);
    ZOOP     += Xt.col(T - 1) * temp(N + nTEMP + 5 + T) * parameters(Ent + T + p);
    
    // beta rho prec (year t in temp)
    initAR1 = AR1(temp(arma::span(N + nTEMP + 6, N + nTEMP + 5 + T)),
      initAR1, na, nb, ga, gb, T);
    
    temp(N + nTEMP) = initAR1(0);
    temp(N + nTEMP + T + 6) = initAR1(1);
    temp(N + nTEMP + 2) = initAR1(2);
    
    // space si (temp)
    for (int iInd = 0; iInd < nTEMP; ++iInd) {
      TRUEtemp -= ITEMP.col(iInd) * temp(N + iInd);
      ZOOP     -= ITEMP.col(iInd) * temp(N + iInd) * parameters(Ent + T + p);
      process = temp(arma::span(N, N + nTEMP - 1));
      process.shed_row(iInd);
      r = RTEMP.row(iInd);
      r.shed_col(iInd);
      delta = 1 / (
        NT1s(iInd) * parameters(Ent + T + p + 3) +
        NT2s(iInd) * pow(alpha(1), 2) * alpha(2) + 
        NY1s(iInd) * pow(parameters(Ent + T + p), 2) * parameters(Ent + T + p + 2) +
        NY2s(iInd) * pow(lambda(1) * parameters(Ent + T + p), 2) * lambda(2) + 
        RTEMP(iInd, iInd) * temp(N + nTEMP + 1));
      chi   = 
        arma::accu((T1(indT1) - TRUEtemp(indT1)) % XsauxT1.col(iInd)) * parameters(Ent + T + p + 3) +
        arma::accu((T2(indT2) - (alpha(0) + alpha(1) * TRUEtemp(indT2))) % XsauxT2.col(iInd)) * alpha(1) * alpha(2) + 
        arma::accu((Y1(indY1) - ZOOP(indY1)) % XsauxY1.col(iInd)) * parameters(Ent + T + p) * parameters(Ent + T + p + 2) +
        arma::accu((Y2(indY2) - (lambda(0) + lambda(1) * ZOOP(indY2))) % XsauxY2.col(iInd)) * lambda(1) * parameters(Ent + T + p) * lambda(2) +
        (- r * process).eval()(0,0) * temp(N + nTEMP + 1);
      temp(N + iInd) = R::rnorm(chi * delta, sqrt(delta));
      TRUEtemp += ITEMP.col(iInd) * temp(N + iInd);
      ZOOP     += ITEMP.col(iInd) * temp(N + iInd) * parameters(Ent + T + p);
    }

    // prec (space s in temp)
    process = temp(arma::span(N, N + nTEMP - 1));
    sumita = (process.t() * RTEMP * process).eval()(0,0); 
    temp(N + nTEMP + 1) = R::rgamma(nTEMP / 2 + ga, 1 / (sumita / 2 + gb));
    
    // temp
    temp(arma::span(0, N - 1)) = TRUEtemp;
    
    // alpha0 alpha1 precT2
    alpha = lmS(
      T2.elem(indT2), temp.elem(indT2), alpha,
      na, nb, ga, gb,
      NT2, 0, 1, 1, 10).t();
    
    // precT1 
    auxN = T1 - temp(arma::span(0, N - 1)); 
    auxNT1 = auxN.elem(indT1); 
    parameters(Ent + T + p + 3) = R::rgamma(NT1 / 2 + ga, 1 / ((auxNT1.t() * auxNT1).eval()(0) / 2 + gb));
    
    // lambda0 lambda1 precY2
    lambda = lmS(
      Y2.elem(indY2), ZOOP.elem(indY2), lambda,
      na, nb, ga, gb,
      NY2, 0, 1, 1, 10).t();
    
    // precY1 OK
    auxN = Y1 - ZOOP; 
    auxNY1 = auxN.elem(indY1); 
    parameters(Ent + T + p + 2) = R::rgamma(NY1 / 2 + ga, 1 / ((auxNY1.t() * auxNY1).eval()(0) / 2 + gb));
    
    // beta's OK
    for (int j = 0; j < p; ++j) {
      ZOOP -= parameters(Ent + T + j) * X.col(j);
      
      auxN   = pow(X.col(j), 2);
      auxNY1 = auxN.elem(indY1); 
      auxN   = pow(lambda(1) * X.col(j), 2);
      auxNY2 = auxN.elem(indY2); 
      
      delta = 1 / (arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(2) + nb);
      
      auxN   = X.col(j) % (Y1 - ZOOP); 
      auxNY1 = auxN.elem(indY1); 
      auxN   = lambda(1) * X.col(j) % (Y2 - (lambda(0) + lambda(1) * ZOOP)); 
      auxNY2 = auxN.elem(indY2); 
      
      chi   = arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(2) + na * nb;
      
      parameters(Ent + T + j) = R::rnorm(chi * delta, sqrt(delta));
      
      ZOOP += parameters(Ent + T + j) * X.col(j);
    }
    
    // beta1 temp OK
    ZOOP -= parameters(Ent + T + p) * temp(arma::span(0, N - 1));
    
    auxN   = pow(temp(arma::span(0, N - 1)), 2);
    auxNY1 = auxN.elem(indY1); 
    auxN   = pow(lambda(1) * temp(arma::span(0, N - 1)), 2);
    auxNY2 = auxN.elem(indY2); 
    
    delta = 1 / (arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(2) + nb);
    
    auxN   = temp(arma::span(0, N - 1)) % (Y1 - ZOOP); 
    auxNY1 = auxN.elem(indY1); 
    auxN   = lambda(1) * temp(arma::span(0, N - 1)) % (Y2 - (lambda(0) + lambda(1) * ZOOP)); 
    auxNY2 = auxN.elem(indY2);
    
    chi   = arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(2) + na * nb;
    
    parameters(Ent + T + p) = R::rnorm(chi * delta, sqrt(delta));
    
    ZOOP += parameters(Ent + T + p) * temp(arma::span(0, N - 1));
    
    // etats
    contador = 0;
    for (int tInd = 0; tInd < T; ++tInd){
      for (int iInd = 0; iInd < nt(tInd); ++iInd){
        
        ZOOP -= Xeta.col(contador) * parameters(contador);
        
        process = parameters(arma::span(CUMSUMnt(tInd), CUMSUMnt(tInd + 1) - 1));
        process.shed_row(iInd);
        r = R(arma::span(iInd, iInd), arma::span(0, nt(tInd) - 1), arma::span(tInd, tInd));
        r.shed_col(iInd);
        
        auxN = Xeta.col(contador);
        auxNY1 = auxN.elem(indY1);
        auxN = pow(lambda(1) * auxN, 2);
        auxNY2 = auxN.elem(indY2);
        
        delta = 1 / (arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(2) + R.slice(tInd)(iInd, iInd) * parameters(Ent + T + p + 4));
        
        auxN   = Xeta.col(contador) % (Y1 - ZOOP); 
        auxNY1 = auxN.elem(indY1); 
        auxN   = lambda(1) * Xeta.col(contador) % (Y2 - (lambda(0) + lambda(1) * ZOOP)); 
        auxNY2 = auxN.elem(indY2); 
        
        chi   = arma::accu(auxNY1) * parameters(Ent + T + p + 2) + arma::accu(auxNY2) * lambda(2) + 
          (parameters(Ent + tInd) * R.slice(tInd)(iInd, iInd) + r * (parameters(Ent + tInd) - process)).eval()(0,0) * parameters(Ent + T + p + 4);
        
        parameters(contador) = R::rnorm(delta * chi, sqrt(delta));
        
        ZOOP += Xeta.col(contador) * parameters(contador);
        
        ++contador;
      }
    }
    
    //// etat OK
    for (int tInd = 0; tInd < T; ++tInd) {
      delta = 1 / (Rsum(tInd) * parameters(Ent + T + p + 4) + parameters(Ent + T + p + 5)); 
      chi   = (oneN(arma::span(0, nt(tInd) - 1)).t() * R.slice(tInd)(arma::span(0, nt(tInd) - 1), arma::span(0, nt(tInd) - 1)) * parameters(arma::span(CUMSUMnt(tInd), CUMSUMnt(tInd + 1) - 1))).eval()(0,0) * parameters(Ent + T + p + 4) + parameters(Ent + T + p + 1) * parameters(Ent + T + p + 5);
      parameters(Ent + tInd) = R::rnorm(chi * delta, sqrt(delta));
    }
    
    ////// beta0 OK
    delta = 1 / (T * parameters(Ent + T + p + 5) + nb);
    chi   = arma::accu(parameters(arma::span(Ent, Ent + T - 1))) * parameters(Ent + T + p + 5) + na * nb;
    parameters(Ent + T + p + 1) = R::rnorm(chi * delta, sqrt(delta));
    
    ////// precEtat OK
    auxT = parameters(arma::span(Ent, Ent + T - 1)) - parameters(Ent + T + p + 1);
    parameters(Ent + T + p + 5) = R::rgamma(T / 2 + ga, 1 / ((auxT.t() * auxT).eval()(0) / 2 + gb));
    
    //// precEtas OK
    process = parameters(arma::span(0, nt(0) - 1)) - parameters(Ent);
    sumita  = (process.t() * R.slice(0)(arma::span(0, nt(0) - 1),arma::span(0, nt(0) - 1)) * process).eval()(0,0);  
    for (int tInd = 1; tInd < T; ++tInd) {
      process = parameters(arma::span(CUMSUMnt(tInd), CUMSUMnt(tInd + 1) - 1)) - parameters(Ent + tInd);
      sumita += (process.t() * R.slice(tInd)(arma::span(0, nt(tInd) - 1),arma::span(0, nt(tInd) - 1)) * process).eval()(0,0);  
    }
    parameters(Ent + T + p + 4) = R::rgamma(Ent / 2 + ga, 1 / (sumita / 2 + gb));
    
    if ((b != 0) && (b % nReport == 0)) {
      Rcpp::Rcout << "The value of b : " << b << "\n";
    }
    
    if ((b > 0) && (b % nThin == 0)) {
      keep(b / nThin - 1, arma::span(0, Ent + T + p + 5)) = parameters.t();
      keep(b / nThin - 1, arma::span(Ent + T + p + 6, (Ent + T + p + 6) + (2))) = alpha.t();
      keep(b / nThin - 1, arma::span((Ent + T + p + 6) + (3), (Ent + T + p + 6) + (3) + (2))) = lambda.t();
      keep(b / nThin - 1, arma::span((Ent + T + p + 6) + (3) + (3), (Ent + T + p + 6) + (3) + (3) + (N + nTEMP + T + 6))) = temp.t();
    }
  }
  
  return keep;
}



