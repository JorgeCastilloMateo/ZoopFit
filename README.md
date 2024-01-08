# ZoopFit
The R package *ZoopFit* is the companion package for the paper Castillo-Mateo et al. (2023) <doi:10.1007/s10651-023-00583-6>. The package includes a data snippet and functions to fit and predict with the double data fusion and calibration models from the paper. 

> **Warning:** Please note that the code provided is for the paper only and a generalization is under development. I do not guarantee that it will work properly outside of the data from the paper.

## Installation
You can install the **development** version from
[GitHub](https://github.com/JorgeCastilloMateo/ZoopFit)

```s
if (!require("remotes")) install.packages("remotes")
remotes::install_github("JorgeCastilloMateo/ZoopFit")
```

## Workflow of the paper
### Simulation Study
To reproduce the simulation study in Section 3 of the paper, see this file: `inst/scripts/00simulation.R`

### Data for Fitting Zooplankton Model
Once the `ZoopFit` package is installed, we progress to the data. 

> Note that the data were collected by and are property of the Center for Coastal Studies. We share a small 200 row subset of these data here to show how the model fitting works. To fully reproduce this analysis, you will need to contact CCS to get access via a data sharing agreement. 

To access the data, simply load in the R data object. 

```s
data("zoop_snippet")
```

This object contains 11 components:

1. `coords` - a 2-column matrix that contains the x and y locations of the sampling stations. These are in projected units. (See below for projection details.)
2. `X` - a three column matrix of covariates: 1) whale presence/absence (0/1); 2) regime 2 (0/1); regime 3  (0/1).
3. `date` - a vector of months of the study period: `"2003-01-06 UTC" "2019-06-09 UTC"`
4. `oneT1` - a binary vector indicating if the CTD measurement of surface temperature was made
5. `oneT2` - a binary vector indicating if the thermistor measurement of surface temperature was made
6. `oneY1` - a binary vector indicating if the surface tow of zooplankton was made
7. `oneY2` - a binary vector indicating if the oblique tow of zooplankton was made
8. `T1` - the measured CTD temperature value (in degrees C)
9. `T1` - the measured thermistor temperature value (in degrees C)
10. `Y1` - the measured total zooplankton recorded (log organisms per m^3) by the surface tow
11. `Y2` - the measured total zooplankton recorded (log organisms per m^3) by the oblique tow

The indicator variables are comprised of 0s and 1s and are all the same length as all the other vectors, whereas the functions require the measured vectors vary in length according to their sampling status, i.e. `sum(oneT1) == length(T1)`. This can be done with `T1 <- T1[oneT1 == 1]`.

#### Projection Details
We use an Albers Equal Area projection

```s
ccb_crs <- paste("+proj=aea +lat_1=40",
                 "+lat_2=45 +lat_0=42 +lon_0=-70",
                 "+x_0=0 +y_0=0",
                 "+ellps=WGS84 +datum=WGS84 +units=m +no_defs")
```

### Data for Fitting the Whale Model
To fit the zooplankton model one of the covariates above `X[, 1]` is a binary vector of whale presence and absence.

> As with the zooplankton data, the right whale presence/absence data are collected by and property of the Center for Coastal Studies. Please contact CCS to get access via an agreement. 

We will need an NARW p/a surface for prediction, and the way that is done is via a temporally dynamic model that is fit separately from the zooplankton model. We detail this below after describing the workflow of the zooplankton model.

To see an example with simulated data fitting the whale model, 
see this file: `inst/scripts/02modelWhales.R`

### Zooplankton Model Fitting
The workhorse function is `GibbsZoop.R`, and will return a matrix with posterior samples from model parameters: rows are iterations and columns are parameters. Here we show how to call the function in order to fit two MCMC chains. We feed all 11 data objects into the function, we choose the model (more below), and we specify MCMC-specific information about the burn-in, the length of the chain, thinning, and how frequently we report the progress. 

```s
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
```

The default model `"corTime"` is the full model used in the manuscript, but the function can flexibly accommodate 5 different formulations: 1) a linear model (`lm`) in the calibration, 2) a lm with a spatially varying intercept `GP`, 3) a lm with spatially varying intercept and slope `GPs`, 4) the same with coregionalization `cor`, and 5) the same with time varying coefs `corTime`. 

After running `GibbsZoop` then we check for scale reduction, i.e., Gelman and Rubin's convergence diagnostic with the `coda` package.

```s
library("coda")
R <- rep(NA, ncol(MODEL[[1]]))
for (i in 1:ncol(MODEL[[1]])) {
  print(i)
  R[i] <- try(gelman.diag(mcmc.list(mcmc(MODEL[[1]][,i]), mcmc(MODEL[[2]][,i])))$psrf[1])
}
max(R[!is.nan(R)])

MODEL <- rbind(MODEL[[1]], MODEL[[2]])
```

We then perform Bayesian Kriging for the Gaussian Processes. This is done with the `GP` function, and the goal is to use coefficients from the locations where we fit the Gaussian Processes (recall the `coords` in the data). We used these parameters to predict the surfaces across a new grid - in this case the entire CCB. (See `inst/scripts/01grid.R` for details of the grid within CCB.)

We note the contents of the function:

```s
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
```

...as well as how its called. Note that since we have fitted yearly GPs, these are called in a loop over year:

```s
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
```

### NARW Model Fitting
This is done in `jags`; see the `02modelWhales.R` script for fitting that model, which needs to be fit prior to prediction. 

You will need to have a current installation of `jags` on your computer (see the `jags` [Sourceforge page](https://mcmc-jags.sourceforge.io) for platform-specific installation information), as well as two R libraries to call `jags`:

```s
library("rjags")
library("jagsUI")
```

Define the model in the `jags` language:

```s
model.JAGS <-"
  model{
    for (i in 1:N) {
      Y[i] ~ dbern(p[i])
      logit(p[i]) <- X[i,] %*% beta + psi[year[i]]
    }
    for (t in 1:T) {
      psi[t] ~ dnorm(beta0, tau_psi)
    }
    beta0 ~ dnorm(0,1.0E-04)
    for (j in 1:k) {
      beta[j] ~ dnorm(0,1.0E-04)
    }
    tau_psi ~ dgamma(2, 1)
 }"
```

We fit the model and then make daily posterior predictions (see the script mentioned to see the arguments used).

```s
model <- autojags(
  data = list('Y' = Y, 'N' = N, 'X' = X, 'year' = year, 'T' = 17, 'k' = 2),
  parameters.to.save = c("beta0", "beta", "psi", "tau_psi"),
  model.file = textConnection(model.JAGS),
  save.all.iter = TRUE,
  n.chains = 2, max.iter = 10000,
  n.burnin = 0, n.thin = 1, Rhat.limit = 1.1)

model <- update(model, parameters.to.save = c("beta0", "beta", "psi", "tau_psi"), 
  n.iter = 10000, n.thin = 10)
```

When we make predictions, we create two objects that are used in the plotting of the zooplankton predictions. These objects are `probWhales` and `predictWhales`. The dimensions of each are the number of stored MCMC iterations (rows) by the prediction days (181 days across 17 years) (columns).

```s
m <- 17 * 181 # number of new data
probWhales    <- matrix(nrow = 2000, ncol = m)
predictWhales <- matrix(nrow = 2000, ncol = m)
p <- rep(NA, m)
d <- rep(1:181, 17)
Xhat <- cbind(sin(2 * pi * d / 365), cos(2 * pi * d / 365))
year <- rep(1:17, each = 181)
for (b in 1:2000) {
  p <- 1 / (1 + exp(-(Xhat %*% model$sims.list$beta[b,] + model$sims.list$psi[b,year])))
  probWhales[b,]    <- p
  predictWhales[b,] <- rbinom(m, size = 1, prob = p)
}
```

Now we have fitted the two models - one for zooplankton and one for the probability of whale presence. Next we turn to prediction and plotting.


### Plotting Results
To reproduce the maps that are presented in the paper, we use the `predictMAPZOOP` function, whose parameters are the fitted `model`, a `date` value, the probability of whale presence (`predictWhales`), which is generated separately from the zooplankton model using the `02modelWhales.R` script (see previous section), the output from the kriged GPs, i.e., the spatial fields for temperature (`gamma0s`) and zooplankton (`etats`), respectively. Three additional arguments are for the spatial grid, a background color, and the spatial extent (these are hard-coded to the extent (in projected units) of Cape Cod Bay). 

The function then predicts and plots the mean zoop, the standard deviation, the exceedance probabilities over two thresholds (1,000 and 1,500 organisms/m^3). The logic is to call a separate function that returns the zooplankton values (`ZOOP`). We then create a data frame that is plotted with `ggplot2` code.

We call the `predictMAPZOOP` function over a series of three different days that are emblematic of the three zooplankton regimes discussed in the paper, e.g.,

```s
predictMAPZOOP(MODEL, "2011-04-29", predictWhales, 
               gamma0s, etats, grid, eeuu, zlim = rbind(c(250, 2650), c(0,1)))
```

We can also examine the yearly average in a similar fashion with a slightly different function (`predictMAPZOOPyear`) whose logic is the same, but accounts for changes within an entire year as opposed to a single day. To make yearly maps:

```s
for (i in 2003:2019) {
  print(i)
  predictMAPZOOPyear(MODEL, i, 1:181, predictWhales, 
                     gamma0s, etats, grid, eeuu)
}
```

To examine the daily time series of zooplankton, we use the `tsZOOP` function. This returns a three element list containing the posterior predictive samples of predicted temperatures, zooplankton and zooplankton exceedance. To summarize the posteriors and make plots, we call the function and then derive the quantiles of the parameter of interest. For example, for zooplankton if we want to extract the daily time series:

```s
my_out <- tsZOOP
Ynew <- my_out$predictZOOP
```

To then visualize it and make the figure from the manuscript:

```s
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
```

Alternatively if we want to mean exceedance, we first extract the different list element and then summarize and plot:

```s
Ynew <- my_out$IndZOOPabove1000
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
```

## How to cite?
To cite *ZoopFit* in publications use:

Castillo-Mateo J, Gelfand AE, Hudak CA, Mayo CA, Schick RS (2023).
“Space-time multi-level modeling for zooplankton abundance employing double data fusion and calibration.” 
*Environmental and Ecological Statistics*, **30**(4), 769--795.
<doi:10.1007/s10651-023-00583-6>.

A BibTeX entry for LaTeX users is

@Article{,  
&nbsp;&nbsp;  title = {Space-time multi-level modeling for zooplankton abundance employing double data fusion and calibration},  
&nbsp;&nbsp;  author = {Jorge Castillo-Mateo and Alan E. Gelfand, Christine A. Hudak, Charles A. Mayo and Robert S. Schick},  
&nbsp;&nbsp;  journal = {Environmental and Ecological Statistics},  
&nbsp;&nbsp;  year = {2023},  
&nbsp;&nbsp;  volume = {30},  
&nbsp;&nbsp;  number = {4},  
&nbsp;&nbsp;  pages = {769--795},  
&nbsp;&nbsp;  doi = {10.1007/s10651-023-00583-6},  
}

  
## References
Castillo-Mateo J, Gelfand AE, Hudak CA, Mayo CA, Schick RS (2023).
“Space-time multi-level modeling for zooplankton abundance employing double data fusion and calibration.” 
*Environmental and Ecological Statistics*, **30**(4), 769--795.
<doi:10.1007/s10651-023-00583-6>.
