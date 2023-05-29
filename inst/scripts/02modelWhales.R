### PAPER EES -- SECTION 5 ###
###      WHALES MODEL      ###
###       USING JAGS       ###

# (This is a code outline to obtain probabilities of whale presence, 
# real whales data are not publicly available, we illustrate with simulation)
set.seed(12345)

# SIMULATED WHALES DATA
N    <- 2000
day  <- sample(181, N, replace = TRUE)
year <- sample(17, N, replace = TRUE)
beta <- c(5.7, 1.25)
X    <- cbind(sin(2 * pi * day / 365), 
              cos(2 * pi * day / 365))
yearEffect <- rnorm(17, -8, 0.84)
Y    <- rbinom(N, size = 1, 
  prob = 1 / (1 + exp(- (X %*% beta + yearEffect[year]))))

# PACKAGES
library("rjags")
library("jagsUI")

# DEFINE THE MODEL
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

#FIT THE MODEL
model <- autojags(
  data = list('Y' = Y, 'N' = N, 'X' = X, 'year' = year, 'T' = 17, 'k' = 2),
  parameters.to.save = c("beta0", "beta", "psi", "tau_psi"),
  model.file = textConnection(model.JAGS),
  save.all.iter = TRUE,
  n.chains = 2, max.iter = 10000,
  n.burnin = 0, n.thin = 1, Rhat.limit = 1.1)

model <- update(model, parameters.to.save = c("beta0", "beta", "psi", "tau_psi"), 
  n.iter = 10000, n.thin = 10)

traceplot(model)

# PREDICT
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
