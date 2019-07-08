### WinBUGS code for the case-control model 
### with contaminated controls fit to Hellbender data
### Rota et al. 2013

# This is WinBUGS model code for the Bayesian model fit to hellbender data.
# b[1] is the intercept parameter, b[2:4] are regression coefficients
model{
  
  # prior distributions for intercept parameter and regression coefficients
  for(j in 1:4){
    b[j] ~  dnorm(0, 0.01)
  }
  # prior distribution for pi
  pi ~ dunif(0,1)
  
  for(i in 1:n){
    eta[i] <- exp(b[1] + log(( n1 / (pi * n0) ) + 1) + b[2] * cover[i] +
                    b[3] * coarse[i] + b[4] * bedrock[i])
    psi[i] <- eta[i] / (1 + eta[i])  # equation 5
    z[i] ~ dbern(psi[i])  # modeling the unobserved 'true' state
    theta[i] <-  n1 / (n1 + pi * n0)  # equation 6
    mu[i] <- theta[i] * z[i] 
    y[i] ~ dbern(mu[i])  # likelihood
  }
}




### Simulating data for the case-control model with contaminated controls

library(boot)  # for the ‘logit’ and ‘inv.logit’ functions
library(R2WinBUGS)  # for data-cloning

b0.true <- 0     # data-generating value for intercept parameter
b1.true <- 3     # data-generating value for continuous covariate
b2.true <- -3    # data-generating value for categorical covariate
N <- 2.8e+6      # total number of resource units
n.sims <- 100    # number of simulations

k <- 10  # number of data 'clones'

# Vectors representing each simulation scenario
# I vary prevalence by changing the mean of the continuous covariate
# pi.scenario corresponds to pi = 0.05, 0.45, 0.75
pi.scenario <- c(-1.59, 0.34, 1.35)
n.scenario <- c(1000, 2000)

# simulating covariates and true occupancy
continuous.cov <- rnorm(n = N, mean = pi.scenario[i])   # cont covariate
categorical.cov <- rbinom(n = N, size = 1, prob = 0.5)  # cat covariate
p.use <- inv.logit(b0.true + b1.true * continuous.cov +
                     b2.true * categorical.cov)  # calculating prob of use
z <- ifelse(test = p.use > runif(n = N), yes = 1, no = 0)  # 1 = used

# simulating sample size
n1 <- n.scenario[j] / 2  # number of used locations sampled
n0 <- n.scenario[j] / 2  # number of available locations sampled

# y = 1 if used, y = 0 if available
y <- c(rep(x = 1, times = n1), rep(x = 0, times = n0))

# randomly sampling indices of used and available locations
used <- sample(x = which(z == 1), size = n1, replace = T)
avail <- sample(x = 1:N, size = n0, replace = T)

# these are the values of covariates at sampled locations
cont <- c(continuous.cov[used], continuous.cov[avail])
cat <- c(categorical.cov[used], categorical.cov[avail])

# data list for export to WinBUGS.  K = 1 clones is a fully Bayesian
# treatment
data <- list(n1 = k * n1,
             n0 = k * n0,
             n = k * (n1 + n0),
             con = rep(x = cont, times = k),
             cat = rep(x = cat, times = k),
             y = rep(x = y, times = k))

# Telling WinBUGS which parameters to monitor
params <- list('b', 'pi')

# Generating random initial values for WinBUGS
inits <- function(){
  list(b = rnorm(n = 3),
       pi = runif(n = 1),
       z = rbinom(n = k * (n1 + n0), size = 1, prob = 0.5))
}

# Fitting the model via R2WinBUGS
fit <- bugs(data = data, inits = inits, parameters.to.save = params,
            model.file = 'Case-Control Simulation.txt', n.chains = 3,
            n.iter = iter, n.burnin = burnin, n.thin = thin, debug = F) 
WinBUGS code for the case-control model with contaminated controls fit to simulated data

model{
  
  # prior distribution for regression coefficients
  for (j in 1:3) {
    b[j] ~  dnorm(0, 0.01)
  }
  # prior distribution for prevalence
  pi ~ dunif(0, 1)
  
  # case control adjustment, eqn. 4
  bstar0 <- log(((n1) / (pi * n0)) + 1) + b[1]
  
  for(i in 1:n){
    z[i] ~ dbern(psi[i])  # process model, eqn. 5
    logit(psi[i]) <- bstar0 + b[2] * con[i] + b[3] * cat[i]
    y[i] ~ dbern(mu[i])  # observation model, eqn. 6
    mu[i] <- (n1 / (n1 + pi * n0)) * z[i]
  }
}
