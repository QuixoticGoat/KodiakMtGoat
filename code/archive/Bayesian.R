#### Bayesian RSF (logistic regression)
#### Author: McCrea Cobb
#### Last edited: 9/19/2017


library(R2WinBUGS)
bugs.dir <- "C:/Program Files/WinBUGS14/"          # Place where your WinBUGS installed


#Load the libraries for parallel processing
library(snow)
library(snowfall)
# Set the number of CPUs to be 2
sfInit(parallel=TRUE, cpus=2)

# Assign the R2WinBUGS library to each CPU
sfLibrary(R2WinBUGS)


####################################################
### Need to deal with covariates containing NAs ####
####################################################



# Bundle data
win.data <- list(Response = d.bug$Response,
                 n = length(d.bug$Response),
                 slope = d.bug$slopeS)



# 5. creating separate directory for each CPU process
folder1 <- paste(getwd(), "/chain1", sep="")
folder2 <- paste(getwd(), "/chain2", sep="")
dir.create(folder1); dir.create(folder2)

# 6. sinking the model into a file in each directory
for (folder in c(folder1, folder2))
{
  sink(paste(folder, "/RSF.txt", sep=""))

  # Specify model in BUGS language
  cat(file = "RSF.txt","
    model {

    # Priors
    alpha ~ dnorm(0, 0.0001)     # intercept
    b.slope ~ dnorm(0, 0.0001)   # betas for slope, normally distributed

    # Likelihood
    for (i in 1:n){
    Response[i] ~ dbern(p[i])
    logit(p[i]) <- alpha + b.slope*slope[i]
    }
    }
    ")
  sink()
}

# 7. Define the function that will run MCMC on each CPU
# Arguments:
# chain - will be 1, 2 or 3
# win.data - the data list
# params - parameters to be monitored
parallel.bugs <- function(chain, win.data, params, seeds)
{
  # 7a. defining directory for each CPU
  sub.folder <- paste(getwd(),"/chain", chain, sep="")
  # Initial values
  inits <- function() list(alpha = 0, b.slope = 1)

  # Parameters monitored
  params <- c("alpha", "b.slope", "p")

  # MCMC settings
  ni <- 1000   ;   nt <- 1   ;   nb <- 200   ;  nc <- 1

  bugs.dir <- "C:/Program Files/WinBUGS14/"          # Place where your WinBUGS installed

  # Call WinBUGS
  bugs(data = win.data, inits = inits, params, model.file = "RSF.txt",
       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
       debug = TRUE, bugs.directory = bugs.dir, working.directory = sub.folder,
       bugs.seed = seeds[chain], codaPkg = T)
}




# 9. Call the sfLapply function that will run
# parallel.bugs on each of the 3 CPUs:
sfLapply(1:2, fun=parallel.bugs, win.data=win.data, params=params, seeds=c(1,2))




# Locating position of each CODA chain:
chain1 <- paste(folder1, "/CODAchain1.txt", sep="")
chain2 <- paste(folder2, "/CODAchain1.txt", sep="")

# Get the results:
a <- read.bugs(c(chain1, chain2))
plot(a)
