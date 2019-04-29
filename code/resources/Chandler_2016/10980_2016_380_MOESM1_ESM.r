
## ----------------------------- Utility functions --------------------------------------

## Function to generate spatial covariate as multivariate Gaussian process with exponential correlation.
## Modified from Andy Royle's code
spcov <- function(R, alpha=2, standardize=TRUE) {
    v <- sqrt(nrow(R))
    D <- as.matrix(dist(R)) #e2dist1(R,R)
    V <- exp(-D/alpha)
    cov1 <- t(chol(V)) %*% rnorm(nrow(R))
    Rd <- as.data.frame(R)
    colnames(Rd) <- c("x", "y")
    if(standardize)
        Rd$z <- as.numeric((cov1 - mean(cov1)) / sd(cov1))
    else
        Rd$z <- as.numeric(cov1)
    return(Rd)
}




## Compute Euclidean distance between two sets of points
e2dist <- function (x, y) {
    i <- sort(rep(1:nrow(y), nrow(x)))
    dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
    matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}






## ----------------------------- Likelihood functions -----------------------------------

## Negative log-likelihood
nll <- function(p, N) {
    sig <- exp(p[1])
    beta0 <- p[2]
    beta1 <- p[3]
    beta2 <- p[4]
    cov.w <- apply(D, 1, function(x) {
        w0 <- exp(-x^2 / (2*sig^2))
        w0[which.min(x)] <- 0
        w <- w0/sum(w0)
        sum(cov$z * w)
    })
    lambda <- exp(beta0 + beta1*cov0 + beta2*cov.w)
    -sum(dpois(N, lambda, log=TRUE))
}


nll2 <- function(p, N) {
    beta0 <- p[1]
    beta1 <- p[2]
    lambda <- exp(beta0 + beta1*cov0)
    -sum(dpois(N, lambda, log=TRUE))
}







## --------------------- Simulate data and fit model under one of the 40 cases -------------------------



nSim <- 1000
simout <- matrix(NA, nSim, 4)
colnames(simout) <- c("sigma", "beta0", "beta1", "beta2")#, "conv", "nll")
##                       "beta0p", "beta1p", "AIC1", "AIC2")
sigma <- 0.03
beta0 <- 2
beta1 <- 1 ## Patch level effect
beta2 <- 0.5 ## Landscape level effect
B <- 50
co <- seq(0, 1, length=B)
Z <- cbind(rep(co, each=B), rep(co, times=B))


## Matrix to hold simulated spatial covariates
covout <- matrix(NA, nrow(Z), nSim)


## Survey locations
l <- 10
cop <- seq(0.2, 0.8, length=l)
X <- cbind(rep(cop, each=l), rep(cop, times=l))


## Distance between survey locations and pixels
D <- e2dist(X, Z)


## Matrix to hold simulated N
Nout <- matrix(NA, nrow(X), nSim)

## List to hold full optim output
olist <- list()

set.seed(90340)
for(i in 1:nSim) {
    cat("doing", i, "\n")
    cat("   ", format(Sys.time()), "\n")
    cov <- spcov(Z, alpha=0.01)
    cov.w <- apply(D, 1, function(x) {
        w0 <- exp(-x^2 / (2*sigma^2))
        ## Ignore habitat in focal site.
        w0[which.min(x)] <- 0
        w <- w0/sum(w0)
        sum(cov$z * w)
    })
    covout[,i] <- cov$z
    cov0 <- rnorm(nrow(D))  ##apply(D, 1, which.min)]
    lambda <- exp(beta0 + beta1*cov0 + beta2*cov.w)
    N <- rpois(length(lambda), lambda)
    Nout[,i] <- N
    o1 <- try(optim(c(log(0.05), 0, 0, 0), nll, N=N, hessian=TRUE,
                    lower=c(-6, -Inf, -Inf, -Inf),
                    upper=c(1, Inf, Inf, Inf)))
    if(identical(class(o1), "try-error")) {
        o1 <- list()
        o1$par <- rep(NA, 4)
        o1$hessian <- matrix(NA, 4, 4)
        o1$convergence <- -99
    }
    simout[i,1:4] <- o1$par
    olist[[i]] <- o1
    cat("     ", round(simout[i,], 2), "\n\n")
}






