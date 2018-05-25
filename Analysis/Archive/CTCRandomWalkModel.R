## Continuous-time correlated random walk model
## McCrea Cobb
## 5/11/2017


## Required packages
library(crawl)


crwMod <- crwMLE(mov.model = ~1, 
                 err.model = NULL, 
                 activity = NULL, 
                 drift = FALSE,
                 data = df 
                 coord = c("Long", "Lat"), 
                 Time.name = "Date"
                 initial.state, 
                 theta, 
                 fixPar,
                 method = "Nelder-Mead", 
                 control = NULL, 
                 constr = list(lower = -Inf, upper = Inf), 
                 prior = NULL, 
                 need.hess = TRUE, 
                 initialSANN = list(maxit = 200), 
                 attempts = 1)



## Example:
