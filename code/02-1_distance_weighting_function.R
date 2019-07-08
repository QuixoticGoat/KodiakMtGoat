################################################################################
# Creates a distance weighting layer as a predictor in an RSF model            #
#                                                                              #
# Author: McCrea Cobb <mccrea_cobb@fws.gov                                     #
# Date created: 05-13-2019                                                     #
#                                                                              #
# Lasted edited by: McCrea Cobb <mccrea_cobb@fws.gov                           #
# Date last edited: 05-13-2019                                                 #
################################################################################



#-------------------------------------------------------------------------------
##----weighting_functions

## From Miguet et al. 2017. How to quantify a distance-dependent landscape
## effect on a biological response. Methods in Ecology and Evolution 8:1717-1724.



# 1] Four families of weighting functions: threshold, negative exponential, 
#    Gaussian, exponential power

# Threshold functions:

# A: one parameter, X: distance
fThreshold_0 = function(A,X) { # Non-standardized version
  Y=X<A
  return(Y)
}
fThreshold <- function(A,X) { # standardized version
  Y=X<A
  Y=Y/sum(Y)
  return(Y)
}
# Negative exponential functions
# A: one parameter, X: distance
fNExp_0 <- function(A,x) { # Non-standardized version
  res=exp(-A*x)
  return(res)
}
fNExp <- function(A,X) { # standardized version
  Y=exp(-A*X)
  Y=Y/sum(Y)
  return(Y)
}
# Gaussian functions
# A: one parameter, X: distance
fGaus_0 = function (A,X) { # Non-standardized version
  Y=exp(-A[1]*(X^2))
  return(Y)
}
fGaus = function (A,X) { # standardized version
  Y=exp(-A[1]*(X^2))
  Y=Y/sum(Y)
  return(Y)
}
# Exponential power functions
# A: two parameters, X: distance
fNExpb_0 = function (A,X) { # Non-standardized version
  Y=exp(-A[1]*(X^A[2])) 
  return(Y)
}
fNExpb = function (A,X) { # standardized version
  Y=exp(-A[1]*(X^A[2]))
  Y=Y/sum(Y)
  return(Y)
}


#-------------------------------------------------------------------------------
## ---- calculate_weighted_variable

## From Miguet et al. 2017. How to quantify a distance-dependent landscape
## effect on a biological response. Methods in Ecology and Evolution 8:1717-1724.



# 3] Calculation of the weighted landscape variable
# Input
# 	dataLandscape: landscape variable value for each ring and each focal sampling unit
# 	dataArea: area of the rings within which the landscape variable is calculated for each sampling unit
# 	sList: list of distances used in the analysis (one distance per ring, the radius of the external radius of each ring)
# 	fW: weighting function
# 	A: parameter(s) for the weighting function
# 	distMax: value above which the landscape weight should be zero
# Output
# 	the weighted landscape variable, one value per focal sampling unit

MakeVarW = function(dataLandscape, dataArea, sList, fW, A, distMax) {
  Ns  = length(sList)
  N = dim(dataLandscape)[1]
  sm1List = c(0, sList[1:(Ns-1)])
  s2List = 2/3*(sList^3 - sm1List^3)/(sList^2 - sm1List^2)
  W = fW(A, s2List)
  W[s2List > distMax] = 0 
  WS = dataArea * t(matrix(W,Ns,N))
  WSSum = rowSums(WS)
  VarW = rowSums(dataLandscape * WS)/WSSum
  return(VarW)
}


#-------------------------------------------------------------------------------
## ---- run_calculate_weighted_variable

## Runs the MakeVarW() function on the Kodiak goat data.






#-------------------------------------------------------------------------------
## Run the buff.extract function

library(raster)
library(velox)  # A faster spatial extract function based on C++

# Load the spatial dataframe of GPS fixes
load("./data/derived_data/geodata/dfSp.ua.Rdata")  

# Source buff.extract()
source("./code/functions/buff.extract.R")

# Read in the habitat covariate rasterstack and convert it to a vx object for speedier processing:
covar <- stack("./data/derived_data/geodata/covar.grd")
covar.vx <- velox::velox(covar)
rm(covar)

# Run buff.extract(), 50 m radius rings out to 1 km:
df1 <- buff.extract(dfSp, covar.vx, 50, 20, n)  
save(df1, file="./data/derived_data/temp/df1.Rdata")
rm(df1)

source("./code/functions/buff.extract.R")
df2 <- buff.extract(dfSp, covar.vx, 50, 20, n)
save(df2, file="./data/derived_data/temp/df2.Rdata")
rm(df2)

source("./code/functions/buff.extract.R")
df3 <- buff.extract(dfSp, covar.vx, 50, 20, n)
save(df3, file="./data/derived_data/temp/df3.Rdata")
rm(df3)


# subset the data
df <- as.data.frame(dfSp)
df <- df[df$id %in% c("IG18", "IG20"), ]

foo2 <- cbind(df, foo1)

#-------------------------------------------------------------------------------

## 37 sec for 1000 fixes buffered to 1000 m (file = 315 Mb)...
## Need to split the dfSp into sections for analysis due to file size and memory limitations:

load("./data/derived_data/geodata/dfSp.ua.Rdata")
dfSp <- dfSp[1:100000,]
system.time(dfSp.buf50.1 <- multi.buffer(dfSp, 50, 20))
save(dfSp.buf50.1, file="./data/derived/geodata/dfSp.buf50.1.Rdata")
rm(dfSp.buf50.1, dfSp)
gc()

load("./data/derived_data/geodata/dfSp.ua.Rdata")
dfSp <- dfSp[100001:200000,]
system.time(dfSp.buf50.2 <- multi.buffer(dfSp, 50, 20))
save(dfSp.buf50.2, file="./data/derived_data/geodata/dfSp.buf50.2.Rdata")
rm(dfSp.buf50.2, dfSp)
gc()

load("./data/derived_data/geodata/dfSp.ua.Rdata")
dfSp <- dfSp[200001:300000,]
system.time(dfSp.buf50.3 <- multi.buffer(dfSp, 50, 20))
save(dfSp.buf50.3, file="./data/derived/geodata/dfSp.buf50.3.Rdata")
rm(dfSp.buf50.3, dfSp)
gc()

load("./data/derived_data/geodata/dfSp.ua.Rdata")
dfSp <- dfSp[300001:400000,]
system.time(dfSp.buf50.4 <- multi.buffer(dfSp, 50, 20))
save(dfSp.buf50.4, file="./data/derived/geodata/dfSp.buf50.4.Rdata")
rm(dfSp.buf50.4, dfSp)
gc()

load("./data/derived_data/geodata/dfSp.ua.Rdata")
dfSp <- dfSp[400001:500000,]
system.time(dfSp.buf50.5 <- multi.buffer(dfSp, 50, 20))
save(dfSp.buf50.5, file="./data/derived/geodata/dfSp.buf50.5.Rdata")
rm(dfSp.buf50.5, dfSp)

load("./data/derived_data/geodata/dfSp.ua.Rdata")
dfSp <- dfSp[500001:600000,]
system.time(dfSp.buf50.6 <- multi.buffer(dfSp, 50, 20))
save(dfSp.buf50.6, file="./data/derived/geodata/dfSp.buf50.6.Rdata")
rm(dfSp.buf50.6, dfSp)
gc()

load("./data/derived_data/geodata/dfSp.ua.Rdata")
dfSp <- dfSp[600001:700000,]
system.time(dfSp.buf50.7 <- multi.buffer(dfSp, 50, 20))
save(dfSp.buf50.7, file="./data/derived/geodata/dfSp.buf50.7.Rdata")
rm(dfSp.buf50.7, dfSp)
gc()

load("./data/derived_data/geodata/dfSp.ua.Rdata")
dfSp <- dfSp[700001:length(dfSp),]
system.time(dfSp.buf50.8 <- multi.buffer(dfSp, 50, 20))
save(dfSp.buf50.8, file="./data/derived/geodata/dfSp.buf50.8.Rdata")
rm(dfSp.buf50.8, dfSp)
gc()




# look at it:
plot(foo[1:100])
plot(foo[1001:1101], add =T)

plot(foo[[1]][1:100])
plot(foo[[2]][1:100], add=T)





#-------------------------------------------------------------------------------
## ---- extract

## Extracts the proportions of habitats within each ring surrounding use/available fixes.

library(raster)
library(velox)  # A faster extract function
library(future.apply)  # For parallel processing apply function

# Read in the habitat covariate rasterstack:
covar <- stack("./data/derived_data/geodata/covar.grd")

# Convert the rasterstack to a vx object of the habitat covariates (for speedier processing):
covar.vx <- velox::velox(covar)



# Read in the Rdata file containing the buffered rings:
system.time(load("./data/derived_data/geodata/dfSp.buf50.1.Rdata"))

# Extract covariate data from the rasterstack to the spatial dataframe (polygons):
system.time(covar.ex <- lapply(1:length(dfSp.buf50.1), function(x) {
  covar.vx$extract(sp = dfSp.buf50.1[[x]], fun = mean)
}))

# Convert output to a dataframe and rename columns using the rasterstack names:
covar.ex <- as.data.frame(covar.ex)
dist <- sort(rep.int(seq(50, 1000, 50), 13), decreasing=F)
colnames(covar.ex) <- paste0(names(covar), "_", dist)

save(covar.ex, file="./data/derived_data/covar.ex.1.Rdata")


#############
function() {
  lapply(1:8, function(x) {
    load(paste0("./data/derived/geodata/dfSp.buf50.", x, ".Rdata"))
    covar.ex <- lapply(1:length(paste0("dfSp.buf50.", x)), function(x) {
      covar.vx$extract(sp = dfSp.buf50.1[[x]], fun = mean)
    })
  })
  
}

#############



# Combine output dataframe with the original dataframe (df.ua):
dfRSF.new <- cbind(df.ua, covar.ex)
