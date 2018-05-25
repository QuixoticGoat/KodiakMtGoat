
################################################################################
# Creates "available" fixes based on step length and turning angles of "used"  #
# GPS collar fixes. Extracts covariates from a rasterstack to the              #
# used/available dataframe.Performs conditional RSF model selection via AICc   #
# to quantify resource selection.                                              #
#                                                                              #
# Author: McCrea Cobb <mccrea.cobb@fws.gov>                                    #
# Last modified: 5/24/2018                                                     #
################################################################################

## ---- 4_ConditionalRSFmodel.R


# Load the GPS collar data
load("./Data/Rdata/df.Rdata")



#-------------------------------------------------------------------------------
## STEP 1: Create a dataset of "available" points. 10 points associated with
## each GPS fix.


# A. Create a traj object for each CollarID. Summarize the movements of 
#    each animal.

get.traj <- function(df) {
  ## Converts a df of GPS fixes to a trajectory object
  ## Step 1: Adds UTM coordinates and ordered by time
  ## Step 2: Creates and outputs a list of trajectory objects, by CollarID

  # required packages:
  require(adehabitatLT)

  # Order the df by CollarID and Date:

  xyConv <- function(df, xy = c('Long', 'Lat'), CRSin = '+proj=longlat',
                     CRSout = '+proj=utm +zone=5') {
    # Function to add UTM coordinates to the dataframe:

    df <- df[complete.cases(df[, xy]), ]
    conv <- SpatialPoints(coordinates(cbind('x' = df[, xy[1]],
                                            'y' = df[, xy[2]])),
                          proj4string = CRS(CRSin))
    conv <- spTransform(conv, CRS(CRSout))
    conv <- data.frame(conv)
    colnames(conv) <- c('x', 'y')
    df <- cbind(df, conv)

    df <- df[with(df, order(CollarID, Date)), ]

    return(df)
  }
  dat <- xyConv(df)


  to_ltraj <- function(dat) {
    ## Function to create a trajectory object

    dat <- dat[complete.cases(dat[, c("x", "y", "Date")]), ]
    traj <- adehabitatLT::as.ltraj(dat[, c("x", "y")], date = dat$Date, id = dat$CollarID,
                                   proj4string = CRS("+proj=utm +zone=5 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
    return(traj)
  }
  df.traj <- to_ltraj(dat)
  return(df.traj)
}


# Create trajectory objects seperately for ATS collars and Telonics collar data:

df.traj.ats <- get.traj(droplevels(subset(df, Collar == "ATS")))
df.traj.tel <- get.traj(droplevels(subset(df, Collar == "Telonics")))


# Load study area mask:
load("./Data/GIS/KodiakBound/StudyArea.RData")


## B. Draw 10 random "available" steps for each "used" step (takes a little time..)
library(hab); library(beepr)
df.ua.ats <- hab::rdSteps(df.traj.ats, nrs = 10, only.others = TRUE, mask = StudyArea) ; beep()
df.ua.tel <- hab::rdSteps(df.traj.tel, nrs = 10, only.others = TRUE, mask = StudyArea) ; beep()

# Merge the ATS and Telonics dfs back together:
df <- rbind(df.ua.ats, df.ua.tel)
rm(df.ua.ats, df.ua.tel)

# Compute "available" UTM coordinates based on the random step lengths and turning angles:
df$utmE <- df$x + df$dx
df$utmN <- df$y + df$dy

# Save the results:
save(data = df, file = "./Data/Rdata/df.ua.Rdata")



#-------------------------------------------------------------------------------
## STEP 2: Extract spatial covariate raster values to used and available points

###### Training data

## Load the data:
load("./Data/RData/df.ua.RData")

## Import the raster stack:
require(rgeos)
library(raster)
Covar.all <- stack("./Data/GIS/RasterStack/Covar.all.grd")

# Create a spatial df of the df.ua dataframe:
df <- SpatialPointsDataFrame(df[ ,c("utmE", "utmN")], df,
                             proj4string = CRS("+proj=utm +zone=5 +datum=NAD83"))

# Extract the covariate data to the used/available spatial dataframe (df.ua):
dfTemp <- extract(Covar.all, df[, c("utmE", "utmN")],
                  method = "simple",
                  df = T) ; beepr::beep(4)
df <- cbind(df, dfTemp) ; rm(dfTemp)

# Save it:
save(df, file="./Data/RData/df.cond.RData", compress="gzip")



#------------------------------------------------------------------------------
## STEP 3. Split data into testing/training datasets


test.train.fn <- function(df, p) {
  ## Function to split df into a list of testing and training dfs
  ## df = a dataframe of GPS collar fixes
  ## p = the propotion of the sample to include in the testing df.

  # 1. Create a sample size and randomly select rows:
  SmpSize <- floor(p * nrow(df))
  TrainInd <- sample(seq_len(nrow(df)), size = SmpSize)
  # Subset the data:
  dfTrain <- df[TrainInd, ]
  dfTest <- df[-TrainInd, ]
  dfTrain <- dfTrain[ order(dfTrain$Date), ]
  dfTest <- dfTest[ order(dfTest$Date), ]
  # Create a list of the dfs to return:
  df <- list("Test" = dfTest, "Train" = dfTrain)
  return(df)
}
# Run it and save:
dfTestTrain <- test.train.fn(df, 0.80)
save(data = dfTestTrain, file = "./Data/Rdata/df.testtrain.Rdata",
     compress = "gzip")



#------------------------------------------------------------------------------
## STEP 4: Run the conditional mixed effects logistic regression 


#library(survival)  # For clogit function - Conditional logistic regression

# Load the data:
load("./Data/RData/df.testtrain.Rdata")

# Create a dataframe of just the training data:
df <- dfTestTrain$Train

# Required library:
library(mclogit)  # For mixed effects conditional logistic regression (mclogit)
source("./Analysis/AIC.mclogit.R")  # For the AIC.mclogit function


################
## A. Summer
################

df <- subset(df, season=="summer")

### STEP 1: GRAIN. Determine the best spatial scale (grain) for vrm and slope covariates
### using univariate models (Lowrey et al. 2017)

mod <- list()
mod[[1]] <- mclogit(cbind(case, strata)~slopeS.30, random=~1|year/id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[2]] <- mclogit(cbind(case, strata)~slopeS.100, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[3]] <- mclogit(cbind(case, strata)~slopeS.500, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[4]] <- mclogit(cbind(case, strata)~slopeS.1000, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[5]] <- mclogit(cbind(case, strata)~vrmS.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[6]] <- mclogit(cbind(case, strata)~vrmS.100, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[7]] <- mclogit(cbind(case, strata)~vrmS.500, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[8]] <- mclogit(cbind(case, strata)~vrmS.1000, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))

Covariate <-  c("slopeS.30", "slopeS.100", "slopeS.500", "slopeS.1000", "vrmS.30", "vrmS.100", "vrmS.500", "vrmS.1000")

require(plyr)
mod[["aic"]] <- ldply(mod, function(x) AIC.mclogit(x))
mod[["aic"]]["AIC"] <- mod[["aic"]]$V1
mod[["aic"]] <- mod[["aic"]]["AIC"]
mod[["aic"]] <- cbind(Covariate, mod[["aic"]])
mod[["aic"]] <- mod[["aic"]][order(mod[["aic"]]["AIC"]), ]
View(mod[["aic"]])

# Save memory by removing the data from the model objects:
for(i in 1:length(mod)) {
  mod[[i]]$data <- c()
}

save(mod, file = "./Data/Rdata/ModG.s.cond.Rdata", compress = "gzip")
# rm(mod, Covariate)




### STEP 2: TERRAIN. Create a candidate list of mixed effects
### logistic models of terrain covariates

mod <- list()
mod[[1]] <- mclogit(cbind(case, strata) ~ elevS.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[2]] <- mclogit(cbind(case, strata) ~ elevS.30 + elevS2.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[3]] <- mclogit(cbind(case, strata) ~ slopeS.500, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[4]] <- mclogit(cbind(case, strata) ~ slopeS.500 + slopeS2.500, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[5]] <- mclogit(cbind(case, strata) ~ aspectS.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[6]] <- mclogit(cbind(case, strata) ~ vrmS.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[7]] <- mclogit(cbind(case, strata) ~ elevS.30 + elevS2.30 + vrmS.30 + aspectS.30 + slopeS.30 + slopeS2.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[8]] <- mclogit(cbind(case, strata) ~ elevS.30 + elevS2.30 + vrmS.30 + aspectS.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[9]] <- mclogit(cbind(case, strata) ~ elevS.30 + elevS2.30 + vrmS.30 + slopeS.30 + slopeS2.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[10]] <- mclogit(cbind(case, strata) ~ elevS.30 + elevS2.30 + aspectS.30 + slopeS.30 + slopeS2.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[11]] <- mclogit(cbind(case, strata) ~ vrmS.30 + aspectS.30 + slopeS.30 + slopeS2.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))


Covariate <-  c("elevS.30", "elevS2.30", "slopeS.500", "slopeS2.500", "aspectS.30", "vrmS.30", "elevS2.30 + slopeS2.500 + aspectS.30 + vrmS.30",
                "elevS2.30 + aspectS.30 + vrmS.30", "elevS2.30 + slopeS2.500 + vrmS.30", "elevS2.30 + slopeS2.500 + aspectS.30",
                "slopeS2.500 + aspectS.30 + vrmS.30")

require(plyr)
mod[["aic"]] <- ldply(mod, function(x) AIC.mclogit(x))
mod[["aic"]]["AIC"] <- mod[["aic"]]$V1
mod[["aic"]] <- mod[["aic"]]["AIC"]
mod[["aic"]] <- cbind(Covariate, mod[["aic"]])
mod[["aic"]] <- mod[["aic"]][order(mod[["aic"]]["AIC"]), ]
View(mod[["aic"]])

# Save memory by removing the data from the model objects:
for(i in 1:length(mod)) {
  mod[[i]]$data <- c()
}

save(mod, file = "./Data/Rdata/ModT.s.cond.Rdata", compress = "gzip")
rm(mod, Covariate)



### STEP 3. HABITAT MODEL: Create a candidate list of mixed effects
###    logistic models of habitat covariates


## A) First, select the best scale (30, 100, 500, 100 m) for each habitat variable:

mod <- list()

mod[[1]] <- mclogit(cbind(case, strata) ~ forest.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[2]] <- mclogit(cbind(case, strata) ~ forest.100, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[3]] <- mclogit(cbind(case, strata) ~ forest.500, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[4]] <- mclogit(cbind(case, strata) ~ forest.1000, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))

mod[[5]] <- mclogit(cbind(case, strata) ~ shrub.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[6]] <- mclogit(cbind(case, strata) ~ shrub.100, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[7]] <- mclogit(cbind(case, strata) ~ shrub.500, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[8]] <- mclogit(cbind(case, strata) ~ shrub.1000, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))

mod[[9]] <- mclogit(cbind(case, strata) ~ tundraHeath.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[10]] <- mclogit(cbind(case, strata) ~ tundraHeath.100, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[11]] <- mclogit(cbind(case, strata) ~ tundraHeath.500, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[12]] <- mclogit(cbind(case, strata) ~ tundraHeath.1000, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))

mod[[13]] <- mclogit(cbind(case, strata) ~ meadow.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[14]] <- mclogit(cbind(case, strata) ~ meadow.100, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[15]] <- mclogit(cbind(case, strata) ~ meadow.500, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[16]] <- mclogit(cbind(case, strata) ~ meadow.1000, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))

mod[[17]] <- mclogit(cbind(case, strata) ~ water.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[18]] <- mclogit(cbind(case, strata) ~ water.100, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[19]] <- mclogit(cbind(case, strata) ~ water.500, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[20]] <- mclogit(cbind(case, strata) ~ water.1000, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))

Covariate <-  c("forest.30", "forest.100", "forest.500", "forest.1000", "shrub.30", "shrub.100", "shrub.500", "shrub.1000",
                "tundra.30", "tundra.100", "tundra.500", "tundra.1000", "meadow.30", "meadow.100", "meadow.500", "meadow.1000",
                "water.30", "water.100", "water.500", "water.1000")

require(plyr)
mod[["aic"]] <- ldply(mod, function(x) AIC.mclogit(x))
mod[["aic"]]["AIC"] <- mod[["aic"]]$V1
mod[["aic"]] <- mod[["aic"]]["AIC"]
mod[["aic"]] <- cbind(Covariate, mod[["aic"]])
mod[["aic"]] <- mod[["aic"]][order(mod[["aic"]]["AIC"]), ]
View(mod[["aic"]])

# Save memory by removing the data from the model objects:
for(i in 1:length(mod)) {
  mod[[i]]$data <- c()
}

save(mod, file = "./Data/Rdata/ModH.s.cond.Rdata", compress = "gzip")
rm(mod, Covariate)


## B) Compare univariate habitat models at their "best" scale:

mod <- list()

mod[[1]] <- mclogit(cbind(case, strata) ~ elevS.30 + elevS2.30 + vrmS.30 + aspectS.30 + slopeS.30 + slopeS2.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[2]] <- mclogit(cbind(case, strata) ~ tundraHeath.100 + elevS.30 + elevS2.30 + vrmS.30 + aspectS.30 + slopeS.30 + slopeS2.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[3]] <- mclogit(cbind(case, strata) ~ shrub.100 + elevS.30 + elevS2.30 + vrmS.30 + aspectS.30 + slopeS.30 + slopeS2.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[4]] <- mclogit(cbind(case, strata) ~ forest.1000 + elevS2.30 + vrmS.30 + aspectS.30 + slopeS.30 + slopeS2.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[5]] <- mclogit(cbind(case, strata) ~ meadow.1000 + elevS.30 + elevS2.30 + vrmS.30 + aspectS.30 + slopeS.30 + slopeS2.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[6]] <- mclogit(cbind(case, strata) ~ water.1000 + elevS.30 + elevS2.30 + vrmS.30 + aspectS.30 + slopeS.30 + slopeS2.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[7]] <- mclogit(cbind(case, strata) ~ tundraHeath.100 + shrub.100 + elevS.30 + elevS2.30 + vrmS.30 + aspectS.30 + slopeS.30 + slopeS2.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[8]] <- mclogit(cbind(case, strata) ~ tundraHeath.100 + shrub.100 + water.1000 + elevS.30 + elevS2.30 + vrmS.30 + aspectS.30 + slopeS.30 + slopeS2.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[9]] <- mclogit(cbind(case, strata) ~ tundraHeath.100 + shrub.100 + water.1000 + elevS.30 + elevS2.30 + vrmS.30 + aspectS.30 + slopeS.30 + slopeS2.30 + sex, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))

Covariate <-  c("null", "tundra.100", "shrub.100", "forest.1000", "meadow.1000", "water.1000", "tundra.100 + shrub.100",
                "tundra.100 + shrub.100 + water.1000", "tundra.100 + shrub.100 + water.1000 + sex")

require(plyr)
mod[["aic"]] <- ldply(mod, function(x) AIC.mclogit(x))
mod[["aic"]]["AIC"] <- mod[["aic"]]$V1
mod[["aic"]] <- mod[["aic"]]["AIC"]
mod[["aic"]] <- cbind(Covariate, mod[["aic"]])
mod[["aic"]] <- mod[["aic"]][order(mod[["aic"]]["AIC"]), ]
View(mod[["aic"]])

# Save memory by removing the data from the model objects:
for(i in 1:length(mod)) {
  mod[[i]]$data <- c()
}

save(mod, file = "./Data/Rdata/ModH.s.cond.Rdata", compress = "gzip")
rm(mod, Covariate)

## Save top model
best <- mclogit(cbind(case, strata) ~ tundraHeath.100 + shrub.100 + water.1000 + elevS.30 + elevS2.30 + vrmS.30 + aspectS.30 + slopeS.30 + slopeS2.30 + sex, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
save(best, file = "./Data/Rdata/mod_summer_final.Rdata", compress = "gzip")

load(file = "./Data/Rdata/mod_summer_final.Rdata")



################
## B. Winter
################

# Load the data:
load("./Data/RData/df.cond.Rdata")
df <- subset(df, season=="winter")

### STEP 1: GRAIN. Determine the best spatial scale (grain) for each covariate
### using univariate models (Lowrey et al. 2017)

mod <- list()
mod[[1]] <- mclogit(cbind(case, strata)~slopeS.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[2]] <- mclogit(cbind(case, strata)~slopeS.100, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[3]] <- mclogit(cbind(case, strata)~slopeS.500, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[4]] <- mclogit(cbind(case, strata)~slopeS.1000, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[5]] <- mclogit(cbind(case, strata)~vrmS.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[6]] <- mclogit(cbind(case, strata)~vrmS.100, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[7]] <- mclogit(cbind(case, strata)~vrmS.500, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[8]] <- mclogit(cbind(case, strata)~vrmS.1000, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))

Covariate <-  c("slopeS.30", "slopeS.100", "slopeS.500", "slopeS.1000", "vrmS.30", "vrmS.100", "vrmS.500", "vrmS.1000")

require(plyr)
mod[["aic"]] <- ldply(mod, function(x) AIC.mclogit(x))
mod[["aic"]]["AIC"] <- mod[["aic"]]$V1
mod[["aic"]] <- mod[["aic"]]["AIC"]
mod[["aic"]] <- cbind(Covariate, mod[["aic"]])
mod[["aic"]] <- mod[["aic"]][order(mod[["aic"]]["AIC"]), ]
View(mod[["aic"]])

# Save memory by removing the data from the model objects:
for(i in 1:length(mod)) {
  mod[[i]]$data <- c()
}

save(mod, file = "./Data/Rdata/ModG.w.cond.Rdata", compress = "gzip")
rm(mod)



### STEP 2: TERRAIN. Create a candidate list of mixed effects
### logistic models of terrain covariates

mod <- list()
mod[[1]] <- mclogit(cbind(case, strata) ~ elevS.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[2]] <- mclogit(cbind(case, strata) ~ elevS.30 + elevS2.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[3]] <- mclogit(cbind(case, strata) ~ slopeS.100, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[4]] <- mclogit(cbind(case, strata) ~ slopeS.100 + slopeS2.100, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[5]] <- mclogit(cbind(case, strata) ~ aspectS.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[6]] <- mclogit(cbind(case, strata) ~ vrmS.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[7]] <- mclogit(cbind(case, strata) ~ elevS.30 + elevS2.30 + vrmS.30 + aspectS.30 + slopeS.100 + slopeS2.100, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[8]] <- mclogit(cbind(case, strata) ~ elevS.30 + elevS2.30 + vrmS.30 + aspectS.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[9]] <- mclogit(cbind(case, strata) ~ elevS.30 + elevS2.30 + vrmS.30 + slopeS.100 + slopeS2.100, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[10]] <- mclogit(cbind(case, strata) ~ elevS.30 + elevS2.30 + aspectS.30 + slopeS.100 + slopeS2.100, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[11]] <- mclogit(cbind(case, strata) ~ vrmS.30 + aspectS.30 + slopeS.100 + slopeS2.100, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))

Covariate <-  c("elevS.30", "elevS2.30", "slopeS.100", "slopeS2.100", "aspectS.30", "vrmS.30", "elevS2.30 + slopeS2.100 + aspectS.30 + vrmS.30",
                "elevS2.30 + aspectS.30 + vrmS.30", "elevS2.30 + slopeS2.100 + vrmS.30", "elevS2.30 + slopeS2.100 + aspectS.30",
                "slopeS2.100 + aspectS.30 + vrmS.30")

require(plyr)
mod[["aic"]] <- ldply(mod, function(x) AIC.mclogit(x))
mod[["aic"]]["AIC"] <- mod[["aic"]]$V1
mod[["aic"]] <- mod[["aic"]]["AIC"]
mod[["aic"]] <- cbind(Covariate, mod[["aic"]])
mod[["aic"]] <- mod[["aic"]][order(mod[["aic"]]["AIC"]), ]
View(mod[["aic"]])

# Save memory by removing the data from the model objects:
for(i in 1:length(mod)) {
  mod[[i]]$data <- c()
}

save(mod, file = "./Data/Rdata/ModT.w.cond.Rdata", compress = "gzip")
rm(mod, Covariate)



## STEP 3. HABITAT MODEL: Create a candidate list of mixed effects
## logistic models of habitat covariates

## A) First, select the best scale (30, 100, 500, 100 m) for each habitat variable:

mod <- list()

mod[[1]] <- mclogit(cbind(case, strata) ~ forest.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[2]] <- mclogit(cbind(case, strata) ~ forest.100, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[3]] <- mclogit(cbind(case, strata) ~ forest.500, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[4]] <- mclogit(cbind(case, strata) ~ forest.1000, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))

mod[[5]] <- mclogit(cbind(case, strata) ~ shrub.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[6]] <- mclogit(cbind(case, strata) ~ shrub.100, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[7]] <- mclogit(cbind(case, strata) ~ shrub.500, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[8]] <- mclogit(cbind(case, strata) ~ shrub.1000, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))

mod[[9]] <- mclogit(cbind(case, strata) ~ tundraHeath.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[10]] <- mclogit(cbind(case, strata) ~ tundraHeath.100, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[11]] <- mclogit(cbind(case, strata) ~ tundraHeath.500, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[12]] <- mclogit(cbind(case, strata) ~ tundraHeath.1000, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))

mod[[13]] <- mclogit(cbind(case, strata) ~ meadow.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[14]] <- mclogit(cbind(case, strata) ~ meadow.100, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[15]] <- mclogit(cbind(case, strata) ~ meadow.500, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[16]] <- mclogit(cbind(case, strata) ~ meadow.1000, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))

mod[[17]] <- mclogit(cbind(case, strata) ~ water.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[18]] <- mclogit(cbind(case, strata) ~ water.100, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[19]] <- mclogit(cbind(case, strata) ~ water.500, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[20]] <- mclogit(cbind(case, strata) ~ water.1000, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))

Covariate <-  c("forest.30", "forest.100", "forest.500", "forest.1000", "shrub.30", "shrub.100", "shrub.500", "shrub.1000",
                "tundra.30", "tundra.100", "tundra.500", "tundra.1000", "meadow.30", "meadow.100", "meadow.500", "meadow.1000",
                "water.30", "water.100", "water.500", "water.1000")

require(plyr)
mod[["aic"]] <- ldply(mod, function(x) AIC.mclogit(x))
mod[["aic"]]["AIC"] <- mod[["aic"]]$V1
mod[["aic"]] <- mod[["aic"]]["AIC"]
mod[["aic"]] <- cbind(Covariate, mod[["aic"]])
mod[["aic"]] <- mod[["aic"]][order(mod[["aic"]]["AIC"]), ]
View(mod[["aic"]])

# Save memory by removing the data from the model objects:
for(i in 1:length(mod)) {
  mod[[i]]$data <- c()
}

save(mod, file = "./Data/Rdata/ModH.s.cond.Rdata", compress = "gzip")
rm(mod, Covariate)


## B) Compare univariate habitat models at their "best" scale:

mod <- list()

mod[[1]] <- mclogit(cbind(case, strata) ~ vrmS.30 + slopeS.100 + slopeS2.100 + aspectS.30 + elevS.30 + elevS2.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[2]] <- mclogit(cbind(case, strata) ~ forest.500 + vrmS.30 + slopeS.100 + slopeS2.100 + aspectS.30 + elevS.30 + elevS2.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[3]] <- mclogit(cbind(case, strata) ~ meadow.100 + vrmS.30 + slopeS.100 + slopeS2.100 + aspectS.30 + elevS.30 + elevS2.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[4]] <- mclogit(cbind(case, strata) ~ tundraHeath.30 + vrmS.30 + slopeS.100 + slopeS2.100 + aspectS.30 + elevS.30 + elevS2.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[5]] <- mclogit(cbind(case, strata) ~ shrub.100 + vrmS.30 + slopeS.100 + slopeS2.100 + aspectS.30 + elevS.30 + elevS2.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[6]] <- mclogit(cbind(case, strata) ~ water.30 + vrmS.30 + slopeS.100 + slopeS2.100 + aspectS.30 + elevS.30 + elevS2.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[7]] <- mclogit(cbind(case, strata) ~ tundraHeath.30 + forest.500 + vrmS.30 + slopeS.100 + slopeS2.100 + aspectS.30 + elevS.30 + elevS2.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[8]] <- mclogit(cbind(case, strata) ~ tundraHeath.30 + forest.500 + shrub.100 + vrmS.30 + slopeS.100 + slopeS2.100 + aspectS.30 + elevS.30 + elevS2.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[9]] <- mclogit(cbind(case, strata) ~ tundraHeath.30 + forest.500 + shrub.100 + water.30 + vrmS.30 + slopeS.100 + slopeS2.100 + aspectS.30 + elevS.30 + elevS2.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
mod[[10]] <- mclogit(cbind(case, strata) ~ tundraHeath.30 + forest.500 + shrub.100 + water.30 + vrmS.30 + slopeS.100 + slopeS2.100 + aspectS.30 + elevS.30 + elevS2.30 + sex, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))

Covariate <-  c("null", "forest.500", "meadow.100", "tundra.30", "shrub.100", "water.30", "tundra.30 + forest.500",
                "tundra.30 + forest.500 + shrub.100", "tundra.30 + forest.500 + shrub.100 + water.30",
                "tundra.30 + forest.500 + shrub.100 + water.30 + sex")

require(plyr)
mod[["aic"]] <- ldply(mod, function(x) AIC.mclogit(x))
mod[["aic"]]["AIC"] <- mod[["aic"]]$V1
mod[["aic"]] <- mod[["aic"]]["AIC"]
mod[["aic"]] <- cbind(Covariate, mod[["aic"]])
mod[["aic"]] <- mod[["aic"]][order(mod[["aic"]]["AIC"]), ]
View(mod[["aic"]])

# Save memory by removing the data from the model objects:
for(i in 1:length(mod)) {
  mod[[i]]$data <- c()
}

save(mod, file = "./Data/Rdata/ModH.w.cond.Rdata", compress = "gzip")
rm(mod, Covariate)

## Save top model
best <- mclogit(cbind(case, strata) ~ tundraHeath.30 + forest.500 + shrub.100 + water.30 + vrmS.30 + slopeS.100 + slopeS2.100 + aspectS.30 + elevS.30 + elevS2.30, random=~1|id, data=df, control = mclogit.control(epsilon = 1e-03, maxit = 100))
save(best, file = "./Data/Rdata/mod_winter_final.Rdata", compress = "gzip")




#------------------------------------------------------------------------------
## STEP 5 Evaluate fit of the best model with the testing data


# Load the testing data:
load(file="./Data/RData/df.ua.RData")

df <- df.ua.RData$Test

# Filter out rows with NAs:
df <- df[complete.cases(df),]

# Predict values:
df$predicted <- predict(best, newdata = df, type = "response")


plot(df$case ~ df$predicted)

