
#################################################################
### RSF Model for GPS collar data
### Author: McCrea Cobb
### Last modified: 9/15/2017
#################################################################

## ---- 3_RSFmodel.R

# Load the used/available data frame:
load("./Data/RData/dfRSF.Rdata")





### First, check covariance among covariates:

cov.fn <- function(dfRSF) {
  dfRSFCov <- data.frame(cor(dfRSF[-c(1:4)],
                             use = "complete.obs",
                             method = "pearson"))
 # dfRSFCov <- dfRSFCov[5:16, 5:16]

  return(dfRSFCov)
  #write.csv(dfRSFCov, "./Output/Tables/Covariance.csv")
}

# Run it for the list of 30, 100, 500, and 100 m extracted data frames:
dfRSF.covs <- cov.fn(dfRSF)







## Split data into testing/training datasets

# Determine the sample size of the training set:

test.train.fn <- function(df) {

  # Create a sample size and randomly select rows:
  SmpSize <- floor(0.80 * nrow(df))
  TrainInd <- sample(seq_len(nrow(df)), size = SmpSize)

  # Subset the data for training:
  dfTrain <- df[TrainInd, ]
  attr(dfTrain, "id") <- "dfTrain"

  # Subset the remaining data for testing:
  dfTest <- df[-TrainInd, ]
  attr(dfTest, "id") <- "dfTest"


  return(list("Train"=dfTrain, "Test"=dfTest))
}

# Run it:
dfTT <- test.train.fn(dfRSF)

dfRSF.t <- as.data.frame(dfTT[[1]])








### STEP 1: GRAIN. Determine the best spatial scale (grain) for each covariate
### using univariate models (Lowrey et al. 2017)

mod.scale.fn <- function(df = dfRSF) {
  require(lme4)
  require(AICcmodavg)

  # Create a df of variables:
  df.scale <- data.frame("CollarID" = df$CollarID,
                         "Response" = df$Response,
                         "slope.30" = df$slopeS.30,
                         "slope.100" = df$slopeS.100,
                         "slope.500" = df$slopeS.500,
                         "slope.1000" = df$slopeS.1000,
                         "slope2.30" = df$slopeS2.30,
                         "slope2.100" = df$slopeS2.100,
                         "slope2.500" = df$slopeS2.500,
                         "slope2.1000" = df$slopeS2.1000,
                         "vrm.30" = df$vrmS.30,
                         "vrm.100" = df$vrmS.100,
                         "vrm.500" = df$vrmS.500,
                         "vrm.1000" = df$vrmS.1000)

  # Create a list of the names of the variables:
  pred <- names(df.scale[, -c(1:2)])

  modG <- lapply(df.scale[pred], function(x) glmer(Response ~ x + (1|CollarID),
                                                             data = df.scale,
                                                             family = binomial(link = "logit")) )

  # Create an AICc table:
  modG[["aic"]] <- aictab(cand.set = modG, modnames = pred, sort = T)
  modG$aic <- as.data.frame(modG$aic)
  modG$aic <- modG$aic[order(modG$aic$Delta_AICc), ]

  # Alert when completed:
  beepr::beep()    # Beep :)

  modG          # Display the results
  return(modG)  # Save the results
}

# Run it:
modG.w <- mod.scale.fn(dfRSF.w)


save(modG, file = "./Data/Rdata/ModG.Rdata", compress = "gzip")
rm(modG)







### STEP 2: TERRAIN. Create a candidate list of mixed effects
### logistic models of terrian covariates

require(lme4)
library(afex)  # for all_fit()

# Create the candidate models:
modT <- list()

modT[[1]] <- glmer(Response ~ elevS.30 + (1|CollarID),   # Elevation, 30m
                   data = dfRSF,
                   family = binomial(link = "logit"))
modT[[2]] <- glmer(Response ~ elevS.30 + elevS2.30 + (1|CollarID),   # Elevation^2, 30m
                  data = dfRSF,
                  family = binomial(link = "logit"))
modT[[3]] <- glmer(Response ~ slopeS.500 + (1|CollarID),           # Slope, 500m
                  data = dfRSF,
                  family = binomial(link = "logit"))
modT[[4]] <- glmer(Response ~ aspectS.30 + (1|CollarID),          # Aspect, 30m
                   data = dfRSF,
                   family = binomial(link = "logit"))
modT[[5]] <- glmer(Response ~ vrmS.500 + (1|CollarID),             # VRM, 500m
                  data = dfRSF,
                  family = binomial(link = "logit"))
modT[[6]] <- glmer(Response ~ 1 + (1|CollarID),                # Null ### Not converging..
                    data = dfRSF,
                   family = binomial(link = "logit"))
modT[[6]] <- all_fit(modT[[6]])  # To refit the null model using a range of optimizers.
modT[[6]] <- modT[[6]]$bobyqa.
modT[[7]] <- glmer(Response ~ slopeS.500 + vrmS.500 + (1|CollarID),           # Slope, 500m and vrm, 1000m
                   data = dfRSF,
                   family = binomial(link = "logit"))
modT[[8]] <- glmer(Response ~ slopeS.500 + vrmS.500 + elevS.30 + elevS2.30 + (1|CollarID),           # Slope, 500m and vrm, 1000m
                   data = dfRSF,
                   family = binomial(link = "logit"))
modT[[9]] <- glmer(Response ~ slopeS.500 + vrmS.500 + elevS.30 + elevS2.30 + aspectS.30 + (1|CollarID),           # Slope, 500m and vrm, 1000m
                   data = dfRSF,
                   family = binomial(link = "logit"))

### Check models for overdispersion:

source("./Analysis/OverdispersionFunction.R")
ModTDisp <- lapply("modT", function(x)
               overdisp_fun(x))
# None are overdispersed.



### Model selection

require(AICcmodavg)

# Name the list of models:
pred <- c("Elev.30", "Elev^2.30", "Slope.500", "Aspect.30", "vrm.1000",
          "(Null)", "Slope.500+vrm.1000","Slope.500+vrm.1000+Elev^2.30",
          "Slope.500+vrm.1000+Elev^2.30+Aspect.30")

# Create a model selection table:
modT$aic <- aictab(cand.set = modT[1:9], modnames = pred, sort = T)
modT$aic <- as.data.frame(modT$aic)
modT$aic <- modT$aic[order(modT$aic$Delta_AICc), ]

  # Top model includes slopeS with 500m grain. AICcWeight = 1


save(modT, file = "./Data/Rdata/ModT.Rdata", compress = "gzip")








### STEP 3. HABITAT MODEL: Create a candidate list of mixed effects
###    logistic models of habitat covariates


    ## A) First, select the best scale (30, 100, 500, 100 m) for each habitat variable:

mod.scale.fn <- function(df) {
  require(lme4)
  require(AICcmodavg)

  df.scale <- data.frame("CollarID" = df$CollarID,
                         "Response" = df$Response,
                         # "forest.30" = df$forest.30,
                         # "forest.100" = df$forest.100,
                         # "forest.500" = df$forest.500,
                         # "forest.1000" = df$forest.1000
                         # "shrub.30" = df$shrub.30,
                         # "shrub.100" = df$shrub.100,
                         # "shrub.500" = df$shrub.500,
                         # "shrub.1000" = df$shrub.1000
                         # "tundraHeath30" = df$tundraHeath.30,
                         # "tundraHeath.100" = df$tundraHeath.100,
                         # "tundraHeath.500" = df$tundraHeath.500,
                         # "tundraHeath.1000" = df$tundraHeath.1000
                         # "meadow.30" = df$meadow.30,
                         # "meadow.100" = df$meadow.100,
                         # "meadow.500" = df$meadow.500,
                         # "meadow.1000" = df$meadow.1000,
                         # "water.30" = df$water.30,
                         # "water.100" = df$water.100,
                         # "water.500" = df$water.500,
                         # "water.1000" = df$water.1000
                         # "snowWater.30" = df$snowWater.30,
                         # "snowWater.100" = df$snowWater.100,
                         # "snowWater.500" = df$snowWater.500,
                         # "snowWater.1000" = df$snowWater.1000),
                         "rock.30" = df$rock.30,
                         "rock.100" = df$rock.100,
                         "rock.500" = df$rock.500,
                         "rock.1000" = df$rock.1000)


  # Create a list of the names of the variables:
  pred <- names(df.scale[, -c(1:2)])

  mod <- lapply(df.scale[pred], function(x) glmer(Response ~ x + (1|CollarID),
                                                  data = df.scale,
                                                  family = binomial(link = "logit")) )

  # Create an AICc table:
  mod[["aic"]] <- aictab(cand.set = mod, modnames = pred, sort = T)
  mod$aic <- as.data.frame(mod$aic)
  mod$aic <- mod$aic[order(mod$aic$Delta_AICc), ]

  # Alert when completed:
  beepr::beep()    # Beep :)

  mod         # Display the results
  return(mod)  # Save the results
}

# Run it:
modG.forest <- mod.scale.fn(dfRSF)
modG.shrub <- mod.scale.fn(dfRSF)
modG.tundra <- mod.scale.fn(dfRSF)
modG.meadow <- mod.scale.fn(dfRSF)
modG.water <- mod.scale.fn(dfRSF)
modG.snowWater <- mod.scale.fn(dfRSF)
modG.rock <- mod.scale.fn(dfRSF)

save(modG.rock,
     file = "./Data/Rdata/ModG.rock.Rdata",
     compress = "gzip")






    ## B) Compare univariate habitat models at their "best" scale:

# Make an empty list to hold the candidate models:
modH <- list()

require(lme4)

# Create the candidate models:

modH[[1]] <- glmer(Response ~ slopeS.500 + vrmS.500 + elevS.30 + elevS2.30 +   # forest, 500m
                     aspectS.30 + forest.500 + (1|CollarID),
                   data = dfRSF, family = binomial(link = "logit"))
modH[[2]] <- glmer(Response ~ slopeS.500 + vrmS.500 + elevS.30 + elevS2.30 +   # shrub, 100m
                     aspectS.30 + shrub.100 + (1|CollarID),
                   data = dfRSF, family = binomial(link = "logit"))
modH[[3]] <- glmer(Response ~ slopeS.500 + vrmS.500 + elevS.30 + elevS2.30 +   # tundra, 100m
                     aspectS.30 + tundraHeath.100 + (1|CollarID),
                   data = dfRSF, family = binomial(link = "logit"))
modH[[4]] <- glmer(Response ~ slopeS.500 + vrmS.500 + elevS.30 + elevS2.30 +   # meadow, 500m
                     aspectS.30 + meadow.500 + (1|CollarID),
                   data = dfRSF, family = binomial(link = "logit"))
modH[[5]] <- glmer(Response ~ slopeS.500 + vrmS.500 + elevS.30 + elevS2.30 +   # water, 1000m
                     aspectS.30 + water.1000 + (1|CollarID),
                   data = dfRSF, family = binomial(link = "logit"))
modH[[6]] <- glmer(Response ~ slopeS.500 + vrmS.500 + elevS.30 + elevS2.30 +   # snow, 1000m
                     aspectS.30 + snowWater.1000 + (1|CollarID),
                   data = dfRSF, family = binomial(link = "logit"))
modH[[7]] <- glmer(Response ~ slopeS.500 + vrmS.500 + elevS.30 + elevS2.30 +   # rock, 1000m
                     aspectS.30 + rock.1000 + (1|CollarID),
                   data = dfRSF, family = binomial(link = "logit"))
modH[[8]] <- glmer(Response ~ slopeS.500 + vrmS.500 + elevS.30 + elevS2.30 +   # null
                     aspectS.30 + (1|CollarID),
                   data = dfRSF, family = binomial(link = "logit"))
modH[[9]] <- glmer(Response ~ slopeS.500 + vrmS.500 + elevS.30 + elevS2.30 +   # tundra + snow
                     aspectS.30 + tundraHeath.100 + snowWater.1000 + (1|CollarID),
                   data = dfRSF, family = binomial(link = "logit"))
modH[[10]] <- glmer(Response ~ slopeS.500 + vrmS.500 + elevS.30 + elevS2.30 +   # tundra + water
                     aspectS.30 + tundraHeath.100 + water.1000 + (1|CollarID),
                   data = dfRSF, family = binomial(link = "logit"))
modH[[11]] <- glmer(Response ~ slopeS.500 + vrmS.500 + elevS.30 + elevS2.30 +   # tundra + snow + water
                      aspectS.30 + tundraHeath.100 + water.1000 + snowWater.1000 + (1|CollarID),
                    data = dfRSF, family = binomial(link = "logit"))
modH[[12]] <- glmer(Response ~ slopeS.500 + vrmS.500 + elevS.30 + elevS2.30 +   # tundra + snow + water
                      aspectS.30 + tundraHeath.100 + water.1000 + snowWater.1000 + shrub.100 + (1|CollarID),
                    data = dfRSF, family = binomial(link = "logit"))
modH[[13]] <- glmer(Response ~ slopeS.500 + vrmS.500 + elevS.30 + elevS2.30 +   # tundra + snow + water
                      aspectS.30 + tundraHeath.100 + water.1000 + snowWater.1000 + shrub.100 + meadow.500 + (1|CollarID),
                    data = dfRSF, family = binomial(link = "logit"))
modH[[14]] <- glmer(Response ~ slopeS.500 + vrmS.500 + elevS.30 + elevS2.30 +   # tundra + snow + water
                      aspectS.30 + tundraHeath.100 + water.1000 + snowWater.1000 + shrub.100 + meadow.500 + forest.500 + (1|CollarID),
                    data = dfRSF, family = binomial(link = "logit"))
modH[[15]] <- glmer(Response ~ slopeS.500 + vrmS.500 + elevS.30 + elevS2.30 +   # tundra + snow + water
                      aspectS.30 + tundraHeath.100 + water.1000 + snowWater.1000 + shrub.100 + meadow.500 + forest.500 + rock.1000 + (1|CollarID),
                    data = dfRSF, family = binomial(link = "logit"))

# Alert when completed:
beepr::beep()    # Beep :)

# Check models for overdispersion:
# source("./Analysis/OverdispersionFunction.R")
# ModHDisp <- lapply(modH, function(x)
#   overdisp_fun(x))
# No dispersion (p = 1)


# Model selection
require(AICcmodavg)

pred.base <- c("slopeS.500 + vrmS.1000 + elevS.30 + elevS2.30 +aspectS.30")
pred <- c("slopeS.500 + vrmS.1000 + elevS.30 + elevS2.30 +aspectS.30 + Forest.500",
          "slopeS.500 + vrmS.1000 + elevS.30 + elevS2.30 + aspectS.30 + Shrub.100",
          "slopeS.500 + vrmS.1000 + elevS.30 + elevS2.30 +aspectS.30 + Tundra.100",
          "slopeS.500 + vrmS.1000 + elevS.30 + elevS2.30 + aspectS.30 + Meadow.500",
          "slopeS.500 + vrmS.1000 + elevS.30 + elevS2.30 +aspectS.30 + Water.1000",
          "slopeS.500 + vrmS.1000 + elevS.30 + elevS2.30 +aspectS.30 + Snow.1000",
          "slopeS.500 + vrmS.1000 + elevS.30 + elevS2.30 +aspectS.30 + Rock.1000",
          "slopeS.500 + vrmS.1000 + elevS.30 + elevS2.30 +aspectS.30",
          "slopeS.500 + vrmS.1000 + elevS.30 + elevS2.30 +aspectS.30 + Tundra.100 + Snow.1000",
          "slopeS.500 + vrmS.1000 + elevS.30 + elevS2.30 +aspectS.30 + Tundra.100 + Water.1000",
          "slopeS.500 + vrmS.1000 + elevS.30 + elevS2.30 +aspectS.30 + Tundra.100 + Snow.1000 + Water.1000",
          "slopeS.500 + vrmS.1000 + elevS.30 + elevS2.30 +aspectS.30 + Tundra.100 + Snow.1000 + Water.1000 + shrub.100",
          "slopeS.500 + vrmS.1000 + elevS.30 + elevS2.30 +aspectS.30 + Tundra.100 + Snow.1000 + Water.1000 + Shrub.100 + Meadow.500",
          "slopeS.500 + vrmS.1000 + elevS.30 + elevS2.30 +aspectS.30 + Tundra.100 + Snow.1000 + Water.1000 + Shrub.100 + Meadow.500 + Forest.500",
          "slopeS.500 + vrmS.1000 + elevS.30 + elevS2.30 +aspectS.30 + Tundra.100 + Snow.1000 + Water.1000 + Shrub.100 + Meadow.500 + Forest.500 + Rock.1000")


# Create a model selection table:
modH$aic <- aictab(cand.set = modH[1:15], modnames = pred, sort = TRUE)
modH$aic <- as.data.frame(modH$aic)
modH$aic <- modH$aic[order(modH$aic$Delta_AICc), ]

# Save the candidate model set:
save(modH, file = "./Data/Rdata/ModH.Rdata", compress = "gzip")

library(MuMIn)
r.squaredGLMM(modH[[14]])

# Extract the best model:
best <- modH[[14]]

# Save the best model:
save(best, file = "./Data/Rdata/best.Rdata", compress = "gzip")





### STEP 4: SUMMARY TABLES from the top model

# Get confidence intervals from the top model to interpret results:
per1_se <- sqrt(diag(vcov(best)))

# table of parameter estimates with 95% CI:
tab_per1 <- cbind(Est = fixef(best), LL = fixef(best) - 1.96 * per1_se,
                  UL = fixef(best) + 1.96 * per1_se)






### STEP 5: GOODNESS OF FIT of the top model

# Add predicted values (based on the best model) to dfRSF:
dfRSF$pred <- exp(predict(best, newdata=dfRSF, re.form=~(1|CollarID), type="link"))

# Evaluate goodness of fit:
library(ResourceSelection)
hoslem.test(dfRSF$Response, Pred, 5)
