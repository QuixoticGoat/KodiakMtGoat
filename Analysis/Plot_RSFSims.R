########################################
## Plot the fixed effects of a RSF model
## Author: McCrea Cobb
## Last updated: 10/20/2017
########################################



library(merTools)

#  Simulate (bootstrapped) fixed effects posterior distributions of the top model
best.fe <- FEsim(best, 1000, oddsRatio = F)

# Plot it:
plotFEsim(data=best.fe, level=0.95, stat="median", sd=T, intercept=F)



# Predicted values (winter)
pred.w <- c("Response", "CollarID", "slopeS.100", "slopeS2.100", "vrmS.30",
            "elevS.30", "elevS2.30", "aspectS.30", "rock.500", "tundraHeath.30",
            "water.100", "forest.500", "shrub.100", "meadow.100")

new <- dfRSF.w[complete.cases(dfRSF.w), ]
new <- new[pred.w]

new$intFit <- predictInterval(best, new, type = "probability")



# Plots of predictions of relative probability of use for the top covariates
# in the final RSF model:
shinyMer(best, simData = new, pos=1)


# Predicted values (summer)
pred.s <- c("Response", "CollarID", "slopeS.500", "slopeS2.500", "vrmS.500",
            "elevS.30", "elevS2.30", "aspectS.30", "rock.500", "tundraHeath.30",
            "water.1000", "forest.500", "shrub.500", "meadow.1000")

new <- dfRSF.s[complete.cases(dfRSF.s), ]
new <- new[pred.s]

intFit <- predictInterval(best, new, type = "probability")
new <- cbind(new, intFit)


## Plots of marginal effects of the best model:
library(sjPlot)
library(sjmisc)
library(sjlabelled)

sjp.glm(best, type = "slope")

library(ggplot2)
ggplot(new, aes(x = slopeS.500, y = fit, ymin = lwr, ymax = upr)) +
  geom_ribbon()


#################################
#################################

library(merTools)

## Creates a dataframe from which a range of values for a chosen covariate,
## while keeping all other covariates at the mean value

# Calculate the mean values for each covariate:
df.p <- colMeans(dfRSF.w[pred.w[3:length(pred.w)]], na.rm = T)
df.p <- as.data.frame(t(df.p))

df <- wiggle(data = df.p, "slopeS.100",
             seq(from = min(dfRSF.w$slopeS.100, na.rm = T),
                 to = max(dfRSF.w$slopeS.100, na.rm = T),
                 by = 0.001))
df$CollarID <- "IG02"

fit <- predictInterval(best, df, type = "probability", which="fixed")

df <- cbind(df, fit)

library(ggplot2)

ggplot(df, aes(x = slopeS.100, y = fit, ymin = lwr, ymax = upr)) +
  geom_line() +
  geom_ribbon(alpha = 0.3)




