
###################################################################
#### Movement modeling of GPS collar data
#### Author: McCrea Cobb
#### Last modified: 5/10/2017
###################################################################

##---- Movement.R


## (First run "ImportFormat" script)


# Required libraries
library(adehabitatLT)
library(movementAnalysis)
library(BBMM)

junk <- split(df, df$CollarID)

## Create "ltraj" object (trajectory object)
# dfIG18 <- subset(df, df$CollarID == "IG18")
DatTraj <- lapply(junk, function(x) {
  as.ltraj(xy = x[ ,c("Long", "Lat")],
                    date = x$Date,
                    id = "CollarID")})

head(DatTraj[[1]])
plot(DatTraj)
DatTraj2 <- ld(DatTraj)  # convert to a data frame

# Check to see if the relocations are regular (i.e. missed fixes)
is.regular(DatTraj)
# Plot time interval between successive relocations (dt):
plot.new()
plotltr(DatTraj, "dt")

# Insert missing values (NAs) to make time interval b/w fixes regular
refda <- strptime("00:00", "%H")
DatTrajReg <- setNA(DatTraj, refda, 12.5,units = "hour")

# Check whether each animal was monitored for the same duration
is.sd(DatTrajReg)
# Plot time interval between successive relocations (dt):
plot.new()
plotltr(DatTrajReg, "dt")

DatTrajRegLim <- set.limits(DatTrajReg, begin = "2015-07-26 02:00",
                         dur = 3360, units = "hour", pattern = "%Y-%m-%d %H:%M")

runsNAltraj(DatTrajRegLim)  # Check whether NAs are randomly distributed over time




# Wald-Wolfowitx test of randomness (test of sequential autocorrelation)
wawotest(DatTrajRegLim)



## Movement functions

move.sum <- function(df.traj) {
  ## Summarize the mean, min and max distances traveled, by CollarID.

  CollarID <- sapply(df.traj, function(x) {
    attr(x, "id")
  })
  DistMean <- sapply(df.traj, function(x) {
    round(mean(x$dist, na.rm = T), 2)
  })
  DistMax <- sapply(df.traj, function(x) {
    round(max(x$dist, na.rm = T), 2)
  })
  DistMin <- sapply(df.traj, function(x) {
    round(min(x$dist, na.rm = T), 2)
  })
  DistSD <- sapply(df.traj, function(x) {
    round(sd(x$dist, na.rm = T), 2)
  })
  DistMed <- sapply(df.traj, function(x) {
    round(median(x$dist, na.rm = T), 2)
  })
  df <- cbind(CollarID, DistMean, DistMed, DistSD, DistMin, DistMax)
  df <- as.data.frame(df)

  return(df)
}

dist.sum <- move.sum(df.traj)


## Plots of the step lengths:
plotEm <- function() {
  library(ggplot2)

  # Convert from list to df:
  df.trajDF <- do.call("rbind", df.traj)

  p <- ggplot(df.traj[[4]], aes(x = dist)) +
    geom_histogram(binwidth = 50,
                   colour = "black",
                   fill = "white")

  # facet_grid(CollarID)
  p
}

plotEm()
