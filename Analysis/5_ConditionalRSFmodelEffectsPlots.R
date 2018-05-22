###### Effects plots


# Create prediction dataset
rm(df.pr)
tundraHeath.100 <- shrub.100 <- water.1000 <- elevS.30 <- elevS2.30 <- vrmS.30 <- aspectS.30 <- slopeS.30 <- slopeS2.30 <- rep(0.5, 1000)

df.pr <- data.frame(tundraHeath.100, shrub.100, water.1000, elevS.30, elevS2.30, vrmS.30, aspectS.30, slopeS.30, slopeS2.30)
df.pr$case <- 1
df.pr$strata <- 1
df.pr$elevS.30 <- 1.1
df.pr$elevS2.30 <- df.pr$elevS.30^2

df.pr$shrub.100 <- seq(from = 0, to = 1, length.out = 1000)


## Elevation
########
df.pr$elevS.30 <- seq(from = -2, to = 3, length.out = 1000)
df.pr$elevS2.30 <- df.pr$elevS.30^2
########


df.pr$predicted <- predict(best, df.pr, type = "response")

plot(df.pr$predicted ~ df.pr$shrub.100)
