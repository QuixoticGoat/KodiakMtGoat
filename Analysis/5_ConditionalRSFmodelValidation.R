######################################################################
### Model Validation

## Use the testing data to validate the accuracy of the
## conditional logistic model by evaluating the correlation
## between the frequency of occurrance of the testing data
## and the relative RSF score using Spearman's Rank correlation
## statistic (Boyce et al. 2002).
#####################################################################


# Required packages
library(raster)


#######################
####### Summer ########
#######################

# Load the testing data
load(file="./Data/RData/df.ua.RData")
df <- df.ua.RData$Test

# Load the binned (10 bins) RSF:
rsf <- raster("./Data/GIS/RSFsurface/rsfBinsummer.grd")
plot(rsf)

# Subset to summer locations:
df <- subset(df, season == "summer")

# Convert df to a spatial dataframe:
df.sp <- SpatialPointsDataFrame(coords = data.frame(df$utmE, df$utmN), data = df,
                                   proj4string = CRS("+proj=utm +zone=5 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))

# clip the binned RSF values to the testing data
rsfValues <- extract(rsf, df.sp, method = "simple", buffer = 0)
rsfValues <- unlist(rsfValues)
df$rsfbin <- rsfValues

# Spearman Rank correlation:
df.u <- subset(df, case == 1) # subset used locations
foo <- data.frame(table(df.u$rsfbin))
foo$Var1 <- as.numeric(foo$Var1)
(rsf.spearman <- cor.test(foo$Var1, foo$Freq, method = "spearman"))

# Save it:
save(rsf.spearman, file = "./Output/SpearmanRanks/SpearmanSummer.RData")

# Plot a histogram
hist(df.a$rsfbin)





####################
##### Winter #######
####################

# Load the testing data
load(file="./Data/RData/df.ua.RData")
df <- df.ua.RData$Test

# Load the binned (10 bins) RSF:
rsf <- raster("./Data/GIS/RSFsurface/rsfBinsummer.grd")
plot(rsf)

# Subset to summer locations:
df <- subset(df, season == "winter")

# Convert df to a spatial dataframe:
df.sp <- SpatialPointsDataFrame(coords = data.frame(df$utmE, df$utmN), data = df,
                                proj4string = CRS("+proj=utm +zone=5 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))

# clip the binned RSF values to the testing data
rsfValues <- extract(rsf, df.sp, method = "simple", buffer = 0)
rsfValues <- unlist(rsfValues)
df$rsfbin <- rsfValues

# Spearman Rank correlation:
df.u <- subset(df, case == 1) # subset used locations
foo <- data.frame(table(df.u$rsfbin))
foo$Var1 <- as.numeric(foo$Var1)
(rsf.spearman <- cor.test(foo$Var1, foo$Freq, method = "spearman"))

# Save it:
save(rsf.spearman, file = "./Output/SpearmanRanks/SpearmanSummer.RData")

# Plot a histogram
hist(df.a$rsfbin)
