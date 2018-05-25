
## Extract the covariate values from the raster to used and available fixes
## McCrea Cobb
## 5/17/2017

## ---- 2_ExtractCov.R

library(rgdal)
library(beepr)
load("./Data/Rdata/df.Rdata")
load("./Data/Rdata/dfSp.Rdata")
load("./Data/GIS/KodiakBound/StudyArea.Rdata")

### STEP 1: Define "available" points at the kernel home range scale

## Create a 99% kernel HR from all "used" fixes
require(adehabitatHR)
require(rasterVis)  # for crs()

crs(dfSp) <- crs(StudyArea)
Kernel <- kernelUD(SpatialPoints(dfSp), h = "href")
crs(Kernel) <- crs(dfSp)
KernelVert <- getverticeshr(Kernel, percent = 99, unin = "m", unout = "km2")
# Crop the kernel to Kodiak Island (i.e., no ocean)
KernelVert <- crop(KernelVert, StudyArea)

## Generate a random sample of "available" points within the kernel:
require(sp)
AvailS <- spsample(KernelVert, type = "random", n = length(df$CollarID)*10)

# Reproject available points to WGS84:
AvailS <- spTransform(AvailS, CRS = CRS("+proj=longlat +datum=WGS84"))

# Convert available points spatial dataframe to a df, add CollarID values:
Avail <- data.frame(Response = 0, Lat = AvailS@coords[,2],
                      Long = AvailS@coords[,1])
Avail$CollarID <- rep(df$CollarID, 10)  # Add CollarID
Avail$Date <- rep(df$Date, 10)  # Add CollarID
Avail$Collar <- rep(df$Collar, 10)  # Add CollarID
Avail$sex <- rep(df$sex, 10)  # Add CollarID
Avail$season <- rep(df$season, 10)  # Add season




## Bind used and available dfs together:
common_cols <- intersect(colnames(df), colnames(Avail))
dfRSF <- rbind(subset(df[common_cols]),
               subset(Avail[common_cols]))

# Create a spatial points dataframe (dfRSFs) of dfRSF:
dfRSFs <- SpatialPointsDataFrame(dfRSF[ ,c("Long", "Lat")], dfRSF,
                                 proj4string = CRS("+proj=longlat +datum=WGS84"))
dfRSFs <- spTransform(dfRSFs, CRS("+proj=utm +zone=5 +datum=NAD83"))

# Clean up:
rm(common_cols, AvailS, Avail, KernelVert, Kernel)




### STEP 2: Extract spatial covariate raster values to used and available points:

require(rgeos)

## Import the raster stack:
Covar.all <- stack("./Data/GIS/RasterStack/Covar.all.grd")

# Extract the covariate data to the used/available spatial dataframe (dfRSFs):

dfTemp <- extract(Covar.all, dfRSFs[, c("Long", "Lat")],
               method = "simple",
               df = T) ; beepr::beep(4)

# Bind the extracted data to the response and lat/long df:
dfRSF <- cbind(dfRSF, dfTemp)
dfRSF$Response <- as.integer(dfRSF$Response)



# Add weights (# fixes/CollarID)
library(dplyr)

junk <- dfRSF %>%
  group_by(CollarID) %>%
  summarise(count = length(CollarID))
range01.fn <- function(x, ...){(x - max(x, ...)) / (min(x, ...) - max(x, ...))}
junk$weight <- range01.fn(junk$count)
junk$count <- NULL

dfRSF <- merge(dfRSF, junk)
rm(junk)


# Save it:
save(dfRSF, file="./Data/Rdata/dfRSF.Rdata", compress=T)





## Subset winter and summer
dfRSF.w <- subset(dfRSF, season == "winter")

# Add weights (# fixes/CollarID)
library(dplyr)

junk <- dfRSF.w %>%
  group_by(CollarID) %>%
  summarise(count = length(CollarID))
range01.fn <- function(x, ...){(x - max(x, ...)) / (min(x, ...) - max(x, ...))}
junk$weight <- range01.fn(junk$count)
junk$count <- NULL

dfRSF.w <- merge(dfRSF.w, junk)
rm(junk)

save(dfRSF.w, file="./Data/Rdata/dfRSF.w.Rdata", compress=T)




dfRSF.s <- subset(dfRSF, season == "summer")

# Add weights (# fixes/CollarID)
library(dplyr)

junk <- dfRSF.s %>%
  group_by(CollarID) %>%
  summarise(count = length(CollarID))
range01.fn <- function(x, ...){(x - max(x, ...)) / (min(x, ...) - max(x, ...))}
junk$weight <- range01.fn(junk$count)
junk$count <- NULL

dfRSF.s <- merge(dfRSF.s, junk)
rm(junk)

save(dfRSF.s, file="./Data/Rdata/dfRSF.s.Rdata", compress=T)



# Create a spatial points dataframe of the results and save as shapefile:
dfRSFs <- SpatialPointsDataFrame(dfRSF[ ,c("Long", "Lat")], dfRSF,
                                 proj4string = CRS("+proj=longlat +datum=WGS84"))
dfRSFs <- spTransform(dfRSFs, CRS("+proj=utm +zone=5 +datum=NAD83"))
writeOGR(obj = dfRSFs,
         dsn = "./Data/GIS/dfRSFs",
         layer = "dfRSFs",
         driver = "ESRI Shapefile",
         overwrite_layer = T)
save(dfRSFs, file="./Data/Rdata/dfRSFs.Rdata", compress=T)

## Clean up
rm(dfTemp, dfRSFs, Avail, AvailS)
