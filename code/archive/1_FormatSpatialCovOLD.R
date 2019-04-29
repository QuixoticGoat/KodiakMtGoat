## Create the spatial covariates for mountain goat habitat selection
## Author: McCrea Cobb
## Created: 5/10/2017

## ---- 1_FormatSpatialCov.R



### STEP 1. Import and format spatial covariate data ###

library(rgdal)
library(maptools)
library(raster)
library(rasterVis)
library(sp)

## Shapefile of study area boundary (Kodiak Island)
StudyArea <- readOGR(dsn = "Data/GIS/KodiakBound/kodiak_island.shp")

## Elevation, slope, aspect and
## VRM (vector ruggedness measure) rasters. Calculate elev and slope squared.
elev <- raster("Data/GIS/Elevation/dem.tif")
elev@data@names <- "Elev"
elev2 <- elev^2
elev2@data@names <- "Elev2"
slope <- raster("Data/GIS/Slope/slope.tif")
slope@data@names <- "Slope"
slope2 <- slope^2
slope2@data@names <- "Slope2"
aspect <- raster("Data/GIS/Aspect/aspect.tif")
aspect@data@names <- "Aspect"
vrm <- raster("Data/GIS/VRM/vrm.tif")
sri <- raster("Data/GIS/SRI/sri.tif")

## Read in land cover classification raster (7 classes):
lcc <- raster("Data/GIS/LandCover/lccNew1.tif")
# Set lcc extent to match elevation raster extent:
ex <- extent(elev)
lcc <- crop(lcc, ex)



# Calculate "southness":
aspect[aspect < 0] <- NA  # set flat values (-1) as NA
aspectS <- abs((abs(aspect-180)/180)-1)  # Calculate index (S = 1, N = 0)
aspectS@data@names <- "aspectS"
writeRaster(x = aspectS, filename = "Data/GIS/Aspect/aspectS.grd",
            format = "raster", progress = "text", overwrite = T)


# Standardize elev^2, slope, and vrm around zero +- 1 SD. Export as grid file.
require(raster)
library(scales)

elevS <- scale(elev)
elevS@data@names <- "elevS"
writeRaster(x = elevS, filename = "Data/GIS/Elevation/elevS.grd",
            format = "raster", progress = "text", overwrite = T)

elev2S <- scale(elev2)
elev2S@data@names <- "elevS2"
writeRaster(x = elev2S, filename = "Data/GIS/Elevation/elev2S.grd",
            format = "raster", progress = "text", overwrite = T)

slopeS <- scale(slope)
slopeS@data@names <- "slopeS"
writeRaster(x = slopeS, filename = "Data/GIS/Slope/slopeS.grd",
            format = "raster", progress = "text", overwrite = T)

slope2S <- scale(slope2)
slope2S@data@names <- "slopeS2"
writeRaster(x = slope2S, filename = "Data/GIS/Slope/slope2S.grd",
            format = "raster", progress = "text", overwrite = T)

vrmS <- scale(vrm)
vrmS@data@names <- "vrmS"
writeRaster(x = vrmS, filename = "Data/GIS/VRM/vrmS.grd",
            format = "raster", progress = "text", overwrite = T)

sriS <- scale(sri)
sriS@data@names <- "sriS"
writeRaster(x = sriS, filename = "Data/GIS/SRI/sriS.grd",
            format = "raster", progress = "text", overwrite = T)

# Resample lcc ("nearest neighbor" method) to elev:
lcc <- resample(lcc, elev, "ngb")  # TAKES A WHILE! (~5-10 min)
writeRaster(x = lcc, filename = "Data/GIS/LandCover/lcc.grd",
            format = "raster", progress = "text", overwrite = T)





### STEP 2: FORMAT SPATIAL COVARIATES AT DIFFERENT SPATIAL SCALES

## Create binomial rasters of each habitat class from the lcc raster,
## based on the mean value within 50m radius of each pixel cell.

# First, need to define focal weight (50 m) for moving window (focal()):

Format.LCC <- function(radius) {

  # Allow for parallel processing to speed things up:
  require(snow)
  require(parallel)
  beginCluster( detectCores() - 1)

  require(raster)

  lcc <- raster("Data/GIS/LandCover/lcc.grd")
  elevS <- raster("Data/GIS/Elevation/elevS.grd")
  elev2S <- raster("Data/GIS/Elevation/elev2S.grd")
  slopeS <- raster("Data/GIS/Slope/slopeS.grd")
  slope2S <- raster("Data/GIS/Slope/slope2S.grd")
  aspectS <- raster("Data/GIS/Aspect/aspectS.grd")
  vrmS <- raster("Data/GIS/VRM/vrmS.grd")
  sriS <- raster("Data/GIS/SRI/sriS.grd")

  # Define the focal weight (radius around fix location to average habitat values)
  fw <- focalWeight(lcc, radius, "circle")

  # Then, create the rasters:

  # Forest:
  message("Working on Forest layer")
  lcc1 <- lcc
  lcc1[lcc1 > 1] <- 0
  lcc1 <- focal(lcc1, w = fw, fun = "mean", na.rm = T) # mean value within "fw"
  lcc1@data@names <- "forest"
  writeRaster(x = lcc1, filename = paste("Data/GIS/LandCover/Forest",
                                         radius,".grd", sep = ""),
              format = "raster",
              progress = "text",
              overwrite = T)

  # Shrub:
  message("Working on Shrub layer")
  lcc2 <- lcc
  lcc2[lcc2 != 2] <- 0
  lcc2[lcc == 2] <- 1
  lcc2 <- focal(lcc2, w = fw, fun = "mean", na.rm = T)
  lcc2@data@names <- "shrub"
  writeRaster(x = lcc2, filename = paste("Data/GIS/LandCover/Shrub",
                                         radius,".grd", sep = ""),
              format = "raster",
              progress = "text",
              overwrite = T)

  # Tundra/Heath:
  message("Working on Tundra/Heath layer")
  lcc3 <- lcc
  lcc3[lcc3 != 3] <- 0
  lcc3[lcc3 == 3] <- 1
  lcc3 <- focal(lcc3, w = fw, fun = "mean", na.rm = T)
  lcc3@data@names <- "tundraHeath"
  writeRaster(x = lcc3, filename = paste("Data/GIS/LandCover/TundraHeath",
                                         radius,".grd", sep = ""),
              format = "raster",
              progress = "text",
              overwrite = T)

  # Meadow:
  message("Working on Meadow layer")
  lcc4 <- lcc
  lcc4[lcc4 != 4] <- 0
  lcc4[lcc4 == 4] <- 1
  lcc4 <- focal(lcc4, w = fw, fun = "mean", na.rm = T)
  lcc4@data@names <- "meadow"
  writeRaster(x = lcc4, filename = paste("Data/GIS/LandCover/Meadow",
                                         radius,".grd", sep = ""),
              format = "raster",
              progress = "text",
              overwrite = T)

  # Water:
  message("Working on Water layer")
  lcc5 <- lcc
  lcc5[lcc5 != 5] <- 0
  lcc5[lcc5 == 5] <- 1
  lcc5 <- focal(lcc5, w = fw, fun = "mean", na.rm = T)
  lcc5@data@names <- "water"
  writeRaster(x = lcc5, filename = paste("Data/GIS/LandCover/Water",
                                         radius,".grd", sep = ""),
              format = "raster",
              progress = "text",
              overwrite = T)

  # Snow/Water:
  message("Working on Snow layer")
  lcc6 <- lcc
  lcc6[lcc6 != 6] <- 0
  lcc6[lcc6 == 6] <- 1
  lcc6 <- focal(lcc6, w = fw, fun = "mean", na.rm = T)
  lcc6@data@names <- "snowWater"
  writeRaster(x = lcc6, filename = paste("Data/GIS/LandCover/SnowWater",
                                         radius,".grd", sep = ""),
              format = "raster",
              progress = "text",
              overwrite = T)

  # Rock:
  message("Working on Rock layer")
  lcc7 <- lcc
  lcc7[lcc7 != 7] <- 0
  lcc7[lcc7 == 7] <- 1
  lcc7 <- focal(lcc7, w = fw, fun = "mean", na.rm = T)
  lcc7@data@names <- "rock"
  writeRaster(x = lcc7, filename = paste("Data/GIS/LandCover/Rock",
                                         radius,".grd", sep = ""),
              format = "raster",
              progress = "text",
              overwrite = T)


  elevS.F <- focal(elevS, w = fw, fun = "mean", na.rm = T)
  elevS.F@data@names <- "Elevation"
  writeRaster(x = elevS.F, filename = paste("Data/GIS/Elevation/ElevS.F",
                                            radius,".grd", sep = ""),
              format = "raster",
              progress = "text",
              overwrite = T)


  # Save the final stacked raster as a grid:
  message("Saving raster stack")
  writeRaster(r, paste("./Data/GIS/RasterStack/LCC", radius, ".grd", sep = ""),
              format = "raster",
              options = c("INTERLEAVE=BAND"),
              progress = "text",
              prj = T,
              overwrite = T)

  # Done. Free up cores:
  endCluster()

  removeTmpFiles(h = 0)
  message("Done!")
}


Format.LCC(50)





## Save spatial covariates at a raster stack
## (One stack for each focal weight area)

## Open the final spatial covariate raster layers:
elevS <- raster("Data/GIS/Elevation/elevS.grd")
elev2S <- raster("Data/GIS/Elevation/elev2S.grd")
slopeS <- raster("Data/GIS/Slope/slopeS.grd")
slope2S <- raster("Data/GIS/Slope/slope2S.grd")
aspectS <- raster("Data/GIS/Aspect/aspectS.grd")
vrmS <- raster("Data/GIS/VRM/vrmS.grd")
sriS <- raster("Data/GIS/SRI/sriS.grd")
forest <- raster("Data/GIS/LandCover/Forest.grd")
shrub <- raster("Data/GIS/LandCover/Shrub.grd")
tundraHeath <- raster("Data/GIS/LandCover/TundraHeath.grd")
meadow <- raster("Data/GIS/LandCover/Meadow.grd")
water <- raster("Data/GIS/LandCover/Water.grd")
snowWater <- raster("Data/GIS/LandCover/SnowWater.grd")
rock <- raster("Data/GIS/LandCover/Rock.grd")

require(rasterVis)
r <- stack(elevS, elev2S, slopeS, slope2S, aspectS, vrmS, forest, shrub,
           tundraHeath, meadow, water, snowWater, rock)

# Save the final stacked raster as a grid:
writeRaster(r, "./Data/GIS/RasterStack/Covar.grd",
            format = "raster",
            options = c("INTERLEAVE=BAND"),
            progress = "text",
            prj = T,
            overwrite = T)



# Import the raster stack containing all the covariates:
Covar <- stack("./Data/GIS/RasterStack/Covar.grd")




## Save a plot of the final raster stack:
pdf("Plots/Maps/Covariates.pdf", width = 7, height = 4, title = "Covariates")
plot(Covar)
dev.off()


# Clean up:
rm(r, elevS, elev2S, slopeS, slope2S, aspectS, vrmS, forest, shrub, tundraHeath, meadow, water,
   snowWater, rock, elev, elev2, aspect, aspect30,aspectCos, lcc, lccJunk, ned, slope,
   slope2, vrm)

# Free up memory:
removeTmpFiles(h = 0)
