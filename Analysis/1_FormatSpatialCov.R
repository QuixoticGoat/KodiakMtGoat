
################################################################################
# Creates spatial covariates for mountain goat habitat selection at Kodiak
# Island, Alaska.
#
# Author: McCrea Cobb
# Created: 5/10/2017
# Last modified: 5/24/2018
################################################################################


## ---- 1_FormatSpatialCov.R


#-------------------------------------------------------------------------------
## STEP 1. Import and format spatial covariate data ###

# #Required packages:
packages <- c("rgdal", "maptools", "raster", "rasterVis", "sp")
lapply(packages, require, character.only = TRUE)


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




#-------------------------------------------------------------------------------
## STEP 2: Format LCC raster as binary present/absent (0/1)

Format.LCC <- function() {

  ## Create binomial rasters of each habitat class from the lcc raster. Save each
  ## as a .grd file.

  require(snow)
  require(parallel)
  require(raster)

  # Allow for parallel processing to speed things up:
  beginCluster(detectCores() - 1)

  lcc <- raster("Data/GIS/LandCover/lcc.grd")
  elevS <- raster("Data/GIS/Elevation/elevS.grd")
  elev2S <- raster("Data/GIS/Elevation/elev2S.grd")
  slopeS <- raster("Data/GIS/Slope/slopeS.grd")
  slope2S <- raster("Data/GIS/Slope/slope2S.grd")
  aspectS <- raster("Data/GIS/Aspect/aspectS.grd")
  vrmS <- raster("Data/GIS/VRM/vrmS.grd")
  sriS <- raster("Data/GIS/SRI/sriS.grd")

  # Create the rasters

  # Forest:
  message("Working on Forest layer")
  lcc1 <- lcc
  lcc1[lcc1 > 1] <- 0
  lcc1@data@names <- "forest"
  writeRaster(x = lcc1, filename = paste("Data/GIS/LandCover/Forest",
                                         ".grd", sep = ""),
              format = "raster",
              progress = "text",
              overwrite = T)
  message("Saved!")

  # Shrub:
  message("Working on Shrub layer")
  lcc2 <- lcc
  lcc2[lcc2 != 2] <- 0
  lcc2[lcc == 2] <- 1
  lcc2@data@names <- "shrub"
  writeRaster(x = lcc2, filename = paste("Data/GIS/LandCover/Shrub",
                                         ".grd", sep = ""),
              format = "raster",
              progress = "text",
              overwrite = T)
  message("Saved!")

  # Tundra/Heath:
  message("Working on Tundra/Heath layer")
  lcc3 <- lcc
  lcc3[lcc3 != 3] <- 0
  lcc3[lcc3 == 3] <- 1
  lcc3@data@names <- "tundraHeath"
  writeRaster(x = lcc3, filename = paste("Data/GIS/LandCover/TundraHeath",
                                         ".grd", sep = ""),
              format = "raster",
              progress = "text",
              overwrite = T)
  message("Saved!")

  # Meadow:
  message("Working on Meadow layer")
  lcc4 <- lcc
  lcc4[lcc4 != 4] <- 0
  lcc4[lcc4 == 4] <- 1
  lcc4@data@names <- "meadow"
  writeRaster(x = lcc4, filename = paste("Data/GIS/LandCover/Meadow",
                                         ".grd", sep = ""),
              format = "raster",
              progress = "text",
              overwrite = T)
  message("Saved!")

  # Water:
  message("Working on Water layer")
  lcc5 <- lcc
  lcc5[lcc5 != 5] <- 0
  lcc5[lcc5 == 5] <- 1
  lcc5@data@names <- "water"
  writeRaster(x = lcc5, filename = paste("Data/GIS/LandCover/Water",
                                         ".grd", sep = ""),
              format = "raster",
              progress = "text",
              overwrite = T)
  message("Saved!")

  # Snow/Water:
  message("Working on Snow layer")
  lcc6 <- lcc
  lcc6[lcc6 != 6] <- 0
  lcc6[lcc6 == 6] <- 1
  lcc6@data@names <- "snowWater"
  writeRaster(x = lcc6, filename = paste("Data/GIS/LandCover/SnowWater",
                                         ".grd", sep = ""),
              format = "raster",
              progress = "text",
              overwrite = T)
  message("Saved!")

  # Rock:
  message("Working on Rock layer")
  lcc7 <- lcc
  lcc7[lcc7 != 7] <- 0
  lcc7[lcc7 == 7] <- 1
  lcc7@data@names <- "rock"
  writeRaster(x = lcc7, filename = paste("Data/GIS/LandCover/Rock",
                                         ".grd", sep = ""),
              format = "raster",
              progress = "text",
              overwrite = T)
  message("Saved!")

  # Done. Free up cores and memory:
  endCluster()
  removeTmpFiles(h = 0)

  message("All done!")
}

# Run it:
Format.LCC()



#-------------------------------------------------------------------------------
## STEP 3: Save spatial covariates together as a raster stack

raster.stack <- function() {

  require(parallel)
  require(raster)
  require(rasterVis)

  # Allow for parallel processing to speed things up:
  beginCluster(detectCores() - 1)

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

  # Create a raster stack:
  r <- stack(elevS, elev2S, slopeS, slope2S, aspectS, vrmS, forest, shrub,
             tundraHeath, meadow, water, snowWater, rock)

  # Save the final raster stack as a grid (.grd) file:
  writeRaster(r, "./Data/GIS/RasterStack/Covar.grd",
              format = "raster",
              options = c("INTERLEAVE=BAND"),
              progress = "text",
              prj = T,
              overwrite = T)

  removeTmpFiles(h = 0)

}

# Run it:
raster.stack()


#-------------------------------------------------------------------------------
## STEP 4. Extract the focal mean pixel values for each covariate by 100, 500 
## and 1000 m.

library(raster)

# Load the covariate raster stack:
Covar.30 <- stack("./Data/GIS/RasterStack/Covar.grd")

# Define the focal weight (change as needed..)
fw <- focalWeight(Covar.30, 100, type="circle")

# Function to run focal() across all layers in a rasterstack:
multiFocal <- function(x, fw, ...) {
  library(raster)

  if(is.character(x)) {
    x <- brick(x)
  }

  # The function to be applied to each individual layer
  fun <- function(ind, x, w, ...){
    focal(x[[ind]], w=fw, na.rm=TRUE)
  }

  n <- seq(nlayers(x))
  list <- lapply(X=n, FUN=fun, x=x, w=fw, ...)

  out <- stack(list)
  return(out)
}

# Run it:
Covar.100  <- multiFocal(x=Covar.30, fw)


# Rename the layers for each:
names <- names(Covar.30)
names(Covar.30) <- paste0(names, ".30")
names(Covar.100) <- paste0(names, ".100")
names(Covar.500) <- paste0(names, ".500")
names(Covar.1000) <- paste0(names, ".1000")

# Save them:
writeRaster(Covar.30, "./Data/GIS/Rasterstack/Covar.30.grd",
            format = "raster",
            options = c("INTERLEAVE=BAND"),
            progress = "text",
            prj = T,
            overwrite = T)
writeRaster(Covar.100, "./Data/GIS/Rasterstack/Covar.100.grd",
            format = "raster",
            options = c("INTERLEAVE=BAND"),
            progress = "text",
            prj = T,
            overwrite = T)
writeRaster(Covar.500, "./Data/GIS/Rasterstack/Covar.500.grd",
            format = "raster",
            options = c("INTERLEAVE=BAND"),
            progress = "text",
            prj = T,
            overwrite = T)
writeRaster(Covar.1000, "./Data/GIS/Rasterstack/Covar.1000.grd",
            format = "raster",
            options = c("INTERLEAVE=BAND"),
            progress = "text",
            prj = T,
            overwrite = T)

# Combine into a single raster stack:
Covar.all <- stack(Covar.30, Covar.100, Covar.500, Covar.1000)

# Save it:
writeRaster(Covar.all, "./Data/GIS/Rasterstack/Covar.all.grd",
            format = "raster",
            options = c("INTERLEAVE=BAND"),
            progress = "text",
            prj = T,
            overwrite = T)


#-------------------------------------------------------------------------------
## STEP 5. Import the raster stack containing all the covariates and save
## a plot of it:

Covar.all <- stack("./Data/GIS/RasterStack/Covar.all.grd")

pdf("Plots/Maps/Covariates.pdf", width = 7, height = 4, title = "Covariates")
plot(Covar)
dev.off()

