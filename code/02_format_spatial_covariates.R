
################################################################################
# Creates spatial covariates for mountain goat habitat selection at Kodiak
# Island, Alaska.
#
# Author: McCrea Cobb
# Created: 5/10/2017
# Last modified: 4/29/2019
################################################################################


## ---- 1_FormatSpatialCov.R


#-------------------------------------------------------------------------------
## STEP 1. Import and format spatial covariate data ###

# #Required packages:
packages <- c("rgdal", "maptools", "raster", "rasterVis", "sp")
lapply(packages, require, character.only = TRUE)


## Shapefile of study area boundary (Kodiak Island)
StudyArea <- readOGR(dsn = "resources/geodata/KodiakBound/kodiak_island.shp")

## Elevation, slope, aspect and
## VRM (vector ruggedness measure) rasters. Calculate elev and slope squared.
elev <- raster("resources/geodata/Elevation/dem.tif")
elev@data@names <- "Elev"
elev2 <- elev^2
elev2@data@names <- "Elev2"
slope <- raster("resources/geodata/Slope/slope.tif")
slope@data@names <- "Slope"
slope2 <- slope^2
slope2@data@names <- "Slope2"
aspect <- raster("resources/geodata/Aspect/aspect.tif")
aspect@data@names <- "Aspect"
vrm <- raster("resources/geodata/VRM/vrm.tif")
sri <- raster("resources/geodata/SRI/sri.tif")

## Read in land cover classification raster (7 classes):
lcc <- raster("resources/geodata/LandCover/lccNew1.tif")
# Set lcc extent to match elevation raster extent:
ex <- extent(elev)
lcc <- crop(lcc, ex)



# Calculate "southness":
aspect[aspect < 0] <- NA  # set flat values (-1) as NA
aspectS <- abs((abs(aspect-180)/180)-1)  # Calculate index (S = 1, N = 0)
aspectS@data@names <- "aspectS"
writeRaster(x = aspectS, filename = "data/derived/geodata/aspectS.grd",
            format = "raster", progress = "text", overwrite = T)


# Standardize elev^2, slope, and vrm around zero +- 1 SD. Export as grid file.
require(raster)
library(scales)

elevS <- scale(elev)
elevS@data@names <- "elevS"
writeRaster(x = elevS, filename = "data/derived/geodata/elevS.grd",
            format = "raster", progress = "text", overwrite = T)

elev2S <- scale(elev2)
elev2S@data@names <- "elevS2"
writeRaster(x = elev2S, filename = "data/derived/geodata/elev2S.grd",
            format = "raster", progress = "text", overwrite = T)

slopeS <- scale(slope)
slopeS@data@names <- "slopeS"
writeRaster(x = slopeS, filename = "data/derived/geodata/slopeS.grd",
            format = "raster", progress = "text", overwrite = T)

slope2S <- scale(slope2)
slope2S@data@names <- "slopeS2"
writeRaster(x = slope2S, filename = "data/derived/geodata/slope2S.grd",
            format = "raster", progress = "text", overwrite = T)

vrmS <- scale(vrm)
vrmS@data@names <- "vrmS"
writeRaster(x = vrmS, filename = "data/derived/geodata/vrmS.grd",
            format = "raster", progress = "text", overwrite = T)

sriS <- scale(sri)
sriS@data@names <- "sriS"
writeRaster(x = sriS, filename = "data/derived/geodata/sriS.grd",
            format = "raster", progress = "text", overwrite = T)

# Resample lcc ("nearest neighbor" method) to elev:
lcc <- resample(lcc, elev, "ngb")  # TAKES A WHILE! (~5-10 min)
writeRaster(x = lcc, filename = "data/derived/geodata/lcc.grd",
            format = "raster", progress = "text", overwrite = T)




#-------------------------------------------------------------------------------
## STEP 2: Format LCC raster as binary present/absent (0/1)

format.lcc <- function() {

  ## Create binomial rasters of each habitat class from the lcc raster. Save each
  ## as a .grd file.

  rm(list=ls())  # Clear working directory to speed things up
  
  require(snow)
  require(parallel)
  require(raster)

  # Allow for parallel processing to speed things up:
  beginCluster(detectCores() - 1)

  lcc <- raster("data/derived/geodata/lcc.grd")
  elevS <- raster("data/derived/geodata/elevS.grd")
  elev2S <- raster("data/derived/geodata/elev2S.grd")
  slopeS <- raster("data/derived/geodata/slopeS.grd")
  slope2S <- raster("data/derived/geodata/slope2S.grd")
  aspectS <- raster("data/derived/geodata/aspectS.grd")
  vrmS <- raster("data/derived/geodata/vrmS.grd")
  sriS <- raster("data/derived/geodata/sriS.grd")

  # Create the rasters

  # Forest:
  message("Working on Forest layer...")
  lcc1 <- lcc
  lcc1[lcc1 > 1] <- 0
  lcc1@data@names <- "forest"
  writeRaster(x = lcc1, filename = paste("data/derived/geodata/forest",
                                         ".grd", sep = ""),
              format = "raster",
              progress = "text",
              overwrite = T)
  message("...Saved!")

  # Shrub:
  message("Working on Shrub layer...")
  lcc2 <- lcc
  lcc2[lcc2 != 2] <- 0
  lcc2[lcc == 2] <- 1
  lcc2@data@names <- "shrub"
  writeRaster(x = lcc2, filename = paste("data/derived/geodata/shrub",
                                         ".grd", sep = ""),
              format = "raster",
              progress = "text",
              overwrite = T)
  message("Saved!")

  # Tundra/Heath:
  message("Working on Tundra/Heath layer...")
  lcc3 <- lcc
  lcc3[lcc3 != 3] <- 0
  lcc3[lcc3 == 3] <- 1
  lcc3@data@names <- "tundraHeath"
  writeRaster(x = lcc3, filename = paste("data/derived/geodata/tundraheath",
                                         ".grd", sep = ""),
              format = "raster",
              progress = "text",
              overwrite = T)
  message("...Saved!")

  # Meadow:
  message("Working on Meadow layer...")
  lcc4 <- lcc
  lcc4[lcc4 != 4] <- 0
  lcc4[lcc4 == 4] <- 1
  lcc4@data@names <- "meadow"
  writeRaster(x = lcc4, filename = paste("data/derived/geodata/meadow",
                                         ".grd", sep = ""),
              format = "raster",
              progress = "text",
              overwrite = T)
  message("...Saved!")

  # Water:
  message("Working on Water layer...")
  lcc5 <- lcc
  lcc5[lcc5 != 5] <- 0
  lcc5[lcc5 == 5] <- 1
  lcc5@data@names <- "water"
  writeRaster(x = lcc5, filename = paste("data/derived/geodata/water",
                                         ".grd", sep = ""),
              format = "raster",
              progress = "text",
              overwrite = T)
  message("...Saved!")

  # Snow/Water:
  message("Working on Snow layer...")
  lcc6 <- lcc
  lcc6[lcc6 != 6] <- 0
  lcc6[lcc6 == 6] <- 1
  lcc6@data@names <- "snowWater"
  writeRaster(x = lcc6, filename = paste("data/derived/geodata/snowWater",
                                         ".grd", sep = ""),
              format = "raster",
              progress = "text",
              overwrite = T)
  message("...Saved!")

  # Rock:
  message("Working on Rock layer...")
  lcc7 <- lcc
  lcc7[lcc7 != 7] <- 0
  lcc7[lcc7 == 7] <- 1
  lcc7@data@names <- "rock"
  writeRaster(x = lcc7, filename = paste("data/derived/geodata/rock",
                                         ".grd", sep = ""),
              format = "raster",
              progress = "text",
              overwrite = T)
  message("...Saved!")

  # Done. Free up cores and memory:
  endCluster()
  removeTmpFiles(h = 0)

  message("All done!")
}

# Run it:
format.lcc()



#-------------------------------------------------------------------------------
## STEP 3: Save spatial covariates together as a raster stack

raster.stack <- function() {

  require(parallel)
  require(raster)
  require(rasterVis)

  # Allow for parallel processing to speed things up:
  beginCluster(detectCores() - 1)

  elevS <- raster("data/derived/geodata/elevS.grd")
  elev2S <- raster("data/derived/geodata/elev2S.grd")
  slopeS <- raster("data/derived/geodata/slopeS.grd")
  slope2S <- raster("data/derived/geodata/slope2S.grd")
  aspectS <- raster("data/derived/geodata/aspectS.grd")
  vrmS <- raster("data/derived/geodata/vrmS.grd")
  sriS <- raster("data/derived/geodata/sriS.grd")
  forest <- raster("data/derived/geodata/forest.grd")
  shrub <- raster("data/derived/geodata/shrub.grd")
  tundraHeath <- raster("data/derived/geodata/tundraheath.grd")
  meadow <- raster("data/derived/geodata/meadow.grd")
  water <- raster("data/derived/geodata/water.grd")
  snowWater <- raster("data/derived/geodata/snowwater.grd")
  rock <- raster("data/derived/geodata/rock.grd")

  # Create a raster stack:
  r <- stack(elevS, elev2S, slopeS, slope2S, aspectS, vrmS, forest, shrub,
             tundraHeath, meadow, water, snowWater, rock)

  # Save the final raster stack as a grid (.grd) file:
  writeRaster(r, "./data/derived/geodata/covar.grd",
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
covar <- stack("./data/derived/geodata/covar.grd")



# Function to run focal() across all layers in a rasterstack:
multifocal <- function(x, focw, ...) {
  library(raster)

  # Define the focal weight (change as needed..)
  fw <- focalWeight(x, focw, type="circle")
  
  if(is.character(x)) {
    x <- brick(x)
  }

  # The function to be applied to each individual layer
  fun <- function(ind, x, w, ...){
    focal(x[[ind]], w=fw, na.rm=TRUE)
  }

  n <- seq(nlayers(x))
  message(paste("Calculating the focal area within a", focw, "m radius..."))
  list <- lapply(X=n, FUN=fun, x=x, w=fw)
  
  names <- names(x)
  names(list) <- paste0(names, ".", focw)
  
  out <- stack(list)
  
  # message("Saving the output as a rasterstack...")
  # writeRaster(out, paste0("./data/derived/geodata/covar.", focw, ".grd"),
  #             format = "raster",
  #             options = c("INTERLEAVE=BAND"),
  #             progress = "text",
  #             prj = T,
  #             overwrite = T)
  message("....Done!")
  return(out)
}

# Run it:
function() {
  fcs <- list(50, 100)
  
  covar.all <- lapply(fcs, function(x) {
    ls <- list()
    ls[x] <- multifocal(covar, x)
  })
  names(covar.all) <- paste0("covar.", fcs)
  covar.all <- stack(covar.all)
  # Save it:
  writeRaster(covar.all, "./data/derived/geodata/covar.all.grd",
              format = "raster",
              options = c("INTERLEAVE=BAND"),
              progress = "text",
              prj = T,
              overwrite = T)
  return(covar.all)
}






#-------------------------------------------------------------------------------
## STEP 5. Import the raster stack containing all the covariates and save
## a plot of it:

covar.all <- stack("./data/derived/geodata/covar.all.grd")

pdf("output/maps/Covariates.pdf", width = 7, height = 4, title = "Covariates")
plot(covar)
dev.off()

