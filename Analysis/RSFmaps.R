#################################################################
#### Create RSF maps from model output
#### Author: McCrea Cobb
#### Last modified: 9/15/2017
#################################################################

##---- RSFmaps.R

library(raster)

###############
# Load the best model and the raster stack of the spatial covariates
###############

Covar.all <- stack("./Data/GIS/RasterStack/Covar.all.grd")

# Load best summer model:
load("./Data/Rdata/mod_summer_final.Rdata")

# Load best winter mode:
load("./Data/Rdata/mod_winter_final.Rdata")


calc.RSF.surface <- function(mod, cov.stack=Covar.all) {
  # Function to calculate the RSF surface raster using the exponential equation
  # and save it.
                      # best = top model

  library(raster)
  library(scales)
  library(leaflet)
  library(rgdal)

  # Get the coefficients and their names from the best model:
  coefs <- summary(mod)$coef[2:nrow(summary(mod)$coef)]
  facts <- rownames(summary(mod)$coef)[2:nrow(summary(mod)$coef)]

  # calculate the RSF surface:
  out.RSF <- exp(sum(coefs * cov.stack[[facts]]))

  # Scale between 0 and 1:
  scaleRSF.fn <- function(out.RSF = out.RSF) {
    RSF <- (out.RSF - cellStats(out.RSF, "min", na.rm = TRUE))/(cellStats(out.RSF, "max", na.rm = TRUE) - cellStats(out.RSF, "min", na.rm = TRUE))
    return(RSF)
    }

  out.RSF.scaled <- scaleRSF.fn(out.RSF)


  return(out.RSF.scaled)
}

##############
# Run it
##############

RSFsurface <- calc.RSF.surface(best)


###############
# Save the summer and winter RSFs
###############

writeRaster(RSFsurface, "./Data/GIS/RSFsurface/rsfwinter.grd",
            format = "raster",
            options = c("INTERLEAVE=BAND"),
            progress = "text",
            prj = T,
            overwrite = T)

writeRaster(RSFsurface, "./Data/GIS/RSFsurface/rsfwinter.tif",
            format = "GTiff",
            options = c("INTERLEAVE=BAND"),
            progress = "text",
            prj = T,
            overwrite = T)




######################
# Bin the RSF surface
######################

quantile.bin <- function(x,n.bins = 10, n.transparent = round(n.bins * 0.6)) {
  q <- quantile(x,
                probs = seq(0.1, 0.9, by = 0.1),
                na.rm = TRUE)
  rcl.mat <- cbind(c(0, q[1:9]),
                   c(q[1:9],
                     cellStats(x,'max') + 0.1),
                   1:10)
  reclass.rast <- reclassify(x, rcl = rcl.mat, include.lowest = T)
  reclass.rast[reclass.rast <= n.transparent] <- NA  # Remove values >0.6
  return(reclass.rast)
}

RSFsurfaceBin <- quantile.bin(RSFsurface)

writeRaster(RSFsurfaceBin, "./Data/GIS/RSFsurface/rsfBinsummer7to10.grd",
            format = "raster",
            options = c("INTERLEAVE=BAND"),
            progress = "text",
            prj = T,
            overwrite = T)

writeRaster(RSFsurfaceBin, "./Data/GIS/RSFsurface/rsfBinsummer7to10.tif",
            format = "GTiff",
            options = c("INTERLEAVE=BAND"),
            progress = "text",
            prj = T,
            overwrite = T)



###########################
## Plot the RSF surface map
###########################
library(raster)
kodiak <- shapefile("./Data/GIS/KodiakBound/kodiak_island")


#### Summer
rsf <- raster("./Data/GIS/RSFsurface/rsfBinsummer7to10.grd")
jpeg(file = "./Output/Maps/RSFsummer.jpg", width = 1000, height = 1000)
plot(kodiak)
plot(rsf, col = rainbow(10, alpha = 1), add=T)
dev.off()

#### Winter
rsf <- raster("./Data/GIS/RSFsurface/rsfBinwinter7to10.grd")
jpeg(file = "./Output/Maps/RSFwinter.jpg", width = 1000, height = 1000)
plot(kodiak)
plot(rsf, col = rainbow(10, alpha = 1), add=T)
dev.off()

