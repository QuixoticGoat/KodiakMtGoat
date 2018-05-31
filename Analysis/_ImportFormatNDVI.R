
################################################################################
# Imports MODIS-based NDVI and snow indices. Formats them as rasterstacks.     #
# Reprojects into UTM NAD83. Adds a date value (z) for each raster layer in    #
# the stack. Saves it as a .grd.                                               #
#                                                                              #
# NDVI/snow data source: <https://lpdaacsvc.cr.usgs.gov/appeears/>             #
#                                                                              #
# Author: McCrea Cobb <mccrea_cobb@fws.gov>                                    #
# Last modified: 5/30/2018                                                     #
################################################################################


# Required packages:
packages <- c("raster", "rgdal", "MODIStsp")
lapply(packages, require, character.only = TRUE)


#-------------------------------------------------------------------------------
## Import NDVI geotiffs and save them as a rasterstack
## Data were acquired on 3/23/18 using USGS AppEEARS. They are rasters of 
## NDVIs for Kodiak Island at a 16 day time interval, from 2013-2018.

ImportNDVI <- function() {
  # Assign a path to the object:
  ndvi_path <- "./Data/GIS/MODIS/ndvi"

  # Create a list of the file names:
  ndvi_files <- list.files(ndvi_path,
                           full.names = TRUE,
                           pattern = ".tif$")

  # Read in the geotiffs into a raster stack (into memory):
  message("Reading it in as a raster stack..")
  ndvi.250 <- stack(ndvi_files)
  message("Converting to a brick..")
  ndvi.250 <- brick(ndvi.250)  # Convert to a brick to speed up reprojection

  # Reproject to UTM NAD83
  message("Reprojecting into UTM NAD83..")
  ndvi.250 <- projectRaster(ndvi.250, res = c(250, 250), 
                            crs = "+proj=utm +zone=5 +datum=NAD83")

  ## Set a z (date) value for each raster in the stack:
  names <- names(ndvi.250)
  names <- substr(names, 35, 41) # extract the date values from the name
  ## Create a POSIXct
  names <- as.Date(names, format = "%Y%j")
  ndvi.250 <- setZ(ndvi.250, names)

  # Save it:
  message("Saving as a grid file..")
  writeRaster(ndvi.250, "./Data/GIS/MODIS/ndvi/Rasterstack/ndvi.250.grd",
              format = "raster",
              options = c("INTERLEAVE=BAND"),
              progress = "text",
              prj = T,
              overwrite = T)
  message("Done.")
  return(ndvi.250)
}

ndvi.250 <- ImportNDVI()




#-------------------------------------------------------------------------------
## Extract the focal mean pixel values for each covariate by 250, 500, 
## and 1,000 m.


 # 1. Create raster stacks for each focal mean buffer distance (500, 1000):

multiFocal <- function(x, fw, ...) {
  # Function to run focal() across all layers in a raster brick
  # x = a raster stack in UTM NAD83, fw = numeric value of buffer distance (m)

  require(raster)

  # Define the focus weight (i.e. buffer):
  fw <- focalWeight(x, fw, type = "circle")

  # Focal weight function:
  fun <- function(ind, x, w, ...) {
    focal(x[[ind]], w = fw, na.rm = TRUE)
  }

  # Apply the focal weight function to each raster in the brick:
  n <- seq(nlayers(x))
  list <- lapply(
    X = n,
    FUN = fun,
    x = x,
    w = fw,
    ...
  )

  # Convert result to a stack:
  out <- stack(list)

  # Match the output raster names with the input:
  names <- names(x)
  names(out) <- paste0(names)

  return(out)
}


# **Did not run at 250 m because the raw raster is 250m pixels.**

# Run it
ndvi.500  <- multiFocal(ndvi.250, 500)

# Save it:
writeRaster(ndvi.500, "./Data/GIS/MODIS/ndvi/Rasterbrick/ndvi.500.grd",
            format = "raster",
            options = c("INTERLEAVE=BAND"),
            progress = "text",
            prj = T,
            overwrite = T)

# Run it:
ndvi.1000  <- multiFocal(ndvi.250, 1000)

# Save it:
writeRaster(ndvi.1000, "./Data/GIS/MODIS/ndvi/Rasterbrick/ndvi.1000.grd",
            format = "raster",
            options = c("INTERLEAVE=BAND"),
            progress = "text",
            prj = T,
            overwrite = T)



#-------------------------------------------------------------------------------
## Get monthly composites of NDVI values (raw data are 16 day intervals)

## Load the ndvi rasterstacks at 250, 500, and 1000 m buffers:
ndvi.250 <- stack("./Data/GIS/MODIS/ndvi/Rasterbrick/ndvi.250.grd")
ndvi.500 <- stack("./Data/GIS/MODIS/ndvi/Rasterbrick/ndvi.500.grd")
ndvi.1000 <- stack("./Data/GIS/MODIS/ndvi/Rasterbrick/ndvi.1000.grd")
ndvi.ls <- list("ndvi.250" = ndvi.250, 
                "ndvi.500" = ndvi.500, 
                "ndvi.1000" = ndvi.1000)

MonthlyNDVI <- function(r) {
  # r = a rasterstack of ndvi values
  # Get the date from the names of the layers and extract the YearMonth
  indices <- format(as.Date(substr(names(r), 35, 41), format = "%Y%j"), 
                    format = "%Y%m")
  indices <- as.numeric(indices)
  
  message("Creating a brick with monthly means...")
  # Create a new rasterbrick with monthly mean NDVI values:
  ndvi.m <- stackApply(r, indices, fun = mean)
  
  # Add names to each raster in the rasterbrick ("ndvi_YearMonth")
  indices <- unique(indices)
  names(ndvi.m) <- paste0("ndvi_", indices)
  
  # Convert to a rasterstack:
  ndvi.m <- stack(ndvi.m)
  
  # Convert the indice to a date value (the 15th of every month) and
  # add as a z-value (not needed)
  # indices <- paste0(as.character(indices), "15")
  # indices <- as.Date(indices, format = "%Y%m%d")
  # ndvi.250.m <- setZ(ndvi.250.m, indices) 
  message("Saving it..")
  filename = paste("./Data/GIS/MODIS/ndvi/Rasterbrick/", deparse(substitute(r)), ".m.grd", sep="")
  # Save it:
  writeRaster(x = ndvi.m, 
              filename = filename,
              format = "raster",
              options = c("INTERLEAVE=BAND"),
              progress = "text",
              prj = T,
              overwrite = T)
  message("Done!")
}


lapply(ndvi.ls, MonthlyNDVI)

MonthlyNDVI(ndvi.250)
MonthlyNDVI(ndvi.500)
MonthlyNDVI(ndvi.1000)

#-------------------------------------------------------------------------------
## Extract data from each matching NDVI to each GPS collar fix


library(MODIStsp)


## Load the ndvi:
ndvi.250.m <- stack("./Data/GIS/MODIS/ndvi/Rasterbrick/ndvi.250.m.grd")


## Load some test GPS collar data (100 rows), create an $ym value (YearMonth)
load("./Data/RData/df.ua.RData")
df <- subset(df.ua, case == 1)
df <- df[101:350, ]
df$ym <- format(df$date, format = "%Y%m")
rm(df.ua)

# Convert to a spatial df:
df <- SpatialPointsDataFrame(df[ ,c("utmE", "utmN")], df,
                             proj4string = CRS("+proj=utm +zone=5 +datum=NAD83"))

# Split the spatial df by yearmonth
df.spl <- split(df, df$ym)


df$ndvi250 <- unsplit(lapply(df.spl, function(x) {
    i <- paste("ndvi_", x$ym, sep = "")
    f <- extract(ndvi.250.m[[i]], x)
    f <- f[,1]
  }), df$ym)



#-------------------------------------------------------------------------------
## SANDBOX ###




#-------------------------------------------------------------------------------
## Calculate the sum of rain values by month
## https://stackoverflow.com/questions/36722492/conditionally-extract-data-from-a-raster-stack-data-based-on-values-in-a-spatial

rain_sum <- unsplit(lapply(df.spl, function(x) {
  # current year
  y <- as.numeric(substr(x$ym, 1, 4))
  # current month
  m <- as.numeric(substr(x$ym, 5, 6))
  # if month is after Oct, start from that year's Oct
  # if month is before Oct, start from previous year's Oct
  if(m < 11) y <- y-1
  start_date <- as.Date(sprintf('%s/10/01', y))
  # if start_date is earlier than first time slice, reset to first time slice
  start_date <- max(min(as.Date(sub('rain', '01', names(rain.250.m)), '%d%Y%m')), 
                    start_date)
  end_date <- as.Date(paste0(x$date, '01'), '%Y%m%d')
  # Sequence of dates to sum over
  i <- format(seq(start_date, end_date, by='month'), 'rain%Y%m')
  # Extract values
  extract(rain.250.m[[i]], x)
}), df$ym)



















#
#
# for (i in 1:nrow(df)) {
#   subset(df, which(getZ(ndvi.250.m) == mydate))
#   if(any(format(df[[i]]$Date, format = "%Y%m") == getZ(ndvi.250.m[[i]])))  {
#     df$ndvi <- extract(ndvi.250.m, df[, c("utmE", "utmN")],
#                       method = "simple",
#                       df = T) ; beepr::beep(4)
#   }
# }


# foo <- MODIStsp_extract(ndvi.250.m, df, FUN = "mean", out_format = "dframe")


## Extract then reshape then merge
library(reshape2)
dfTemp <- extract(ndvi.250.m, df[, c("utmE", "utmN")],
                  method = "simple",
                  df = T) ; beepr::beep(4)
dfTemp <- melt(dfTemp, id = "ID")

# Create a start and end date:
dfTemp$DateStart <- paste0(dfTemp$variable, "01")
dfTemp$DateStart <- as.Date(substr(dfTemp$DateStart,  6, 16), format = "%Y%m%d")
dfTemp$DateEnd <- dfTemp$DateStart + days_in_month(dfTemp$DateStart) - 1


