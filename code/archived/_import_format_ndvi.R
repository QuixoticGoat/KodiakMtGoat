
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
  ndvi_path <- "./resources/geodata/MODIS/ndvi"

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
  writeRaster(ndvi.250, "./data/derived/geodata/ndvi.250.grd",
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
writeRaster(ndvi.500, "./data/derived/geodata/ndvi.500.grd",
            format = "raster",
            options = c("INTERLEAVE=BAND"),
            progress = "text",
            prj = T,
            overwrite = T)

# Run it:
ndvi.1000  <- multiFocal(ndvi.250, 1000)

# Save it:
writeRaster(ndvi.1000, "./data/derived/geodata/ndvi.1000.grd",
            format = "raster",
            options = c("INTERLEAVE=BAND"),
            progress = "text",
            prj = T,
            overwrite = T)



#-------------------------------------------------------------------------------
## Get monthly composites of NDVI values (raw data are 16 day intervals)

## Load the ndvi rasterstacks at 250, 500, and 1000 m buffers:
ndvi.250 <- stack("./data/derived/geodata/ndvi.250.grd")
ndvi.500 <- stack("./data/derived/geodata/ndvi.500.grd")
ndvi.1000 <- stack("./data/derived/geodata/ndvi.1000.grd")

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
  filename = paste("./data/derived/geodata/", deparse(substitute(r)), ".m.grd", sep="")
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

MonthlyNDVI(ndvi.250)
MonthlyNDVI(ndvi.500)
MonthlyNDVI(ndvi.1000)


#-------------------------------------------------------------------------------
## Extract data from each matching NDVI to each GPS collar fix

## Load the ndvi:
ndvi.250.m <- stack("./data/derived/geodata/ndvi.250.m.grd")
ndvi.500.m <- stack("./data/derived/geodata/ndvi.500.m.grd")
ndvi.1000.m <- stack("./data/derived/geodata/ndvi.1000.m.grd")


loadit <- function() {
  ## Loads GPS collar data, creates an $ym value (YearMonth), cleans it, and
  ## 
  load("./data/derived/df.ua.RData")
  df <- df.ua
  df$ym <- format(df$date, format = "%Y%m")
  df <- subset(df, date > "2013-07-01 00:00:00") # remove June '13 fixes
  # df <- df[1:100,]  # Take just the first 100 rows...
  rm(df.ua)
  
  df <- df[, c("utmE", "utmN", "ym")]
  # Convert to a spatial df:
  df <- SpatialPointsDataFrame(df[ ,c("utmE", "utmN")], df,
                               proj4string = CRS("+proj=utm +zone=5 +datum=NAD83"))
  df.spl <- split(df, df$ym)  # Split the spatial df by yearmonth
  rm(df)
  message("Loaded")
  dat <- list("df" = df, "df.spl" = df.spl)
  return(dat)
}

## NDVI.250
start.time <- Sys.time() # Start time
dat <- loadit()  # Load the list of spatial dataframes
# Extract the values to each fix (TAKES A WHILE..):
ndvi.250 <- data.frame("ndvi.250" = unsplit(lapply(dat$df.spl, function(x) {
  i <- paste("ndvi_", x$ym, sep = "")
  f <- raster::extract(ndvi.250.m[[i]], x)
}), dat$df$ym))
end.time <- Sys.time()  # End time
end.time-start.time  # Duration
save(ndvi.250, file = "./data/derived/df.ndvi.250.RData")

## NDVI.500
start.time <- Sys.time() # Start time
df.spl <- loadit()  # Load the list of spatial dataframes
ndvi.500 <- data.frame("ndvi.500" = unsplit(lapply(df.spl, function(x) {
  i <- paste("ndvi_", x$ym, sep = "")
  f <- extract(ndvi.500.m[[i]], x)
  f <- f[,1]
}), df$ym))
end.time <- Sys.time()  # End time
end.time-start.time  # Duration

save(ndvi.500, file = "./data/derived/df.ndvi.500.RData")

## NDVI.1000
start.time <- Sys.time() # Start time
df.spl <- loadit()  # Load the list of spatial dataframes
ndvi.1000 <- unsplit(lapply(df.spl, function(x) {
  i <- paste("ndvi_", x$ym, sep = "")
  f <- extract(ndvi.1000.m[[i]], x)
  f <- f[,1]
}), df$ym)
save(ndvi.1000, file = "./data/derived/df.ndvi.1000.RData")
end.time <- Sys.time()  # End time
end.time-start.time  # Duration

load("./data/derived/df.ndvi.250.RData")
load("./data/derived/df.ndvi.500.RData")
load("./data/derived/df.ndvi.1000.RData")
ndvi <- cbind(ndvi.250, ndvi.500, ndvi.1000)  # Combine them

