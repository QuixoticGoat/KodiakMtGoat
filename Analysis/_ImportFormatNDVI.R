##################################################
#### Import format NDVI and snow index time series data
##################################################


require(raster)
require(rgdal)



#################
## Import the geotiffs and save as a raster stack
## Data were acquired on 3/23/18 using USGS AppEEARS.
#################

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
  ndvi.250 <- projectRaster(ndvi.250, res = c(250, 250), crs = "+proj=utm +zone=5 +datum=NAD83")

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




########
### Extract the focal mean pixel values for each covariate by 250, 500 and 1000 m.
########

 # 1. Create raster stacks for each focal mean buffer distance (250, 500, 1000):

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

# Run it:
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



#######
### Get monthly composites
#######

#get the date from the names of the layers and extract the YearMonth
indices <- format(as.Date(substr(names(ndvi.250), 35, 41), format = "%Y%j"), format = "%Y%m")
indices <- as.numeric(indices)

#sum layers
ndvi.250.m <- stackApply(ndvi.250, indices, fun = mean)

indices <- unique(indices)
names(ndvi.250.m) <- paste0("ndvi_", indices)


ndvi.250.m <- stack(ndvi.250.m)

indices <- paste0(as.character(indices), "15")
indices <- as.Date(indices, format = "%Y%m%d")

ndvi.250.m <- setZ(ndvi.250.m, indices)

# Save it:
writeRaster(ndvi.250.m, "./Data/GIS/MODIS/ndvi/Rasterbrick/ndvi.250.m.grd",
            format = "raster",
            options = c("INTERLEAVE=BAND"),
            progress = "text",
            prj = T,
            overwrite = T)


##################################
### Extract data from NDVI to fixes
#################################

library(MODIStsp)


## Load the ndvi:
ndvi.250.m <- stack("./Data/GIS/MODIS/ndvi/Rasterbrick/ndvi.250.m.grd")


## Load some data:
load("./Data/RData/df.ua.RData")
df <- subset(df.ua, case == 1)
df <- df[1:100, ]
df$ym <- format(df$date, format = "%Y%m")
rm(df.ua)

# Create a spatial df of the df.ua dataframe:
df <- SpatialPointsDataFrame(df[ ,c("utmE", "utmN")], df,
                             proj4string = CRS("+proj=utm +zone=5 +datum=NAD83"))

# Split the spatial df by yearmonth
df.spl <- split(df, df$ym)

##


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
  start_date <- max(min(as.Date(sub('rain', '01', names(ndvi.250.m)), '%d%Y%m')), start_date)
  end_date <- as.Date(paste0(x$date, '01'), '%Y%m%d')
  # Sequence of dates to sum over
  i <- format(seq(start_date, end_date, by='month'), 'rain%Y%m')
  # Extract values
  extract(ndvi.250.m[[i]], x)
}), df$ym)



foo <- for(i in 1:length(df.spl)) {
  for(i in 1:nlayers(ndvi.250.m)) {
    if(as.numeric(substr(names(ndvi.250.m[[1]]), 6, 12)) == as.numeric(names(df.spl[i]))) {
      extract(ndvi.250.m[[i]], df.spl[[i]])
    }
  }
}



substr(names(ndvi.250.m[[1]]), 6, 12)












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


