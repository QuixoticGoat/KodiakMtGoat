#### HOME RANGES  ######
#### McCrea Cobb  ######
#### 12/30/2015   ######



## Load required packages
library(adehabitatHR)

## Format Data (convert from data frame to SpatialPointsDataFrame) ----
dfSp <- SpatialPointsDataFrame(df[ ,6:5], df, 
                                proj4string = CRS("+proj=longlat +datum=WGS84"))

# Reproject into UTM NAD83, Zone 5N
dfSp <- spTransform(dfSp, CRS("+proj=utm +zone=5 +datum=NAD83"))
head(dfSp) # Check the headers to make sure it worked
plot(dfSp, col = dfSp$CollarID) # plot the points by CollarID




## MINIMUM CONVEX POLYGONS (MCPs) ----

require(adehabitatHR)

# Create 95% MCPs:
MCPs <- mcp(dfSp[ ,11], percent = 95, unout = "km2")
# Calculate the areas within the MCPs
MCPArea <- mcp.area(dfSp[, 11], percent = seq(50, 100, by = 5),
                    unout = "km2", plotit = FALSE)
MCPsArea # Look at it

# Convert the tbl of 95% MCPs areas to a dataframe
MCPdf <- data.frame(MCPs)  
row.names(MCPdf) <- NULL  # Remove row names
colnames(MCPdf) <- c("CollarID", "Area")  # Rename the columns

mean(MCPdf$area)  # Calculate mean 95% MCP area (km2)
sd(MCPdf$area)/sqrt(length(MCPdf$area))  # Calculate SE of all 95% MCP areas

#  Create a map of 95% MCPs for all collars:
library(ggmap)
library(rgdal)
# transform to WGS84
MCPs <- spTransform(MCPs, CRS("+proj=longlat +datum=WGS84"))
#  Get the map data (requires network connection):
kodiak <- get_map("kodiak island", zoom = 8, maptype = "terrain")
MapMcps <- ggmap(kodiak, extent = "panel") + 
  geom_polygon(data = MCPs, aes(x = long, y = lat, group = id, fill = id), 
               alpha =  0.5) +
  geom_path(data = MCPs, aes(x = long, y = lat, group = id, fill = id)) +
  geom_point(aes(x = Long, y = Lat), data = df, size = 1, pch = 20)
MapMcps
ggsave("Plots/Maps/HRs/MCPs/MapMcps.jpg", MapMcps, width = 10, height = 10)




## FIXED KERNEL HRs ----

require(adehabitatHR)

# Create raster kernels:
Kernel <- kernelUD(dfSp, h = "href")
image(Kernel)  # Plot them
# Create contour lines for first animal, plot result:
KernelContour <- as.image.SpatialGridDataFrame(Kernel[[1]])
contour(KernelContour, add = TRUE)  

# Plot 95% fixed kernel polygons
KernelVert <- getverticeshr(Kernel, percent = 99, unin = "m", unout = "km2")
plot(KernelVert)
crs(dfSp) <- crs(StudyArea)

# Crop the kernel to Kodiak Island
KernelVert <- crop(KernelVert, StudyArea)

# get means and SDs (first transform to WGS84)
KernelVert <- spTransform(KernelVert, CRS("+proj=longlat +datum=WGS84"))
KernelVertdf <- data.frame(KernelVert)  # convert to df
mean(KernelVertdf$area)  # Calculate mean 95% kernel hr area (km2) 
sd(KernelVert$area)/sqrt(length(KernelVert$area))  # Calculate SE

# Map kernel polygons
library(ggmap)
library(rgdal)

KernelVert <- spTransform(KernelVert, CRS("+proj=longlat +datum=WGS84"))

#  Get the map data (requires network connection):
kodiak <- get_map("kodiak island", zoom = 8, maptype = "terrain")
MapKernel <- ggmap(kodiak, extent = "panel") + 
  geom_polygon(data = KernelVert, 
               aes(x = long, y = lat, group = id, fill = id), 
               alpha =  0.5) +
  geom_path(data = KernelVert, 
            aes(x = long, y = lat, group = group, fill = id)) +
  geom_point(aes(x = Long, y = Lat), data = df, size = 1, pch = 20)
MapKernel
ggsave("Plots/Maps/HRs/Kernels/MapKernel.jpg", MapKernel, width = 8, height = 8)  # save it




## KERNEL BROWNIAN BRIDGE HOME RANGES ----

require(adehabitatHR)

# First, create a list of trajectory objects for each animal:

traj <- function(d) {
  ## Function to create a list of trajectory objects for each animal. 
 
   require(adehabitatHR)
 
   dat = split(d, d$CollarID)  # Split df by CollarId into a list
  
   # Create a list of trajectory objects
  datTraj = lapply(dat, function(x) as.ltraj(xy = x[ ,c("Long", "Lat")], 
                    date = x$Date,
                    id = "CollarID"))
  
  refda <- strptime("21:00", "%H")  # reference date for SetNA function below
  # Set to interval between consequetive fixes at 1 fix/12.5 hrs, and add
  # NAs into times with no fix success.
  datTrajReg = lapply(datTraj, function(x) setNA(x, x[[1]]$date[1], 12.5, units = "hour"))
  
  return(datTrajReg)
}


## Run the traj() function and store list of trajectories:
dfTrajLs <- traj(df)

# Convert trajectories to a list of dataframes:
dfTraj = lapply(dfTrajLs, function(x) ld(x))
  
# Estimate sig1 for BB kernels:
Sig1Est <- lapply(dfTrajLs, function(x) liker(x, sig2 = 0.00001,  # To estimate sig1
                                        rangesig1 = c(-1, 1), 
                                        plotit = FALSE,
                                        byburst = FALSE) )
# shows that sig1 = 0.00501

# Create a list of BB kernel HRs:
KernelBB <- lapply(dfTrajLs, function(x) kernelbb(x, sig1 = 0, sig2 = .001, grid = 100, 
                     extent = .4))

# Look at a plot of the BB kernel and the 95% contour:
image(KernelBB[[1]])
plot(getverticeshr(KernelBB[[1]], 95), 
     add = TRUE,
     lwd = 2)

# Create contours from the BB kernel and add it to the previous image:
KernelBBContour <- as.image.SpatialGridDataFrame(KernelBB[[1]])
contour(KernelBBContour, add = TRUE) 

# Create the 95% vertice contour of the BB kernel and add to the plot as well:
KernelBBVert <- getverticeshr(KernelBB[[1]], percent = 95, unin = "m", unout = "km2")
plot(KernelBBVert, add = TRUE)




## Display and save pretty maps of the BB kernels

# Required package:
library(ggmap)

# Project the BB kernel vertice contour to WGS84:
proj4string(KernelBBContour) <- "+proj=longlat +datum=WGS84"

#  Get the map data (requires network connection):
kodiak <- get_googlemap(center = c(lon = mean(KernelBBVert@polygons[[1]]@Polygons[[1]]@coords[,1]), 
                                   lat = mean(KernelBBVert@polygons[[1]]@Polygons[[1]]@coords[,2])),
                  zoom = 11, 
                  maptype = "terrain")

# Create the map object:
MapKernelBB <- ggmap(kodiak, extent = "panel") + 
  geom_polygon(data = KernelBBVert, 
               aes(x = long, y = lat, group = id, fill = id), 
               alpha =  0.5) +
  geom_path(data = KernelBBVert, 
            aes(x = long, y = lat, group = id, fill = id))
 geom_point(aes(x = Long, y = Lat), data = df, size = 1, pch = 20)
MapKernelBB
ggsave("MapKernel.jpg", MapKernelBB, width = 8, height = 8)  # save it

