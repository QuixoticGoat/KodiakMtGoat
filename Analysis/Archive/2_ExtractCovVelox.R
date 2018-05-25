
#########################################################################
#### Extract covariate information to a used/available data frame
#### Author: McCrea Cobb
#### Last modified: 9/2017
#########################################################################


## ---- ExtractCovVelox.R


### STEP 1: Define "available" points at the kernel home range scale


avail.fn <- function(df = dfSp, StudyArea = StudyArea) {

  require(sp)
  require(snow)
  require(rgeos)
  require(adehabitatHR)
  require(rasterVis)  # for crs()

  ## Create a 99% kernel HR from all "used" fixes
  Kernel <- kernelUD(dfSp, h = "href")
  KernelVert <- getverticeshr(Kernel, percent = 99, unin = "m", unout = "km2")
  #crs(dfSp) <- crs(StudyArea)
  # Crop the kernel to Kodiak Island (i.e., no ocean)
  KernelVert <- crop(KernelVert, StudyArea)
  # plot(KernelVert)

  ## Generate a random sample of "available" points within the kernel:
  AvailS <- spsample(KernelVert, type = "random", n = length(df$CollarID)*10)

  # Take a look at a plot of available points:
  # plot(StudyArea)
  # plot(KernelVert, add = T)
  # plot(AvailS, add = T)

  # Reproject available points to WGS84:
  AvailS <- spTransform(AvailS, CRS = CRS("+proj=longlat +datum=WGS84"))

  # Convert available points spatial data frame to a data frame:
  Avail <- data.frame(Response = 0, Lat = AvailS@coords[,2],
                      Long = AvailS@coords[,1])
  Avail$CollarID <- rep(df$CollarID, 10)  # Add CollarID

  ## Bind "used" and "available" dfs:
  common_cols <- intersect(colnames(df), colnames(Avail))
  dfRSF <- rbind(subset(df[common_cols]),
                 subset(Avail[common_cols]))

  return(dfRSF)
}

# Run it:
dfRSF <- avail.fn()
save(dfRSF, file="./Data/Rdata/dfRSF.Rdata", compress="gzip")


### STEP 2: Extract the covariate data to used and available points:

extract.fn <- function(df = dfRSF, w) {

  require(velox)   # velox
  require(rgeos)   # rgeos()

  # Convert dfRSF to a SpatialPoint called dfRSFSp:
  dfRSFSp <- SpatialPoints(df[ ,2:1],
                           proj4string = CRS("+proj=longlat +datum=WGS84"))

  # Reproject dfRSFSp into UTM NAD83, Zone 5N:
  dfRSFSp <- spTransform(dfRSFSp, CRS("+proj=utm +zone=5 +datum=NAD83"))

  ## Create buffered spatial polygons from dfRSFSp:
  dfSpBuff <- gBuffer(dfRSFSp, byid = T, width = w)

  # Load the habitat covariate rasters:
  Covar <- stack("./Data/GIS/RasterStack/Covar.grd")

  # Create the vx object:
  covar.vx <- velox(Covar)

  # Extract data to the spatial polygons:
  covar.ex <- covar.vx$extract(sp = dfSpBuff, fun = mean)

  # Convert to a df and name columns:
  covar.ex <- as.data.frame(covar.ex)
  colnames(covar.ex) <- paste0(names(Covar), ".", w)

  # Combine with the original df:
  dfRSF.new <- cbind(dfRSF, covar.ex)

  removeTmpFiles(h = 0)

  return(dfRSF.new)
}

# Run it for 30, 100, 500, and 1000 m buffers:
dfRSF30 <- extract.fn( ,30)   # 30 m buffer
dfRSF100 <- extract.fn( ,100)      # 100 m buffer
dfRSF500 <- extract.fn( ,500)      # 500 m buffer
dfRSF1000 <- extract.fn( ,1000)    # 1000 m buffer

# Save them:
save(dfRSF30, file = "./Data/Rdata/dfRSF30.Rdata", compress = "gzip")
save(dfRSF100, file = "./Data/Rdata/dfRSF100.Rdata", compress = "gzip")
save(dfRSF500, file = "./Data/Rdata/dfRSF500.Rdata", compress = "gzip")
save(dfRSF1000, file = "./Data/Rdata/dfRSF1000.Rdata", compress = "gzip")


### STEP 3. Combine them into a single df

dfRSF <- cbind(dfRSF30, dfRSF100[5:17], dfRSF500[5:17], dfRSF1000[5:17])

# Save the combined df:
save(dfRSF, file = "./Data/Rdata/dfRSF.Rdata", compress = "gzip")

# clean up:
rm(dfRSF30, dfRSF100, dfRSF500, dfRSF1000)


