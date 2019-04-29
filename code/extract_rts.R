### Using rst package to extract ndvi values to gps fixes


library(rts)
library(tidyverse)
library(sf)
library(raster)

# Load the df:
load("./Data/RData/df.ua.RData")
df <- df.ua
df$ym <- format(df$date, format = "%Y%m")
df <- subset(df, date > "2013-07-01 00:00:00") # remove June '13 fixes
df$time <- as.yearmon(df$date)

df <- df[, c("utmE", "utmN")]
# df <- df[1:1000,]  # Take just the first 100 rows...
rm(df.ua)

# df <- st_as_sf(x = df, 
#                 coords = c("utmE", "utmN"), 
#                 crs = "+proj=utm +zone=5 +datum=NAD83")

df <- SpatialPoints(df[ ,c("utmE", "utmN")],
                             proj4string = CRS("+proj=utm +zone=5 +datum=NAD83"))

# time <- substr(names(ndvi.250.m), 6, 13) 
# time <- paste(substr(time, 1, 4), "-", substr(time, 5, 7), sep = "")
# time <- as.yearmon(time)
# 
# ndvi <- rts(ndvi.250.m, time)
# 
# timels <- df$time
# df$foo <- NA
# 
# for (i in 1:length(timels))
#   df[foo, i]

start.time <- Sys.time()
foo <- raster::extract(ndvi.1000.m, df)
end.time <- Sys.time()
end.time - start.time



foo <- df %>% 
  group_by(time) %>%
  i <- df$time %>%
  extract(ndvi.250.m[[i]], df)


library(raster)

# Vector of dates
dates <- format(seq(as.Date('1999/1/11'), as.Date('2016/1/10'), by='month'), '%Y%m')

# RasterStack with random data
s <- setNames(stack(replicate(length(dates), raster(matrix(runif(10000), 100)))), 
              paste0('rain', dates))

# Create a SpatialPointsDataFrame with some random dates and coords
d <- data.frame(x=runif(10000), y=runif(10000), date=sample(dates, 10000, replace=T))
coordinates(d) <- ~x+y

# Split the spdf by date
d_by_date <- split(d, d$date)


foo <- data.frame(unsplit(lapply(d_by_date, function(x) {
  # i <- x$date
  i <- paste("rain", x$date, sep = "")
  
  f <- raster::extract(s[[i]], x)
}), d$date))



raster::raster(100, 100, rnorm(1, 100))
