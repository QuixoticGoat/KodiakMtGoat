
################################################################################
# mapAll() is a function to produce a map of all individual (in this case      #
# mountain goats on Kodiak) in a dataset. Requires internet connectivity to    #
# download basemaps from Google.  Mountain Goat GPS Collar Maps, Kodiak        #
# Island, Alaska                                                               #
# mapCollars() is a function that creates seperate maps (one/page) for each    #
# individual collared animal in a dataset. Adds a title with with the collar   #
# ID and the date of last fix                                                  #
#                                                                              #
# Author: McCrea Cobb <mccrea_cobb@fws.gov>                                    #
# Date created: 8/3/2017                                                       #
# Date last edited: 4/17/2018                                                  #
################################################################################


# ------------------------------------------------------------------------------
## ---- MapAll

mapAll <- function(dat) {

  # Required packages:
  library(ggmap)
  library(ggsn)
  library(rgdal)

  basemap = get_googlemap("Kodiak Island, Alaska",
                          zoom = 8,
                          size = c(600,640),
                          scale = 2,
                          maptype = "terrain")
  bb = attr(basemap, "bb")
  # Create basic map object:
  Map = ggmap(basemap, extent = "panel", legend = "topleft") +
    geom_point(aes(x = Long, y = Lat),
               dat = dat,
               size = .5,
               pch = 20) +
    scalebar(data = NULL,
             location = "bottomright",
             dist = 25,
             height = 0.01,
             dd2km = TRUE,
             model = "WGS84",
             st.bottom = FALSE,
             st.size = 3,
             y.min = (bb[1,1]) + 0.1,
             x.min = bb[1,2],
             y.max = bb[1,3],
             x.max = (bb[1,4]) - 0.15)
  Map
}

# Run and save it:
mapAll(df)
ggsave(filename = paste("Plots/Maps/_MapAll", ".pdf", sep = ""))


# -----------------------------------------------------------------------------
## ---- MapCollars

mapCollars <- function(dat) {
  # Function to create individual maps for each collared goat and
  # add a title with the collar ID and the date of last fix.

  require(ggmap)
  library(ggsn)
  library(rgdal)

  dat = split(dat, dat$CollarID)
  CollarID = names(dat)

  Basemap <-  lapply(dat, function(x) {
    Sys.sleep(10)  # add a delay to make google happy..
    get_googlemap(center = c(lon = mean(x$Long), lat = mean(x$Lat)),
                  zoom = 10,
                  size = c(600,640),
                  maptype = "terrain",
                  messaging = FALSE)
  })

  Maps = lapply(CollarID, function(i) {

    bb = attr(Basemap[[i]], "bb")

    dfLastFix = data.frame(dat[[i]][dat[[i]]$Date==max(dat[[i]]$Date), ])

    ggmap(Basemap[[i]], extent = "panel", legend = "topleft") +
      geom_path(aes(x = Long, y = Lat),  # Path
                data = dat[[i]],
                size = 0.5) +
      geom_point(aes(x = Long, y = Lat), # All fixes
                 data = dat[[i]],
                 size = 1) +
      geom_point(aes(x = Long, y = Lat),  # Last fix
                 data = dfLastFix,
                 size = 2,
                 color = "red") +
      scalebar(data = NULL,
               location = "bottomright",
               dist = 5,
               height = 0.01,
               dd2km = TRUE,
               model = "WGS84",
               st.bottom = FALSE,
               st.size = 3,
               y.min = (bb[1,1]) + 0.01,
               x.min = bb[1,2],
               y.max = bb[1,3],
               x.max = (bb[1,4]) - 0.05) +
      ggtitle(paste(i,
                    " ",
                    "(",
                    "Last fix: ",
                    as.Date(dfLastFix$Date),
                    ")",
                    sep = ""))
  })

  names(Maps) <- CollarID

  return(Maps)
}

#Run function and save maps as pdfs:
Maps <- mapCollars(df)

lapply(names(Maps), function(i) {
  ggsave(filename = paste("./Output/Maps/Collars/", i, "Map.pdf", sep = ""),
         plot = Maps[[i]])
})
