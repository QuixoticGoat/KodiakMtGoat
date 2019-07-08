### Function to create maps of last fixes of Kodiak mountain goats
### McCrea Cobb
### 10/4/16

## First run the "ImportFormat" script on the raw .txt file from ATS ###


MortMaps <- function(id) {
  
  # Load the required packages:
  library(ggmap)
  library(ggsn)
  library(rgdal)

  # Subset the dataframe:
  dfLastFix = subset(df, CollarID == id)
  dfLastFix = subset(dfLastFix, Date == max(Date))
  dfLastFix = droplevels(dfLastFix)
  
  # Load a basemap from google:
  basemap = get_googlemap(center = c(lon = dfLastFix$Long,
                                     lat = dfLastFix$Lat), 
                          zoom = 13, 
                          size = c(600,640),
                          scale = 2,
                          maptype = "hybrid")
  
  bb = attr(basemap, "bb")  # Create bounding box for scalebar
  
  # Create map object:
  Map = ggmap(basemap, extent = "panel", legend = "topleft") + 
    geom_point(aes(x = Long, y = Lat), 
               data = dfLastFix, 
               size = 6,
               col = "yellow",
               pch = 13) +
    scalebar(data = NULL,
             location = "bottomright",
             dist = 2,
             height = 0.01,
             dd2km = TRUE,
             model = "WGS84",
             st.bottom = FALSE,
             st.size = 3,
             y.min = (bb[1,1]) + 0.001,
             x.min = bb[1,2],
             y.max = bb[1,3],
             x.max = (bb[1,4]) - 0.005) +
    ggtitle(paste("IG27",
                  "\n",
                  "Last fix: ",
                  as.Date(dfLastFix$Date),
                  "\n",
                  dfLastFix$Lat,
                  " N     ",
                  dfLastFix$Long,
                  " W",
                  sep = ""))
  
  # Save it:
  ggsave(filename = paste("./Plots/Maps/", dfLastFix$CollarID, ".pdf"),
       plot = junkmap)

  Map  # Show the map
}

# Run the function:
MortMaps(id = "IG34")

