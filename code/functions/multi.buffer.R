
#' multi.buffer
#'
#' @description Takes a spatial dataframe (df) and outputs a list of spatial 
#'    polygons containing buffered rings (donuts) extending from the points in the 
#'    spatial dataframe.
#'
#' @param df A spatial dataframe of point locations, projected in UTMs
#' @param d Distance buffered between individual rings. If projects UTMs, then this is meters.
#' @param nr Number of rings/donuts to output.
#'
#' @return A list of buffered spatial polygons.
#' @export
#'
#' @examples multi.buffer(dfSp, 50, 20)

multi.buffer <- function(df, d, nr) {
  
  library(raster)
  library(rgeos)
  library(rlist)
  
  # Create a list of buffered spatial polygons:
  b <- lapply(1:nr, function(x) {
    buffer(df, width=d*x, dissolve=F)
  })
  
  # Erase the inner portions of the buffered points:
  f <- lapply(2:length(b), function(i) {
    #foo <- list()
    ls <- 1:length(b[[1]])
    f1 <- lapply(ls, function(x) {  # Create the ring for a single buffered distance
      gDifference(b[[i]][x,], b[[i-1]][x,])
    })
    f2 <- do.call(bind, f1)  # Combine the indiv. sp. polys for each pt. into a sp. poly for each ring
  })
  
  f <- rlist::list.append(f, b[[1]])  # add in the central buffered circle

  return(f)
}
