#' buff.extract
#'
#' @description Takes a spatial dataframe of GPS collar fixes and returns a 
#' dataframe of extracted habitat values (from the input raster stack) at
#' buffered rings surrounding each GPS collar fix. Output can be used to
#' produce an RSF. 
#'
#' @param dfsp A spatial dataframe of GPS collar fixes, projected in UTM NAD83.
#' @param covar.vx A velox raster layer of habitat covariates. See the velox package for details on how to create this layer from a raster stack.
#' @param d Distance (m) buffered between individual rings.
#' @param nr Number of rings/donuts to output.
#' @param n A character string of names of the habitat variables in the raster
#'
#' @return A dataframe of mean habitat values (from the input raster stack) extracted from buffered rings surrounding each GPS collar fix.
#' @export
#'
#' @examples buff.extract()


buff.extract <- function(dfSp, covar.vx, d, nr, n) {
  
  start_time <- Sys.time()  # For benchmarking code
  
  source("./code/functions/multi.buffer.R")  # source multi.buffer()
  
  dfSp.split <- split(dfSp, dfSp$id)  # split gps collar data into lists by collar
  
  # 1:number of collars
  out <- lapply(31:40, function(x) {  # loop across 1:n collars (here's where you can subset/split the collar dataset for increased speed)
    
    b <- multi.buffer(dfSp.split[[x]], d, nr)  # run multi.buffer() a collar, then:
    
    b.ext <- lapply(1:length(b), function(y) {  # for each buffered distance,
      temp <- as.data.frame(covar.vx$extract(sp = b[[y]], fun = mean))  # extract habitat data from each ring
      colnames(temp) <- paste0(n, "_", y*d)  # name each column with habitat info using habitat name from covar raster and the distance info
      return(temp)
    })
    b.ext <- as.data.frame(b.ext)  # merge the lists for each distance into a single dataframe for each collar
    return(b.ext)
  })
  
  names(out) <- names(dfSp.split[1:2])  # name each collar
  
  # convert list to a dataframe
  out1 <- do.call(rbind.data.frame, out)
  out1$id <- as.factor(rep(names(out), sapply(out, nrow)))
  rownames(out1) <- NULL
  
  end_time <- Sys.time()
  total_time <- round((end_time - start_time)/60, 3)
  message(paste("It took", total_time, "minutes to run."))
  
  return(out1)
}
