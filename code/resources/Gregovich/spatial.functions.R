#buffer an spatialpointsdataframe by radius, output a spatialpolygonsdataframe with data intact
buffer.points<-function(in.points,radius){
	buffs<-lapply(1:nrow(in.points),function(x){
			a<-as(disc(radius,c(coordinates(in.points)[x,1],coordinates(in.points)[x,2])),"SpatialPolygons")
			slot(a@polygons[[1]],"ID")<-as.character(x)
			return(a)
			
		})
	buffs2<-do.call(rbind,buffs)
	out.polys<-SpatialPolygonsDataFrame(buffs2,data=slot(in.points,"data"))
	proj4string(out.polys)<-proj4string(in.points)
	return(out.polys)
}


#calc.RSF.surface takes RSF coefficients and corresponding factor names (as they are named in an input stack of covariate rasters)
#and calculates an RSF surface raster using the exponential equation.
calc.RSF.surface <- function(coefs,facts,cov.stack){
					out.RSF <- exp(sum(coefs*cov.stack[[facts]]))
				}

#convert telonics format of datetime to as.POSIXct
#x being a single or character vector of date/times as they are formatted in telonics collars.
format.telonics.time<-function(x){
		as.POSIXct(strptime(x,'%m/%d/%Y %H:%M:%S'))
	}


#gme.kern.poly produces an img, line shapefile, and poly shapefile from animal locations.
#and returns a SpatialPolygonsDataFrame of the 95% isopleth
#in.pts is the input points file
#out.* are the output files
#cell.size is the cell size of the .img output file.

gme.kern.poly<-function(in.pts,out.img,bw = 'LSCV',out.line,out.poly, q = 0.95, cell.size = 30){
        #kde portion...
        parameters <- c('in','out','bandwidth','cellsize','kernel')
        kernel <- 'GAUSSIAN'
        vals<-list(in.pts,out.img,bw,cell.size,kernel)
        R.to.GME('kde',parameters,vals,close.GME=T)
        #isopleth portion...
        parameters <- c('in','out','quantiles','poly')
        vals<-list(out.img, out.line,q, out.poly)
        R.to.GME('isopleth',parameters,vals,close.GME=T)
        #read poly back in...
        t.split <- strsplit(out.poly,'/')[[1]]
	split.dsn<-t.split[1:(length(t.split) - 1)]
	t.dsn<-paste(split.dsn,collapse='/')
	t.lay<-strsplit(t.split[length(t.split)],'.shp')[[1]][1]
        gme.poly <- readOGR(t.dsn, t.lay)
        return(gme.poly)
}



#hr.kern creates a list:
#[[1]] a raster of the output kernel density as computed by the 'kde' function of the ks package
#[[2]] a polygon of the specified percent isopleth
#[[3]] a 4X4 matrix of [1,1]= bandwidth, [1,2] and [2,1] = xy covariance, [2,2] = y variance
#loc.coords, 2 column matrix of x, y coordinates of animal locations
#cell.size is the desired res of the output kernel raster
#prob is a vector of length 1 that defines the output polygon isopleth value
hr.kern<-function(loc.coords,  h = NULL, out.res = c(30,30), probs = 0.95){
	conts <- probs*100
	grid.dim<-grid.size(loc.coords,cell.size=out.res[1])
	if(is.null(h)){
		bw.mat<-Hlscv(loc.coords)
		k<-kde(loc.coords,H=bw.mat,gridsize = grid.dim)
		h<-bw.mat
	} else {
		k<-kde(loc.coords,h=h,gridsize = grid.dim)
	}
	spkde <- image2Grid(list(x = k$eval.points[[1]], y = k$eval.points[[2]], z = k$estimate),digits=6) 
	kr<-raster(spkde)
	new.rast<-raster(extent(spkde),res = out.res)
	kr<-resample(kr,new.rast)
	threshs <- contourLevels(k, cont = conts)
	out.polys<-list()
	for(i in 1:length(threshs)){
		bin.kr <- kr
		bin.kr[bin.kr < threshs[i]]<-NA
		bin.kr[bin.kr > 0]<-1
		k.poly<-rasterToPolygons(bin.kr,dissolve=T)
		out.polys[[i]]<-k.poly
		#x11();plot(k.poly)
	}
	
	for(i.id in 1:length(out.polys)){
		slot(slot(out.polys[[i.id]],'polygons')[[1]],'ID')<-as.character(i.id)
	}

	all.polys<-do.call(rbind,out.polys)
	proj4string(all.polys)<-proj4string(kr)<-p4s.akspz1
	out<-list(kr,all.polys,h)
	return(out)
}





#grid.size computes the desired n for an n X n grid, given a location coordinates and desired cell.size.
#loc.coords is a 2-column matrix of location coordinates
#cell.size is the desired cell size
grid.size<-function(loc.coords,cell.size){
	ranges<-apply(apply(loc.coords,2,range),2,function(x){x[2] - x[1]})
	t.range<-ranges[which.max(ranges)]
	n.cells<-round(t.range/cell.size)
	return(n.cells)
}

#poly.area.sqkm calculates the area of a SpatialPolygonsDataFrame in square kilometers
#sp.poly is a SpatialPolygonsDataFrame
#watch out if there are (large) holes, as they are included in the calculation.
poly.area.sqkm<-function(sp.poly){
		sum(sapply(sp.poly@polygons[[1]]@Polygons,function(x){x@area}))/1000000
}

#clip.poly clips a polygon layer with another, either inside or outside the clipping polygon layer
clip.poly<-function(poly.to.clip,clipper,inside=T){
	if(inside){
		gI <- gIntersects(poly.to.clip, clipper, byid=TRUE)
		clip.frame.inside<-poly.to.clip@data[gI,]
		out.inside <- vector(mode="list", length=length(which(gI)))
		ii <- 1
		for (i in seq(along=gI)){
			if (gI[i]){
				out.inside[[ii]] <- gIntersection(poly.to.clip[i,], clipper)
				row.names(out.inside[[ii]]) <- row.names(poly.to.clip)[i]
				ii <- ii+1
			}
		}
		out.inside.comb <- do.call("rbind", out.inside)
		clip.frame.inside<-as.data.frame(clip.frame.inside)
		rownames(clip.frame.inside)<-sapply(out.inside.comb@polygons,function(x){x@ID})
		clipped<-SpatialPolygonsDataFrame(out.inside.comb,data=clip.frame.inside)
	} else {
		gD <- gDisjoint(poly.to.clip, clipper, byid=TRUE)
		gOver <- gOverlaps(poly.to.clip,clipper,byid=TRUE)
		clip.frame.overlap<-poly.to.clip@data[gOver,]
		out.outside <- vector(mode="list", length=length(which(gD)))
		out.overlap <- vector(mode="list", length=length(which(gOver)))
		ii <- 1
		for (i in seq(along=gD)){
			if (gD[i]){
				out.outside[[ii]] <- poly.to.clip[i,]
				row.names(out.outside[[ii]]) <- row.names(poly.to.clip)[i]
				ii <- ii+1
			}	
		}
		ii<-1
		for (i in seq(along=gOver)){
			if (gOver[i]){
				out.overlap[[ii]] <- gDifference(poly.to.clip[i,],clipper)
				row.names(out.overlap[[ii]]) <- row.names(poly.to.clip)[i]
				ii <- ii+1
			}	
		}
		outside.frame<-do.call("rbind",out.outside)
		overlap <- do.call("rbind", out.overlap)
		clip.frame.overlap<-as.data.frame(clip.frame.overlap)
		rownames(clip.frame.overlap)<-sapply(outside.frame@polygons,function(x){x@ID})
		overlap.frame<-SpatialPolygonsDataFrame(overlap,data=clip.frame.overlap)
		clipped<-spRbind(outside.frame,overlap.frame)
		return(clipped)
	}
}


#proj4strings for State plane alaska zone 1, state plane alaska zone 1 with 7 parameter transformation, and wgs84

p4s.akspz1<-'+proj=omerc +lat_0=57 +lonc=-133.6666666666667 +alpha=-36.86989764583333 +k=0.9999 +x_0=5000000 +y_0=-5000000 +datum=NAD83 +units=m +no_uoff +no_defs +ellps=GRS80 +towgs84=0,0,0'
p4s.z1w7parm<-'+proj=omerc +lat_0=57 +lonc=-133.6666666666667 +alpha=-36.86989764583333 +k=0.9999 +x_0=5000000 +y_0=-5000000 +datum=NAD83 +units=m +no_uoff +no_defs +ellps=GRS80 +towgs84=-0.995600,1.901300,0.521500,0.025915,0.009426,0.011599,0.000620'
p4s.wgs84<-'+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'


#quantile.bin bins a continuous raster into quantiles. By default, the bottom %60 of bins
#are changed to NAs to cause transparency in plotting (10 bins is the default number, x is the input raster object)
quantile.bin<- function(x,n.bins = 10,n.transparent = round(n.bins *0.6)){
			q <- quantile(x,probs = seq(0.1,0.9,by=0.1),na.rm=TRUE)
			rcl.mat <- cbind(c(0,q[1:9]),c(q[1:9],cellStats(x,'max')+0.1),1:10)
			reclass.rast<-reclassify(x,rcl=rcl.mat,include.lowest=T)
			reclass.rast[reclass.rast <= n.transparent] <- NA
			return(reclass.rast)
		}


#function to call GME, supplying tool name, parameters, and values
R.to.GME<-function(gme.function,parameters,values,script.file=file.path(getwd(),"GME.script.txt"),
				script.append=F,close.GME=F,run=T,inv=F,exe.file="C:/Program Files (x86)/SpatialEcology/GME/SEGME.exe"){
	op<-options()
	options(warn = -1)
	values<-lapply(1:length(values),function(x){if(!is.numeric(values[[x]][1])){gsub("/","\\\\",values[[x]])} else {values[[x]]}})
	sep.char<-lapply(1:length(values),function(x){
			sapply(1:length(values[[x]]),function(y){
				ifelse(is.na(as.numeric(values[[x]][y])),'"',"")
			})
		})
	statement=""
	for(i.param in 1:length(parameters)){
		if(i.param==1){
			statement<-paste(statement,gme.function,"(",sep="")
		}
		for(i.elem in 1:length(values[[i.param]])){
			t.sep<-sep.char[[i.param]][i.elem]
			starter<-paste(parameters[i.param],"= c(",sep="")
			alt.sep.mat<-matrix(c(starter,",","",")","",",","",")"),4,2,byrow=T)
			if(i.elem==1){ind<-1}else{if(i.elem==length(values[[i.param]])){ind<-4}}
			if(i.elem>1 & length(values[[i.param]])==2){ind<-2}
			if(i.elem>1 & i.elem!=length(values[[i.param]])){ind<-3}	
			alt.seps<-alt.sep.mat[ind,]
			if(length(values[[i.param]])==1){
				statement<-paste(statement,sub("c\\(","",starter),t.sep,values[[i.param]][i.elem],t.sep,sep="")
				} else {
					statement<-paste(statement,alt.seps[1],t.sep,values[[i.param]][i.elem],
							t.sep,alt.seps[2],sep="")
			}
		}		
	end.char<-ifelse(i.param==length(parameters),");",",")
	statement<-paste(statement,end.char,sep="")
	}
	write(statement,file=script.file,append=script.append)
	options(op)
	if(run){
		script.file=gsub("/","\\\\",script.file)
		argument<-paste("run(in=",shQuote(shQuote(script.file)),");",sep="")
		if(close.GME){argument<-paste("-c",argument)}
		system2(command=exe.file,args=argument,invisible=inv)
	}	
}


#reformat.garmin.datetimes
#is a function that takes the very specific format of datetimes in garmin files
#downloaded with gpsbabel, and converts them to ISO datetimes
#datetimes is a vector of characters (usually a field in a data.frame) representing datetimes like this: '27-FEB-15 1:10:07PM'
#notice the way hour ('%I') and the PM indicator ('$p') are handled.
reformat.garmin.datetimes <- function(datetimes){
						strptime(datetimes,'%d-%b-%y %I:%M:%S%p')
					}


#sampleEvery
#from Barry Rowlingson via r-sig-geo
#' sample every distance interval along a line
#'
#' @param cc a two-column matrix of x and y coords for the line
#' @param dist distance separating sample points
#'
#' @return points along the line, separated by the given distance
#measured along the line
#'
sampleEvery <- function (cc, dist){
    lengths = LineLength(cc, longlat = FALSE, sum = FALSE)
    if (any(abs(lengths) < .Machine$double.eps)) {
        wl <- which(abs(lengths) < .Machine$double.eps)
        cc <- cc[-(wl), ]
        lengths <- lengths[-(wl)]
    }
    csl = c(0, cumsum(lengths))
    maxl = csl[length(csl)]
    pts = seq(0,sum(lengths),by=dist)
    int = findInterval(pts, csl, all.inside = TRUE)
    where = (pts - csl[int])/diff(csl)[int]
    xy = cc[int, , drop = FALSE] + where * (cc[int + 1, , drop = FALSE] -
        cc[int, , drop = FALSE])
    if (nrow(xy) < 1)
        return(NULL)
    return(xy)
}


#tpiw is a function to perform a topographic position index analysis
#at varying focal windows, where x is the elevation raster, and w the length (in pixels) of a side of the focal window.
tpiw <- function(x, w) {
	m <- matrix(1/(w^2-1), nc=w, nr=w)
	m[ceiling(0.5 * length(m))] <- 0
	f <- focal(x, m)
	x - f
}

#scale a raster or raster stack.
scale.rast <- function(x){
			(x - cellStats(x, 'mean'))/cellStats(x,'sd')
		}

#ArcGIS sometimes has problems constructing .prj files from proj4string information.
#this script writes a shapefile from an R spatial feature object, and then ensures that the 
#AK State plane zone 1 projection is specified as being of natural origin, which is not the case if 'writeOGR' is used without this wrapper.
write.5001<-function(obj,dsn,out.name){
       	writeOGR(obj,dsn,out.name,overwrite_layer=T,driver="ESRI Shapefile")
        prj.path<-paste(file.path(dsn,out.name),".prj",sep="")
        prj<-scan(prj.path,what=character())
       	prj<-sub("Hotine_Oblique_Mercator_Azimuth_Center","NAD_1983_StatePlane_Alaska_1_FIPS_5001",prj)
	prj<-sub("Hotine_Oblique_Mercator_Azimuth_Center","Hotine_Oblique_Mercator_Azimuth_Natural_Origin",prj)
        write(prj,prj.path)
}


