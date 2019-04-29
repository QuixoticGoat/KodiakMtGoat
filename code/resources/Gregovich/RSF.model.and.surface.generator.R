rm(list=ls())
#load required packages, note that you will need to install the twostageRSF package from a local zip file.
#within RStudio the way to install a local, zip package is to go to the “Tools” menu, select “Install Packages” 
#and then select “Package Archive File (.zip; .tar.gz)” from the “Install from:” dropdown box.
library(twostageRSF)
#and then the rest of the required packages if you haven't already installed them...
packagesRequired <- c("raster", "rgdal", "plotKML", "dismo")
packagesNeeded <- packagesRequired[!(packagesRequired %in% installed.packages()[,"Package"])]
if(length(packagesNeeded)) install.packages(packagesNeeded)
remove(packagesNeeded, packagesRequired)
library(raster);library(rgdal);library(plotKML);library(dismo)
#_________INPUT: enter root folder, where you extracted the 'RSF.workshop' folder, here_____________
#eg. root.folder<- 'C:/Joesdata'
#_________End INPUT_________________________________________________________________________________
root.folder<- 'C:/dave/temp'
main.folder <- file.path(root.folder,'RSF.workshop')
#create a subfolder that points to all geographic information regarding the target locale
#in this case, the Takhin Ridge near Haines, Alaska ('takhin')
locale <- 'takhin'
locale.folder<-file.path(main.folder,locale)
#obtain some important raster processing functions supplied with the workshop materials...
source(file.path(main.folder,'code','spatial.functions.R'))

#load all modelling inputs into your R session....
load(file.path(locale.folder,'available/available.Rdata'))
load(file.path(locale.folder,'locations/rdata/locations.Rdata'))
load(file.path(locale.folder,'means.and.sds/means.and.sds.Rdata'))
load(file.path(locale.folder,'covariates/covariate.names.Rdata'))
load(file.path(locale.folder,'background.image/bg.image.Rdata'))
covariate.stack <- stack(file.path(locale.folder,'covariates','covariate.stack.tif'))
names(covariate.stack)<-covariate.names

#___________MORE INPUTS__________________________________________________________________________________
#which field identifies animals in your input used and available files?
ID.field<-'animal_ID'

#what is the name of this model? (output files and folders will be named via this character string)
mod.name <- 'test'


#what terms do you want to put into your model?
#full complement of terms
model.terms <- c('dem_s','dem_sq','slope_s','slope_sq','solrad_s','solrad_sq','edist_s','edist_sq','conifer','otherveg','nonveg')
#a subset...
#model.terms <- c('dem_s','dem_sq','slope_s','edist_s','edist_sq','conifer','otherveg','nonveg')

#build output location from the name of your model
out.fold<-file.path(locale.folder,'output')
kml.fold<-file.path(out.fold,mod.name,'kml.rasters');if(!file.exists(kml.fold)){dir.create(kml.fold)}
tif.fold<-file.path(out.fold,mod.name,'tif.rasters');if(!file.exists(tif.fold)){dir.create(tif.fold)}


#run the RSF model!!!
test<-RSF(used=locations,available=available,mod.terms = model.terms,model.name=mod.name,
		output.folder = out.fold,ID.field='animal_ID',means.and.sds=means.and.sds,out.plot.type='jpeg')
#take a look at the components of the output from the model 'test'
test

#build output raster folders from the name of your model
kml.fold<-file.path(out.fold,mod.name,'kml.rasters');if(!file.exists(kml.fold)){dir.create(kml.fold)}
tif.fold<-file.path(out.fold,mod.name,'tif.rasters');if(!file.exists(tif.fold)){dir.create(tif.fold)}

#obtain the model coefficients and their names
model.coefs <- test[[1]][,1]; model.factors <- names(model.coefs)

#calculate the output RSF surface
out.RSF <- calc.RSF.surface(model.coefs,model.factors,covariate.stack)

#reclassify the raster into quantile bins
reclass.rast <- quantile.bin(out.RSF,n.bins = 10,n.transparent = 6)

#supply transformation parameters, and project output raster to WGS84
proj4string(reclass.rast)<-p4s.z1w7parm
reclass.wgs<- projectRaster(reclass.rast,crs='+init=epsg:4326',method='ngb')

#define colors of output raster
t.colors <- colorRampPalette(colors = c('yellow','orange','red'))(10)
t.colors <- paste(t.colors,c('1A','33','4D','73','8C','99','B3','CC','E6','FF'),sep='')

#make a google earth file
KML(reclass.wgs,filename=file.path(kml.fold,'bin_rast.kml'),blur = 3,
		col=t.colors,overwrite=T)

#plot the raster on top of a 2-D google earth image.
x11(10,10)
plotRGB(bg.image)
plot(reclass.wgs,add=T,legend=F,col=t.colors)

#identify raster values.
#click(reclass.wgs)
























