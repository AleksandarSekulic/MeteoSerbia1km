library(devtools)
devtools::install_github("AleksandarSekulic/Rmeteo")
# install.packages("meteo", repos="http://R-Forge.R-project.org")
library(meteo)

library(sp)
library(spacetime)
library(raster)
library(rgdal)
library(CAST)
library(doParallel)
library(ranger)
library(gstat)
library(plyr)
library(ggsn)
library(ggpubr)
library(ggplot2)
library(RSAGA)
library(fields)
library(caret)
library(grid)
library(gridExtra)
library(lemon)
library(GISTools)
library(maps)  
#library(devtools)
#devtools::install_github("eliocamp/ggnewscale")
library(ggnewscale)

# repeat code for each variable!
var = "tmax"
# var = "tmin"
# var = "tmean"
# var = "slp"
# var = "prcp"

wd=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

borders <- readOGR("borders/osm_nominatim/Polygon.shp")
# projection for Serbia UTM zone 34
wgs84 <- CRS("+proj=longlat +datum=WGS84")
utm34 <- CRS('+proj=utm +zone=34 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')

border.buffer <- buffer(borders, 0.8, dissolve=T)
border.buffer@polygons[[1]]@Polygons[[1]]@coords[, 1] = border.buffer@polygons[[1]]@Polygons[[1]]@coords[, 1]-0.1

###### Functions ######
scale_x_longitude <- function(xmin=-180, xmax=180, step=1, limits=NA, ...) {
  xbreaks <- seq(xmin,xmax,step)
  xlabels <- unlist(sapply(xbreaks, function(x) ifelse(x < 0, parse(text=paste0(x,"^o", "*W")), ifelse(x > 0, parse(text=paste0(x,"^o", "*E")),x))))
  if (is.na(limits)){
    return(scale_x_continuous("Longitude", breaks = xbreaks, labels = xlabels, expand = c(0, 0), ...))
  } else {
    return(scale_x_continuous("Longitude", breaks = xbreaks, labels = xlabels, expand = c(0, 0), limits=limits, ...))
  }
}
scale_y_latitude <- function(ymin=-90, ymax=90, step=0.5, limits=NA, ...) {
  ybreaks <- seq(ymin,ymax,step)
  ylabels <- unlist(sapply(ybreaks, function(x) ifelse(x < 0, parse(text=paste0(x,"^o", "*S")), ifelse(x > 0, parse(text=paste0(x,"^o", "*N")),x))))
  if (is.na(limits)){
    return(scale_y_continuous("Latitude", breaks = ybreaks, labels = ylabels, expand = c(0, 0), ...))
  } else {
    return(scale_y_continuous("Latitude", breaks = ybreaks, labels = ylabels, expand = c(0, 0), limits=limits, ...))
  }
}

# GTT function
temp_geom4 <- function(day,fi,a, b) {
  f=ifelse(fi==0,1e-10,fi)
  costeta= cos( (day-18 )*pi/182.5 +2^(1-sign(fi) ) *pi) 
  cosfi = cos(fi*pi/180 )
  A=cosfi
  B= (1-costeta ) * abs(sin(fi*pi/180 ) )
  x=a*A + b*B 
  return(x)}
# GTT coefficients
if(var=='tmin'){
  a <- 24.16453
  b <-  -15.71751
} else if(var=='tmax'){
  a <- 37.03043
  b <- -15.43029
} else if(var=='tmean'){
  a <- 30.419375
  b <- -15.539232
}

### Prepare data ####################################################################

load(file="ogimet/ogimet_serbia.rda") # 591898
load(file="ogimet/stations_serbia.rda")
serbia = point.in.polygon(stations_serbia$lon, stations_serbia$lat,
                          border.buffer@polygons[[1]]@Polygons[[1]]@coords[,1],
                          border.buffer@polygons[[1]]@Polygons[[1]]@coords[,2])
summary(as.factor(serbia))
# c(18,23.8,23.8,18),
# c(41.4,41.4,47,47))
stations_serbia = stations_serbia[serbia != 0, ] # 93 stations
stations_serbia$wmo_id <- as.numeric(as.character(stations_serbia$wmo_id))
stations_serbia$station_names <- as.character(stations_serbia$station_names)
names(stations_serbia)[1] <- "station_ID"
ogimet_serbia <- ogimet_serbia[ogimet_serbia$station_ID %in% stations_serbia$station_ID, ]
ogimet_serbia <- ogimet_serbia[order(ogimet_serbia$Date, ogimet_serbia$station_ID), ]

ogimet_serbia <- join(ogimet_serbia, stations_serbia, by="station_ID")

# save(ogimet_serbia, file="ogimet/ogimet_serbia08.rda")
load(file="ogimet/ogimet_serbia08.rda")

names(ogimet_serbia)
el <- c('station_ID','station_names','lat', 'lon', 'alt', 'Date')

col.names <- c('staid', 'name', 'lat', 'lon', 'elevation', 'date')

# get variable specific columns
if (var == "tmax"){
  # summary(ogimet_serbia$TemperatureCMax)
  ogimet_serbia <- ogimet_serbia[((ogimet_serbia$TemperatureCAvg < ogimet_serbia$TemperatureCMax) & (ogimet_serbia$TemperatureCAvg > ogimet_serbia$TemperatureCMin)), ]
  ogimet_serbia <- ogimet_serbia[ogimet_serbia$TemperatureCMax > ogimet_serbia$TemperatureCMin, ] # 400057 (5 were out)
  ogimet_serbia <- ogimet_serbia[, c(el, 'TemperatureCMax')]
} else if (var == "tmin"){
  # summary(ogimet_serbia$TemperatureCMin)
  ogimet_serbia <- ogimet_serbia[((ogimet_serbia$TemperatureCAvg < ogimet_serbia$TemperatureCMax) & (ogimet_serbia$TemperatureCAvg > ogimet_serbia$TemperatureCMin)), ]
  ogimet_serbia <- ogimet_serbia[ogimet_serbia$TemperatureCMax > ogimet_serbia$TemperatureCMin, ]
  ogimet_serbia <- ogimet_serbia[, c(el, 'TemperatureCMin')]
} else if (var == "tmean"){
  # summary(ogimet_serbia$TemperatureCAvg)
  ogimet_serbia <- ogimet_serbia[((ogimet_serbia$TemperatureCAvg < ogimet_serbia$TemperatureCMax) & (ogimet_serbia$TemperatureCAvg > ogimet_serbia$TemperatureCMin)), ]
  ogimet_serbia <- ogimet_serbia[ogimet_serbia$TemperatureCMax > ogimet_serbia$TemperatureCMin, ]
  ogimet_serbia <- ogimet_serbia[, c(el, 'TemperatureCAvg')]
} else if (var == "slp") {
  # summary(ogimet_serbia$PresslevHp)
  ogimet_serbia <- ogimet_serbia[, c(el, 'PresslevHp')]
} else if (var == "prcp") {
  # summary(ogimet_serbia$Precmm)
  # summary(as.factor(ifelse(ogimet_serbia$Precmm>0, 1, 0)))
  ogimet_serbia <- ogimet_serbia[, c(el, 'Precmm')]
}
names(ogimet_serbia) <- c(col.names, var)

# remove NA
ogimet_serbia <- ogimet_serbia[!is.na(ogimet_serbia[, var]), ]

# dem, twi, cdate, doy
ogimet_sp <- ogimet_serbia
coordinates(ogimet_sp) <-~ lon + lat
dem <- raster('dem_twi/dem_buff.tif')
ogimet_serbia$dem <- raster::extract(dem, ogimet_sp)
twi <- raster('dem_twi/twi_buff.tif')
ogimet_serbia$twi <- raster::extract(twi, ogimet_sp)
ogimet_serbia$cdate = floor(unclass(as.POSIXct(as.POSIXct(paste(ogimet_serbia$date), format="%Y-%m-%d")))/86400)
ogimet_serbia$doy = as.integer(strftime(as.POSIXct(paste(ogimet_serbia$date), format="%Y-%m-%d"), format = "%j"))

# ### GTT fit locally #####
# # f=ifelse(fi==0,1e-10,fi)
# # costeta= cos( (day-18 )*pi/182.5 +2^(1-sign(fi) ) *pi) 
# # cosfi = cos(fi*pi/180 )
# # A=cosfi
# # B= (1-costeta ) * abs(sin(fi*pi/180 ) )
# # x=a*A + b*B 
# 
# a.values <- cos(ogimet_serbia$lat*pi/180)
# summary(a.values)
# teta <- (ogimet_serbia$doy-18)*pi/182.5 + 2^(1-sign(ogimet_serbia$lat)) *pi
# summary(teta)
# b.values <- (1-cos(teta)) * abs(sin(ogimet_serbia$lat*pi/180))
# summary(b.values)
# x <- cbind(a.values, b.values)
# head(x)
# y <- ogimet_serbia[, var]
# mnk <- lsfit(x, y, wt = NULL, intercept = FALSE, tolerance = 1e-07,
#       yname = NULL)
# mnk$coefficients
# a <- mnk$coefficients["a.values"]
# b <- mnk$coefficients["b.values"]
# # tmax
# # a.values  b.values 
# # 41.39717 -17.50594 
### Other covariates ##############################

ogimet_serbia <- ogimet_serbia[order(ogimet_serbia$date, ogimet_serbia$staid), ]
# overlay with predictors
if (var %in% c("tmax", "tmin", "tmean")){
  # GTT, DEM, TWI, CDATE, DOY
  ogimet_serbia$gtt <- temp_geom4(ogimet_serbia$doy,ogimet_serbia$lat, a, b )
} else if (var == "slp") {
  # DEM, CDATE, DOY
} else if (var == "prcp") {
  # normals, CDATE, DOY ???
  ### TMAX and TMIN overlay ##########
  ogimet_prcp <- ogimet_serbia
  
  # tmin
  load(file = paste("ogimet/ogimet_serbia08_", "tmin", ".rda", sep=""))
  load(file=paste('ogimet/model_', "tmin", '.rda', sep=''))
  
  a <- 24.16453
  b <-  -15.71751
  ogimet_prcp$gtt <- temp_geom4(ogimet_prcp$doy,ogimet_prcp$lat, a, b)
  
  pre <- pred.rfsi(model=final.model, # RFSI model iz rfsi ili tune rfsi funkcije
                   data=ogimet_serbia, # data.frame(x,y,obs,time) | STFDF - with covariates | SpatialPointsDataFrame | SpatialPixelsDataFrame
                   zcol="tmin",
                   data.staid.x.y.time = c(1,4,3,6), # if data.frame
                   newdata=ogimet_prcp, # data.frame(x,y,time,ec1,ec2,...) | STFDF - with covariates | SpatialPointsDataFrame | SpatialPixelsDataFrame
                   newdata.staid.x.y.time = c(1,4,3,6), # if data.frame
                   output.format = "data.frame",
                   zero.tol=0,
                   # n.obs=10, # nearest obs 3 vidi iz modela
                   # time.nmax, # use all if not specified
                   s.crs=wgs84,
                   t.crs=utm34,
                   newdata.s.crs=wgs84,
                   newdata.t.crs=utm34,
                   # parallel.processing = FALSE, # doParallel - videti zbog ranger-a
                   cpus=detectCores()-1,
                   progress=F)
  
  pre$pred <- round(pre$pred, 1)
  ogimet_prcp <- join(ogimet_prcp, pre[, c(1,4,7)])
  # ogimet_prcp$tmin <- ogimet_prcp$pred
  # ogimet_prcp$pred <- NULL
  names(ogimet_prcp)[13] <- "tmin"
  
  # tmax
  load(file = paste("ogimet/ogimet_serbia08_", "tmax", ".rda", sep=""))
  load(file=paste('ogimet/model_', "tmax", '.rda', sep=''))
  
  a <- 37.03043
  b <- -15.43029
  ogimet_prcp$gtt <- temp_geom4(ogimet_prcp$doy,ogimet_prcp$lat, a, b)
  
  pre <- pred.rfsi(model=final.model, # RFSI model iz rfsi ili tune rfsi funkcije
                   data=ogimet_serbia, # data.frame(x,y,obs,time) | STFDF - with covariates | SpatialPointsDataFrame | SpatialPixelsDataFrame
                   zcol="tmax",
                   data.staid.x.y.time = c(1,4,3,6), # if data.frame
                   newdata=ogimet_prcp, # data.frame(x,y,time,ec1,ec2,...) | STFDF - with covariates | SpatialPointsDataFrame | SpatialPixelsDataFrame
                   newdata.staid.x.y.time = c(1,4,3,6), # if data.frame
                   output.format = "data.frame",
                   zero.tol=0,
                   # n.obs=10, # nearest obs 3 vidi iz modela
                   # time.nmax, # use all if not specified
                   s.crs=wgs84,
                   t.crs=utm34,
                   newdata.s.crs=wgs84,
                   newdata.t.crs=utm34,
                   # parallel.processing = FALSE, # doParallel - videti zbog ranger-a
                   cpus=detectCores()-1,
                   progress=F)
  
  pre$pred <- round(pre$pred, 1)
  ogimet_prcp <- join(ogimet_prcp, pre[, c(1,4,7)])
  # ogimet_prcp$tmax <- ogimet_prcp$pred
  # ogimet_prcp$pred <- NULL
  names(ogimet_prcp)[14] <- "tmax"
  
  ogimet_prcp$gtt <- NULL
  
  # slp
  load(file = paste("ogimet/ogimet_serbia08_", "slp", ".rda", sep=""))
  load(file=paste('ogimet/model_', "slp", '.rda', sep=''))
  
  pre <- pred.rfsi(model=final.model, # RFSI model iz rfsi ili tune rfsi funkcije
                   data=ogimet_serbia, # data.frame(x,y,obs,time) | STFDF - with covariates | SpatialPointsDataFrame | SpatialPixelsDataFrame
                   zcol="slp",
                   data.staid.x.y.time = c(1,4,3,6), # if data.frame
                   newdata=ogimet_prcp, # data.frame(x,y,time,ec1,ec2,...) | STFDF - with covariates | SpatialPointsDataFrame | SpatialPixelsDataFrame
                   newdata.staid.x.y.time = c(1,4,3,6), # if data.frame
                   output.format = "data.frame",
                   zero.tol=0,
                   # n.obs=10, # nearest obs 3 vidi iz modela
                   # time.nmax, # use all if not specified
                   s.crs=wgs84,
                   t.crs=utm34,
                   newdata.s.crs=wgs84,
                   newdata.t.crs=utm34,
                   # parallel.processing = FALSE, # doParallel - videti zbog ranger-a
                   cpus=detectCores()-1,
                   progress=F)
  
  pre$pred <- round(pre$pred, 1)
  ogimet_prcp <- join(ogimet_prcp, pre[, c(1,4,7)])
  # ogimet_prcp$slp <- ogimet_prcp$pred
  # ogimet_prcp$pred <- NULL
  names(ogimet_prcp)[14] <- "slp"
  
  ogimet_serbia <- ogimet_prcp
  rm(ogimet_prcp)
  
  # IMERG
  # time=index(stfdf@time)
  # daysNum = length(time)
  # days=gsub("-","",time,fixed=TRUE)
  dates <- sort(unique(ogimet_serbia$date))
  
  registerDoParallel(cores=detectCores()-1)
  imerg = foreach (date = dates) %dopar% { 
    ogimet_day <- ogimet_serbia[ogimet_serbia$date==date, ]
    coordinates(ogimet_day) <-~ lon+lat
    
    version <- 'V06B'
    day<-gsub("-","",date,fixed=TRUE)
    n_30 = sprintf("%04d", as.numeric(strftime(date, format = "%j")) * 30 - 30)
    r <- try(raster(paste("imerg/3B-DAY-GIS.MS.MRG.3IMERG.", day, "-S000000-E235959.", n_30, ".", version, ".total.accum.tif", sep = "")))
    if(inherits(r, "try-error")) {
      return(rep(NA, length(ogimet_day)))
    }
    r@file@nodatavalue = 29999
    
    # resample to dem 1 km
    im_res <- resample(r, dem, method="bilinear")
    
    e <- raster::extract(im_res,ogimet_day)
    return(e / 10)
  }
  stopImplicitCluster()
  # imerg <- do.call("cbind", imerg)
  # imerg <- as.vector(imerg)
  imerg = unlist(imerg)
  # imerg = imerg * 2.4 #/ 10 * 24
  imerg[imerg==2999.9] = NA
  ogimet_serbia$imerg = imerg
  
}
summary(ogimet_serbia)
save(ogimet_serbia, file = paste("ogimet/ogimet_serbia08_", var, ".rda", sep=""))

### EOBS #####

load(file = paste("ogimet/ogimet_serbia08_", var, ".rda", sep=""))

dates <- sort(unique(ogimet_serbia$date))

registerDoParallel(cores=detectCores())
eobs <- foreach(date = dates, .packages = c("raster","spacetime","gstat","rgdal","raster","doParallel","meteo")) %dopar% {
# for (date in dates) {
  # day<-gsub("-","",date,fixed=TRUE)
  band = as.integer(date - as.Date("1950-01-01")) + 1
  r <- raster(paste("eobs/", var, "_ens_mean_0.1deg_reg_v21.0e.nc", sep=""), band=band)
  # r@z[[1]]
  # get ogimet_serbia date
  ogimet_day <- ogimet_serbia[ogimet_serbia$date==date, ]
  coordinates(ogimet_day) <-~ lon+lat
  return(raster::extract(r, ogimet_day))
}
stopImplicitCluster()
eobs <- unlist(eobs)
summary(eobs)
ogimet_serbia$eobs <- eobs

# ogimet_serbia <- ogimet_serbia[complete.cases(ogimet_serbia), ]
# save(ogimet_serbia, file = paste("ogimet/ogimet_serbia08_", var, ".rda", sep=""))


### TPS ##############################

load(file = paste("ogimet/ogimet_serbia08_", var, ".rda", sep=""))

dates <- as.character(sort(unique(ogimet_serbia$date)))
r <- raster("dem_twi/dem_buff.tif")

tps.df <- 13

registerDoParallel(cores=detectCores())
tps <- foreach(date = dates, .packages = c("RSAGA", "raster","spacetime","gstat","rgdal","raster","doParallel","meteo")) %dopar% {
  day.ogimet <- ogimet_serbia[ogimet_serbia$date == date, ]
  coordinates(day.ogimet) <- c("lon", "lat")
  day.ogimet@proj4string <- wgs84
  if ((nrow(day.ogimet)-3) <= tps.df) {
    tps.df.day <- nrow(day.ogimet)-3
  } else {
    tps.df.day <- tps.df
  }
  m <- Tps(coordinates(day.ogimet), as.data.frame(day.ogimet)[, var], lon.lat = T,
           # lambda = 20)
           df=tps.df.day)
  tps <- as.vector(m$fitted.values)
  return(tps)
}
stopImplicitCluster()

tps = unlist(tps)
ogimet_serbia$tps <- tps

save(ogimet_serbia, file = paste("ogimet/ogimet_serbia08_", var, ".rda", sep=""))

### Outliers detection ##############################################
# not for SLP #

load(file = paste("ogimet/ogimet_serbia08_", var, ".rda", sep=""))

if (var=="prcp"){
  
  # remove prcp > 200 mm - impossible (max 211.1 mm in 10 October 1955 in Negotin)
  # ogimet_serbia <- ogimet_serbia[ogimet_serbia$prcp<200, ]
  
  summary(ogimet_serbia[, var])
  ogimet_serbia_test <- ogimet_serbia[ogimet_serbia$prcp>40, c("staid","name","date",var,"eobs", "imerg", "lon", "lat")]
  nrow(ogimet_serbia_test)
  treshold = 4
} else if (var=="tmax" | var=="tmin" | var=="tmean") {
  summary(ogimet_serbia[, var])
  ogimet_serbia_test <- ogimet_serbia[abs(ogimet_serbia[, var] - ogimet_serbia$eobs) > 5, c("staid","name","date",var,"eobs","tps", "lon", "lat")]
  nrow(ogimet_serbia_test)
  treshold = 10
}

br = 0
out <- c()

for(o in 1:nrow(ogimet_serbia_test)) {
  om <- ogimet_serbia_test[o, ]
  coordinates(om) <- c("lon", "lat")
  om@proj4string <- wgs84
  ogimet_date <- ogimet_serbia[ogimet_serbia$date==om$date, ] # c("staid","date",var,"eobs", "imerg", "lon", "lat")]
  coordinates(ogimet_date) <- c("lon", "lat")
  ogimet_date@proj4string <- wgs84
  ogimet_date$dist <- spDists(ogimet_date, om, longlat = T)[, 1]
  ogimet_date <- ogimet_date[order(ogimet_date$dist), ]

  ogimet_date <- ogimet_date[ogimet_date$dist>1&ogimet_date$dist<100, ]
  if (nrow(ogimet_date)==0) {
    next
  }
  
  if (var=="prcp"){
    if ((om$prcp/max(ogimet_date$prcp)) > treshold & # or mean
        (om$prcp/om$eobs) > treshold # & (om$prcp/om$imerg) > treshold
    ) {
      print(paste(om$staid, om$date, om$prcp, om$eobs, sep="  "))
      print(ogimet_date$prcp)
      print("")
      br = br+1
      out <- cbind(out, c(om$staid, om$date))
      ogimet_serbia <- ogimet_serbia[!(ogimet_serbia$staid==om$staid & ogimet_serbia$date==om$date), ]
    }
  } else if (var=="tmax" | var=="tmin" | var=="tmean") {
    if (abs(om@data[, var] - mean(ogimet_date@data[, var], na.rm=T)) > treshold & # or mean
        abs(om@data[, var] - om@data$eobs) > treshold & # & (om$prcp/om$imerg) > treshold
        abs(om@data[, var] - om@data$tps) > treshold 
    ) {
      print(paste(om@data$staid, om@data$date, om@data[, var], om@data$eobs, om@data$tps, sep="  "))
      print(ogimet_date@data[, var])
      print("")
      br = br+1
      out <- cbind(out, c(om@data$staid, om@data$date))
      ogimet_serbia <- ogimet_serbia[!(ogimet_serbia$staid==om@data$staid & ogimet_serbia$date==om@data$date), ]
    }
  }
  
}

if (var=="prcp") {
  ogimet_serbia <- ogimet_serbia[ogimet_serbia$prcp!=99, ]
  ogimet_serbia[which.max(ogimet_serbia$prcp), ]
  ogimet_serbia[ogimet_serbia$prcp==99, ]
  ogimet_serbia[ogimet_serbia$date==as.Date("2009-06-24"), ]
  summary(ogimet_serbia$prcp)
}

br # prcp 110 out # tmax 32 out # tmin 37 out # tmean 1 out
ogimet_serbia[which.max(ogimet_serbia[, var]-ogimet_serbia$tps), ]
# ogimet_serbia[which.max(ogimet_serbia[, var]), ]

### Save the data ###

if (var!="slp") {
  ogimet_serbia <- ogimet_serbia[complete.cases(ogimet_serbia), ]
}

cor(ogimet_serbia[, var], ogimet_serbia$eobs) # gtt # imerg

if (var=="prcp") {
  ogimet_serbia$prcp_cl <- as.factor(ifelse(ogimet_serbia$prcp == 0, 0, 1))
  # summary(ogimet_serbia$prcp_cl)
}

save(ogimet_serbia, file = paste("ogimet/ogimet_serbia08_", var, ".rda", sep=""))

#####################################################################################

### Tune and CV #####################################################################

# get variable specific columns
if (var == "tmax" | var == "tmin" | var == "tmean"){
  frm <- paste(var, ' ~ dem + twi + gtt + doy', sep='') # + cdate 
} else if (var == "slp") {
  frm <- paste(var, ' ~ dem + doy', sep='') # + cdate
} else if (var == "prcp") {
  frm <- paste(var, ' ~ dem + doy + imerg + tmax + tmin + slp', sep='') # + cdate
  frm.cl <- "prcp_cl ~ dem + doy + imerg + tmax + tmin + slp" # + cdate
}

load(file = paste("ogimet/ogimet_serbia08_", var, ".rda", sep=""))

# check number of stations by days
n.st <- c()
dates <- as.character(sort(unique(ogimet_serbia$date)))
rm.dates <- c()
for (t in 1:length(dates)){
  sub.ogimet <- ogimet_serbia[ogimet_serbia$date == dates[t], ]
  n.st[t] <- nrow(sub.ogimet)
  if (nrow(sub.ogimet) < 17) rm.dates <- c(rm.dates, dates[t])
}
summary(n.st)
summary(n.st<17) # we will tune it from 4:12
summary(n.st<30)
summary(n.st>50)
dates[n.st==57]

# remove days with less than 17 observations
ogimet_serbia <- ogimet_serbia[!(ogimet_serbia$date %in% as.Date(rm.dates)), ]
###

### prcp classification ####
if (var=="prcp") {
  # tune mtry and n.obs
  nos <- seq(6, 10, 1)
  min.node.size <- 10
  sample.fraction <- 0.95
  splitrule <- "gini"
  ntree <- 250
  mtry <- 3:(2+2*max(nos))
  tgrid = expand.grid(min.node.size=min.node.size, num.trees=ntree,
                      mtry=mtry, n.obs=nos, sample.fraction=sample.fraction) 
  
  # ogimet_serbia <- ogimet_serbia[order(ogimet_serbia$date, ogimet_serbia$staid), ]
  
  # chek time for all data for 1 combination
  start <- Sys.time()
  tuned.model.cl <- tune.rfsi(formula = frm.cl, # without nearest obs
                           data = ogimet_serbia,
                           data.staid.x.y.time = c(1,4,3,6),
                           zero.tol=0,
                           s.crs=wgs84,
                           t.crs=utm34,
                           use.idw=FALSE,
                           tgrid = tgrid,
                           tgrid.n = 50,
                           tune.type = "LLO",
                           k = 5, # number of folds
                           acc.metric = "Kappa",
                           cpus=detectCores()-1,
                           progress=F,
                           fit.final.model=FALSE,
                           importance = 'none', # speed up rf
                           seed = 42,
                           oob.error = FALSE, # speed up rf
                           verbose = FALSE)
  Sys.time() - start # 5.667484 hours
  tuned.model.cl
  #   min.node.size num.trees mtry n.obs sample.fraction     RMSE
  # 1            10       250   10     8            0.95 1.041758
  
  # save(tuned.model.cl, file=paste('ogimet/tuned1_', var, '_cl.rda', sep=''))
  load(file=paste('ogimet/tuned1_', var, '_cl.rda', sep=''))
  
  # tune min.node.size, sample.fraction - 30 combinations
  # tune mtry and n.obs
  nos <- tuned.model.cl$tuned.parameters$n.obs
  min.node.size <- 1:15
  sample.fraction <- seq(1, 0.65, -0.05)
  splitrule <- "gini"
  ntree <- 250
  mtry <- tuned.model.cl$tuned.parameters$mtry
  tgrid = expand.grid(min.node.size=min.node.size, num.trees=ntree,
                      mtry=mtry, n.obs=nos, sample.fraction=sample.fraction) 
  
  # chek time for all data for 1 combination
  start <- Sys.time()
  tuned.model.cl <- tune.rfsi(formula = frm.cl, # without nearest obs
                           data = ogimet_serbia,
                           data.staid.x.y.time = c(1,4,3,6),
                           zero.tol=0,
                           s.crs=wgs84,
                           t.crs=utm34,
                           use.idw=FALSE,
                           tgrid = tgrid,
                           tgrid.n = 30,
                           tune.type = "LLO",
                           k = 5, # number of folds
                           acc_param = "Kappa",
                           cpus=detectCores()-1,
                           progress=F,
                           fit.final.model=FALSE,
                           importance = 'none', # speed up rf
                           seed = 42,
                           oob.error = FALSE, # speed up rf
                           verbose = FALSE)
  Sys.time() - start # 2.941284 hours
  tuned.model.cl
  # tmax
  #     min.node.size num.trees mtry n.obs sample.fraction     RMSE
  # 102             7       250    7    10            0.91 1.805666
  
  # save(tuned.model.cl, file=paste('ogimet/tuned2_', var, '_cl.rda', sep=''))
  load(file=paste('ogimet/tuned2_', var, '_cl.rda', sep=''))
  params <- tuned.model.cl$tuned.parameters
}
#####

### Tune 1 (n.obs, mtry) ####
# tune mtry and n.obs
nos <- seq(6, 10, 1)
min.node.size <- 10
sample.fraction <- 0.95
splitrule <- "variance"
ntree <- 250
mtry <- 3:(2+2*max(nos))
idw.p <- 3
tgrid = expand.grid(min.node.size=min.node.size, num.trees=ntree,
                    mtry=mtry, n.obs=nos, sample.fraction=sample.fraction,
                    idw.p=idw.p) 

# ogimet_serbia <- ogimet_serbia[order(ogimet_serbia$date, ogimet_serbia$staid), ]

# chek time for all data for 1 combination
start <- Sys.time()
tuned.model <- tune.rfsi(formula = frm, # without nearest obs
                         data = ogimet_serbia,
                         data.staid.x.y.time = c(1,4,3,6),
                         zero.tol=0,
                         s.crs=wgs84,
                         t.crs=utm34,
                         use.idw = TRUE,
                         tgrid = tgrid,
                         tgrid.n = 50,
                         tune.type = "LLO",
                         k = 5, # number of folds
                         acc.metric = "RMSE",
                         cpus=detectCores()-1,
                         progress=F,
                         fit.final.model=FALSE,
                         importance = 'none', # speed up rf
                         seed = 42,
                         oob.error = FALSE, # speed up rf
                         verbose = FALSE)
Sys.time() - start # 6.557313 hours
tuned.model
#   min.node.size num.trees mtry n.obs sample.fraction     RMSE
# 1            10       250   10     8            0.95 1.041758

save(tuned.model, file=paste('ogimet/tuned1_', var, '.rda', sep=''))

#####

### Tune 2 (min.node.size, sample.fraction) ####
load(file=paste('ogimet/tuned1_', var, '.rda', sep=''))

# tune min.node.size, sample.fraction - 30 combinations
# tune mtry and n.obs
nos <- tuned.model$tuned.parameters$n.obs
min.node.size <- 5:15
sample.fraction <- seq(1, 0.9, -0.01)# seq(1, 0.632, -0.05) # 0.632 without / 1 with replacement
splitrule <- "variance"
ntree <- 250 # 500
mtry <- tuned.model$tuned.parameters$mtry
idw.p <- 3 # 2.5
tgrid = expand.grid(min.node.size=min.node.size, num.trees=ntree,
                    mtry=mtry, n.obs=nos, sample.fraction=sample.fraction,
                    idw.p=idw.p) 

# chek time for all data for 1 combination
start <- Sys.time()
tuned.model <- tune.rfsi(formula = frm, # without nearest obs
                         data = ogimet_serbia,
                         data.staid.x.y.time = c(1,4,3,6),
                         zero.tol=0,
                         s.crs=wgs84,
                         t.crs=utm34,
                         use.idw = T,
                         tgrid = tgrid,
                         tgrid.n = 30,
                         tune.type = "LLO",
                         k = 5, # number of folds
                         acc_param = "RMSE",
                         cpus=detectCores()-1,
                         progress=F,
                         fit.final.model=FALSE,
                         importance = 'none', # speed up rf
                         seed = 42,
                         oob.error = FALSE, # speed up rf
                         verbose = FALSE)
Sys.time() - start # 2.941284 hours
tuned.model
# tmax
#     min.node.size num.trees mtry n.obs sample.fraction     RMSE
# 102             7       250    7    10            0.91 1.805666
save(tuned.model, file=paste('ogimet/tuned2_', var, '.rda', sep=''))

#####

### Tune 3 (IDW) ####
load(file=paste('ogimet/tuned2_', var, '.rda', sep=''))

# tune idw.p - 10 combinations
nos <- tuned.model$tuned.parameters$n.obs
min.node.size <- tuned.model$tuned.parameters$min.node.size
sample.fraction <- tuned.model$tuned.parameters$sample.fraction
splitrule <- "variance"
ntree <- 250
mtry <- tuned.model$tuned.parameters$mtry
# idw.p <- seq(1,4,0.5)
idw.p <- seq(2.5,3.5,0.1)
tgrid = expand.grid(min.node.size=min.node.size, num.trees=ntree,
                    mtry=mtry, n.obs=nos, sample.fraction=sample.fraction,
                    idw.p=idw.p) 

# chek time for all data for 1 combination
start <- Sys.time()
tuned.model <- tune.rfsi(formula = frm, # without nearest obs
                         data = ogimet_serbia,
                         data.staid.x.y.time = c(1,4,3,6),
                         zero.tol=0,
                         s.crs=wgs84,
                         t.crs=utm34,
                         use.idw=TRUE,
                         tgrid = tgrid,
                         tgrid.n = 10,
                         tune.type = "LLO",
                         k = 5, # number of folds
                         acc_param = "RMSE",
                         cpus=detectCores()-1,
                         progress=F,
                         fit.final.model=FALSE,
                         importance = 'none', # speed up rf
                         seed = 42,
                         oob.error = FALSE, # speed up rf
                         verbose = FALSE)
Sys.time() - start # 2.941284 hours
tuned.model
# tmax
#     min.node.size num.trees mtry n.obs sample.fraction     RMSE
# 102             7       250    7    10            0.91 1.805666
save(tuned.model, file=paste('ogimet/tuned3_', var, '.rda', sep=''))

#####
#####################################################################################

### nested 5-fold LLOCV ####
### prcp classification CV ####
if (var=="prcp") {
  load(file=paste('ogimet/tuned2_', var, '_cl.rda', sep=''))
  params <- tuned.model.cl$tuned.parameters
  
  # cross validation on whole dataset based on tuned parameters +-2 mtry and n.obs (just a few combinations - 10)
  nos <- (params$n.obs-1):(params$n.obs+1)
  min.node.size <- (params$min.node.size-2):(params$min.node.size+2)
  sample.fraction <- seq(params$sample.fraction-0.1,params$sample.fraction+0.1, 0.05)
  splitrule <- "gini"
  ntree <- 250
  mtry <- (params$mtry-2):(params$mtry+2)
  tgrid = expand.grid(min.node.size=min.node.size, num.trees=ntree,
                      mtry=mtry, n.obs=nos, sample.fraction=sample.fraction)
  
  start <- Sys.time()
  cv.model.cl <- cv.rfsi(formula = frm.cl, # without nearest obs
                      data = ogimet_serbia,
                      data.staid.x.y.time = c(1,4,3,6),
                      output.format = "data.frame",
                      zero.tol=0,
                      s.crs=wgs84,
                      t.crs=utm34,
                      use.idw = FALSE,
                      tgrid = tgrid,
                      tgrid.n = 20,
                      tune.type = "LLO",
                      k = 5, # number of folds
                      acc_param = "Kappa",
                      cpus=detectCores()-1,
                      progress=F,
                      importance = 'none', # speed up rf
                      seed = 42,
                      oob.error = FALSE, # speed up rf
                      verbose = FALSE)
  Sys.time() - start # 3.063333 hours
  
  # save(cv.model.cl, file=paste('ogimet/cv_', var, '_cl.rda', sep=''))
  load(file=paste('ogimet/cv_', var, '_cl.rda', sep=''))
}
#####

### CV ####
load(file=paste('ogimet/tuned3_', var, '.rda', sep=''))
params <- tuned.model$tuned.parameters

# cross validation on whole dataset based on tuned parameters +-2 mtry and n.obs (just a few combinations - 10)
nos <- (params$n.obs-1):(params$n.obs+1)
min.node.size <- (params$min.node.size-2):(params$min.node.size+2) 
sample.fraction <- seq(params$sample.fraction-0.02,params$sample.fraction+0.02, 0.01)
# sample.fraction <- seq(params$sample.fraction-0.04,params$sample.fraction, 0.01)
splitrule <- "variance"
ntree <- 250
mtry <- (params$mtry-2):(params$mtry+2)
idw.p <- seq(params$idw.p-0.2,params$idw.p+0.2, 0.1)

tgrid = expand.grid(min.node.size=min.node.size, num.trees=ntree,
                    mtry=mtry, n.obs=nos, sample.fraction=sample.fraction,
                    idw.p=idw.p)

start <- Sys.time()
cv.model <- cv.rfsi(formula = frm, # without nearest obs
                    data = ogimet_serbia,
                    data.staid.x.y.time = c(1,4,3,6),
                    output.format = "data.frame",
                    zero.tol=0,
                    s.crs=wgs84,
                    t.crs=utm34,
                    use.idw = TRUE,
                    tgrid = tgrid,
                    tgrid.n = 20,
                    tune.type = "LLO",
                    k = 5, # number of folds
                    acc_param = "RMSE",
                    cpus=detectCores()-1,
                    progress=F,
                    importance = 'none', # speed up rf
                    seed = 42,
                    oob.error = FALSE, # speed up rf
                    verbose = FALSE)
Sys.time() - start # 7.144371 hours

save(cv.model, file=paste('ogimet/cv_', var, '.rda', sep=''))

### CV results ####

load(file=paste('ogimet/cv_', var, '.rda', sep=''))

if (var=="prcp") {
  load(file=paste('ogimet/cv_', var, '_cl.rda', sep=''))
  # join cv prcp results
  names(cv.model.cl)[5:7] <- c("folds_cl", "obs_cl", "pred_cl")
  cv.model <- join(cv.model.cl, cv.model, by = c("staid", "lon", "lat", "date"))
  identical(cv.model$folds_cl, cv.model$folds) # check
  cv.model$pred <- ifelse(cv.model$pred_cl==0, 0, cv.model$pred)
}

# cv.model1 <- cv.model[order(cv.model$date, cv.model$staid), ]
# identical(cv.model1$obs, ogimet_serbia$tmax)

# serbia only
serbia = point.in.polygon(cv.model$lon, cv.model$lat,
                          borders@polygons[[1]]@Polygons[[1]]@coords[,1],
                          borders@polygons[[1]]@Polygons[[1]]@coords[,2])
summary(as.factor(serbia))
cv.model <- cv.model[serbia != 0, ]

if (var=="prcp") {
  library(caret)
  cm <- confusionMatrix(cv.model$pred_cl, cv.model$obs_cl)
  cm
  # Prediction      0      1
  #          0 189974  19643
  #          1  14623  97150
  # Kappa : 0.7674
  
  # summary(cv.model$obs_cl==0 & cv.model$pred_cl==1)
  # summary(cv.model$pred[cv.model$obs_cl==0 & cv.model$pred_cl==1])
  length(cv.model$pred[cv.model$obs_cl==0 & cv.model$pred > 1])
  round(length(cv.model$pred[cv.model$obs_cl==0 & cv.model$pred > 1])/
    length(cv.model$obs_cl[cv.model$obs_cl==0]) * 100, 2)
  
  # summary(cv.model$obs_cl==1 & cv.model$pred_cl==0)
  # summary(cv.model$obs[cv.model$obs_cl==1 & cv.model$pred_cl==0])
  round(length(cv.model$pred[cv.model$pred_cl==0 & cv.model$obs > 1])/
          length(cv.model$obs_cl[cv.model$obs_cl==1]) * 100, 2)
}

## RMSE
sqrt(mean((cv.model$obs - cv.model$pred)^2, na.rm = T))
# 1.667313
## CCC
DescTools::CCC(cv.model$obs, cv.model$pred, ci = "z-transform", conf.level = 0.95, na.rm=TRUE)$rho.c
# 0.9869776
## MAE
mean(abs(cv.model$obs - cv.model$pred), na.rm=TRUE)
# 1.091398
## R2
1 - (t(cv.model$obs - cv.model$pred) %*% (cv.model$obs - cv.model$pred)) / (t(cv.model$obs - mean(cv.model$obs)) %*% (cv.model$obs - mean(cv.model$obs)))
# 0.974496
# ## ME
# mean((cv.model$obs - cv.model$pred), na.rm=TRUE)

#####################################################################################

### Create model ####################################################################

load(file = paste("ogimet/ogimet_serbia08_", var, ".rda", sep=""))

# get variable specific columns
if (var == "tmax" | var == "tmin" | var == "tmean"){
  frm <- paste(var, ' ~ dem + twi + gtt + doy', sep='') # + cdate 
} else if (var == "slp") {
  frm <- paste(var, ' ~ dem + doy', sep='') # + cdate
} else if (var == "prcp") {
  frm <- paste(var, ' ~ dem + doy + imerg + tmax + tmin + slp', sep='') # + cdate
  frm.cl <- "prcp_cl ~ dem + doy + imerg + tmax + tmin + slp" # + cdate
}

if (var=="prcp") {
  load(file=paste('ogimet/tuned2_', var, '_cl.rda', sep=''))
  params <- tuned.model.cl$tuned.parameters
  final.model.cl <- rfsi(formula = frm.cl, # without nearest obs
                         data = ogimet_serbia,
                         data.staid.x.y.time = c(1,4,3,6),
                         zero.tol=0,
                         n.obs = params$n.obs,
                         mtry=params$mtry,
                         sample.fraction=params$sample.fraction,
                         min.node.size=params$min.node.size,
                         num.trees=params$num.trees,
                         s.crs=wgs84,
                         t.crs=utm34,
                         cpus=detectCores()-1,
                         progress=TRUE,
                         importance = 'impurity',
                         seed = 42)
  final.model.cl
  final.model.cl$prediction.error
  
  xlP.g <- as.list(round(final.model.cl$variable.importance))
  df = t(data.frame(xlP.g[order(unlist(xlP.g), decreasing=TRUE)])) # [1:15,]/100000
  print(df)
  
  # save(final.model.cl, file=paste('ogimet/model_', var, '_cl.rda', sep=''))
  load(file=paste('ogimet/model_', var, '_cl.rda', sep=''))
}

load(file=paste('ogimet/tuned3_', var, '.rda', sep=''))
params <- tuned.model$tuned.parameters
# final model based on the upper tuned parameters
final.model <- rfsi(formula = frm, # without nearest obs
                    data = ogimet_serbia,
                    data.staid.x.y.time = c(1,4,3,6),
                    zero.tol=0,
                    n.obs= params$n.obs,
                    mtry=params$mtry,
                    sample.fraction=params$sample.fraction,
                    min.node.size=params$min.node.size,
                    num.trees=250,
                    s.crs=wgs84,
                    t.crs=utm34,
                    use.idw = TRUE,
                    idw.p = params$idw.p,
                    cpus=detectCores(),
                    progress=TRUE,
                    importance = 'impurity',
                    seed = 42)

save(final.model, file=paste('ogimet/model_', var, '.rda', sep=''))
# load(file=paste('ogimet/model_', var, '.rda', sep=''))

final.model
final.model$r.squared
#      tmax      prcp
# 0.9852171 0.6711121
sqrt(final.model$prediction.error)
#    tmax      prcp
# 1.248015 3.361527

xlP.g <- as.list(round(final.model$variable.importance))
df = t(data.frame(xlP.g[order(unlist(xlP.g), decreasing=TRUE)])) # [1:15,]/100000
print(df)

#####################################################################################

### Daily predictions - all ######################################################################

dir.create("data")
dir.create("data/ltm_wgs84")
dir.create("data/ltm_utm34")
dir.create("data/ann_wgs84")
dir.create("data/ann_utm34")
dir.create("data/mon_wgs84")
dir.create("data/mon_utm34")

for(year in 2000:2019) {
  dir.create(paste("data/day_", year, "_wgs84", sep=""))
  dir.create(paste("data/day_", year, "_utm34", sep=""))
}

for (var in c("tmax", "tmin", "tmean", "slp", "prcp")) {
  load(file = paste("ogimet/ogimet_serbia08_", var, ".rda", sep=""))
  load(file=paste('ogimet/model_', var, '.rda', sep=''))
  assign(paste("ogimet_", var, sep=""), ogimet_serbia)
  assign(paste("final_model_", var, sep=""), final.model)
}

load(file=paste('ogimet/model_', var, '_cl.rda', sep=''))
assign(paste("final_model_", var, "_cl", sep=""), final.model.cl)

dates <- sort(unique(ogimet_tmax$date))
# dates <- sort(unique(ogimet_prcp$date))

dem <- readGDAL("dem_twi/dem.tif")
r = raster(dem)

names(dem) = 'dem'
dem$twi<- readGDAL('dem_twi/twi.tif')$band1
gg <- as(dem, "SpatialPixelsDataFrame")
gg@proj4string = CRS("+proj=longlat +datum=WGS84")
rm(dem) 
gg$staid <- 1:length(gg)
gg_plot=gg
gg_raster <- raster(gg_plot)
gg_utm34 <- as(projectRaster(from=gg_raster, crs=utm34@projargs), "SpatialPixelsDataFrame")

# dates <- c('2000-01-01', '2005-03-01', '2010-06-01', '2015-09-01', '2019-11-01')
# dates <- dates[1583:7305]
# dates <- c('2000-06-03', '2000-06-07')

start = Sys.time()
registerDoParallel(cores=detectCores())
foreach(date = dates, .packages = c("raster","spacetime","gstat","rgdal","raster","doParallel","meteo")) %dopar% {
  # for (date in dates) {
  day<-gsub("-","",date,fixed=TRUE)
  year<-substr(date, 1, 4)
  
  df <- gg
  df <- as.data.frame(df)
  lon.lat <- 4:5
  names(df)[lon.lat] <- c("lon", "lat")
  
  df$doy <- as.integer(strftime(as.POSIXct(paste(date), format="%Y-%m-%d"), format = "%j"))
  df$cdate <- floor(unclass(as.POSIXct(as.POSIXct(paste(date), format="%Y-%m-%d")))/86400)
  
  for (var in c("tmax", "tmin", "tmean", "slp", "prcp")){
    ogimet_serbia <- get(paste("ogimet_", var, sep=""))
    obs_st <- ogimet_serbia[ogimet_serbia$date==date, ]
    obs_st <-  obs_st[!is.na(obs_st[, var]),]
    # if (nrow(obs_st) == 0) {
    if (nrow(obs_st) < 10) { # because of prcp n.obs = 9 (one day have < 11)
      next
    }
    
    if (var == "tmax" | var == "tmin" | var == "tmean"){
      if(var=='tmin'){
        a <- 24.16453
        b <-  -15.71751
      } else if(var=='tmax'){
        a <- 37.03043
        b <- -15.43029
      } else if(var=='tmean'){
        a <- 30.419375
        b <- -15.539232
      }
      df$gtt <- temp_geom4(df$doy,gg@coords[,2], a, b)
    }
    
    if (var=="prcp"){
      version <- 'V06B'
      n_30 = sprintf("%04d", as.numeric(strftime(date, format = "%j")) * 30 - 30)
      imerg <- raster(paste("imerg/3B-DAY-GIS.MS.MRG.3IMERG.", day, "-S000000-E235959.", n_30, ".", version, ".total.accum.tif", sep=""))
      im_res <- resample(imerg, r, method="bilinear")
      df$imerg <- raster::extract(im_res, gg_plot)/10
      # df$imerg <- raster::extract(imerg, gg_plot)/10
      # prediction
      pre <- pred.rfsi(model=get(paste("final_model_", var, "_cl", sep="")),
                       data=obs_st,
                       zcol="prcp_cl",
                       data.staid.x.y.time = c(1,4,3,NA),
                       newdata=df,
                       newdata.staid.x.y.time = c(3,lon.lat,NA),
                       output.format = "data.frame",
                       zero.tol=0,
                       s.crs=wgs84,
                       t.crs=utm34,
                       newdata.s.crs=wgs84,
                       cpus=detectCores()-1,
                       progress=F)
      df[, "prcp_cl"] <- pre$pred
      
    }
    
    # prediction
    pre <- pred.rfsi(model=get(paste("final_model_", var, sep="")),
                     data=obs_st,
                     zcol=var,
                     data.staid.x.y.time = c(1,4,3,NA),
                     newdata=df,
                     newdata.staid.x.y.time = c(3,lon.lat,NA),
                     output.format = "data.frame",
                     zero.tol=0,
                     s.crs=wgs84,
                     t.crs=utm34,
                     newdata.s.crs=wgs84,
                     cpus=detectCores()-1,
                     progress=F)
    
    pre <- round(pre$pred * 10)
    if (var=="prcp") {
      pre <- ifelse(df$prcp_cl==0, 0, pre)
    }
    df[, var] <- pre/10 # because it was trained with one decimal!
    
    gg_plot@data$pred=pre
    p <- try( raster( as(gg_plot["pred"], "SpatialPixelsDataFrame")) )
    if(inherits(p, "try-error")) {
      p <- rasterize( gg_plot, r, "pred")
    }
    # plot(p)
    # day_2000_wgs84/tmax_day_20000101_wgs84.tif
    writeRaster(p, paste("data/day_", year, "_wgs84/",
                         var, "_day_", day, "_wgs84.tif", sep = ""), "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
    p_utm34 <- projectRaster(from=p, crs=utm34)
    # plot(p_utm34)
    writeRaster(p_utm34, paste("data/day_", year, "_utm34/",
                               var, "_day_", day, "_utm34.tif", sep = ""), "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
  }
  
}
stopImplicitCluster()
Sys.time()-start # 15 hours

#####################################################################################

### E-OBS comparison calculation ######################################################################

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

for (var in c("tmax", "tmin", "tmean", "slp", "prcp")) {
  
  if (var == "prcp") {
    dates <- seq(as.Date("2000-06-01"), as.Date("2019-12-31"), by="day")
  } else {
    dates <- seq(as.Date("2000-01-01"), as.Date("2019-12-31"), by="day")
  }
  daysNum <- length(dates)
  
  registerDoParallel(cores=detectCores())
  eobs.list <- foreach(d = 1:daysNum,
                     .packages = c("raster","spacetime","gstat","rgdal","raster","doParallel","meteo"),
                     .combine='comb', .multicombine=TRUE,
                     .init=list(list(), list())) %dopar% {
  # for(d in 1:daysNum) {
    date <- dates[d]
    year<-substr(date, 1, 4)
    day<-gsub("-","",date,fixed=TRUE)
    r <- raster(paste("data/day_", year, "_wgs84/",
                      var, "_day_", day, "_wgs84.tif", sep = ""))
    band = as.integer(date - as.Date("1950-01-01")) + 1
    eobs <- raster(paste("eobs/", var, "_ens_mean_0.1deg_reg_v21.0e.nc", sep=""), band=band)
    # eobs@z
    eobs <- crop(eobs, r)
    eobs <- mask(eobs, borders)
    # resampling
    r <- resample(r, eobs)
    r <- mask(r, borders)
    return(list(values(eobs), values(r) / 10))
  }
  stopImplicitCluster()
  eobs.df <- eobs.list[[1]]
  eobs.df <- do.call("rbind", eobs.df)
  resample.df <- eobs.list[[2]]
  resample.df <- do.call("rbind", resample.df)
  
  row.names(eobs.df) <- as.character(dates)
  row.names(resample.df) <- as.character(dates)
  save(eobs.df, file = paste("ogimet/", var, "_eobs_values.rda", sep=""))
  save(resample.df, file = paste("ogimet/", var, "_res_values.rda", sep=""))
}

### ASV validation ######################################################################

load("ogimet/asv_df.rda")
asv.df <- asv.df[order(asv.df$date), ]
dates <- sort(unique(asv.df$date))

for (var in c("tmax", "tmin", "tmean", "slp", "prcp")) {
  print("###################")
  print(var)
  
  registerDoParallel(cores=detectCores())
  meteo.serbia <- foreach(date = dates, .packages = c("raster","spacetime","gstat","rgdal","raster","doParallel","meteo")) %dopar% {
    day<-gsub("-","",date,fixed=TRUE)
    year<-substr(date, 1, 4)
    r <- raster(paste("data/day_", year, "_wgs84/",
                      var, "_day_", day, "_wgs84.tif", sep = ""))
    # r@z[[1]]
    asv.day <- asv.df[asv.df$date==date, ]
    coordinates(asv.day) <-~ lon+lat
    return(raster::extract(r, asv.day)/10)
  }
  stopImplicitCluster()
  meteo.serbia <- unlist(meteo.serbia)
  asv.df[, paste(var, "_ms", sep="")] <- meteo.serbia
  
}
# save(asv.df, file = "ogimet/asv_df_ms.rda")

load(file = "ogimet/asv_df_ms.rda")
for (var in c("tmax", "tmin", "tmean", "prcp")) {
  print("###################")
  print(var)
  print(acc_param_fun(asv.df[, var], asv.df[, paste(var, "_ms", sep="")], "RMSE"))
}

### Agg monthly, yearly, LTM day ######################################################################

dates.sub <- substr(seq(as.Date("2000-01-01"), as.Date("2000-12-31"), by="day"), 6, 10)

for (var in c("tmax", "tmin", "tmean", "slp", "prcp")) {
  print("###################")
  print(var)
  
  if (var == "prcp") {
    dates <- seq(as.Date("2000-06-01"), as.Date("2019-12-31"), by="day")
  } else {
    dates <- seq(as.Date("2000-01-01"), as.Date("2019-12-31"), by="day")
  }
  daysNum <- length(dates)
  
  registerDoParallel(cores=detectCores())
  all.days.list <- foreach(d = 1:daysNum, .packages = c("raster","spacetime","gstat","rgdal","raster","doParallel","meteo")) %dopar% {
    
    date <- dates[d]
    year<-substr(date, 1, 4)
    day<-gsub("-","",date,fixed=TRUE)
    r <- raster(paste("data/day_", year, "_wgs84/",
                      var, "_day_", day, "_wgs84.tif", sep = ""))
    return(values(r) / 10)
    
  }
  stopImplicitCluster()
  all.days.df <- do.call("rbind", all.days.list)
  rm(all.days.list)
  gc();gc()
  
  row.names(all.days.df) <- as.character(dates)
  # save(all.days.df, file = paste("ogimet/", var, "_all_days_values.rda", sep=""))
  
  r <- raster(paste("data/day_", 2019, "_wgs84/",
                    var, "_day_", 20190101, "_wgs84.tif", sep = ""))
  
  ### annual
  ### tmax_ann_2000.tif
  years <- 2000:2019
  for (year in years) {
    temp.year <- all.days.df[as.integer(substr(row.names(all.days.df), 1, 4)) == year, ]
    if (var == "prcp") {
      fnc <- "sum"
    } else {
      fnc <- "mean"
    }
    year.mean <- apply(temp.year, 2, fnc)
    # summary(year.mean)
    values(r) <- round(year.mean*10)
    
    # tmax_ann_2000_wgs84.tif
    writeRaster(r, paste("data/ann_wgs84/",
                         var, "_ann_", year, "_wgs84.tif", sep = ""),
                "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
    p_utm34 <- projectRaster(from=r, crs=utm34)
    # plot(p_utm34)
    writeRaster(p_utm34, paste("data/ann_utm34/",
                               var, "_ann_", year, "_utm34.tif", sep = ""),
                "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)

    ### monthly
    ### tmax_mon_200001.tif
    if (var=="prcp" & year==2000) {
      months <- sprintf("%02d", 6:12)
    } else {
      months <- sprintf("%02d", 1:12)
    }
    
    for (month in months) {
      temp.month <- temp.year[substr(row.names(temp.year), 6, 7) == month, ]

      month.mean <- apply(temp.month, 2, fnc)
      # summary(month.mean)
      values(r) <- round(month.mean*10)
      
      # tmax_mon_200001_wgs84.tif
      writeRaster(r, paste("data/mon_wgs84/",
                           var, "_mon_", year, month, "_wgs84.tif", sep = ""),
                  "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
      p_utm34 <- projectRaster(from=r, crs=utm34)
      # plot(p_utm34)
      writeRaster(p_utm34, paste("data/mon_utm34/",
                                 var, "_mon_", year, month, "_utm34.tif", sep = ""),
                  "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)

    }
  }
  
  ### ltm day ###
  # tmax_ltm_day_20000101.tif

  for (day in dates.sub) {
    
    temp.day <- all.days.df[substr(row.names(all.days.df), 6, 10) == day, ]
    fnc <- "mean"
    day.mean <- apply(temp.day, 2, fnc)
    # summary(day.mean)
    values(r) <- round(day.mean*10)
    
    # tmax_ltm_day_0101_wgs84.tif
    writeRaster(r, paste("data/ltm_wgs84/",
                         var, "_ltm_day_", gsub("-","",day,fixed=TRUE), "_wgs84.tif", sep = ""),
                "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
    p_utm34 <- projectRaster(from=r, crs=utm34)
    # plot(p_utm34)
    writeRaster(p_utm34, paste("data/ltm_utm34/",
                               var, "_ltm_day_", gsub("-","",day,fixed=TRUE), "_utm34.tif", sep = ""),
                "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
    
  }
  
  rm(all.days.df)
  gc();gc()
}

### LTM mon, ann ######################################################################

months <- sprintf("%02d", 1:12)

for (var in c("tmax", "tmin", "tmean", "slp", "prcp")) {
  print("###################")
  print(var)
  
  r <- raster(paste("data/day_", 2019, "_wgs84/",
                    var, "_day_", 20190101, "_wgs84.tif", sep = ""))
  
  # ltm ann
  # prcp_ann_2000_wgs84.tif
  ann.files <- list.files("data/ann_wgs84", full.names = T)
  year.files <- ann.files[sapply(ann.files, function(x) strsplit(strsplit(x, "/")[[1]][3], "_")[[1]][1]) == var]
  year.df <- c()
  for (yf in year.files) {
    year.df <- rbind(year.df, values(raster(yf))/10)
  }
  fnc <- "mean"
  year.mean <- apply(year.df, 2, fnc)
  # summary(day.mean)
  values(r) <- round(year.mean*10)
  
  # tmax_ltm_ann_wgs84.tif
  writeRaster(r, paste("data/ltm_wgs84/",
                       var, "_ltm_ann_wgs84.tif", sep = ""),
              "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
  p_utm34 <- projectRaster(from=r, crs=utm34)
  # plot(p_utm34)
  writeRaster(p_utm34, paste("data/ltm_utm34/",
                             var, "_ltm_ann_utm34.tif", sep = ""),
              "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
  
  for (month in months) {
    # ltm mon
    # tmax_ltm_mon_01.tif
    mon.files <- list.files("data/mon_wgs84", full.names = T)
    month.files <- mon.files[sapply(mon.files, function(x) strsplit(strsplit(x, "/")[[1]][3], "_")[[1]][1]) == var &
                               sapply(mon.files, function(x) substr(strsplit(strsplit(x, "/")[[1]][3], "_")[[1]][3], 5, 6)) == month]
    month.df <- c()
    for (mf in month.files) {
      month.df <- rbind(month.df, values(raster(mf))/10)
    }
    month.mean <- apply(month.df, 2, fnc)
    # summary(day.mean)
    values(r) <- round(month.mean*10)
    
    # tmax_ltm_mon_01_wgs84.tif
    writeRaster(r, paste("data/ltm_wgs84/",
                         var, "_ltm_mon_", month, "_wgs84.tif", sep = ""),
                "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
    p_utm34 <- projectRaster(from=r, crs=utm34)
    # plot(p_utm34)
    writeRaster(p_utm34, paste("data/ltm_utm34/",
                               var, "_ltm_mon_", month, "_utm34.tif", sep = ""),
                "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
    
  }

}

### ZIP folders ######################################################################

setwd("data/")
zip.dirs <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE)[-1]
for (zd in zip.dirs) {
  # zd=zip.dirs[1]
  files2zip <- dir(zd, full.names = TRUE)
  zip.name <- strsplit(zd, "/")[[1]][2]
  zip(zipfile = zip.name, files = files2zip)
}

######
### Plots ######################################################################
### stations with DEM ######################################################################

neighbours <- readOGR(dsn="borders/neighbours.shp", layer="neighbours") #, driver = "GeoJSON")
neighbours <- crop(neighbours, extent(18, 24, 41, 47))
load("ogimet/ogimet_serbia08_tmax.rda")
stations <- unique(ogimet_serbia[, 1:5])
st <- stations
coordinates(st) <- c("lon", "lat")
# bbox(st)
# plot(st)
min(stations$lat); max(stations$lat);
min(stations$lon); max(stations$lon);
dem_raster <- raster("dem_twi/dem_buff.tif") 
# hillshade
slope = terrain(dem_raster, opt='slope') # nagib
aspect = terrain(dem_raster, opt='aspect') # orjentacija
hill_raster = hillShade(slope, aspect, 45, 270) # senka deklinacija 45 , azimut 270

r <- as(dem_raster, "SpatialPixelsDataFrame")
hill <- as(hill_raster, "SpatialPixelsDataFrame")

theme = theme_set(theme_minimal())
sta_dem_plot <- ggplot() + # watch out for attribute name color order
  geom_raster(data=as.data.frame(hill), aes(x=x, y=y, fill=layer), alpha=1) +
  scale_fill_gradientn(colours=grey(0:100/100),breaks=layer,guide="none") +
  new_scale_fill() + 
  geom_raster(data=as.data.frame(r), aes(x=x, y=y, fill=dem_buff), alpha=0.35) +
  scale_fill_gradientn(colours = terrain.colors(12, alpha=0.35), name = "DEM [m]") +
  geom_polygon(data = borders, aes(x = long, y = lat, group = group), alpha = 0, color = "black", fill=NA, size = 0.5) +
  geom_polygon(data = neighbours, aes(x = long, y = lat, group = group), alpha = 0, color = "black", fill=NA, size = 0.1) +
  geom_point(data = stations, aes(x = lon, y = lat, color = "red", shape = as.factor("staid")), size = 0.8) +
  theme(plot.title = element_text(hjust = 8),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        text = element_text(size = 8),
        # legend.key.size= unit(0.2, "cm"),
        # legend.margin = unit(0, "cm"),
        legend.title = element_text(size=8, face="bold"),
        legend.text=element_text(size=8),
        legend.position = c(1.25,.4),
        plot.margin=unit(c(5.5,80,5.5,5.5),"points")) +
  labs(x = "Longitude", y = "Latitude") +
  scale_colour_manual(name = "Stations",
                      labels = c("SYNOP"),
                      values = c("red"="red")) +
  scale_shape_manual(name = "Stations",
                     labels = c("SYNOP"),
                     values = c(17)) +
  scalebar(x.min = 18, x.max = 24,
           y.min = 41, y.max = 47,
           st.size = 2, location="bottomleft", st.dist=0.03, border.size=0.3,
           dist = 75, dist_unit = "km",
           transform = T, model = "WGS84",
           anchor=c(x=18.3, y=41.6)) +
  north(x.min = 18, x.max = 24,
        y.min = 41, y.max = 47,
        scale = 0.2, symbol = 3, location="topright",
        anchor=c(x=23.8, y=46.5)) +
  scale_x_longitude(xmin=18, xmax=24, step=1, limits=c(18, 24)) +
  scale_y_latitude(ymin=41, ymax=47, step=1, limits=c(41, 47)) +
  # xlim(18, 23) + ylim(41, 47) +
  coord_fixed()
sta_dem_plot

### Fig1 ###
# tiff("plot/Fig1.tiff", width = 100, height = 70, units = 'mm', res = 600, compression = "lzw")
jpeg("plot/Fig1.jpeg", width = 100, height = 70, units = 'mm', res = 600)
sta_dem_plot
dev.off()

### daily predictions 2014-07-27 ######################################################################

bbox(borders)

tmax <- raster("data/day_2014_wgs84/tmax_day_20140727_wgs84.tif")
values(tmax) <- values(tmax)/10
tmax <- as(tmax, "SpatialPixelsDataFrame")
names(tmax) <- "tmax"

tmax_plot <- ggplot() + # watch out for attribute name color order
  geom_raster(data=as.data.frame(tmax), aes(x=x, y=y, fill=tmax), alpha=0.8) +
  scale_fill_gradientn(colours = colorRampPalette(c("white", "orange", "red"))(30), name = expression(paste("Temp [",degree,"C]", sep = ""))) +
  geom_polygon(data = borders, aes(x = long, y = lat, group = group), alpha = 0, color = "black", fill=NA, size = 0.1) +
  theme(plot.title = element_text(hjust = 8, face="bold"),
        axis.text = element_text(size = 8),
        axis.title = element_blank(),#element_text(size = 8),
        text = element_text(size = 8),
        legend.position = "bottom",
        legend.direction = "horizontal",
        # legend.key.size= unit(0.2, "cm"),
        # legend.margin = unit(0, "cm"),
        legend.title = element_text(size=8), #, face="bold"),
        legend.text=element_text(size=8)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text = element_text(size = 8)) + # changes axis labels
  theme(axis.title = element_text(size = 8)) + # change axis titles
  # theme(text = element_text(size = 8))+
  ggtitle("Tmax") +
  labs(x = "Longitude", y = "Latitude") +
  coord_fixed() +
  scale_x_longitude(xmin=19, xmax=23, step=1, limits=c(18.81, 23.01)) +
  scale_y_latitude(ymin=42, ymax=46, step=1, limits=c(41.85, 46.2))

tmean <- raster("data/day_2014_wgs84/tmean_day_20140727_wgs84.tif")
values(tmean) <- values(tmean)/10
tmean <- as(tmean, "SpatialPixelsDataFrame")
names(tmean) <- "tmean"

tmean_plot <- ggplot() + # watch out for attribute name color order
  geom_raster(data=as.data.frame(tmean), aes(x=x, y=y, fill=tmean), alpha=0.8) +
  scale_fill_gradientn(colours = colorRampPalette(c("white", "orange", "red"))(30), name = expression(paste("Temp [",degree,"C]", sep = ""))) +
  geom_polygon(data = borders, aes(x = long, y = lat, group = group), alpha = 0, color = "black", fill=NA, size = 0.1) +
  theme(plot.title = element_text(hjust = 8, face="bold"),
        axis.text = element_text(size = 8),
        axis.title = element_blank(),#element_text(size = 8),
        text = element_text(size = 8),
        legend.position = "bottom",
        legend.direction = "horizontal",
        # legend.key.size= unit(0.2, "cm"),
        # legend.margin = unit(0, "cm"),
        legend.title = element_text(size=8), #, face="bold"),
        legend.text=element_text(size=8)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text = element_text(size = 8)) + # changes axis labels
  theme(axis.title = element_text(size = 8)) + # change axis titles
  # theme(text = element_text(size = 8))+
  ggtitle("Tmean") +
  labs(x = "Longitude", y = "Latitude") +
  coord_fixed() +
  scale_x_longitude(xmin=19, xmax=23, step=1, limits=c(18.81, 23.01)) +
  scale_y_latitude(ymin=42, ymax=46, step=1, limits=c(41.85, 46.2))

slp <- raster("data/day_2014_wgs84/slp_day_20140727_wgs84.tif")
values(slp) <- values(slp)/10
slp <- as(slp, "SpatialPixelsDataFrame")
names(slp) <- "slp"

slp_plot <- ggplot() + # watch out for attribute name color order
  geom_raster(data=as.data.frame(slp), aes(x=x, y=y, fill=slp), alpha=0.8) +
  scale_fill_continuous(name = "SLP [mbar]", high = "#132B43", low = "#56B1F7") +
  geom_polygon(data = borders, aes(x = long, y = lat, group = group), alpha = 0, color = "black", fill=NA, size = 0.1) +
  theme(plot.title = element_text(hjust = 8, face="bold"),
        axis.text = element_text(size = 8),
        axis.title = element_blank(),#element_text(size = 8),
        text = element_text(size = 8),
        legend.position = "bottom",
        legend.direction = "horizontal",
        # legend.key.size= unit(0.2, "cm"),
        # legend.margin = unit(0, "cm"),
        legend.title = element_text(size=8), #, face="bold"),
        legend.text=element_text(size=8)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text = element_text(size = 8)) + # changes axis labels
  theme(axis.title = element_text(size = 8)) + # change axis titles
  # theme(text = element_text(size = 8))+
  ggtitle("SLP") +
  labs(x = "Longitude", y = "Latitude") +
  coord_fixed() +
  scale_x_longitude(xmin=19, xmax=23, step=1, limits=c(18.81, 23.01)) +
  scale_y_latitude(ymin=42, ymax=46, step=1, limits=c(41.85, 46.2))

tmin <- raster("data/day_2014_wgs84/tmin_day_20140727_wgs84.tif")
values(tmin) <- values(tmin)/10
tmin <- as(tmin, "SpatialPixelsDataFrame")
names(tmin) <- "tmin"
# plot(tmin)

theme = theme_set(theme_minimal())
tmin_plot <- ggplot() + # watch out for attribute name color order
  geom_raster(data=as.data.frame(tmin), aes(x=x, y=y, fill=tmin), alpha=0.8) +
  scale_fill_gradientn(colours = colorRampPalette(c("white", "orange", "red"))(30), name = expression(paste("Temp [",degree,"C]", sep = ""))) +
  geom_polygon(data = borders, aes(x = long, y = lat, group = group), alpha = 0, color = "black", fill=NA, size = 0.1) +
  theme(plot.title = element_text(hjust = 8, face="bold"),
        axis.text = element_text(size = 8),
        axis.title = element_blank(),#element_text(size = 8),
        text = element_text(size = 8),
        legend.position = "bottom",
        legend.direction = "horizontal",
        # legend.key.size= unit(0.2, "cm"),
        # legend.margin = unit(0, "cm"),
        legend.title = element_text(size=8), #, face="bold"),
        legend.text=element_text(size=8)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text = element_text(size = 8)) + # changes axis labels
  theme(axis.title = element_text(size = 8)) + # change axis titles
  # theme(text = element_text(size = 8))+
  ggtitle("Tmin") +
  labs(x = "Longitude", y = "Latitude") +
  coord_fixed() +
  scale_x_longitude(xmin=19, xmax=23, step=1, limits=c(18.81, 23.01)) +
  scale_y_latitude(ymin=42, ymax=46, step=1, limits=c(41.85, 46.2)) +
  scalebar(x.min = 19, x.max = 23,
           y.min = 42, y.max = 46,
           st.size = 3, location="bottomleft", st.dist=0.03, border.size=0.3,
           dist = 50, dist_unit = "km",
           transform = T, model = "WGS84",
           anchor=c(x=19, y=42.1)) +
  north(x.min = 19, x.max = 23,
        y.min = 42, y.max = 46,
        scale = 0.25, symbol = 3, location="topright",
        anchor=c(x=22.8, y=46))
tmin_plot

prcp <- raster("data/day_2014_wgs84/prcp_day_20140727_wgs84.tif")
values(prcp) <- values(prcp)/10
prcp <- as(prcp, "SpatialPixelsDataFrame")
names(prcp) <- "prcp"
# plot(prcp)

prcp_plot <- ggplot() + # watch out for attribute name color order
  geom_raster(data=as.data.frame(prcp), aes(x=x, y=y, fill=prcp), alpha=0.8) +
  scale_fill_gradientn(colours = colorRampPalette(c("white", "#56B1F7", "#132B43"))(30), name = "PRCP [mm/day]") +
  geom_polygon(data = borders, aes(x = long, y = lat, group = group), alpha = 0, color = "black", fill=NA, size = 0.1) +
  theme(plot.title = element_text(hjust = 8, face="bold"),
        axis.text = element_text(size = 8),
        axis.title = element_blank(),#element_text(size = 8),
        text = element_text(size = 8),
        legend.position = "bottom",
        legend.direction = "horizontal",
        # legend.key.size= unit(0.2, "cm"),
        # legend.margin = unit(0, "cm"),
        legend.title = element_text(size=8), #, face="bold"),
        legend.text=element_text(size=8)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text = element_text(size = 8)) + # changes axis labels
  theme(axis.title = element_text(size = 8)) + # change axis titles
  # theme(text = element_text(size = 8))+
  ggtitle("PRCP") +
  labs(x = "Longitude", y = "Latitude") +
  coord_fixed() +
  scale_x_longitude(xmin=19, xmax=23, step=1, limits=c(18.81, 23.01)) +
  scale_y_latitude(ymin=42, ymax=46, step=1, limits=c(41.85, 46.2)) # +

### Fig3 ###
# tiff("plot/Fig3.tiff", width = 232, height = 180, units = 'mm', res = 600, compression = "lzw")
jpeg("plot/Fig3.jpeg", width = 232, height = 180, units = 'mm', res = 600)
ggarrange(ggarrange(tmax_plot, tmean_plot, tmin_plot, ncol=3, nrow=1, common.legend = T, legend="bottom"),
                # labels = c("1 Jan 2016", "2 Jan 2016", "3 Jan 2016", "4 Jan 2016")),
                ggarrange(slp_plot, prcp_plot, ncol=2, nrow=1, common.legend = F, legend="bottom"),
                ncol=1, nrow=2)
dev.off()

### CV bubble plots ######################################################################

for (var in c("tmax", "tmin", "tmean", "slp", "prcp")) {
  print(var)
  load(paste("ogimet/cv_", var, ".rda", sep=""))
  stNum = length(unique(cv.model$staid))
  results = matrix(nrow = stNum, ncol = 4)

  for (st in 1:stNum){
    staid = unique(cv.model$staid)[st]
    lon = cv.model[cv.model$staid==staid, ]$lon[1]
    lat = cv.model[cv.model$staid==staid, ]$lat[1]
    rmse = sqrt(mean((cv.model[cv.model$staid==staid, ]$obs - cv.model[cv.model$staid==staid, ]$pred)^2, na.rm = T))
    results[st,] = c(staid, lon, lat, rmse)
  }
  results <- as.data.frame(results)
  names(results) <- c('staid', 'lon', 'lat', 'rmse')

  load(paste("ogimet/ogimet_serbia08_", var, ".rda", sep=""))
  stations <- unique(ogimet_serbia[, 1:5])
  
  results <- join(results, stations, type = "left", match="first")
  results_order = results[order(results$rmse, decreasing = T), ]

  head(results_order, 20)
  tail(results_order, 20)

  if (var == "tmax") {
    title = "Tmax"
    yax <- expression(paste("Tmax [",degree,"C]", sep = ""))
  } else if (var == "tmin") {
    title = "Tmin"
    yax <- expression(paste("Tmin [",degree,"C]", sep = ""))
  } else if (var == "tmean") {
    title = "Tmean"
    yax <- expression(paste("Tmean [",degree,"C]", sep = ""))
  } else if (var == "slp") {
    title = "SLP"
    yax <- "SLP [mbar]"
  } else if (var == "prcp") {
    title = "PRCP"
    yax <- "PRCP [mm]"
  }
  
  coordinates(results_order) <-~ lon+lat
  assign(paste(var, "_bubble", sep=""),
         ggplot(as.data.frame(results_order), aes(x = lon, y = lat, size = rmse)) +
           geom_polygon(data = borders, aes(x = long, y = lat, group = group), alpha = 0.8, color = "black", fill="white", size = 0.1) +
           geom_point(shape = 21, colour = "mediumvioletred", fill = "springgreen", alpha = 0.8) +
           theme(plot.title = element_text(hjust = 0.5, face="bold"),
                 axis.text = element_text(size = 8),
                 axis.title = element_text(size = 8),
                 text = element_text(size = 8),
                 legend.position = "bottom",
                 legend.direction = "horizontal",
                 legend.key.size= unit(0.2, "cm"),
                 legend.margin = unit(0, "cm"),
                 legend.title = element_text(size=8, face="bold")) +
           labs(x = "Longitude", y = "Latitude") + labs(size = "RMSE") +
           ggtitle(title) +
           coord_fixed() +
           scale_size_identity(guide="legend", breaks = c(1, 2, 3, 4, 5, 6, 7, 8), limits = c(1, 8)) +
           scale_x_longitude(xmin=18, xmax=23, step=1, limits=c(18, 23.6)) +
           scale_y_latitude(ymin=41, ymax=47, step=1, limits=c(40.8, 47.2))
  )
  
  if (var=="tmin"){
    tmin_bubble <- tmin_bubble +
      scalebar(x.min = 18, x.max = 23.6,
               y.min = 40.8, y.max = 47.2,
               st.size = 3, location="bottomleft", st.dist=0.03, border.size=0.3,
               dist = 75, dist_unit = "km",
               transform = T, model = "WGS84",
               anchor=c(x=18.3, y=41.5)) +
      north(x.min = 18, x.max = 23.6,
            y.min = 40.8, y.max = 47.2,
            scale = 0.22, symbol = 3, location="bottomleft",
            anchor=c(x=22, y=45.6))
  }
  
  # time series for year 2014 and station Belgrade
  # stations[stations$name=="Beograd ",] # 13274
  Sys.setlocale("LC_TIME", "C")
  ts <- cv.model[cv.model$staid==13274 & substr(cv.model$date, 1, 4)=="2014", ]
  assign(paste(var, "_ts", sep=""),
         ggplot() +
           geom_line(data=ts, aes(x = date, y = obs, color = "black")) +
           geom_line(data=ts, aes(x = date, y = pred, color = "red")) +
           theme(plot.title = element_text(hjust = 0.5, face="bold"),
                 axis.text = element_text(size = 8),
                 axis.title = element_text(size = 8),
                 text = element_text(size = 8),
                 legend.position = "bottom",
                 legend.direction = "horizontal",
                 legend.key.size= unit(0.2, "cm"),
                 legend.margin = unit(0, "cm"),
                 legend.title = element_text(size=8, face="bold")) +
           labs(x = "", y = yax) +
           scale_colour_manual(name = "",
                               labels = c("Observation", "LLOCV prediction"),
                               values = c("black"="black", "red"="red"))
  )
  
  if (var=="prcp") {
    prcp_ts <- prcp_ts + labs(x = "Date")
  }
  
}

### Fig4 ###
# tiff("plot/Fig4.tiff", width = 232, height = 174, units = 'mm', res = 600, compression = "lzw")
jpeg("plot/Fig4.jpeg", width = 232, height = 174, units = 'mm', res = 600)

legend <- g_legend(prcp_bubble + theme(legend.position='bottom'))
grid.arrange(tmax_bubble+theme(legend.position='hidden'), tmean_bubble+theme(legend.position='hidden'), tmin_bubble+theme(legend.position='hidden'), 
             slp_bubble+theme(legend.position='hidden'), prcp_bubble+theme(legend.position='hidden'), legend,
             ncol=3, nrow=2)
dev.off()

### Fig5 ###
# tiff("plot/Fig5.tiff", width = 232, height = 250, units = 'mm', res = 600, compression = "lzw")
jpeg("plot/Fig5.jpeg", width = 232, height = 250, units = 'mm', res = 600)
ggarrange(tmax_ts, tmean_ts, tmin_ts, slp_ts, prcp_ts,
          ncol=1, nrow=5, common.legend = T, legend="bottom")
dev.off()

### E-OBS correlation ######################################################################

# # day <- "20140727"
# r <- raster("data/day_2014_wgs84/tmax_day_20140727_wgs84.tif")
# # eobs
# band = as.integer(date - as.Date("1950-01-01")) + 1
# eobs <- raster(paste("eobs/", var, "_ens_mean_0.1deg_reg_v21.0e.nc", sep=""), band=band)
# # eobs@z
# eobs <- crop(eobs, r)
# r <- mask(eobs, borders)
# writeRaster(r, "eobs_ref.tif", "GTiff", NAflag= -32767, datatype='INT2S', overwrite=T)
r <- raster("eobs_ref.tif")

for (var in c("tmax", "tmin", "tmean", "slp", "prcp")) {
  print(var)
  load(file = paste("ogimet/", var, "_eobs_values.rda", sep=""))
  load(file = paste("ogimet/", var, "_res_values.rda", sep=""))

  eobs <- as.vector(eobs.df)
  eobs <- eobs[!is.na(eobs)]
  srb <- as.vector(resample.df) # /10
  srb <- srb[!is.na(srb)]
  
  print(cor.test(eobs, srb))
  # 0.9923896
  ## RMSE
  print(sqrt(mean((eobs - srb)^2, na.rm = T)))
  # 1.253691
  ## R2
  print(1 - (t(eobs - srb) %*% (eobs - srb)) / (t(eobs - mean(eobs)) %*% (eobs - mean(eobs))))
  # 0.9845831
  
  corel <- c()
  for (pix in 1:dim(eobs.df)[2]) {
    a <- as.vector(eobs.df[, pix])
    b <- as.vector(resample.df[, pix]) # /10
    if (all(is.na(a))) {
      corel[pix] <- NA
    } else {
      corel[pix] <- cor.test(a, b)$estimate
    }
  }
  
  summary(corel)
  values(r) <- corel
  
  spdf <- as(r, "SpatialPixelsDataFrame")
  names(spdf) <-var
  
  if (var == "tmax") {
    title = "Tmax"
    brks = c(0.96,0.97,0.98, 0.99)
  } else if (var == "tmin") {
    title = "Tmin"
    brks = c(0.96,0.97,0.98, 0.99)
  } else if (var == "tmean") {
    title = "Tmean"
    brks = c(0.96,0.97,0.98, 0.99)
  } else if (var == "slp") {
    title = "SLP"
    brks = c(0.80,0.85,0.90, 0.95)
  } else if (var == "prcp") {
    title = "PRCP"
    brks = c(0.50,0.60,0.70)
  }
  
  theme = theme_set(theme_minimal())
  # plot(r)
  assign(paste(var, "_eobs", sep=""),
         ggplot() + # watch out for attribute name color order
           geom_raster(data=as.data.frame(spdf), aes_string(x="x", y="y", fill=var), alpha=0.8) +
           scale_fill_gradientn(colours = colorRampPalette(c("red", "yellow", "blue"))(30), name = "PCC", breaks=brks) +
           geom_polygon(data = borders, aes(x = long, y = lat, group = group), alpha = 0, color = "black", fill=NA, size = 0.1) +
           theme(plot.title = element_text(hjust = 8, face="bold"),
                 axis.text = element_text(size = 8),
                 axis.title = element_blank(),#element_text(size = 8),
                 text = element_text(size = 8),
                 legend.position = "bottom",
                 legend.direction = "horizontal",
                 # legend.key.size= unit(0.2, "cm"),
                 # legend.margin = unit(0, "cm"),
                 legend.title = element_text(size=8), #, face="bold"),
                 legend.text=element_text(size=8)) +
           theme(plot.title = element_text(hjust = 0.5)) +
           theme(axis.text = element_text(size = 8)) + # changes axis labels
           theme(axis.title = element_text(size = 8)) + # change axis titles
           # theme(text = element_text(size = 8))+
           ggtitle(title) +
           labs(x = "Longitude", y = "Latitude") +
           coord_fixed() +
           scale_x_longitude(xmin=19, xmax=23, step=1, limits=c(18.81, 23.01)) +
           scale_y_latitude(ymin=42, ymax=46, step=1, limits=c(41.85, 46.2))
         )
  
  if (var=="tmin") {
    tmin_eobs <- tmin_eobs +
      scalebar(x.min = 19, x.max = 23,
             y.min = 42, y.max = 46,
             st.size = 3, location="bottomleft", st.dist=0.03, border.size=0.3,
             dist = 50, dist_unit = "km",
             transform = T, model = "WGS84",
             anchor=c(x=19, y=42.1)) +
      north(x.min = 19, x.max = 23,
            y.min = 42, y.max = 46,
            scale = 0.25, symbol = 3, location="topright",
            anchor=c(x=22.8, y=46))
  }
  
}

### Fig6 ###
# tiff("plot/Fig6.tiff", width = 232, height = 180, units = 'mm', res = 600, compression = "lzw")
jpeg("plot/Fig6.jpeg", width = 232, height = 180, units = 'mm', res = 600)
ggarrange(ggarrange(tmax_eobs, tmean_eobs, tmin_eobs, ncol=3, nrow=1, common.legend = T, legend="bottom"),
          ggarrange(slp_eobs, prcp_eobs, ncol=2, nrow=1, common.legend = F, legend="bottom"),
          ncol=1, nrow=2, common.legend = F, legend="bottom")
dev.off()

#########################