# library(remotes)
# install_github("bczernecki/climate")

library(sp)
library(spacetime)
library(raster)
library(rgdal)
library(climate)
library(plyr)

wd=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

borders <- readOGR("borders/osm_nominatim/Polygon.shp")
border.buffer <- buffer(borders, 2, dissolve=T) # 1.5 degree -> 150 km

# # dem and twi #
# dem = crop(raster('../meteo/dem_twi/dem.sdat'), borders)
# dem = mask(dem, borders)
# writeRaster(dem, '../dem_twi/dem.tif', "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
# 
# twi = crop(raster('../meteo/dem_twi/twi.sdat'), borders)
# twi = mask(twi, borders)
# writeRaster(twi, '../dem_twi/twi.tif', "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
# 
# dem = crop(raster('../meteo/dem_twi/dem.sdat'), border.buffer)
# # dem = mask(dem, border.buffer)
# writeRaster(dem, 'dem_twi/dem_buff.tif', "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
# 
# twi = crop(raster('../meteo/dem_twi/twi.sdat'), border.buffer)
# # twi = mask(twi, border.buffer)
# writeRaster(twi, 'dem_twi/twi_buff.tif', "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)


years <- 2000:2019

countries = c("Serbia", "Hungary", "Romania", "Bulgaria", "Macedonia", "Albania",
              "Montenegro", "Croatia", "Bosnia", "Slovenia", "Austria", "Slovakia") # , "Kosovo"
stations_serbia <- c()
for (country in countries) {
  stations_serbia <- rbind(stations_serbia, nearest_stations_ogimet(country = country,
                                                                    date = Sys.Date(),
                                                                    add_map = F,
                                                                    point = c(20.806635, 44.082873),
                                                                    no_of_stations = 200)
  )
}

serbia = point.in.polygon(stations_serbia$lon, stations_serbia$lat,
                          border.buffer@polygons[[1]]@Polygons[[1]]@coords[,1],
                          border.buffer@polygons[[1]]@Polygons[[1]]@coords[,2])
# c(18,23.8,23.8,18),
# c(41.4,41.4,47,47))
stations_serbia = stations_serbia[serbia != 0, ] # 93 stations
stations_serbia$wmo_id <- as.numeric(as.character(stations_serbia$wmo_id))
stations_serbia$station_names <- as.character(stations_serbia$station_names)
# save(stations_serbia, file = "ogimet/stations_serbia.rda")
load(file = "ogimet/stations_serbia.rda")

st <- stations_serbia
coordinates(st) <-~ lon + lat
# plot(stfdf@sp)
plot(st) #, pch=0, add=T)
plot(borders, add=T)
plot(border.buffer, add=T)

# ### by year #################
# ogimet_serbia <- c()
# wmo_year_missing <- as.data.frame(matrix(nrow = 0, ncol=2))
# for (wmo in stations_serbia$wmo_id) {
#   for (year in years){
#     wmo_year <- try(ogimet_daily(date = c(paste(year,"-01-02", sep=""), paste(year+1,"-01-01", sep="")),
#                                  coords = FALSE,
#                                  station = wmo,# c(12326, 12330),
#                                  hour = 0))
#     if(inherits(wmo_year, "try-error")) {
#       wmo_year_missing <- rbind(wmo_year_missing, c(wmo, year))
#     } else {
#       ogimet_serbia <- rbind(ogimet_serbia, wmo_year)
#     }
#   }
# }
# # save(ogimet_serbia, file = "ogimet/ogimet_serbia.rda")
# load(file = "ogimet/ogimet_serbia.rda")
# 
# wmo_year_missing <- as.data.frame(matrix(nrow = 0, ncol=2))
# for (wmo in stations_serbia$wmo_id) {
#   df_years <- as.integer(unique(substr(ogimet_serbia[ogimet_serbia$station_ID==wmo, "Date"],1,4)))
#   dif_years <- setdiff(years, df_years)
#   if(length(dif_years) !=0) {
#     # wmo_year_missing <- rbind(wmo_year_missing, data.frame(wmo, dif_years))
#   }
# }
# # save(wmo_year_missing, file="ogimet/wmo_year_missing.rda")
# load(file="ogimet/wmo_year_missing.rda")
# 
# unique(wmo_year_missing$wmo)

# ### by month ############################

# station_ID
# Date 
# TemperatureCAvg - average air temperature at 2 metres above ground level. Values given in Celsius degrees
# TemperatureCMax - maximum air temperature at 2 metres above ground level. Values given in Celsius degrees
# TemperatureCMin - minimum air temperature at 2 metres above ground level. Values given in Celsius degrees
# TdAvgC - average dew point temperature at 2 metres above ground level. Values given in Celsius degrees
# HrAvg - average relative humidity. Values given in %
# WindkmhDir - wind direction
# WindkmhInt - wind speed in km/h
# WindkmhGust - wind gust in km/h
# PresslevHp - Sea level pressure in hPa
# Precmm - precipitation totals in mm
# TotClOct - total cloudiness in octants
# lowClOct - cloudiness by low level clouds in octants
# SunD1h - sunshine duration in hours
# PreselevHp - atmospheric pressure measured at altitude of station in hPa
# SnowDepcm - depth of snow cover in centimetres

ogimet_serbia_month <- as.data.frame(matrix(nrow = 0, ncol=18))
wmo_month_missing <- as.data.frame(matrix(nrow = 0, ncol=3))
for (wmo in stations_serbia$wmo_id[1:50]) {
  for (year in years){
    for (month in sprintf("%02d", 1:12)) {
      if (month != "12") {
        range.dates = c(paste(year,"-", month, "-02", sep=""), paste(year,"-", sprintf("%02d", as.integer(month)+1), "-01", sep=""))
      } else {
        range.dates = c(paste(year,"-", month, "-02", sep=""), paste(year+1,"-01-01", sep=""))
      }
      wmo_month <- try(ogimet_daily(date = range.dates,
                                    coords = FALSE,
                                    station = wmo,# c(12326, 12330),
                                    hour = 0))
      if(inherits(wmo_month, "try-error")) {
        wmo_month_missing <- rbind(wmo_month_missing, c(wmo, year, month))
        if (is.factor(wmo_month_missing[, 1])){
          wmo_month_missing[, 1] <- as.character(wmo_month_missing[, 1])
          wmo_month_missing[, 2] <- as.character(wmo_month_missing[, 2])
          wmo_month_missing[, 3] <- as.character(wmo_month_missing[, 3])
        }
      } else {
        if(length(ogimet_serbia_month) == length(wmo_month)) {
          ogimet_serbia_month <- rbind(ogimet_serbia_month, wmo_month)
        } else {
          print("NECE!!!")
        }
        
      }
    }
  }
}

# save(ogimet_serbia_month, file = "ogimet/ogimet_serbia_month.rda")
load(file = "ogimet/ogimet_serbia_month.rda")

ogimet_serbia <- ogimet_serbia_month
# save(ogimet_serbia, file = "ogimet/ogimet_serbia.rda")
load(file = "ogimet/ogimet_serbia.rda")

ogimet_serbia <- ogimet_serbia[!duplicated(ogimet_serbia), ]

summary(ogimet_serbia$Date)
ogimet_serbia$Date <- ogimet_serbia$Date - 1 # because of 00 hour - previous day
summary(ogimet_serbia$Date)

# save(ogimet_serbia, file = "ogimet/ogimet_serbia.rda")
load(file = "ogimet/ogimet_serbia.rda")

# # remove 12-31 and do it again
# nrow(ogimet_serbia[substr(ogimet_serbia$Date, 6, 10) == "12-31", ])
# summary(ogimet_serbia[substr(ogimet_serbia$Date, 6, 10) == "12-31", ])
# ogimet_serbia <- ogimet_serbia[substr(ogimet_serbia$Date, 6, 10) != "12-31", ]

### remove all na
# summary(is.na(ogimet_serbia$TemperatureCAvg))
summary(!rowSums(is.na(ogimet_serbia[, c(3:5,14,16)])) == ncol(ogimet_serbia[, c(3:5,14,16)]))
ogimet_serbia <- ogimet_serbia[!(rowSums(is.na(ogimet_serbia[, c(3:5,14,16)])) == ncol(ogimet_serbia[, c(3:5,14,16)])), ] # 533747

summary(duplicated(ogimet_serbia))
summary(duplicated(ogimet_serbia[, 1:2]))
duplicates <- ogimet_serbia[duplicated(ogimet_serbia[, 1:2]), ]
head(duplicates)
summary(duplicates)
ogimet_serbia <- ogimet_serbia[!duplicated(ogimet_serbia[, 1:2]), ] # 533595
# save(ogimet_serbia, file = "ogimet/ogimet_serbia.rda")
load(file = "ogimet/ogimet_serbia.rda")

###
# try by day
all.dates <- as.character(seq(as.Date("2000-01-01"), as.Date("2019-12-31"), by="day"))
wmo_dates_missing <- as.data.frame(matrix(nrow = 0, ncol=2))
for (wmo in stations_serbia$wmo_id) {
  df_dates <- unique(as.character(ogimet_serbia[ogimet_serbia$station_ID==wmo, "Date"]))
  dif_dates <- setdiff(all.dates, df_dates)
  if(length(dif_dates) !=0) {
    wmo_dates_missing <- rbind(wmo_dates_missing, data.frame(wmo, dif_dates))
    if (is.factor(wmo_dates_missing[, 1])){
      wmo_dates_missing[, 1] <- as.character(wmo_dates_missing[, 1])
      wmo_dates_missing[, 2] <- as.character(wmo_dates_missing[, 2])
      wmo_dates_missing[, 3] <- as.character(wmo_dates_missing[, 3])
    }
  }
}
# save(wmo_dates_missing, file="ogimet/wmo_dates_missing.rda")
load(file="ogimet/wmo_dates_missing.rda")

### missing days ###

ogimet_serbia_dates <- as.data.frame(matrix(nrow = 0, ncol=18))
for (wmo in unique(wmo_dates_missing$wmo)) {
  for (date in as.character(wmo_dates_missing[wmo_dates_missing$wmo==wmo, ]$dif_dates)){
    if (substr(date, 6, 10) == "12-31") {
      hour=23
      # range.dates = date
      add=0
    } else {
      hour=0
      # range.dates = as.Date(date)+1
      add=1
    }
    wmo_date <- try(ogimet_daily(date = as.Date(date)+add,
                                 coords = FALSE,
                                 station = wmo,# c(12326, 12330),
                                 hour = hour))
    if(inherits(wmo_date, "try-error")) {
      # wmo_month_missing <- rbind(wmo_month_missing, c(wmo, year, month))
    } else {
      if (add==1){
        wmo_date$Date <- wmo_date$Date-1
      }
      if(length(ogimet_serbia_dates) == length(wmo_date)) {
        ogimet_serbia_dates <- rbind(ogimet_serbia_dates, wmo_date)
      } else {
        print("It doesn't work!!!")
      }
    }
  }
}
# save(ogimet_serbia_dates, file = "ogimet/ogimet_serbia_dates.rda")
load(file = "ogimet/ogimet_serbia_dates.rda")
summary(duplicated(ogimet_serbia_dates[, 1:2]))

load(file = "ogimet/ogimet_serbia.rda")