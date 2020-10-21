library(httr)
library(RJSONIO)

wd=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

library(rgdal)
borders <- readOGR("borders/osm_nominatim/Polygon.shp")

users <- c("bt", "vs", "vš", "ns", "zr", "ki", "pa", "ru", "sa", "sm", "so", "su")
stations <- matrix(nrow=0, ncol=8)

for (user in users) {
  body = paste('{"grant_type": "password",
	"client_id": "FieldclimateNG",
	"client_secret": "618a5baf48287eecbdfc754e9c933a",
	"username": "', user, '",
	"password": "pis"}', sep="")
  login <- POST("https://oauth.fieldclimate.com/token", content_type_json(), body = body)
  token <- content(login)$access_token
  
  station_request <- GET("https://api.fieldclimate.com/v2/user/stations",
                  add_headers("Authorization" = paste("Bearer", token, sep = " ")))
                  # authenticate(user, "pis"))
  # content(station_request, "text")
  # headers(stations)
  station_request <- fromJSON(content(station_request, "text"))
  for(st in 1:length(station_request)) {
    st.obs <- station_request[st][[1]]
    # c("id", "name", "lon", "lat", "from", "altitude", "from", "to", "user")
    stations <- rbind(stations, c(st.obs$name[1],
                                  st.obs$name[2],
                                  st.obs$position$geo$coordinates[1],
                                  st.obs$position$geo$coordinates[2],
                                  st.obs$position$altitude,
                                  st.obs$dates[1],
                                  st.obs$dates[2],
                                  user
                                  ))
  }
}

stations <- as.data.frame(stations)
names(stations) <- c("id", "name", "lon", "lat", "altitude", "from", "to", "user")
stations$id <- as.character(stations$id)
stations$name <- as.character(stations$name)
stations$lon <- as.numeric(as.character(stations$lon))
stations$lat <- as.numeric(as.character(stations$lat))
stations$altitude <- as.numeric(as.character(stations$altitude))
stations$from <- as.Date(as.character(stations$from))+1
stations$to <- as.Date(as.character(stations$to))-1
stations$user <- as.character(stations$user)
summary(stations)
# save(stations, file="ogimet/asv_stations.rda")

load(file="ogimet/asv_stations.rda")

aaa <- stations
coordinates(aaa) <-~ lon+lat
plot(aaa)
plot(borders, add=T)


asv.df <- matrix(nrow=0, ncol=6)
for (st in 1:nrow(stations)) {
  print(st)
  station <- stations[st, ]
  print(station$id)
  print("###########")
  
  body = paste('{"grant_type": "password",
	"client_id": "FieldclimateNG",
	"client_secret": "618a5baf48287eecbdfc754e9c933a",
	"username": "', station$user, '",
	"password": "pis"}', sep="")
  login <- POST("https://oauth.fieldclimate.com/token", content_type_json(), body = body)
  token <- content(login)$access_token
  
  station_data <- GET(paste("https://api.fieldclimate.com/v2/fc/",
                            station$id, "/daily/from/946684800/to/1577750400", sep=""),
                      add_headers("Authorization" = paste("Bearer", token, sep = " ")))
  # authenticate(user, "pis"))
  content(station_data, "text")
  # headers(stations)
  station_data <- fromJSON(content(station_data, "text"))
  length(station_data[[3]]$data)
  station.df <- lapply(station_data[[3]]$data, function(x) unlist(x)[c("datetime",
                                                                       "sensor_x_x_22_506_a",
                                                                       "sensor_x_x_22_506_mx",
                                                                       "sensor_x_x_22_506_mn",
                                                                       "sensor_x_x_5_6_s")])
  station.df <- do.call("rbind", station.df)
  # head(station.df)
  
  # c("id", "date", "tmean", "tmax", "tmin", "prcp")
  asv.df <- rbind(asv.df, cbind(station$id, station.df))
  # head(asv.df)
}

asv.df <- as.data.frame(asv.df)
names(asv.df) <- c("id", "date", "tmean", "tmax", "tmin", "prcp")
asv.df$id  <- as.character(asv.df$id)
asv.df$date  <- as.Date(as.character(asv.df$date))
asv.df$tmean  <- as.numeric(as.character(asv.df$tmean))
asv.df$tmax  <- as.numeric(as.character(asv.df$tmax))
asv.df$tmin  <- as.numeric(as.character(asv.df$tmin))
asv.df$prcp  <- as.numeric(as.character(asv.df$prcp))
summary(asv.df)
# save(asv.df, file="ogimet/asv_obs.rda")

load(file="ogimet/asv_obs.rda") # 182858
asv.df <- asv.df[complete.cases(asv.df),] # 171688
asv.df <- asv.df[asv.df$tmax > asv.df$tmean, ] # 171634
asv.df <- asv.df[asv.df$tmax > asv.df$tmin, ] # 171634
asv.df <- asv.df[asv.df$tmean > asv.df$tmin, ] # 171601

summary(asv.df)

# save(asv.df, file="ogimet/asv_obs.rda")

### Outliers detection ##############################################

library(plyr)
library(sp)
load(file="ogimet/asv_obs.rda")
load(file="ogimet/asv_stations.rda")
asv.df <- join(asv.df, stations)
summary(asv.df)

asv.df_test <- asv.df[asv.df$prcp>40 | asv.df$tmax>30 | asv.df$tmin<(-15), ]
nrow(asv.df_test)

# remove prcp > 200 mm - impossible (max 211.1 mm in 10 October 1955 in Negotin)
prcp.treshold = 4
temp.treshold = 10

br = 0
out <- c()

wgs84 <- CRS("+proj=longlat +datum=WGS84")

for(o in 1:nrow(asv.df_test)) {

  om <- asv.df_test[o, ]
  
  coordinates(om) <- c("lon", "lat")
  om@proj4string <- wgs84
  asv_date <- asv.df[asv.df$date==om$date, ] # c("staid","date",var,"eobs", "imerg", "lon", "lat")]
  coordinates(asv_date) <- c("lon", "lat")
  asv_date@proj4string <- wgs84
  asv_date$dist <- spDists(asv_date, om, longlat = T)[, 1]
  asv_date <- asv_date[order(asv_date$dist), ]
  
  asv_date <- asv_date[asv_date$dist>1&asv_date$dist<50, ]
  if (nrow(asv_date)==0) {
    next
  }
  
  for (var in c("tmax", "tmin", "tmean", "prcp")) {
    
    if (var=="prcp"){
      if (max(asv_date$prcp) == 0) {
        cond <- (om$prcp-max(asv_date$prcp)) > temp.treshold
      } else {
        cond <- (om$prcp/max(asv_date$prcp)) > prcp.treshold
      }
      if (cond) {
        print(paste(om$id, om$date, om$prcp, sep="  "))
        print(asv_date$prcp)
        print("")
        br = br+1
        out <- cbind(out, c(om$staid, om$date))
        asv.df <- asv.df[!(asv.df$id==om$id & asv.df$date==om$date), ]
        break
      }
    } else if (var=="tmax" | var=="tmin" | var=="tmean") {
      if (abs(om@data[, var] - mean(asv_date@data[, var], na.rm=T)) > temp.treshold) {
        print(paste(om@data$id, om@data$date, om@data[, var], sep="  "))
        print(asv_date@data[, var])
        print("")
        br = br+1
        out <- cbind(out, c(om@data$id, om@data$date))
        asv.df <- asv.df[!(asv.df$id==om@data$id & asv.df$date==om@data$date), ]
        break
      }
    }
  }
}

br # 634 out

#############

asv.df <- asv.df[asv.df$tmax < 45, ] # 171504
asv.df <- asv.df[asv.df$tmin > -40, ] # 171504
asv.df <- asv.df[asv.df$prcp < 211.1, ] # 171462
summary(asv.df)

asv.df[which.max(asv.df$tmax), ]
asv.df[asv.df$date==asv.df[which.max(asv.df$tmax), "date"], ]
stations[stations$id==asv.df[which.max(asv.df$tmax), "id"], ]

save(asv.df, file="ogimet/asv_df.rda")

# # https://oauth.fieldclimate.com/token
# # https://api.fieldclimate.com/v2/user/stations
# # do 2019-12-31 -> 1577750400
# as.Date("2019-12-31") - as.Date("2005-03-08") # 5411
# 1577750400-1110240000 # 467510400
# 467510400/5411 # 86400 # 60*60*24
# 1577750400-(as.Date("2019-12-31")-as.Date("2005-03-08"))*60*60*24
# # 2000-01-01 -> 946684800
# as.Date("2019-12-31") - as.Date("2000-01-01") # 7304
# 1577750400-7304*60*60*24
# # https://api.fieldclimate.com/v2/view/fc/00000A1E/from/1110240000/to/1577750400
# # https://api.fieldclimate.com/v2/view/fc/00000A1E/from/946684800/to/1577750400