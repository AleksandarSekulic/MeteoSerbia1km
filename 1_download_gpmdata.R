library(R.utils)
library(doParallel)
#library(RCurl)

# year=2015

wd = '/media/geocomp/060C0FE30C0FCC9D/meteo_serbia/'
# wd = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)
dir.create("imerg")
setwd("imerg")
years <- 2000:2019

time=seq(as.Date(paste(years[1],"-06-01", sep="")), as.Date(paste(years[length(years)],"-12-31", sep="")), by="day")
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

version <- 'V06B' # c('06B', '05B', '04B')

missing_days <- c()

# ftp://arthurhou.pps.eosdis.nasa.gov/gpmdata/2000/06/01/gis/
# 3B-DAY-GIS.MS.MRG.3IMERG.20000601-S000000-E235959.4560.V06B.zip

registerDoParallel(cores=detectCores())
foreach(i = 1:daysNum, .packages = c('R.utils')) %dopar% {
# for(i in c(1:daysNum)) {
    #Sys.sleep(5)
  n_30 = sprintf("%04d", as.numeric(strftime(time[i], format = "%j")) * 30 - 30)
  url = paste('ftp://aleksandarsale.sekulic%40gmail.com:aleksandarsale.sekulic%40gmail.com@arthurhou.pps.eosdis.nasa.gov/gpmallversions/',
              substr(version, 1,3), '/', substr(days[i], 1, 4), '/', substr(days[i], 5, 6), '/', substr(days[i], 7, 8), '/gis/',
              '3B-DAY-GIS.MS.MRG.3IMERG.', days[i], '-S000000-E235959.', n_30, '.', version,'.zip', sep = '')
  # url = paste('ftp://aleksandarsale.sekulic%40gmail.com:aleksandarsale.sekulic%40gmail.com@jsimpson.pps.eosdis.nasa.gov/data/imerg/gis/'
  #             ,substr(days[i], 1, 4), '/', substr(days[i], 5, 6), '/', 
  #             '3B-HHR-L.MS.MRG.3IMERG.', days[i], '-S233000-E235959.1410.V06A.1day.zip', sep = '')
  destfile = paste(days[i], '.zip', sep = '')
  d <- try(download.file(url, destfile))
  # v=1
  if(inherits(d, "try-error")) {
    cat(time[i])
    missing_days <- c(missing_days, time[i])
    # for (v in 1:length(version)) {
    #   url = paste('ftp://aleksandarsale.sekulic%40gmail.com:aleksandarsale.sekulic%40gmail.com@jsimpson.pps.eosdis.nasa.gov/data/imerg/gis/'
    #               ,substr(days[i], 1, 4), '/', substr(days[i], 5, 6), '/', 
    #               # V04/
    #               '3B-HHR-L.MS.MRG.3IMERG.', days[i], '-S233000-E235959.1410.V', version[v], '.1day.zip', sep = '')
    #   d <- try(download.file(url, destfile))
    #   if(inherits(d, "try-error")) {
    #     next
    #   } else {
    #     break
    #   }
    # }
  }
  unzip(destfile)
  unlink(destfile)
}
stopImplicitCluster()

files_remove <- list.files(pattern = "liquidPercent")
unlink(files_remove)
files_remove <- list.files(pattern = "liquid.accum")
unlink(files_remove)
files_remove <- list.files(pattern = "ice.accum")
unlink(files_remove)
files_remove <- list.files(pattern = "total.rate")
unlink(files_remove)
files_remove <- list.files(pattern = "liquid.rate")
unlink(files_remove)
files_remove <- list.files(pattern = "ice.rate")
unlink(files_remove)
