# MeteoSerbia1km

This repository contains scripts and data used in the manuscript "A high-resolution daily gridded meteorological dataset for Serbia made by Random Forest Spatial Interpolation" by Sekulić A., Kilibarda M., Protić D. and Bajat B., that was submitted to the Scientific Data journal in October 2020.

MeteoSerbia1km is the first daily gridded meteorological dataset at a 1-km spatial resolution across Serbia for the 2000–2019 period. The dataset consists of five daily variables: maximum, minimum and mean temperature, mean sea level pressure, and total precipitation. Besides daily summaries, it contains monthly and annual summaries, daily, monthly, and annual long term means (LTM).

## Data download

Data can be downloaded from [ZENODO](https://doi.org/10.5281/zenodo.4058167)

## Methodology

Daily gridded data were interpolated using the Random Forest Spatial Interpolation ([RFSI](https://github.com/AleksandarSekulic/RFSI)) (Sekulić et al. 2020) methodology based on Random Forest (Breiman 2001) and using nearest observations and distances to them as spatial covariates, together with environmental covariates.

## Scripts and data description

The scripts are:
- [1_download_gpmdata.R](1_download_gpmdata.R) - a script for downloading the IMERG data (Huffman et al. 2014),
- [2_download_ogimet.R](2_download_ogimet.R) - a script for downloading the [OGIMET](https://www.ogimet.com/) data,
- [3_asv_download.R](3_asv_download.R) - a script for downloading the automated meteorological stations network in Vojvodina ([ASV](http://www.pisvojvodina.com/Shared%20Documents/AMS_pristup.aspx)) data,
- [4_final_script.R](4_final_script.R) - a script for data preparation, tuning, fitting, cross-validation, independent validation, and prediction for RFSI models of all daily meterological variables for Serbia, and making plots.

*Note that the number in the script name refers to the order in which the scripts should be run.*

The data are:
- [borders](borders/) - a folder containing Serbian border (ESRI Shapefile),
- [dem_twi](dem_twi/) - a folder containing digital elevation model (DEM) and topographic wetness index (TWI) for Serbia (*_buff files are DEM and TWI for 200 km - buffer area from Serbian border) (GeoTIFF),
- [eobs_ref.tif](eobs_ref.tif) - E-OBS (Cornes et al. 2018) reference raster cropped by Serbian border,
- [ogimet](ogimet/) - a folder containing all of the data (Rdata). The data are:
    - asv_*.rda - data used for validation with ASV,
    - cv_*.rda - nested 5-fold leave-location-out cross-validation results,
    - stations_serbia.rda - OGIMET stations,
    - ogimet_serbia*.rda - OGIMET daily meteorlogical observations,
    - *_eobs_values.rda - E-OBS raster values,
    - *_res_values.rda - MeteoSerbia1km resampled raster to E-OBS values,
    - tuned1_*.rda - mtry and n.obs tuning results,
    - tuned2_*.rda - sample.fraction and min.node.size tuning results,
    - tuned3_*.rda - idw.p tuning results,

## References

* Breiman, L. (2001). Random Forests. Machine Learning, 45(1), 5–32. https://doi.org/10.1023/A:1010933404324
* Cornes, R. C., van der Schrier, G., van den Besselaar, E. J. M., & Jones, P. D. (2018). An Ensemble Version of the E-OBS Temperature and Precipitation Data Sets. Journal of Geophysical Research: Atmospheres, 123(17), 9391–9409. https://doi.org/10.1029/2017JD028200
* Huffman, G. J., Bolvin, D. T. & Nelkin, E. J. Integrated Multi-satellitE Retrievals for GPM (IMERG), Final Run, versionV06B. ftp://arthurhou.pps.eosdis.nasa.gov/gpmdata/ (2014). Accessed: 31 July, 2019.
* Sekulić, A., Kilibarda, M., Heuvelink, G. B., Nikolić, M. & Bajat, B. Random Forest Spatial Interpolation.Remote. Sens. 12, 1687, https://doi.org/10.3390/rs12101687 (2020b).