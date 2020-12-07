#################################
# Merge DEM raster tiles into one raster
#################################
# Data source: https://dwtkns.com/srtm30m/ (need free NASA account first)
# Code source: https://stackoverflow.com/questions/15876591/merging-multiple-rasters-in-r; https://stackoverflow.com/questions/47471117/create-a-list-of-all-loaded-rasters-in-a-for-loop-in-r
# For the WPPP, I downloaded the following tiles, 17 in total:
#   N49W122 to N49W117
#   N48W122 to N48W117
#   N47W122 to N47W117

# clear workspace
rm(list = ls())

# automatically set the working directory as the location of the .r file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# check the working directory
getwd()

# load libraries
library(raster)
library(gdata)
library(rgdal)

# load the tiles
tiles <- list.files(pattern = "\\.hgt$")
tiles #18 files total

# read in all the .hgt raster files and name them by their tile name (same as file name)
for (i in 1:(length(tiles))) {
  temp.rast <- raster(tiles[i])
  temp.name <- NULL
  temp.name[i] <- temp.rast@data@names
  mv(from = "temp.rast", to = temp.name[i])
  print(temp.name[i])
}

# check by plotting one
plot(N49W120)

# compile the files into a list
tiles.list <- lapply(tiles, raster)

# mosaic into one raster
WPPP_DEM <- do.call(merge, tiles.list)
plot(WPPP_DEM)

# write the combined raster
writeRaster(WPPP_DEM, file="WPPP_DEM_30m.tif", format="GTiff")
