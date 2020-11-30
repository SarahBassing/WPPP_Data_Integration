##################################
# Resample DEM to Study Area Grid
##################################

# This script resamples elevation to a pre-specified grid

# Some sources:
# https://stackoverflow.com/questions/32278825/how-to-change-the-resolution-of-a-raster-layer-in-r
# https://rdrr.io/cran/raster/man/resample.html

# clear workspace
rm(list = ls())

# automatically set the working directory as the location of the .r file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# check the workding directory
getwd()

# load libraries
library(raster)
library(rgdal)
library(rgeos)

# load DEM layer and grid
DEM <- raster("Washington_DEM.tif")
grid_NE <- raster("NE_grid_4k.img") #928 cells
NE <- shapefile("NE_covariate_area_2855.shp")

#check projection
crs(DEM)
crs(grid_NE)
crs(NE)

# project needed files to EPSG 2855 and check
DEM <- projectRaster(DEM, crs = "+init=epsg:2855")
crs(DEM)
DEM

# crop the larger DEM to the NE study area boundary
DEM_NE <- crop(DEM, NE)
plot(DEM_NE)
plot(NE, add=T)

# resample the WA DEM from 30 x 30 m to 4 x 4 km
DEM_4k_NE <- resample(DEM, grid_NE, method="bilinear")
plot(DEM_4k_NE)
plot(NE, add=T)
DEM_4k_NE #928 cells
summary(DEM_4k_NE)

# write the raster file
writeRaster(DEM_4k_NE,"NE_DEM_4k.tif")
