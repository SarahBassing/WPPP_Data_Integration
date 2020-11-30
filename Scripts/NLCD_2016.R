##################################
# Resample NLCD 2016 to Study Area Grid
##################################

# This script resamples NLCD 2016 to a pre-specified grid

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

# load NLCD 2016 layer and grids
NLCD <- raster("NLCD2016_NEp1.tif")
grid_NE <- raster("NE_grid_4k.img") #928 cells
NE <- shapefile("NE_covariate_area_2855.shp")

#check projection
crs(NLCD)
crs(grid_NE)
crs(NE)

# project the WPPP grid to the crs of the NLCD layer and clip
NLCD <- crop(NLCD, grid_WPPP)

# project the cropped NLCD raster to 2855
NLCD <- projectRaster(NLCD, crs = "+init=epsg:2855")

#check projection
crs(NLCD)
crs(grid_NE)
crs(NE)

# investigate the NLCD layer
plot(NLCD)
plot(grid_NE, add=T)
nlcd_val <- getValues(NLCD)
summary(nlcd_val)

#reclassify the NLCD raster into fewer groups
#create a vector of values, each set of 3 in reclass is a "from","to","reclassify as"
reclass_df <- c(10,12,1, #water
                20,22,2, #low intensity development
                22,24,3, #high intensity development
                30,31,4, #barren
                40,43,5, #forest
                50,52,6, #shrubland
                70,74,7, #grassland
                80,82,8, #fields/crops
                89,95,9) #wetlands

#create into a matrix
reclass_m <- matrix(reclass_df, ncol=3, byrow=TRUE)

#reclassify the nlcd raster
NLCD <- reclassify(NLCD, reclass_m)
summary(getValues(NLCD))

# crop the larger DEM to the NE study area boundary
NLCD_NE <- crop(NLCD, grid_NE)
plot(NLCD_NE)
plot(NE, add=T)

# resample the WA DEM from 30 x 30 m to 4 x 4 km
NLCD_4k_NE <- resample(NLCD, grid_NE, method="ngb") #use nearest neighbor resample for categorical raster
plot(NLCD_4k_NE)
plot(NE, add=T)
NLCD_4k_NE #928 cells
summary(NLCD_4k_NE)

# write the raster file
writeRaster(NLCD_4k_NE,"NE_NLCD_4k.tif")
