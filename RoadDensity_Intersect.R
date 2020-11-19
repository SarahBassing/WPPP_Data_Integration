#################################################################################
# ROAD DENSITY FROM GRID - This script calculates road density within a preset raster grid #
#################################################################################

# Source code adapted from: https://gis.stackexchange.com/questions/119993/convert-line-shapefile-to-raster-value-total-length-of-lines-within-cell

# clear workspace
rm(list = ls())

# automatically set the working directory as the location of the .r file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# check the workding directory
getwd()

library(raster)
library(rgdal)
library(rgeos)

# First a toy example with roads data from Tanzania
# Link to Tanzania roads file: http://www.diva-gis.org/gdata

roads <- shapefile("TZA_roads.shp")
roads <- spTransform(roads, CRS("+proj=utm +zone=37 +south +datum=WGS84"))
rs <- raster(extent(roads), crs=projection(roads))
rs[] <- 1:ncell(rs)

# Intersect lines with raster "polygons" and add length to new lines segments
rsp <- rasterToPolygons(rs)
rp <- intersect(roads, rsp)
rp$length <- gLength(rp, byid=TRUE) / 1000 # to convert m to km
x <- tapply(rp$length, rp$layer, sum)
r <- raster(rs)
r[as.integer(names(x))] <- x

# Visualize Results
plot(r)
plot(roads, add = T)

############################################################################

# Now calculate road density with our WPPP study area data on a preset grid 

#load shapefiles and grid raster
roads <- shapefile("WPPProads_Cascades.shp")
NE <- shapefile("NE_covariate_area.shp")
WPPP <- shapefile("WPPP_CovariateBoundary.shp")
grid <- raster("ref_grid_4k.img")
grid_NE <- raster("NE_grid_4k.img")

#check projection
crs(roads)
crs(NE)
crs(WPPP)
crs(grid)
crs(grid_NE)

# project needed files to EPSG 2855 and check
roads <- spTransform(roads, CRS("+init=epsg:2855"))
NE <- spTransform(NE, CRS("+init=epsg:2855"))
WPPP <- spTransform(WPPP, CRS("+init=epsg:2855"))

# visualize the data
plot(grid)
plot(roads)
plot(WPPP, add=T)
plot(NE, add=T)


# crop full roads layer to grid boundary (since they are off otherwise)
roads_NE <- crop(roads, grid_NE)
plot(grid_NE)
plot(roads_NE, add=T)
plot(NE, add=T)

# create empty cells within the grid
grid_NE[] <- 1:ncell(grid_NE)
grid_NE
plot(grid_NE)
plot(roads_NE, add=T)

# Intersect lines with raster "polygons" and add length to new lines segments
grid_NEp <- rasterToPolygons(grid_NE)
rp <- intersect(roads_NE, grid_NEp)
summary(rp)
rp$length <- gLength(rp, byid=TRUE) / 1000 # to convert m to km
x <- tapply(rp$length, rp$NE_grid_4k, sum)
r <- raster(grid_NE)
r[as.integer(names(x))] <- x

# Visualize Results
plot(r)
plot(roads_NE, add = T)
writeRaster(r,"NE_RoadDensity_4k.tif")

####################################################
# now copy the steps above for the entire study area
# create empty cells within the grid
grid[] <- 1:ncell(grid)
grid
plot(grid)
plot(roads, add=T)

# Intersect lines with raster "polygons" and add length to new lines segments
gridp <- rasterToPolygons(grid)
rp_all <- intersect(roads, gridp)
summary(rp_all)
rp_all$length <- gLength(rp_all, byid=TRUE) / 1000 # to convert m to km
x_all <- tapply(rp_all$length, rp_all$ref_grid_4k, sum)
r_all <- raster(grid)
r_all[as.integer(names(x_all))] <- x_all

# Visualize Results
plot(r_all)
plot(roads, add = T)
writeRaster(r_all,"WPPP_RoadDensity_4k.tif")
