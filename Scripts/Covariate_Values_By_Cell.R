###########################################
# THIS SCRIPT EXTRACTS COVARIATE VALUES BY 
# CELL AND CREATES A COVARIATE DATA TABLE
##########################################

# clear workspace
rm(list = ls())

# automatically set the working directory as the location of the .r file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# check the working directory
getwd()

# load libraries
library(raster)

# set projection
crs <- CRS("+init=epsg:2855")

# read in the final covariate rasters and NE grid
DEM <- raster("NE_DEM_4k.tif")
roads <- raster("NE_RoadDensity_4k.tif")
NLCD <- raster("NE_NLCD_4k.tif")
grid_NE <- raster("NE_grid_4k.img")
NE <- shapefile("NE_covariate_area_2855.shp")


#check projection
crs(DEM)
crs(roads)
crs(NLCD)
crs(grid_NE)
crs(NE)


# extract values by grid cell
grid_cell_NE <- getValues(grid_NE) 
DEM_val <- getValues(DEM)
roads_val <- getValues(roads)
NLCD_val <- getValues(NLCD)
index <- as.vector(as.numeric(seq_along(DEM_val)))
table(NLCD_val)
levels(NLCD_val)

# pair the NLCD numeric values with their categorical values

# first a toy example
# testvec <- c(1,1,1,2,3,4,5,5,6,7,6,7,9, NA)
# testvec
# testout <- factor(testvec, levels = c(1:9), labels = c("water","lowdevelop","highdevelop","barren","forest","shrubland","grassland","fields_crops","wetlands"))
# testout

NLCD_label <- factor(NLCD_val, level = c(1:9), labels = c("water","lowdevelop","highdevelop","barren","forest","shrubland","grassland","fields_crops","wetlands"))

# combine into a dataframe
covars <- data.frame(index, grid_cell_NE, DEM_val, roads_val, NLCD_val, NLCD_label, row.names = T)
head(covars)
View(covars)

# write .csv file
write.csv(covars, "Covariates_by_cell.csv", row.names = FALSE)


