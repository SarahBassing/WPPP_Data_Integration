  #'  WPPP Data Integration Project
  #'  SEFS 521 class project
  #'  
  #'  Sarah Bassing
  #'  November 2020
  #'  ------------------------------
  #'  This script generates a grid to guide data sampling process. 
  #'  ------------------------------


  #'  Load libraries
  library(sp)
  library(rgdal)
  library(rgeos)
  library(raster) 
  
  #'  Read spatial data
  OK_SA <- readOGR("./Shapefiles/fwdstudyareamaps", layer = "METHOW_SA") 
  NE_SA <- readOGR("./Shapefiles/fwdstudyareamaps", layer = "NE_SA") 
  NE_box <- readOGR("./Shapefiles/NE_covariate_area_2855", layer = "NE_covariate_area_2855")
  extra_SA <- readOGR("./Shapefiles/WPPP_CovariateBoundary", layer = "WPPP_CovariateBoundary")
  NE_box2 <- readOGR("./Shapefiles/NE_covariate_area", layer = "NE_covariate_area")
  
  #'  Get everything in same projected coordinate system (UTMs)
  #'  NE_box already in correct project in using NE_covariate_area_2855.shp
  sa_proj <- projection("+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs ")
  extra_SA <- spTransform(extra_SA, sa_proj)
  NE_SA <- spTransform(NE_SA, sa_proj)
  OK_SA <- spTransform(OK_SA, sa_proj)
  NE_box2 <- spTransform(NE_box2, sa_proj)
  
  plot(extra_SA)
  plot(NE_SA, add = T)
  plot(NE_box, add = T)
  plot(OK_SA, add = T)
  
  #'  Set bounding box for full extent, includes both study areas and then some
  print(bb_extra <- bbox(extra_SA))

  #'  Create reference raster to sample from
  #'  Set grid cell resolution (units = meters)
  #'  Smallest is 1 km x 1 km ==> 1000 m x 1000 m ==> 1e+6 m^2
  #ref <- raster(extent(bb_extra), crs = projection(sa_proj), res = 1000)
  #'  For computational reasons, starting with a larger grid cell (4 km^2)
  ref <- raster(extent(bb_extra), crs = projection(sa_proj), res = 4000)
  
  #'  Transform raster into polygons for visualization purposes
  #'  And add grid cell numbers to keep track of everything
  poly_ref <- rasterToPolygons(ref)
  cells <- 1:ncell(poly_ref)
  poly_ref@data$cell_ID <- c(1:length(cells))
  
  #'  Rasterize grid back into a raster for later data extraction
  ref_rast <- rasterize(poly_ref, ref)
  
  #'  Mask extended NE study area to the larger reference raster
  NE_mask <- mask(ref_rast, NE_box)
  #'  Crop masked raster to just the extent of the NE study area
  NE_crop <- crop(NE_mask, extent(NE_box))
  #'  Convert to a polygon
  NE_grid <- rasterToPolygons(NE_crop)
  #'  Count the number of grid cells in this new raster
  NE_cells <- ncell(NE_grid)
  #'  Append the new grid cell count to the NE raster so each cell has a unique
  #'  sequential number that's different from the much larger reference raster
  NE_grid@data$NE_cell_ID <- c(1:length(NE_cells))
  #'  Rasterize NE grid back into a raster
  NE_rast <- rasterize(NE_grid, NE_crop)

  #'  Check it out
  area(ref_rast)
  area(NE_rast)
  summary(ref_rast)
  summary(NE_rast)
  #'  What are these NA's? Are they cells that don't fully overlap NE_box?

  #  Plot the new rasters
  plot(ref_rast)
  plot(NE_box, add = T)
  plot(NE_SA, add = T, col = "blue")
  
  plot(NE_rast)
  plot(poly_ref, add = T)
  plot(NE_SA, add = T)
  plot(NE_box, add = T)
  
  #'  Save for other analyses
  #'  Use ".img" and "HFA" if planning to ever use this file in ArcGIS
  #'  Otherwise ".grid" and "raster" are better for working in R 
  writeRaster(ref_rast, filename = "./Shapefiles/ref_grid_4k.img", format = 'HFA', overwrite = T)
  writeRaster(NE_rast, filename = "./Shapefiles/NE_grid_4k.img", format = 'HFA', overwrite = T)
  