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
  extra_SA <- readOGR("./Shapefiles/WPPP_CovariateBoundary", layer = "WPPP_CovariateBoundary")
  
  #'  Get everything in same projected coordinate system (UTMs)
  print(sa_proj <- projection(NE_SA))
  extra_SA <- spTransform(extra_SA, sa_proj)
  
  plot(extra_SA)
  plot(NE_SA, add = T)
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
  ref_grid <- rasterize(poly_ref, ref)
  
  #'  Save for other analyses
  #'  Use ".img" and "HFA" if planning to ever use this file in ArcGIS
  #'  Otherwise ".grid" and "raster" are better for working in R 
  writeRaster(ref_grid, filename = "./Shapefiles/ref_grid_4k.img", format = 'HFA', overwrite = T)