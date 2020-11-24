  #'  Counting camera detections by cell
  #'  SEFS 521 Group Project
  #'  
  #'  Sarah Bassing
  #'  November 2020
  #'  ----------------------------------
  #'  Script takes photo capture data from WPPP camera traps, extracts detections
  #'  from user-defined independent events, and counts the number of independent
  #'  detections of a given species per grid cell.
  #'  
  #'  A detection is considered an independent event if a 30 minute interval passes  
  #'  between subsequent camera triggers on the same species at a camare site.
  #'  ----------------------------------
  
  #'  Load libraries
  library(lubridate)
  library(overlap)
  library(tidyverse)
  library(sf)
  library(rgdal)
  library(raster)
  
  #'  Read in and format camera detection data
  camdat <- read.csv("SEFS521_camdata.csv") %>%
    mutate(
      DateTime = as.POSIXct(DateTime, 
                            format="%Y-%m-%d %H:%M:%S",tz="America/Los_Angeles")
    ) %>%
    dplyr::select(-X)
  
  #'  Read in spatial data
  NE_grid2 <- raster("./Shapefiles/NE_grid_4k.img", layer = 'grid', crs = 2855)
  NE_SA <- readOGR("./Shapefiles/fwdstudyareamaps", layer = "NE_SA") %>%
    spTransform(crs(NE_grid))
  
  #'  Identify projection as NAD83 (harn) Washington North EPSG:2855
  print(proj <- proj4string(NE_grid))
  
  
  #' #'  Turns out this radian time step is not necessary, but helpful code to have!
  #' #'  -----------------------------------------
  #' #'  Identify independent detection events based on time
  #' #'  Requires converting times to radians since working with circular time data
  #' #'  To do this- change every unit of time to minutes => hr*60 & sec/60
  #' #'  1 minute time = 1 turn/1440; 1 minute of time = pi/720 rad
  #' #'  1 Radian = the angle subtended at the center of a circle by an arc whose 
  #' #'  length is equal to the circle's radius; 1 full turn = 2pi Radians
  #' #'  ------------------------------------------
  #' 
  #' #'  Pull out lat/long of each detection for sunTime conversion
  #' coords <- SpatialPoints(cbind(camdat$Camera_Long, camdat$Camera_Lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
  #' 
  #' #'  Add column where time is converted to radians 
  #' #'  Use minutes as unit of time then convert to rad (1 min = 1 turn/1440)
  #' #'  And add column where radians are converted to suntime to adjust for seasonality
  #' #'  in daylight... don't actually need this for SEFS251  project
  #' raddat <- camdat %>%
  #'   mutate(
  #'     radTime = (((hour(DateTime)*60 + minute(DateTime) + second(DateTime)/(60))/1440)*2*pi),
  #'     sunTime = sunTime(radTime, Dates = DateTime, Coords = coords)
  #'   ) %>%
  #'   relocate(radTime, .after = Time) %>%
  #'   #  Make sure data are in some relevant order based on site & detection time
  #'   arrange(CameraLocation, DateTime)
  
  
  #'  Extract independent detections
  #'  Create a column identifying whether each image is an "independent" event
  #'  If camera site is diff from previous row then give unique value. If not then...
  #'  If species detected is diff from previous row at same site then give unique value. If not then...
  #'  If DateTime is >30 min from previous DateTime at same site for same species then give unique value. If not then...
  #'  Capture value is the same as that in the previous row.
  dat <- arrange(camdat, CameraLocation, DateTime)
  caps <- c()
  caps[1] <- 1
  for (i in 2:nrow(dat)){
    if (dat$CameraLocation[i-1] != dat$CameraLocation[i]) caps[i] = i
    else (if (dat$Species[i-1] != dat$Species[i]) caps[i] = i
          else (if (difftime(dat$DateTime[i], dat$DateTime[i-1], units = c("mins")) > 30) caps[i] = i
                else caps[i] = caps[i-1]))
  }

  caps <- as.factor(caps)
  
  #'  Add new column to larger data set
  capdata <- cbind(as.data.frame(dat), caps)
  
  #'  Retain only the first iamge from each unique detection event  
  detections <- capdata %>% 
    group_by(caps) %>% 
    slice(1L) %>%
    ungroup()
  
  #'  Separate by species
  coug_cam <- as.data.frame(detections[detections$Species == "Cougar",])
  elk_cam <- as.data.frame(detections[detections$Species == "Elk",])
  
  #'  Winnow down to date range of interest (Dec. 1 2018 - March 31, 2019)
  coug_winter <- coug_cam %>%
    filter(DateTime >= "2018-12-1" & DateTime <= "2019-03-31") 
  elk_winter <- elk_cam %>%
    filter(DateTime >= "2018-12-1" & DateTime <= "2019-03-31") 
  
  #'  Make detections spatial and reproject to match grid
  cam_proj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  # coug_spdf <- SpatialPointsDataFrame(coords = coug_winter[,25:26], coug_winter, proj4string = cam_proj)
  # coug_proj <- spTransform(coug_spdf, proj)
  # elk_spdf <- SpatialPointsDataFrame(coords = elk_winter[,25:26], elk_winter, proj4string = cam_proj)
  # elk_proj <- spTransform(elk_spdf, proj)
  coug_sf <- st_as_sf(coug_winter, coords = c("Camera_Long", "Camera_Lat"), crs = cam_proj) %>%
    st_transform(crs(proj))
  elk_sf <- st_as_sf(elk_winter, coords = c("Camera_Long", "Camera_Lat"), crs = cam_proj) %>%
    st_transform(crs(proj))
  
  plot(NE_grid, axes = TRUE)
  plot(NE_SA, add = TRUE)
  plot(elk_sf$geometry, add = TRUE, col = as.factor(elk_sf$CameraLocation), pch = 19)
  plot(coug_sf$geometry, add = TRUE, col = as.factor(coug_sf$CameraLocation), pch = 19)
  
    
  #'  Extract count of independent camera locations in each grid cell
  #'  Convert sf objects to sp objects, required for rasterization step
  coug_sp <- as_Spatial(coug_sf)
  elk_sp <- as_Spatial(elk_sf)
  #'  Rasterize camera detections with grid
  coug_rast <- rasterize(coug_sp, NE_grid, field = coug_sp$caps, fun = "count")
  elk_rast <- rasterize(elk_sp, NE_grid, field = elk_sp$caps, fun = "count")

  #'  Overwrite NAs with zero in raster cells that were not sampled
  #'  
  ####  THIS IS NOT DOING WHAT I WANT YET! NEED TO ONLY TURN NON-DETECTION CAMERA LOCATIONS TO 0  ####
 
   coug_rast[is.na(coug_rast)] <- 0
  elk_rast[is.na(elk_rast)] <- 0
  #'  Store values as a new data frame
  coug_det <- as.data.frame(values(coug_rast))
  colnames(coug_det) <- "Cougar Detections"
  elk_det <- as.data.frame(values(elk_rast))
  colnames(elk_det) <- "Elk Detections"
  #'  Append grid cell number to detection data frames
  cell <- values(NE_grid)
  detections <- cbind(cell, coug_det, elk_det)
  
