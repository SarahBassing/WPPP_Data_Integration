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
  camdat <- read.csv("full_camdata.csv") %>%
    mutate(
      DateTime = as.POSIXct(DateTime, 
                            format="%Y-%m-%d %H:%M:%S",tz="America/Los_Angeles")
    ) %>%
    dplyr::select(-X) %>%
    #'  Only keep images of animals (for now)
    filter(Animal == "TRUE" | Animal == "true") %>%
    #'  Drop images marked as "Animal" but with no Species classification
    filter(!is.na(Species))

  #  Double check I have the right number of cameras in the dataset
  cams <- unique(camdat$CameraLocation)
  length(cams)
  
  #'  Read in spatial data
  NE_grid <- raster("./Shapefiles/NE_grid_4k.img", layer = 'grid', crs = 2855)
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
  
  #'  Summarize detection data
  #'  Number of cameras that detected species of interest
  length(droplevels(unique(as.factor(coug_winter$CameraLocation))))
  droplevels(unique(as.factor(coug_winter$CameraLocation)))
  length(droplevels(unique(as.factor(elk_winter$CameraLocation))))
  droplevels(unique(as.factor(elk_winter$CameraLocation)))
  #'  Number of independent detection events
  length(droplevels(unique(as.factor(coug_winter$caps))))
  length(droplevels(unique(as.factor(elk_winter$caps))))
  
  #'  Make detections spatial and reproject to match grid
  cam_proj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  coug_sf <- st_as_sf(coug_winter, coords = c("Camera_Long", "Camera_Lat"), crs = cam_proj) %>%
    st_transform(crs(proj))
  elk_sf <- st_as_sf(elk_winter, coords = c("Camera_Long", "Camera_Lat"), crs = cam_proj) %>%
    st_transform(crs(proj))
  
  #  Make all camera locations spatial
  cam_locs <- camdat[!duplicated(camdat$CameraLocation),] %>%
    dplyr::select(CameraLocation, Year, Study_Area, Camera_Long, Camera_Lat) %>%
    st_as_sf(., coords = c("Camera_Long", "Camera_Lat"), crs = cam_proj) %>%
    st_transform(crs(proj))
  
  #png(file = "./Elk_camera_detections.png", width = 600, height = 600)
  plot(NE_grid, axes = TRUE, cex.lab = 1.2, cex.axis = 1.2, main = "Elk Detections on Camera")
  plot(NE_SA, add = TRUE)
  plot(cam_locs, add = TRUE, col = "black", pch = 3, cex = 1.75)
  plot(elk_sf$geometry, add = TRUE, col = "blue", pch = 19, cex = 1.5) # col = as.factor(elk_sf$CameraLocation)
  #dev.off()
  
  #png(file = "./Cougar_camera_detections.png", width = 600, height = 600)
  plot(NE_grid, axes = TRUE, cex.lab = 1.2, cex.axis = 1.2, main = "Cougar Detections on Camera")
  plot(NE_SA, add = TRUE)
  plot(cam_locs, add = TRUE, col = "black", pch = 3, cex = 1.75)
  plot(coug_sf$geometry, add = TRUE, col = "dark red", pch = 19, cex = 1.5) # col = as.factor(coug_sf$CameraLocation)
  #dev.off()
  
    
  #'  Extract count of independent camera locations in each grid cell
  #'  Convert sf objects to sp objects, required for rasterization step
  coug_sp <- as_Spatial(coug_sf)
  elk_sp <- as_Spatial(elk_sf)
  cam_sp <- as_Spatial(cam_locs)
  #'  Rasterize camera detections with grid
  coug_rast <- rasterize(coug_sp, NE_grid, field = coug_sp$caps, fun = "count")
  elk_rast <- rasterize(elk_sp, NE_grid, field = elk_sp$caps, fun = "count")
  cam_rast <- rasterize(cam_sp, NE_grid, field = cam_sp$CameraLocation, fun = "count")

  #'  Overwrite NAs with zero in raster cells that were not sampled
  #'  
  ####  THIS IS NOT DOING WHAT I WANT YET! NEED TO ONLY TURN NON-DETECTION CAMERA LOCATIONS TO 0  ####
 
  # coug_rast[is.na(coug_rast)] <- 0
  # elk_rast[is.na(elk_rast)] <- 0
  
  #'  Store values as a new data frame
  coug_det <- as.data.frame(values(coug_rast))
  colnames(coug_det) <- "coug_det"
  elk_det <- as.data.frame(values(elk_rast))
  colnames(elk_det) <- "elk_det"
  #'  Count number of cameras per grid cell
  cam_samp <- as.data.frame(values(cam_rast))
  colnames(cam_samp) <- "Camera_Sampled"
  
  #'  Combine grid cells with detection data
  #'  Change NA values in detection data to 0 if camera sampled site but no detections occured
  cell <- values(NE_grid)
  detections <- cbind(cell, cam_samp, coug_det, elk_det) %>%
    mutate(
      Cougar_Detections = ifelse(Camera_Sampled >= 1 & is.na(coug_det), 0, coug_det),
      Elk_Detections = ifelse(Camera_Sampled >= 1 & is.na(elk_det), 0, elk_det)
    ) %>%
    dplyr::select(-c(coug_det, elk_det))
  
  #'  Make sure you get back what you put in (ignore the cell total count)
  colSums(detections, na.rm = TRUE)
  
  #'  Final data set is a data frame with a row for each grid cell, the number of
  #'  camera traps sampling that grid cell ("Camera_Sampled"), and the number of
  #'  independent detection events of cougars and elk in each grid cell, respectively.
  write.csv(detections, "Camera_detections.csv")
