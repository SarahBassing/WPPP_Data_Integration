  #'  Data viz
  #'  SEFS 521
  #'  December 2021
  #'  --------------------------------------------
  #'  Get a feel for how these different data sets overlap and what the covariates
  #'  look like in these grid cells.
  #'  -------------------------------------------

  #'  Clean out the environment
  rm(list = ls())

  #'  Libraries
  library(sf)
  library(raster)
  library(ggplot2)
  library(tidyverse)
  
  #'  Read in telemetry & camera count data
  elk <- read.csv("./Elk_4hrFix_Cell_Count_withNAs 2020-11-25.csv") %>%
    dplyr::select(-X)
  coug <- read.csv("./Cougar_4hrFix_Cell_Count 2020-11-28.csv")
  cams <- read.csv("./Camera_detections.csv") %>%
    dplyr::select(-X)
  #'  Read in covariate data
  cov <- read.csv("./Covariates_by_cell.csv")
  
  #'  Read in spatial data
  grid <- raster("./Shapefiles/NE_grid_4k.img")
  NE_sa <- st_read("./Shapefiles/fwdstudyareamaps", layer = "NE_SA") %>%
    st_transform(., crs = st_crs(grid))
  NE_box <- st_read("./Shapefiles/NE_covariate_area_2855", layer = "NE_covariate_area_2855")
  
  #'  Make count data and MCP info spatial
  #'  Extract coordinate data from grid (coordinates from center of each grid cell)
  coords <- coordinates(grid)
  #'  Append coordinates to the telemetry, camera, & covariate data
  elk <- cbind(elk, coords)
  coug <- cbind(coug, coords)
  cams <- cbind(cams, coords)
  cov <- cbind(cov, coords)
  
  #'  Focus only on grid cells within MCPs
  all <- full_join(elk, coug, by = c("cell", "x", "y")) %>%
    full_join(cams, by = c("cell", "x", "y")) %>%
    mutate(
      total_MCP = ifelse(MCP.x == 1 | MCP.y == 1, 1, 0)
    )
  mcp_elk <- elk[all$total_MCP == 1,]
  mcp_coug <- coug[all$total_MCP == 1,]
  mcp_cov <- cov[all$total_MCP == 1,]
  mcp_cams <- cams[all$total_MCP == 1,]
  #'  Double check that we don't lose any cameras that fall oustide total MCP
  #'  Should have 55 cameras
  colSums(mcp_cams, na.rm = TRUE)
  
  #'  Sum count data for each grid cell
  #'  Drop extra columns so they are excluded from sum function
  thin_obs <- dplyr::select(all, -c(cell, MCP.x, MCP.y, total_MCP, x, y, 
                                    Camera_Sampled, Cougar_Detections, Elk_Detections))
  
  elkSum <- thin_obs %>% 
    dplyr::select(contains("X")) %>% 
    rowSums()
  cougSum <- thin_obs %>%
    dplyr::select(contains("NE")) %>%
    rowSums()
  
  #'  Sum across rows (all observations in a grid cell, cougars & elk combined)
  sumAll <- as.data.frame(rowSums(thin_obs))
  sumElk <- as.data.frame(elkSum)
  sumCoug <- as.data.frame(cougSum)
  colnames(sumAll) <- "sumAll"
  colnames(sumElk) <- "sumElk"
  colnames(sumCoug) <- "sumCoug"
  #'  Tack that total count at end of data frame
  all <- cbind(all, sumElk, sumCoug, sumAll)
  
  #'  Same thing but for MCP grid cells only
  sumElk_mcp <- as.data.frame(sumElk[all$total_MCP == 1,])
  sumCoug_mcp <- as.data.frame(sumCoug[all$total_MCP == 1,])
  colnames(sumElk_mcp) <- "sumElk"
  colnames(sumCoug_mcp) <- "sumCoug"
  mcp_elk <- cbind(mcp_elk, sumElk_mcp)
  mcp_coug <- cbind(mcp_coug, sumCoug_mcp)
  
  # summary(all)
  
  #'  Explore the observation and covariate data
  hist(all$sumAll, xlim = c(0, 2500), ylim = c(0, 1000), breaks = 25)
  hist(all$sumElk, xlim = c(0, 2500), ylim = c(0, 1000), breaks = 25)
  hist(all$sumCoug, xlim = c(0, 500), ylim = c(0, 1000), breaks = 25)
  #'  Really hard to see what's going on with so many cells with 0 observations
  #'  What about only grid cells within the larger MCP?
  mcp_obs <- all[all$total_MCP == 1,]
  hist(mcp_obs$sumAll, xlim = c(0, 2500), ylim =  c(0, 300), breaks = 25)
  #'  Now just the frequency of grid cells with 1+ observation
  #'  By default they are all within the larger MCP
  obs_only <- all$sumAll[all$sumAll > 0]
  hist(obs_only, xlim = c(0, 2500), breaks = 25)
  elk_only <- all$sumElk[all$sumElk > 0]
  hist(elk_only, xlim = c(0, 2500), breaks = 25)
  coug_only <- all$sumCoug[all$sumCoug > 0]
  hist(coug_only, xlim = c(0, 500), ylim = c(0, 80), breaks = 25)
  
  #'  Check out the covariate data too!
  #'  Grid cells with telemetry observations
  hist(mcp_cov$DEM_val, ylim = c(0, 80))
  hist(mcp_cov$roads_val)
  landcov <- table(mcp_cov$NLCD_label)
  barplot(landcov[order(landcov, decreasing = TRUE)], ylim = c(0, 300))
  #'  Grid cells with camera observations
  hist(mcp_cov$DEM_val[mcp_cams$Camera_Sampled > 0], xlim = c(0, 1600))
  hist(mcp_cov$roads_val[mcp_cams$Camera_Sampled > 0], xlim = c(0, 50), ylim = c(0, 15))
  landcov <- table(mcp_cov$NLCD_label[mcp_cams$Camera_Sampled > 0])
  barplot(landcov[order(landcov, decreasing = TRUE)])
  
  #'  Map out the spatial distribution of these data!
  #'  Make observation data spatial
  elk_cells <- st_as_sf(mcp_elk, coords = c("x", "y"))
  coug_cells <- st_as_sf(mcp_coug, coords = c("x", "y"))
  cam_cells <- mcp_cams[!is.na(mcp_cams$Camera_Sampled),]
  cam_cells <- st_as_sf(cam_cells, coords = c("x", "y"))
  cov_cells <- st_as_sf(mcp_cov, coords = c("x", "y"))

  #'  Rasterize spatial observations to match grid cell raster
  elk_rast <- rasterize(elk_cells, grid)
  coug_rast <- rasterize(coug_cells, grid)
  cam_rast <- rasterize(cam_cells, grid)
  #cov_rast <- rasterize(cov_cells, grid)
  
  #'  Map out the variation in telemetry and camera observations
  plot(elk_rast$sumElk)
  plot(coug_rast$sumCoug)
  plot(cam_rast$Elk_Detections)
  plot(cam_rast$Cougar_Detections)
  
  #'  Plotting raster of raw location data in ggplot
  #'  Pull out single raster from larger rasterBrick
  elk_count <- elk_rast$sumElk
  coug_count <- coug_rast$sumCoug
  elk_det <- cam_rast$Elk_Detections
  coug_det <- cam_rast$Cougar_Detections
  #'  Convert raster to a spatial points data frame
  elk_count_spdf <- rasterToPoints(elk_count, spatial = TRUE)
  coug_count_spdf <- rasterToPoints(coug_count, spatial = TRUE)
  elk_det_spdf <- rasterToPoints(elk_det, spatial = TRUE)
  coug_det_spdf <- rasterToPoints(coug_det, spatial = TRUE)
  #'  Drop the spatial component so it's just a data frame with x, y data
  elk_count_df <- data.frame(elk_count_spdf) %>%
    dplyr::mutate(
      sumElk = ifelse(sumElk < 1, NA, sumElk)
    )
  coug_count_df <- data.frame(coug_count_spdf)%>%
    dplyr::mutate(
      sumCoug = ifelse(sumCoug < 1, NA, sumCoug)
    )
  elk_det_df <- data.frame(elk_det_spdf)  # Keep the 0's since they're real!
  coug_det_df <- data.frame(coug_det_spdf)  # Keep the 0's since they're real!
  #'  Remove the spatial data so there's no confusion
  rm(elk_count_spdf, elk_count)
  rm(coug_count_spdf, coug_count)
  rm(elk_det_spdf, elk_det)
  rm(coug_det_spdf, coug_det)
 
  require(scales)

  ggplot() +
    geom_sf(data = NE_sa) +
    geom_raster(data = elk_count_df, aes(x = x, y = y, fill = sumElk)) +
    scale_fill_distiller(palette = "YlGn", na.value = "transparent", trans = "reverse") +
    labs(title = "Total number of elk telemetry locations per grid cell")

  ggplot() +
    geom_sf(data = NE_sa) +
    geom_raster(data = coug_count_df, aes(x = x, y = y, fill = sumCoug)) +
    scale_fill_distiller(palette = "OrRd", na.value = "transparent", trans = "reverse") +
    labs(title = "Total number of cougar telemetry locations per grid cell")
  
  ggplot() +
    geom_sf(data = NE_sa) +
    geom_raster(data = elk_det_df, aes(x = x, y = y, fill = Elk_Detections)) +
    scale_fill_distiller(palette = "YlGn", trans = "reverse") +
    labs(title = "Total number of elk detections on camera per grid cell")
    
  ggplot() +
    geom_sf(data = NE_sa) +
    geom_raster(data = coug_det_df, aes(x = x, y = y, fill = Cougar_Detections)) +
    scale_fill_distiller(palette = "OrRd", trans = "reverse") +
    labs(title = "Total number of cougar detections on camera per grid cell")
  