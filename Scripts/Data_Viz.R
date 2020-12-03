  #'  Data viz
  #'  SEFS 521
  #'  December 2021
  #'  --------------------------------------------
  #'  Get a feel for how these different data sets overlap and what the covariates
  #'  look like in these grid cells.
  #'  -------------------------------------------
  #'  Libraries
  library(sf)
  library(raster)
  library(tidyverse)
  
  #'  Read in telemetry & camera count data
  elk <- read.csv("./Elk_4hrFix_Cell_Count_withNAs 2020-11-25.csv") %>%
    dplyr::select(-X)
  coug <- read.csv("./Cougar_4hrFix_Cell_Count 2020-11-28.csv")
  cams <- read.csv("./Camera_detections.csv") %>%
    dplyr::select(-X)
  #'  Read in covariate data
  cov <- read.csv("./Covariates_by_cell_120220.csv")
  
  #'  Read in spatial data
  grid <- raster("./Shapefiles/NE_grid_4k.img")
  NE_sa <- st_read("./Shapefiles/fwdstudyareamaps", layer = "NE_SA") %>%
    st_transform(., crs = st_crs(grid))
  NE_box <- st_read("./Shapefiles/NE_covariate_area_2855", layer = "NE_covariate_area_2855")
  
  #'  Focus only on grid cells within MCPs
  mcp_elk <- elk[elk$MCP == 1,]
  mcp_coug <- coug[coug$MCP == 1,]
  all <- full_join(elk, coug, by = "cell") %>%
    mutate(
      total_MCP = ifelse(MCP.x == 1 | MCP.y == 1, 1, 0)
    )
  mcp_cams <- cams[all$total_MCP == 1,]
  #'  Double check that we don't lose any cameras that fall oustide total MCP
  #'  Should have 55 cameras
  colSums(mcp_cams, na.rm = TRUE)
  
  #'  Sum count data for each grid cell
  #'  Drop extra columns so they are excluded from sum function
  thin_obs <- dplyr::select(all, -c(cell, MCP.x, MCP.y, total_MCP))
  #'  Sum across rows (all observations in a grid cell, cougars & elk combined)
  sum_obs <- as.data.frame(rowSums(thin_obs))
  colnames(sum_obs) <- "sum_obs"
  #'  Tack that total count at end of data frame
  all <- cbind(all, sum_obs)
  summary(all)
  
  #'  Explore the observation and covariate data
  hist(all$sum_obs, xlim = c(0, 2500), ylim = c(0, 1000))
  #'  Really hard to see what's going on with so many cells with 0 observations
  #'  What about only grid cells within the larger MCP?
  mcp_obs <- all[all$total_MCP == 1,]
  hist(mcp_obs$sum_obs, xlim = c(0, 2500), ylim =  c(0, 300))
  #'  Now just the frequency of grid cells with 1+ observation
  #'  By default they are all within the larger MCP
  obs_only <- all$sum_obs[all$sum_obs > 0]
  hist(obs_only, xlim = c(0, 2500))
  
  #'  Check out the covariate data too!
  hist(cov$DEM_val, xlim = c(0, 2000), ylim = c(0, 200))
  hist(cov$roads_val)
  landcov <- table(cov$NLCD_label)
  barplot(landcov[order(landcov, decreasing = TRUE)], ylim = c(0, 800))
