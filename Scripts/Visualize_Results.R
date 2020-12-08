  #'  Visualize Integrated RSF Results
  #'  SEFS 521
  #'  December 2020
  #'  --------------------------------------------
  #'  Plot lambda or probability of use across study area. Overlay animal locations
  #'  and detections to visualize results.
  #'  --------------------------------------------

  #'  Load libraries
  library(sf)
  library(ggplot2)
  library(raster)
  
  #'  Read in and pull out relevant data from integrated RSF model
  load("./Output/combo_road_output.RData")
  load("./Output/prob_use_tel.RData")
  load("./Output/prob_use_cam.RData")
  
  #'  Lambda estimate for each grid cell
  #'  Telemetry half of model
  mu_lam <- combo_road_output$BUGSoutput$mean$mu_lam # posterior mean
  mu_grid1 <- mean(combo_DEM_output$BUGSoutput$sims.list$mu_lam[,1]) # mean of all iterations for 1st grid cell
  #'  Camera half of model
  lamc <- combo_road_output$BUGSoutput$mean$lamc
  
  #'  Probability of use for each grid cell
  tel_prob <- prob_use_tel 
  cam_prob <- prob_use_cam
  
  #'  Merge into a single data frame for each data type
  #'  Important to retain the original grid cell numbers here!
  telem <- cbind(tel_prob, mu_lam)
  cam <- cbind(cam_prob, lam_cam)
  
  #'  Beta coefficient
  b_dem <- combo_road_output$BUGSoutput$mean$b1
  
  #'  Predict probability of use and intensity of use across MCP grid cell
  #'  using intercept from camera half of model and beta coefficient from joint model
  #'  Need to average across iterations & camera sites to get the mean intercept
  mu.a0 <- mean(combo_road_output$BUGSoutput$sims.list$a_cam)
  #'  Average across iterations to get mean beta coefficient
  b1_dem <- mean(combo_road_output$BUGSoutput$sims.list$b1)
  #'  Gather covariate value for each grid cell we want to predict over
  require(tidyverse)
  cov <- read.csv("./Covariates_by_cell.csv") %>%
    mutate(
      zDEM = scale(DEM_val, center = TRUE, scale = TRUE),
      zRoads = scale(roads_val, center = TRUE, scale = TRUE)
    )
  #'  Predict across all grid cells based on mean intercept and effect of covariate
  #'  Have to exponentiate linear model so predicted lambdas are on a meaningful scale
  #'  And convert lambda to a probability of use with 1-exp(-pred_lam)
  #'  Should end up with 1 predicted value per grid cell
  pred_lam <- c()
  pred_prob <- c()
  for(i in 1:nrow(cov)) {
    pred_lam[i] <- exp(mu.a0 + b1_dem*cov$zRoads[i])
    pred_prob[i] <- 1-exp(-pred_lam[i])
  }
  
  pred_lam <- as.data.frame(pred_lam)
  pred_prob <- as.data.frame(pred_prob)
  
  #'  Now plot these results!
  #'  Load spatial data
  grid <- raster("./Shapefiles/NE_grid_4k.img")
  NE_sa <- st_read("./Shapefiles/fwdstudyareamaps", layer = "NE_SA") %>%
    st_transform(., crs = st_crs(grid))
  NE_box <- st_read("./Shapefiles/NE_covariate_area_2855", layer = "NE_covariate_area_2855")
  
  #'  Inspect the results
  #'  Make sure to adjust the xlim & ylim depending on species!
  hist(mu_lam)
  hist(tel_prob$tel_prob)
  hist(lamc)
  hist(cam_prob$det_prob)
  hist(pred_lam$pred_lam)
  hist(pred_prob$pred_prob)
  #'  QUESTION: why is the model predicting no grid cells with 0 telemetry locations?
  #'  Is this a product of using the Poisson distribution (positive, non-zero)
  #'  Or a product of the intercept and thus lambda not converging?
  
  #'  Make model results spatial so they can be mapped
  #'  Grab coordinates of grid cells included in the MCP
  coords <- raster::coordinates(grid)
  cells <- 1:length(grid)
  mcp_coords <- coords[telem$cell,]
  cam_coords <- coords[cam$cell,]
  
  #'  Append coordinates to model results
  telem_out <- cbind(telem, mcp_coords)
  cam_out <- cbind(cam, cam_coords)
  lam_out <- cbind(cells, pred_lam, coords)
  prob_out <- cbind(cells, pred_prob, coords)
  
  #'  Turn model results into spatial objects
  telem_cells <- st_as_sf(telem_out, coords = c("x", "y"))
  cam_cells <- st_as_sf(cam_out, coords = c("x", "y"))
  lam_cells <- st_as_sf(lam_out, coords = c("x", "y"))
  prob_cells <- st_as_sf(prob_out, coords = c("x", "y"))
  
  #'  Rasterize spatial observations to match grid cell raster
  telem_rast <- rasterize(telem_cells, grid)
  cam_rast <- rasterize(cam_cells, grid)
  lam_rast <- rasterize(lam_cells, grid)
  prob_rast <- rasterize(prob_cells, grid)
  
  #'  Plotting raster of raw location data in ggplot
  #'  Pull out single raster from larger rasterBrick
  telem_prob <- telem_rast$tel_prob
  telem_lambda <- telem_rast$mu_lam
  cam_prob <- cam_rast$det_prob
  cam_lambda <- cam_rast$lamc
  pred_lambda <- lam_rast$pred_lam
  pred_prob <- prob_rast$pred_prob 
  
  #'  Convert raster to a spatial points data frame
  telem_prob_spdf <- rasterToPoints(telem_prob, spatial = TRUE)
  telem_lam_spdf <- rasterToPoints(telem_lambda, spatial = TRUE)
  cam_prob_spdf <- rasterToPoints(cam_prob, spatial = TRUE)
  cam_lam_spdf <- rasterToPoints(cam_lambda, spatial = TRUE)
  pred_lam_spdf <- rasterToPoints(pred_lambda, spatial = TRUE)
  pred_prob_spdf <- rasterToPoints(pred_prob, spatial = TRUE)

  #'  Drop the spatial component so it's just a data frame with x, y data
  telem_prob_df <- data.frame(telem_prob_spdf)
  telem_lam_df <- data.frame(telem_lam_spdf)
  cam_prob_df <- data.frame(cam_prob_spdf)
  cam_lam_df <- data.frame(cam_lam_spdf)
  pred_lam_df <- data.frame(pred_lam_spdf)
  pred_prob_df <- data.frame(pred_prob_spdf)

  require(scales)
  
  ggplot() +
    geom_sf(data = NE_sa) +
    geom_raster(data = telem_prob_df, aes(x = x, y = y, fill = tel_prob)) +
    scale_fill_distiller(palette = "YlGn", trans = "reverse") +
    labs(title = "Probability of use of grid cell by at least 1 telemetered animal")  
  
  ggplot() +
    geom_sf(data = NE_sa) +
    geom_raster(data = telem_lam_df, aes(x = x, y = y, fill = mu_lam)) +
    scale_fill_distiller(palette = "YlGn", trans = "reverse") +
    labs(title = "Intensity of use of grid cell by at least 1 telemetered animal")  
  
  ggplot() +
    geom_sf(data = NE_sa) +
    geom_raster(data = cam_prob_df, aes(x = x, y = y, fill = det_prob)) +
    scale_fill_distiller(palette = "OrRd", trans = "reverse") +
    labs(title = "Probability of use of camera site in a grid cell")  
  
  ggplot() +
    geom_sf(data = NE_sa) +
    geom_raster(data = cam_lam_df, aes(x = x, y = y, fill = lamc)) +
    scale_fill_distiller(palette = "OrRd", trans = "reverse") +
    labs(title = "Intensity of use of camera site in a grid cell")  
  

  #'  Then plot results to get a continuous surface of predicted use
  
  ggplot() +
    geom_sf(data = NE_sa) +
    geom_raster(data = pred_lam_df, aes(x = x, y = y, fill = pred_lam)) +
    scale_fill_distiller(palette = "OrRd", trans = "reverse") +
    labs(title = "Predicted intensity of use in a grid cell")  
  
  ggplot() +
    geom_sf(data = NE_sa) +
    geom_raster(data = pred_prob_df, aes(x = x, y = y, fill = pred_prob)) +
    scale_fill_distiller(palette = "OrRd", trans = "reverse") +
    labs(title = "Predicted Probability of use in a grid cell")  
  