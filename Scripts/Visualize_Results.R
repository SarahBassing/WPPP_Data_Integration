  #'  Visualize Integrated RSF Results
  #'  SEFS 521
  #'  December 2020
  #'  --------------------------------------------
  #'  Plot lambda or probability of use across study area. Overlay animal locations
  #'  and detections to visualize results.
  #'  --------------------------------------------

  #'  Clean out the environment
  rm(list = ls())

  #'  Load libraries
  library(sf)
  library(ggplot2)
  library(raster)
  
  #'  Read in integrated model results for each combination of covariate & species
  load("./Output/elk_combo_elev_output.RData")
  load("./Output/elk_combo_road_output.RData")
  load("./Output/elk_combo_elev_road_output.RData")
  load("./Output/elk_combo_nlcd_output.RData")
  load("./Output/coug_combo_elev_output.RData")
  load("./Output/coug_combo_road_output.RData")
  load("./Output/coug_combo_elev_road_output.RData")
  load("./Output/coug_combo_nlcd_output.RData")
  
  #'  Read in telemetry- & camera-only model results for elev + road density model
  #'  for each species for comparison to integrated model results
  load("./Output/elk_telem_elev_road_output.RData")
  load("./Output/elk_cam_elev_road_output.RData")
  load("./Output/coug_telem_elev_road_output.RData")
  load("./Output/coug_cam_elev_road_output.RData")
  
  #'  Extract beta means
  belev_elk <- elk_combo_elev_output$BUGSoutput$mean$b_elev
  broad_elk <- elk_combo_road_output$BUGSoutput$mean$b_road
  belev2_elk <- elk_combo_elev_road_output$BUGSoutput$mean$b_elev
  broad2_elk <- elk_combo_elev_road_output$BUGSoutput$mean$b_road
  belev_coug <- coug_combo_elev_output$BUGSoutput$mean$b_elev
  broad_coug <- coug_combo_road_output$BUGSoutput$mean$b_road
  belev2_coug <- coug_combo_elev_road_output$BUGSoutput$mean$b_elev
  broad2_coug <- coug_combo_elev_road_output$BUGSoutput$mean$b_road
  #'  Extract 95% CI's of betas
  belev_elk_ci <- quantile(elk_combo_elev_output$BUGSoutput$sims.list$b_elev, c(0.025, 0.975))
  broad_elk_ci <- quantile(elk_combo_road_output$BUGSoutput$sims.list$b_road, c(0.025, 0.975))
  belev2_elk_ci <- quantile(elk_combo_elev_road_output$BUGSoutput$sims.list$b_elev, c(0.025, 0.975))
  broad2_elk_ci <- quantile(elk_combo_elev_road_output$BUGSoutput$sims.list$b_road, c(0.025, 0.975))
  belev_coug_ci <- quantile(coug_combo_elev_output$BUGSoutput$sims.list$b_elev, c(0.025, 0.975))
  broad_coug_ci <- quantile(coug_combo_road_output$BUGSoutput$sims.list$b_road, c(0.025, 0.975))
  belev2_coug_ci <- quantile(coug_combo_elev_road_output$BUGSoutput$sims.list$b_elev, c(0.025, 0.975))
  broad2_coug_ci <- quantile(coug_combo_elev_road_output$BUGSoutput$sims.list$b_road, c(0.025, 0.975))
  
  #'  Build table with integrated model results for both species
  species <- as.data.frame(c("elk", "elk", "elk", "elk", "cougar", "cougar", "cougar", "cougar"))
  mod <- rep(c("elev", "road", "elev + road", "elev + road"), 2)
  variable <- rep(c("elev", "road"), 4)
  beta <- round(c(belev_elk, broad_elk, belev2_elk, broad2_elk, belev_coug, 
                  broad_coug, belev2_coug, broad2_coug), 3)
  CI_2.5 <- round(c(belev_elk_ci[1], broad_elk_ci[1], belev2_elk_ci[1], broad2_elk_ci[1], 
                belev_coug_ci[1], broad_coug_ci[1], belev2_coug_ci[1], broad2_coug_ci[1]), 3)
  CI_97.5 <- round(c(belev_elk_ci[2], broad_elk_ci[2], belev2_elk_ci[2], broad2_elk_ci[2], 
                    belev_coug_ci[2], broad_coug_ci[2], belev2_coug_ci[2], broad2_coug_ci[2]), 3)
  tbl <- cbind(species, mod, variable, beta, CI_2.5, CI_97.5)
  colnames(tbl) <- c("Species", "Model Name", "Variable", "Estimate", "Lower CI", "Upper CI")
  
  #'  Extract intercepts for integrated & camera-only models for each species
  #'  Remember- no intercept for the telemetry-only model
  mu_a0_elk <- mean(elk_combo_elev_road_output$BUGSoutput$sims.list$a_cam)
  mu_a0_coug <- mean(coug_combo_elev_road_output$BUGSoutput$sims.list$a_cam)
  mu_a0_elk_cam <- mean(elk_cam_elev_road_output$BUGSoutput$sims.list$a_cam)
  mu_a0_coug_cam <- mean(coug_cam_elev_road_output$BUGSoutput$sims.list$a_cam)
  #'  And 95% CI on mean intercept for each species and model
  mu_a0_elk_ci <- quantile(elk_combo_elev_road_output$BUGSoutput$sims.list$a_cam, c(0.025, 0.975))
  mu_a0_coug_ci <- quantile(coug_combo_elev_road_output$BUGSoutput$sims.list$a_cam, c(0.025, 0.975))
  mu_a0_elk_cam_ci <- quantile(elk_cam_elev_road_output$BUGSoutput$sims.list$a_cam, c(0.025, 0.975))
  mu_a0_coug_cam_ci <- quantile(coug_cam_elev_road_output$BUGSoutput$sims.list$a_cam, c(0.025, 0.975))
  
  
  #'  Extract beta estimates & 95% CIs from telemetry-only and camera-only models
  #'  for each species 
  belev2_elk_tel <- elk_telem_elev_road_output$BUGSoutput$mean$b_elev
  broad2_elk_tel <- elk_telem_elev_road_output$BUGSoutput$mean$b_road
  belev2_elk_cam <- elk_cam_elev_road_output$BUGSoutput$mean$b_elev
  broad2_elk_cam <- elk_cam_elev_road_output$BUGSoutput$mean$b_road
  
  belev2_coug_tel <- coug_telem_elev_road_output$BUGSoutput$mean$b_elev
  broad2_coug_tel <- coug_telem_elev_road_output$BUGSoutput$mean$b_road
  belev2_coug_cam <- coug_cam_elev_road_output$BUGSoutput$mean$b_elev
  broad2_coug_cam <- coug_cam_elev_road_output$BUGSoutput$mean$b_road
  
  belev2_elk_tel_ci <- quantile(elk_telem_elev_road_output$BUGSoutput$sims.list$b_elev, c(0.025, 0.975))
  broad2_elk_tel_ci <- quantile(elk_telem_elev_road_output$BUGSoutput$sims.list$b_road, c(0.025, 0.975))
  belev2_elk_cam_ci <- quantile(elk_cam_elev_road_output$BUGSoutput$sims.list$b_elev, c(0.025, 0.975))
  broad2_elk_cam_ci <- quantile(elk_cam_elev_road_output$BUGSoutput$sims.list$b_road, c(0.025, 0.975))
  
  belev2_coug_tel_ci <- quantile(coug_telem_elev_road_output$BUGSoutput$sims.list$b_elev, c(0.025, 0.975))
  broad2_coug_tel_ci <- quantile(coug_telem_elev_road_output$BUGSoutput$sims.list$b_road, c(0.025, 0.975))
  belev2_coug_cam_ci <- quantile(coug_cam_elev_road_output$BUGSoutput$sims.list$b_elev, c(0.025, 0.975))
  broad2_coug_cam_ci <- quantile(coug_cam_elev_road_output$BUGSoutput$sims.list$b_road, c(0.025, 0.975))
  
  #'  Build table with results from all three models
  species <- as.data.frame(c("elk", "elk", "elk", "elk", "elk", "elk", "elk", 
                             "elk", "elk", "cougar", "cougar", "cougar", "cougar", 
                             "cougar", "cougar", "cougar", "cougar", "cougar"))
  mod <- rep(c("telem", "telem", "telem", "cam", "cam", "cam", "integrated", 
               "integrated", "integrated"), 2)
  param <- rep(c("Intercept", "B_elev", "B_road"), 6)
  est <- round(as.numeric(c("NA", belev2_elk_tel, broad2_elk_tel, mu_a0_elk_cam, belev2_elk_cam, 
                 broad2_elk_cam, mu_a0_elk, belev2_elk, broad2_elk, "NA", 
                 belev2_coug_tel, broad2_coug_tel, mu_a0_coug_cam, belev2_coug_cam, 
                 broad2_coug_cam, mu_a0_coug, belev2_coug, broad2_coug)), 3)
  CI_2.5 <- round(as.numeric(c("NA", belev2_elk_tel_ci[1], broad2_elk_tel_ci[1], mu_a0_elk_cam_ci[1],  
              belev2_elk_cam_ci[1], broad2_elk_cam_ci[1], mu_a0_elk_ci[1],
              belev2_elk_ci[1], broad2_elk_ci[1], "NA", belev2_coug_tel_ci[1], 
              broad2_coug_tel_ci[1], mu_a0_coug_cam_ci[1], belev2_coug_cam_ci[1], 
              broad2_coug_cam_ci[1], mu_a0_coug_ci[1], belev2_coug_ci[1], broad2_coug_ci[1])), 3)
  CI_97.5 <- round(as.numeric(c("NA", belev2_elk_tel_ci[2], broad2_elk_tel_ci[2], mu_a0_elk_cam_ci[2],  
               belev2_elk_cam_ci[2], broad2_elk_cam_ci[2], mu_a0_elk_ci[2],
               belev2_elk_ci[2], broad2_elk_ci[2], "NA", belev2_coug_tel_ci[2], 
               broad2_coug_tel_ci[2], mu_a0_coug_cam_ci[2], belev2_coug_cam_ci[2], 
               broad2_coug_cam_ci[2], mu_a0_coug_ci[2], belev2_coug_ci[2], broad2_coug_ci[2])), 3)
  tbl <- cbind(species, mod, param, est, CI_2.5, CI_97.5) 
  colnames(tbl) <- c("Species", "Model Name", "Parameter", "Estimate", 
                     "Lower 95% CI", "Upper 95% CI")
  
  #'  Save for sharing
  # write.csv(tbl, "./Output/Elev+Road_telem-cam-integrated_Results.csv")
  
  
  #'  Predict probability of use and intensity of use across MCP grid cell
  #'  using intercept from camera half of model and beta coefficient from joint model
  #'  Need to average across iterations & camera sites to get the mean intercept
  #'  Remember- no intercept for the telemetry-only model
  mu_a0_elk <- mean(elk_combo_elev_road_output$BUGSoutput$sims.list$a_cam)
  mu_a0_coug <- mean(coug_combo_elev_road_output$BUGSoutput$sims.list$a_cam)
  mu_a0_elk_cam <- mean(elk_cam_elev_road_output$BUGSoutput$sims.list$a_cam)
  mu_a0_coug_cam <- mean(coug_cam_elev_road_output$BUGSoutput$sims.list$a_cam)
  
  #'  Average across iterations to get mean beta coefficient
  b_elev_elk <- mean(elk_combo_elev_road_output$BUGSoutput$sims.list$b_elev)
  b_road_elk <- mean(elk_combo_elev_road_output$BUGSoutput$sims.list$b_road)
  b_elev_coug <- mean(coug_combo_elev_road_output$BUGSoutput$sims.list$b_elev)
  b_road_coug <- mean(coug_combo_elev_road_output$BUGSoutput$sims.list$b_road)
  b_elev_elk_tel <- mean(elk_telem_elev_road_output$BUGSoutput$sims.list$b_elev)
  b_road_elk_tel <- mean(elk_telem_elev_road_output$BUGSoutput$sims.list$b_road)
  b_elev_coug_tel <- mean(coug_telem_elev_road_output$BUGSoutput$sims.list$b_elev)
  b_road_coug_tel <- mean(coug_telem_elev_road_output$BUGSoutput$sims.list$b_road)
  b_elev_elk_cam <- mean(elk_cam_elev_road_output$BUGSoutput$sims.list$b_elev)
  b_road_elk_cam <- mean(elk_cam_elev_road_output$BUGSoutput$sims.list$b_road)
  b_elev_coug_cam <- mean(coug_cam_elev_road_output$BUGSoutput$sims.list$b_elev)
  b_road_coug_cam <- mean(coug_cam_elev_road_output$BUGSoutput$sims.list$b_road)
  
  
  #'  Gather covariate value for each grid cell we want to predict over
  require(tidyverse)
  cov <- read.csv("./Covariates_by_cell.csv") %>%
    mutate(
      zDEM = scale(DEM_val, center = TRUE, scale = TRUE),
      zRoads = scale(roads_val, center = TRUE, scale = TRUE)
    )
  
  #'  Function to predict across all grid cells based on mean intercept and 
  #'  effect of a given covariate.
  #'  Exponentiate linear model so predicted lambdas are on a meaningful scale,
  #'  Then convert lambda to a probability of use with 1-exp(-pred_lam)
  #'  Should end up with 1 predicted value per grid cell
  predict_use <- function(int, beta1, beta2, cov1, cov2) {
    predict_lam <- c()
    predict_prob <- c()
    for(i in 1:nrow(cov)) {
      predict_lam[i] <- exp(int + beta1*cov1[i] + beta2*cov2[i])
      #predict_lam[i] <- exp(mu_a0_elk + b_elev_elk*cov$zDEM[i])
      predict_prob[i] <- 1-exp(-predict_lam[i])
    }
    predict_prob <- as.data.frame(predict_prob)
    return(predict_prob)
  }
 
  #'  Extract desired covariate
  elev <- cov$zDEM
  road <- cov$zRoads
  
  #'  Run estimates from linear model through function to predict probability of use
  # elk_combo_el_rd_predict_elev <- predict_use(mu_a0_elk, b_elev_elk, elev)
  # elk_combo_el_rd_predict_road <- predict_use(mu_a0_elk, b_road_elk, road) 
  elk_combo_el_rd_predict_elev_road <- predict_use(mu_a0_elk, b_elev_elk, b_road_elk, elev, road)
  # elk_cam_el_rd_predict_elev <- predict_use(mu_a0_elk_cam, b_elev_elk, elev)  
  # elk_cam_el_rd_predict_road <- predict_use(mu_a0_elk_cam, b_road_elk, road)
  elk_cam_el_rd_predict_elev_road <- predict_use(mu_a0_elk_cam, b_elev_elk, b_road_elk, elev, road)
  #'  Feeding it 0 for intercept of telemetry-only model since no intercept estimated
  # elk_telem_el_rd_predict_elev <- predict_use(0, b_elev_elk, elev)
  # elk_telem_el_rd_predict_road <- predict_use(0, b_road_elk, road)  
  elk_telem_el_rd_predict_elev_road <- predict_use(0, b_elev_elk, b_road_elk, elev, road)
  
  # coug_combo_el_rd_predict_elev <- predict_use(mu_a0_coug, b_elev_coug, elev)  
  # coug_combo_el_rd_predict_road <- predict_use(mu_a0_coug, b_road_coug, road)
  coug_combo_el_rd_predict_elev_road <- predict_use(mu_a0_coug, b_elev_coug, b_road_coug, elev, road)
  # coug_cam_el_rd_predict_elev <- predict_use(mu_a0_coug_cam, b_elev_coug_cam, elev)  
  # coug_cam_el_rd_predict_road <- predict_use(mu_a0_coug_cam, b_road_coug_cam, road)
  coug_cam_el_rd_predict_elev_road <- predict_use(mu_a0_coug_cam, b_elev_coug, b_road_coug, elev, road)
  #'  Feeding it 0 for intercept of telemetry-only model since no intercept estimated
  # coug_telem_el_rd_predict_elev <- predict_use(0, b_elev_coug_tel, elev)  
  # coug_telem_el_rd_predict_road <- predict_use(0, b_road_coug_tel, road)
  coug_telem_el_rd_predict_elev_road <- predict_use(0, b_elev_coug, b_road_coug, elev, road)
  

  #'  Now plot these results!
  #'  Load spatial data
  grid <- raster("./Shapefiles/NE_grid_4k.img")
  NE_sa <- st_read("./Shapefiles/fwdstudyareamaps", layer = "NE_SA") %>%
    st_transform(., crs = st_crs(grid))
  NE_box <- st_read("./Shapefiles/NE_covariate_area_2855", layer = "NE_covariate_area_2855")
  
  #'  Inspect the results
  hist(elk_combo_el_rd_predict_elev$predict_prob)
  hist(coug_combo_el_rd_predict_elev$predict_prob)
  hist(elk_cam_el_rd_predict_elev$predict_prob)
  hist(coug_cam_el_rd_predict_elev$predict_prob)
  hist(elk_telem_el_rd_predict_elev$predict_prob)
  hist(coug_telem_el_rd_predict_elev$predict_prob)

  
  #'  Make model results spatial so they can be mapped
  #'  Grab coordinates of grid cells
  coords <- raster::coordinates(grid)
  cells <- 1:length(grid)
  
  #'  Function to merge grid cell locations and predicted probability of use, 
  #'  make spatial and rasterize, then strip spatial aspect so it can
  #'  be used to map with ggplot
  rasterize_predictions <- function(cells, predictions, coords, grid) {
    #'  Append coordinates to predictions
    bind_pred_coord <- cbind(cells, predictions, coords)
    #'  Turn predictions into spatial objects
    predicted_sf <- st_as_sf(bind_pred_coord, coords = c("x", "y"))
    #'  Rasterize spatial observations to match grid cell raster
    predicted_rast <- rasterize(predicted_sf, grid)
    #'  Pull out single raster from larger rasterBrick
    predict_prob <- predicted_rast$predict_prob 
    #'  Convert raster to a spatial points data frame
    predict_prob_spdf <- rasterToPoints(predict_prob, spatial = TRUE)
    #'  Drop the spatial component so it's just a data frame with x, y data
    predict_prob_df <- data.frame(predict_prob_spdf)
    
    return(predict_prob_df)
  }
  
  #'  Run predictions through function to rasterize and prepare for ggplot
  # elk_combo_df_elev <- rasterize_predictions(cells, elk_combo_el_rd_predict_elev, coords, grid)
  # elk_combo_df_road <- rasterize_predictions(cells, elk_combo_el_rd_predict_road, coords, grid)
  elk_combo_df_elev_road <- rasterize_predictions(cells, elk_combo_el_rd_predict_elev_road, coords, grid)
  # elk_telem_df_elev <- rasterize_predictions(cells, elk_telem_el_rd_predict_elev, coords, grid)
  # elk_telem_df_road <- rasterize_predictions(cells, elk_telem_el_rd_predict_road, coords, grid)
  elk_telem_df_elev_road <- rasterize_predictions(cells, elk_telem_el_rd_predict_elev_road, coords, grid)
  # elk_cam_df_elev <- rasterize_predictions(cells, elk_cam_el_rd_predict_elev, coords, grid)
  # elk_cam_df_road <- rasterize_predictions(cells, elk_cam_el_rd_predict_road, coords, grid)
  elk_cam_df_elev_road <- rasterize_predictions(cells, elk_cam_el_rd_predict_elev_road, coords, grid)

  # coug_combo_df_elev <- rasterize_predictions(cells, coug_combo_el_rd_predict_elev, coords, grid)
  # coug_combo_df_road <- rasterize_predictions(cells, coug_combo_el_rd_predict_road, coords, grid)
  coug_combo_df_elev_road <- rasterize_predictions(cells, coug_combo_el_rd_predict_elev_road, coords, grid)
  # coug_telem_df_elev <- rasterize_predictions(cells, coug_telem_el_rd_predict_elev, coords, grid)
  # coug_telem_df_road <- rasterize_predictions(cells, coug_telem_el_rd_predict_road, coords, grid)
  coug_telem_df_elev_road <- rasterize_predictions(cells, coug_telem_el_rd_predict_elev_road, coords, grid)
  # coug_cam_df_elev <- rasterize_predictions(cells, coug_cam_el_rd_predict_elev, coords, grid)
  # coug_cam_df_road <- rasterize_predictions(cells, coug_cam_el_rd_predict_road, coords, grid)
  coug_cam_df_elev_road <- rasterize_predictions(cells, coug_cam_el_rd_predict_elev_road, coords, grid)
  
 
  #'  Function to plot predicted probability of use across study area based on
  #'  model outputs
  require(scales)
  

  #'  Function to plot ELK maps
  plot_elk_use <- function(studyarea, prob_df_combo, prob_df_tel, prob_df_cam) {
    combo <- ggplot() +
      geom_sf(data = studyarea) +
      geom_raster(data = prob_df_combo, aes(x = x, y = y, fill = predict_prob)) +
      geom_sf(data = NE_sa, fill = NA, color = "black") +
      scale_fill_distiller(palette = "YlGn", trans = "reverse") +
      labs(title = "Probability of Use by Elk, Integrated Data Model")
    
    telem <- ggplot() +
      geom_sf(data = studyarea) +
      geom_raster(data = prob_df_tel, aes(x = x, y = y, fill = predict_prob)) +
      geom_sf(data = NE_sa, fill = NA, color = "black") +
      scale_fill_distiller(palette = "YlGn", trans = "reverse") +
      labs(title = "Relative Probability of Use by Elk, Telemetry-Only Model")
    
    cam <- ggplot() +
      geom_sf(data = studyarea) +
      geom_raster(data = prob_df_cam, aes(x = x, y = y, fill = predict_prob)) +
      geom_sf(data = NE_sa, fill = NA, color = "black") +
      scale_fill_distiller(palette = "YlGn", trans = "reverse") +
      labs(title = "Probability of Use by Elk, Camera-Only Model")
    
    plot(combo)
    plot(telem)
    plot(cam)

  }
  
  # pdf("./Output/Map_Elk_Probability_Use_Elevation.pdf")
  # plot_elk_use(NE_sa, elk_combo_df_elev, elk_telem_df_elev, elk_cam_df_elev)
  # dev.off()
  # pdf("./Output/Map_Elk_Probability_Use_RoadDensity.pdf")
  # plot_elk_use(NE_sa, elk_combo_df_road, elk_telem_df_road, elk_cam_df_road)
  # dev.off()
  pdf("./Output/Map_Elk_Probability_Use_Elevation&RoadDensity.pdf")
  plot_elk_use(NE_sa, elk_combo_df_elev_road, elk_telem_df_elev_road, elk_cam_df_elev_road)
  dev.off()

  #'  Function to plot COUGAR maps
  plot_coug_use <- function(studyarea, prob_df_combo, prob_df_tel, prob_df_cam) {
    combo <- ggplot() +
      geom_sf(data = studyarea) +
      geom_raster(data = prob_df_combo, aes(x = x, y = y, fill = predict_prob)) +
      geom_sf(data = NE_sa, fill = NA, color = "black") +
      scale_fill_distiller(palette = "OrRd", trans = "reverse") +
      labs(title = "Probability of Use by Cougar, Integrated Data Model")
    
    telem <- ggplot() +
      geom_sf(data = studyarea) +
      geom_raster(data = prob_df_tel, aes(x = x, y = y, fill = predict_prob)) +
      geom_sf(data = NE_sa, fill = NA, color = "black") +
      scale_fill_distiller(palette = "OrRd", trans = "reverse") +
      labs(title = "Relative Probability of Use by Cougar, Telemetry-Only Model")
    
    cam <- ggplot() +
      geom_sf(data = studyarea) +
      geom_raster(data = prob_df_cam, aes(x = x, y = y, fill = predict_prob)) +
      geom_sf(data = NE_sa, fill = NA, color = "black") +
      scale_fill_distiller(palette = "OrRd", trans = "reverse") +
      labs(title = "Probability of Use by Cougar, Camera-Only Model")
    
    plot(combo)
    plot(telem)
    plot(cam)
  }
  
  # pdf("./Output/Map_Cougar_Probability_Use_Elevation.pdf")
  # plot_coug_use(NE_sa, coug_combo_df_elev, coug_telem_df_elev, coug_cam_df_elev)
  # dev.off()
  # pdf("./Output/Map_Cougar_Probability_Use_RoadDensity.pdf")
  # plot_coug_use(NE_sa, coug_combo_df_road, coug_telem_df_road, coug_cam_df_road)
  # dev.off()
  pdf("./Output/Map_Cougar_Probability_Use_Elevation&RoadDensity.pdf")
  plot_coug_use(NE_sa, coug_combo_df_elev_road, coug_telem_df_elev_road, coug_cam_df_elev_road)
  dev.off()
  
  
  #'  Testing how to plot results to get a continuous surface of predicted use
  
  tst <- ggplot() +
    geom_sf(data = NE_sa) +
    geom_raster(data = elk_combo_df_elev, aes(x = x, y = y, fill = predict_prob)) +
    geom_sf(data = NE_sa, fill = NA, color = "black") +
    scale_fill_distiller(palette = "OrRd", trans = "reverse") +
    labs(title = "Predicted Probability of use in a grid cell")  
  
  ggsave("test.png", plot = tst)