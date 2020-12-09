  #'  Data Integration RSF Model
  #'  SEFS 521
  #'  November 2020
  #'  --------------------------------------------
  #'  Combines observations of cougar and elk, respectively, in the WPPP Northeast 
  #'  study area from Dec 1, 2018 - Mar 31, 2019. Observations collected from 55  
  #'  camera traps, 44 elk GPS collars & 24 cougar GPS collars under one habitat 
  #'  selection model. Script runs a series of models for each species and varies
  #'  the covariates included in the linear model. Each model series includes
  #'  an integrated model, a telemetry-only model, and a camera-only model to 
  #'  evaluate how results change across data source and when they are combined.
  #'  
  #'  Model written by Dr. Beth Gardner
  #'  --------------------------------------------

  #'  Clean out the environment
  rm(list = ls())

  #'  Source input data wrangling script so data are formatted for model correctly
  source("./Scripts/Data_Integration_Input_Wrangling.R")
  
  #'  Load libraries
  library(R2jags)
  load.module("glm")
  library(mcmcplots)
  library(tidyverse)
  
  #'  ---------------------------------------------------
  #'  Helpful definitions for input data and estimated/derived parameters
  #'  DATA
  #'  ngrid = number of grid cells
  #'  ncam = number of grid cells with camera traps*
  #'  n = number of telemetered animals
  #'  R = total number of telemetry locations for an individual animal
  #'  M = number of telemetry locations per individual per grid cell
  #'  y = total number of species-specific detections in a grid cell
  #'  
  #'  PATAMETERS
  #'  a_telem = random effect for each telemetered animal
  #'  pi = a categorical probability that each grid cell is used by an individual
  #'       animal; calculated by dividing the site-specific estimates of lambda
  #'       for a given individual by the mean lambda for that individual 
  #'       (averaged across all sites)
  #'  lam_telem = (lambda) the intensity or rate parameter representing the estimated  
  #'        mean number of telemetry locations in a grid cell for an individual 
  #'        animal during the study period
  #'  a_cam = random effect for each camera site
  #'  lam_cam = (lambda) the intensity or rate parameter representing the estimated
  #'         mean number of independent detection events of a given species at
  #'         a camera site during the study period
  #'  b1 = beta coefficient for the effect of a given covariate on the number of 
  #'       animal observations; NOTE: estimating the same b1 for both types 
  #'       of data to connect the two observation processes under one model
  #'  ---------------------------------------------------
  #'  *will want to update this eventually so it's the true number of cameras,
  #'  not grid cells with cameras

  ####  ELK  ####
  ####  Elk Elevation Model  ####
  #'  Combined model for both telemetry and camera data
  cat("
  model{
  
    #'  Telemetry half of the model
    for(i in 1:n){
    
      #'  Random intercept for each individual
      a_telem[i] ~ dnorm(int_telem, tau_telem)
      
      #'  The number of observed telemetry locations in a grid cell arises from a 
      #'  categorical distribution based on the probability that a given cell is  
      #'  used by that individual and the total number of its telemetry locations
      M[i,1:ngrid] ~ dmulti(pi[i,1:ngrid], R[i])
    
    
      for(j in 1:ngrid){
      
        #'  Estimate the site- and individual-specific values of lambda
        log(lam_telem[i,j]) = a_telem[i] + b_elev*telev[j]
        
        #'  Calculating pi based on lambda estimates
        pi[i,j] = lam_telem[i,j]/sum(lam_telem[i,1:ngrid])
      
      }
    }
    
    #'  Camera trap half of the model
    for(k in 1:ncam){
      
      #'  Random effect for each camera
      a_cam[k] ~ dnorm(int_cam, tau_cam)
      
      #'  The number of observed detections at a camera site arises from a 
      #'  Poisson distribution based on the intensity or mean number of independent
      #'  detection events for a given species during the study period 
      y[k] ~ dpois(lam_cam[k])
      
      #'  Estimate the camera-site specific values of lambda 
      log(lam_cam[k]) <- a_cam[k] + b_elev*celev[k]
    }
    
    
    #'  Priors
    int_telem ~ dnorm(0, 0.1)
    tau_telem ~ dunif(0, 10)
    int_cam ~ dnorm(0, 0.1)
    tau_cam ~ dunif(0, 4)
    b_elev ~ dnorm(0, 0.1)
    
    
    #'  Derived parameters

    #'  Averaged across telemetered animals
    #'  Expected number of total locations per grid cell
    for(j in 1:ngrid) {
      mu_lam[j] <- mean(lam_telem[,j])
    }
    
  }
  
  ", fill = TRUE, file = "combo.txt")
  
  
  #'  Arguments for jags
  data <- list(M = elk_telem, R = sumelk, telev = as.vector(mcp_cov$zDEM),
               celev = as.vector(elk_covs$zDEM), ngrid = ngrid, n = nelk,
               ncam = ncam, y = elk_cams)

  parameters = c('a_telem', 'a_cam', 'int_telem', 'tau_telem', 'int_cam', 'tau_cam',
                'mu_lam', 'lam_cam', 'b_elev')
  
  inits = function() {list(b_elev = rnorm(1))}

  
  #'  Call to jags
  out <- jags(data, inits, parameters, "combo.txt", 
              n.chains = 3, n.thin = 1, n.iter = 6000, n.burnin = 3000)
  print(out, dig = 3)
  # mcmcplot(out)
  
  
  #'  Hold on to model output
  elk_combo_elev_output <- out
  save(elk_combo_elev_output, file = "./Output/elk_combo_elev_output.RData")
  
  #'  Put everything on the probability scale to estimate Probability of Use, which
  #'  represents the probability that a cell is used by an individual (telemetry)
  #'  and the probability that a camera site is used in a grid cell (camera)
  #'  during the study period.
  
  #'  Extract all iterations of telemetry- and camera-based lambdas
  mu_lam <- elk_combo_elev_output$BUGSoutput$sims.list$mu_lam
  lam_cam <- elk_combo_elev_output$BUGSoutput$sims.list$lam_cam
  
  #'  Loop through each iteration for each grid cell and calculate probability
  #'  Telemetry results
  tel_prob <- matrix(nrow = nrow(mu_lam), ncol = ngrid)
  for(i in 1:nrow(mu_lam)){
    for(j in 1:ngrid){
      tel_prob[i,j] <- 1-exp(-(mu_lam[i,j]))
    }
  }
  #'  Camera results
  det_prob <- matrix(nrow = nrow(lam_cam), ncol = ncam)
  for(i in 1:nrow(lam_cam)){
    for(k in 1:ncam){
      det_prob[i,k] <- 1-exp(-(lam_cam[i,k]))
    }
  }
  #'  Transpose so rows = grid cells & save as a data frame
  tel_prob <- as.data.frame(t(tel_prob))
  det_prob <- as.data.frame(t(det_prob))
  
  #' Save estimated probability of use based on each data source
  #' Pair probability of use with the original grid cell number it was estimated for
  cell <- mcp_elk$cell
  elk_pr_use_elev_tel <- cbind(cell, tel_prob)
  cell <- cam_and_covs$cell[!is.na(cam_and_covs$Elk_Detections)]
  elk_pr_use_elev_cam <- cbind(cell, det_prob)
  
  save(elk_pr_use_elev_tel, file = "./Output/elk_pr_use_elev_tel.RData")
  save(elk_pr_use_elev_cam, file = "./Output/elk_pr_use_elev_cam.RData")
  
  #'-------------------------------------------
  
  #'  RSF model with telemetry data only
  cat("
  model{
    
    for(i in 1:n){
    
      a_telem[i] ~ dnorm(int_telem, tau_telem)
      M[i,1:ngrid] ~ dmulti(pi[i,1:ngrid], R[i])
      
      for(j in 1:ngrid){
      
        log(lam_telem[i,j]) = a_telem[i] + b_elev*telev[j]
        pi[i,j] = lam_telem[i,j]/sum(lam_telem[i,1:ngrid])
      
      }
    }
    
    #'  Priors
    int_telem ~ dnorm(0, 0.1)
    tau_telem ~ dunif(0, 10)
    b_elev ~ dnorm(0, 0.1)

    
    #'  Derived parameters
    #'  Averaged across telemetered animals
    #'  Expected number of total locations per grid cell
    for(j in 1:ngrid) {
      mu_lam[j] <- mean(lam_telem[,j])
    }
  
  }
  
  ", fill = TRUE, file = "telem.txt")
  
  
  #'  Arguments for jags
  data <- list(M = elk_telem, R = sumelk, telev = as.vector(mcp_cov$zDEM),
               ngrid = ngrid, n = nelk)

  parameters = c('a_telem', 'int_telem', 'tau_telem', 'mu_lam', 'b_elev')

  inits = function() {list(b_elev = rnorm(1))}
  
  
  # call to jags
  out <- jags(data, inits, parameters, "telem.txt", 
              n.chains = 3, n.thin = 1, n.iter = 6000, n.burnin = 3000)
  print(out, dig = 3)
  # mcmcplot(out)
  
  #'  Hold on to model output
  elk_telem_elev_output <- out
  save(elk_telem_elev_output, file = "./Output/elk_telem_elev_output.RData")
  
  #'------------------------------------------
  
  #'  Camera data only model
  cat("
    model{
    
    for(k in 1:ncam){
      a_cam[k] ~ dnorm(int_cam, tau_cam)
      y[k] ~ dpois(lam_cam[k])
      log(lam_cam[k]) <- a_cam[k] + b_elev*celev[k]
    }
    
    #'  Priors
    b_elev ~ dnorm(0, 0.1)
    int_cam ~ dnorm(0, 0.1)
    tau_cam ~ dunif(0, 4)
    
  }
  
  ", fill = TRUE, file = "cam.txt")
  
  #'  Arguments for jags
  data <- list(ncam = ncam, y = elk_cams, celev = as.vector(elk_covs$zDEM))
 
  parameters = c('a_cam', 'int_cam', 'tau_cam', 'lam_cam', 'b_elev') 
  
  inits = function() {list(b_elev = rnorm(1))}
  
  # call to jags
  out <- jags(data, inits, parameters, "cam.txt", 
              n.chains = 3, n.thin = 1, n.iter = 6000, n.burnin = 3000)
  print(out, dig = 2)
  # mcmcplot(out)
  
  #'  Hold on to model output
  elk_cam_elev_output <- out
  save(elk_cam_elev_output, file = "./Output/elk_cam_elev_output.RData")
  
  
  
  ####  Elk Road Density Model  ####
  cat("
  model{
  
    #'  Telemetry half of the model
    for(i in 1:n){
    
      #'  Random intercept for each individual
      a_telem[i] ~ dnorm(int_telem, tau_telem)
      
      #'  The number of observed telemetry locations in a grid cell arises from a 
      #'  categorical distribution based on the probability that a given cell is  
      #'  used by that individual and the total number of its telemetry locations
      M[i,1:ngrid] ~ dmulti(pi[i,1:ngrid], R[i])
    
    
      for(j in 1:ngrid){
      
        #'  Estimate the site- and individual-specific values of lambda
        log(lam_telem[i,j]) = a_telem[i] + b_road*troad[j]
        
        #'  Calculating pi based on lambda estimates
        pi[i,j] = lam_telem[i,j]/sum(lam_telem[i,1:ngrid])
      
      }
    }
    
    #'  Camera trap half of the model
    for(k in 1:ncam){
      
      #'  Random effect for each camera
      a_cam[k] ~ dnorm(int_cam, tau_cam)
      
      #'  The number of observed detections at a camera site arises from a 
      #'  Poisson distribution based on the intensity or mean number of independent
      #'  detection events for a given species during the study period 
      y[k] ~ dpois(lam_cam[k])
      
      #'  Estimate the camera-site specific values of lambda 
      log(lam_cam[k]) <- a_cam[k] + b_road*croad[k]
    }
    
    
    #'  Priors
    int_telem ~ dnorm(0, 0.1)
    tau_telem ~ dunif(0, 10)
    int_cam ~ dnorm(0, 0.1)
    tau_cam ~ dunif(0, 4)
    b_road ~ dnorm(0, 0.1)
    
    
    #'  Derived parameters

    #'  Averaged across telemetered animals
    #'  Expected number of total locations per grid cell
    for(j in 1:ngrid) {
      mu_lam[j] <- mean(lam_telem[,j])
    }
    
  }
  
  ", fill = TRUE, file = "combo.txt")
  
  
  #'  Arguments for jags
  data <- list(M = elk_telem, R = sumelk, ngrid = ngrid, n = nelk, ncam = ncam, 
               y = elk_cams, troad = as.vector(mcp_cov$zRoads),
               croad = as.vector(elk_covs$zRoads))
  
  parameters = c('a_telem', 'a_cam', 'int_telem', 'tau_telem', 'int_cam', 'tau_cam',
                 'mu_lam', 'lam_cam', 'b_road')
  
  inits = function() {list(b_road = rnorm(1))}
  
  
  #'  Call to jags
  out <- jags(data, inits, parameters, "combo.txt", 
              n.chains = 3, n.thin = 1, n.iter = 6000, n.burnin = 3000)
  print(out, dig = 3)
  # mcmcplot(out)
  
  
  #'  Hold on to model output
  elk_combo_road_output <- out
  save(elk_combo_road_output, file = "./Output/elk_combo_road_output.RData")
  
  #'  Put everything on the probability scale to estimate Probability of Use, which
  #'  represents the probability that a cell is used by an individual (telemetry)
  #'  and the probability that a camera site is used in a grid cell (camera)
  #'  during the study period.
  
  #'  Extract all iterations of telemetry- and camera-based lambdas
  mu_lam <- elk_combo_road_output$BUGSoutput$sims.list$mu_lam
  lam_cam <- elk_combo_road_output$BUGSoutput$sims.list$lam_cam
  
  #'  Loop through each iteration for each grid cell and calculate probability
  #'  Telemetry results
  tel_prob <- matrix(nrow = nrow(mu_lam), ncol = ngrid)
  for(i in 1:nrow(mu_lam)){
    for(j in 1:ngrid){
      tel_prob[i,j] <- 1-exp(-(mu_lam[i,j]))
    }
  }
  #'  Camera results
  det_prob <- matrix(nrow = nrow(lam_cam), ncol = ncam)
  for(i in 1:nrow(lam_cam)){
    for(k in 1:ncam){
      det_prob[i,k] <- 1-exp(-(lam_cam[i,k]))
    }
  }
  #'  Transpose so rows = grid cells & save as a data frame
  tel_prob <- as.data.frame(t(tel_prob))
  det_prob <- as.data.frame(t(det_prob))
  
  #' Save estimated probability of use based on each data source
  #' Pair probability of use with the original grid cell number it was estimated for
  cell <- mcp_elk$cell
  elk_pr_use_road_tel <- cbind(cell, tel_prob)
  cell <- cam_and_covs$cell[!is.na(cam_and_covs$Elk_Detections)]
  elk_pr_use_road_cam <- cbind(cell, det_prob)
  
  save(elk_pr_use_road_tel, file = "./Output/elk_pr_use_road_tel.RData")
  save(elk_pr_use_road_cam, file = "./Output/elk_pr_use_road_cam.RData")
  
  #'-------------------------------------------
  
  #'  RSF model with telemetry data only
  cat("
  model{
    
    for(i in 1:n){
    
      a_telem[i] ~ dnorm(int_telem, tau_telem)
      M[i,1:ngrid] ~ dmulti(pi[i,1:ngrid], R[i])
      
      for(j in 1:ngrid){
      
        log(lam_telem[i,j]) = a_telem[i] + b_road*troad[j]
        pi[i,j] = lam_telem[i,j]/sum(lam_telem[i,1:ngrid])
      
      }
    }
    
    #'  Priors
    int_telem ~ dnorm(0, 0.1)
    tau_telem ~ dunif(0, 10)
    b_road ~ dnorm(0, 0.1)

    
    #'  Derived parameters
    #'  Averaged across telemetered animals
    #'  Expected number of total locations per grid cell
    for(j in 1:ngrid) {
      mu_lam[j] <- mean(lam_telem[,j])
    }
  
  }
  
  ", fill = TRUE, file = "telem.txt")
  
  
  #'  Arguments for jags
  data <- list(M = elk_telem, R = sumelk, troad = as.vector(mcp_cov$zRoads),
               ngrid = ngrid, n = nelk)
  
  parameters = c('a_telem', 'int_telem', 'tau_telem', 'mu_lam', 'b_road')
  
  inits = function() {list(b_road = rnorm(1))}
  
  
  # call to jags
  out <- jags(data, inits, parameters, "telem.txt", 
              n.chains = 3, n.thin = 1, n.iter = 6000, n.burnin = 3000)
  print(out, dig = 3)
  # mcmcplot(out)
  
  #'  Hold on to model output
  elk_telem_road_output <- out
  save(elk_telem_road_output, file = "./Output/elk_telem_road_output.RData")
  
  #'------------------------------------------
  
  #'  Camera data only model
  cat("
    model{
    
    for(k in 1:ncam){
      a_cam[k] ~ dnorm(int_cam, tau_cam)
      y[k] ~ dpois(lam_cam[k])
      log(lam_cam[k]) <- a_cam[k] + b_road*croad[k]
    }
    
    #'  Priors
    b_road ~ dnorm(0, 0.1)
    int_cam ~ dnorm(0, 0.1)
    tau_cam ~ dunif(0, 4)
    
  }
  
  ", fill = TRUE, file = "cam.txt")

  
  #'  Arguments for jags
  data <- list(ncam = ncam, y = elk_cams, croad = as.vector(elk_covs$zRoads))
  
  parameters = c('a_cam', 'int_cam', 'tau_cam', 'lam_cam', 'b_road') 
  
  inits = function() {list(b_road = rnorm(1))}
  
  # call to jags
  out <- jags(data, inits, parameters, "cam.txt", 
              n.chains = 3, n.thin = 1, n.iter = 6000, n.burnin = 3000)
  print(out, dig = 2)
  # mcmcplot(out)
  
  #'  Hold on to model output
  elk_cam_road_output <- out
  save(elk_cam_road_output, file = "./Output/elk_cam_road_output.RData")
  

  
  
  ####  Elk Elevation & Road Density Model  ####
  cat("
  model{
  
    #'  Telemetry half of the model
    for(i in 1:n){
    
      #'  Random intercept for each individual
      a_telem[i] ~ dnorm(int_telem, tau_telem)
      
      #'  The number of observed telemetry locations in a grid cell arises from a 
      #'  categorical distribution based on the probability that a given cell is  
      #'  used by that individual and the total number of its telemetry locations
      M[i,1:ngrid] ~ dmulti(pi[i,1:ngrid], R[i])
    
    
      for(j in 1:ngrid){
      
        #'  Estimate the site- and individual-specific values of lambda
        log(lam_telem[i,j]) = a_telem[i] + b_elev*telev[j] + b_road*troad[j]
        
        #'  Calculating pi based on lambda estimates
        pi[i,j] = lam_telem[i,j]/sum(lam_telem[i,1:ngrid])
      
      }
    }
    
    #'  Camera trap half of the model
    for(k in 1:ncam){
      
      #'  Random effect for each camera
      a_cam[k] ~ dnorm(int_cam, tau_cam)
      
      #'  The number of observed detections at a camera site arises from a 
      #'  Poisson distribution based on the intensity or mean number of independent
      #'  detection events for a given species during the study period 
      y[k] ~ dpois(lam_cam[k])
      
      #'  Estimate the camera-site specific values of lambda 
      log(lam_cam[k]) <- a_cam[k] + b_elev*celev[k] + b_road*croad[k]
    }
    
    
    #'  Priors
    int_telem ~ dnorm(0, 0.1)
    tau_telem ~ dunif(0, 10)
    int_cam ~ dnorm(0, 0.1)
    tau_cam ~ dunif(0, 4)
    b_elev ~ dnorm(0, 0.1)
    b_road ~ dnorm(0, 0.1)
    
    
    #'  Derived parameters

    #'  Averaged across telemetered animals
    #'  Expected number of total locations per grid cell
    for(j in 1:ngrid) {
      mu_lam[j] <- mean(lam_telem[,j])
    }
    
  }
  
  ", fill = TRUE, file = "combo.txt")
  
  
  #'  Arguments for jags
  data <- list(M = elk_telem, R = sumelk, telev = as.vector(mcp_cov$zDEM),
               celev = as.vector(elk_covs$zDEM), troad = as.vector(mcp_cov$zRoads), 
               croad = as.vector(elk_covs$zRoads), ngrid = ngrid, n = nelk,
               ncam = ncam, y = elk_cams)
  
  parameters = c('a_telem', 'a_cam', 'int_telem', 'tau_telem', 'int_cam', 'tau_cam',
                 'mu_lam', 'lam_cam', 'b_elev', 'b_road')
  
  inits = function() {list(b_elev = rnorm(1), b_road = rnorm(1))}
  
  
  #'  Call to jags
  out <- jags(data, inits, parameters, "combo.txt", 
              n.chains = 3, n.thin = 1, n.iter = 6000, n.burnin = 3000)
  print(out, dig = 3)
  # mcmcplot(out)
  
  
  #'  Hold on to model output
  elk_combo_elev_road_output <- out
  save(elk_combo_elev_road_output, file = "./Output/elk_combo_elev_road_output.RData")
  
  #'  Put everything on the probability scale to estimate Probability of Use, which
  #'  represents the probability that a cell is used by an individual (telemetry)
  #'  and the probability that a camera site is used in a grid cell (camera)
  #'  during the study period.
  
  #'  Extract all iterations of telemetry- and camera-based lambdas
  mu_lam <- elk_combo_elev_road_output$BUGSoutput$sims.list$mu_lam
  lam_cam <- elk_combo_elev_road_output$BUGSoutput$sims.list$lam_cam
  
  #'  Loop through each iteration for each grid cell and calculate probability
  #'  Telemetry results
  tel_prob <- matrix(nrow = nrow(mu_lam), ncol = ngrid)
  for(i in 1:nrow(mu_lam)){
    for(j in 1:ngrid){
      tel_prob[i,j] <- 1-exp(-(mu_lam[i,j]))
    }
  }
  #'  Camera results
  det_prob <- matrix(nrow = nrow(lam_cam), ncol = ncam)
  for(i in 1:nrow(lam_cam)){
    for(k in 1:ncam){
      det_prob[i,k] <- 1-exp(-(lam_cam[i,k]))
    }
  }
  #'  Transpose so rows = grid cells & save as a data frame
  tel_prob <- as.data.frame(t(tel_prob))
  det_prob <- as.data.frame(t(det_prob))
  
  #' Save estimated probability of use based on each data source
  #' Pair probability of use with the original grid cell number it was estimated for
  cell <- mcp_elk$cell
  elk_pr_use_elev_road_tel <- cbind(cell, tel_prob)
  cell <- cam_and_covs$cell[!is.na(cam_and_covs$Elk_Detections)]
  elk_pr_use_elev_road_cam <- cbind(cell, det_prob)
  
  save(elk_pr_use_elev_road_tel, file = "./Output/elk_pr_use_elev_road_tel.RData")
  save(elk_pr_use_elev_road_cam, file = "./Output/elk_pr_use_elev_road_cam.RData")
  
  #'-------------------------------------------
  
  #'  RSF model with telemetry data only
  cat("
  model{
    
    for(i in 1:n){
    
      a_telem[i] ~ dnorm(int_telem, tau_telem)
      M[i,1:ngrid] ~ dmulti(pi[i,1:ngrid], R[i])
      
      for(j in 1:ngrid){
      
        log(lam_telem[i,j]) = a_telem[i] + b_elev*telev[j] + b_road*troad[j]
        pi[i,j] = lam_telem[i,j]/sum(lam_telem[i,1:ngrid])
      
      }
    }
    
    #'  Priors
    int_telem ~ dnorm(0, 0.1)
    tau_telem ~ dunif(0, 10)
    b_elev ~ dnorm(0, 0.1)
    b_road ~ dnorm(0, 0.1)

    
    #'  Derived parameters
    #'  Averaged across telemetered animals
    #'  Expected number of total locations per grid cell
    for(j in 1:ngrid) {
      mu_lam[j] <- mean(lam_telem[,j])
    }
  
  }
  
  ", fill = TRUE, file = "telem.txt")
  
  
  #'  Arguments for jags
  data <- list(M = elk_telem, R = sumelk, telev = as.vector(mcp_cov$zDEM), 
               troad = as.vector(mcp_cov$zRoads), ngrid = ngrid, n = nelk)
  
  parameters = c('a_telem', 'int_telem', 'tau_telem', 'mu_lam', 'b_elev', 'b_road')
  
  inits = function() {list(b_elev = rnorm(1), b_road = rnorm(1))}
  
  
  # call to jags
  out <- jags(data, inits, parameters, "telem.txt", 
              n.chains = 3, n.thin = 1, n.iter = 6000, n.burnin = 3000)
  print(out, dig = 3)
  # mcmcplot(out)
  
  #'  Hold on to model output
  elk_telem_elev_road_output <- out
  save(elk_telem_elev_road_output, file = "./Output/elk_telem_elev_road_output.RData")
  
  #'------------------------------------------
  
  #'  Camera data only model
  cat("
    model{
    
    for(k in 1:ncam){
      a_cam[k] ~ dnorm(int_cam, tau_cam)
      y[k] ~ dpois(lam_cam[k])
      log(lam_cam[k]) <- a_cam[k] + b_elev*celev[k] + b_road*croad[k]
    }
    
    #'  Priors
    b_elev ~ dnorm(0, 0.1)
    b_road ~ dnorm(0, 0.1)
    int_cam ~ dnorm(0, 0.1)
    tau_cam ~ dunif(0, 4)
    
  }
  
  ", fill = TRUE, file = "cam.txt")
  
  
  
  #'  Arguments for jags
  data <- list(ncam = ncam, y = elk_cams, celev = as.vector(elk_covs$zDEM), 
               croad = as.vector(elk_covs$zRoads))
  
  parameters = c('a_cam', 'int_cam', 'tau_cam', 'lam_cam', 'b_elev', 'b_road') 
  
  inits = function() {list(b_elev = rnorm(1), b_road = rnorm(1))}
  
  # call to jags
  out <- jags(data, inits, parameters, "cam.txt", 
              n.chains = 3, n.thin = 1, n.iter = 6000, n.burnin = 3000)
  print(out, dig = 2)
  # mcmcplot(out)
  
  #'  Hold on to model output
  elk_cam_elev_road_output <- out
  save(elk_cam_elev_road_output, file = "./Output/elk_cam_elev_road_output.RData")
  
  
  
  ####  Elk Land Cover Model  ####
  cat("
  model{
  
    #'  Telemetry half of the model
    for(i in 1:n){
    
      #'  Random intercept for each individual
      a_telem[i] ~ dnorm(int_telem, tau_telem)
      
      #'  The number of observed telemetry locations in a grid cell arises from a 
      #'  categorical distribution based on the probability that a given cell is  
      #'  used by that individual and the total number of its telemetry locations
      M[i,1:ngrid] ~ dmulti(pi[i,1:ngrid], R[i])
    
    
      for(j in 1:ngrid){
      
        #'  Estimate the site- and individual-specific values of lambda
        #'  Forest is reference variable for NLCD categorical covariate
        log(lam_telem[i,j]) = a_telem[i] + b_shrub*tshrub[j] + b_crop*tcrop[j] + 
            b_water*twater[j] + b_grass*tgrass[j] + b_other*tother[j]

        
        #'  Calculating pi based on lambda estimates
        pi[i,j] = lam_telem[i,j]/sum(lam_telem[i,1:ngrid])
      
      }
    }
    
    #'  Camera trap half of the model
    for(k in 1:ncam){
      
      #'  Random effect for each camera
      a_cam[k] ~ dnorm(int_cam, tau_cam)
      
      #'  The number of observed detections at a camera site arises from a 
      #'  Poisson distribution based on the intensity or mean number of independent
      #'  detection events for a given species during the study period 
      y[k] ~ dpois(lam_cam[k])
      
      #'  Estimate the camera-site specific values of lambda 
      #'  Forest is reference variable for NLCD categorical covariate
      log(lam_cam[k]) <- a_cam[k] + b_shrub*cshrub[k] + b_crop*ccrop[k]
    }
    
    
    #'  Priors
    int_telem ~ dnorm(0, 0.1)
    tau_telem ~ dunif(0, 10)
    int_cam ~ dnorm(0, 0.1)
    tau_cam ~ dunif(0, 4)
    b_shrub ~ dnorm(0, 0.1)
    b_crop ~ dnorm(0, 0.1)
    b_water ~ dnorm(0, 0.1)
    b_grass ~ dnorm(0, 0.1)
    b_other ~ dnorm(0, 0.1)
    
    
    #'  Derived parameters

    #'  Averaged across telemetered animals
    #'  Expected number of total locations per grid cell
    for(j in 1:ngrid) {
      mu_lam[j] <- mean(lam_telem[,j])
    }
    
  }
  
  ", fill = TRUE, file = "combo.txt")
  
  
  #'  Arguments for jags
  data <- list(M = elk_telem, R = sumelk, ngrid = ngrid, n = nelk, ncam = ncam, y = elk_cams,
               tshrub = as.vector(mcp_cov$shrub), tcrop = as.vector(mcp_cov$crops),
               twater = as.vector(mcp_cov$water), tgrass = as.vector(mcp_cov$grass),
               tother = as.vector(mcp_cov$other), cshrub = as.vector(elk_covs$shrub),
               ccrop = as.vector(elk_covs$crops))
  
  parameters = c('a_telem', 'a_cam', 'int_telem', 'tau_telem', 'int_cam', 'tau_cam',
                 'mu_lam', 'lam_cam', 'b_shrub', 'b_crop', 'b_water', 'b_grass', 'b_other')
  
  inits = function() {list(b_shrub = rnorm(1), b_crop = rnorm(1), b_water = rnorm(1),
                           b_grass = rnorm(1), b_other = rnorm(1))}
  
  
  #'  Call to jags
  out <- jags(data, inits, parameters, "combo.txt", 
              n.chains = 3, n.thin = 1, n.iter = 6000, n.burnin = 3000)
  print(out, dig = 3)
  # mcmcplot(out)
  
  
  #'  Hold on to model output
  elk_combo_nlcd_output <- out
  save(elk_combo_nlcd_output, file = "./Output/elk_combo_nlcd_output.RData")
  
  #'  Put everything on the probability scale to estimate Probability of Use, which
  #'  represents the probability that a cell is used by an individual (telemetry)
  #'  and the probability that a camera site is used in a grid cell (camera)
  #'  during the study period.
  
  #'  Extract all iterations of telemetry- and camera-based lambdas
  mu_lam <- elk_combo_nlcd_output$BUGSoutput$sims.list$mu_lam
  lam_cam <- elk_combo_nlcd_output$BUGSoutput$sims.list$lam_cam
  
  #'  Loop through each iteration for each grid cell and calculate probability
  #'  Telemetry results
  tel_prob <- matrix(nrow = nrow(mu_lam), ncol = ngrid)
  for(i in 1:nrow(mu_lam)){
    for(j in 1:ngrid){
      tel_prob[i,j] <- 1-exp(-(mu_lam[i,j]))
    }
  }
  #'  Camera results
  det_prob <- matrix(nrow = nrow(lam_cam), ncol = ncam)
  for(i in 1:nrow(lam_cam)){
    for(k in 1:ncam){
      det_prob[i,k] <- 1-exp(-(lam_cam[i,k]))
    }
  }
  #'  Transpose so rows = grid cells & save as a data frame
  tel_prob <- as.data.frame(t(tel_prob))
  det_prob <- as.data.frame(t(det_prob))
  
  #' Save estimated probability of use based on each data source
  #' Pair probability of use with the original grid cell number it was estimated for
  cell <- mcp_elk$cell
  elk_pr_use_nlcd_tel <- cbind(cell, tel_prob)
  cell <- cam_and_covs$cell[!is.na(cam_and_covs$Elk_Detections)]
  elk_pr_use_nlcd_cam <- cbind(cell, det_prob)
  
  save(elk_pr_use_nlcd_tel, file = "./Output/elk_pr_use_nlcd_tel.RData")
  save(elk_pr_use_nlcd_cam, file = "./Output/elk_pr_use_nlcd_cam.RData")
  
  #'-------------------------------------------
  
  #'  RSF model with telemetry data only
  cat("
  model{
    
    for(i in 1:n){
    
      a_telem[i] ~ dnorm(int_telem, tau_telem)
      M[i,1:ngrid] ~ dmulti(pi[i,1:ngrid], R[i])
      
      for(j in 1:ngrid){
      
        #'  Forest is reference variable for NLCD categorical covariate
        log(lam_telem[i,j]) = a_telem[i] + b_shrub*tshrub[j] + b_crop*tcrop[j] + 
            b_water*twater[j] + b_grass*tgrass[j] + b_other*tother[j]

        pi[i,j] = lam_telem[i,j]/sum(lam_telem[i,1:ngrid])
      
      }
    }
    
    #'  Priors
    int_telem ~ dnorm(0, 0.1)
    tau_telem ~ dunif(0, 10)
    b_shrub ~ dnorm(0, 0.1)
    b_crop ~ dnorm(0, 0.1)
    b_water ~ dnorm(0, 0.1)
    b_grass ~ dnorm(0, 0.1)
    b_other ~ dnorm(0, 0.1)

    
    #'  Derived parameters
    #'  Averaged across telemetered animals
    #'  Expected number of total locations per grid cell
    for(j in 1:ngrid) {
      mu_lam[j] <- mean(lam_telem[,j])
    }
  
  }
  
  ", fill = TRUE, file = "telem.txt")
  
  
  #'  Arguments for jags
  data <- list(M = elk_telem, R = sumelk, ngrid = ngrid, n = nelk, 
               tshrub = as.vector(mcp_cov$shrub), tcrop = as.vector(mcp_cov$crops),
               twater = as.vector(mcp_cov$water), tgrass = as.vector(mcp_cov$grass),
               tother = as.vector(mcp_cov$other))
  
  parameters = c('a_telem', 'int_telem', 'tau_telem', 'mu_lam', 
                 'b_shrub', 'b_crop', 'b_water', 'b_grass', 'b_other')
  
  inits = function() {list(b_shrub = rnorm(1), b_crop = rnorm(1), b_water = rnorm(1),
                           b_grass = rnorm(1), b_other = rnorm(1))}
  
  
  # call to jags
  out <- jags(data, inits, parameters, "telem.txt", 
              n.chains = 3, n.thin = 1, n.iter = 6000, n.burnin = 3000)
  print(out, dig = 3)
  # mcmcplot(out)
  
  #'  Hold on to model output
  elk_telem_nlcd_output <- out
  save(elk_telem_nlcd_output, file = "./Output/elk_telem_nlcd_output.RData")
  
  #'------------------------------------------
  
  #'  Camera data only model
  cat("
    model{
    
    for(k in 1:ncam){
    
      a_cam[k] ~ dnorm(int_cam, tau_cam)
      
      y[k] ~ dpois(lam_cam[k])
            
      #'  Forest is reference variable for NLCD categorical covariate
      log(lam_cam[k]) <- a_cam[k] + b_shrub*cshrub[k] + b_crop*ccrop[k]
    }
    
    #'  Priors
    int_cam ~ dnorm(0, 0.1)
    tau_cam ~ dunif(0, 4)
    b_shrub ~ dnorm(0, 0.1)
    b_crop ~ dnorm(0, 0.1)
    
  }
  
  ", fill = TRUE, file = "cam.txt")
  
  
  
  #'  Arguments for jags
  data <- list(ncam = ncam, y = elk_cams, cshrub = as.vector(elk_covs$shrub), 
               ccrop = as.vector(elk_covs$crops))
  
  parameters = c('a_cam', 'int_cam', 'tau_cam', 'lam_cam', 'b_shrub', 'b_crop') 
  
  inits = function() {list(b_shrub = rnorm(1), b_crop = rnorm(1))}
  
  # call to jags
  out <- jags(data, inits, parameters, "cam.txt", 
              n.chains = 3, n.thin = 1, n.iter = 6000, n.burnin = 3000)
  print(out, dig = 2)
  # mcmcplot(out)
  
  #'  Hold on to model output
  elk_cam_nlcd_output <- out
  save(elk_cam_nlcd_output, file = "./Output/elk_cam_nlcd_output.RData")
  
  
  
  ####  Elk Elevation, Road Density, & Land Cover Model  ####
  
  
  
  ####  COUGAR  ####
  ####  Cougar Elevation Model  ####
  #'  Combined model for both telemetry and camera data
  cat("
  model{
  
    #'  Telemetry half of the model
    for(i in 1:n){
    
      #'  Random intercept for each individual
      a_telem[i] ~ dnorm(int_telem, tau_telem)
      
      #'  The number of observed telemetry locations in a grid cell arises from a 
      #'  categorical distribution based on the probability that a given cell is  
      #'  used by that individual and the total number of its telemetry locations
      M[i,1:ngrid] ~ dmulti(pi[i,1:ngrid], R[i])
    
    
      for(j in 1:ngrid){
      
        #'  Estimate the site- and individual-specific values of lambda
        log(lam_telem[i,j]) = a_telem[i] + b_elev*telev[j]
        
        #'  Calculating pi based on lambda estimates
        pi[i,j] = lam_telem[i,j]/sum(lam_telem[i,1:ngrid])
      
      }
    }
    
    #'  Camera trap half of the model
    for(k in 1:ncam){
      
      #'  Random effect for each camera
      a_cam[k] ~ dnorm(int_cam, tau_cam)
      
      #'  The number of observed detections at a camera site arises from a 
      #'  Poisson distribution based on the intensity or mean number of independent
      #'  detection events for a given species during the study period 
      y[k] ~ dpois(lam_cam[k])
      
      #'  Estimate the camera-site specific values of lambda 
      log(lam_cam[k]) <- a_cam[k] + b_elev*celev[k]
    }
    
    
    #'  Priors
    int_telem ~ dnorm(0, 0.1)
    tau_telem ~ dunif(0, 10)
    int_cam ~ dnorm(0, 0.1)
    tau_cam ~ dunif(0, 4)
    b_elev ~ dnorm(0, 0.1)
    
    
    #'  Derived parameters

    #'  Averaged across telemetered animals
    #'  Expected number of total locations per grid cell
    for(j in 1:ngrid) {
      mu_lam[j] <- mean(lam_telem[,j])
    }
    
  }
  
  ", fill = TRUE, file = "combo.txt")
  
  
  #'  Arguments for jags
  data <- list(M = coug_telem, R = sumcoug, telev = as.vector(mcp_cov$zDEM),
               celev = as.vector(coug_covs$zDEM), ngrid = ngrid, n = ncoug,
               ncam = ncam, y = coug_cams)
  
  parameters = c('a_telem', 'a_cam', 'int_telem', 'tau_telem', 'int_cam', 'tau_cam',
                 'mu_lam', 'lam_cam', 'b_elev')
  
  inits = function() {list(b_elev = rnorm(1))}
  
  
  #'  Call to jags
  out <- jags(data, inits, parameters, "combo.txt", 
              n.chains = 3, n.thin = 1, n.iter = 6000, n.burnin = 3000)
  print(out, dig = 3)
  # mcmcplot(out)
  
  
  #'  Hold on to model output
  coug_combo_elev_output <- out
  save(coug_combo_elev_output, file = "./Output/coug_combo_elev_output.RData")
  
  #'  Put everything on the probability scale to estimate Probability of Use, which
  #'  represents the probability that a cell is used by an individual (telemetry)
  #'  and the probability that a camera site is used in a grid cell (camera)
  #'  during the study period.
  
  #'  Extract all iterations of telemetry- and camera-based lambdas
  mu_lam <- coug_combo_elev_output$BUGSoutput$sims.list$mu_lam
  lam_cam <- coug_combo_elev_output$BUGSoutput$sims.list$lam_cam
  
  #'  Loop through each iteration for each grid cell and calculate probability
  #'  Telemetry results
  tel_prob <- matrix(nrow = nrow(mu_lam), ncol = ngrid)
  for(i in 1:nrow(mu_lam)){
    for(j in 1:ngrid){
      tel_prob[i,j] <- 1-exp(-(mu_lam[i,j]))
    }
  }
  #'  Camera results
  det_prob <- matrix(nrow = nrow(lam_cam), ncol = ncam)
  for(i in 1:nrow(lam_cam)){
    for(k in 1:ncam){
      det_prob[i,k] <- 1-exp(-(lam_cam[i,k]))
    }
  }
  #'  Transpose so rows = grid cells & save as a data frame
  tel_prob <- as.data.frame(t(tel_prob))
  det_prob <- as.data.frame(t(det_prob))
  
  #' Save estimated probability of use based on each data source
  #' Pair probability of use with the original grid cell number it was estimated for
  cell <- mcp_coug$cell
  coug_pr_use_elev_tel <- cbind(cell, tel_prob)
  cell <- cam_and_covs$cell[!is.na(cam_and_covs$Cougar_Detections)]
  coug_pr_use_elev_cam <- cbind(cell, det_prob)
  
  save(coug_pr_use_elev_tel, file = "./Output/coug_pr_use_elev_tel.RData")
  save(coug_pr_use_elev_cam, file = "./Output/coug_pr_use_elev_cam.RData")
  
  #'-------------------------------------------
  
  #'  RSF model with telemetry data only
  cat("
  model{
    
    for(i in 1:n){
    
      a_telem[i] ~ dnorm(int_telem, tau_telem)
      M[i,1:ngrid] ~ dmulti(pi[i,1:ngrid], R[i])
      
      for(j in 1:ngrid){
      
        log(lam_telem[i,j]) = a_telem[i] + b_elev*telev[j]
        pi[i,j] = lam_telem[i,j]/sum(lam_telem[i,1:ngrid])
      
      }
    }
    
    #'  Priors
    int_telem ~ dnorm(0, 0.1)
    tau_telem ~ dunif(0, 10)
    b_elev ~ dnorm(0, 0.1)

    
    #'  Derived parameters
    #'  Averaged across telemetered animals
    #'  Expected number of total locations per grid cell
    for(j in 1:ngrid) {
      mu_lam[j] <- mean(lam_telem[,j])
    }
  
  }
  
  ", fill = TRUE, file = "telem.txt")
  
  
  #'  Arguments for jags
  data <- list(M = coug_telem, R = sumcoug, telev = as.vector(mcp_cov$zDEM),
               ngrid = ngrid, n = ncoug)
  
  parameters = c('a_telem', 'int_telem', 'tau_telem', 'mu_lam', 'b_elev')
  
  inits = function() {list(b_elev = rnorm(1))}
  
  
  # call to jags
  out <- jags(data, inits, parameters, "telem.txt", 
              n.chains = 3, n.thin = 1, n.iter = 6000, n.burnin = 3000)
  print(out, dig = 3)
  # mcmcplot(out)
  
  #'  Hold on to model output
  coug_telem_elev_output <- out
  save(coug_telem_elev_output, file = "./Output/coug_telem_elev_output.RData")
  
  #'------------------------------------------
  
  #'  Camera data only model
  cat("
    model{
    
    for(k in 1:ncam){
      a_cam[k] ~ dnorm(int_cam, tau_cam)
      y[k] ~ dpois(lam_cam[k])
      log(lam_cam[k]) <- a_cam[k] + b_elev*celev[k]
    }
    
    #'  Priors
    b_elev ~ dnorm(0, 0.1)
    int_cam ~ dnorm(0, 0.1)
    tau_cam ~ dunif(0, 4)
    
  }
  
  ", fill = TRUE, file = "cam.txt")
  
  #'  Arguments for jags
  data <- list(ncam = ncam, y = coug_cams, celev = as.vector(coug_covs$zDEM))
  
  parameters = c('a_cam', 'int_cam', 'tau_cam', 'lam_cam', 'b_elev') 
  
  inits = function() {list(b_elev = rnorm(1))}
  
  # call to jags
  out <- jags(data, inits, parameters, "cam.txt", 
              n.chains = 3, n.thin = 1, n.iter = 6000, n.burnin = 3000)
  print(out, dig = 2)
  # mcmcplot(out)
  
  #'  Hold on to model output
  coug_cam_elev_output <- out
  save(coug_cam_elev_output, file = "./Output/coug_cam_elev_output.RData")
  
  
  
  
  ####  Cougar Road Density Model  ####
  cat("
  model{
  
    #'  Telemetry half of the model
    for(i in 1:n){
    
      #'  Random intercept for each individual
      a_telem[i] ~ dnorm(int_telem, tau_telem)
      
      #'  The number of observed telemetry locations in a grid cell arises from a 
      #'  categorical distribution based on the probability that a given cell is  
      #'  used by that individual and the total number of its telemetry locations
      M[i,1:ngrid] ~ dmulti(pi[i,1:ngrid], R[i])
    
    
      for(j in 1:ngrid){
      
        #'  Estimate the site- and individual-specific values of lambda
        log(lam_telem[i,j]) = a_telem[i] + b_road*troad[j]
        
        #'  Calculating pi based on lambda estimates
        pi[i,j] = lam_telem[i,j]/sum(lam_telem[i,1:ngrid])
      
      }
    }
    
    #'  Camera trap half of the model
    for(k in 1:ncam){
      
      #'  Random effect for each camera
      a_cam[k] ~ dnorm(int_cam, tau_cam)
      
      #'  The number of observed detections at a camera site arises from a 
      #'  Poisson distribution based on the intensity or mean number of independent
      #'  detection events for a given species during the study period 
      y[k] ~ dpois(lam_cam[k])
      
      #'  Estimate the camera-site specific values of lambda 
      log(lam_cam[k]) <- a_cam[k] + b_road*croad[k]
    }
    
    
    #'  Priors
    int_telem ~ dnorm(0, 0.1)
    tau_telem ~ dunif(0, 10)
    int_cam ~ dnorm(0, 0.1)
    tau_cam ~ dunif(0, 4)
    b_road ~ dnorm(0, 0.1)
    
    
    #'  Derived parameters

    #'  Averaged across telemetered animals
    #'  Expected number of total locations per grid cell
    for(j in 1:ngrid) {
      mu_lam[j] <- mean(lam_telem[,j])
    }
    
  }
  
  ", fill = TRUE, file = "combo.txt")
  
  
  #'  Arguments for jags
  data <- list(M = coug_telem, R = sumcoug, troad = as.vector(mcp_cov$zRoads),
               croad = as.vector(coug_covs$zRoads), ngrid = ngrid, n = ncoug,
               ncam = ncam, y = coug_cams)
  
  parameters = c('a_telem', 'a_cam', 'int_telem', 'tau_telem', 'int_cam', 'tau_cam',
                 'mu_lam', 'lam_cam', 'b_road')
  
  inits = function() {list(b_road = rnorm(1))}
  
  
  #'  Call to jags
  out <- jags(data, inits, parameters, "combo.txt", 
              n.chains = 3, n.thin = 1, n.iter = 6000, n.burnin = 3000)
  print(out, dig = 3)
  # mcmcplot(out)
  
  
  #'  Hold on to model output
  coug_combo_road_output <- out
  save(coug_combo_road_output, file = "./Output/coug_combo_road_output.RData")
  
  #'  Put everything on the probability scale to estimate Probability of Use, which
  #'  represents the probability that a cell is used by an individual (telemetry)
  #'  and the probability that a camera site is used in a grid cell (camera)
  #'  during the study period.
  
  #'  Extract all iterations of telemetry- and camera-based lambdas
  mu_lam <- coug_combo_road_output$BUGSoutput$sims.list$mu_lam
  lam_cam <- coug_combo_road_output$BUGSoutput$sims.list$lam_cam
  
  #'  Loop through each iteration for each grid cell and calculate probability
  #'  Telemetry results
  tel_prob <- matrix(nrow = nrow(mu_lam), ncol = ngrid)
  for(i in 1:nrow(mu_lam)){
    for(j in 1:ngrid){
      tel_prob[i,j] <- 1-exp(-(mu_lam[i,j]))
    }
  }
  #'  Camera results
  det_prob <- matrix(nrow = nrow(lam_cam), ncol = ncam)
  for(i in 1:nrow(lam_cam)){
    for(k in 1:ncam){
      det_prob[i,k] <- 1-exp(-(lam_cam[i,k]))
    }
  }
  #'  Transpose so rows = grid cells & save as a data frame
  tel_prob <- as.data.frame(t(tel_prob))
  det_prob <- as.data.frame(t(det_prob))
  
  #' Save estimated probability of use based on each data source
  #' Pair probability of use with the original grid cell number it was estimated for
  cell <- mcp_coug$cell
  coug_pr_use_road_tel <- cbind(cell, tel_prob)
  cell <- cam_and_covs$cell[!is.na(cam_and_covs$Cougar_Detections)]
  coug_pr_use_road_cam <- cbind(cell, det_prob)
  
  save(coug_pr_use_road_tel, file = "./Output/coug_pr_use_road_tel.RData")
  save(coug_pr_use_road_cam, file = "./Output/coug_pr_use_road_cam.RData")
  
  #'-------------------------------------------
  
  #'  RSF model with telemetry data only
  cat("
  model{
    
    for(i in 1:n){
    
      a_telem[i] ~ dnorm(int_telem, tau_telem)
      M[i,1:ngrid] ~ dmulti(pi[i,1:ngrid], R[i])
      
      for(j in 1:ngrid){
      
        log(lam_telem[i,j]) = a_telem[i] + b_road*troad[j]
        pi[i,j] = lam_telem[i,j]/sum(lam_telem[i,1:ngrid])
      
      }
    }
    
    #'  Priors
    int_telem ~ dnorm(0, 0.1)
    tau_telem ~ dunif(0, 10)
    b_road ~ dnorm(0, 0.1)

    
    #'  Derived parameters
    #'  Averaged across telemetered animals
    #'  Expected number of total locations per grid cell
    for(j in 1:ngrid) {
      mu_lam[j] <- mean(lam_telem[,j])
    }
  
  }
  
  ", fill = TRUE, file = "telem.txt")
  
  
  #'  Arguments for jags
  data <- list(M = coug_telem, R = sumcoug, troad = as.vector(mcp_cov$zRoads),
               ngrid = ngrid, n = ncoug)
  
  parameters = c('a_telem', 'int_telem', 'tau_telem', 'mu_lam', 'b_road')
  
  inits = function() {list(b_road = rnorm(1))}
  
  
  # call to jags
  out <- jags(data, inits, parameters, "telem.txt", 
              n.chains = 3, n.thin = 1, n.iter = 6000, n.burnin = 3000)
  print(out, dig = 3)
  # mcmcplot(out)
  
  #'  Hold on to model output
  coug_telem_road_output <- out
  save(coug_telem_road_output, file = "./Output/coug_telem_road_output.RData")
  
  #'------------------------------------------
  
  #'  Camera data only model
  cat("
    model{
    
    for(k in 1:ncam){
      a_cam[k] ~ dnorm(int_cam, tau_cam)
      y[k] ~ dpois(lam_cam[k])
      log(lam_cam[k]) <- a_cam[k] + b_road*croad[k]
    }
    
    #'  Priors
    b_road ~ dnorm(0, 0.1)
    int_cam ~ dnorm(0, 0.1)
    tau_cam ~ dunif(0, 4)
    
  }
  
  ", fill = TRUE, file = "cam.txt")
  
  
  #'  Arguments for jags
  data <- list(ncam = ncam, y = coug_cams, croad = as.vector(coug_covs$zRoads))
  
  parameters = c('a_cam', 'int_cam', 'tau_cam', 'lam_cam', 'b_road') 
  
  inits = function() {list(b_road = rnorm(1))}
  
  # call to jags
  out <- jags(data, inits, parameters, "cam.txt", 
              n.chains = 3, n.thin = 1, n.iter = 6000, n.burnin = 3000)
  print(out, dig = 2)
  # mcmcplot(out)
  
  #'  Hold on to model output
  coug_cam_road_output <- out
  save(coug_cam_road_output, file = "./Output/coug_cam_road_output.RData")
  
  
  
  ####  Cougar Elevation & Road Density Model  ####
  cat("
  model{
  
    #'  Telemetry half of the model
    for(i in 1:n){
    
      #'  Random intercept for each individual
      a_telem[i] ~ dnorm(int_telem, tau_telem)
      
      #'  The number of observed telemetry locations in a grid cell arises from a 
      #'  categorical distribution based on the probability that a given cell is  
      #'  used by that individual and the total number of its telemetry locations
      M[i,1:ngrid] ~ dmulti(pi[i,1:ngrid], R[i])
    
    
      for(j in 1:ngrid){
      
        #'  Estimate the site- and individual-specific values of lambda
        log(lam_telem[i,j]) = a_telem[i] + b_elev*telev[j] + b_road*troad[j]
        
        #'  Calculating pi based on lambda estimates
        pi[i,j] = lam_telem[i,j]/sum(lam_telem[i,1:ngrid])
      
      }
    }
    
    #'  Camera trap half of the model
    for(k in 1:ncam){
      
      #'  Random effect for each camera
      a_cam[k] ~ dnorm(int_cam, tau_cam)
      
      #'  The number of observed detections at a camera site arises from a 
      #'  Poisson distribution based on the intensity or mean number of independent
      #'  detection events for a given species during the study period 
      y[k] ~ dpois(lam_cam[k])
      
      #'  Estimate the camera-site specific values of lambda 
      log(lam_cam[k]) <- a_cam[k] + b_elev*celev[k] + b_road*croad[k]
    }
    
    
    #'  Priors
    int_telem ~ dnorm(0, 0.1)
    tau_telem ~ dunif(0, 10)
    int_cam ~ dnorm(0, 0.1)
    tau_cam ~ dunif(0, 4)
    b_elev ~ dnorm(0, 0.1)
    b_road ~ dnorm(0, 0.1)
    
    
    #'  Derived parameters

    #'  Averaged across telemetered animals
    #'  Expected number of total locations per grid cell
    for(j in 1:ngrid) {
      mu_lam[j] <- mean(lam_telem[,j])
    }
    
  }
  
  ", fill = TRUE, file = "combo.txt")
  
  
  #'  Arguments for jags
  data <- list(M = coug_telem, R = sumcoug, telev = as.vector(mcp_cov$zDEM),
               celev = as.vector(coug_covs$zDEM), troad = as.vector(mcp_cov$zRoads), 
               croad = as.vector(coug_covs$zRoads), ngrid = ngrid, n = ncoug,
               ncam = ncam, y = coug_cams)
  
  parameters = c('a_telem', 'a_cam', 'int_telem', 'tau_telem', 'int_cam', 'tau_cam',
                 'mu_lam', 'lam_cam', 'b_elev', 'b_road')
  
  inits = function() {list(b_elev = rnorm(1), b_road = rnorm(1))}
  
  
  #'  Call to jags
  out <- jags(data, inits, parameters, "combo.txt", 
              n.chains = 3, n.thin = 1, n.iter = 6000, n.burnin = 3000)
  print(out, dig = 3)
  # mcmcplot(out)
  
  
  #'  Hold on to model output
  coug_combo_elev_road_output <- out
  save(coug_combo_elev_road_output, file = "./Output/coug_combo_elev_road_output.RData")
  
  #'  Put everything on the probability scale to estimate Probability of Use, which
  #'  represents the probability that a cell is used by an individual (telemetry)
  #'  and the probability that a camera site is used in a grid cell (camera)
  #'  during the study period.
  
  #'  Extract all iterations of telemetry- and camera-based lambdas
  mu_lam <- coug_combo_elev_road_output$BUGSoutput$sims.list$mu_lam
  lam_cam <- coug_combo_elev_road_output$BUGSoutput$sims.list$lam_cam
  
  #'  Loop through each iteration for each grid cell and calculate probability
  #'  Telemetry results
  tel_prob <- matrix(nrow = nrow(mu_lam), ncol = ngrid)
  for(i in 1:nrow(mu_lam)){
    for(j in 1:ngrid){
      tel_prob[i,j] <- 1-exp(-(mu_lam[i,j]))
    }
  }
  #'  Camera results
  det_prob <- matrix(nrow = nrow(lam_cam), ncol = ncam)
  for(i in 1:nrow(lam_cam)){
    for(k in 1:ncam){
      det_prob[i,k] <- 1-exp(-(lam_cam[i,k]))
    }
  }
  #'  Transpose so rows = grid cells & save as a data frame
  tel_prob <- as.data.frame(t(tel_prob))
  det_prob <- as.data.frame(t(det_prob))
  
  #' Save estimated probability of use based on each data source
  #' Pair probability of use with the original grid cell number it was estimated for
  cell <- mcp_coug$cell
  coug_pr_use_elev_road_tel <- cbind(cell, tel_prob)
  cell <- cam_and_covs$cell[!is.na(cam_and_covs$Cougar_Detections)]
  coug_pr_use_elev_road_cam <- cbind(cell, det_prob)
  
  save(coug_pr_use_elev_road_tel, file = "./Output/coug_pr_use_elev_road_tel.RData")
  save(coug_pr_use_elev_road_cam, file = "./Output/coug_pr_use_elev_road_cam.RData")
  
  #'-------------------------------------------
  
  #'  RSF model with telemetry data only
  cat("
  model{
    
    for(i in 1:n){
    
      a_telem[i] ~ dnorm(int_telem, tau_telem)
      M[i,1:ngrid] ~ dmulti(pi[i,1:ngrid], R[i])
      
      for(j in 1:ngrid){
      
        log(lam_telem[i,j]) = a_telem[i] + b_elev*telev[j] + b_road*troad[j]
        pi[i,j] = lam_telem[i,j]/sum(lam_telem[i,1:ngrid])
      
      }
    }
    
    #'  Priors
    int_telem ~ dnorm(0, 0.1)
    tau_telem ~ dunif(0, 10)
    b_elev ~ dnorm(0, 0.1)
    b_road ~ dnorm(0, 0.1)

    
    #'  Derived parameters
    #'  Averaged across telemetered animals
    #'  Expected number of total locations per grid cell
    for(j in 1:ngrid) {
      mu_lam[j] <- mean(lam_telem[,j])
    }
  
  }
  
  ", fill = TRUE, file = "telem.txt")
  
  
  #'  Arguments for jags
  data <- list(M = coug_telem, R = sumcoug, telev = as.vector(mcp_cov$zDEM), 
               troad = as.vector(mcp_cov$zRoads), ngrid = ngrid, n = ncoug)
  
  parameters = c('a_telem', 'int_telem', 'tau_telem', 'mu_lam', 'b_elev', 'b_road')
  
  inits = function() {list(b_elev = rnorm(1), b_road = rnorm(1))}
  
  
  # call to jags
  out <- jags(data, inits, parameters, "telem.txt", 
              n.chains = 3, n.thin = 1, n.iter = 6000, n.burnin = 3000)
  print(out, dig = 3)
  # mcmcplot(out)
  
  #'  Hold on to model output
  coug_telem_elev_road_output <- out
  save(coug_telem_elev_road_output, file = "./Output/coug_telem_elev_road_output.RData")
  
  #'------------------------------------------
  
  #'  Camera data only model
  cat("
    model{
    
    for(k in 1:ncam){
      a_cam[k] ~ dnorm(int_cam, tau_cam)
      y[k] ~ dpois(lam_cam[k])
      log(lam_cam[k]) <- a_cam[k] + b_elev*celev[k] + b_road*croad[k]
    }
    
    #'  Priors
    b_elev ~ dnorm(0, 0.1)
    b_road ~ dnorm(0, 0.1)
    int_cam ~ dnorm(0, 0.1)
    tau_cam ~ dunif(0, 4)
    
  }
  
  ", fill = TRUE, file = "cam.txt")
  
  
  
  #'  Arguments for jags
  data <- list(ncam = ncam, y = coug_cams, celev = as.vector(coug_covs$zDEM), 
               croad = as.vector(coug_covs$zRoads))
  
  parameters = c('a_cam', 'int_cam', 'tau_cam', 'lam_cam', 'b_elev', 'b_road') 
  
  inits = function() {list(b_elev = rnorm(1), b_road = rnorm(1))}
  
  # call to jags
  out <- jags(data, inits, parameters, "cam.txt", 
              n.chains = 3, n.thin = 1, n.iter = 6000, n.burnin = 3000)
  print(out, dig = 2)
  # mcmcplot(out)
  
  #'  Hold on to model output
  coug_cam_elev_road_output <- out
  save(coug_cam_elev_road_output, file = "./Output/coug_cam_elev_road_output.RData")
  
  
  
  ####  Cougar Land Cover Model  ####
  cat("
  model{
  
    #'  Telemetry half of the model
    for(i in 1:n){
    
      #'  Random intercept for each individual
      a_telem[i] ~ dnorm(int_telem, tau_telem)
      
      #'  The number of observed telemetry locations in a grid cell arises from a 
      #'  categorical distribution based on the probability that a given cell is  
      #'  used by that individual and the total number of its telemetry locations
      M[i,1:ngrid] ~ dmulti(pi[i,1:ngrid], R[i])
    
    
      for(j in 1:ngrid){
      
        #'  Estimate the site- and individual-specific values of lambda
        #'  Forest is reference variable for NLCD categorical covariate
        log(lam_telem[i,j]) = a_telem[i] + b_shrub*tshrub[j] + b_crop*tcrop[j] + 
            b_water*twater[j] + b_grass*tgrass[j] + b_other*tother[j]

        
        #'  Calculating pi based on lambda estimates
        pi[i,j] = lam_telem[i,j]/sum(lam_telem[i,1:ngrid])
      
      }
    }
    
    #'  Camera trap half of the model
    for(k in 1:ncam){
      
      #'  Random effect for each camera
      a_cam[k] ~ dnorm(int_cam, tau_cam)
      
      #'  The number of observed detections at a camera site arises from a 
      #'  Poisson distribution based on the intensity or mean number of independent
      #'  detection events for a given species during the study period 
      y[k] ~ dpois(lam_cam[k])
      
      #'  Estimate the camera-site specific values of lambda 
      #'  Forest is reference variable for NLCD categorical covariate
      log(lam_cam[k]) <- a_cam[k] + b_shrub*cshrub[k] + b_crop*ccrop[k]
    }
    
    
    #'  Priors
    int_telem ~ dnorm(0, 0.1)
    tau_telem ~ dunif(0, 10)
    int_cam ~ dnorm(0, 0.1)
    tau_cam ~ dunif(0, 4)
    b_shrub ~ dnorm(0, 0.1)
    b_crop ~ dnorm(0, 0.1)
    b_water ~ dnorm(0, 0.1)
    b_grass ~ dnorm(0, 0.1)
    b_other ~ dnorm(0, 0.1)
    
    
    #'  Derived parameters

    #'  Averaged across telemetered animals
    #'  Expected number of total locations per grid cell
    for(j in 1:ngrid) {
      mu_lam[j] <- mean(lam_telem[,j])
    }
    
  }
  
  ", fill = TRUE, file = "combo.txt")
  
  
  #'  Arguments for jags
  data <- list(M = coug_telem, R = sumcoug, ngrid = ngrid, n = ncoug, ncam = ncam, y = coug_cams,
               tshrub = as.vector(mcp_cov$shrub), tcrop = as.vector(mcp_cov$crops),
               twater = as.vector(mcp_cov$water), tgrass = as.vector(mcp_cov$grass),
               tother = as.vector(mcp_cov$other), cshrub = as.vector(coug_covs$shrub),
               ccrop = as.vector(coug_covs$crops))
  
  parameters = c('a_telem', 'a_cam', 'int_telem', 'tau_telem', 'int_cam', 'tau_cam',
                 'mu_lam', 'lam_cam', 'b_shrub', 'b_crop', 'b_water', 'b_grass', 'b_other')
  
  inits = function() {list(b_shrub = rnorm(1), b_crop = rnorm(1), b_water = rnorm(1),
                           b_grass = rnorm(1), b_other = rnorm(1))}
  
  
  #'  Call to jags
  out <- jags(data, inits, parameters, "combo.txt", 
              n.chains = 3, n.thin = 1, n.iter = 6000, n.burnin = 3000)
  print(out, dig = 3)
  # mcmcplot(out)
  
  
  #'  Hold on to model output
  coug_combo_nlcd_output <- out
  save(coug_combo_nlcd_output, file = "./Output/coug_combo_nlcd_output.RData")
  
  #'  Put everything on the probability scale to estimate Probability of Use, which
  #'  represents the probability that a cell is used by an individual (telemetry)
  #'  and the probability that a camera site is used in a grid cell (camera)
  #'  during the study period.
  
  #'  Extract all iterations of telemetry- and camera-based lambdas
  mu_lam <- coug_combo_nlcd_output$BUGSoutput$sims.list$mu_lam
  lam_cam <- coug_combo_nlcd_output$BUGSoutput$sims.list$lam_cam
  
  #'  Loop through each iteration for each grid cell and calculate probability
  #'  Telemetry results
  tel_prob <- matrix(nrow = nrow(mu_lam), ncol = ngrid)
  for(i in 1:nrow(mu_lam)){
    for(j in 1:ngrid){
      tel_prob[i,j] <- 1-exp(-(mu_lam[i,j]))
    }
  }
  #'  Camera results
  det_prob <- matrix(nrow = nrow(lam_cam), ncol = ncam)
  for(i in 1:nrow(lam_cam)){
    for(k in 1:ncam){
      det_prob[i,k] <- 1-exp(-(lam_cam[i,k]))
    }
  }
  #'  Transpose so rows = grid cells & save as a data frame
  tel_prob <- as.data.frame(t(tel_prob))
  det_prob <- as.data.frame(t(det_prob))
  
  #' Save estimated probability of use based on each data source
  #' Pair probability of use with the original grid cell number it was estimated for
  cell <- mcp_coug$cell
  coug_pr_use_nlcd_tel <- cbind(cell, tel_prob)
  cell <- cam_and_covs$cell[!is.na(cam_and_covs$Cougar_Detections)]
  coug_pr_use_nlcd_cam <- cbind(cell, det_prob)
  
  save(coug_pr_use_nlcd_tel, file = "./Output/coug_pr_use_nlcd_tel.RData")
  save(coug_pr_use_nlcd_cam, file = "./Output/coug_pr_use_nlcd_cam.RData")
  
  #'-------------------------------------------
  
  #'  RSF model with telemetry data only
  cat("
  model{
    
    for(i in 1:n){
    
      a_telem[i] ~ dnorm(int_telem, tau_telem)
      M[i,1:ngrid] ~ dmulti(pi[i,1:ngrid], R[i])
      
      for(j in 1:ngrid){
      
        #'  Forest is reference variable for NLCD categorical covariate
        log(lam_telem[i,j]) = a_telem[i] + b_shrub*tshrub[j] + b_crop*tcrop[j] + 
            b_water*twater[j] + b_grass*tgrass[j] + b_other*tother[j]

        pi[i,j] = lam_telem[i,j]/sum(lam_telem[i,1:ngrid])
      
      }
    }
    
    #'  Priors
    int_telem ~ dnorm(0, 0.1)
    tau_telem ~ dunif(0, 10)
    b_shrub ~ dnorm(0, 0.1)
    b_crop ~ dnorm(0, 0.1)
    b_water ~ dnorm(0, 0.1)
    b_grass ~ dnorm(0, 0.1)
    b_other ~ dnorm(0, 0.1)

    
    #'  Derived parameters
    #'  Averaged across telemetered animals
    #'  Expected number of total locations per grid cell
    for(j in 1:ngrid) {
      mu_lam[j] <- mean(lam_telem[,j])
    }
  
  }
  
  ", fill = TRUE, file = "telem.txt")
  
  
  #'  Arguments for jags
  data <- list(M = coug_telem, R = sumcoug, ngrid = ngrid, n = ncoug, 
               tshrub = as.vector(mcp_cov$shrub), tcrop = as.vector(mcp_cov$crops),
               twater = as.vector(mcp_cov$water), tgrass = as.vector(mcp_cov$grass),
               tother = as.vector(mcp_cov$other))
  
  parameters = c('a_telem', 'int_telem', 'tau_telem', 'mu_lam', 
                 'b_shrub', 'b_crop', 'b_water', 'b_grass', 'b_other')
  
  inits = function() {list(b_shrub = rnorm(1), b_crop = rnorm(1), b_water = rnorm(1),
                           b_grass = rnorm(1), b_other = rnorm(1))}
  
  
  # call to jags
  out <- jags(data, inits, parameters, "telem.txt", 
              n.chains = 3, n.thin = 1, n.iter = 6000, n.burnin = 3000)
  print(out, dig = 3)
  # mcmcplot(out)
  
  #'  Hold on to model output
  coug_telem_nlcd_output <- out
  save(coug_telem_nlcd_output, file = "./Output/coug_telem_nlcd_output.RData")
  
  #'------------------------------------------
  
  #'  Camera data only model
  cat("
    model{
    
    for(k in 1:ncam){
    
      a_cam[k] ~ dnorm(int_cam, tau_cam)
      
      y[k] ~ dpois(lam_cam[k])
            
      #'  Forest is reference variable for NLCD categorical covariate
      log(lam_cam[k]) <- a_cam[k] + b_shrub*cshrub[k] + b_crop*ccrop[k]
    }
    
    #'  Priors
    int_cam ~ dnorm(0, 0.1)
    tau_cam ~ dunif(0, 4)
    b_shrub ~ dnorm(0, 0.1)
    b_crop ~ dnorm(0, 0.1)
    
  }
  
  ", fill = TRUE, file = "cam.txt")
  
  
  
  #'  Arguments for jags
  data <- list(ncam = ncam, y = coug_cams, cshrub = as.vector(coug_covs$shrub),
                ccrop = as.vector(coug_covs$crops))
  
  parameters = c('a_cam', 'int_cam', 'tau_cam', 'lam_cam', 'b_shrub', 'b_crop') 
  
  inits = function() {list(b_shrub = rnorm(1), b_crop = rnorm(1))}
  
  # call to jags
  out <- jags(data, inits, parameters, "cam.txt", 
              n.chains = 3, n.thin = 1, n.iter = 6000, n.burnin = 3000)
  print(out, dig = 2)
  # mcmcplot(out)
  
  #'  Hold on to model output
  coug_cam_nlcd_output <- out
  save(coug_cam_nlcd_output, file = "./Output/coug_cam_nlcd_output.RData")

  
  
  ####  Cougar Elevation, Road Density, & Land Cover Model  ####

  