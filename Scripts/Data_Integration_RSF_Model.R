  #'  Data Integration RSF Model
  #'  SEFS 521
  #'  November 2020
  #'  --------------------------------------------
  #'  Combines observations of cougar and elk, respectively, in the WPPP Northeast 
  #'  study area from Dec 1, 2018 - Mar 31, 2019. Observations collected from 55  
  #'  camera traps, 44 elk GPS collars & 24 cougar GPS collars under one habitat 
  #'  selection model.
  #'  
  #'  Model written by Dr. Beth Gardner
  #'  --------------------------------------------

  #'  Load libraries
  library(rjags)
  library(mcmcplots)
  library(tidyverse)
  
  #'  Read in data and format for model
  #'  Species specific telemetry data
  elk <- read.csv("./Elk_4hrFix_Cell_Count_withNAs 2020-11-25.csv") %>%
    select(-X)
  coug <- read.csv("./Cougar_4hrFix_Cell_Count 2020-11-28.csv")
  #'  Camera trap detection data
  cams <- read.csv("./Camera_detections.csv") %>%
    select(-X)
  #'  Covariate data
  cov <- read.csv("./Covariates_by_cell.csv")
  
  #'  Remove cells from all data where covariates have NAs
  elk <- elk[!is.na(cov$DEM_val) & !is.na(cov$roads_val),]
  coug <- coug[!is.na(cov$DEM_val) & !is.na(cov$roads_val),]
  cams <- cams[!is.na(cov$DEM_val) & !is.na(cov$roads_val),] # lose cell 677 with 1 cam (no detections tho)
  cov <- cov[!is.na(cov$DEM_val) & !is.na(cov$roads_val),]
  
  #'  Telemetry observation data only
  elk_telem <- as.matrix(select(elk, -c(cell, MCP)))
  elk_telem <- t(elk_telem) # transposing matrix so rows = animal, column = grid cell
  coug_telem <- as.matrix(select(coug, -c(cell, MCP)))
  coug_telem <- t(coug_telem) # transposing matrix so rows = animal, column = grid cell
  
  #'  Species-specific camera detection data
  #'  Remove un-sampled grid cells
  #'  Change to vector instead of data frame
  elk_cams <- select(cams, Elk_Detections) %>%
    filter(!is.na(.))
  elk_cams <- elk_cams[["Elk_Detections"]]
  coug_cams <- select(cams, Cougar_Detections) %>%
    filter(!is.na(.))
  coug_cams <- coug_cams[["Cougar_Detections"]]

  #'  Covariate data
  cov <- cov %>%
    select(-grid_cell_NE) %>%
    #'  Center and scale covariate data
    mutate(
      zDEM = scale(DEM_val, center = TRUE, scale = TRUE),
      zRoads = scale(roads_val, center = TRUE, scale = TRUE)
    )
  #'  Add camera sampling effort & MCP boundary to covariate data
  covs <- cbind(cov, cams$Camera_Sampled, elk$MCP, coug$MCP) %>%
    relocate(zDEM, .after = DEM_val) %>%
    relocate(zRoads, .after = roads_val)
  #'  Rename columns for easier use in model  
  colnames(covs) <- c("DEM", "zDEM", "Roads", "zRoads", "NLCD", "NLCD_names", 
                      "Camera_survey", "Elk_MCP", "Coug_MCP") 
  #'  Covariate data for only grid cells with cameras
  cam_and_covs <- cbind(cams, covs)
  elk_covs <- cam_and_covs[!is.na(cam_and_covs$Elk_Detections),] %>%
    select(-c(cell, Camera_Sampled, Cougar_Detections, Elk_Detections))
  coug_covs <- cam_and_covs[!is.na(cam_and_covs$Cougar_Detections),] %>%
    select(-c(cell, Camera_Sampled, Cougar_Detections, Elk_Detections))
  
  # fake <- rnorm(nrow(covs), 0, 1)
  # fakec <- sample(fake, length(elk_cams), replace = TRUE)
  
  #'  Look over the data & create a few summary statistics
  #'  Keep in mind telemetry data are arranged as individual grid cells (rows) by
  #'  number of individual locations per grid cell (columns)
  dim(elk_telem)
  dim(elk_cams)
  dim(coug_telem)
  dim(coug_cams)
  dim(covs)
  
  head(elk_telem[,1:6])
  head(coug_telem[,1:6])
  
  #'  Sum total number of telemetry locations per individual animal
  #sumelk <- colSums(elk_telem)
  #sumcoug <- colSums(coug_telem)
  sumelk <- rowSums(elk_telem)  # if matrix is transposed
  sumcoug <- rowSums(coug_telem)  # if matrix is transposed
  
  #'  Number of grid cells, cameras, and telemetered animals
  #'  Keep in mind multiple cameras fall within one grid cell so think about 
  #'  whether this should be number of cameras total (55) or number of grids with
  #'  cameras in them (<55)
  ngrid <- nrow(covs)
  #ncam <- sum(cams$Camera_Sampled, na.rm = TRUE) # Number of cameras total
  ncam <- sum(!is.na(cams$Camera_Sampled)) # Number of cells with cameras
  # nelk <- length(elk_telem)
  # ncoug <- length(coug_telem)
  nelk <- nrow(elk_telem)  # if matrix is transposed
  ncoug <- nrow(coug_telem)  # if matrix is transposed

  #'  ---------------------------------------------------
  #'  Helpful definitions for input data and estimated/derived parameters
  #'  DATA
  #'  ngrid = number of grid cells
  #'  ncam = number of camera traps
  #'  n = number of telemetered animals
  #'  R = total number of telemetry locations for an individual animal
  #'  M = number of telemetry locations per individual per grid cell
  #'  y = total number of species-specific detections in a grid cell
  #'  
  #'  PATAMETERS
  #'  alpha = random effect for each telemetered animal
  #'  pi = a categorical probability that each grid cell is used by an individual
  #'       animal; calculated by dividing the site-specific estimates of lambda
  #'       for a given individual by the mean lambda for that individual 
  #'       (averaged across all sites)
  #'  lam = (lambda) the intensity or rate parameter representing the estimated  
  #'        mean number of telemetry locations in a grid cell for an individual 
  #'        animal during the study period
  #'  a0 = random effect for each camera site
  #'  lamc = (lambda) the intensity or rate parameter representing the estimated
  #'         mean number of independent detection events of a given species at
  #'         a camera site during the study period
  #'  b1 = beta coefficient for the effect of a given covariate on the number of 
  #'       animal observations; NOTE: estimating the same b1 for both types 
  #'       of data to connect the two observation processes under one model
  #'  ---------------------------------------------------

  #'  Combined model for both gps and camera data
  cat("
  model{
  
    #'  Telemetry half of the model
    for(i in 1:n){
    
      #'  Random intercept for each individual
      alpha[i] ~ dnorm(int, tau)
      
      #'  The number of observed telemetry locations in a grid cell arises from a 
      #'  categorical distribution based on the probability that a given cell is  
      #'  used by that individual and the total number of its telemetry locations
      M[i,1:ngrid] ~ dmulti(pi[i,1:ngrid], R[i])
    
    
      for(j in 1:ngrid){
      
        #'  Estimate the site- and individual-specific values of lambda
        log(lam[i,j]) = alpha[i] + b1*X[j]
        
        #'  Calculating pi based on lambda estimates
        pi[i,j] = lam[i,j]/sum(lam[i,1:ngrid])
      
      }
    }
    
    #'  Camera trap half of the model
    for(k in 1:ncam){
      
      #'  Random effect for each camera
      a0[k] ~ dnorm(intc, tauc)
      
      #'  The number of observed detections at a camera site arises from a 
      #'  Poisson distribution based on the intensity or mean number of independent
      #'  detection events for a given species during the study period 
      y[k] ~ dpois(lamc[k])
      
      #'  Estimate the camera-site specific values of lambda 
      log(lamc[k]) <- a0[k] + b1*Xc[k]
    }
    
    
    #'  Priors
    int ~ dnorm(0, 0.1)
    tau ~ dunif(0, 10)
    b1 ~ dnorm(0, 0.1)
    
    intc ~ dnorm(0, 0.1)
    tauc ~ dunif(0, 4)
    
    
    #'  Derived parameters
    #'  Back-transform results from log scale to original scale
    
    #'  Telemetry model
    for(i in 1:n){
      for(j in 1:ngrid){
        exp_lam[i,j] <- exp(alpha[i] + b1*X[j])
      }
    }
    #'  Averaged across telemetered animals
    #'  Estimated total number of locations per grid cell
    for(j in 1:ngrid) {
      mu_lam[j] <- mean(exp_lam[,j])
    }
    
    #'  Camera model
    #'  Estimated number of independent detections per grid cell with a camera
    for(k in 1:ncam){
      exp_lamc[k] <- exp(a0[k] + b1*Xc[k])
    }
    
    #'  Put everything on the probability scale to estimate Probability of Use
    #'  In reality, it's the probability of 1 or more telemetry locations/camera
    #'  detections occuring in a given grid cell during the study period.

    for(j in 1:ngrid){
      tel_prob[j] <- 1-exp(-mu_lam[j])
    }
    
    for(k in 1:ncam) {
      det_prob[k] <- 1-exp(-exp_lamc)
    }
    #' This only gives prob. for sites with cameras--- how do we merge with 
    #' telemetry data to make a single 'use' surface?
  }
  
  ", fill = TRUE, file = "combo.txt")
  
  
  
  #'  Arguments for jags
  #'  Elk data   # X = as.vector(covs$zDEM), Xc = as.vector(elk_covs$zDEM),
  data <- list(M = elk_telem, R = sumelk, X = as.vector(covs$zDEM), Xc = as.vector(elk_covs$zDEM), ngrid = ngrid, n = nelk, ncam = ncam, y = elk_cams)
  #data <- list(M = elk_telem, R = sumelk, X = as.vector(fake), Xc = as.vector(fakec), ngrid = ngrid, n = nelk, ncam = ncam, y = elk_cams)
  
  #'  Cougar data
  #data <- list(M = coug_telem, R = sumcoug, X = as.vector(fake), Xc = as.vector(fakec), ngrid = ngrid, n = ncoug, ncam = ncam, y = coug_cams)
  
  parameters = c('alpha', 'a0', 'b1', 'int', 'tau', 'intc', 'tauc', 'tel_prob', 'det_prob')
  
  inits = function() {list(b1 = rnorm(1))}
  
  
  #'  Call to jags
  mod <- jags.model("combo.txt", data, inits, n.chains = 3, n.adapt = 100)
  fit <- coda.samples(mod, parameters, n.iter = 1000)
  summary(fit)
  mcmcplot(fit)
  #plot(fit)
  
  #'  Hold on to model output

  
  
  #################################
  
  #'  RSF with telemetry data only
  cat("
  model{
    
    for(i in 1:n){
    
      alpha[i] ~ dnorm(int, tau)
      M[i,1:ngrid] ~ dmulti(pi[i,1:ngrid], R[i])
      
      for(j in 1:ngrid){
      
        lam[i,j] = alpha[i] + b1*X[j]
        hold[i,j] <- exp(lam[i,j])
        pi[i,j] = hold[i,j]/sum(hold[i,1:ngrid])
      
      }
    }
    
    #'  Priors
    int ~ dnorm(0,.1)
    tau ~ dunif(0,10)
    b1 ~ dnorm(0, .1)
  
  }
  
  ", fill = TRUE, file = "rsf.txt")
  
  
  #'  Arguments for jags
  data <- list(M = elk_telem, R = sumelk, X = as.vector(fake), ngrid = ngrid, n = nelk)
  #data <- list(M = coug_telem, R = sumcoug, X = as.vector(fake), ngrid = ngrid, n = ncoug)
  parameters = c('alpha','b1', 'int','tau')
  
  inits = function() {list(b1=rnorm(1))}
  
  
  # call to jags
  mod <- jags.model("rsf.txt", data, inits, n.chains = 3, n.adapt = 500)
  fit <- coda.samples(mod,parameters, n.iter = 5000)
  summary(fit)
  mcmcplot(fit)
  #plot(fit)
  
  ############################################
  
  #'  Camera data only
  cat("
    model{
    
    for(k in 1:ncam){
      a0[k] ~ dnorm(intc, tauc)
      y[k] ~ dpois(lamc[k])
      log(lamc[k]) <- a0[k] + b1*Xc[k]
    }
    
    #'  Priors
    b1~dnorm(0, .1)
    intc~dnorm(0,.1)
    tauc~dunif(0,4)
    
  }
  
  ", fill = TRUE, file = "cam.txt")
  
  
  
  #'  Arguments for jags
  data = list(Xc = as.vector(fakec), ncam = ncam, y = elk_cams)
  #data = list(Xc = as.vector(fakec), ncam = ncam, y = coug_cams)
  parameters = c('a0', 'b1', 'intc', 'tauc')
  
  inits = function() {list(b1=rnorm(1))}
  
  
  # call to jags
  mod <- jags.model("cam.txt", data, inits, n.chains = 3, n.adapt = 500)
  fit <- coda.samples(mod, parameters, n.iter = 15000)
  summary(fit)
  mcmcplot(fit)
  #plot(fit)
  