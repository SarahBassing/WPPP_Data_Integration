  #'  Data Integration RSF Model
  #'  SEFS 521- Final Presentation Example
  #'  December 11, 2020
  #'  --------------------------------------------
  #'  Combines observations of cougar and elk, respectively, in the WPPP Northeast 
  #'  study area from Dec 1, 2018 - Mar 31, 2019. Observations collected from 55  
  #'  camera traps, 44 elk GPS collars & 24 cougar GPS collars under one habitat 
  #'  selection model. Script runs a series of models that includes an integrated 
  #'  model, a telemetry-only model, and a camera-only model to evaluate how
  #'  results change across data sources and when they are combined.
  #'  
  #'  Original model written by Dr. Beth Gardner
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
  #'  y = total number of species-specific detections in a grid cell- number of
  #'      unique animals observed in frame during an independent detection event
  #'  
  #'  PARAMETERS
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
  
  
  ####  Elk Elevation & Road Density Model  ####
  
  #'  Integrated model!
  cat("
    model{
    
      #'  Telemetry sub-model
      for(i in 1:n){
        
        #'  The number of observed telemetry locations in a grid cell arises from a 
        #'  categorical distribution based on the probability that a given cell is  
        #'  used by that individual and the total number of its telemetry locations
        M[i,1:ngrid] ~ dmulti(pi[i,1:ngrid], R[i])
      
      
        for(j in 1:ngrid){
        
          #'  Calculating pi based on lambda estimates
          pi[i,j] = lam_telem[i,j]/sum(lam_telem[i,1:ngrid])
          
          #'  Estimate the site- and individual-specific values of lambda
          log(lam_telem[i,j]) = b_elev*telev[j] + b_road*troad[j]
        
        }
      }
      
      #'  Camera trap sub-model
      for(k in 1:ncam){
        
        #'  The number of observed detections at a camera site arises from a 
        #'  Poisson distribution based on the intensity or mean number of independent
        #'  detection events for a given species during the study period 
        y[k] ~ dpois(lam_cam[k])
        
        #'  Estimate the camera-site specific values of lambda 
        log(lam_cam[k]) <- a_cam[k] + b_elev*celev[k] + b_road*croad[k]
        
        #'  Random effect for each camera
        a_cam[k] ~ dnorm(int_cam, tau_cam)
      }
      
      
      #'  Priors
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
  
  parameters = c('a_cam', 'int_cam', 'tau_cam', 'mu_lam', 'lam_cam', 'b_elev', 'b_road')
  
  inits = function() {list(b_elev = rnorm(1), b_road = rnorm(1))}
  
  
  #'  Call to jags
  out <- jags(data, inits, parameters, "combo.txt", 
              n.chains = 3, n.thin = 1, n.iter = 10000, n.burnin = 5000)
  print(out, dig = 3)
  # mcmcplot(out)
  
  
  #'  Hold on to model output
  elk_combo_elev_road_output <- out

  #'-------------------------------------------
  
  #'  RSF model with telemetry data only
  cat("
    model{
      
      for(i in 1:n){
      
        M[i,1:ngrid] ~ dmulti(pi[i,1:ngrid], R[i])
        
        for(j in 1:ngrid){
        
          log(lam_telem[i,j]) = b_elev*telev[j] + b_road*troad[j]
          
          pi[i,j] = lam_telem[i,j]/sum(lam_telem[i,1:ngrid])
        
        }
      }
      
      #'  Priors
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
  
  parameters = c('mu_lam', 'b_elev', 'b_road')
  
  inits = function() {list(b_elev = rnorm(1), b_road = rnorm(1))}
  
  
  # call to jags
  out <- jags(data, inits, parameters, "telem.txt", 
              n.chains = 3, n.thin = 1, n.iter = 10000, n.burnin = 5000)
  print(out, dig = 3)
  # mcmcplot(out)
  
  #'  Hold on to model output
  elk_telem_elev_road_output <- out
 
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
              n.chains = 3, n.thin = 1, n.iter = 10000, n.burnin = 5000)
  print(out, dig = 2)
  # mcmcplot(out)
  
  #'  Hold on to model output
  elk_cam_elev_road_output <- out

