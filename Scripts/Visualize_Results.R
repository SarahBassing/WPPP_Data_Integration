  #'  Visualize Integrated RSF Results
  #'  SEFS 521
  #'  December 2020
  #'  --------------------------------------------
  #'  Plot lambda or probability of use across study area. Overlay animal locations
  #'  and detections to visualize results.
  #'  --------------------------------------------

  #'  Load libraries
  
  #'  Read in and pull out relevant data
  out <- load("./Output/combo_DEM_output.RData")
  
  mu_lam <- combo_DEM_output$BUGSoutput$mean$mu_lam
  mu_lam_sd <- combo_DEM_output$BUGSoutput$sd$mu_lam
  tel_prob <- combo_DEM_output$BUGSoutput$mean$tel_prob
  tel_prob_sd <- combo_DEM_output$BUGSoutput$sd$tel_prob
  exp_lamc <- combo_DEM_output$BUGSoutput$mean$exp_lamc
  exp_lamc_sd <- combo_DEM_output$BUGSoutput$sd$exp_lamc
  det_prob <- combo_DEM_output$BUGSoutput$mean$det_prob
  det_prob_sd <- combo_DEM_output$BUGSoutput$sd$det_prob
  
  #'  Load spatial data