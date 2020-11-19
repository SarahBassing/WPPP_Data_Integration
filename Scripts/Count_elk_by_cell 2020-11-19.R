#####
### Taylor Ganz
### November 18, 2020
### SEFS 521 Group Project

### Prep Elk telemetry data for integrated camera trap - telemetry model of species use

library(tidyverse)
library(lubridate)
library(sf)
library(sp)
library(rgdal)
library(rgeos)
library(raster)

options(scipen = 999) # turn off scientific notation
options(digits = 15) # set digits to 15 to ensure GPS coordinates aren't truncated

rm(list = ls())
#dev.off()

NE_grid <- raster("NE_grid_4k.img", layer = 'grid', crs = 2855)
NE_SA <- readOGR("NE/Neast.shp") %>%
  spTransform(crs(NE_grid))

#set the date range for the analysis

start <- lubridate::ymd_hms('2018-12-01 00:00:00', tz = 'Etc/GMT+8') 
end <- lubridate::ymd_hms('2019-03-31 23:59:59', tz = 'Etc/GMT+8') 
  
##########################################################################
#
#           BIG NOTE: need to to deal with fix rate                     #
#
##########################################################################

elk <- read.csv("elk_skinny 2020-11-17.csv")%>%
  as.data.frame() %>%
  na.omit() %>%
  mutate(ts = lubridate::ymd_hms(Floordt, tz = "Etc/GMT+8")) %>% 
  filter(VEC_Height <= 2000 & VEC_Height >= 0) %>%
  dplyr::select(ID, CollarID, Latitude, Longitude, ts, Floordt) %>%
  filter(ts >= start & ts <= end)


# Make the elk locations spatial
elk_sf <- st_as_sf(elk, coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs(NE_grid))

# Check that it looks ok 
plot(NE_grid, axes = TRUE)
plot(NE_SA, add = TRUE)
plot(elk_sf$geometry, add = TRUE, col = as.factor(elk_sf$ID), pch = 3)


elkids <- unique(elk$ID) 
nelk <- length(elkids)

#set up the first column of the data frame
cell <- values(NE_grid)

# now iterate by elk
for (i in 1:nelk){
 (who <- elkids[i]) # chose your animal
 cow <- filter(elk_sf, ID == who)
 
 cow_sp <- as_Spatial(cow)  #convert to sp, required for the rasterize function
 #will throw some warnings but it works
 
 cowrast <- rasterize(cow_sp, NE_grid, fun = 'count') # This is where the magic happens
 cow_ct <-cowrast$ID #the other metrics are all duplicates, we only need one
 cow_ct[is.na(cow_ct)] = 0 # overwrite the NA with a zero
 cow_id <- values(cow_ct) #store the values in a string
 
 ############################################################################
 ### Some code to do a little visusalizing if you want
 ### Just set i to what you want (or set who)
 #
 # Start by making the background easier to look at
 # NE_zg <- NE_grid
 # values(NE_zg) <- 0 
 #
 # plot(NE_zg, axes=TRUE, legend = FALSE) # baseline, everything set to zero
 # plot(cow_ct, add=TRUE) # see the values for the cells 
 # plot(NE_SA, add = TRUE)
 # plot(cow$geometry, add = TRUE, pch = 20, cex = .3) # see the points
 ###########################################################################
 
 if(i == 1){
  elk_cell_count <- as.data.frame(cbind(cell, cow_id)) #start to build out the data frame
  names(elk_cell_count)[names(elk_cell_count) == "cow_id"] <- who  # rename the new col with animal ID
 } else {
   elk_cell_count <- cbind(elk_cell_count, cow_id) #tack on the new record
   names(elk_cell_count)[names(elk_cell_count) == "cow_id"] <- who # name it right
 }

}
 
#you'll get a lot of warnings refering to line 69 but it seems to work ok

write.csv(elk_cell_count, paste0('Elk_Cell_Count ', Sys.Date(), '.csv'))
