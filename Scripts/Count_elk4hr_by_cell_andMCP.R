#####
### Taylor Ganz
### November 18, 2020
### SEFS 521 Group Project

###############################################################################
#
# This script prepares telemetry data for a model integrating camera trap and collar data
# It also claculates the fix rate and filters it to 4 hours 
#
###############################################################################

library(tidyverse)
library(lubridate)
library(sf)
library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(amt)
library(adehabitatHR)

options(scipen = 999) # turn off scientific notation
options(digits = 15) # set digits to 15 to ensure GPS coordinates aren't truncated

rm(list = ls())


##############################################################################
#--------------------------Bring in your data-------------------------------
NE_grid <- raster("NE_grid_4k.img", layer = 'grid', crs = 2855)
NE_SA <- readOGR("NE/Neast.shp") %>%
  spTransform(crs(NE_grid))

start <- lubridate::ymd_hms('2018-12-01 00:00:00', tz = 'Etc/GMT-8') 
end <- lubridate::ymd_hms('2019-03-31 23:59:59', tz = 'Etc/GMT-8') 

elk <- read.csv("elk_skinny 2020-11-17.csv")%>%
  as.data.frame() %>%
  na.omit() %>%
  mutate(ts = lubridate::ymd_hms(Floordt, tz = "Etc/GMT-8")) %>% 
  filter(VEC_Height <= 2000 & VEC_Height >= 0) %>%
  dplyr::select(ID, CollarID, Latitude, Longitude, ts) %>%
  filter(ts >= start & ts <= end)
#-----------------------------------------------------------------------------
##############################################################################



#############################################################################
#----Thin telemetry data to a 4 hr fix rate, regardless of schedule --------

#start by getting set up for the loops
elkids <- unique(elk$ID) 
ne <- length(elkids)
collars <- unique(elk$CollarID)


# I shouldn't have to do this in a loop but I'm not getting id to work from the amt package 
for(i in 1:ne) {

elk1 <- filter(elk, ID == elkids[i])
elk_amt <- mk_track(elk1, .x=Longitude, .y=Latitude, .t=ts, crs=sp::CRS("+init=epsg:4326")) %>%
    track_resample(rate = hours(4), tolerance = minutes(1)) %>%
    base::as.data.frame() %>%
    transmute(Longitude = x_,
           Latitude = y_,
           ts = lubridate::ymd_hms(t_, tz = "Etc/GMT-8"))

elk_amt$ID <- elkids[i]
elk_amt$CollarID <- collars[i]


if (i == 1) {elk_4hr <- elk_amt} else {elk_4hr <- rbind(elk_4hr, elk_amt)} #reconstruct the dataframe
}
 
elk_4hr <- elk_4hr %>% # just organize it a little nicer
  relocate(ID, .before = Longitude) %>%
  relocate(CollarID, .after = ID)

table(elk_4hr$ID) # see how many fixes we have per animal
# Do we want to consider eliminating those (3696EA17) with really low level of fixes?

#-----------------------------------------------------------------------------------
################################################################################


################################################################
#                                                              # 
#                 THIS STEP IS NOT NECESSARY                   #
#                                                              # 
#   A check to manually determine and examine the fix rate     #
#   You can compare to before the thinning by using the 'elk'  #      
#   data frame instead of 'elk_4hr' in lines 101 and 115       #
#   I also used this to validate the thinning method above     #
#                                                              # 
#--------------------------------------------------------------

for (j in 1:ne){
  who <- elkids[j] # chose your animal
  cow <- filter(elk_4hr, ID == who)

  # now go through and figure out the time between fixes for each animal
  cow$fix_rate[1] <- NA # first value should be NA since we don't have a reference
  nr <- nrow(cow)  
  
    for (k in 2:nr){
    int <- lubridate::interval(cow$ts[k-1], cow$ts[k])
    cow$fix_rate[k] <- time_length(int, 'hour')
    }
  
  if (j == 1) {fr <- cow$fix_rate} else {fr <- c(fr, cow$fix_rate)} # make a list
  #Chug throught and plot the fix rates by elk to make sure all >= 4
  plot(cow$fix_rate, main = who)
}
elk_4hr$fix_rate <- fr # add the list to the elk data frame


# END UNECESSARTY STEP
#---------------------------------------------------------------#
#################################################################



###################################################################
#                                                                 #
#               NOW WE GET TO COUNT THE DATA                      #
#        Make the data  spatial and do some visualizing           #
#------------------------------------------------------------------

elk_sf <- st_as_sf(elk_4hr, coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs(NE_grid))

# Check that it looks ok 
plot(NE_grid, axes = TRUE)
plot(NE_SA, add = TRUE)
plot(elk_sf$geometry, add = TRUE, col = as.factor(elk_sf$ID), pch = 3)

# Now start a plot to visualize the counts as they come through
#NE_na <- NE_grid
#values(NE_na) <- NA
#plot(NE_na, axes=TRUE, legend = FALSE) # baseline, everything set to zero

###=== Count the elk locations by grid cell ====####
#set up the first column of the new data frame
cell <- values(NE_grid)

dev.off()
# iterate by elk
for (i in 1:ne){
 (who <- elkids[i]) # chose your animal
 cow <- filter(elk_sf, ID == who) 
 cow_sp <- as_Spatial(cow)  #convert to sp, required for the rasterize function

 cowrast <- rasterize(cow_sp, NE_grid, field = cow_sp$ID, fun = 'count') # This is where the magic happens
 
 ############################################################################
 ### Some code to do a littl visusalizing if you want - non essential to analysis
 if (i==1) {plot(NE_SA, axes=TRUE)} else {}
 plot(cowrast, add=TRUE, legend = FALSE)
 plot(cow$geometry, add = TRUE, pch = 20, cex = .3) # and see the points
 ###########################################################################
 
 cowrast[is.na(cowrast)] = 0 # overwrite the NA with a zero
 cow_id <- values(cowrast) #store the values in a string
 
 if(i == 1){
  elk_cell_count <- as.data.frame(cbind(cell, cow_id)) #start to build out the data frame
  names(elk_cell_count)[names(elk_cell_count) == "cow_id"] <- who  # rename the new col with animal ID
      } else {
  elk_cell_count <- cbind(elk_cell_count, cow_id) #tack on the new record
  names(elk_cell_count)[names(elk_cell_count) == "cow_id"] <- who # name it right
  }

}

#you'll get a lot of warnings refering to line 153 but it still works
View(elk_cell_count)  

#-------------------- End elk counts by cell --------------------------------#
##############################################################################



############################################################################
#------------------- Now make the MCP where all the elk are -------------#
elk_sp <- as_Spatial(elk_sf) #required format
elk_mcp <- mcp(elk_sp, percent = 100) #make the MCP

mcp_line = as(elk_mcp, "SpatialLinesDataFrame") # we have to deal with the boudnary seperately

mcprast <- rasterize(elk_mcp, NE_grid) #rasterize inside the poly
line_rast <- rasterize(mcp_line, NE_grid) #rasterize the MCP line
mcprast[is.na(mcprast)] = 0 # overwrite the NA with a zero
line_rast[is.na(line_rast)] = 0
full_mcp <- (mcprast + line_rast) # now add these rasters together
full_mcp <- reclassify(full_mcp, c(1.9, 2.1, 1)) # and reclassify to deal with raster values of 2

# make sure it looks good
plot(full_mcp)
plot(elk_mcp, add = TRUE) # add it to our map before

plot(elk_sf$geometry, add = TRUE, col = as.factor(elk_sf$ID), pch = 20, cex = .3)
lines(NE_SA)


elk_cell_count$MCP <- values(full_mcp) # add a new col to the data frame for the MCP

View(elk_cell_count)


#### JUST NOTE THAT THIS WILL CHANGE HOW THE DATA ARE INDEXED IN THE MODEL
elk_cell_count <- elk_cell_count %>% # this is important, we remove all the lines that are not part of the original grid 
  filter(cell != 'NA')

View(elk_cell_count)
write.csv(elk_cell_count, paste0('Elk_4hrFix_Cell_Count ', Sys.Date(), '.csv'))

