#### Create a covariate layer box for the NE ####
# Taylor Ganz #
# November 16, 2020

library(lubridate)
library(tidyverse)
library(sf)
library(sp)
library(rgdal)
library(rgeos)
library(raster)

rm(list = ls())

options(scipen = 999) # turn off scientific notation
options(digits = 15) # set digits to 15 to ensure GPS coordinates aren't truncated

sa_proj <- projection("+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs ")

NE_SA <- readOGR("NE/Neast.shp") %>%
  spTransform(sa_proj)

elk <- read.csv("telem/elk_skinny 2020-11-17.csv")
wtd <- read.csv("telem/wtd_skinny 2020-10-28.csv")

ung <- rbind(elk, wtd) %>%
  as.data.frame() %>%
  filter(VEC_Height <= 2000 & VEC_Height >= 0) %>%
  mutate(ObsDayTimePST = lubridate::mdy_hms(ObservationDateTimePST)) %>%
  dplyr::select(ID, CollarID, Longitude, Latitude, VEC_DOP, VEC_Height, 
                VEC_FixType, ObsDayTimePST, TransmissionDateTimePST, Species) %>%
  na.omit() #%>%

ung_sf <- st_as_sf(ung, coords = c("Longitude", "Latitude"), crs = 4326) %>% #make it spatial
  st_transform(sa_proj)

plot(ung_sf$geometry, axes = TRUE, col = 2, pch = 20, cex = .3)
plot(NE_SA, add = TRUE)


min(ung$Latitude) #find the extreme values
  # 47.37379296
plot(ung$Latitude)

true_ung <- filter(ung, Latitude > 47.8) #remove the outliers

dim(ung) - dim(true_ung) # lost 6 rows, seems reasonable

true_ung_sf <- st_as_sf(true_ung, coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(sa_proj)

st_crs(NE_SA) == st_crs(true_ung_sf)
st_crs(NE_SA) == st_crs(ung_sf) # should get TRUE

plot(true_ung_sf$geometry, axes = TRUE, col = 2, pch = 20, cex = .3)
plot(NE_SA, add = TRUE)

box_sfc <- st_as_sfc(st_bbox(true_ung_sf$geometry)) # create a polygon of the bbox
st_crs(box_sfc) == st_crs(NE_SA)
plot(box_sfc, add = TRUE)


box_poly <- st_polygonize(st_bbox(true_ung_sf$geometry))

buff_box <- st_buffer(true_ung_sf$geometry, 10000)
buff_sfc <- st_as_sfc(buff_box)

#### manually create bbox ####
(xmin <- min(true_ung$Longitude) - .5)
(xmax <- max(true_ung$Longitude) + .5)
(ymin <- min(true_ung$Latitude) - .2)
(ymax <- max(true_ung$Latitude) + .2) 

box <- matrix(c(xmin, ymin,
                xmax, ymin,
                xmax, ymax,
                xmin, ymax),  ncol = 2, byrow = TRUE)

poly_box = Polygon(box)
NE_box = SpatialPolygons(list(Polygons(list(poly_box), ID = "NE_box"), proj4string = st_crs(NE_SA))) 

%>%
  
  spTransform(st_crs(NE_SA))

plot(NE_box, axes = TRUE)
plot(true_ung_sf$geometry, add=TRUE)
plot(NE_SA, add = TRUE)

     

#now make a bigger buffered box     
(bxmin <- min(true_ung$Longitude) - .25)
(bxmax <- max(true_ung$Longitude) + .25)
(bymin <- min(true_ung$Latitude) - .05)
(bymax <- max(true_ung$Latitude) + .05) 
     
buff <- matrix(c(bxmin, bymin,
                     bxmax, bymin,
                     bxmax, bymax,
                     bxmin, bymax),  ncol = 2, byrow = TRUE)
     
buff_box = Polygon(buff)
NE_buff = SpatialPolygons(list(Polygons(list(buff_box), ID = "NE_buff")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

plot(NE_buff, axes = TRUE, border = 'blue')
#plot(NE_box, add=TRUE)

plot(true_ung_sf$geometry, add=TRUE, pch = 19, cex = .5, col = as.factor(true_ung_sf$Species))
plot(NE_SA, add = TRUE)

NE_buff_sp <- as(NE_buff, "SpatialPolygonsDataFrame")
writeOGR(NE_buff_sp, dsn = '.', layer = 'NE_covariate_area', driver = "ESRI Shapefile")
