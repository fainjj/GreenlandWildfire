require(tidyverse)
require(magrittr)
require(rgdal)
require(rgeos)
require(sf)
require(sp)
require(gdalUtils)
require(raster)
require(rasterVis)
require(RColorBrewer)
require(h5)

## Greenland ####
# Get the entire area of Greenland
grl_all <- getData('GADM', level = 0, country = "Greenland")
# Limit this to a polygon around our study area
study_area <- readOGR('bbox_grl.kml')
grl <- crop(grl_all, study_area)
# Also create a larger version to include the Arctic Circle Trail
with_trail <- readOGR('bbox_grl_trail.kml')
grl_trail <- crop(grl_all, with_trail)
#


## Roads ####
# OSM roads/paths
roads <- readOGR('roads.shp')
#


## SMAP ####
# Soil Moisture via SMAP, 36km
sm_files <- list.files(pattern = 'SMAP')
sm <- get_subdatasets(sm_files)
#


## MODIS Burned Area ####
# Burned area for August
ba_files_8 <- list.files(pattern = 'MCD64.*213.*\\.hdf') %>% as.list()
ba_sds_8 <- lapply(ba_files_8, get_subdatasets) %>% lapply(as.list)
# Burned area for September
ba_files_9 <- list.files(pattern = 'MCD64.*244.*\\.hdf') %>% as.list()
ba_sds_9 <- lapply(ba_files_9, get_subdatasets) %>% lapply(as.list)
#

# Reproject Tiles
for(i in seq_len(length(ba_sds_8))){
  for(j in seq_len(length(ba_sds_8[[i]]))){
    gdal_translate(ba_sds_8[[i]][[j]],
                   paste('ba8_', i, j, '.tif', sep = ''))
  }
}
# Merge reprojected tiles
ba8 <- list.files(pattern = "ba8.*\\.tif$") %>%
  as.list() %>%
  lapply(raster)
ba8 <- do.call(merge, ba8)
#

# Reproject tiles
for(i in seq_len(length(ba_sds_9))){
  for(j in seq_len(length(ba_sds_9[[i]]))){
    gdal_translate(ba_sds_9[[i]][[j]],
                   paste('ba9_', i, j, '.tif', sep = ''))
  }
}
# Merge reprojected tiles
ba9 <- list.files(pattern = "ba9.*\\.tif$") %>%
  as.list() %>%
  lapply(raster)
ba9 <- do.call(merge, ba9)
#

# Both Months
ba <- stack(ba8, ba9)
ba <- projectRaster(ba, crs = CRS('+proj=longlat +datum=WGS84 +no_defs
                                   +ellps=WGS84 +towgs84=0,0,0'))
#


## Study Burned Area ####
fire <- unzip('GreenlandWildfireExtent.kml.zip')
fire <- readOGR(fire)
#


## Soil Carbon ####
# Read in soil carbon
soil_c_files <- list.files(pattern = 'OCS.*\\.tif$') %>% as.list()
soil_c <- lapply(soil_c_files, raster)
# Align all rasters
collate <- function(x){extent(x) <- extent(soil_c[[1]]); return(x)}
soil_c <- lapply(soil_c, collate)
soil_c <- brick(soil_c)
# List of standard depths
sd <- c('0–5 cm',	'5–15 cm', '15–30 cm', '30–60 cm', '60–100 cm',	'100–200 cm')
#


## Soil Type ####
# Read in soil type
soil_t <- list.files(pattern = 'TAXN.*\\.tif$') %>% raster()
# Give a description for each soil type code
soil_t <- as.factor(soil_t)
rat <- levels(soil_t)[[1]]
rat[["type"]] <- c('Haplic Chernozems','Acric Ferrasols',
                   'Plinthic Gleysols', 'Plinthic Luvisols',
                   'Humic Nitosols', 'Haplic Xerosols')
levels(soil_t) <- rat
# Plot Soil type
levelplot(soil_t, col.regions=brewer.pal(6, 'Dark2'), xlab="", ylab="",
          main = "Soil Types")
#


## Extraction ####
pull <- function(x){a <- mask(x, fire); return(a)}
punch <- function(x){a <- crop(x, fire); return(a)}
#

# Study area Soil Carbon
fire_c <- soil_c %>%
  as.list() %>%
  lapply(pull) %>%
  lapply(collate) %>%
  lapply(punch) %>%
  brick() # Missing values at sd5. Interpolate from sd4 and sd6.
sd5_interp <- (fire_c[[4]]+fire_c[[6]])/2
fire_c[[5]] <- merge(fire_c[[5]], sd5_interp)
names(fire_c[[5]]) <- str_extract(names(soil_c[[5]]), 'OCS.*m\\.tif$')
#

# Study area Soil Type
fire_t <- soil_t %>% pull() %>% punch()
levelplot(fire_t, col.regions=brewer.pal(6, 'Dark2'), xlab="", ylab="",
          main = "Soil Types")
#

# Greenland Soil Carbon
pull_grl <- function(x){a <- mask(x, grl); return(a)}
punch_grl <- function(x){a <- crop(x, grl); return(a)}
grl_c <- lapply(as.list(soil_c), pull_grl) %>%
  lapply(collate) %>%
  lapply(punch_grl) %>%
  brick()
#

## Fuel ####
fuel <- raster('Global_fuelbeds_map_Tile2.tif')
fuel %<>% pull_grl() %>% punch_grl()
fuel_params <- readxl::read_xlsx('Global_fuelbeds_parameters.xlsx',
                                 'Fuelbeds_metric')
fuel <- as.factor(fuel)
rat <- levels(fuel)[[1]]
rat[["type"]] <- c("Mosaic trees-shrubs/grasses", "Shrubs", "Grasses",
                   "Sparse vegetation", "Shrubs-grasses regularly flooded")
levels(fuel) <- rat
levelplot(fuel, col.regions=brewer.pal(5, 'Dark2'), xlab="", ylab="",
          main = "Fuelbeds")
#

## Carbon Comparison ####
par(mfrow = c(3, 4), bg = '#9ab7c3', mar = c(3,3,2,1))
for(i in seq_len(6)){
  plot(fire_c[[i]],
       breaks = seq(0, 900, length.out = 275),
       legend = FALSE,
       main = paste('Soil Carbon at Depth', sd[i], sep = ' '))
  plot(fire, add = T)
  breaks <- hist(fire_c[[6]], plot=FALSE)$breaks
  hist(fire_c[[i]],
       xlim = c(0,800),
       main = 'Histogram of Values',
       xlab = 'tonnes/ha',
       breaks = seq(0, 850, length.out = 17),
       col = rev(terrain.colors(length(breaks)-1)))
}
#

## Arctic Circle Trail ####
par(mfrow = c(1,1), mar = c(0,0,2,0))
plot(grl_trail, col = 'aliceblue', border = NA)
plot(roads, add = T, col = "blue")
plot(fire, add = T, col = "red", border = NA)
title(main = 'Burned Area and Arctic Circle Trail', line = -0.8)
