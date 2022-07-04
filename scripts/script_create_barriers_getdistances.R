# Script to generate 50 different scenarios of % forest cover
# Generates distance matrices among sites considering forest areas as forbidden terrain

#loading libraries

library(here)
library(sf)
library(rmapshaper)
library(raster)
library(ggplot2)
library(ggpubr)
library(igraph)
library(patchwork)
library(leaflet)
library(fasterize)
library(gdistance)
library(dplyr)

#reading and verifying data structure

# Loading sampled points and their coordinates
data<-read.csv2(here("src", "distance_polygon.csv"), header=TRUE)

# Creating sf data frame with point id and coordinates
sites_coords<-data.frame(id=data$localidade, utmx=data$utmx_point, utmy=data$utmy_point)
sites_coords<-st_as_sf(sites_coords, coords=c("utmx", "utmy"))

# Loading polygons
shape_original<-st_read(here("src/shape", "Vegetacao.shp"))
# Transforming shape file into an sf dataframe
shape<-st_as_sf(shape_original)

#Setting crs system to sites coordinates
sites_coords<-sites_coords %>% st_set_crs(st_crs(shape))

#Setting amount of increases/decreases of polygons of original shape
increases<-seq(from=-171, to=270, by=9)

# Simplifying polygons for better performance
shape<-ms_simplify(shape, keep=0.05, keep_shapes=TRUE)

# Creating different scenarios varying % of forest cover
shape_buffer<-lapply(increases, function(x, shape){return(st_buffer(shape, x))}, shape=shape)
# Cleaning empty geometries
shape_buffer<-lapply(shape_buffer, function(x){d<-x; d<-d[!st_is_empty(d),,drop=FALSE]; return(d)})
# Calculating percentage forest cover for each scenario
shape_landcover<-lapply(shape_buffer, function(x){box<-st_as_sfc(st_bbox(x)); return(sum(st_area(x))/st_area(box))})

# Creating 500x500 grid rasters for each scenario of forest cover
shapes_rasters<-lapply(shape_buffer, function(x, res){
  a<-st_as_sf(st_union(x))
  r<-raster(extent(a), nrow=res, ncol=res);
  r<-fasterize(a, r);
  return(r)
}, res=500)

# Inverting raster cell values so that forest areas have infinite resistance
rasters_inverted<-lapply(shapes_rasters, function(x){
  #make all vegetation = -999
  y<-x
  y[is.na(y)] <- -999
  #this turns all vegetations to 0, e.g. they become forbidden terrain (infinite resistance) because the cost function uses 1/value
  y[y>-999] <- 0
  #assign unit cost to all grid cells in land
  y[y==-999] <- 1
  return(y)
})

clean_units <- function(x){
  attr(x,"units") <- NULL
  class(x) <- setdiff(class(x),"units")
  x
}

# Creating distance matrices with the assumption that points can only reach each other going around forest areas
D_matrices<-lapply(rasters_inverted, function(x, points_coords){
  r_tr <- transition(x, transitionFunction = mean, directions = 16)
  r_tr <- geoCorrection(r_tr, type = "c")
  dist <- costDistance(r_tr,fromCoords = as(as_Spatial(points_coords), "SpatialPoints"),
                       toCoords = as(as_Spatial(points_coords), "SpatialPoints"))
  return(clean_units(dist))
}, points_coords=sites_coords)