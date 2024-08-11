library(tidyverse)
library(sf)
library(tmap)

data("World")

cat_points <- readRDS("data/species_data_prionailurus_bengalensis.rds")
nrow(cat_points)
sum(cat_points$value)

plot(cat_points$longitude, cat_points$latitude)
plot(World$geometry, add = TRUE)
cat_points %>% 
    filter(value == 1) %>% 
    dplyr::select(longitude, latitude) %>% 
    points(col = "red")
abline(v = 65)
abline(v = 145)


pts <- sf::st_as_sf(cat_points, coords = c("longitude", "latitude"), crs = 4326)

tmap_mode("plot")


tm_shape(pts) +
    tm_dots(size = 0.2, alpha = 0.5) +
    tm_shape(World) +
    tm_borders("gray", lwd = .5) +
    tm_shape(pts[pts$value == 1, ]) +
    tm_dots(col = "red", size = 0.2, alpha = 0.5)



plant_points <- readRDS("data/species_data_zamia_polymorpha.rds")
nrow(plant_points)
sum(plant_points$value)


plot(
    plant_points$longitude,
    plant_points$latitude
)
plot(World$geometry, add = TRUE)
plant_points %>% 
    filter(value == 1) %>% 
    dplyr::select(longitude, latitude) %>% 
    points(col = "red")

