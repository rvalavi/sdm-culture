library(tidyverse)
library(spocc)
library(rgbif)
library(terra)

sf_bbox <- function(x) {
  return(
    terra::ext(x) %>%
      as.polygons() %>% 
      sf::st_as_sf()
  )
}


# Download Data From GBIF -------------------------------------------------
## Indian leopard cat
## Prionailurus bengalensis
## Distinct from Prionailurus javanensis since 2017. Remove any of these

cat_ext <- c(50, 150, -10, 50)

cat_data <- occ(query = "Prionailurus bengalensis",
                from = "gbif",
                limit = 10000,
                geometry = sf_bbox(cat_ext))

cat_background <- occ(query = "Felis",
                      from = "gbif",
                      limit = 50000,
                      geometry = sf_bbox(cat_ext))

## Sinkhole Cycad
## Zamia prasina
## Previously used incorrectly for Z. decumbens. Remove any of these
## Synonym: Z. polymorpha

# plant_ext <- c(-100, -80, 12, 24)
plant_ext <- c(-100, -65, 6, 24)

plant_data <- occ(query = "Zamia prasina",
                  from = "gbif",
                  limit = 10000,
                  geometry = sf_bbox(plant_ext))

plant_data2 <- occ(query = "Zamia polymorpha",
                   from = "gbif",
                   limit = 10000,
                   geometry = sf_bbox(plant_ext))

plant_background <- occ(query = "Cycadopsida",
                        from = "gbif",
                        limit = 50000,
                        geometry = sf_bbox(plant_ext))


# Clean Data --------------------------------------------------------------
## Indian leopard cat
cat_df <- cat_data$gbif$data[[1]] %>%
  as_tibble() %>%
  filter(name %in% c("Prionailurus bengalensis (Kerr, 1792)",
                     "Prionailurus bengalensis bengalensis",
                     "Prionailurus bengalensis euptilurus (Elliot, 1871)",
                     "Felis bengalensis Kerr, 1792",
                     "Prionailurus bengalensis chinensis (Gray, 1837)",
                     "Prionailurus bengalensis borneoensis Brongersma, 1936",
                     "Prionailurus bengalensis heaneyi Groves, 1997",
                     "Prionailurus bengalensis javanensis (Desmarest, 1816)",
                     "Prionailurus bengalensis sumatranus (Horsfield, 1821)",
                     "Prionailurus bengalensis trevelyani Pocock, 1939",
                     "Prionailurus bengalensis horsfieldii (Gray, 1842)"),
         !is.na(longitude),
         !is.na(latitude),
         year >= 1981,
         coordinateUncertaintyInMeters <= 5000 | is.na(coordinateUncertaintyInMeters))

cat_bg <- cat_background$gbif$data[[1]] %>%
  as_tibble() %>%
  bind_rows(cat_df) %>%
  filter(!is.na(longitude),
         !is.na(latitude),
         year >= 1981,
         coordinateUncertaintyInMeters <= 5000 | is.na(coordinateUncertaintyInMeters))


plot(cat_bg$latitude ~ cat_bg$longitude,
     xlim = c(60, 140),
     ylim = c(-10, 50))

points(cat_df$latitude ~ cat_df$longitude,
     xlim = c(60,140),
     ylim = c(-10, 50),
     col = "red")

cat_points <- bind_rows(cat_df %>%
                          select(longitude,
                                 latitude) %>%
                          mutate(species = "Prionailurus bengalensis",
                                 value = 1),
                        cat_bg %>%
                          select(longitude,
                                 latitude) %>%
                          mutate(species = "Prionailurus bengalensis",
                                 value = 0)) %>%
  relocate(species, 
           .before = longitude)

## Sinkhole cycad

plant_df <- bind_rows(plant_data$gbif$data[[1]],
                      plant_data2$gbif$data[[1]]) %>%
  as_tibble() %>%
  filter(!is.na(longitude),
         !is.na(latitude),
         year >= 1981,
         coordinateUncertaintyInMeters <= 5000 | is.na(coordinateUncertaintyInMeters))
  
plant_bg <- plant_background$gbif$data[[1]] %>%
  as_tibble() %>%
  filter(!is.na(longitude),
         !is.na(latitude),
         year >= 1981,
         coordinateUncertaintyInMeters <= 5000 | is.na(coordinateUncertaintyInMeters))
  
plot(plant_bg$latitude ~ plant_bg$longitude)

points(plant_df$latitude ~ plant_df$longitude,
     # xlim = c(-100, -80),
     # ylim = c(12, 24),
     col = "red")

plant_points <- bind_rows(plant_df %>%
                          select(longitude,
                                 latitude) %>%
                          mutate(species = "Zamia polymorpha",
                                 value = 1),
                        plant_bg %>%
                          select(longitude,
                                 latitude) %>%
                          mutate(species = "Zamia polymorpha",
                                 value = 0)) %>%
  relocate(species, 
           .before = longitude)


# Write To File -----------------------------------------------------------
saveRDS(cat_points, "data/species_data_prionailurus_bengalensis.rds")

saveRDS(plant_points, "data/species_data_zamia_polymorpha.rds")

