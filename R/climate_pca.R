library(terra)
library(geodata)

world_map <- geodata::world(resolution = 4, path = "data")

covar_rast <- terra::rast(
    list.files(
        path = "data/CHELSA_data/1981-2010/", 
        pattern = "_cat.tif$", 
        full.names = TRUE
    )[1:19]
) %>% 
    terra::mask(world_map)

plot(covar_rast[["bio1"]])
plot(world_map, add = TRUE)

plot(covar_rast[[paste0("bio", 12:19)]])

terra::global(covar_rast[["bio12"]], fun = "min")
plot(log(covar_rast[["bio12"]]))



# PCA model ---------------------------------------------------------------
tm <- Sys.time()
set.seed(3010)
pca <- terra::prcomp(
    x = covar_rast, 
    center = TRUE, 
    scale. = TRUE,
    maxcell = 1e7
)
Sys.time() - tm


pca$center
pca$scale

pca_var <- pca$sdev^2 / sum(pca$sdev^2)
num_pcs <- seq_along(pca_var)
plot(num_pcs, pca_var, type = "b", xlab = "bio", ylab = "Explained variance")


num_pcs[cumsum(pca_var) < 0.99]


# PCA prediction (all times) ----------------------------------------------

pc_rast <- predict(covar_rast, pca, index = 1:9)
plot(pc_rast)

terra::writeRaster(
    pc_rast,
    filename = paste0("data/CHELSA_data/PCA/current/pc", 1:9)
)


