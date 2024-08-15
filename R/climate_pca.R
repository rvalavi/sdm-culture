library(terra)
library(geodata)

world_map <- geodata::world(resolution = 4, path = "data")

species_list <- c("cat", "plant")
periods <- c(
    "1981-2010",
    "2041-2070_gfdl-esm4_ssp585",
    "2071-2100_gfdl-esm4_ssp585"
)


for (species in species_list) {
    for (period in periods) {
        cat("Species:", species, "for", period, "\n")
        
        covar_rast <- terra::rast(
            list.files(
                path = sprintf("data/CHELSA_data/%s", period), 
                pattern = sprintf("_%s.tif$", species), 
                full.names = TRUE
            )[1:19]
        ) %>% 
            terra::mask(world_map)
        
        plot(covar_rast[["bio1"]])
        plot(world_map, add = TRUE)
        
        # add one to Aisa dataset to avoid log of zero
        if (species == "cat") {
            covar_rast[["bio14"]] <- covar_rast[["bio14"]] + 1
        }
        # log transform all the rainfall data
        covar_rast[[paste0("bio", 12:19)]] <- log(covar_rast[[paste0("bio", 12:19)]])
        
        # PCA model ---------------------------------------------------------------
        if (period == "1981-2010") {
            tm <- Sys.time()
            set.seed(3010)
            pca <- terra::prcomp(
                x = covar_rast, 
                center = TRUE, 
                scale. = TRUE,
                maxcell = 1e7
            )
            print(Sys.time() - tm)
            
            saveRDS(pca, file = sprintf("data/CHELSA_data/PCA/%s_pca.rds", species))
            
            pca_var <- pca$sdev^2 / sum(pca$sdev^2)
            num_pcs <- seq_along(pca_var)
            plot(
                num_pcs,
                pca_var,
                type = "b",
                xlab = "PC",
                ylab = "Explained variance"
            )
            
            selected_pcs <- num_pcs[cumsum(pca_var) <= 0.995]
            print(
                cumsum(pca_var)
            )
            print(
                selected_pcs
            )
            
        }
        
        # PCA prediction (all times) ----------------------------------------------
        pc_rast <- predict(covar_rast, pca, index = selected_pcs)
        plot(pc_rast)
        
        out_path <- paste0("data/CHELSA_data/PCA/", period)
        if (!dir.exists(out_path)) {
            dir.create(out_path, recursive = TRUE)
        }
        
        terra::writeRaster(
            pc_rast,
            filename = sprintf("%s/pc%s_%s.tif", out_path, selected_pcs, species)
        )
    }
}

