library(geodata)

# plant_ext <- terra::ext(c(-100, -65, 11, 24)) 
# get the extent from PCA to exactly match
plant_ext <- terra::ext(
    terra::rast("data/CHELSA_data/PCA/1981-2010/pc1_plant.tif")
) 

out_dir <- "data/Soil"
soil_vars <- c("sand", "silt", "clay", "bdod")

for (sv in soil_vars) {
    
    out_pl <- sprintf("%s/%s_%s.tif", out_dir, sv, "plant")
    
    if (file.exists(out_pl)) {
        cat("The", sv, "for plant species exists!\n")
    } else {
        global_rast <- geodata::soil_world(var = sv, depth = 5, path = out_dir)
        
        terra::crop(
            x = global_rast,
            y = plant_ext,
            filename = out_pl
        )
    
        # remove the main file
        unlink(global_rast)
        print(global_rast)
    }
}

