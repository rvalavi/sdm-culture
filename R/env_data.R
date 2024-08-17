library(geodata)

plant_ext <- terra::ext(c(-100, -65, 11, 24)) 

out_dir <- "data/soil"
soil_vars <- c("phh2o", "cec", "soc")

for (sv in soil_vars) {
    
    out_ca <- sprintf("%s/%s_%s.tif", out_dir, sv, "cat")
    out_pl <- sprintf("%s/%s_%s.tif", out_dir, sv, "plant")
    
    dwn_name <- geodata::soil_world(var = sv, depth = 5, path = "data/soil")
    
    if (file.exists(out_pl)) {
        cat("The", sv, "for plant species exists!\n")
    } else {
        terra::crop(
            x = terra::rast(dwn_name),
            y = plant_ext,
            filename = out_pl
        )
    }
    
    # remove the main file
    unlink(dwn_name)
    print(dwn_name)
}

