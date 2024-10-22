library(magrittr)
library(terra)

cat_ext <- terra::ext(c(65, 145, -10, 50)) 
plant_ext <- terra::ext(c(-100, -65, 11, 24)) 

# periods to download
periods <- c(
    "1981-2010",
    "2041-2070/GFDL-ESM4/ssp585",
    "2071-2100/GFDL-ESM4/ssp585"
)
# variable of interest
vars <- c(
    paste0("bio", 1:19),
    # "scd"
    # "gdd5",
    # "pet_penman_mean",
    # "pet_penman_range",
    # "npp"
)

clean_name <- function(x) {
    return(
        tolower(gsub("/", "_", x))
    )
}


dir_base <- "data/CHELSA_data"
url_base <- "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/%s/bio"

for (tp in periods) {
    for (vr in vars) {
        # get the full url for download
        url_period <- sprintf(url_base, tp)
        url_full <- sprintf("%s/CHELSA_%s_%s_V.2.1.tif", url_period, vr, clean_name(tp))
        # NOTE: for download just use the parent dir; then delete original
        dwn_name <- sprintf("%s/%s.tif", dir_base, vr)
        
        # output dir base on the input period
        out_dir <- file.path(dir_base, clean_name(tp))
        if (!dir.exists(out_dir)) 
            dir.create(out_dir, recursive = TRUE)
        
        out_ca <- sprintf("%s/%s_%s.tif", out_dir, vr, "cat")
        out_pl <- sprintf("%s/%s_%s.tif", out_dir, vr, "plant")
        
        # do not download the file if already available
        if (any(!file.exists(c(out_ca, out_pl)))) {
            tryCatch(
                {
                    download.file(url = url_full, destfile = dwn_name, method = "curl")
                },
                error = function(cond) {
                    message("Failed download: ", vr, " of ", tp)
                    next
                }
            )
        }
        
        # read and mask the data
        if (file.exists(out_ca)) {
            cat("The", vr, "for cat species exists!\n")
        } else {
            terra::crop(
                x = terra::rast(dwn_name),
                y = cat_ext,
                filename = out_ca
            )
        }
        
        if (file.exists(out_pl)) {
            cat("The", vr, "for plant species exists!\n")
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
}

