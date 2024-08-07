library(magrittr)
library(terra)

cat_ext <- terra::ext(c(65, 145, -10, 50)) 
plant_ext <- terra::ext(c(-100, -80, 12, 24)) 

# periods to download
periods <- c("1981-2010", "2041-2070/GFDL-ESM4/ssp585")
# variable of interest
vars <- c(
    paste0("bio", 1:19),
    "gdd10",
    
)

clean_name <- function(x) {
    return(
        tolower(gsub("/", "_", x))
    )
}

clean_name("1981-2010")
clean_name("2041-2070/GFDL-ESM4/ssp585")

"CHELSA_bio11_1981-2010_V.2.1.tif"
"CHELSA_bio11_2041-2070_gfdl-esm4_ssp585_V.2.1.tif" 

out_path <- "data/CHELSA_data"
url_base <- "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/%s/bio"


for (tp in periods) {
    for (vr in vars) {
        
        # get the full url for download
        url_full <- sprintf("%s/CHELSA_%s_%s_V.2.1.tif", url_base, vr, tp)
        # NOTE: for download just use the parent dir; then delete original
        dwn_name <- sprintf("%s/%s.tif", out_path, vr)
        
        tryCatch(
            {
                download.file(url = url_full, destfile = dwn_name, method = "curl")
            },
            error = function(cond) {
                message("Failed download: ", vr, " of ", tp)
                next
            }
        )
        
        # output dir base on the input period
        out_dir <- file.path(out_path, clean_name(tp))
        if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
        out_name <- file.path(out_dir, vr)
        
        # read and mask the data
        terra::crop(
            x = terra::rast(dwn_name), 
            y = cat_ext,
            filename = sprintf("%s/%s_%s.tif", out_dir, vr, "cat")
        )
        
        terra::crop(
            x = terra::rast(dwn_name), 
            y = plant_ext,
            filename = sprintf("%s/%s_%s.tif", out_dir, vr, "plant")
        )
        
        # remove the main file
        unlink(out_name)
        print(out_name)
    }
}
