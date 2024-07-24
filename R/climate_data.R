library(magrittr)
library(terra)


months <- 1:12
years <- 2003:2017
vars <- c("pr", "tas", "tasmax", "tasmin")

out_path <- "data/CHELSA_data/"
url_base <- "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/monthly"

# australia <- terra::vect("data/bioregions/Australia_border.shp")
plot(australia)

for(vr in vars){
  for(yr in years){
    for(mn in months){
      
      name <- sprintf("CHELSA_%s_%02d_%s_V.2.1.tif", vr, mn, yr)
      out_name <- paste0(out_path, name)
      # get the full url for download
      url_full <- file.path(url_base, vr, name)
      
      tryCatch(
        {
          download.file(url = url_full, destfile = out_name, method = "curl")
        },
        error = function(cond) {
          message("This data failed:", vy, yr, mn)
          next
        }
      )
      
      # read and mask the data
      r <- terra::rast(out_name) %>% 
        terra::crop(australia) %>% 
        terra::mask(
          mask = australia,
          filename = file.path(out_path, vr, name)
        )
      
      # remove the main file
      unlink(out_name)
      print(out_name)
    }
  }
}

