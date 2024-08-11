library(tidyverse)
library(terra)

# read models code
source("R/models.R")

extract_value <- function(r, x, drop_na = TRUE) {
  vals <- cbind(
    x,
    terra::extract(r, x[, c("x", "y")], ID = FALSE)
  )
  
  return(
    if (drop_na) {
      tidyr::drop_na(vals)
    } else {
      vals
    }
  )
}

# get species data --------------------------------------------------------
occ_xy <- readRDS("data/species_data_zamia_polymorpha.rds") %>% 
  dplyr::select(longitude, latitude, value) %>% 
  setNames(c("x", "y", "occ"))

# get climate data --------------------------------------------------------
covar_rast <- terra::rast(
  list.files("data/CHELSA_data/PCA/1981-2010/", pattern = "_plant.tif$", full.names = TRUE)
)
plot(covar_rast[[1:4]])

model_data <- extract_value(covar_rast, occ_xy)
head(model_data)
table(model_data$occ)

# spatial cv --------------------------------------------------------------
library(blockCV)
library(sf)

data_sf <- sf::st_as_sf(model_data, coords = c("x", "y"), crs = 4326)

# # cv_spatial_autocor(r = rr)
# sac <- cv_spatial_autocor(x = data_sf, column = "prod")
# sac$range

set.seed(3010)
scv <- cv_spatial(
  x = data_sf,
  column = "occ",
  r = covar_rast,
  k = 5,
  # size = 360000,
  selection = "random", 
  progress = TRUE, 
  max_pixel = 2e6
) 

# cv_plot(cv = scv, x = data_sf)

# model fitting -----------------------------------------------------------
# select only modelling columns
training_data <- dplyr::select(model_data, -x, -y) #%>% 
str(training_data)
anyNA(training_data)

# fitting the with spatial CV model tuning
tm <- Sys.time()
model <- ensemble(
  x = training_data,
  y = "occ", 
  fold_ids = scv$folds_ids, 
  models = c("GLM", "GAM", "GBM", "RF")
)
Sys.time() - tm

print(model)

# model evaluation --------------------------------------------------------


#
# predicting rasters ------------------------------------------------------
# predict to raster layers
tm <- Sys.time()
pred_current <- terra::predict(
  object = covar_rast,
  model = model,
  type = "response",
  cpkgs = c(
    "ranger",
    "dismo",
    "gbm",
    "mgcv",
    "glmnet"
  ),
  na.rm = TRUE
  # filename = "outputs/pred_current.tif"
)
Sys.time() - tm

plot(pred_current)

