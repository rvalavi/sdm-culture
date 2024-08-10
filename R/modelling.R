library(tidyverse)
library(terra)

# read models code
source("R/models.R")

# get climate data --------------------------------------------------------
cover_rast <- terra::rast(
  list.files("data/CHELSA_data/1981-2010/", pattern = "_cat.tif$", full.names = TRUE)
)
plot(cover_rast[[1:4]])


# spatial cv --------------------------------------------------------------
library(blockCV)
library(sf)

data_sf <- sf::st_as_sf(model_data_xy, coords = c("x", "y"), crs = 4326)

# # cv_spatial_autocor(r = rr)
# sac <- cv_spatial_autocor(x = data_sf, column = "prod")
# sac$range

set.seed(101)
scv <- cv_spatial(
  x = data_sf,
  r = rr,
  k = 5,
  # size = 360000,
  selection = "random", 
  progress = TRUE, 
  max_pixel = 2e6
) 

cv_plot(cv = scv, x = data_sf)

# model fitting -----------------------------------------------------------
# select only modelling columns
model_data <- dplyr::select(model_data_xy, -x, -y) #%>% 
# mutate(prod = log(prod))
str(model_data)


# fitting the with spatial CV model tuning
tm <- Sys.time()
scv_mod <- ensemble(
  x = model_data,
  y = "prod", 
  fold_ids = scv$folds_ids
)
Sys.time() - tm


# model evaluation --------------------------------------------------------


#
# predicting rasters ------------------------------------------------------
# predict to raster layers
tm <- Sys.time()
pred_current <- predict(
  object = rr,
  model = rcv_mod,
  type = "response",
  cpkgs = c(
    "ranger",
    "dismo",
    "gbm",
    "mgcv",
    "glmnet"
  ),
  na.rm = TRUE,
  filename = "outputs/pred_current.tif"
)
Sys.time() - tm

plot(pred_current)

