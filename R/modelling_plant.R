library(tidyverse)
library(terra)
library(blockCV)
library(sf)

# read models and other functions
source("R/models.R")
source("R/helper_functions.R")

#
# get climate data --------------------------------------------------------
world_map <- geodata::world(resolution = 4, path = "data")
the_ext <- terra::ext(c(-100, -65, 11, 24))

covar_rast <- terra::rast(
    list.files(
        path = "data/CHELSA_data/PCA/1981-2010/", 
        pattern = "_plant.tif$",
        full.names = TRUE
    )
) %>% 
    terra::crop(the_ext)

plot(covar_rast[[1]])
plot(world_map, add = TRUE)

# get and clean species data ----------------------------------------------
occ_xy <- readRDS("data/species_data_zamia_polymorpha.rds") %>% 
    dplyr::select(longitude, latitude, value) %>% 
    setNames(c("x", "y", "occ")) %>% 
    rm_duplicates(r = covar_rast, column = "occ")

table(occ_xy$occ)

# extract values ----------------------------------------------------------
model_data <- extract_value(covar_rast, occ_xy, drop_na = TRUE)
head(model_data)
table(model_data$occ)

# spatial cv --------------------------------------------------------------
data_sf <- sf::st_as_sf(model_data, coords = c("x", "y"), crs = 4326)

sac <- blockCV::cv_spatial_autocor(r = covar_rast)
sac$range

set.seed(3010)
scv <- blockCV::cv_spatial(
    x = data_sf,
    column = "occ",
    r = covar_rast,
    k = 5,
    size = 210000,
    selection = "random", 
    progress = TRUE, 
    max_pixel = 3e6
) 

# cv_plot(cv = scv, x = data_sf)

# model evaluation --------------------------------------------------------



# final model fitting -----------------------------------------------------
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
    models = c("GLM", "GAM", "GBM", "RF", "Maxent")
)
Sys.time() - tm

print(model)

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
    # filename = "outputs/plant/pred_current.tif"
)
Sys.time() - tm

plot(pred_current)


# projections -------------------------------------------------------------

f1_rast <- terra::rast(
    list.files(
        path = "data/CHELSA_data/PCA/2041-2070_gfdl-esm4_ssp585/",
        pattern = "_plant.tif$", 
        full.names = TRUE
    )
)

tm <- Sys.time()
pred_f1 <- terra::predict(
    object = f1_rast,
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
    # filename = "outputs/plant/pred_f1.tif"
)
Sys.time() - tm

plot(pred_f1)



f2_rast <- terra::rast(
    list.files(
        path = "data/CHELSA_data/PCA/2071-2100_gfdl-esm4_ssp585/",
        pattern = "_plant.tif$", 
        full.names = TRUE
    )
)

tm <- Sys.time()
pred_f2 <- terra::predict(
    object = f2_rast,
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
    # filename = "outputs/plant/pred_f2.tif"
)
Sys.time() - tm

plot(pred_f2)

