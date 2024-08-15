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

covar_rast <- terra::rast(
    list.files(
        path = "data/CHELSA_data/PCA/1981-2010/", 
        pattern = "_cat.tif$",
        full.names = TRUE
    )
) %>% 
    terra::subset(paste0("bio", c(1, 5, 12, 15))) %>% 
    terra::mask(world_map)

plot(covar_rast)

plot(covar_rast[[1]])
plot(world_map, add = TRUE)

# get and clean species data ----------------------------------------------
occ_xy <- readRDS("data/species_data_prionailurus_bengalensis.rds") %>% 
    dplyr::select(longitude, latitude, value) %>% 
    setNames(c("x", "y", "occ")) %>% 
    rm_duplicates(r = covar_rast, column = "occ")

table(occ_xy$occ)
points(occ_xy$x, occ_xy$y)
points(occ_xy$x[occ_xy$occ == 1], occ_xy$y[occ_xy$occ == 1], col = "red")

# extract values ----------------------------------------------------------
species_data <- extract_value(covar_rast, occ_xy, drop_na = TRUE)
# head(model_data)
table(species_data$occ)


# select only modelling columns
model_data <- dplyr::select(species_data, -x, -y) #%>% 
str(model_data)
anyNA(model_data)


# spatial cv --------------------------------------------------------------
data_sf <- sf::st_as_sf(species_data, coords = c("x", "y"), crs = 4326)

# sac <- blockCV::cv_spatial_autocor(r = covar_rast)
# sac$range

set.seed(3010)
scv <- blockCV::cv_spatial(
    x = data_sf,
    column = "occ",
    r = covar_rast,
    k = 5,
    # size = 210000,
    selection = "random", 
    progress = TRUE, 
    max_pixel = 3e6
) 

blockCV::cv_similarity(scv, data_sf, covar_rast)

# cv_plot(cv = scv, x = data_sf)

# model evaluation --------------------------------------------------------
folds <- scv$folds_list

AUCs <- c()

for(k in seq_len(length(folds))){
    train_set <- unlist(folds[[k]][1]) 
    test_set <- unlist(folds[[k]][2])
    
    mod <- ensemble(
        x = model_data[train_set, ],
        y = "occ", 
        fold_ids = scv$folds_ids,
        models = c("GLM", "GAM", "GBM", "RF", "Maxent")
    )
    
    preds <- predict(mod, model_data[test_set, ], type = "response")
    AUCs[k] <- precrec::auc(evalmod(scores = pred, labels = model_data$occ[test_set]))[1,4]
}

print(AUCs)
mean(AUCs)
sd(AUCs)

#
# final model fitting -----------------------------------------------------
# fitting the with spatial CV model tuning
tm <- Sys.time()
model <- ensemble(
    x = model_data,
    y = "occ", 
    fold_ids = scv$folds_ids,
    models = c("GLM", "GAM", "GBM", "RF", "Maxent")
)
Sys.time() - tm

print(model)

# predicting rasters ------------------------------------------------------
terra::terraOptions(steps = 30)

covar_rast <- terra::aggregate(covar_rast, fact = 10)

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
    na.rm = TRUE,
    filename = "outputs/cat/pred_current.tif"
)
Sys.time() - tm

plot(pred_current)


# projections -------------------------------------------------------------
f1_rast <- terra::rast(
    list.files(
        path = "data/CHELSA_data/PCA/2041-2070_gfdl-esm4_ssp585/",
        pattern = "_cat.tif$", 
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
    na.rm = TRUE,
    filename = "outputs/cat/pred_f1.tif"
)
Sys.time() - tm

plot(pred_f1)



f2_rast <- terra::rast(
    list.files(
        path = "data/CHELSA_data/PCA/2071-2100_gfdl-esm4_ssp585/",
        pattern = "_cat.tif$", 
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
    na.rm = TRUE,
    filename = "outputs/cat/pred_f2.tif"
)
Sys.time() - tm

plot(pred_f2)

