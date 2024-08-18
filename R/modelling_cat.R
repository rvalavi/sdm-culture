# remotes::install_github("rvalavi/myspatial")
library(tidyverse)
library(terra)
library(blockCV)
library(sf)

# read models and other functions
source("R/models.R")
source("R/helper_functions.R")

#
# get climate data --------------------------------------------------------
exc_countries <- c("MNG", "KAZ") # "KGZ"
world_map <- geodata::world(resolution = 4, path = "data")
world_map <- world_map[!world_map$GID_0 %in% exc_countries]

the_ext <- terra::ext(c(65, 145, -10, 48))

covar_rast <- terra::rast(
    c(
        list.files(
            path = "data/CHELSA_data/PCA/1981-2010/", 
            pattern = "_cat.tif$",
            full.names = TRUE
        ),
        "data/Topo/MSTPI_cat.tif"
    )
) %>% 
    terra::crop(the_ext) %>% 
    terra::mask(world_map)

# plot(covar_rast)
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

#
# spatial cv --------------------------------------------------------------
data_sf <- sf::st_as_sf(species_data, coords = c("x", "y"), crs = 4326)

set.seed(3010)
scv <- blockCV::cv_spatial(
    x = data_sf,
    column = "occ",
    r = covar_rast,
    k = 5,
    size = 400000,
    selection = "random", 
    progress = TRUE, 
    max_pixel = 3e6
) 

blockCV::cv_similarity(scv, data_sf, covar_rast)

# cv_plot(cv = scv, x = data_sf)

# model evaluation --------------------------------------------------------
folds <- scv$folds_list

AUCs <- c()
PRCs <- c()
Boyce <- c()
models <- list()

for(k in seq_len(length(folds))) {
    train_set <- unlist(folds[[k]][1]) 
    test_set <- unlist(folds[[k]][2])
    
    mod <- ensemble(
        x = model_data[train_set, ],
        y = "occ", 
        models = c("GLM", "GAM", "GBM", "RF", "Maxent")
    )
    
    preds <- predict(mod, model_data[test_set, -1], type = "response")
    AUCs[k] <- calc_auc(preds, model_data$occ[test_set])
    PRCs[k] <- calc_prc(preds, model_data$occ[test_set])
    Boyce[k] <- calc_boyce(preds, preds[which(model_data[test_set, ]$occ == 1)])
    models[[k]] <- mod
}

print(AUCs)
print(mean(AUCs))
print(sd(AUCs))

# check the response curves of with CI of CV
myspatial::ggResponse2(
    models = models, 
    covariates = model_data[, -1], 
    type = "response"
)

#
# final model fitting -----------------------------------------------------
# fitting the with model tuning
tm <- Sys.time()
model <- ensemble(
    x = model_data,
    y = "occ", 
    fold_ids = NULL, # random-cv tuning
    models = c("GLM", "GAM", "GBM", "RF", "Maxent")
)
Sys.time() - tm

print(model)

# check the response curves
myspatial::ggResponse(
    models = model, 
    covariates = model_data[, -1], 
    type = "response"
)

#
# predicting rasters ------------------------------------------------------
terra::terraOptions(steps = 30)

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
    filename = "outputs/cat/pred_current.tif",
    wopt = list(names = "current"),
    overwrite = TRUE
)
Sys.time() - tm

plot(pred_current, range = c(0, 1))


# projections -------------------------------------------------------------
f1_rast <- terra::rast(
    c(
        list.files(
            path = "data/CHELSA_data/PCA/2041-2070_gfdl-esm4_ssp585/",
            pattern = "_cat.tif$", 
            full.names = TRUE
        ),
        "data/Topo/MSTPI_cat.tif"
    )
) %>% 
    terra::crop(the_ext) %>% 
    terra::mask(world_map)

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
    filename = "outputs/cat/pred_f1.tif",
    wopt = list(names = "proj-2041-2070"),
    overwrite = TRUE
)
Sys.time() - tm

plot(pred_f1, range = c(0, 1))



f2_rast <- terra::rast(
    c(
        list.files(
            path = "data/CHELSA_data/PCA/2071-2100_gfdl-esm4_ssp585/",
            pattern = "_cat.tif$", 
            full.names = TRUE
        ),
        "data/Topo/MSTPI_cat.tif"
    )
) %>% 
    terra::crop(the_ext) %>% 
    terra::mask(world_map)

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
    filename = "outputs/cat/pred_f2.tif",
    wopt = list(names = "proj-2071-2100"),
    overwrite = TRUE
)
Sys.time() - tm

plot(pred_f2, range = c(0, 1))


terra::panel(
    x = c(pred_current, pred_f1, pred_f2), 
    range = c(0, 1),
    nr = 2
)

