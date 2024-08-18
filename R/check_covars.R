library(terra)

# preparing topo data for cat ---------------------------------------------
cat_rast <- terra::rast("data/CHELSA_data/PCA/1981-2010/pc1_cat.tif")

r <- terra::resample(
    x = terra::rast("data/Topo/MSTPI.tif"),
    y = cat_rast,
    method = "near",
    filename = "data/Topo/MSTPI_cat.tif",
    names = "TPI"
)

plot(
    c(cat_rast, r)
)



# check colinearity -------------------------------------------------------
covars <- terra::rast(
    c(
        list.files(
            path = "data/CHELSA_data/PCA/1981-2010/", 
            pattern = "_cat.tif$",
            full.names = TRUE
        ),
        "data/CHELSA_data/1981-2010/scd_cat.tif",
        "data/Topo/MSTPI_cat.tif"
    )
)

plot(covars)

pairs(covars)
