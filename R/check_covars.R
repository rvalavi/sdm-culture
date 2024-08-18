library(terra)

# preparing topo data for cat ---------------------------------------------
the_rast <- terra::rast("data/CHELSA_data/PCA/1981-2010/pc1_plant.tif")

r <- terra::resample(
    x = terra::rast("data/Topo/MSTPI_pl.tif"),
    y = the_rast,
    method = "near",
    filename = "data/Topo/MSTPI_plant.tif",
    names = "TPI"
)

plot(
    c(the_rast, r)
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
