# . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ .
# SINGULARITY
# . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ .

# Set a working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(terra) # Processing of spatial data
library(parallel)
library(foreach)

# Import GIS data
grid <- vect("../data/grid/elrs10000sqkm.shp")
countries <- vect("../data/World_Countries/World_Countries.shp")
CHELSA_present <- rast(list.files("../data/present_variables/", pattern = "\\.tif$", full.names = TRUE))
CHELSA_LGM_original <- rast(list.files("../data/glacial_variables/", pattern = "\\.tif$", full.names = TRUE))

altitude <- rast("../data/wc2.1_10m_elev.tif")

# ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ .
# 1. PREDICTORS EDITION
# ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ .

# Delete zero cell values
# CHELSA_LGM_nozero <- clamp(CHELSA_LGM_original, lower = 0)

# Change CHELSA LGM units
CHELSA_LGM <- CHELSA_LGM_original / 10
  object_names <- names(CHELSA_LGM)
  new_object_names <- sub("^CHELSA_PMIP_CCSM4_", "PMIP_", object_names)
  names(CHELSA_LGM) <- new_object_names
CHELSA_LGM_temperatureK <- CHELSA_LGM[[c(1:11)]]
CHELSA_LGM_temperatureC <- CHELSA_LGM_temperatureK - 273
CHELSA_LGM_precipitation <- CHELSA_LGM[[c(12:19)]]
CHELSA_LGM_modified <- c(CHELSA_LGM_temperatureC, CHELSA_LGM_precipitation)

# Mask land data an aggregate CHELSA climate and altitude into reference grid
grid_land <- mask(grid, countries)


dir.create("../CHELSA_LGM_modified")
for (l in 1:nlyr(CHELSA_LGM_modified)) {
    writeRaster(CHELSA_LGM_modified[[l]], paste0("../CHELSA_LGM_modified/", names(CHELSA_LGM_modified)[l], ".tif"))
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Alternative: load CHELSA_modified layers:
CHELSA_LGM_modified <- rast(list.files("../data/CHELSA_LGM_modified/", pattern = "\\.tif$", full.names = TRUE)) 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

CHELSA_present_grid <- as.data.frame(matrix(ncol = nlyr(CHELSA_present) + 1, nrow = grid_cells))
time0 <- Sys.time()
for (i in 1:length(grid_land$ID)) {
  CHELSA_present_grid[i,] <- as.vector(extract(CHELSA_present, grid_land[i], fun = mean, method = "simple"))
  
}
time1 <- Sys.time()- time0
time1

CHELSA_present_grid$ID_cell <- grid_land$ID
present_climate <- CHELSA_present_grid[, c(21, 2:20)]
  colnames(present_climate) <- c("ID_cells", names(CHELSA_present))


write.csv(present_climate, "../CHELSA_PRESENT_GRID.csv")


grid_cells <- length(grid_land$ID)
CHELSA_LGM_grid <- as.data.frame(matrix(ncol = nlyr(CHELSA_LGM_modified) + 1, nrow = grid_cells))

time0 <- Sys.time()
for (i in 1:length(grid_land$ID)) {
  CHELSA_LGM_grid[i,] <- as.vector(extract(CHELSA_LGM_modified, grid_land[i], fun = mean, method = "simple"))
  
}
time1 <- Sys.time()- time0
time
  CHELSA_LGM_grid$ID_cell <- grid_land$ID
  colnames(CHELSA_LGM_grid)<- c("X", names(CHELSA_LGM_modified), "ID_cells")
  LGM_climate <- CHELSA_LGM_grid[, c(21, 2:20)]
write.csv(LGM_climate, "../CHELSA_PRESENT_LGM_GRID.csv")




altitude_grid <- as.data.frame(matrix(ncol = nlyr(altitude) + 1, nrow = grid_cells))

time0 <- Sys.time()
for (i in 1:length(grid_land$ID)) {
  altitude_grid[i,] <- as.vector(extract(altitude, grid_land[i], fun = mean, method = "simple"))
  
}
time1 <- Sys.time()- time0
time

  altitude_grid$ID_dells <-  grid_land$ID
  alitude_modified <- altitude_grid[, c(3,2)]
write.csv(alitude_modified, "../ALTITUDE_GRID.csv")




# Aggregate CHELSA climate and altitude into reference grid

### Parallelization Test ###

# Detect number of cores in computer
n_cores <- parallel::detectCores() - 1
# Define cluster with number of cores - 1 to avoid overloading computer
my_cluster <- parallel::makeCluster(
    n_cores,
    type = "PSOCK"
)

print(my_cluster)

doParallel::registerDoParallel(cl = my_cluster, cores = n_cores)
parallel::clusterCall(my_cluster, function(x) .libPaths(x), .libPaths())

# CHELSA_present <- rast(list.files("../data/present_variables/", pattern = "\\.tif$", full.names = TRUE))

grid_cells <- length(grid_land$ID)

CHELSA_present_grid <- foreach(
    i = 1:2,
    .combine = "rbind",
    .packages = "terra"
    ) %dopar% {
      library(terra)
      CHELSA_present <- terra::rast(list.files("../data/present_variables/", pattern = "\\.tif$", full.names = TRUE))
      grid <- terra::vect("../data/grid/elrs10000sqkm.shp")
      grid_land <- terra::mask(grid, countries)
      terra::extract(CHELSA_present, grid_land[i], fun = mean, method = "simple")
      
      }

time0 <- Sys.time()
extract(CHELSA_present, grid_land[i], fun = mean, method = "simple")
time1 <- Sys.time()- time0
time1


save.image("preliminar_result.RData")

dir.create("./CHELSA_present_grid")
for (l in 1:nlyr(CHELSA_present_grid)) {
  writeRaster(CHELSA_present_grid[[l]], paste0("./CHELSA_present_grid/", names(CHELSA_present_grid)[l], ".tif"))
}

CHELSA_LGM_grid <- list()
for (l in 1:nlyr(CHELSA_LGM_modified)) {
    print(paste("Procesando capa", i, "de", nlyr(CHELSA_LGM_modified), "..."))
    layer <- CHELSA_LGM_modified[[l]]
    grid_rows <- max(grid$FILA)
    CHELSA_LGM_grid[[names(CHELSA_LGM_modified[[l]])]] <- rast(ext(CHELSA_LGM_modified[[l]]), resolution = res(CHELSA_LGM_modified[[l]]))
    CHELSALGM_grid[[names(CHELSA_LGM_modified[[l]])]] <- foreach(
        i = 1:grid_rows,
        .combine = "c",
        .packages = c("terra")
        ) %dopar% {
            grid_subset <- grid[grid$FILA == i]
            extract(layer, grid_subset, mean)
        }
}    

altitude_grid <- rast(ext(altitude), resolution = res(altitude))
    altitude_grid <- foreach(
        i = 1:grid_rows,
        .combine = "c",
        .packages = c("terra")
        ) %dopar% {
            grid_subset <- grid[grid$FILA == i]
            extract(altitude, grid_subset, mean)
        }



parallel::stopCluster(cl = my_cluster)

### End Test ###


# Export modified rasters
dir.create("./CHELSA_present_grid")
for (l in 1:nlyr(CHELSA_present_grid)) {
    writeRaster(CHELSA_present_grid[[l]], paste0("./CHELSA_present_grid/", names(CHELSA_present_grid)[l], ".tif"))
}

dir.create("./CHELSA_LGM_grid")
for (l in 1:nlyr(CCHELSA_LGM_grid)) {
    writeRaster(CHELSA_LGM_grid[[l]], paste0("./CHELSA_LGM_grid/", names(CHELSA_LGM_grid)[l], ".tif"))
}

writeRaster(altitude_grid, "altitude_grig.tiff")