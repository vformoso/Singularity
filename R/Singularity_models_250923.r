# . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ .
# SINGULARITY
# . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ .

# Set a working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(terra) # Processing of spatial data
library(fuzzySim)
library(MASS)

# ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ .
# 0. IMPORT DATA
# ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ .

# Import singularity data
amph_sing <- read.table("../data/amph-singularity-data.txt", sep = " ", header = TRUE)
bird_sing <- read.table("../data/bird-singularity-data.txt", sep = " ", header = TRUE)
mamm_sing <- read.table("../data/mamm-singularity-data.txt", sep = " ", header = TRUE)
rept_sing <- read.table("../data/rept-singularity-data.txt", sep = " ", header = TRUE)

# Import GIS data
grid <- vect("../data/grid/elrs10000sqkm.shp")
countries <- vect("../data/World_Countries/World_Countries.shp")
CHELSA_present <- rast(list.files("../data/present_variables/", pattern = "\\.tif$", full.names = TRUE))
CHELSA_LGM_original <- rast(list.files("../data/glacial_variables/", pattern = "\\.tif$", full.names = TRUE))

altitude_original <- rast("../data/wc2.1_10m_elev.tif")

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


dir.create("../CHELSA_LGM_modified")
for (l in 1:nlyr(CHELSA_LGM_modified)) {
  writeRaster(CHELSA_LGM_modified[[l]], paste0("../CHELSA_LGM_modified/", names(CHELSA_LGM_modified)[l], ".tif"))
}


# Mask land data an aggregate CHELSA climate and altitude into reference grid

grid_land <- mask(grid, countries)
grid_cells <- length(grid_land$ID)

CHELSA_present_grid <- as.data.frame(matrix(ncol = nlyr(CHELSA_present) + 1, nrow = grid_cells))

for (i in 1:length(grid_land$ID)) {
  CHELSA_present_grid[i,] <- as.vector(extract(CHELSA_present, grid_land[i], fun = mean, method = "simple"))
  
}
  CHELSA_present_grid$ID_cell <- grid_land$ID
  
present_climate <- CHELSA_present_grid[, c(21, 2:20)]
  colnames(present_climate) <- c("ID_cells", names(CHELSA_present))
  write.csv(present_climate, "../CHELSA_PRESENT_GRID.csv")


CHELSA_LGM_grid <- as.data.frame(matrix(ncol = nlyr(CHELSA_LGM_modified) + 1, nrow = grid_cells))

for (i in 1:length(grid_land$ID)) {
  CHELSA_LGM_grid[i,] <- as.vector(extract(CHELSA_LGM_modified, grid_land[i], fun = mean, method = "simple"))
  
}
  CHELSA_LGM_grid$ID_cell <- grid_land$ID
  colnames(CHELSA_LGM_grid)<- c("X", names(CHELSA_LGM_modified), "ID_cells")

LGM_climate <- CHELSA_LGM_grid[, c(21, 2:20)]
  write.csv(LGM_climate, "../CHELSA_PRESENT_LGM_GRID.csv")


altitude_grid <- as.data.frame(matrix(ncol = nlyr(altitude_original) + 1, nrow = grid_cells))
altitude_max <- as.data.frame(matrix(ncol = nlyr(altitude_original) + 1, nrow = grid_cells))
altitude_min <- as.data.frame(matrix(ncol = nlyr(altitude_original) + 1, nrow = grid_cells))


for (i in 1:length(grid_land$ID)) {
  altitude_grid[i,] <- as.vector(extract(altitude_original, grid_land[i], fun = mean, method = "simple"))
  altitude_max[i,] <- as.vector(extract(altitude_original, grid_land[i], fun = max, method = "simple"))
  altitude_min[i,] <- as.vector(extract(altitude_original, grid_land[i], fun = min, method = "simple"))
  altitude_range <- altitude_max - altitude_min
  
}
  altitude_grid$ID_dells <-  grid_land$ID
  altitude_range$ID_dells <-  grid_land$ID

altitude <- cbind(altitude_grid[ , c(3,2)], altitude_range[ , 2])
  colnames(altitude)<- c("ID_cells", "mean_altitude", "altitude_range")
  write.csv(altitude, "../ALTITUDE_GRID.csv")

  
# ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ .
# 2. PREPARE DATA FOR MODELLING
# ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ .

all_taxa <- list(amph_sing, rept_sing, bird_sing, mamm_sing)
selected_taxa <- menu(c("Amphibians", "Reptiles", "Birds", "Mammals"), title = "Please select a taxa to work with:")
  
taxa_df <- all_taxa[[selected_taxa]]

taxa_climate <- cbind(taxa_df, present_climate[match(taxa_df$sites, present_climate$ID_cell), ])
  taxa_climate <- cbind(taxa_climate, altitude$mean_altitude[match(taxa_df$sites, altitude$ID_cells)], altitude$altitude_range[match(taxa_df$sites, altitude$ID_cells)])
  names_to_change <- colnames(taxa_climate)
  new_names <- c(names_to_change[1:35], "mean_altitude", "range_altitude")
  colnames(taxa_climate) <- new_names
  taxa_climate$mean_altitude[is.nan(taxa_climate$mean_altitude)] <- 0
  taxa_climate$range_altitude[is.nan(taxa_climate$range_altitude)] <- 0
  taxa_climate$dif_temp <- present_climate$bio1[match(taxa_df$sites, present_climate$ID_cell)] - LGM_climate$PMIP_BIO_01[match(taxa_df$sites, LGM_climate$ID_cells)]
  taxa_formodel <- na.omit(taxa_climate)
  

# ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ .
# 3. MODELLING
# ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ .
  
response_variables <- c("mean.sing.value", "a.exp", "b.exp")
selected_response <- menu(c("Singularity", "Intercept", "Slope"), title = "Please select a response variable to work with:")


vars <- colnames(taxa_formodel[ ,c(17:38)]) # Identify climatic variables in the dataframe
  
# Identify weakly correlated variables (R < 0.8) and select the more relevant for the target species, according to VIF
uncor_vars <- corSelect(data = taxa_formodel, sp.cols = response_variables[[selected_response]], var.cols = vars, cor.thresh = 0.8, select = "VIF")

# Select predictor variables for GLM using a stepwise procedure
null_formula <- formula(paste(response_variables[[selected_response]], "~ 1"))
glm_null <- glm(formula = null_formula, data = taxa_formodel, family = gaussian(link = "identity"))
glm_maximal <- formula(paste(response_variables[[selected_response]], "~ ", paste(uncor_vars$selected.vars, collapse = " + ")))
glm_step <- step(glm_null, glm_maximal, data = taxa_formodel, direction = "both", trace = TRUE)
glm_step  

