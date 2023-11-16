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

for (i in 1:length(grid_land$ID)) {
  altitude_grid[i,] <- as.vector(extract(altitude_original, grid_land[i], fun = mean, method = "simple"))
  
}
  altitude_grid$ID_dells <-  grid_land$ID

alitude_ <- altitude_grid[, c(3,2)]
  write.csv(alitude_modified, "../ALTITUDE_GRID.csv")

  
# ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ .
# 2. PREPARE DATA FOR MODELLING
# ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ .

amph_climate <- cbind(amph_sing, present_climate[match(amph_sing$sites, present_climate$ID_cell), ])
  amph_climate <- cbind(amph_climate, altitude$V2[match(amph_sing$sites, altitude$ID_dells)])
  names_to_change <- colnames(amph_climate)
  new_names <- c(names_to_change[1:35], "altitude")
  colnames(amph_climate) <- new_names
  amph_climate$altitude[is.nan(amph_climate$altitude)] <- 0
  amph_climate$dif_temp <- present_climate$bio1[match(amph_sing$sites, present_climate$ID_cell)] - LGM_climate$PMIP_BIO_01[match(amph_sing$sites, LGM_climate$ID_cells)]
  
bird_climate <- cbind(bird_sing, present_climate[match(bird_sing$sites, present_climate$ID_cell), ])
  bird_climate <- cbind(bird_climate, altitude$V2[match(bird_sing$sites, altitude$ID_dells)])
  names_to_change <- colnames(bird_climate)
  new_names <- c(names_to_change[1:35], "altitude")
  colnames(bird_climate) <- new_names
  bird_climate$altitude[is.nan(bird_climate$altitude)] <- 0
  bird_climate$dif_temp <- present_climate$bio1[match(bird_sing$sites, present_climate$ID_cell)] - LGM_climate$PMIP_BIO_01[match(bird_sing$sites, LGM_climate$ID_cells)]
  
rept_climate <- cbind(rept_sing, present_climate[match(rept_sing$sites, present_climate$ID_cell), ])
  rept_climate <- cbind(rept_climate, altitude$V2[match(rept_sing$sites, altitude$ID_dells)])
  names_to_change <- colnames(rept_climate)
  new_names <- c(names_to_change[1:35], "altitude")
  colnames(rept_climate) <- new_names
  rept_climate$altitude[is.nan(rept_climate$altitude)] <- 0
  rept_climate$dif_temp <- present_climate$bio1[match(rept_sing$sites, present_climate$ID_cell)] - LGM_climate$PMIP_BIO_01[match(rept_sing$sites, LGM_climate$ID_cells)] 
  
mamm_climate <- cbind(mamm_sing, present_climate[match(mamm_sing$sites, present_climate$ID_cell), ])
  mamm_climate <- cbind(mamm_climate, altitude$V2[match(mamm_sing$sites, altitude$ID_dells)])
  names_to_change <- colnames(mamm_climate)
  new_names <- c(names_to_change[1:35], "altitude")
  colnames(mamm_climate) <- new_names
  mamm_climate$altitude[is.nan(mamm_climate$altitude)] <- 0
  mamm_climate$dif_temp <- present_climate$bio1[match(mamm_sing$sites, present_climate$ID_cell)] - LGM_climate$PMIP_BIO_01[match(mamm_sing$sites, LGM_climate$ID_cells)]
  
  
# . . . . . . . . . . . . . 
# 2.2. SELECT PREDICTORS
# . . . . . . . . . . . . . 

# Amphibians
amph_subset <- amph_climate[ , c(14, 17:37)]
vars <- colnames(amph_subset[ ,c(2:22)]) # Identify climatic variables in the dataframe
  
# Identify weakly correlated variables (R < 0.8) and select the more relevant for the target species, according to VIF
uncor_vars <- corSelect(data = amph_subset, sp.cols = "mean.sing.value", var.cols = vars, cor.thresh = 0.8, select = "VIF")

# Select predictor variables for GLM using a stepwise procedure
glm_null <- glm(formula = mean.sing.value ~ 1, data = amph_subset, family = gaussian(link = "identity"))
glm_maximal <- formula(paste("mean.sing.value ~", paste(uncor_vars$selected.vars, collapse = " + ")))
glm_step <- step(glm_null, glm_maximal, data = amph_subset, direction = "both", trace = TRUE)
  