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


# ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ .
# 3. COVARIANCE ANALYSIS
# ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ . ~ .

# 3.1. ONLY CLIMATIC VARIABLES

# Create subsets of  grouped-by-variables data
singularity_data <- taxa_formodel[ , c(1:16)]
names_singularity_data <- colnames(singularity_data)
current_climatic_data <- cbind(singularity_data, taxa_formodel[ , grep("bio", names(taxa_formodel))])
current_climatic_vars  <- colnames(current_climatic_data[ , c(17:35)])
altitude_data <- cbind(singularity_data, taxa_formodel[ , grep("altitude", names(taxa_formodel))])
altitude_vars <- colnames(altitude_data[ , c(17:18)])
dif_temp <- taxa_formodel[ , "dif_temp"]
LGM_data <- cbind(singularity_data, dif_temp)
LGM_vars <- colnames(LGM_data)[17]

# Select variable response to work with
response_variables <- c("mean.sing.value", "a.exp", "b.exp")
selected_response <- menu(c("Singularity", "Intercept", "Slope"), title = "Please select a response variable to work with:")


# Identify weakly correlated variables (R < 0.8) and select the more relevant for the target species, according to VIF
uncor_current_vars <- corSelect(data = current_climatic_data, sp.cols = response_variables[[selected_response]], var.cols = current_climatic_vars, cor.thresh = 0.8, select = "VIF")

uncor_altitude_vars <- corSelect(data = altitude_data, sp.cols = response_variables[[selected_response]], var.cols = altitude_vars, cor.thresh = 0.8, select = "VIF")

# Select predictor variables for GLM using a stepwise procedure
null_formula <- formula(paste(response_variables[[selected_response]], "~ 1"))
glm_null <- glm(formula = null_formula, data = taxa_formodel, family = gaussian(link = "identity"))

glm_maximal_current <- formula(paste(response_variables[[selected_response]], "~ ", paste(uncor_current_vars$selected.vars, collapse = " + ")))
glm_maximal_altitude <- formula(paste(response_variables[[selected_response]], "~ ", paste(uncor_altitude_vars$selected.vars, collapse = " + ")))
glm_maximal_LGM<- formula(paste(response_variables[[selected_response]], "~ ", LGM_vars))

glm_step_current <- step(glm_null, glm_maximal_current, data = current_climatic_data, direction = "both", trace = TRUE)
  glm_step_current
  variance_glm_current <- (glm_step_current$null.deviance - glm_step_current$deviance) / glm_step_current$null.deviance
glm_step_altitude <- step(glm_null, glm_maximal_altitude, data = altitude_data, direction = "both", trace = TRUE)
  glm_step_altitude  
  variance_glm_altitude <- (glm_step_altitude$null.deviance - glm_step_altitude$deviance) / glm_step_altitude$null.deviance
glm_step_LGM <- step(glm_null, glm_maximal_LGM, data = LGM_data, direction = "both", trace = TRUE)
  glm_step_LGM  
  variance_glm_LGM <- (glm_step_LGM$null.deviance - glm_step_LGM$deviance) / glm_step_LGM$null.deviance
  

# Compare by pairs:
formula_current_altitude <- formula(paste(response_variables[[selected_response]], "~ ", paste(glm_step_current$formula, glm_step_altitude$formula, sep = " + ")[3]))
current_altitude_data <- cbind(singularity_data, taxa_formodel[ , grep("bio", names(taxa_formodel))], taxa_formodel[ , grep("altitude", names(taxa_formodel))])

glm_current_altitude <- glm(formula = formula_current_altitude, data = current_altitude_data, family = gaussian(link = "identity"))
  glm_current_altitude
  variance_glm_current_altitude <- (glm_current_altitude$null.deviance - glm_current_altitude$deviance) / glm_current_altitude$null.deviance
  
  
formula_current_LGM <- formula(paste(response_variables[[selected_response]], "~ ", paste(glm_step_current$formula, glm_step_LGM$formula, sep = " + ")[3]))
current_LGM_data <- cbind(singularity_data, taxa_formodel[ , grep("bio", names(taxa_formodel))], dif_temp)

glm_current_LGM <- glm(formula = formula_current_LGM, data = current_LGM_data, family = gaussian(link = "identity"))
  glm_current_LGM
  variance_glm_current_LGM <- (glm_current_LGM$null.deviance - glm_current_LGM$deviance) / glm_current_LGM$null.deviance


formula_altitude_LGM <- formula(paste(response_variables[[selected_response]], "~ ", paste(glm_step_altitude$formula, glm_step_LGM$formula, sep = " + ")[3]))
altitude_LGM_data <- cbind(singularity_data, taxa_formodel[ , grep("altitude", names(taxa_formodel))], dif_temp)

glm_altitude_LGM <- glm(formula = formula_altitude_LGM, data = altitude_LGM_data, family = gaussian(link = "identity"))
  glm_altitude_LGM
  variance_glm_altitude_LGM <- (glm_altitude_LGM$null.deviance - glm_altitude_LGM$deviance) / glm_altitude_LGM$null.deviance
  
# Global model
formula_global <- formula(paste(response_variables[[selected_response]], "~ ", paste(glm_step_current$formula, glm_step_altitude$formula, glm_step_LGM$formula, sep = " + ")[3]))

glm_global<- glm(formula = formula_global, data = taxa_formodel, family = gaussian(link = "identity"))
  glm_global
  variance_glm_global <- (glm_global$null.deviance - glm_global$deviance) / glm_global$null.deviance
  
variances <- c(variance_glm_current, variance_glm_altitude, variance_glm_LGM, variance_glm_current_altitude, variance_glm_current_LGM, variance_glm_altitude_LGM, variance_glm_global)
  names(variances) <- c("variance_glm_current", "variance_glm_altitude", "variance_glm_LGM", "variance_glm_current_altitude", "variance_glm_current_LGM", "variance_glm_altitude_LGM", "variance_glm_global")
  
taxa <- c("Amphibians", "Reptiles", "Birds", "Mammals")
  write.csv(variances, paste("../variance", response_variables[[selected_response]], taxa[selected_taxa], ".csv", sep = "_" ))
    
        
# 3.2. CLIMATIC AND SPATIAL VARIABLES 

# Create subsets of  grouped-by-variables data
spatial_data <- singularity_data

# Select predictor variables for GLM using a stepwise procedure
glm_maximal_spatial <- formula(paste(response_variables[[selected_response]], "~ long + lat + I(long * lat) + I(long ^2) + I(lat^2)"))


glm_step_spatial <- step(glm_null, glm_maximal_spatial, data = spatial_data, direction = "both", trace = TRUE)
  glm_step_spatial
  variance_glm_spatial <- (glm_step_spatial$null.deviance - glm_step_spatial$deviance) / glm_step_spatial$null.deviance


# Compare by pairs:
formula_current_spatial<- formula(paste(response_variables[[selected_response]], "~ ", paste(glm_step_current$formula, glm_step_spatial$formula, sep = " + ")[3]))
current_spatial_data <- cbind(singularity_data, taxa_formodel[ , grep("bio", names(taxa_formodel))])

glm_current_spatial <- glm(formula = formula_current_spatial, data = current_spatial_data, family = gaussian(link = "identity"))
  glm_current_spatial
  variance_glm_current_spatial <- (glm_current_spatial$null.deviance - glm_current_spatial$deviance) / glm_current_spatial$null.deviance

formula_altitude_spatial<- formula(paste(response_variables[[selected_response]], "~ ", paste(glm_step_altitude$formula, glm_step_spatial$formula, sep = " + ")[3]))
altitude_spatial_data <- cbind(singularity_data, taxa_formodel[ , grep("altitude", names(taxa_formodel))])

glm_altitude_spatial <- glm(formula = formula_altitude_spatial, data = altitude_spatial_data, family = gaussian(link = "identity"))
  glm_altitude_spatial
  variance_glm_altitude_spatial <- (glm_altitude_spatial$null.deviance - glm_altitude_spatial$deviance) / glm_altitude_spatial$null.deviance
  
formula_LGM_spatial<- formula(paste(response_variables[[selected_response]], "~ ", paste(glm_step_LGM$formula, glm_step_spatial$formula, sep = " + ")[3]))
  LGM_spatial_data <- cbind(singularity_data, dif_temp)
  
  glm_LGM_spatial <- glm(formula = formula_LGM_spatial, data = LGM_spatial_data, family = gaussian(link = "identity"))
  glm_LGM_spatial
  variance_glm_LGM_spatial <- (glm_LGM_spatial$null.deviance - glm_LGM_spatial$deviance) / glm_LGM_spatial$null.deviance

# Compare by 3
formula_current_altitude_LGM<- formula(paste(response_variables[[selected_response]], "~ ", paste(glm_step_current$formula, glm_step_altitude$formula, glm_step_LGM$formula, sep = " + ")[3]))
  current_altitude_LGM_data <- cbind(singularity_data, taxa_formodel[ , grep("bio", names(taxa_formodel))], taxa_formodel[ , grep("altitude", names(taxa_formodel))], dif_temp)
  
glm_current_altitude_LGM <- glm(formula = formula_current_altitude_LGM, data = current_altitude_LGM_data, family = gaussian(link = "identity"))
  glm_current_altitude_LGM
  variance_glm_current_altitude_LGM <- (glm_current_altitude_LGM$null.deviance - glm_current_altitude_LGM$deviance) / glm_current_altitude_LGM$null.deviance

formula_current_altitude_spatial<- formula(paste(response_variables[[selected_response]], "~ ", paste(glm_step_current$formula, glm_step_altitude$formula, glm_step_spatial$formula, sep = " + ")[3]))
  current_altitude_spatial_data <- cbind(singularity_data, taxa_formodel[ , grep("bio", names(taxa_formodel))], taxa_formodel[ , grep("altitude", names(taxa_formodel))])
  
glm_current_altitude_spatial <- glm(formula = formula_current_altitude_spatial, data = current_altitude_spatial_data, family = gaussian(link = "identity"))
  glm_current_altitude_spatial
  variance_glm_current_altitude_spatial <- (glm_current_altitude_spatial$null.deviance - glm_current_altitude_spatial$deviance) / glm_current_altitude_spatial$null.deviance

formula_current_spatial_LGM<- formula(paste(response_variables[[selected_response]], "~ ", paste(glm_step_current$formula, glm_step_spatial$formula, glm_step_LGM$formula, sep = " + ")[3]))
  current_spatial_LGM_data <- cbind(singularity_data, taxa_formodel[ , grep("bio", names(taxa_formodel))],  dif_temp)
  
glm_current_spatial_LGM <- glm(formula = formula_current_spatial_LGM, data = current_spatial_LGM_data, family = gaussian(link = "identity"))
  glm_current_spatial_LGM
  variance_glm_current_spatial_LGM <- (glm_current_spatial_LGM$null.deviance - glm_current_spatial_LGM$deviance) / glm_current_spatial_LGM$null.deviance
  
formula_LGM_altitude_spatial<- formula(paste(response_variables[[selected_response]], "~ ", paste(glm_step_LGM$formula, glm_step_altitude$formula, glm_step_spatial$formula, sep = " + ")[3]))
  LGM_altitude_spatial_data <- cbind(singularity_data, dif_temp, taxa_formodel[ , grep("altitude", names(taxa_formodel))])
  
glm_LGM_altitude_spatial <- glm(formula = formula_LGM_altitude_spatial, data = LGM_altitude_spatial_data, family = gaussian(link = "identity"))
  glm_LGM_altitude_spatial
  variance_glm_LGM_altitude_spatial <- (glm_LGM_altitude_spatial$null.deviance - glm_LGM_altitude_spatial$deviance) / glm_LGM_altitude_spatial$null.deviance
  
  
# Global model
formula_global <- formula(paste(response_variables[[selected_response]], "~ ", paste(glm_step_current$formula, glm_step_altitude$formula, glm_step_LGM$formula, glm_step_spatial$formula, sep = " + ")[3]))

glm_global<- glm(formula = formula_global, data = taxa_formodel, family = gaussian(link = "identity"))
  glm_global
  variance_glm_global <- (glm_global$null.deviance - glm_global$deviance) / glm_global$null.deviance
  
variances <- c(variance_glm_current, variance_glm_altitude, variance_glm_LGM, variance_glm_spatial, variance_glm_current_altitude, variance_glm_current_LGM, variance_glm_altitude_LGM, variance_glm_current_spatial, variance_glm_altitude_spatial, variance_glm_LGM_spatial, variance_glm_current_altitude_LGM, variance_glm_current_altitude_spatial, variance_glm_current_spatial_LGM, variance_glm_LGM_altitude_spatial, variance_glm_global)
  names(variances) <- c("variance_glm_current", "variance_glm_altitude", "variance_glm_LGM", "variance_glm_spatial",  "variance_glm_current_altitude", "variance_glm_current_LGM", "variance_glm_altitude_LGM", "variance_glm_current_spatial", "variance_glm_altitude_spatial", "variance_glm_LGM_spatial", "variance_glm_current_altitude_LGM", "variance_glm_current_altitude_spatial", "variance_glm_current_spatial_LGM", "variance_glm_LGM_altitude_spatial", "variance_glm_global")
  
taxa <- c("Amphibians", "Reptiles", "Birds", "Mammals")
  write.csv(variances, paste("../variance_spatial", response_variables[[selected_response]], taxa[selected_taxa], ".csv", sep = "_" ))
            