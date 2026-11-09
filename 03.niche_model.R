# Loading required libraries
library(terra)
library(sf)
library(dplyr)
library(readr)
library(ENMeval)

#########################################
# Settings
periods <- paste0("yr", 1:5)

dir.create("data/bg", showWarnings = FALSE)
dir.create("output/sdm/gbif", showWarnings = FALSE, recursive = TRUE)
dir.create("output/sdm/combined", showWarnings = FALSE, recursive = TRUE)


#########################################
# Creating a common calibration area
# Use all combined occurrence records to define the study area
com <- read_csv("data/occ/com.csv", show_col_types = FALSE)

# Environmental raster used only as the spatial template
env_overall <- rast("data/env/env_overall.tif")

# Occurrence coordinates
occs_all <- com %>%
  dplyr::select(longitude, latitude) %>%
  distinct()

# Convert occurrences to sf
occs_sf <- st_as_sf(occs_all,
                    coords = c("longitude", "latitude"),
                    crs = crs(env_overall))

# Equal-area projection for buffering
eckertIV <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# Create 500-km buffer around all occurrences
calibration_area <- occs_sf %>%
  st_transform(eckertIV) %>%
  st_buffer(dist = 1000000) %>%
  st_union() %>%
  st_sf() %>%
  st_transform(crs(env_overall))

calibration_vect <- vect(calibration_area)

#########################################
# Importing background points

bg <- read_csv("data/bg.csv", show_col_types = FALSE)


#########################################
# Function for running SDMs

run_sdm <- function(data_file, dataset_name) {
  
  # Import occurrence data
  occ_data <- read_csv(data_file, show_col_types = FALSE)
  
  for (i in periods) {
    
    print(paste(dataset_name, i))
    
    
    # Import environmental predictors for this year interval
    envs <- rast(paste0("data/env/env_", i, ".tif"))
    
    
    # Restrict environmental predictors to calibration area
    envs_bg <- crop(envs, calibration_vect, mask = TRUE)
    
    
    # Filter occurrence data by year interval
    occs <- occ_data %>%
      filter(year_group == i) %>%
      dplyr::select(longitude, latitude)
    
    
    # Remove occurrences without complete environmental data
    occ_values <- terra::extract(
      envs_bg,
      as.matrix(occs[, c("longitude", "latitude")])
    )
    
    occs <- occs[complete.cases(occ_values), ]
    
    
    # Keep background points with complete environmental data
    bg_values <- terra::extract(
      envs_bg,
      as.matrix(bg[, c("longitude", "latitude")])
    )
    
    bg_i <- bg[complete.cases(bg_values), ]
    
    
    # Run Maxent models
    model <- ENMevaluate(
      occs = occs,
      envs = envs_bg,
      bg = bg_i,
      partitions = "checkerboard",
      partition.settings = list(
        aggregation.factor = c(2, 2)
      ),
      tune.args = list(
        fc = c("L", "LQ", "LQP", "LQPT"), # Hinge caused errors during cross validation
        rm = seq(0.5, 4, 0.5)
      ),
      algorithm = "maxnet",
      raster.preds = FALSE,
      rmm = list(),
      parallel = FALSE
    )
    
    
    # Model evaluation results
    res <- eval.results(model)
    
    
    # Select model with lowest AICc
    opt_aicc <- res %>%
      filter(!is.na(AICc)) %>%
      slice_min(AICc, n = 1, with_ties = FALSE)
    
    
    # Extract fitted models
    models <- eval.models(model)
    
    
    # Identify selected model
    best_name <- as.character(opt_aicc$tune.args)
    best_model <- models[[best_name]]
    
    
    # Make suitability raster for selected model only
    r <- ENMeval::maxnet.predictRaster(
      best_model,
      envs_bg
    )
    
    
    # Extract suitability values at occurrence locations
    occ_predictions <- terra::extract(
      r,
      as.matrix(occs[, c("longitude", "latitude")])
    )
    
    
    # Calculate 10th-percentile occurrence threshold
    threshold <- quantile(
      occ_predictions[, 1],
      probs = 0.10,
      na.rm = TRUE
    )
    
    
    # Add threshold to model summary
    opt_aicc$threshold_10p <- threshold
    
    
    # Export model summary
    write_csv(
      opt_aicc,
      paste0(
        "output/sdm/", dataset_name, "/",
        dataset_name, "_model_sum_", i, ".csv"
      )
    )
    
    
    # Save complete ENMeval object
    saveRDS(
      model,
      paste0(
        "output/sdm/", dataset_name, "/",
        dataset_name, "_model_", i, ".rds"
      )
    )
    
    
    # Export suitability raster
    writeRaster(
      r,
      paste0(
        "output/sdm/", dataset_name, "/",
        dataset_name, "_map_", i, ".tif"
      ),
      NAflag = -9999,
      overwrite = TRUE
    )
    
    
    # Create binary suitability raster
    r_bin <- ifel(
      r >= threshold,
      1,
      0
    )
    
    
    # Export binary raster
    writeRaster(
      r_bin,
      paste0(
        "output/sdm/", dataset_name, "/",
        dataset_name, "_bin_", i, ".tif"
      ),
      NAflag = -9999,
      overwrite = TRUE
    )
  }
}

#########################################
# Run GBIF-only models
run_sdm(data_file = "data/occ/gbif.csv", dataset_name = "gbif")

#########################################
# Run combined-data models
run_sdm(data_file = "data/occ/com.csv", dataset_name = "combined")

