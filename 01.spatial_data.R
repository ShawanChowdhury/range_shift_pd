# Loading required libraries
library(terra)
library(dplyr)
library(readr)

#########################################
# Settings
# TerraClimate predictor variables
pred_vars <- c("aet", "def", "PDSI", "pet", "ppt", "q", "soil", "tmax", "tmin", "ws")

# Five 5-year periods
period_years <- list(
  yr1 = 2001:2005,
  yr2 = 2006:2010,
  yr3 = 2011:2015,
  yr4 = 2016:2020,
  yr5 = 2021:2025
)

# Create output folders
dir.create("data/pred_mean", showWarnings = FALSE, recursive = TRUE)
dir.create("data/pred_mean/ov", showWarnings = FALSE, recursive = TRUE)

#########################################
# Calculating mean rasters by year groups

# Variables that represent monthly totals
sum_vars <- c("aet", "def", "pet", "ppt", "q")

for (period in names(period_years)) {
  
  print(period)
  years <- period_years[[period]]
  
  for (v in pred_vars) {
    
    print(v)
    
    # Select the five yearly NetCDF files
    files <- file.path("data/climate",
                       paste0("TerraClimate_", v, "_", years, ".nc"))
    
    if (!all(file.exists(files))) {
      stop("Missing files for ", period, " / ", v)
    }
    
    # Summarise each year's 12 monthly layers
    annual_rasters <- lapply(files, function(f) {
      
      r <- rast(f)
      
      if (nlyr(r) != 12) {
        stop("Expected 12 monthly layers in ", basename(f))
      }
      
      if (v %in% sum_vars) {
        sum(r, na.rm = TRUE)
      } else {
        mean(r, na.rm = TRUE)
      }
    })
    
    # Combine the five annual rasters
    r <- do.call(c, annual_rasters)
    
    # Calculate mean across the five years
    mean_ras <- mean(r, na.rm = TRUE)
    
    # Export raster
    writeRaster(mean_ras,
                paste0("data/pred_mean/", period, "_", v, ".tif"),
                NAflag = -9999,
                overwrite = TRUE)
  }
}


#########################################
# Calculating overall mean rasters for 2001-2025

for (v in pred_vars) {
  
  print(v)
  
  # Read the five period rasters
  files <- file.path("data/pred_mean",
                     paste0(names(period_years), "_", v, ".tif"))
  
  r <- rast(files)
  
  # Calculate overall mean
  mean_ras <- mean(r, na.rm = TRUE)
  
  # Export raster
  writeRaster(mean_ras,
              paste0("data/pred_mean/ov/", v, ".tif"),
              NAflag = -9999,
              overwrite = TRUE)
}


#########################################
# Resampling elevation raster

# Use one of the climate rasters as the template
template <- rast("data/pred_mean/yr1_aet.tif")

# Import elevation raster
elev <- rast("data/climate/elev/elev.tif")

# Match elevation to climate raster
if (!same.crs(elev, template)) {
  elev_re <- project(elev, template, method = "bilinear")
} else {
  elev_re <- resample(elev, template, method = "bilinear")
}

# Apply the same spatial mask
elev_re <- mask(elev_re, template)

# Save the same elevation layer for each year interval
for (period in names(period_years)) {
  
  writeRaster(elev_re,
              paste0("data/pred_mean/", period, "_elev.tif"),
              NAflag = -9999,
              overwrite = TRUE)
}

# Also save elevation with the overall 2001-2025 predictors
writeRaster(elev_re,
            "data/pred_mean/ov/elev.tif",
            NAflag = -9999,
            overwrite = TRUE)


#########################################
# Predictor variable collinearity

# All candidate predictor variables
vars <- c("aet", "def", "PDSI", "pet", "ppt", "q", "soil",
          "tmax", "tmin", "ws", "elev")

# Attach environmental rasters
env_all <- rast(file.path("data/pred_mean/ov", paste0(vars, ".tif")))
names(env_all) <- vars

# Randomly sample environmental values
set.seed(123)

env_sample <- spatSample(env_all,
                         size = 20000,
                         method = "random",
                         na.rm = TRUE,
                         as.df = TRUE)

# Correlation matrix
cor_mat <- cor(env_sample, use = "pairwise.complete.obs")

# View correlations
round(cor_mat, 2)

# Plot relationships between predictors
pairs(env_sample)

# Identify correlations > 0.75
cor_index <- which(abs(cor_mat) > 0.75 & upper.tri(cor_mat), arr.ind = TRUE)

correlated_pairs <- data.frame(
  variable1 = rownames(cor_mat)[cor_index[, 1]],
  variable2 = colnames(cor_mat)[cor_index[, 2]],
  correlation = cor_mat[cor_index]
) %>%
  arrange(desc(abs(correlation)))

correlated_pairs


#########################################
# Final selected predictor variables
vars_final <- c("PDSI", "ppt", "soil", "tmin", "ws", "elev")

# Create folder for final environmental predictors
dir.create("data/env", showWarnings = FALSE, recursive = TRUE)

#########################################
# Overall environmental predictors: 2001-2025
env <- rast(file.path("data/pred_mean/ov",
                      paste0(vars_final, ".tif")))

names(env) <- vars_final

# Save as one multi-layer raster
writeRaster(env,
            "data/env/env_overall.tif",
            NAflag = -9999,
            overwrite = TRUE)


#########################################
# Environmental predictors for each 5-year interval
for (period in names(period_years)) {
  
  print(period)
  
  # Select only the final predictor variables
  files <- file.path("data/pred_mean",
                     paste0(period, "_", vars_final, ".tif"))
  
  # Check that all files exist
  if (!all(file.exists(files))) {
    stop("Missing predictor files for ", period)
  }
  
  # Combine selected predictors
  env_period <- rast(files)
  names(env_period) <- vars_final
  
  # Save as one multi-layer raster
  writeRaster(env_period,
              paste0("data/env/env_", period, ".tif"),
              NAflag = -9999,
              overwrite = TRUE)
}


#########################################
# Creating background points for the distribution model
# Create a raster containing only cells where all predictors have data
valid_cells <- ifel(sum(is.na(env)) == 0, 1, NA)

# Create 20,000 random background points
set.seed(123)

bg <- spatSample(valid_cells, size = 20000, method = "random", na.rm = TRUE, 
                 xy = TRUE, values = FALSE, as.df = TRUE)

# Convert to data frame
bg <- as.data.frame(bg)

# Rename coordinates
colnames(bg) <- c("longitude", "latitude")

# Check
head(bg)
class(bg)

# Export background points
write_csv(bg, "data/bg.csv")
