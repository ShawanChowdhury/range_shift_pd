# Loading libraries
library(terra)

# Full-extent template
template <- rast("data/env/env_overall.tif")[[1]]

# Assign explicit WGS84 longitude-latitude CRS
crs(template) <- "OGC:CRS84"

# List binary rasters
tifs <- list.files(
  path = "output/sdm/",
  pattern = "_bin_.*\\.tif$",
  recursive = TRUE,
  full.names = TRUE
)

# Create output folder
dir.create("output/bin_maps_full", showWarnings = FALSE, recursive = TRUE)

for (i in tifs) {
  
  print(i)
  
  # File name
  file_name <- tools::file_path_sans_ext(basename(i))
  
  # Read binary raster
  ras <- rast(i)
  
  # Standardise CRS metadata
  crs(ras) <- "OGC:CRS84"
  
  # Resample onto full-extent grid using nearest neighbour
  ras_full <- resample(
    ras,
    template,
    method = "near"
  )
  
  # Export raster
  writeRaster(
    ras_full,
    paste0("output/bin_maps_full/", file_name, ".tif"),
    NAflag = -9999,
    overwrite = TRUE
  )
}
