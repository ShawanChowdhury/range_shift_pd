# Loading libraries
library(terra)

# Full-extent template
template <- rast("data/env/env_overall.tif")[[1]]

# Assign explicit WGS84 longitude-latitude CRS
crs(template) <- "OGC:CRS84"

# Five-year windows (>=3 years)
annual_windows <- list(yr1 = 2001:2005, yr2 = 2006:2010, yr3 = 2011:2015,
                       yr4 = 2016:2020, yr5 = 2021:2025)
dir.create("output/sdm/annual5", showWarnings = FALSE, recursive = TRUE)
for (p in names(annual_windows)) {
  yrs <- terra::rast(sprintf("output/sdm/annual/annual_bin_%d.tif", annual_windows[[p]]))
  writeRaster(ifel(sum(yrs) >= 3, 1, 0),
              paste0("output/sdm/annual5/annual5_bin_", p, ".tif"),
              NAflag = -9999, overwrite = TRUE)
}

# List period binary rasters
tifs <- list.files(
  path = "output/sdm/",
  pattern = "_bin_yr[1-5]\\.tif$",
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

  # Resample to full extent
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
