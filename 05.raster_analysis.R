# Loading libraries
library(terra)
library(geosphere)
library(rnaturalearth)
library(sf)

years <- c("yr1", "yr2", "yr3", "yr4", "yr5")
sources <- c("gbif", "combined")

old_year <- c("yr1", "yr2", "yr3", "yr4", "yr1")
new_year <- c("yr2", "yr3", "yr4", "yr5", "yr5")

#####################################
# Suitable area, centroid and mean elevation
#####################################
period_metrics <- data.frame()

# Suitable area, centroid and mean elevation
#####################################

period_metrics <- data.frame()

for (source in sources) {
  
  centroid_yr1 <- NULL
  
  for (i in years) {
    
    print(i)
    
    if (source == "gbif") {
      r <- rast(paste0("output/bin_maps_full/gbif_bin_", i, ".tif"))
    } else {
      r <- rast(paste0("output/bin_maps_full/combined_bin_", i, ".tif"))
    }
    
    env <- rast(paste0("data/env/env_", i, ".tif"))
    
    # Make CRS metadata consistent
    crs(env) <- crs(r)
    
    elev <- env[["elev"]]
    
    # Cell area in km2
    area <- cellSize(r, unit = "km")
    
    # Suitable area
    suitable_area <- global(
      ifel(r == 1, area, NA),
      "sum",
      na.rm = TRUE
    )[1, 1]
    
    # Centroid
    cells <- which(values(r)[, 1] == 1)
    centroid <- colMeans(xyFromCell(r, cells))
    
    if (i == "yr1") {
      centroid_yr1 <- centroid
    }
    
    centroid_shift <- distHaversine(
      centroid_yr1,
      centroid
    ) / 1000
    
    # Mean elevation of suitable cells
    elev_suitable <- ifel(r == 1, elev, NA)
    
    mean_elev <- global(
      elev_suitable,
      "mean",
      na.rm = TRUE
    )[1, 1]
    
    period_metrics <- rbind(
      period_metrics,
      data.frame(
        source = source,
        year = i,
        suitable_area_km2 = suitable_area,
        centroid_lon = centroid[1],
        centroid_lat = centroid[2],
        centroid_shift_from_yr1_km = centroid_shift,
        mean_elevation_m = mean_elev
      )
    )
  }
}

write.csv(period_metrics, "output/period_metrics.csv", row.names = FALSE)

#####################################
# Temporal gain, loss and persistence
#####################################
years <- c("yr1", "yr2", "yr3", "yr4", "yr5")
sources <- c("gbif", "combined")

old_year <- c("yr1", "yr2", "yr3", "yr4", "yr1")
new_year <- c("yr2", "yr3", "yr4", "yr5", "yr5")

temporal_change <- data.frame()
country_change <- data.frame()
period_metrics <- read.csv("output/period_metrics.csv")
period_metrics$source[
  period_metrics$source == "com"
] <- "combined"

# Cell area only needs to be calculated once
template <- rast("output/bin_maps_full/gbif_bin_yr1.tif")
area <- cellSize(template, unit = "km")

# Country raster only needs to be created once
countries <- ne_countries(
  scale = "medium",
  returnclass = "sv"
)

# Same geographic CRS, make metadata consistent
crs(countries) <- crs(template)

countries$id <- 1:nrow(countries)

country_raster <- rasterize(
  countries,
  template,
  field = "id"
)

for (source in sources) {
  
  for (j in seq_along(old_year)) {
    
    print(paste(source, old_year[j], "to", new_year[j]))
    
    # Load binary maps
    
    if (source == "gbif") {
      
      old <- rast(
        paste0(
          "output/bin_maps_full/gbif_bin_",
          old_year[j],
          ".tif"
        )
      )
      
      new <- rast(
        paste0(
          "output/bin_maps_full/gbif_bin_",
          new_year[j],
          ".tif"
        )
      )
      
    } else {
      
      old <- rast(
        paste0(
          "output/bin_maps_full/combined_bin_",
          old_year[j],
          ".tif"
        )
      )
      
      new <- rast(
        paste0(
          "output/bin_maps_full/combined_bin_",
          new_year[j],
          ".tif"
        )
      )
    }
    
    # Temporal change raster
    # 0 = unsuitable in both
    # 1 = loss
    # 2 = gain
    # 3 = persistent suitability
    #
    # NA remains NA if either period is NA
    
    change <- old + 2 * new
    
    
    writeRaster(
      change,
      paste0(
        "output/",
        source,
        "_change_",
        old_year[j],
        "_",
        new_year[j],
        ".tif"
      ),
      overwrite = TRUE
    )
    
    # Area of loss, gain and persistence
    loss_area <- mask(
      area,
      change == 1,
      maskvalues = 0
    )
    
    gain_area <- mask(
      area,
      change == 2,
      maskvalues = 0
    )
    
    persistence_area <- mask(
      area,
      change == 3,
      maskvalues = 0
    )
    
    
    change_area <- c(
      loss_area,
      gain_area,
      persistence_area
    )
    
    names(change_area) <- c(
      "loss_km2",
      "gain_km2",
      "persistence_km2"
    )
    
    
    totals <- global(
      change_area,
      "sum",
      na.rm = TRUE
    )
    
    loss <- totals["loss_km2", 1]
    gain <- totals["gain_km2", 1]
    persistence <- totals["persistence_km2", 1]
    
    # Centroid and elevational change
    old_metrics <- period_metrics[
      period_metrics$source == source &
        period_metrics$year == old_year[j],
    ]
    
    new_metrics <- period_metrics[
      period_metrics$source == source &
        period_metrics$year == new_year[j],
    ]
    
    
    centroid_shift <- distHaversine(
      c(
        old_metrics$centroid_lon,
        old_metrics$centroid_lat
      ),
      c(
        new_metrics$centroid_lon,
        new_metrics$centroid_lat
      )
    ) / 1000
    
    
    elevational_change <-
      new_metrics$mean_elevation_m -
      old_metrics$mean_elevation_m
    
    # Save overall temporal results
    temporal_change <- rbind(
      temporal_change,
      data.frame(
        source = source,
        transition = paste0(
          old_year[j],
          "_to_",
          new_year[j]
        ),
        old_area_km2 = old_metrics$suitable_area_km2,
        new_area_km2 = new_metrics$suitable_area_km2,
        gain_km2 = gain,
        loss_km2 = loss,
        persistence_km2 = persistence,
        net_change_km2 = gain - loss,
        centroid_shift_km = centroid_shift,
        elevational_change_m = elevational_change
      )
    )
    
    # Country-wise gain and loss
    country <- zonal(
      change_area,
      country_raster,
      fun = "sum",
      na.rm = TRUE
    )
    
    names(country)[1] <- "id"
    
    
    country$country <- countries$admin[
      match(country$id, countries$id)
    ]
    
    country$iso3 <- countries$iso_a3[
      match(country$id, countries$id)
    ]
    
    country$source <- source
    
    country$transition <- paste0(
      old_year[j],
      "_to_",
      new_year[j]
    )
    
    country$net_change_km2 <-
      country$gain_km2 -
      country$loss_km2
    
    
    country <- country[
      ,
      c(
        "country",
        "iso3",
        "source",
        "transition",
        "gain_km2",
        "loss_km2",
        "persistence_km2",
        "net_change_km2"
      )
    ]
    
    
    country_change <- rbind(
      country_change,
      country
    )
  }
}

# Export output
write.csv(temporal_change, "output/temporal_range_change.csv", row.names = FALSE)
write.csv(country_change, "output/country_range_change.csv", row.names = FALSE)

#####################################
# GBIF vs combined maps
#####################################
#####################################
# GBIF vs combined maps
#####################################
gbif_combined <- data.frame()
combined_only_environment <- data.frame()

# Cell area only needs to be calculated once
template <- rast("output/bin_maps_full/gbif_bin_yr1.tif")
area <- cellSize(template, unit = "km")

for (i in years) {
  
  print(i)
  
  gbif <- rast(
    paste0(
      "output/bin_maps_full/gbif_bin_",
      i,
      ".tif"
    )
  )
  
  com <- rast(
    paste0(
      "output/bin_maps_full/combined_bin_",
      i,
      ".tif"
    )
  )
  
  # Compare GBIF and combined maps
  # 0 = unsuitable in both
  # 1 = GBIF only
  # 2 = combined only
  # 3 = suitable in both
  
  comparison <- gbif + 2 * com
  
  
  writeRaster(
    comparison,
    paste0(
      "output/gbif_combined_",
      i,
      ".tif"
    ),
    overwrite = TRUE
  )
  
  # Area in each category
  gbif_only_area <- mask(
    area,
    comparison == 1,
    maskvalues = 0
  )
  
  com_only_area <- mask(
    area,
    comparison == 2,
    maskvalues = 0
  )
  
  both_area <- mask(
    area,
    comparison == 3,
    maskvalues = 0
  )
  
  
  comparison_area <- c(
    gbif_only_area,
    com_only_area,
    both_area
  )
  
  names(comparison_area) <- c(
    "gbif_only_km2",
    "combined_only_km2",
    "overlap_km2"
  )
  
  
  totals <- global(
    comparison_area,
    "sum",
    na.rm = TRUE
  )
  
  gbif_only <- totals["gbif_only_km2", 1]
  com_only <- totals["combined_only_km2", 1]
  both <- totals["overlap_km2", 1]
  
  
  # Jaccard geographic overlap
  geographic_overlap <-
    both / (gbif_only + com_only + both)
  
  gbif_combined <- rbind(
    gbif_combined,
    data.frame(
      year = i,
      gbif_only_km2 = gbif_only,
      combined_only_km2 = com_only,
      overlap_km2 = both,
      geographic_overlap_jaccard = geographic_overlap
    )
  )
  
  # Environment in combined-only cells
  env <- rast(
    paste0(
      "data/env/env_",
      i,
      ".tif"
    )
  )
  
  # Equivalent CRS, make metadata consistent
  crs(env) <- crs(comparison)
  
  
  # Keep only cells unique to combined model
  env_unique <- mask(
    env,
    comparison == 2,
    maskvalues = 0
  )
  
  
  env_mean <- global(
    env_unique,
    "mean",
    na.rm = TRUE
  )
  
  env_sd <- global(
    env_unique,
    "sd",
    na.rm = TRUE
  )
  
  
  combined_only_environment <- rbind(
    combined_only_environment,
    data.frame(
      year = i,
      variable = names(env),
      mean = env_mean[, 1],
      sd = env_sd[, 1]
    )
  )
}

write.csv(gbif_combined, "output/gbif_combined_comparison.csv", row.names = FALSE)
write.csv(combined_only_environment, "output/combined_only_environment.csv", row.names = FALSE)