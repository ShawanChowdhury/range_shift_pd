# Load libraries
library(terra)
library(geosphere)
library(rnaturalearth)
library(sf)
library(dplyr)

# Settings and shared spatial layers
periods <- paste0("yr", 1:5)
# Model sets to analyse
sources <- strsplit(Sys.getenv("SDM_SOURCES", "gbif,combined,unified"), ",")[[1]]
old_year <- c("yr1", "yr2", "yr3", "yr4", "yr1")
new_year <- c("yr2", "yr3", "yr4", "yr5", "yr5")

# dir.create("output", showWarnings = FALSE, recursive = TRUE)

template <- terra::rast(
  sprintf("output/bin_maps_full/%s_bin_yr1.tif", sources[1])
)

# Land layer as in 07
land <- rnaturalearth::ne_download(
  scale = 50,
  type = "land",
  category = "physical",
  returnclass = "sf"
)
land_v <- terra::vect(land)

if (!terra::same.crs(land_v, template)) {
  land_v <- terra::project(land_v, terra::crs(template))
}

land_v$land_flag <- 1L
land_mask <- terra::rasterize(
  land_v,
  template,
  field = "land_flag",
  background = 0
)

# Rasterise countries on grid
countries <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sv"
)

if (!terra::same.crs(countries, template)) {
  countries <- terra::project(countries, terra::crs(template))
}

countries$id <- seq_len(nrow(countries))
country_raster <- terra::rasterize(
  countries,
  template,
  field = "id"
)
country_raster <- terra::mask(
  country_raster,
  land_mask,
  maskvalues = 0
)

cell_area <- terra::cellSize(template, unit = "km")

# Static elevation
elev <- terra::rast("data/env/env_yr1.tif")[["elev"]]

if (!terra::compareGeom(
  elev, template, stopOnError = FALSE
)) {
  elev <- terra::project(elev, template, method = "bilinear")
}

stopifnot(terra::compareGeom(
  elev, template, stopOnError = FALSE
))

# Helper functions
read_binary <- function(source, period) {
  file <- sprintf(
    "output/bin_maps_full/%s_bin_%s.tif",
    source, period
  )
  if (!file.exists(file)) stop("Missing raster: ", file)

  r <- terra::rast(file)
  if (!terra::compareGeom(r, template, stopOnError = FALSE)) {
    stop("Raster geometry does not match template: ", file)
  }

  # Keep terrestrial cells
  terra::mask(r, land_mask, maskvalues = 0)
}

area_table <- function(r, zones = NULL) {
  terra::expanse(
    r, unit = "km", byValue = TRUE, zones = zones
  )
}

area_for <- function(tab, code) {
  sum(tab$area[tab$value == code], na.rm = TRUE)
}

suitable_area <- function(r) {
  area_for(area_table(r), 1)
}

# Check matching land cells
assert_valid_extent <- function(old, new, label) {
  old_valid <- !is.na(old)
  new_valid <- !is.na(new)
  mismatch <- terra::global(
    old_valid != new_valid, "sum", na.rm = TRUE
  )[1, 1]
  if (!is.finite(mismatch) || mismatch > 0) {
    stop("Valid-data masks differ in: ", label,
         ". Reconcile masks before comparing areas.")
  }
}

# Area, centroid and elevation
period_rows <- list()
k <- 0L

for (source in sources) {
  for (period in periods) {
    message("Period metrics: ", source, " / ", period)
    r <- read_binary(source, period)

    # Area per raster class
    area_km2 <- suitable_area(r)

    # Area-weighted spherical centroid
    cells <- which(terra::values(r, mat = FALSE) == 1)
    if (!length(cells)) {
      stop("No suitable land cells: ", source, " / ", period)
    }

    xy <- terra::xyFromCell(r, cells)
    w <- terra::extract(cell_area, cells)[, 1]
    rad <- xy * pi / 180

    x <- sum(w * cos(rad[, 2]) * cos(rad[, 1]))
    y <- sum(w * cos(rad[, 2]) * sin(rad[, 1]))
    z <- sum(w * sin(rad[, 2]))

    centroid_lon <- atan2(y, x) * 180 / pi
    centroid_lat <- atan2(z, sqrt(x^2 + y^2)) * 180 / pi

    elevation <- terra::extract(elev, cells)[, 1]
    valid <- is.finite(elevation) & is.finite(w) & w > 0
    mean_elev <- if (any(valid)) {
      weighted.mean(elevation[valid], w[valid])
    } else {
      NA_real_
    }

    k <- k + 1L
    period_rows[[k]] <- data.frame(
      source = source,
      year = period,
      suitable_area_km2 = area_km2,
      centroid_lon = centroid_lon,
      centroid_lat = centroid_lat,
      centroid_shift_from_yr1_km = NA_real_,
      mean_elevation_m = mean_elev
    )
  }
}

period_metrics <- do.call(rbind, period_rows)

for (source in sources) {
  ids <- which(period_metrics$source == source)
  first_xy <- as.numeric(period_metrics[
    ids[1], c("centroid_lon", "centroid_lat")
  ])

  period_metrics$centroid_shift_from_yr1_km[ids] <- vapply(
    ids,
    function(i) {
      new_xy <- as.numeric(period_metrics[
        i, c("centroid_lon", "centroid_lat")
      ])
      geosphere::distGeo(first_xy, new_xy) / 1000
    },
    numeric(1)
  )
}

write.csv(
  period_metrics,
  "output/period_metrics.csv",
  row.names = FALSE
)

# 4. Temporal and country summaries
change_rows <- list()
country_rows <- list()
k <- 0L

for (source in sources) {
  for (j in seq_along(old_year)) {
    old_p <- old_year[j]
    new_p <- new_year[j]
    transition <- paste0(old_p, "_to_", new_p)
    message("Temporal change: ", source, " / ", transition)

    old <- read_binary(source, old_p)
    new <- read_binary(source, new_p)
    assert_valid_extent(old, new, transition)

    # 0 neither,1 loss,2 gain,3 persistence
    change <- old + 2 * new
    tab <- area_table(change)

    if (any(!tab$value %in% 0:3)) {
      stop("Unexpected binary change values: ", transition)
    }

    loss <- area_for(tab, 1)
    gain <- area_for(tab, 2)
    persistence <- area_for(tab, 3)
    old_area <- loss + persistence
    new_area <- gain + persistence

    # Validate independently calculated areas.
    old_check <- suitable_area(old)
    new_check <- suitable_area(new)
    if (abs(old_area - old_check) > 1 ||
        abs(new_area - new_check) > 1) {
      stop("Area reconciliation failed: ", transition)
    }

    old_metrics <- period_metrics[
      period_metrics$source == source &
        period_metrics$year == old_p, , drop = FALSE
    ]
    new_metrics <- period_metrics[
      period_metrics$source == source &
        period_metrics$year == new_p, , drop = FALSE
    ]

    centroid_shift <- geosphere::distGeo(
      c(old_metrics$centroid_lon, old_metrics$centroid_lat),
      c(new_metrics$centroid_lon, new_metrics$centroid_lat)
    ) / 1000

    elev_change <- new_metrics$mean_elevation_m -
      old_metrics$mean_elevation_m

    k <- k + 1L
    change_rows[[k]] <- data.frame(
      source = source,
      transition = transition,
      old_area_km2 = old_area,
      new_area_km2 = new_area,
      gain_km2 = gain,
      loss_km2 = loss,
      persistence_km2 = persistence,
      net_change_km2 = gain - loss,
      centroid_shift_km = centroid_shift,
      elevational_change_m = elev_change
    )

    writeRaster(
      change,
      sprintf("output/%s_change_%s_%s.tif",
              source, old_p, new_p),
      overwrite = TRUE
    )

    # Category areas by country
    tab_country <- area_table(change, zones = country_raster)
    country <- data.frame(
      id = countries$id,
      country = countries$admin,
      iso3 = countries$iso_a3
    )

    for (code in 1:3) {
      field <- c("loss_km2", "gain_km2",
                 "persistence_km2")[code]
      z <- tab_country[
        tab_country$value == code,
        c("zone", "area"), drop = FALSE
      ]
      country[[field]] <- z$area[match(country$id, z$zone)]
      country[[field]][is.na(country[[field]])] <- 0
    }

    country$source <- source
    country$transition <- transition
    country$net_change_km2 <- country$gain_km2 -
      country$loss_km2

    country_rows[[k]] <- country[
      , c("country", "iso3", "source", "transition",
          "gain_km2", "loss_km2", "persistence_km2",
          "net_change_km2")
    ]
  }
}

temporal_change <- do.call(rbind, change_rows)
country_change <- do.call(rbind, country_rows)

write.csv(
  temporal_change,
  "output/temporal_range_change.csv",
  row.names = FALSE
)
write.csv(
  country_change,
  "output/country_range_change.csv",
  row.names = FALSE
)

# GBIF-only vs combined maps
run_source_comparison <- all(c("gbif", "combined") %in% sources)
if (run_source_comparison) {
gbif_rows <- list()
env_rows <- list()

for (i in seq_along(periods)) {
  period <- periods[i]
  message("Source comparison: ", period)

  gbif <- read_binary("gbif", period)
  combined <- read_binary("combined", period)
  assert_valid_extent(gbif, combined, period)

  # 0 neither,1 GBIF,2 combined,3 both
  comparison <- gbif + 2 * combined
  tab <- area_table(comparison)

  if (any(!tab$value %in% 0:3)) {
    stop("Unexpected source comparison values: ", period)
  }

  gbif_only <- area_for(tab, 1)
  combined_only <- area_for(tab, 2)
  shared <- area_for(tab, 3)
  union <- gbif_only + combined_only + shared

  gbif_area <- gbif_only + shared
  combined_area <- combined_only + shared

  if (abs(gbif_area - suitable_area(gbif)) > 1 ||
      abs(combined_area - suitable_area(combined)) > 1) {
    stop("Source comparison failed area check: ", period)
  }

  gbif_rows[[i]] <- data.frame(
    year = period,
    gbif_only_km2 = gbif_only,
    combined_only_km2 = combined_only,
    overlap_km2 = shared,
    geographic_overlap_jaccard = if (union > 0) {
      shared / union
    } else {
      NA_real_
    }
  )

  writeRaster(
    comparison,
    paste0("output/gbif_combined_", period, ".tif"),
    overwrite = TRUE
  )

  # Align predictors to SDM grid
  env <- terra::rast(paste0("data/env/env_", period, ".tif"))
  if (!terra::compareGeom(
    env[[1]], template, stopOnError = FALSE
  )) {
    if (terra::same.crs(env, template)) {
      env <- terra::resample(env, template, method = "bilinear")
    } else {
      env <- terra::project(env, template, method = "bilinear")
    }
  }

  stopifnot(terra::compareGeom(
    env[[1]], template, stopOnError = FALSE
  ))

  # Combined-only cell environments
  combined_unique <- terra::ifel(comparison == 2, 1, NA)
  env_unique <- terra::mask(env, combined_unique)

  env_mean <- terra::global(
    env_unique, "mean", na.rm = TRUE
  )[, 1]
  env_sd <- terra::global(
    env_unique, "sd", na.rm = TRUE
  )[, 1]

  env_rows[[i]] <- data.frame(
    year = period,
    variable = names(env),
    mean = env_mean,
    sd = env_sd
  )
}

gbif_combined <- do.call(rbind, gbif_rows)
combined_only_environment <- do.call(rbind, env_rows)

write.csv(
  gbif_combined,
  "output/gbif_combined_comparison.csv",
  row.names = FALSE
)
write.csv(
  combined_only_environment,
  "output/combined_only_environment.csv",
  row.names = FALSE
)
}

# Final checks and summaries
stopifnot(all(abs(
  temporal_change$net_change_km2 -
    (temporal_change$new_area_km2 -
       temporal_change$old_area_km2)
) < 1))

print(period_metrics)
print(temporal_change)
if (run_source_comparison) print(gbif_combined)

#######################
# Largest net increases and decreases
country_change %>%
  filter(
    source == "combined",
    transition == "yr1_to_yr5"
  ) %>%
  arrange(desc(net_change_km2)) %>%
  select(
    country,
    gain_km2,
    loss_km2,
    persistence_km2,
    net_change_km2
  )

# Australia and China transitions
country_change %>%
  filter(
    source == "combined",
    country %in% c("Australia", "China")
  ) %>%
  select(
    country,
    transition,
    net_change_km2
  ) %>%
  arrange(country, transition)

country_change %>%
  filter(
    source == "combined",
    transition == "yr1_to_yr5"
  ) %>%
  arrange(net_change_km2) %>%
  select(country, net_change_km2) %>%
  head(10)

country_change %>%
  group_by(source, transition) %>%
  summarise(
    country_net_km2 = sum(net_change_km2),
    .groups = "drop"
  ) %>%
  left_join(
    temporal_change %>%
      select(source, transition,
             global_net_km2 = net_change_km2),
    by = c("source", "transition")
  ) %>%
  mutate(
    difference_km2 = country_net_km2 - global_net_km2
  )
