# Run from project root

# Loading required libraries
library(terra)
library(dplyr)
library(readr)
library(countrycode)
library(CoordinateCleaner)
library(dplyr)
library(tidyr)
library(sf)

#########################################
# Clean GBIF data
gbif_data <- read_delim("data/gbif/occurrence.txt", delim = "\t")

# Keep precise coordinates
coord_txt <- read_delim("data/gbif/occurrence.txt", delim = "\t",
                        col_select = c(gbifID, decimalLatitude, decimalLongitude),
                        col_types = cols(.default = col_character()))
n_dec <- function(x) ifelse(grepl(".", x, fixed = TRUE), nchar(sub("^[^.]*[.]", "", x)), 0L)
coord_txt <- coord_txt %>%
  transmute(gbifID = as.character(gbifID),
            coord_decimals = pmin(n_dec(decimalLatitude), n_dec(decimalLongitude)))

gbif_data <- gbif_data %>%
  mutate(gbifID = as.character(gbifID)) %>%
  left_join(coord_txt, by = "gbifID") %>%
  filter((!is.na(coordinateUncertaintyInMeters) & coordinateUncertaintyInMeters < 10000) |
           (is.na(coordinateUncertaintyInMeters) & coord_decimals >= 3)) %>%
  # Drop captive or escaped records
  filter(basisOfRecord != "LIVING_SPECIMEN",
         countryCode != "NL" | is.na(countryCode)) %>%
  dplyr::select(-coord_decimals)

# Remove incomplete and duplicates
gbif_data <- gbif_data %>%
  filter(!is.na(species),
         species != "",
         !is.na(decimalLongitude),
         !is.na(decimalLatitude)) %>%
  distinct()

# ISO2 to ISO3 codes
gbif_data$countryCode <- countrycode(gbif_data$countryCode,
                                     origin = "iso2c",
                                     destination = "iso3c")

# Flag problematic coordinates
flags <- clean_coordinates(
  x = as.data.frame(gbif_data),
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  countries = "countryCode",
  species = "species",
  tests = c("centroids", "gbif", "zeros", "countries")
)

# Keep passing records
gbif_cl <- gbif_data[flags$.summary, ]

# Standardise GBIF columns
gbif_cl <- gbif_cl %>%
  transmute(
    longitude = decimalLongitude,
    latitude = decimalLatitude,
    year = year,
    source = "GBIF"
  )

# Export cleaned GBIF data
write_csv(gbif_cl, "data/gbif/gbif_clean.csv")

#########################################
# Keep screened Flickr records
rescreen_flickr <- FALSE

if (rescreen_flickr) {
# Import and clean Flickr records
flickr <- read_csv("data/flickr/flickr_data.csv") %>%
  filter(trimws(bioclip_label) == "Papilio demoleus") %>%
  mutate(year = as.integer(substr(date, 1, 4))) %>%
  separate(
    geotag,
    into = c("latitude", "longitude"),
    sep = ",",
    convert = TRUE
  ) %>%
  filter(
    between(year, 2001, 2025),
    between(latitude, -90, 90),
    between(longitude, -180, 180)
  ) %>%
  distinct(image_url, .keep_all = TRUE)

# Load V1 suitability maps
maps <- rast(paste0(
  "output_V1/bin_maps_full/combined_bin_yr",
  1:5, ".tif"
))

# Cells ever suitable
suitable <- app(maps, fun = "max", na.rm = TRUE)

# Suitability at Flickr records
xy <- as.matrix(flickr[, c("longitude", "latitude")])
flickr$suitable <- terra::extract(suitable, xy)[[1]]

# Save records outside suitability
write_csv(
  filter(flickr, is.na(suitable) | suitable != 1),
  "data/flickr/flickr_outside_suitability.csv"
)

# Keep records in suitable areas
flickr <- flickr %>%
  filter(suitable == 1)

# Export screened records
write_csv(flickr, "data/flickr/flickr_screened.csv")

# Prepare Flickr records
flickr <- flickr %>%
  transmute(
    longitude,
    latitude,
    year,
    source = "Flickr"
  )

write_csv(flickr, "data/flickr/flickr_clean.csv")
}

#########################################
# Import clean datasets
gbif <- read_csv("data/gbif/gbif_clean.csv")
facebook <- read_csv("data/facebook/facebook.csv") %>%
  dplyr::select(lon, lat, year) %>%
  dplyr::rename(
    longitude = lon,
    latitude = lat
  ) %>%
  mutate(source = "Facebook")

lit_review <- read_csv("data/lit_review/lit_review.csv") %>%
  dplyr::select(lon, lat, year) %>%
  dplyr::rename(
    longitude = lon,
    latitude = lat
  ) %>%
  mutate(source = "Literature")
flickr <- read_csv("data/flickr/flickr_clean.csv")

#########################################
# Lohman et al. 2025 records
lohman <- read_csv(
  "data_add/tropical_asian_butterflies/occ_tropical_asian_butterflies_v1.1.csv",
  col_types = cols(.default = col_character())
) %>%
  filter(lohman_final_genus_species == "Papilio demoleus" |
           source_genus_species == "Papilio demoleus",
         is.na(gbifID) | gbifID %in% c("", "NA"),
         is.na(locality_type) | locality_type != "district") %>%
  transmute(
    longitude = as.numeric(decimalLon),
    latitude = as.numeric(decimalLat),
    year = suppressWarnings(as.integer(as.numeric(year))),
    source = "Lohman et al. 2025",
    reference = reference
  ) %>%
  filter(!is.na(longitude), !is.na(latitude), !is.na(year),
         between(year, 2001, 2025)) %>%
  # Fix longitude typo (Cambodia)
  mutate(longitude = if_else(abs(longitude - 7.0469) < 1e-6 & abs(latitude - 13.7475) < 1e-6,
                             107.0469, longitude)) %>%
  distinct(longitude, latitude, year, .keep_all = TRUE)

# Remove duplicates of existing records
existing_keys <- bind_rows(mutate(facebook, year = as.integer(year)),
                           mutate(lit_review, year = suppressWarnings(as.integer(year)))) %>%
  transmute(key = paste(round(longitude, 4), round(latitude, 4), year))
lohman <- lohman %>%
  filter(!paste(round(longitude, 4), round(latitude, 4), year) %in% existing_keys$key)

dir.create("data/lohman", showWarnings = FALSE)
write_csv(lohman, "data/lohman/lohman_clean.csv")
lohman <- dplyr::select(lohman, -reference)

#########################################
# Combine all occurrence sources
gbif <- gbif %>%
  mutate(year = as.integer(year)) %>%
  filter(!is.na(longitude),
         !is.na(latitude),
         !is.na(year))

facebook <- facebook %>%
  mutate(year = as.integer(year)) %>%
  filter(!is.na(longitude),
         !is.na(latitude),
         !is.na(year))

lohman <- lohman %>%
  mutate(year = as.integer(year))

lit_review <- lit_review %>%
  mutate(year = as.integer(year)) %>%
  filter(!is.na(longitude),
         !is.na(latitude),
         !is.na(year))

flickr <- flickr %>%
  mutate(year = as.integer(year)) %>%
  filter(!is.na(longitude),
         !is.na(latitude),
         !is.na(year))

# Combine datasets
sp_data <- bind_rows(gbif, facebook, lit_review, flickr, lohman)

# Remove incomplete records
sp_data <- sp_data %>%
  mutate(year = as.integer(year)) %>%
  filter(!is.na(longitude),
         !is.na(latitude),
         !is.na(year),
         year >= 2001,
         year <= 2025) %>%
  distinct()

#########################################
# Assigning records to 5-year intervals
sp_data <- sp_data %>%
  mutate(year_group = case_when(
    year >= 2001 & year <= 2005 ~ "yr1",
    year >= 2006 & year <= 2010 ~ "yr2",
    year >= 2011 & year <= 2015 ~ "yr3",
    year >= 2016 & year <= 2020 ~ "yr4",
    year >= 2021 & year <= 2025 ~ "yr5"
  ))

# Export data
write_csv(sp_data, "data/com_data.csv")

#########################################
# Deduplicate by raster cell

# Spatial template
env <- rast("data/env/env_overall.tif")

# Assign occurrences to cells
xy <- as.matrix(sp_data[, c("longitude", "latitude")])
sp_data$cell <- cellFromXY(env[[1]], xy)

# Extract environmental values
env_values <- terra::extract(env, xy)

# Complete occurrences only
sp_data <- sp_data[complete.cases(env_values), ]

# Sources per cell
cell_sources <- sp_data %>%
  group_by(year_group, cell) %>%
  summarise(
    sources = paste(sort(unique(source)), collapse = "; "),
    .groups = "drop"
  )

# One record per cell-period
sp_data <- sp_data %>%
  distinct(year_group, cell, .keep_all = TRUE) %>%
  left_join(cell_sources, by = c("year_group", "cell"))

#########################################
# Export combined occurrence datasets
periods <- paste0("yr", 1:5)

for (period in periods) {
  occ_period <- sp_data %>%
    filter(year_group == period)

  write_csv(occ_period, paste0("data/occ/", period, "_com.csv"))
}

write_csv(sp_data, "data/occ/com.csv")

# Create GBIF-only datasets
gbif <- gbif %>%
  mutate(
    year = as.integer(year),
    year_group = case_when(
      between(year, 2001, 2005) ~ "yr1",
      between(year, 2006, 2010) ~ "yr2",
      between(year, 2011, 2015) ~ "yr3",
      between(year, 2016, 2020) ~ "yr4",
      between(year, 2021, 2025) ~ "yr5"
    )
  ) %>%
  filter(!is.na(year_group))

# Assign GBIF records to cells
gbif_xy <- as.matrix(gbif[, c("longitude", "latitude")])
gbif$cell <- terra::cellFromXY(env[[1]], gbif_xy)

# Remove incomplete records
gbif_env <- terra::extract(env, gbif_xy)
gbif <- gbif[complete.cases(gbif_env), ] %>%
  filter(!is.na(cell)) %>%
  distinct(year_group, cell, .keep_all = TRUE)

# Export GBIF-only datasets
for (period in periods) {
  gbif_period <- gbif %>%
    filter(year_group == period)

  write_csv(gbif_period, paste0("data/occ/", period, "_gbif.csv"))
}

write_csv(gbif, "data/occ/gbif.csv")

######################################
# Check final occurrence data
com_grid <- read_csv("data/occ/com.csv")
gbif_grid <- read_csv("data/occ/gbif.csv")
occ_data <- read_csv("data/com_data.csv")

###############
# Total occurrence records by period
occ_data %>%
  count(year_group, name = "combined_records")

# Total occurrence records by source
occ_data %>%
  count(source, name = "records")

# Total records
nrow(occ_data) # 14646

# Contribution of supplementary occurrence records
occ_data %>%
  mutate(source_group = if_else(
    source == "GBIF", "GBIF", "Supplementary"
  )) %>%
  count(year_group, source_group) %>%
  group_by(year_group) %>%
  mutate(percentage = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  print(n = Inf)

###############
# Total unique cells
NROW(unique(com_grid$cell)) # 6565

# Unique cells in each dataset
com_grid %>%
  count(year_group, name = "combined_cells")

gbif_grid %>%
  count(year_group, name = "gbif_cells")

# Unique cells by source
com_grid %>%
  tidyr::separate_rows(sources, sep = ";\\s*") %>%
  count(year_group, sources, name = "cells") %>%
  print(n = Inf)

# Unique cells by source combination
com_grid %>%
  count(year_group, sources, name = "cells") %>%
  print(n = Inf)

# Total cell-period observations by source
com_grid %>%
  tidyr::separate_rows(sources, sep = ";\\s*") %>%
  count(sources, name = "cell_periods")

# Unique cells by source
com_grid %>%
  tidyr::separate_rows(sources, sep = ";\\s*") %>%
  distinct(sources, cell) %>%
  count(sources, name = "unique_cells")
