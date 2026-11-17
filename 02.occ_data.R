# Set working directory
setwd("S:/Sci-Biol/ShawanChowdhury-Lab/HonoursProjects/Sian Zahid/Data/")

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
# Cleaning GBIF occurrence data
# Import raw GBIF data
gbif_data <- read_delim("data/gbif/occurrence.txt", delim = "\t")

# Remove records with missing species or coordinates and exact duplicates
gbif_data <- gbif_data %>%
  filter(!is.na(species),
         species != "",
         !is.na(decimalLongitude),
         !is.na(decimalLatitude)) %>%
  distinct()

# Convert country codes from ISO2 to ISO3
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

# Keep records that pass all tests
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

# Load existing suitability maps for all five periods
maps <- rast(paste0(
  "output/bin_maps_full/combined_bin_yr",
  1:5, ".tif"
))

# Retain cells suitable in at least one period
suitable <- app(maps, fun = "max", na.rm = TRUE)

# Extract predicted suitability at Flickr locations
xy <- as.matrix(flickr[, c("longitude", "latitude")])
flickr$suitable <- terra::extract(suitable, xy)[[1]]

# Save records outside previously predicted suitable areas
write_csv(
  filter(flickr, is.na(suitable) | suitable != 1),
  "data/flickr/flickr_outside_suitability.csv"
)

# Retain records within predicted suitable areas
flickr <- flickr %>%
  filter(suitable == 1)

# Export screened records with image URLs
write_csv(flickr, "data/flickr/flickr_screened.csv")

# Prepare Flickr records for the combined dataset
flickr <- flickr %>%
  transmute(
    longitude,
    latitude,
    year,
    source = "Flickr"
  )

write_csv(flickr, "data/flickr/flickr_clean.csv")

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
# Combining all occurrence sources
# Clean year and remove incomplete records
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
sp_data <- bind_rows(gbif, facebook, lit_review, flickr)

# Clean year and remove incomplete records
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
# Removing duplicate occurrence records by raster cell

# Use final environmental raster as spatial template
env <- rast("data/env/env_overall.tif")

# Assign each occurrence to a raster cell
xy <- as.matrix(sp_data[, c("longitude", "latitude")])
sp_data$cell <- cellFromXY(env[[1]], xy)

# Extract environmental values
env_values <- terra::extract(env, xy)

# Keep only occurrences with complete environmental data
sp_data <- sp_data[complete.cases(env_values), ]

# Record all sources represented in each raster cell
cell_sources <- sp_data %>%
  group_by(year_group, cell) %>%
  summarise(
    sources = paste(sort(unique(source)), collapse = "; "),
    .groups = "drop"
  )

# Retain one occurrence per cell and period
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

# Assign GBIF records to environmental raster cells
gbif_xy <- as.matrix(gbif[, c("longitude", "latitude")])
gbif$cell <- terra::cellFromXY(env[[1]], gbif_xy)

# Remove records without complete environmental data
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

# Check final occurrence data
# Unique cells in each dataset
sp_data %>%
  count(year_group, name = "combined_cells")

gbif %>%
  count(year_group, name = "gbif_cells")

# Unique cells represented by each source
sp_data %>%
  tidyr::separate_rows(sources, sep = ";\\s*") %>%
  count(year_group, sources, name = "cells") %>%
  print(n = Inf)

# Unique cells represented by each source combination
sp_data %>%
  count(year_group, sources, name = "cells") %>%
  print(n = Inf)

# Total unique cells by source
sp_data %>%
  tidyr::separate_rows(sources, sep = ";\\s*") %>%
  count(sources, name = "cells")

# Facebook: 371, Flickr: 831, GBIF: 6611, lit: 146