# Set working directory
setwd("S:/Sci-Biol/ShawanChowdhury-Lab/HonoursProjects/Sian Zahid/Data/")

# Loading required libraries
library(terra)
library(dplyr)
library(readr)
library(countrycode)
library(CoordinateCleaner)


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
gbif <- gbif_cl %>%
  transmute(
    longitude = decimalLongitude,
    latitude = decimalLatitude,
    year = year,
    source = "GBIF"
  )

# Export cleaned GBIF data
write_csv(gbif_cl, "data/gbif/gbif_clean.csv")

#########################################
# Importing Facebook and literature records
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

# Combine datasets
sp_data <- bind_rows(gbif, facebook, lit_review)

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

# Remove duplicate cells within each 5-year interval
sp_data <- sp_data %>%
  distinct(year_group, cell, .keep_all = TRUE)

#########################################
# Exporting combined occurrence data by interval
periods <- paste0("yr", 1:5)

for (period in periods) {
  
  occ_period <- sp_data %>%
    filter(year_group == period)
  
  write_csv(occ_period,
            paste0("data/occ/", period, "_com.csv"))
}


#########################################
# Exporting combined occurrence dataset
write_csv(sp_data, "data/occ/com.csv")

#########################################
# Creating GBIF-only datasets

# Assign raster cells to cleaned GBIF records separately
gbif_xy <- as.matrix(gbif[, c("longitude", "latitude")])
gbif$cell <- cellFromXY(env[[1]], gbif_xy)

# Extract environmental values
gbif_env <- terra::extract(env, gbif_xy)

# Keep records with complete environmental data
gbif <- gbif[complete.cases(gbif_env), ]

# Assign year intervals
gbif <- gbif %>%
  mutate(
    year = as.integer(year),
    year_group = case_when(
      year >= 2001 & year <= 2005 ~ "yr1",
      year >= 2006 & year <= 2010 ~ "yr2",
      year >= 2011 & year <= 2015 ~ "yr3",
      year >= 2016 & year <= 2020 ~ "yr4",
      year >= 2021 & year <= 2025 ~ "yr5"
    )
  ) %>%
  filter(!is.na(year_group)) %>%
  distinct(year_group, cell, .keep_all = TRUE)

# Export GBIF datasets by interval
for (period in periods) {
  
  gbif_period <- gbif %>%
    filter(year_group == period)
  
  write_csv(gbif_period,
            paste0("data/occ/", period, "_gbif.csv"))
}

# Export all cleaned GBIF records
write_csv(gbif, "data/occ/gbif.csv")


#########################################
# Checking final occurrence data

# Number of combined records per interval
sp_data %>%
  count(year_group)

# Number of records by interval and source
sp_data %>%
  count(year_group, source)

# Total records by source
sp_data %>%
  count(source)
