# Load libraries
library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(rnaturalearth)
library(scales)
library(geosphere)

############################################
# SHARED SETTINGS AND DATA
############################################
periods <- paste0("yr", 1:5)
period_labels <- c(yr1 = "2001–2005", yr2 = "2006–2010", yr3 = "2011–2015",
                   yr4 = "2016–2020", yr5 = "2021–2025")
period_levels <- unname(period_labels)
short_labels <- c("2001–05", "2006–10", "2011–15", "2016–20", "2021–25")

source_cols <- c("GBIF" = "orange2", "Supplementary" = "darkorchid3")
change_cols <- c("Stable all periods" = "#16A6A3", "Gain" = "#9633CF",
                 "Loss" = "#F4A300", "Intermittent" = "#9A9A9A")
comparison_cols <- c("GBIF only" = "#FFA500", "Combined only" = "#9B21D5",
                     "Suitable in both" = "#21D5E3")
model_cols <- c("GBIF" = "#FFA500", "Combined" = "#9B21D5")

MAP_FACT <- 3L  # display aggregation only
N_PER_CATEGORY <- 6000L  # cells sampled for Fig. 3c
REBUILD_LAT3 <- TRUE  # rebuild cached samples
lat3_cache <- "output/figure/figure3_lat_elev_samples.rds"
dir.create("output/figure", recursive = TRUE, showWarnings = FALSE)

land <- rnaturalearth::ne_download(scale = 50, type = "land", category = "physical", returnclass = "sf")
occ <- read.csv("data/occ/com.csv")
xlim_map <- range(occ$longitude, na.rm = TRUE) + c(-5, 5)
ylim_map <- pmax(-90, pmin(90, range(occ$latitude, na.rm = TRUE) + c(-5, 5)))
map_extent <- terra::ext(xlim_map[1], xlim_map[2], ylim_map[1], ylim_map[2])

read_stack <- function(source) {
  files <- sprintf("output/bin_maps_full/%s_bin_%s.tif", source, periods)
  if (!all(file.exists(files))) stop("Missing SDM raster(s) for source: ", source)
  x <- terra::rast(files); names(x) <- periods
  x
}

# Figure settings
fig2_sources <- strsplit(Sys.getenv("FIG2_SOURCES", "annual5"), ",")[[1]]
fig2_annual  <- as.logical(Sys.getenv("FIG2_ANNUAL", "TRUE"))  # annual panels (b)-(d)
run_fig13 <- as.logical(Sys.getenv("RUN_FIG13", "FALSE"))
r_stack <- read_stack(if (run_fig13) "combined" else fig2_sources[1])
gbif_stack <- if (run_fig13) read_stack("gbif") else r_stack
ref <- r_stack[[1]]
if (!terra::is.lonlat(ref)) stop("The SDM raster should have longitude/latitude coordinates")
for (i in seq_along(periods)) {
  if (!terra::compareGeom(ref, r_stack[[i]], stopOnError = FALSE) ||
      !terra::compareGeom(ref, gbif_stack[[i]], stopOnError = FALSE)) stop("Raster grids do not match: ", periods[i])
}

land_v <- terra::vect(land)
if (!terra::same.crs(land_v, ref)) land_v <- terra::project(land_v, terra::crs(ref))
land_v$land_flag <- 1L
land_mask <- terra::rasterize(land_v, ref, field = "land_flag", background = 0)
cell_area <- terra::cellSize(ref, unit = "km")

# Static elevation (from 01)
elev <- terra::rast("data/env/env_yr1.tif")[["elev"]]
if (!terra::compareGeom(ref, elev, crs = FALSE, stopOnError = FALSE)) {
  elev <- terra::project(elev, ref, method = "bilinear")
}
if (!terra::compareGeom(ref, elev, crs = FALSE, stopOnError = FALSE)) stop("Elevation grid is not aligned")

# Raster to plot data
map_to_df <- function(r, codes, labels, period = NULL) {
  x <- terra::crop(r, map_extent)
  x <- terra::ifel(x %in% codes, x, NA)
  if (MAP_FACT > 1L) x <- terra::aggregate(x, fact = MAP_FACT, fun = "modal", na.rm = TRUE)
  d <- as.data.frame(x, xy = TRUE, na.rm = TRUE)
  names(d) <- c("longitude", "latitude", "code")
  d$category <- factor(d$code, levels = codes, labels = labels)
  if (!is.null(period)) d$Period <- factor(period_labels[[period]], levels = period_levels)
  d
}

wquantile <- function(x, w, p = 0.95) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]; w <- w[ok]; o <- order(x)
  x[o][which(cumsum(w[o]) >= p * sum(w))[1]]
}

panel_theme <- theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8), axis.title = element_text(size = 9),
        plot.margin = margin(6, 8, 6, 6))

if (run_fig13) {
############################################
# FIGURE 1
############################################
# Occurrences for Figure 1
gbif_occ <- read.csv("data/occ/gbif.csv") %>%
  transmute(cell, year_group, source_group = "GBIF")

supp_occ <- read.csv("data/occ/com.csv") %>%
  filter(grepl("Facebook|Literature|Flickr|Lohman", sources,
               ignore.case = TRUE)) %>%
  transmute(cell, year_group, source_group = "Supplementary")

# Combine grid-cell coverage
occ_all <- bind_rows(gbif_occ, supp_occ) %>%
  filter(!is.na(cell), !is.na(year_group))

# Use raster-cell centres for mapping
template <- terra::rast("data/env/env_overall.tif")[[1]]
xy <- terra::xyFromCell(template, occ_all$cell)

occ_all$longitude <- xy[, 1]
occ_all$latitude <- xy[, 2]

fig1_map_sf <- occ_all %>% distinct(source_group, cell, .keep_all = TRUE) %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

p1a <- ggplot() +
  geom_sf(data = land, fill = "white", colour = "grey75", linewidth = 0.25) +
  geom_sf(data = fig1_map_sf, aes(colour = source_group), alpha = 0.8, size = 0.7) +
  scale_colour_manual(values = source_cols) +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 85), expand = FALSE, datum = NA) +
  labs(colour = NULL) + theme_void() +
  theme(legend.position = "bottom", plot.margin = margin(5, 5, 5, 5))

# Source-specific shares
fig1_bar_dat <- occ_all %>% distinct(year_group, source_group, cell) %>%
  count(year_group, source_group, name = "n_cells") %>%
  group_by(year_group) %>% mutate(prop = n_cells / sum(n_cells)) %>% ungroup() %>%
  mutate(year_group = factor(period_labels[year_group], levels = period_levels),
         source_group = factor(source_group, levels = names(source_cols)))

p1b <- ggplot(fig1_bar_dat, aes(year_group, prop, fill = source_group)) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = source_cols) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Source-specific grid coverage (%)", fill = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(5, 5, 5, 5))

figure1 <- p1a + p1b + plot_layout(widths = c(1.7, 1))
print(figure1)
ggsave("output/figure/data_coverage.png", figure1, width = 12.5, height = 5.8, dpi = 400, bg = "white")
}

for (fig2_source in fig2_sources) {
fig_stack <- read_stack(fig2_source)
sfx <- if (fig2_source == "combined") "" else paste0("_", fig2_source)
############################################
# FIGURE 2
############################################
# Codes: stable, gain, loss, intermittent
n_suitable <- terra::app(fig_stack, sum, na.rm = FALSE)
change_map <- terra::ifel(n_suitable == 0, 0,
                          terra::ifel(n_suitable == 5, 1,
                                      terra::ifel(fig_stack[[1]] == 0 & fig_stack[[5]] == 1, 2,
                                                  terra::ifel(fig_stack[[1]] == 1 & fig_stack[[5]] == 0, 3, 4))))
change_map <- terra::ifel(land_mask == 1, change_map, NA)
map2 <- map_to_df(change_map, 1:4, names(change_cols))

p2a <- ggplot() +
  geom_sf(data = land, fill = "white", colour = "grey85", linewidth = 0.2) +
  geom_raster(data = map2, aes(longitude, latitude, fill = category)) +
  scale_fill_manual(values = change_cols, drop = FALSE) +
  coord_sf(xlim = xlim_map, ylim = ylim_map, expand = FALSE, datum = NA) +
  labs(fill = NULL) + theme_void(base_size = 10) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8),
        plot.margin = margin(4, 4, 4, 4))

# Area from land-masked rasters
area_dat <- bind_rows(lapply(seq_along(periods), function(i) {
  x <- terra::ifel(land_mask == 1, fig_stack[[i]], NA)
  a <- terra::expanse(x, unit = "km", byValue = TRUE)
  data.frame(Period = period_labels[[periods[i]]],
             area_mkm2 = sum(a$area[a$value == 1], na.rm = TRUE) / 1e6,
             modelled_mkm2 = sum(a$area, na.rm = TRUE) / 1e6)
})) %>% mutate(Period = factor(Period, levels = period_levels))
write.csv(area_dat, paste0("output/figure/figure2_suitable_area", sfx, ".csv"), row.names = FALSE)

p2b <- ggplot(area_dat, aes(Period, area_mkm2, group = 1)) +
  geom_line(colour = "black", linewidth = 0.7) + geom_point(colour = "black", size = 2.2) +
  scale_x_discrete(labels = short_labels) +
  labs(x = NULL, y = "Suitable habitat (million km²)") + panel_theme

get_shift_metrics <- function(period) {
  r <- fig_stack[[match(period, periods)]]
  cells <- which(terra::values(terra::ifel(land_mask == 1, r, NA), mat = FALSE) == 1)
  if (!length(cells)) stop("No suitable land cells in: ", period)
  latitude <- terra::xyFromCell(r, cells)[, 2]
  area <- terra::extract(cell_area, cells)[, 1]
  z <- terra::extract(elev, cells)[, 1]
  north <- latitude >= 0; south <- latitude < 0
  ok <- is.finite(z) & is.finite(area) & area > 0
  data.frame(period = period, Period = period_labels[[period]],
             north95_deg = wquantile(latitude[north], area[north]),
             south95_deg = wquantile(-latitude[south], area[south]),
             elev95_m = wquantile(z, area), mean_elev_m = if (any(ok)) weighted.mean(z[ok], area[ok]) else NA_real_)
}

metrics <- bind_rows(lapply(periods, get_shift_metrics)) %>%
  mutate(Period = factor(Period, levels = period_levels),
         north_shift_km = (north95_deg - first(north95_deg)) * 111.2,
         south_shift_km = (south95_deg - first(south95_deg)) * 111.2,
         elev95_shift_m = elev95_m - first(elev95_m),
         mean_elev_shift_m = mean_elev_m - first(mean_elev_m))
write.csv(metrics, paste0("output/latitudinal_elevational_shifts", sfx, ".csv"), row.names = FALSE)

latitude_plot <- metrics %>% select(Period, north_shift_km, south_shift_km) %>%
  pivot_longer(-Period, names_to = "Hemisphere", values_to = "shift_km") %>%
  mutate(Hemisphere = recode(Hemisphere, north_shift_km = "Northern", south_shift_km = "Southern"))

p2c <- ggplot(latitude_plot, aes(Period, shift_km, group = Hemisphere,
                                 linetype = Hemisphere, shape = Hemisphere)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.35) +
  geom_line(colour = "black", linewidth = 0.7) + geom_point(colour = "black", size = 2) +
  scale_linetype_manual(values = c(Northern = "solid", Southern = "dashed")) +
  scale_shape_manual(values = c(Northern = 16, Southern = 17)) +
  scale_x_discrete(labels = short_labels) +
  labs(x = NULL, y = "Poleward-edge shift (km)\n(95th percentile)", linetype = NULL, shape = NULL) +
  panel_theme + theme(legend.position = "bottom", legend.text = element_text(size = 8))

p2d <- ggplot(metrics, aes(Period, elev95_shift_m, group = 1)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.35) +
  geom_line(colour = "black", linewidth = 0.7) + geom_point(colour = "black", size = 2.2) +
  scale_x_discrete(labels = short_labels) +
  labs(x = NULL, y = "Elevational-edge shift (m)\n(95th percentile)") + panel_theme

# Annual panels (b)-(d)
if (fig2_annual) {
  sfx <- "_annual"
  ann <- read.csv("output/annual/annual_metrics.csv") %>%
    mutate(north_shift_km = (north95_deg - first(north95_deg)) * 111.2,
           south_shift_km = (south95_deg - first(south95_deg)) * 111.2,
           elev95_shift_m = elev95_m - first(elev95_m))
  year_axis <- scale_x_continuous(breaks = c(2001, 2005, 2010, 2015, 2020, 2025))

  p2b <- ggplot(ann, aes(year, area_mkm2)) +
    geom_line(colour = "black", linewidth = 0.7) + geom_point(colour = "black", size = 2.2) +
    year_axis + labs(x = NULL, y = "Suitable habitat (million km²)") + panel_theme

  ann_lat <- ann %>% select(year, north_shift_km, south_shift_km) %>%
    pivot_longer(-year, names_to = "Hemisphere", values_to = "shift_km") %>%
    mutate(Hemisphere = recode(Hemisphere, north_shift_km = "Northern", south_shift_km = "Southern"))
  # Coloured lines as Fig. 3c
  hemi_cols <- c(Northern = "#9633CF", Southern = "#F4A300")
  p2c <- ggplot(ann_lat, aes(year, shift_km, group = Hemisphere,
                             colour = Hemisphere, shape = Hemisphere)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.35) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
    scale_colour_manual(values = hemi_cols) +
    scale_shape_manual(values = c(Northern = 16, Southern = 17)) +
    year_axis + labs(x = NULL, y = "Poleward-edge shift (km)\n(95th percentile)", colour = NULL, shape = NULL) +
    panel_theme + theme(legend.position = "bottom", legend.text = element_text(size = 8))

  p2d <- ggplot(ann, aes(year, elev95_shift_m)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.35) +
    geom_line(colour = "black", linewidth = 0.7) + geom_point(colour = "black", size = 2.2) +
    year_axis + labs(x = NULL, y = "Elevational-edge shift (m)\n(95th percentile)") + panel_theme
}

figure2 <- (p2a / (p2b | p2c | p2d)) +
  plot_layout(heights = c(1.6, 1), widths = c(1.2, 1, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold", size = 12))
print(figure2)
ggsave(paste0("output/figure/figure2_temporal_range_dynamics", sfx, ".png"), figure2,
       width = 11.69, height = 7.2, dpi = 400, bg = "white")
}

if (run_fig13) {
############################################
# FIGURE 3
############################################
# 0 neither,1 GBIF,2 combined,3 both
comparison <- setNames(lapply(seq_along(periods), function(i) {
  terra::ifel(land_mask == 1, gbif_stack[[i]] + 2 * r_stack[[i]], NA)
}), periods)
map3 <- bind_rows(lapply(periods, function(p) {
  map_to_df(comparison[[p]], 1:3, names(comparison_cols), period = p)
}))

p3a <- ggplot() +
  geom_sf(data = land, fill = "white", colour = "grey90", linewidth = 0.15) +
  geom_raster(data = map3, aes(longitude, latitude, fill = category)) +
  facet_wrap(~Period, ncol = 3) +
  scale_fill_manual(values = comparison_cols, drop = FALSE) +
  coord_sf(xlim = xlim_map, ylim = ylim_map, datum = NA, expand = FALSE) +
  labs(fill = NULL) + theme_void(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 9),
        legend.position = "bottom", plot.margin = margin(2, 2, 2, 2),
        panel.spacing = grid::unit(0.25, "lines"))

############################################
# FIGURE 3a | Geographic comparison
############################################
MAP_FACT <- 2L
map3 <- bind_rows(lapply(periods, function(p) {

  r <- terra::crop(comparison[[p]], map_extent)
  r <- terra::ifel(r %in% 1:3, r, NA)

  if (MAP_FACT > 1L) r <- terra::aggregate(r, fact = MAP_FACT, fun = "modal", na.rm = TRUE)

  d <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(d) <- c("longitude", "latitude", "code")
  d$category <- factor(d$code, levels = 1:3, labels = names(comparison_cols))
  d$Period <- factor(period_labels[[p]], levels = period_levels)
  d
}))

p3a <- ggplot() +
  geom_sf(data = land, fill = "white", colour = "grey85", linewidth = 0.2) +
  geom_raster(data = map3, aes(longitude, latitude, fill = category)) +
  facet_wrap(~Period, ncol = 3) +
  scale_fill_manual(values = comparison_cols, drop = FALSE) +
  coord_sf(xlim = xlim_map, ylim = ylim_map, datum = NA, expand = FALSE) +
  labs(fill = NULL) +
  theme_void(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        legend.position = "bottom",
        legend.text = element_text(size = 9),
        panel.spacing = grid::unit(0.5, "lines"),
        plot.margin = margin(5, 5, 5, 5))

p3a <- p3a +
  theme(
    panel.spacing = grid::unit(0.1, "lines"),
    strip.text = element_text(face = "bold", size = 11),
    plot.margin = margin(0, 0, 0, 0)
  )

get_area_comparison <- function(period) {
  a <- terra::expanse(comparison[[period]], unit = "km", byValue = TRUE)
  km2 <- setNames(rep(0, 4), as.character(0:3))
  km2[as.character(a$value)] <- a$area
  shared <- unname(km2["3"]); gbif_only <- unname(km2["1"]); combined_only <- unname(km2["2"])
  union <- gbif_only + combined_only + shared
  data.frame(period = period, Period = period_labels[[period]],
             gbif_only_km2 = gbif_only, combined_only_km2 = combined_only, shared_km2 = shared,
             modelled_km2 = sum(km2), gbif_km2 = gbif_only + shared,
             combined_km2 = combined_only + shared, difference_km2 = combined_only - gbif_only,
             jaccard = if (union > 0) shared / union else NA_real_)
}

area3 <- bind_rows(lapply(periods, get_area_comparison)) %>%
  mutate(Period = factor(Period, levels = period_levels), difference_mkm2 = difference_km2 / 1e6)
write.csv(area3, "output/figure/area_comparison.csv", row.names = FALSE)
if (diff(range(area3$modelled_km2)) > 1e4) warning("Modelled land area differs among periods")

p3b <- ggplot(area3, aes(Period, difference_mkm2, group = 1)) +
  geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.4) +
  geom_line(colour = "black", linewidth = 0.75) + geom_point(colour = "black", size = 2.2) +
  scale_x_discrete(labels = short_labels) +
  labs(x = NULL, y = "Combined − GBIF (million km²)") + panel_theme

# Sample suitable cells (Fig. 3c)
sample_comparison <- function(period) {
  r <- comparison[[period]]; v <- terra::values(r, mat = FALSE)
  cells <- unlist(lapply(1:3, function(k) {
    idx <- which(v == k)
    if (length(idx) > N_PER_CATEGORY) sample(idx, N_PER_CATEGORY) else idx
  }), use.names = FALSE)
  if (!length(cells)) return(data.frame())
  xy <- terra::xyFromCell(r, cells)
  data.frame(period = period, latitude = xy[, 2],
             elevation = terra::extract(elev, cells)[, 1],
             category = factor(v[cells], levels = 1:3, labels = names(comparison_cols)),
             cell_area = terra::extract(cell_area, cells)[, 1]) %>%
    filter(is.finite(latitude), is.finite(elevation), is.finite(cell_area))
}

if (!REBUILD_LAT3 && file.exists(lat3_cache)) {
  lat3 <- readRDS(lat3_cache)
} else {
  set.seed(123)
  lat3 <- bind_rows(lapply(periods, sample_comparison))
  saveRDS(lat3, lat3_cache)
}

lat3_summary <- lat3 %>% mutate(lat_band = floor(latitude / 10) * 10 + 5) %>%
  group_by(category, lat_band) %>%
  summarise(median_elev = wquantile(elevation, cell_area, p = 0.5), n = n(), .groups = "drop") %>%
  filter(n >= 20)
write.csv(lat3_summary, "output/figure/latitude_elevation_summary.csv", row.names = FALSE)

p3c <- ggplot(lat3_summary, aes(lat_band, median_elev, colour = category, group = category)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
  scale_colour_manual(values = comparison_cols, guide = "none") +
  labs(x = "Latitude (°)", y = "Median elevation (m)") + panel_theme

# Area-weighted centroids
get_centroid <- function(i, source) {
  r <- if (source == "GBIF") gbif_stack[[i]] else r_stack[[i]]
  cells <- which(terra::values(terra::ifel(land_mask == 1, r, NA), mat = FALSE) == 1)
  if (!length(cells)) stop("No predicted suitable cells: ", periods[i], "/", source)
  xy <- terra::xyFromCell(r, cells) * pi / 180
  w <- terra::extract(cell_area, cells)[, 1]
  x <- sum(w * cos(xy[, 2]) * cos(xy[, 1]))
  y <- sum(w * cos(xy[, 2]) * sin(xy[, 1]))
  z <- sum(w * sin(xy[, 2]))
  data.frame(period = periods[i], Period = period_labels[[periods[i]]], Source = source,
             lon = atan2(y, x) * 180 / pi, lat = atan2(z, sqrt(x^2 + y^2)) * 180 / pi)
}

centroids3 <- bind_rows(lapply(names(model_cols), function(source) {
  d <- bind_rows(lapply(seq_along(periods), get_centroid, source = source)) %>%
    mutate(Period = factor(Period, levels = period_levels))
  base <- as.numeric(d[1, c("lon", "lat")])
  d$centroid_shift_km <- vapply(seq_len(nrow(d)), function(i) {
    geosphere::distGeo(base, as.numeric(d[i, c("lon", "lat")])) / 1000
  }, numeric(1))
  d
}))
write.csv(centroids3, "output/figure/centroid_displacement.csv", row.names = FALSE)

p3d <- ggplot(centroids3, aes(Period, centroid_shift_km, group = Source,
                              colour = Source, shape = Source)) +
  geom_line(linewidth = 0.75) + geom_point(size = 2) +
  scale_colour_manual(values = model_cols) +
  scale_shape_manual(values = c(GBIF = 16, Combined = 17)) +
  scale_x_discrete(labels = short_labels) +
  labs(x = NULL, y = "Centroid displacement (km)", colour = NULL, shape = NULL) +
  panel_theme + theme(legend.position = "bottom", legend.text = element_text(size = 8))

# Export plot
figure3 <- (p3a / (p3b | p3c | p3d)) +
  plot_layout(heights = c(2.5, 1), guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = "bottom",
        plot.tag = element_text(face = "bold", size = 12))

print(figure3)

ggsave("output/figure/figure3_data_source_comparison.png",
       figure3, width = 11.69, height = 8.27, dpi = 400, bg = "white")
}

############################################
# Annual Figures 1 and 3
############################################
fig_annual_all <- as.logical(Sys.getenv("FIG13_ANNUAL", "TRUE"))
if (fig_annual_all) {
years_all <- 2001:2025

# Figure 1: yearly coverage
com_occ <- read.csv("data/occ/com.csv")
occ_y <- bind_rows(
  com_occ %>% filter(grepl("GBIF", sources)) %>% transmute(cell, year, source_group = "GBIF"),
  com_occ %>% filter(grepl("Facebook|Literature|Flickr|Lohman", sources, ignore.case = TRUE)) %>%
    transmute(cell, year, source_group = "Supplementary"))
template <- terra::rast("data/env/env_overall.tif")[[1]]
xy <- terra::xyFromCell(template, occ_y$cell)
occ_y$longitude <- xy[, 1]; occ_y$latitude <- xy[, 2]
fig1_map_sf <- occ_y %>% distinct(source_group, cell, .keep_all = TRUE) %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
p1a <- ggplot() +
  geom_sf(data = land, fill = "white", colour = "grey75", linewidth = 0.25) +
  geom_sf(data = fig1_map_sf, aes(colour = source_group), alpha = 0.8, size = 0.7) +
  scale_colour_manual(values = source_cols) +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 85), expand = FALSE, datum = NA) +
  labs(colour = NULL) + theme_void() +
  theme(legend.position = "bottom", plot.margin = margin(5, 5, 5, 5))
fig1_bar_y <- occ_y %>% distinct(year, source_group, cell) %>%
  count(year, source_group, name = "n_cells") %>%
  group_by(year) %>% mutate(prop = n_cells / sum(n_cells)) %>% ungroup() %>%
  mutate(source_group = factor(source_group, levels = names(source_cols)))
write.csv(fig1_bar_y, "output/figure/figure1_coverage_annual.csv", row.names = FALSE)
p1b <- ggplot(fig1_bar_y, aes(factor(year), prop, fill = source_group)) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = source_cols) +
  scale_x_discrete(breaks = as.character(seq(2001, 2025, 4)), expand = expansion(add = c(2, 0.6))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Source-specific grid coverage (%)", fill = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(5, 5, 5, 5))
figure1 <- p1a + p1b + plot_layout(widths = c(1.7, 1))
ggsave("output/figure/data_coverage_annual.png", figure1, width = 12.5, height = 5.8, dpi = 400, bg = "white")

# Figure 3: GBIF vs combined
read_year <- function(ds, y) {
  r <- terra::rast(sprintf("output/sdm/%s/%s_bin_%d.tif", ds, ds, y))
  r <- terra::extend(r, ref); r <- terra::resample(r, ref, method = "near")
  terra::ifel(is.na(r), 0, r)
}
cmp_year <- setNames(lapply(years_all, function(y) {
  terra::ifel(land_mask == 1, read_year("annual_gbif", y) + 2 * read_year("annual", y), NA)
}), years_all)

# (a) five-year windows
win <- split(years_all, rep(periods, each = 5))
comparison <- setNames(lapply(periods, function(p) {
  g <- sum(terra::rast(lapply(as.character(win[[p]]), function(y) cmp_year[[y]] %in% c(1, 3)))) >= 3
  c <- sum(terra::rast(lapply(as.character(win[[p]]), function(y) cmp_year[[y]] %in% c(2, 3)))) >= 3
  terra::ifel(land_mask == 1, g + 2 * c, NA)
}), periods)
MAP_FACT <- 2L
map3 <- bind_rows(lapply(periods, function(p) map_to_df(comparison[[p]], 1:3, names(comparison_cols), period = p)))
p3a <- ggplot() +
  geom_sf(data = land, fill = "white", colour = "grey85", linewidth = 0.2) +
  geom_raster(data = map3, aes(longitude, latitude, fill = category)) +
  facet_wrap(~Period, ncol = 3) +
  scale_fill_manual(values = comparison_cols, drop = FALSE) +
  coord_sf(xlim = xlim_map, ylim = ylim_map, datum = NA, expand = FALSE) +
  labs(fill = NULL) + theme_void(base_size = 11) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 11),
        legend.position = "bottom", legend.text = element_text(size = 9),
        panel.spacing = grid::unit(0.1, "lines"), plot.margin = margin(0, 0, 0, 0))

# (b) yearly area difference
area3y <- bind_rows(lapply(years_all, function(y) {
  a <- terra::expanse(cmp_year[[as.character(y)]], unit = "km", byValue = TRUE)
  km2 <- setNames(rep(0, 4), as.character(0:3)); km2[as.character(a$value)] <- a$area
  u <- km2["1"] + km2["2"] + km2["3"]
  data.frame(year = y, gbif_only_km2 = km2["1"], combined_only_km2 = km2["2"], shared_km2 = km2["3"],
             gbif_km2 = km2["1"] + km2["3"], combined_km2 = km2["2"] + km2["3"],
             difference_mkm2 = (km2["2"] - km2["1"]) / 1e6, jaccard = km2["3"] / u)
}))
write.csv(area3y, "output/figure/area_comparison_annual.csv", row.names = FALSE)
yr_x <- scale_x_continuous(breaks = seq(2001, 2025, 4))
p3b <- ggplot(area3y, aes(year, difference_mkm2)) +
  geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.4) +
  geom_line(colour = "black", linewidth = 0.75) + geom_point(colour = "black", size = 2.2) +
  yr_x + labs(x = NULL, y = "Combined − GBIF (million km²)") + panel_theme

# (c) latitude-elevation profile
set.seed(123)
lat3 <- bind_rows(lapply(periods, sample_comparison <- function(period) {
  r <- comparison[[period]]; v <- terra::values(r, mat = FALSE)
  cells <- unlist(lapply(1:3, function(k) { idx <- which(v == k)
    if (length(idx) > N_PER_CATEGORY) sample(idx, N_PER_CATEGORY) else idx }), use.names = FALSE)
  if (!length(cells)) return(data.frame())
  xy <- terra::xyFromCell(r, cells)
  data.frame(period = period, latitude = xy[, 2], elevation = terra::extract(elev, cells)[, 1],
             category = factor(v[cells], levels = 1:3, labels = names(comparison_cols)),
             cell_area = terra::extract(cell_area, cells)[, 1]) %>%
    filter(is.finite(latitude), is.finite(elevation), is.finite(cell_area))
}))
lat3_summary <- lat3 %>% mutate(lat_band = floor(latitude / 10) * 10 + 5) %>%
  group_by(category, lat_band) %>%
  summarise(median_elev = wquantile(elevation, cell_area, p = 0.5), n = n(), .groups = "drop") %>%
  filter(n >= 20)
write.csv(lat3_summary, "output/figure/latitude_elevation_summary_annual.csv", row.names = FALSE)
p3c <- ggplot(lat3_summary, aes(lat_band, median_elev, colour = category, group = category)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
  scale_colour_manual(values = comparison_cols, guide = "none") +
  labs(x = "Latitude (°)", y = "Median elevation (m)") + panel_theme

# (d) yearly centroid displacement
cent_y <- function(r) {
  cells <- which(terra::values(r, mat = FALSE) == 1)
  xy <- terra::xyFromCell(r, cells) * pi / 180; w <- terra::extract(cell_area, cells)[, 1]
  x <- sum(w * cos(xy[, 2]) * cos(xy[, 1])); y <- sum(w * cos(xy[, 2]) * sin(xy[, 1])); z <- sum(w * sin(xy[, 2]))
  c(atan2(y, x), atan2(z, sqrt(x^2 + y^2))) * 180 / pi
}
centroids3y <- bind_rows(lapply(c("GBIF", "Combined"), function(s) {
  m <- t(sapply(as.character(years_all), function(y) {
    v <- cmp_year[[y]]; cent_y(if (s == "GBIF") v %in% c(1, 3) else v %in% c(2, 3)) }))
  data.frame(year = years_all, Source = s, lon = m[, 1], lat = m[, 2],
             centroid_shift_km = geosphere::distGeo(m[1, ], m) / 1000)
}))
write.csv(centroids3y, "output/figure/centroid_displacement_annual.csv", row.names = FALSE)
# (d) GBIF-combined centroid distance
cdist <- centroids3y %>% select(year, Source, lon, lat) %>%
  tidyr::pivot_wider(names_from = Source, values_from = c(lon, lat)) %>%
  mutate(centroid_distance_km = geosphere::distGeo(cbind(lon_GBIF, lat_GBIF), cbind(lon_Combined, lat_Combined)) / 1000)
write.csv(cdist, "output/figure/centroid_distance_annual.csv", row.names = FALSE)
p3d <- ggplot(cdist, aes(year, centroid_distance_km)) +
  geom_line(colour = "black", linewidth = 0.75) + geom_point(colour = "black", size = 2.2) +
  scale_y_continuous(limits = c(0, NA)) +
  yr_x + labs(x = NULL, y = "Centroid displacement (km)") + panel_theme

figure3 <- (p3a / (p3b | p3c | p3d)) +
  plot_layout(heights = c(2.5, 1), guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = "bottom", plot.tag = element_text(face = "bold", size = 12))
ggsave("output/figure/figure3_data_source_comparison_annual.png",
       figure3, width = 11.69, height = 8.27, dpi = 400, bg = "white")
}

############################################
# GCB export (180 mm)
############################################
# Keep original proportions
dir.create("output/figure/gcb", showWarnings = FALSE)
gcb_save <- function(p, name, h_mm, scaling = 0.8) {
  ggsave(file.path("output/figure/gcb", paste0(name, ".tif")), p, width = 180, height = h_mm, units = "mm",
         dpi = 600, device = ragg::agg_tiff, compression = "lzw", scaling = scaling, bg = "white")
}
tag_theme <- theme(plot.tag = element_text(face = "bold", size = 12))

if (fig_annual_all) {
  g1a <- p1a + guides(colour = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
    theme(legend.position = "top", legend.justification = "right",
          legend.direction = "horizontal", legend.text = element_text(size = 11))
  g1b <- p1b + theme(legend.position = "none")
  gcb_save((free(g1a) + g1b + plot_layout(widths = c(1.7, 1))) +
             plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") & tag_theme,
           "Figure_1", 78)
}
gcb_save(figure2, "Figure_2", 125)
if (fig_annual_all) {
  top_room <- theme(plot.margin = margin(16, 8, 6, 6))
  # One combined legend
  draw_key_sqline <- function(data, params, size) {
    col <- data$colour
    sq <- grid::unit(1, "snpc")  # Square as other legends
    grid::grobTree(
      grid::rectGrob(x = sq * 0.5, width = sq, height = sq, gp = grid::gpar(fill = col, col = NA)),
      grid::segmentsGrob(sq + grid::unit(3, "pt"), grid::unit(0.5, "npc"), grid::unit(1, "npc"), grid::unit(0.5, "npc"),
                         gp = grid::gpar(col = col, lwd = 0.7 * .pt)),
      grid::pointsGrob((sq + grid::unit(3, "pt") + grid::unit(1, "npc")) * 0.5, grid::unit(0.5, "npc"), pch = 19,
                       size = grid::unit(1.4 * .pt, "pt"), gp = grid::gpar(col = col)))
  }
  g3c <- ggplot(lat3_summary, aes(lat_band, median_elev, colour = category, group = category)) +
    geom_line(linewidth = 0.7, key_glyph = draw_key_sqline) + geom_point(size = 1.4, show.legend = FALSE) +
    scale_colour_manual(values = comparison_cols, name = NULL, drop = FALSE) +
    labs(x = "Latitude (°)", y = "Median elevation (m)") + panel_theme + top_room +
    theme(legend.key.width = grid::unit(3.2, "lines"), legend.key.height = grid::unit(1.2, "lines"), legend.text = element_text(size = 9))
  g3 <- ((p3a + guides(fill = "none")) /
         ((p3b + labs(y = "Combined − GBIF\n(million km²)") + top_room) | g3c | (p3d + labs(y = "Centroid\ndisplacement (km)") + top_room))) +
    plot_layout(heights = c(1.9, 1), guides = "collect") +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    theme(legend.position = "bottom", legend.justification = "center",
          plot.tag = element_text(face = "bold", size = 12))
  gcb_save(g3, "Figure_3", 118)
}
