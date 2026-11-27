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
fig2_sources <- strsplit(Sys.getenv("FIG2_SOURCES", "combined,unified"), ",")[[1]]
run_fig13 <- as.logical(Sys.getenv("RUN_FIG13", "TRUE"))
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
