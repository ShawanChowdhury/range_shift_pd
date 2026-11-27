# Load libraries
library(terra)
library(ade4)
library(ecospat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

set.seed(123)
# Yearly niche comparison
periods <- as.character(2001:2025)
period_labels <- setNames(periods, periods)
MIN_REC <- 5  # Years below this are skipped
period_levels <- unname(period_labels)
vars <- c("PDSI", "ppt", "soil", "tmin", "ws", "elev")
categories <- c("GBIF only", "Supplementary only", "Both")
thresholds <- c(0.02, 0.05, 0.10)
REBUILD_GRIDS <- as.logical(Sys.getenv("REBUILD_GRIDS", "TRUE"))  # rebuild cached grids

###############################################
# ENVIRONMENTAL BACKGROUND AND SHARED PCA
env_list <- setNames(lapply(periods, function(p) {
  f <- paste0("data/env_annual/env_", p, ".tif")
  if (!file.exists(f)) stop("Missing environmental raster: ", f)
  r <- terra::rast(f)
  if (!all(vars %in% names(r))) stop("Missing environmental predictors in: ", f)
  r[[vars]]
}), periods)

bg <- read.csv("data/bg.csv")
bg_xy <- as.matrix(bg[, c("longitude", "latitude")])
background <- setNames(lapply(periods, function(p) {
  z <- terra::extract(env_list[[p]], bg_xy)
  z <- z[complete.cases(z[, vars]), vars, drop = FALSE]
  message("Background ", p, ": ", nrow(z), " valid points")
  z
}), periods)

pooled_background <- dplyr::bind_rows(background)
pca_shared <- ade4::dudi.pca(pooled_background, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)
pooled_scores <- as.data.frame(pca_shared$li[, 1:2, drop = FALSE])
background_scores <- lapply(background, function(z) {
  as.data.frame(ade4::suprow(pca_shared, z)$lisup[, 1:2, drop = FALSE])
})
pc_var <- 100 * pca_shared$eig[1:2] / sum(pca_shared$eig)

write.csv(data.frame(variable = rownames(pca_shared$co), PC1 = pca_shared$co[, 1], PC2 = pca_shared$co[, 2]),
          "output/niche_annual_PCA_loadings.csv", row.names = FALSE)
write.csv(data.frame(axis = c("PC1", "PC2"), variance_percent = as.numeric(pc_var)),
          "output/niche_annual_PCA_variance.csv", row.names = FALSE)

###############################################
# OCCURRENCE SCORES AND NICHE GRIDS
com_occ <- read.csv("data/occ/com.csv")
get_scores <- function(period, source) {
  d <- com_occ[com_occ$year == as.integer(period), ]
  # Cells by source
  if (source == "gbif") d <- d[grepl("GBIF", d$sources), ]
  if (source == "supp") d <- d[grepl("Facebook|Literature|Flickr|Lohman", d$sources, ignore.case = TRUE), ]
  if (nrow(d) < MIN_REC) return(NULL)

  z <- terra::extract(
    env_list[[period]],
    as.matrix(d[, c("longitude", "latitude")])
  )

  z <- z[complete.cases(z[, vars]), vars, drop = FALSE]

  if (nrow(z) < MIN_REC) return(NULL)

  as.data.frame(
    ade4::suprow(pca_shared, z)$lisup[, 1:2, drop = FALSE]
  )
}

occ_scores <- setNames(lapply(periods, function(p) {
  setNames(lapply(c("gbif", "supp"), function(s) get_scores(p, s)), c("gbif", "supp"))
}), periods)

record_counts <- dplyr::bind_rows(lapply(periods, function(p) {
  data.frame(period = period_labels[[p]], source = c("GBIF", "Supplementary"),
             n = vapply(occ_scores[[p]], function(z) if (is.null(z)) 0L else nrow(z), integer(1)))
}))
write.csv(record_counts, "output/niche_annual_record_counts.csv", row.names = FALSE)

get_grid <- function(period, source) {
  sp <- occ_scores[[period]][[source]]
  if (is.null(sp)) return(NULL)
  lim <- apply(pooled_scores, 2, range)
  outside <- sp[, 1] < lim[1, 1] | sp[, 1] > lim[2, 1] |
             sp[, 2] < lim[1, 2] | sp[, 2] > lim[2, 2]
  message(period, " | ", source, " | Records: ", nrow(sp), " | Outside background: ", sum(outside))
  ecospat::ecospat.grid.clim.dyn(glob = pooled_scores, glob1 = background_scores[[period]],
                                 sp = sp, R = 300, th.sp = NULL, th.env = NULL)
}

read_or_make_grids <- function(source, file) {
  if (!REBUILD_GRIDS && file.exists(file)) {
    message("Loading saved grids: ", file, " (only valid if the shared PCA is unchanged)")
    z <- readRDS(file)
    if (!identical(names(z), periods)) stop("Saved grid periods do not match: ", file)
    return(z)
  }
  z <- setNames(lapply(periods, get_grid, source = source), periods)
  saveRDS(z, file)
  z
}

grid_gbif <- read_or_make_grids("gbif", "output/niche_annual_grids_gbif.rds")
grid_supp <- read_or_make_grids("supp", "output/niche_annual_grids_supplementary.rds")
ok_years <- periods[!vapply(grid_gbif, is.null, TRUE) & !vapply(grid_supp, is.null, TRUE)]
message("Years used: ", paste(ok_years, collapse = ", "))

###############################################
# Rebuild niche grids

# ENVIRONMENTAL COVERAGE AND SENSITIVITY
classify_env_space <- function(gbif_grid, supp_grid, period, thr) {
  a <- gbif_grid$z.uncor; b <- supp_grid$z.uncor
  if (!terra::compareGeom(a, b, stopOnError = FALSE)) stop("Niche grids are not aligned: ", period)
  g <- terra::values(a, mat = FALSE); s <- terra::values(b, mat = FALSE)
  gmax <- max(g, na.rm = TRUE); smax <- max(s, na.rm = TRUE)
  if (!is.finite(gmax) || !is.finite(smax) || gmax <= 0 || smax <= 0) stop("Empty niche density: ", period)
  gp <- is.finite(g) & g >= thr * gmax
  sp <- is.finite(s) & s >= thr * smax
  cls <- rep(NA_character_, length(g))
  cls[gp & !sp] <- "GBIF only"; cls[!gp & sp] <- "Supplementary only"; cls[gp & sp] <- "Both"
  cells <- which(!is.na(cls)); xy <- terra::xyFromCell(a, cells)
  data.frame(PC1 = xy[, 1], PC2 = xy[, 2], period = factor(period_labels[[period]], levels = period_levels),
             category = factor(cls[cells], levels = categories))
}

coverage_by_threshold <- lapply(thresholds, function(th) {
  dplyr::bind_rows(lapply(ok_years, function(p) classify_env_space(grid_gbif[[p]], grid_supp[[p]], p, th)))
})
names(coverage_by_threshold) <- as.character(thresholds)
env_class_df <- coverage_by_threshold[[as.character(0.05)]]

coverage_summary <- dplyr::bind_rows(lapply(seq_along(thresholds), function(i) {
  tab <- coverage_by_threshold[[i]] %>% count(period, category, name = "cells")
  tidyr::expand_grid(period = factor(period_levels, levels = period_levels),
                     category = factor(categories, levels = categories)) %>%
    filter(period %in% ok_years) %>%
    left_join(tab, by = c("period", "category")) %>%
    mutate(cells = tidyr::replace_na(cells, 0L)) %>%
    group_by(period) %>% mutate(percentage = 100 * cells / sum(cells)) %>% ungroup() %>%
    mutate(threshold = thresholds[i])
}))
write.csv(coverage_summary, "output/niche_annual_coverage_threshold_sensitivity.csv", row.names = FALSE)
supp_pct <- coverage_summary %>% filter(threshold == 0.05, category == "Supplementary only")
print(supp_pct[, c("period", "percentage")])

saveRDS(list(
    env_class_df = env_class_df,
    supp_pct = supp_pct,
    pc_var = pc_var
  ),"output/niche_annual_plotting_data.rds"
)

#####################################################
# MAPS AND BAR CHART
fig4_data <- readRDS("output/niche_annual_plotting_data.rds")

env_class_df <- fig4_data$env_class_df
supp_pct <- fig4_data$supp_pct
pc_var <- fig4_data$pc_var

# MAPS AND BAR CHART
p1 <- ggplot(env_class_df, aes(PC1, PC2, fill = category)) +
  geom_raster() + facet_wrap(~period, ncol = 5, drop = FALSE) +
  scale_fill_manual(values = c("GBIF only" = "darkgoldenrod1", "Supplementary only" = "darkorchid3",
                               "Both" = "skyblue"), drop = FALSE) +
  scale_x_continuous(breaks = c(-6, -3, 0)) + coord_equal(expand = FALSE) +
  labs(x = sprintf("Environmental PC1 (%.1f%%)", pc_var[1]),
       y = sprintf("Environmental PC2 (%.1f%%)", pc_var[2]), fill = NULL) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
        panel.spacing = grid::unit(0.3, "lines"), plot.margin = margin(5, 5, 5, 5))

# (b) yearly supplementary-only share
supp_pct$year <- as.integer(as.character(supp_pct$period))
p2 <- ggplot(supp_pct, aes(year, percentage)) +
  geom_line(colour = "darkorchid3", linewidth = 0.6) + geom_point(colour = "darkorchid3", size = 1.8) +
  scale_x_continuous(breaks = seq(2001, 2025, 4), limits = c(2001, 2025)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Supplementary-only space (%)") + theme_classic(base_size = 11) +
  theme(panel.grid.major.y = element_line(colour = "grey90"), plot.margin = margin(5, 5, 5, 5))

p_final <- (p1 / p2) +
  patchwork::plot_layout(heights = c(4, 1), guides = "collect") +
  patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = "bottom", plot.tag = element_text(face = "bold", size = 12))

print(p_final)

ggsave("output/figure/environmental_coverage_annual.png", p_final,
       width = 11.69, height = 13, dpi = 400, bg = "white")

#####################################################
# Five-year means of annual results
win_labels <- c(yr1 = "2001–2005", yr2 = "2006–2010", yr3 = "2011–2015", yr4 = "2016–2020", yr5 = "2021–2025")
win_years <- split(periods, rep(names(win_labels), each = 5))
mean_density <- function(grids, yrs) {
  z <- lapply(yrs, function(y) { r <- grids[[y]]$z.uncor; r / max(terra::values(r), na.rm = TRUE) })
  terra::app(terra::rast(z), mean, na.rm = TRUE)
}
classify_mean <- function(w, thr) {
  yrs <- intersect(win_years[[w]], ok_years)
  g <- terra::values(mean_density(grid_gbif, yrs), mat = FALSE)
  s <- terra::values(mean_density(grid_supp, yrs), mat = FALSE)
  gp <- is.finite(g) & g >= thr * max(g, na.rm = TRUE)
  sp <- is.finite(s) & s >= thr * max(s, na.rm = TRUE)
  cls <- rep(NA_character_, length(g))
  cls[gp & !sp] <- "GBIF only"; cls[!gp & sp] <- "Supplementary only"; cls[gp & sp] <- "Both"
  cells <- which(!is.na(cls)); xy <- terra::xyFromCell(grid_gbif[[yrs[1]]]$z.uncor, cells)
  data.frame(PC1 = xy[, 1], PC2 = xy[, 2], period = factor(win_labels[[w]], levels = win_labels),
             category = factor(cls[cells], levels = categories))
}
env_class_df <- dplyr::bind_rows(lapply(names(win_labels), classify_mean, thr = 0.05))
supp_pct <- coverage_summary %>% filter(category == "Supplementary only") %>%
  mutate(year = as.character(period),
         window = win_labels[vapply(year, function(y) names(win_years)[sapply(win_years, function(v) y %in% v)], "")]) %>%
  group_by(threshold, window) %>%
  summarise(percentage = mean(percentage), n_years = n(), .groups = "drop") %>%
  rename(period = window)
write.csv(supp_pct, "output/niche_annual5_supp_only_mean.csv", row.names = FALSE)
supp_pct <- supp_pct %>% filter(threshold == 0.05)
print(supp_pct)
period_levels <- unname(win_labels)

# Maps and bar chart
p1 <- ggplot(env_class_df, aes(PC1, PC2, fill = category)) +
  geom_raster() + facet_wrap(~period, nrow = 1) +
  scale_fill_manual(values = c("GBIF only" = "darkgoldenrod1", "Supplementary only" = "darkorchid3",
                               "Both" = "skyblue"), drop = FALSE) +
  scale_x_continuous(breaks = c(-6, -3, 0)) + coord_equal(expand = FALSE) +
  labs(x = sprintf("Environmental PC1 (%.1f%%)", pc_var[1]),
       y = sprintf("Environmental PC2 (%.1f%%)", pc_var[2]), fill = NULL) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
        panel.spacing = grid::unit(0.3, "lines"), plot.margin = margin(5, 5, 5, 5))

xmax <- max(50, ceiling(max(supp_pct$percentage, na.rm = TRUE) / 10) * 10)
p2 <- ggplot(supp_pct %>% mutate(strip = " "), aes(percentage, factor(period, levels = rev(period_levels)))) +
  geom_col(fill = "darkorchid3", width = 0.6) + facet_wrap(~strip, nrow = 1) +
  scale_x_continuous(limits = c(0, xmax), breaks = seq(0, xmax, by = 10),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Supplementary-only space (%)", y = NULL) + theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), axis.text.y = element_text(size = 9),
        panel.grid.major.x = element_line(colour = "grey90"), plot.margin = margin(5, 5, 5, 5))

p_final <- (p1 | p2) +
  patchwork::plot_layout(widths = c(5, 1.5), guides = "collect") &
  theme(legend.position = "bottom")
ggsave("output/figure/environmental_coverage_annual5.png", p_final,
       width = 11.69, height = 4.3, dpi = 400, bg = "white")

# GCB export (180 mm)
saveRDS(list(env_class_df = env_class_df, supp_pct = supp_pct, pc_var = pc_var), "output/niche_annual5_plotting_data.rds")
g4 <- ((p1 + theme(legend.position = "top", legend.justification = "right")) |
       (p2 + theme(legend.position = "none", plot.margin = margin(5, 14, 5, 5)))) +
  patchwork::plot_layout(widths = c(5, 1.6)) +
  patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold", size = 12))
dir.create("output/figure/gcb", showWarnings = FALSE)
ggsave("output/figure/gcb/Figure_4.tif", g4, width = 180, height = 60, units = "mm", dpi = 600,
       device = ragg::agg_tiff, compression = "lzw", scaling = 0.8, bg = "white")
