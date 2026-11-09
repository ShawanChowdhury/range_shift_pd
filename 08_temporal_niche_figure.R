# Load libraries
library(terra)
library(ade4)
library(ecospat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

set.seed(123)
periods <- paste0("yr", 1:5)
period_labels <- setNames(c("2001–2005", "2006–2010", "2011–2015", "2016–2020", "2021–2025"), periods)
period_levels <- unname(period_labels)
vars <- c("PDSI", "ppt", "soil", "tmin", "ws", "elev")
categories <- c("GBIF only", "Supplementary only", "Both")
thresholds <- c(0.02, 0.05, 0.10)
REBUILD_GRIDS <- FALSE  # Set TRUE if environmental data, PCA or grid parameters have changed.

###############################################
# ENVIRONMENTAL BACKGROUND AND SHARED PCA
# Use the existing background coordinates in all five intervals.
env_list <- setNames(lapply(periods, function(p) {
  f <- paste0("data/env/env_", p, ".tif")
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
          "output/niche_shared_PCA_loadings.csv", row.names = FALSE)
write.csv(data.frame(axis = c("PC1", "PC2"), variance_percent = as.numeric(pc_var)),
          "output/niche_shared_PCA_variance.csv", row.names = FALSE)

###############################################
# OCCURRENCE SCORES AND CACHED NICHE GRIDS 
# com = GBIF + Facebook + literature; supp = Facebook + literature only.
get_scores <- function(period, source) {
  f <- paste0("data/occ/", period, "_", if (source == "supp") "com" else source, ".csv")
  if (!file.exists(f)) stop("Missing occurrence file: ", f)
  d <- read.csv(f)
  if (source == "supp") d <- d[tolower(trimws(d$source)) %in% c("facebook", "literature"), , drop = FALSE]
  if (nrow(d) < 5) stop("Too few records for ", period, " / ", source)
  z <- terra::extract(env_list[[period]], as.matrix(d[, c("longitude", "latitude")]))
  z <- z[complete.cases(z[, vars]), vars, drop = FALSE]
  if (nrow(z) < 5) stop("Too few valid environmental records for ", period, " / ", source)
  as.data.frame(ade4::suprow(pca_shared, z)$lisup[, 1:2, drop = FALSE])
}

occ_scores <- setNames(lapply(periods, function(p) {
  setNames(lapply(c("gbif", "com", "supp"), function(s) get_scores(p, s)), c("gbif", "com", "supp"))
}), periods)

record_counts <- dplyr::bind_rows(lapply(periods, function(p) {
  data.frame(period = period_labels[[p]], source = c("GBIF", "Combined", "Supplementary"),
             n = vapply(occ_scores[[p]], nrow, integer(1)))
}))
write.csv(record_counts, "output/niche_record_counts.csv", row.names = FALSE)

get_grid <- function(period, source) {
  sp <- occ_scores[[period]][[source]]
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

grid_gbif <- read_or_make_grids("gbif", "output/niche_grids_gbif.rds")
grid_com  <- read_or_make_grids("com",  "output/niche_grids_combined.rds")
grid_supp <- read_or_make_grids("supp", "output/niche_grids_supplementary.rds")

###############################################
# Rebuild the niche grids in the current R session
grid_com  <- setNames(lapply(periods, get_grid, source = "com"), periods)
grid_gbif <- setNames(lapply(periods, get_grid, source = "gbif"), periods)
grid_supp <- setNames(lapply(periods, get_grid, source = "supp"), periods)

# THREE-CATEGORY ENVIRONMENTAL COVERAGE AND THRESHOLD SENSITIVITY 
# The percentages describe thresholded PCA grid-cell coverage, not geographic area or niche breadth.
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

env_class_df <- dplyr::bind_rows(lapply(periods, function(p) {
  classify_env_space(grid_gbif[[p]], grid_supp[[p]], p, thr = 0.05)
}))

coverage_by_threshold <- lapply(thresholds, function(th) {
  dplyr::bind_rows(lapply(periods, function(p) classify_env_space(grid_gbif[[p]], grid_supp[[p]], p, th)))
})
names(coverage_by_threshold) <- as.character(thresholds)
env_class_df <- coverage_by_threshold[[as.character(0.05)]]

coverage_summary <- dplyr::bind_rows(lapply(seq_along(thresholds), function(i) {
  tab <- coverage_by_threshold[[i]] %>% count(period, category, name = "cells")
  tidyr::expand_grid(period = factor(period_levels, levels = period_levels),
                     category = factor(categories, levels = categories)) %>%
    left_join(tab, by = c("period", "category")) %>%
    mutate(cells = tidyr::replace_na(cells, 0L)) %>%
    group_by(period) %>% mutate(percentage = 100 * cells / sum(cells)) %>% ungroup() %>%
    mutate(threshold = thresholds[i])
}))
write.csv(coverage_summary, "output/niche_coverage_threshold_sensitivity.csv", row.names = FALSE)
supp_pct <- coverage_summary %>% filter(threshold == 0.05, category == "Supplementary only")
print(supp_pct[, c("period", "percentage")])

saveRDS(list(
    env_class_df = env_class_df,
    supp_pct = supp_pct,
    pc_var = pc_var
  ),"output/niche_space_plotting_data.rds"
)

#####################################################
# EQUAL-AXIS MAPS + VERTICALLY ALIGNED BAR CHART
# Starting from plotting
# Load saved plotting data
fig4_data <- readRDS("output/niche_space_plotting_data.rds")

env_class_df <- fig4_data$env_class_df
supp_pct <- fig4_data$supp_pct
pc_var <- fig4_data$pc_var

period_levels <- c("2001–2005", "2006–2010", "2011–2015",
                   "2016–2020", "2021–2025")

# EQUAL-AXIS MAPS + VERTICALLY ALIGNED BAR CHART
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
        panel.grid.major.x = element_line(colour = "grey90"), plot.margin = margin(0, 5, 0, 5))

# Match margins between the two panels
p1 <- p1 + theme(plot.margin = margin(5, 5, 5, 5))
p2 <- p2 + theme(plot.margin = margin(5, 5, 5, 5))

# Combine panels without vertical spacers
p_final <- (p1 | p2) +
  patchwork::plot_layout(widths = c(5, 1.5), guides = "collect") &
  theme(legend.position = "bottom")

print(p_final)

ggsave("output/figure/environmental_coverage.png", p_final,
       width = 11.69, height = 4.3, dpi = 400, bg = "white")
