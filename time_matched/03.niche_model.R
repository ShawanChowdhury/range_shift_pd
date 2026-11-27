# Loading required libraries
library(terra)
library(sf)
library(dplyr)
library(readr)
library(ENMeval)

#########################################
# Settings
periods <- paste0("yr", 1:5)

# Which parts to run
run_per_period <- as.logical(Sys.getenv("RUN_PER_PERIOD", "FALSE"))
run_annual     <- as.logical(Sys.getenv("RUN_ANNUAL_MODEL", "TRUE"))

dir.create("data/bg", showWarnings = FALSE)
dir.create("output/sdm/gbif", showWarnings = FALSE, recursive = TRUE)
dir.create("output/sdm/combined", showWarnings = FALSE, recursive = TRUE)

#########################################
# Common calibration area
com <- read_csv("data/occ/com.csv", show_col_types = FALSE)

# Spatial template
env_overall <- rast("data/env/env_overall.tif")

# Occurrence coordinates
occs_all <- com %>%
  dplyr::select(longitude, latitude) %>%
  distinct()

# Convert occurrences to sf
occs_sf <- st_as_sf(occs_all,
                    coords = c("longitude", "latitude"),
                    crs = crs(env_overall))

# Equal-area projection for buffering
eckertIV <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# 1,000-km buffer around occurrences
calibration_area <- occs_sf %>%
  st_transform(eckertIV) %>%
  st_buffer(dist = 1000000) %>%
  st_union() %>%
  st_sf() %>%
  st_transform(crs(env_overall))

calibration_vect <- vect(calibration_area)

#########################################
# Importing background points
bg <- read_csv("data/bg.csv", show_col_types = FALSE)

#########################################
# Function for running SDMs
run_sdm <- function(data_file, dataset_name) {

  # Import occurrence data
  occ_data <- read_csv(data_file, show_col_types = FALSE)

  for (i in periods) {

    print(paste(dataset_name, i))

    # Predictors for this period
    envs <- rast(paste0("data/env/env_", i, ".tif"))

    # Crop to calibration area
    envs_bg <- crop(envs, calibration_vect, mask = TRUE)

    # Filter records by period
    occs <- occ_data %>%
      filter(year_group == i) %>%
      dplyr::select(longitude, latitude)

    # Remove incomplete occurrences
    occ_values <- terra::extract(
      envs_bg,
      as.matrix(occs[, c("longitude", "latitude")])
    )

    occs <- occs[complete.cases(occ_values), ]

    # Complete background points
    bg_values <- terra::extract(
      envs_bg,
      as.matrix(bg[, c("longitude", "latitude")])
    )

    bg_i <- bg[complete.cases(bg_values), ]

    # Run Maxent models
    model <- ENMevaluate(
      occs = occs,
      envs = envs_bg,
      bg = bg_i,
      partitions = "checkerboard",
      partition.settings = list(
        aggregation.factor = c(2, 2)
      ),
      tune.args = list(
        fc = c("L", "LQ", "LQP", "LQPT"),  # hinge fails in CV
        rm = seq(0.5, 4, 0.5)
      ),
      algorithm = "maxnet",
      raster.preds = FALSE,
      rmm = list(),
      parallel = FALSE
    )

    # Model evaluation results
    res <- eval.results(model)

    # Select model with lowest AICc
    opt_aicc <- res %>%
      filter(!is.na(AICc)) %>%
      slice_min(AICc, n = 1, with_ties = FALSE)

    # Extract fitted models
    models <- eval.models(model)

    # Identify selected model
    best_name <- as.character(opt_aicc$tune.args)
    best_model <- models[[best_name]]

    # Selected model suitability raster
    r <- ENMeval::maxnet.predictRaster(
      best_model,
      envs_bg
    )

    # Suitability at occurrences
    occ_predictions <- terra::extract(
      r,
      as.matrix(occs[, c("longitude", "latitude")])
    )

    # Calculate 10th-percentile occurrence threshold
    threshold <- quantile(
      occ_predictions[, 1],
      probs = 0.10,
      na.rm = TRUE
    )

    # Add threshold to model summary
    opt_aicc$threshold_10p <- threshold

    # Export model summary
    write_csv(
      opt_aicc,
      paste0(
        "output/sdm/", dataset_name, "/",
        dataset_name, "_model_sum_", i, ".csv"
      )
    )

    # Save complete ENMeval object
    saveRDS(
      model,
      paste0(
        "output/sdm/", dataset_name, "/",
        dataset_name, "_model_", i, ".rds"
      )
    )

    # Export suitability raster
    writeRaster(
      r,
      paste0(
        "output/sdm/", dataset_name, "/",
        dataset_name, "_map_", i, ".tif"
      ),
      NAflag = -9999,
      overwrite = TRUE
    )

    # Create binary suitability raster
    r_bin <- ifel(
      r >= threshold,
      1,
      0
    )

    # Export binary raster
    writeRaster(
      r_bin,
      paste0(
        "output/sdm/", dataset_name, "/",
        dataset_name, "_bin_", i, ".tif"
      ),
      NAflag = -9999,
      overwrite = TRUE
    )
  }
}

#########################################
# Run GBIF-only models
if (run_per_period) run_sdm(data_file = "data/occ/gbif.csv", dataset_name = "gbif")

#########################################
# Run combined-data models
if (run_per_period) run_sdm(data_file = "data/occ/com.csv", dataset_name = "combined")

#########################################
# Held-out evaluation and fold sizes
eval_table <- list()
if (run_per_period) for (ds in c("gbif", "combined")) for (i in periods) {
  m   <- readRDS(paste0("output/sdm/", ds, "/", ds, "_model_", i, ".rds"))
  sel <- read_csv(paste0("output/sdm/", ds, "/", ds, "_model_sum_", i, ".csv"),
                  show_col_types = FALSE)
  eval_table[[length(eval_table) + 1]] <- tibble(
    dataset = ds, period = i,
    n_occ = nrow(m@occs), n_bg = nrow(m@bg),
    occ_per_fold = paste(as.integer(table(m@occs.grp)), collapse = "/"),
    bg_per_fold  = paste(as.integer(table(m@bg.grp)), collapse = "/"),
    selected = sel$tune.args, ncoef = sel$ncoef,
    AUC_train = sel$auc.train, AUC_test = sel$auc.val.avg, AUC_test_sd = sel$auc.val.sd,
    AUC_diff = sel$auc.diff.avg, CBI_test = sel$cbi.val.avg, CBI_test_sd = sel$cbi.val.sd,
    OR10_test = sel$or.10p.avg, OR10_test_sd = sel$or.10p.sd, ORmtp_test = sel$or.mtp.avg,
    threshold_10p = sel$threshold_10p)
}
if (run_per_period) write_csv(bind_rows(eval_table), "output/sdm/model_evaluation_table.csv")

#########################################
# Time-matched annual model
years_all <- 2001:2025

run_annual_sdm <- function(dataset_name = "annual", bg_years_per_point = 4, keep_sources = NULL) {

  out_dir <- paste0("output/sdm/", dataset_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  envs_list <- lapply(years_all, function(y) {
    crop(rast(paste0("data/env_annual/env_", y, ".tif")), calibration_vect, mask = TRUE)
  })
  names(envs_list) <- years_all

  # One record per cell-year
  occ_all <- read_csv("data/com_data.csv", show_col_types = FALSE) %>%
    filter(between(year, min(years_all), max(years_all)))
  # Optional source subset (GBIF-only model)
  if (!is.null(keep_sources)) occ_all <- filter(occ_all, source %in% keep_sources)
  occ_all$cell <- cellFromXY(envs_list[[1]][[1]], as.matrix(occ_all[, c("longitude", "latitude")]))
  occ_all <- occ_all %>%
    filter(!is.na(cell)) %>%
    distinct(year, cell, .keep_all = TRUE)

  occs_swd <- bind_rows(lapply(years_all, function(y) {
    o <- occ_all %>% filter(year == y)
    if (nrow(o) == 0) return(NULL)
    v <- terra::extract(envs_list[[as.character(y)]], as.matrix(o[, c("longitude", "latitude")]))
    bind_cols(o[, c("longitude", "latitude")], v, year = y, source = o$source)
  })) %>% filter(complete.cases(.))

  # Background with random years
  set.seed(2026)
  bg_t <- bg[rep(seq_len(nrow(bg)), each = bg_years_per_point), c("longitude", "latitude")]
  bg_t$year <- sample(years_all, nrow(bg_t), replace = TRUE)
  bg_swd <- bind_rows(lapply(years_all, function(y) {
    b <- bg_t %>% filter(year == y)
    v <- terra::extract(envs_list[[as.character(y)]], as.matrix(b[, c("longitude", "latitude")]))
    bind_cols(b[, c("longitude", "latitude")], v, year = y)
  })) %>% filter(complete.cases(.)) %>% distinct()

  write_csv(count(occs_swd, year, name = "n_cells"), paste0(out_dir, "/", dataset_name, "_occ_per_year.csv"))
  print(table(occs_swd$year))

  # Checkerboard folds by location
  occ_xy <- distinct(occs_swd[, c("longitude", "latitude")])
  bg_xy  <- distinct(bg_swd[, c("longitude", "latitude")])
  grp_xy <- get.checkerboard(occs = occ_xy, envs = envs_list[[1]], bg = bg_xy,
                             aggregation.factor = c(2, 2))
  occ_xy$grp <- grp_xy$occs.grp
  bg_xy$grp  <- grp_xy$bg.grp
  grp <- list(
    occs.grp = left_join(occs_swd[, c("longitude", "latitude")], occ_xy, by = c("longitude", "latitude"))$grp,
    bg.grp   = left_join(bg_swd[, c("longitude", "latitude")], bg_xy, by = c("longitude", "latitude"))$grp
  )
  stopifnot(length(grp$occs.grp) == nrow(occs_swd), !anyNA(grp$occs.grp),
            length(grp$bg.grp) == nrow(bg_swd), !anyNA(grp$bg.grp))

  env_names <- names(envs_list[[1]])
  model <- ENMevaluate(
    occs = as.data.frame(occs_swd[, c("longitude", "latitude", env_names)]),
    bg = as.data.frame(bg_swd[, c("longitude", "latitude", env_names)]),
    partitions = "user",
    user.grp = grp,
    tune.args = list(fc = c("L", "LQ", "LQP", "LQPT"), rm = seq(0.5, 4, 0.5)),
    algorithm = "maxnet",
    parallel = TRUE,
    numCores = 8
  )

  res <- eval.results(model)
  opt_aicc <- res %>% filter(!is.na(AICc)) %>% slice_min(AICc, n = 1, with_ties = FALSE)
  best_model <- eval.models(model)[[as.character(opt_aicc$tune.args)]]

  # Common 10th-percentile threshold
  occ_pred <- predict(best_model, as.data.frame(occs_swd[, env_names]), type = "cloglog", clamp = TRUE)
  threshold <- quantile(occ_pred, probs = 0.10, na.rm = TRUE)

  # Held-out TSS
  bg_pred <- predict(best_model, as.data.frame(bg_swd[, env_names]), type = "cloglog", clamp = TRUE)
  tss <- sapply(sort(unique(grp$occs.grp)), function(k) {
    m_k <- maxnet::maxnet(
      p = c(rep(1, sum(grp$occs.grp != k)), rep(0, sum(grp$bg.grp != k))),
      data = rbind(as.data.frame(occs_swd[grp$occs.grp != k, env_names]),
                   as.data.frame(bg_swd[grp$bg.grp != k, env_names])),
      f = maxnet::maxnet.formula(
        c(rep(1, sum(grp$occs.grp != k)), rep(0, sum(grp$bg.grp != k))),
        rbind(as.data.frame(occs_swd[grp$occs.grp != k, env_names]),
              as.data.frame(bg_swd[grp$bg.grp != k, env_names])),
        classes = tolower(as.character(opt_aicc$fc))),
      regmult = as.numeric(as.character(opt_aicc$rm)))
    tr <- predict(m_k, as.data.frame(occs_swd[grp$occs.grp != k, env_names]), type = "cloglog", clamp = TRUE)
    thr_k <- quantile(tr, 0.10, na.rm = TRUE)
    po <- predict(m_k, as.data.frame(occs_swd[grp$occs.grp == k, env_names]), type = "cloglog", clamp = TRUE)
    pb <- predict(m_k, as.data.frame(bg_swd[grp$bg.grp == k, env_names]), type = "cloglog", clamp = TRUE)
    mean(po >= thr_k) + mean(pb < thr_k) - 1
  })

  opt_aicc$threshold_10p <- threshold
  opt_aicc$TSS_test_avg  <- mean(tss)
  opt_aicc$TSS_test_sd   <- sd(tss)
  opt_aicc$n_occ <- nrow(occs_swd)
  opt_aicc$n_bg  <- nrow(bg_swd)
  opt_aicc$occ_per_fold <- paste(as.integer(table(grp$occs.grp)), collapse = "/")
  write_csv(opt_aicc, paste0(out_dir, "/", dataset_name, "_model_sum.csv"))
  write_csv(res, paste0(out_dir, "/", dataset_name, "_tuning_results.csv"))
  saveRDS(model, paste0(out_dir, "/", dataset_name, "_model.rds"))

  # Annual projections
  for (y in years_all) {
    r <- ENMeval::maxnet.predictRaster(best_model, envs_list[[as.character(y)]])
    writeRaster(r, paste0(out_dir, "/", dataset_name, "_map_", y, ".tif"),
                NAflag = -9999, overwrite = TRUE)
    writeRaster(ifel(r >= threshold, 1, 0), paste0(out_dir, "/", dataset_name, "_bin_", y, ".tif"),
                NAflag = -9999, overwrite = TRUE)
  }
}

# Annual model data sets
annual_sets <- strsplit(Sys.getenv("ANNUAL_DATASETS", "annual"), ",")[[1]]
if (run_annual) {
  if ("annual" %in% annual_sets) run_annual_sdm()
  if ("annual_gbif" %in% annual_sets) run_annual_sdm("annual_gbif", keep_sources = "GBIF")
}
