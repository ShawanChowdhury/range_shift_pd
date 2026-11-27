# Loading required libraries
library(ecospat)
library(terra)
library(ade4)

set.seed(123)

# Input and output folders
occ_dir <- Sys.getenv("OCC_DIR", "data/occ")
out_dir <- Sys.getenv("OUT_DIR", "output")
dir.create(file.path(out_dir, "niche_ckpt"), showWarnings = FALSE, recursive = TRUE)

# Niche grids, background density reused
bg_density <- function(glob, R = 100) {
  glob <- as.matrix(glob)
  xmin <- apply(glob, 2, min, na.rm = TRUE)
  xmax <- apply(glob, 2, max, na.rm = TRUE)
  ext <- c(xmin[1], xmax[1], xmin[2], xmax[2])
  dens <- ecospat:::ecospat.kd(x = glob, ext = ext, method = "adehabitat", th = 0, R = R)
  Z <- dens * nrow(glob) / terra::global(dens, "sum")$sum
  list(glob = glob, ext = ext, dens = dens, Z = Z, R = R)
}

grid_fast <- function(bg, sp) {
  sp <- as.matrix(sp)
  R <- bg$R
  ext <- bg$ext
  sp.dens <- ecospat:::ecospat.kd(x = sp, ext = ext, method = "adehabitat", th = 0, R = R,
                                  env.mask = bg$dens > 0)
  z <- sp.dens * nrow(sp) / terra::global(sp.dens, "sum")$sum
  z.uncor <- z / terra::global(z, "max")$max
  z.cor <- z / bg$Z
  z.cor[is.na(z.cor)] <- 0
  z.cor <- z.cor / terra::global(z.cor, "max")$max
  w <- z.uncor
  w[w > 0] <- 1
  l <- list()
  l$y <- seq(from = ext[3], to = ext[4], length.out = R)
  l$x <- seq(from = ext[1], to = ext[2], length.out = R)
  l$z <- z
  l$z.uncor <- z.uncor
  l$z.cor <- z.cor
  l$Z <- bg$Z
  l$glob <- bg$glob
  l$glob1 <- bg$glob
  l$sp <- sp
  l$w <- w
  l
}

# Equivalency test, reused background
equivalency_fast <- function(z1, z2, bg, rep) {
  obs.o <- c(ecospat.niche.overlap(z1, z2, cor = TRUE),
             ecospat.niche.dyn.index(z1, z2, intersection = 0)$dynamic.index.w)
  occ.pool <- rbind(z1$sp, z2$sp)
  row.names(occ.pool) <- c()
  sim <- t(sapply(1:rep, function(k) {
    rand.row <- sample(1:nrow(occ.pool), nrow(z1$sp))
    z1.sim <- grid_fast(bg, occ.pool[rand.row, , drop = FALSE])
    z2.sim <- grid_fast(bg, occ.pool[-rand.row, , drop = FALSE])
    o <- ecospat.niche.overlap(z1.sim, z2.sim, cor = TRUE)
    c(D = o$D, I = o$I)
  }))
  list(sim = as.data.frame(sim), obs = obs.o)
}

years <- c("yr1", "yr2", "yr3", "yr4", "yr5")

vars_final <- c("PDSI", "ppt", "soil", "tmin", "ws", "elev")

# Find occurrence files
occ_list <- list.files(
  "data/occ",
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

# Background = SDM calibration area
library(sf)
occs_all <- unique(read.csv("data/occ/com.csv")[, c("longitude", "latitude")])
env_template <- rast("data/env/env_overall.tif")
calibration_vect <- vect(
  st_as_sf(occs_all, coords = c("longitude", "latitude"), crs = crs(env_template)) |>
    st_transform("+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") |>
    st_buffer(dist = 1000000) |>
    st_union() |>
    st_sf() |>
    st_transform(crs(env_template))
)

# Empty result tables
niche_overlap <- data.frame()
pca_loadings <- data.frame()
pca_variance <- data.frame()

# Objects for figures
niche_objects <- list()

for (i in years) {
  print(i)

  # Resume from checkpoint
  ckpt <- file.path(out_dir, "niche_ckpt", paste0(i, ".rds"))
  if (file.exists(ckpt)) {
    ck <- readRDS(ckpt)
    niche_overlap <- rbind(niche_overlap, ck$niche_overlap)
    pca_loadings  <- rbind(pca_loadings, ck$loadings)
    pca_variance  <- rbind(pca_variance, ck$variance)
    niche_objects[[i]] <- ck$objects
    assign(".Random.seed", ck$seed, envir = globalenv())
    next
  }

  gbif <- read.csv(
    paste0(occ_dir, "/", i, "_gbif.csv")  # period files from 02
  )

  combined <- read.csv(
    paste0(occ_dir, "/", i, "_com.csv")
  )

  # Supplementary-only records (no GBIF)
  combined <- combined[!grepl("GBIF", combined$sources), ]

  gbif <- gbif[, c("longitude", "latitude", "year")]
  combined <- combined[, c("longitude", "latitude", "year")]

  # Environmental data
  pred <- rast(
    paste0(
      "data/env/env_",
      i,
      ".tif"
    )
  )

  pred <- pred[[vars_final]]
  pred <- crop(pred, calibration_vect, mask = TRUE)  # calibration-area background

  # Climate of record's year
  extract_own_year <- function(pts) {
    out <- list()
    for (y in sort(unique(pts$year))) {
      sub <- pts[pts$year == y, ]
      env_y <- rast(paste0("data/env_annual/env_", y, ".tif"))[[vars_final]]
      out[[as.character(y)]] <- extract(env_y, sub[, c("longitude", "latitude")], ID = FALSE)
    }
    do.call(rbind, out)
  }

  pred_gbif <- extract_own_year(gbif)
  pred_combined <- extract_own_year(combined)

  pred_gbif <- pred_gbif[
    complete.cases(pred_gbif),
  ]

  pred_combined <- pred_combined[
    complete.cases(pred_combined),
  ]

  # Environmental background
  globalEnv <- values(
    pred,
    na.rm = TRUE
  )

  # Add occurrences to background
  globalEnv <- rbind(globalEnv, as.matrix(pred_gbif[, colnames(globalEnv)]),
                     as.matrix(pred_combined[, colnames(globalEnv)]))

  # PCA-env
  pca.clim <- dudi.pca(
    globalEnv,
    center = TRUE,
    scale = TRUE,
    scannf = FALSE,
    nf = 2
  )

  # Background scores
  global.scores <- pca.clim$li

  # Occurrence scores
  gbif.scores <- suprow(
    pca.clim,
    pred_gbif[, colnames(globalEnv)]
  )$lisup

  combined.scores <- suprow(
    pca.clim,
    pred_combined[, colnames(globalEnv)]
  )$lisup

  # Niche grids
  bg <- bg_density(global.scores, R = 100)
  gbif.grid <- grid_fast(bg, gbif.scores)
  combined.grid <- grid_fast(bg, combined.scores)

  # Schoener's D
  D <- ecospat.niche.overlap(
    combined.grid,
    gbif.grid,
    cor = TRUE
  )$D

  # Uncorrected overlap
  D_uncor <- ecospat.niche.overlap(
    combined.grid,
    gbif.grid,
    cor = FALSE
  )$D

  # Niche equivalency test
  equiv.test <- equivalency_fast(
    gbif.grid,
    combined.grid,
    bg,
    rep = 100
  )

  obs.D <- equiv.test$obs["D"]
  sim.D <- equiv.test$sim[, "D"]

  # Proportion of simulations >= observed
  equiv.pval <- sum(
    sim.D >= obs.D
  ) / length(sim.D)

  # Niche similarity tests GBIF -> combined
  sim.test.gbif <- ecospat.niche.similarity.test(
    gbif.grid,
    combined.grid,
    rep = 100
  )

  obs.gbif.D <- sim.test.gbif$obs["D"]
  sim.gbif.D <- sim.test.gbif$sim[, "D"]

  sim.pval.gbif <- sum(
    sim.gbif.D >= obs.gbif.D
  ) / length(sim.gbif.D)

  # Combined -> GBIF
  sim.test.combined <- ecospat.niche.similarity.test(
    combined.grid,
    gbif.grid,
    rep = 100
  )

  obs.combined.D <- sim.test.combined$obs["D"]
  sim.combined.D <- sim.test.combined$sim[, "D"]

  sim.pval.combined <- sum(
    sim.combined.D >= obs.combined.D
  ) / length(sim.combined.D)

  # Save niche results
  niche_overlap <- rbind(
    niche_overlap,
    data.frame(
      year = i,
      n_gbif = nrow(pred_gbif),
      n_supplementary = nrow(pred_combined),
      schoener_D = D,
      schoener_D_uncorrected = D_uncor,
      equivalency_p = equiv.pval,
      similarity_gbif_to_supplementary_p = sim.pval.gbif,
      similarity_supplementary_to_gbif_p = sim.pval.combined
    )
  )

  # PCA loadings
  loadings <- data.frame(
    year = i,
    variable = rownames(pca.clim$co),
    PC1 = pca.clim$co[, 1],
    PC2 = pca.clim$co[, 2]
  )

  pca_loadings <- rbind(
    pca_loadings,
    loadings
  )

  # Variance explained by PCA
  pc1 <- (
    pca.clim$eig[1] /
      sum(pca.clim$eig)
  ) * 100

  pc2 <- (
    pca.clim$eig[2] /
      sum(pca.clim$eig)
  ) * 100

  pca_variance <- rbind(
    pca_variance,
    data.frame(
      year = i,
      PC1_percent = pc1,
      PC2_percent = pc2,
      PC1_PC2_percent = pc1 + pc2
    )
  )

  # Save objects for figure script
  niche_objects[[i]] <- list(
    pca = pca.clim,
    global_scores = global.scores,
    gbif_scores = gbif.scores,
    combined_scores = combined.scores,
    gbif_grid = gbif.grid,
    combined_grid = combined.grid,
    equivalency = equiv.test,
    similarity_gbif = sim.test.gbif,
    similarity_combined = sim.test.combined
  )

  saveRDS(list(niche_overlap = niche_overlap[niche_overlap$year == i, ],
               loadings = loadings,
               variance = pca_variance[pca_variance$year == i, ],
               objects = niche_objects[[i]],
               seed = .Random.seed), paste0(ckpt, ".tmp"))
  file.rename(paste0(ckpt, ".tmp"), ckpt)  # keep complete checkpoints only
}

# Export results
write.csv(
  niche_overlap,
  file.path(out_dir, "niche_overlap.csv"),
  row.names = FALSE
)

write.csv(
  pca_loadings,
  file.path(out_dir, "pca_loadings.csv"),
  row.names = FALSE
)

write.csv(
  pca_variance,
  file.path(out_dir, "pca_variance.csv"),
  row.names = FALSE
)

saveRDS(
  niche_objects,
  file.path(out_dir, "niche_objects.rds")
)

# Remove checkpoints
unlink(file.path(out_dir, "niche_ckpt"), recursive = TRUE)
