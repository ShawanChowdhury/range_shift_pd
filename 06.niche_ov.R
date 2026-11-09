# Loading required libraries
library(ecospat)
library(terra)
library(ade4)

set.seed(123)

years <- c("yr1", "yr2", "yr3", "yr4", "yr5")

vars_final <- c("PDSI", "ppt", "soil", "tmin", "ws", "elev")

# Find occurrence files
occ_list <- list.files(
  "data/occ",
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

# Empty result tables
niche_overlap <- data.frame()
pca_loadings <- data.frame()
pca_variance <- data.frame()

# Save objects needed later for figures
niche_objects <- list()

for (i in years) {
  print(i)
  
  gbif <- read.csv(
    paste0("data/occ/yr/", i, "_gbif.csv")
  )
  
  combined <- read.csv(
    paste0("data/occ/yr/", i, "_com.csv")
  )
  
  gbif <- gbif[, c("longitude", "latitude")]
  combined <- combined[, c("longitude", "latitude")]
  
  # Environmental data
  pred <- rast(
    paste0(
      "data/env/env_",
      i,
      ".tif"
    )
  )
  
  pred <- pred[[vars_final]]
  
  # Extract environment at occurrences
  pred_gbif <- extract(
    pred,
    gbif,
    ID = FALSE
  )
  
  pred_combined <- extract(
    pred,
    combined,
    ID = FALSE
  )
  
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
  gbif.grid <- ecospat.grid.clim.dyn(
    global.scores,
    global.scores,
    gbif.scores,
    R = 100
  )
  
  combined.grid <- ecospat.grid.clim.dyn(
    global.scores,
    global.scores,
    combined.scores,
    R = 100
  )
  
  # Schoener's D
  D <- ecospat.niche.overlap(
    combined.grid,
    gbif.grid,
    cor = TRUE
  )$D
  
  # Niche equivalency test
  equiv.test <- ecospat.niche.equivalency.test(
    gbif.grid,
    combined.grid,
    rep = 100
  )
  
  obs.D <- equiv.test$obs["D"]
  sim.D <- equiv.test$sim[, "D"]
  
  # Same calculation used in the previous analysis/paper
  equiv.pval <- sum(
    sim.D >= obs.D
  ) / length(sim.D)
  
  # Niche similarity tests
  # GBIF -> combined
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
      n_combined = nrow(pred_combined),
      schoener_D = D,
      equivalency_p = equiv.pval,
      similarity_gbif_to_combined_p = sim.pval.gbif,
      similarity_combined_to_gbif_p = sim.pval.combined
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
}

# Export results
write.csv(
  niche_overlap,
  "output/niche_overlap.csv",
  row.names = FALSE
)

write.csv(
  pca_loadings,
  "output/pca_loadings.csv",
  row.names = FALSE
)

write.csv(
  pca_variance,
  "output/pca_variance.csv",
  row.names = FALSE
)

saveRDS(
  niche_objects,
  "output/niche_objects.rds"
)
