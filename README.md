# The complex range dynamics of a globally expanding butterfly

R code for the analyses of *Papilio demoleus* range dynamics, 2001–2025.

## Structure

| Folder | Contents |
|---|---|
| `scripts/` | Main workflow: per-period and unified models |
| `time_matched/` | Time-matched annual model |

Main workflow (`scripts/`):

| Script | Step |
|---|---|
| `01.spatial_data.R` | Environmental predictors |
| `02.occ_data.R` | Occurrence data (GBIF, Facebook, Flickr, literature, Lohman et al. 2025) |
| `03.niche_model.R` | MaxEnt models (per-period; unified with `RUN_UNIFIED=TRUE`) |
| `04.full_extent.R` | Binary maps to full extent |
| `05.raster_analysis.R` | Range dynamics |
| `06.niche_ov.R` | Niche overlap (ecospat) |
| `07.figures_up.R` | Figures 1–3 |
| `08_temporal_niche_figure.R` | Figure 4 |

The `time_matched/` scripts follow the same steps for the time-matched annual model.

## Running

Run the scripts in order from the project root, with input data in `data/` and outputs written to `output/`.
Options are set with environment variables:

| Script | Variable | Default |
|---|---|---|
| `03` | `RUN_PER_PERIOD`, `RUN_UNIFIED` | `TRUE`, `FALSE` |
| `05` | `SDM_SOURCES` | `gbif,combined,unified` |
| `06` | `OCC_DIR`, `OUT_DIR` | `data/occ`, `output` |
| `07` | `FIG2_SOURCES`, `RUN_FIG13` | `combined,unified`, `TRUE` |

## Data

- GBIF download: https://doi.org/10.15468/dl.7hvqew
- Tropical Asian butterfly occurrences (Lohman et al. 2025): https://doi.org/10.1038/s41597-025-05333-w
- Climate data: TerraClimate
