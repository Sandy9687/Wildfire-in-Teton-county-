# =============================================================================
# Wildfire–Climate Spatial Analysis — Teton County, Wyoming
# GIST 5220 | Spring 2026
# -----------------------------------------------------------------------------
# Data sources:
#   VIIRS  : NASA FIRMS shapefiles (SUOMI, J1, J2)
#   PRISM  : Gridded .bil rasters via prism R package
#   NDVI   : GEE export GeoTIFF (Landsat 7/8/9 composite 2012–2025)
#   SPEI   : GEE export GeoTIFFs (1, 3, 6, 12-month scales)
#   Daymet : Downloaded via daymetr package
#   Boundary: tigris R package
# -----------------------------------------------------------------------------
# Units : SI (metric) throughout
# CRS   : EPSG:32612 (UTM Zone 12N — appropriate for Teton County)
# Fire season: May–October (winter excluded — snow cover makes fire impossible)
# =============================================================================

# --- Install packages (run once) ---------------------------------------------
# install.packages(c("sf", "terra", "dplyr", "readr", "purrr", "lubridate",
#                    "tidyr", "spdep", "GWmodel", "tmap", "ggplot2",
#                    "tigris", "daymetr", "prism"))

library(sf)        # vector spatial operations
library(terra)     # raster operations
library(dplyr)     # data manipulation
library(readr)     # CSV reading
library(purrr)     # functional iteration
library(lubridate) # date parsing
library(tidyr)     # pivot_wider
library(spdep)     # spatial weights, Moran's I, Getis-Ord Gi*
library(GWmodel)   # Geographically Weighted Regression / Logistic
library(tmap)      # map output
library(ggplot2)   # time series plots
library(tigris)    # Census boundary
library(daymetr)   # Daymet solar radiation
library(prism)     # PRISM climate data download

# --- Global Settings ---------------------------------------------------------
TARGET_CRS         <- "EPSG:32612"      # UTM Zone 12N — metres
TARGET_RES         <- 800               # metres — matches PRISM resolution
FIRE_SEASON_MONTHS <- c(5,6,7,8,9,10)  # May–October fire season
PEAK_MONTHS        <- c(6,7,8,9)        # June–September — climate aggregation
GI_DIST_M          <- 20000             # 20 km Gi* neighbourhood
KDE_RADIUS         <- 20000             # 20 km KDE search radius


# --- Output Directories ------------------------------------------------------
dir.create("output/maps",     recursive = TRUE, showWarnings = FALSE)
dir.create("output/rasters",  recursive = TRUE, showWarnings = FALSE)
dir.create("output/vectors",  recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",   recursive = TRUE, showWarnings = FALSE)
dir.create("output/plots",    recursive = TRUE, showWarnings = FALSE)
dir.create("intermediate",    showWarnings = FALSE)
dir.create("data/Daymet",     recursive = TRUE, showWarnings = FALSE)
#dir.create("data/PRISM_clean",recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# PRISM DOWNLOAD — Run once, then comment out after data is on disk
# =============================================================================
# Downloads monthly temperature and precipitation rasters (800 m) for
# 2012–2025. Files saved to data/PRISM_clean/ as .bil rasters.
# Uncomment the block below on first run only.

# prism_set_dl_dir("data/PRISM_clean")
#
# get_prism_monthlys(type = "tmean", years = 2012:2025, keepZip = FALSE)
# get_prism_monthlys(type = "ppt",   years = 2012:2025, keepZip = FALSE)
#
# cat("PRISM download complete.\n")


# =============================================================================
# STAGE 1 — DATA INGESTION
# =============================================================================
cat("\n===== STAGE 1: DATA INGESTION =====\n")

# -----------------------------------------------------------------------------
# 1.1  Study Boundary — Teton County, Wyoming
#      Source : US Census TIGER via tigris package
#      CRS    : EPSG:4326 on download; reprojected to EPSG:32612 in Stage 2
# -----------------------------------------------------------------------------
cat("  Loading Teton County boundary ...\n")
options(tigris_use_cache = TRUE)

study_area <- counties(state = "WY", cb = TRUE, year = 2022) %>%
  filter(NAME == "Teton") %>%
  st_transform("EPSG:4326")

cat(sprintf("  Teton County loaded. Area: %.1f km²\n",
            as.numeric(st_area(study_area)) / 1e6))


# Teton County bounding box for Daymet [top_lat, left_lon, bot_lat, right_lon]
# After loading study_area (Teton County polygon)
teton_bbox <- st_bbox(study_area)

# Convert to Daymet format: c(lat_max, lon_min, lat_min, lon_max)
DAYMET_BBOX <- c(teton_bbox$ymax,
                 teton_bbox$xmin,
                 teton_bbox$ymin,
                 teton_bbox$xmax)


# -----------------------------------------------------------------------------
# 1.2  VIIRS Active Fire Data — NASA FIRMS Shapefiles
#      Source   : firms.modaps.eosdis.nasa.gov
#      Files    : Three shapefiles — one per platform:
#                   SUOMI VIIRS C2 → S-NPP
#                   J1 VIIRS C2   → NOAA-20
#                   J2 VIIRS C2   → NOAA-21
#      Place files in: data/VIIRS/FIRMS/
#      FIRMS shapefile column names (lowercase after standardisation):
#                   latitude, longitude → geometry
#                   acq_date  : detection date (YYYY-MM-DD)
#                   acq_time  : UTC time (HHMM)
#                   confidence: "l" / "n" / "h" (lowercase in shapefiles)
#                   frp       : fire radiative power (MW)
#                   type      : 0=vegetation, 1=volcano, 2=static, 3=offshore
#                   daynight  : "D" or "N"
#      CRS      : EPSG:4326 (geographic WGS84)
#      Units    : FRP in megawatts (MW)
# -----------------------------------------------------------------------------
load_firms_shp <- function(filepath, platform_label) {
  cat(sprintf("  Loading %s ...\n", platform_label))

  if (!file.exists(filepath))
    stop(sprintf("Shapefile not found: %s", filepath))

  shp <- st_read(filepath, quiet = TRUE)

  # Standardise column names to lowercase
  names(shp) <- tolower(names(shp))

  # Verify essential columns
  required <- c("acq_date", "acq_time", "confidence", "frp", "type")
  missing  <- setdiff(required, names(shp))
  if (length(missing) > 0)
    stop(sprintf("Missing columns in %s: %s",
                 platform_label, paste(missing, collapse = ", ")))

  # Standardise fields
  shp <- shp %>%
    mutate(
      platform   = platform_label,
      acq_date   = ymd(acq_date),           # parse YYYY-MM-DD to Date
      year       = year(acq_date),
      month      = month(acq_date),
      # FIRMS shapefiles use lowercase confidence: "l"/"n"/"h" → uppercase
      confidence = toupper(as.character(confidence)),
      frp        = as.numeric(frp),
      type       = as.numeric(type)
    )

  # Ensure geometry is in EPSG:4326
  if (is.na(st_crs(shp))) {
    shp <- st_set_crs(shp, "EPSG:4326")
    cat(sprintf("    CRS was missing — set to EPSG:4326 for %s\n",
                platform_label))
  } else if (st_crs(shp)$epsg != 4326) {
    shp <- st_transform(shp, "EPSG:4326")
  }

  cat(sprintf("    %s: %d raw detections loaded.\n",
              platform_label, nrow(shp)))
  shp
}

# Load all three platform shapefiles
# Adjust filenames below to match your exact downloaded file names
# NEW — auto-finds the .shp file inside each platform subfolder
# No renaming needed — code searches the folder for whatever .shp is there

find_shp <- function(folder) {
  f <- list.files(folder, pattern = "\\.shp$", full.names = TRUE)
  if (length(f) == 0) stop(paste("No .shp file found in:", folder))
  if (length(f) > 1)  warning(paste("Multiple .shp files found — using first:", f[1]))
  f[1]
}

snpp_raw   <- load_firms_shp(
  find_shp("data/VIIRS/C2/SUOMI_VIIRS_C2_USA"), "SNPP")
noaa20_raw <- load_firms_shp(
  find_shp("data/VIIRS/C2/J1_VIIRS_C2_USA"),    "NOAA20")
noaa21_raw <- load_firms_shp(
  find_shp("data/VIIRS/C2/J2_VIIRS_C2_USA"),    "NOAA21")

# Merge all platforms into one layer
fire_raw <- bind_rows(snpp_raw, noaa20_raw, noaa21_raw)
cat(sprintf("\n  Total raw detections (all platforms): %d\n", nrow(fire_raw)))

# Verify FRP is numeric
stopifnot(is.numeric(fire_raw$frp))
cat("  Field verification passed.\n")

# -----------------------------------------------------------------------------
# 1.3  PRISM Climate Rasters
#      Source     : prism.oregonstate.edu via prism R package
#      Variables  : tmean (°C), ppt (mm) — monthly, 800 m resolution
#      CRS        : EPSG:4269 (NAD83 geographic) — reprojected in Stage 2
#      File format: .bil rasters stored in data/PRISM_clean/
#      Naming convention from prism package:
#        PRISM_tmean_stable_4kmM3_YYYYMM_bil.bil  (4 km product)
#        PRISM_ppt_stable_4kmM3_YYYYMM_bil.bil
#        or 800 m: PRISM_tmean_stable_800mM3_YYYYMM_bil.bil
# -----------------------------------------------------------------------------
cat("  Loading PRISM rasters ...\n")

# List all temperature .bil files and sort by date
prism_temp_files <- list.files(
  "data/PRISM_clean",
  pattern    = "PRISM_tmean.*\\.bil$",
  full.names = TRUE,
  recursive  = TRUE       # search subdirectories (prism package nests files)
) %>% sort()

# List all precipitation .bil files
prism_precip_files <- list.files(
  "data/PRISM_clean",
  pattern    = "PRISM_ppt.*\\.bil$",
  full.names = TRUE,
  recursive  = TRUE
) %>% sort()

if (length(prism_temp_files) == 0)
  stop("No PRISM tmean .bil files found in data/PRISM_clean/
       Run the PRISM download block above first.")
if (length(prism_precip_files) == 0)
  stop("No PRISM ppt .bil files found in data/PRISM_clean/")

# Stack into multi-band rasters — one band per month
# Bands ordered chronologically (Jan 2012, Feb 2012, ..., Dec 2025)
prism_temp   <- rast(prism_temp_files)
prism_precip <- rast(prism_precip_files)

stopifnot(nlyr(prism_temp) == nlyr(prism_precip))
cat(sprintf("  PRISM: %d monthly temperature bands loaded.\n",
            nlyr(prism_temp)))
cat(sprintf("  PRISM: %d monthly precipitation bands loaded.\n",
            nlyr(prism_precip)))

# -----------------------------------------------------------------------------
# 1.4  NDVI — GEE Export GeoTIFF
#      Source     : Google Earth Engine export (Landsat 7/8/9 composite)
#      Sensors    : L7 (2012), L8 (2013–2021), L9 (2021–2025)
#      Season     : June–September fire season median composite
#      Period     : 2012–2025
#      Resolution : 30 m | CRS: EPSG:32612 (exported at this CRS from GEE)
#      Units      : Unitless (-1 to 1)
#      File       : data/Landsat/NDVI_Teton_FireSeason.tif
# -----------------------------------------------------------------------------
cat("  Loading NDVI GeoTIFF from GEE export ...\n")

ndvi_gee <- rast("data/Landsat/NDVI_Teton_FireSeason.tif")

# Verify CRS — GEE was set to export at EPSG:32612
if (!grepl("32612", crs(ndvi_gee, describe = TRUE)$code)) {
  warning("NDVI CRS is not EPSG:32612 — will reproject in Stage 2")
}
cat(sprintf("  NDVI loaded: resolution = %.0f m, CRS = %s\n",
            res(ndvi_gee)[1],
            crs(ndvi_gee, describe = TRUE)$code))

# -----------------------------------------------------------------------------
# 1.5  Daymet Solar Radiation
#      Source     : daymet.ornl.gov via daymetr package
#      Variable   : srad — surface downwelling shortwave radiation
#      Units      : W/m² (daily average; monthly averages downloaded)
#      Resolution : 1 km | CRS: Daymet Lambert Conformal Conic
#      Coverage   : Teton County bounding box
#      Note       : Download skipped if files already exist
# -----------------------------------------------------------------------------
daymet_files_existing <- list.files("data/Daymet",
                                     pattern = "srad.*\\.nc$",
                                     full.names = TRUE)

if (length(daymet_files_existing) == 0) {
  cat("  Downloading Daymet solar radiation for Teton County ...\n")
  download_daymet_ncss(
    location  = TETON_BBOX,
    start     = 2012,
    end       = 2025,
    param     = "srad",
    frequency = "monthly",
    path      = "data/Daymet"
  )
  daymet_files_existing <- list.files("data/Daymet",
                                       pattern = "srad.*\\.nc$",
                                       full.names = TRUE)
} else {
  cat(sprintf("  Daymet files already exist (%d files) — skipping download.\n",
              length(daymet_files_existing)))
}

daymet_srad <- rast(daymet_files_existing)
cat(sprintf("  Daymet srad loaded: %d bands.\n", nlyr(daymet_srad)))

# -----------------------------------------------------------------------------
# 1.6  SPEI Drought Indices — GEE Export GeoTIFFs
#      Source     : GEE (CSIC/SPEI/2_10) exported to Google Drive
#      Scales     : 1, 3, 6, 12 months
#      Period     : January 2012 – December 2023
#      Units      : Unitless (z-values / standardised index)
#      Format     : Multi-band GeoTIFF — 144 bands (one per month)
#      CRS        : EPSG:4326 (exported at this CRS from GEE)
#      Files      : data/SPEI/SPEI_*month_Teton.tif
# -----------------------------------------------------------------------------
cat("  Loading SPEI GeoTIFFs from GEE export ...\n")

spei_1m  <- rast("data/SPEI/SPEI_1month_Teton.tif")
spei_3m  <- rast("data/SPEI/SPEI_3month_Teton.tif")
spei_6m  <- rast("data/SPEI/SPEI_6month_Teton.tif")
spei_12m <- rast("data/SPEI/SPEI_12month_Teton.tif")

cat(sprintf("  SPEI loaded: 1m=%d bands | 3m=%d bands | 6m=%d bands | 12m=%d bands\n",
            nlyr(spei_1m), nlyr(spei_3m),
            nlyr(spei_6m), nlyr(spei_12m)))
# Each should have 144 bands (Jan 2012 – Dec 2023)

cat("Stage 1 complete.\n")


# =============================================================================
# STAGE 2 — PREPROCESSING
# =============================================================================
cat("\n===== STAGE 2: PREPROCESSING =====\n")

# -----------------------------------------------------------------------------
# 2.1  VIIRS Quality Filtering
# -----------------------------------------------------------------------------

# Step A: Fire season filter — May through October only
# Excludes December–April: snow cover makes wildland fire essentially
# impossible in Teton County during these months. Off-season VIIRS
# detections here are more likely to be residential heating anomalies
# or sensor artefacts than vegetation fires.
fire_season <- fire_raw %>%
  filter(month %in% FIRE_SEASON_MONTHS)
cat(sprintf("  Fire season filter (May–Oct): %d detections\n",
            nrow(fire_season)))

# Step B: Effective start date for S-NPP
# First day (19 Jan 2012) was during sensor calibration stabilisation
# — User Guide §4 recommends starting from 20 January 2012
fire_season <- fire_season %>%
  filter(!(platform == "SNPP" & acq_date < ymd("2012-01-20")))

# Step C: Confidence filter — retain Nominal ("N") + High ("H") only
# Low ("L") = class 7 = sun glint / <15 K temperature anomaly — spurious
# High ("H") = class 9 = saturated pixels — valid but intense fires
# Nominal ("N") = class 8 = strong valid detections (>15 K anomaly)
# Do NOT filter to "H" only — that removes all non-saturated valid fires
fire_conf <- fire_season %>%
  filter(confidence %in% c("N", "H"))
cat(sprintf("  Confidence filter (N+H retained): %d detections\n",
            nrow(fire_conf)))

# Step D: Vegetation fire type = 0 only
# Removes: volcanoes (1), industrial/static sources (2), offshore (3)
# This filter is reliable ONLY with FIRMS Collection 2 data
fire_veg <- fire_conf %>% filter(type == 0)
cat(sprintf("  Type=0 filter (vegetation fires only): %d detections\n",
            nrow(fire_veg)))

# Step E: Flag null FRP values (= 0 MW)
# Occur under extreme I4 channel saturation (User Guide §3.1.3)
# Retained for count-based analysis; excluded from FRP-weighted analysis
fire_veg <- fire_veg %>%
  mutate(frp_valid = frp > 0 & !is.na(frp))
cat(sprintf("  FRP-valid detections: %d / %d\n",
            sum(fire_veg$frp_valid), nrow(fire_veg)))

# Step F: Bowtie deduplication
# Co-located 375-m pixels near swath edges can share FRP values (§3.1.3)
# Within each platform + date + overpass group, retain the highest-FRP
# detection among any pair of points closer than 375 m apart
dedup_viirs <- function(sf_pts, dist_m = 375) {
  # Project to metric CRS for accurate distance calculation
  pts_prj <- st_transform(sf_pts, TARGET_CRS)

  pts_prj <- pts_prj %>%
    mutate(group_key = paste(platform,
                             format(acq_date, "%Y%m%d"),
                             acq_time,
                             sep = "_"))

  groups  <- group_split(pts_prj, group_key)

  cleaned <- map_dfr(groups, function(grp) {
    if (nrow(grp) == 1) return(grp)
    d    <- st_distance(grp)        # distances in metres
    diag(d) <- Inf
    keep <- rep(TRUE, nrow(grp))
    for (i in seq_len(nrow(grp))) {
      if (!keep[i]) next
      too_close <- which(as.numeric(d[i, ]) < dist_m)
      if (length(too_close) == 0) next
      candidates <- c(i, too_close)
      best <- candidates[which.max(grp$frp[candidates])]
      keep[setdiff(candidates, best)] <- FALSE
    }
    grp[keep, ]
  })
  st_transform(cleaned, "EPSG:4326")   # return to geographic
}

fire_dedup <- dedup_viirs(fire_veg, dist_m = 375)
cat(sprintf("  After bowtie deduplication: %d detections\n",
            nrow(fire_dedup)))

# Document platform temporal coverage
platform_coverage <- data.frame(
  platform = c("SNPP",          "NOAA20",          "NOAA21"),
  start    = c("2012-01-20",    "2017-01-01",       "2022-01-01"),
  end      = c("present",       "present",          "present"),
  note     = c("Effective start per VIIRS User Guide §4",
               "Pre-2017 record relies on S-NPP only",
               "Added 2022 — three platforms from this year")
)
write_csv(platform_coverage, "output/tables/platform_coverage.csv")

write_csv(st_drop_geometry(fire_dedup),
          "intermediate/fire_filtered_all_platforms.csv")

# -----------------------------------------------------------------------------
# 2.2  Reproject All Layers to EPSG:32612 (UTM Zone 12N)
#      All distance and area calculations now in metres
# -----------------------------------------------------------------------------
cat("  Reprojecting all layers to EPSG:32612 ...\n")

study_area_prj   <- st_transform(study_area,  TARGET_CRS)
fire_points_prj  <- st_transform(fire_dedup,  TARGET_CRS)
prism_temp_prj   <- project(prism_temp,   TARGET_CRS, method = "bilinear")
prism_precip_prj <- project(prism_precip, TARGET_CRS, method = "bilinear")
daymet_srad_prj  <- project(daymet_srad,  TARGET_CRS, method = "bilinear")
spei_1m_prj      <- project(spei_1m,      TARGET_CRS, method = "bilinear")
spei_3m_prj      <- project(spei_3m,      TARGET_CRS, method = "bilinear")
spei_6m_prj      <- project(spei_6m,      TARGET_CRS, method = "bilinear")
spei_12m_prj     <- project(spei_12m,     TARGET_CRS, method = "bilinear")

# NDVI — already exported at EPSG:32612 from GEE; verify then skip if correct
if (!grepl("32612", crs(ndvi_gee, describe = TRUE)$code)) {
  ndvi_prj <- project(ndvi_gee, TARGET_CRS, method = "bilinear")
  cat("  NDVI reprojected to EPSG:32612.\n")
} else {
  ndvi_prj <- ndvi_gee
  cat("  NDVI already in EPSG:32612 — no reprojection needed.\n")
}

# Verify vector CRS
stopifnot(st_crs(study_area_prj)$epsg  == 32612)
stopifnot(st_crs(fire_points_prj)$epsg == 32612)
cat("  CRS verification passed for all layers.\n")

# -----------------------------------------------------------------------------
# 2.3  Clip All Layers to Teton County
# -----------------------------------------------------------------------------
cat("  Clipping all layers to Teton County boundary ...\n")

wy_vect <- vect(study_area_prj)   # terra SpatVector for raster operations

# Clip fire points to county (spatial intersection)
fire_points_clip <- st_intersection(fire_points_prj,
                                    st_union(study_area_prj))
cat(sprintf("  Fire detections within Teton County: %d\n",
            nrow(fire_points_clip)))

# Clip all rasters
prism_temp_clip   <- mask(crop(prism_temp_prj,   wy_vect), wy_vect)
prism_precip_clip <- mask(crop(prism_precip_prj, wy_vect), wy_vect)
ndvi_clip         <- mask(crop(ndvi_prj,         wy_vect), wy_vect)
daymet_srad_clip  <- mask(crop(daymet_srad_prj,  wy_vect), wy_vect)
spei_1m_clip      <- mask(crop(spei_1m_prj,     wy_vect), wy_vect)
spei_3m_clip      <- mask(crop(spei_3m_prj,     wy_vect), wy_vect)
spei_6m_clip      <- mask(crop(spei_6m_prj,     wy_vect), wy_vect)
spei_12m_clip     <- mask(crop(spei_12m_prj,    wy_vect), wy_vect)

# -----------------------------------------------------------------------------
# 2.4  Resample All Rasters to Reference Grid (PRISM 800 m)
# -----------------------------------------------------------------------------
cat("  Resampling all rasters to 800 m reference grid ...\n")

reference_grid <- prism_temp_clip[[1]]   # spatial template — 800 m EPSG:32612

# NDVI: 30 m → 800 m (mean aggregation preserves spatial average)
ndvi_res      <- resample(ndvi_clip,        reference_grid, method = "bilinear")

# Daymet: 1 km → 800 m
daymet_res    <- resample(daymet_srad_clip, reference_grid, method = "bilinear")

# SPEI: ~55 km → 800 m (upscaling by interpolation)
# Note: SPEI at 0.5° resolution has limited spatial variation across Teton Co.
# Temporal variation (month to month) is the primary signal being captured.
spei_1m_res   <- resample(spei_1m_clip,    reference_grid, method = "bilinear")
spei_3m_res   <- resample(spei_3m_clip,    reference_grid, method = "bilinear")
spei_6m_res   <- resample(spei_6m_clip,    reference_grid, method = "bilinear")
spei_12m_res  <- resample(spei_12m_clip,   reference_grid, method = "bilinear")

# Verify all rasters at 800 m
all_rasters <- list(prism_temp_clip, prism_precip_clip,
                    ndvi_res, daymet_res,
                    spei_1m_res, spei_3m_res, spei_6m_res, spei_12m_res)
res_check <- sapply(all_rasters, function(r) round(res(r)[1]))
if (!all(res_check == TARGET_RES))
  warning(sprintf("Resolution mismatch detected: %s",
                  paste(res_check, collapse = ", ")))
cat(sprintf("  All rasters verified at %d m resolution.\n", TARGET_RES))

# -----------------------------------------------------------------------------
# 2.5  Platform Normalisation
#      Detection counts increase as platforms are added (2017, 2022).
#      Assign weight = 1 / n_active_platforms per year to allow
#      fair inter-annual trend comparisons.
# -----------------------------------------------------------------------------
platform_counts <- st_drop_geometry(fire_points_clip) %>%
  count(year, platform) %>%
  pivot_wider(names_from  = platform,
              values_from = n,
              values_fill = 0) %>%
  mutate(n_platforms = rowSums(select(., -year) > 0))

fire_points_clip <- fire_points_clip %>%
  left_join(select(platform_counts, year, n_platforms), by = "year") %>%
  mutate(norm_weight = 1 / n_platforms)

write_csv(platform_counts,
          "output/tables/platform_detection_counts_by_year.csv")
cat("  Platform normalisation weights assigned.\n")
cat("Stage 2 complete.\n")


# =============================================================================
# STAGE 3 — ANALYSIS
# =============================================================================

# =============================================================================
# MODULE 3.1 — WILDFIRE FREQUENCY & HOTSPOT ANALYSIS
# =============================================================================
cat("\n===== MODULE 3.1: HOTSPOT ANALYSIS =====\n")

# --- 3.1.1  Split detections by year and season -----------------------------
fire_annual   <- split(fire_points_clip, fire_points_clip$year)
fire_seasonal <- list(
  fire_season = fire_points_clip,   # full May–Oct season
  peak        = fire_points_clip %>% filter(month %in% PEAK_MONTHS),
  shoulder    = fire_points_clip %>%
                  filter(month %in% c(5, 10))  # May and October
)
years <- sort(unique(fire_points_clip$year))
cat(sprintf("  Detections across %d years (May–Oct fire season).\n",
            length(years)))

# --- 3.1.2  FRP Density Surfaces (MW/km², kernel density weighted by FRP) ---
cat("  Computing FRP density surfaces ...\n")

compute_frp_density <- function(pts, ref_raster, bw_m = KDE_RADIUS) {
  pts_valid <- pts %>% filter(frp_valid == TRUE)
  if (nrow(pts_valid) == 0) return(NULL)

  # Rasterize: sum FRP (MW) per 800 m cell
  r_frp <- rasterize(vect(pts_valid), ref_raster,
                     field = "frp", fun = "sum", background = 0)

  # Convert to MW/km² (cell area = 800 × 800 m = 0.64 km²)
  cell_area_km2 <- (TARGET_RES / 1000)^2
  r_density     <- r_frp / cell_area_km2

  # Gaussian smoothing — sigma in cell units
  sigma_cells <- round(bw_m / TARGET_RES)
  fw <- focalMat(ref_raster, sigma_cells, type = "Gauss")
  focal(r_density, w = fw, fun = "sum", na.rm = TRUE)
}

frp_density <- list()
for (yr in years) {
  pts_yr <- fire_annual[[as.character(yr)]]
  d_r    <- compute_frp_density(pts_yr, reference_grid)
  frp_density[[as.character(yr)]] <- d_r
  if (!is.null(d_r))
    writeRaster(d_r,
                paste0("output/rasters/FRP_density_", yr, ".tif"),
                overwrite = TRUE)
}
cat(sprintf("  FRP density surfaces (MW/km²) saved for %d years.\n",
            length(years)))

# --- 3.1.3  Getis-Ord Gi* Hotspot Analysis ----------------------------------
# Magnitude: FRP (MW) — weights by fire intensity not just occurrence
# Fallback : FRP = 0 → x = 1 (count contribution only)
# Neighbourhood: fixed distance band of 20 km (GI_DIST_M metres)

compute_hotspot_gi <- function(pts, dist_band_m = GI_DIST_M,
                               label = "output") {
  if (nrow(pts) < 10) {
    message(sprintf("  Skipping Gi* — fewer than 10 points: %s", label))
    return(NULL)
  }

  pts <- pts %>%
    mutate(x_mag = ifelse(frp_valid & frp > 0, frp, 1))

  coords <- st_coordinates(pts)

  # Fixed-distance spatial weights in metres (EPSG:32612)
  nb <- dnearneigh(coords, d1 = 0, d2 = dist_band_m)
  lw <- nb2listw(nb, style = "B", zero.policy = TRUE)

  gi  <- localG(pts$x_mag, lw, zero.policy = TRUE)

  pts_out <- pts %>%
    mutate(
      gi_star       = as.numeric(gi),
      p_value       = 2 * pnorm(-abs(gi_star)),
      hotspot_class = case_when(
        gi_star >=  2.576 & p_value < 0.01 ~ "hotspot_99%",
        gi_star >=  1.960 & p_value < 0.05 ~ "hotspot_95%",
        gi_star >=  1.645 & p_value < 0.10 ~ "hotspot_90%",
        gi_star <= -1.645                   ~ "cold_spot",
        TRUE                                ~ "not_significant"
      )
    )

  st_write(pts_out,
           paste0("output/vectors/hotspot_", label, ".gpkg"),
           delete_dsn = TRUE, quiet = TRUE)
  cat(sprintf("  Hotspot saved: %s\n", label))
  pts_out
}

# Annual hotspot layers
hotspot_annual <- list()
for (yr in years) {
  hotspot_annual[[as.character(yr)]] <-
    compute_hotspot_gi(fire_annual[[as.character(yr)]],
                       GI_DIST_M, label = paste0("annual_", yr))
}

# Full fire season pooled hotspot
hotspot_season <- compute_hotspot_gi(
  fire_seasonal$fire_season, GI_DIST_M, label = "fire_season_pooled")

# Peak months (Jun–Sep) hotspot
hotspot_peak <- compute_hotspot_gi(
  fire_seasonal$peak, GI_DIST_M, label = "peak_JunSep")

# --- 3.1.4  Temporal Hotspot Shift Map --------------------------------------
cat("  Computing temporal hotspot shift ...\n")
yr_first <- as.character(min(years))
yr_last  <- as.character(max(years))

if (!is.null(hotspot_annual[[yr_first]]) &&
    !is.null(hotspot_annual[[yr_last]])) {

  hs_first <- hotspot_annual[[yr_first]] %>%
    select(geometry, hotspot_class) %>%
    rename(class_first = hotspot_class)
  hs_last  <- hotspot_annual[[yr_last]]  %>%
    select(geometry, hotspot_class) %>%
    rename(class_last  = hotspot_class)

  hs_joined <- st_join(hs_first, hs_last, join = st_nearest_feature)

  hs_shift <- hs_joined %>%
    mutate(shift_class = case_when(
      grepl("hotspot", class_first) &
        class_last == "not_significant"  ~ "declining_hotspot",
      grepl("hotspot", class_last)  &
        class_first == "not_significant" ~ "emerging_hotspot",
      grepl("hotspot", class_first) &
        grepl("hotspot", class_last)     ~ "persistent_hotspot",
      TRUE                               ~ "no_significant_shift"
    ))

  st_write(hs_shift, "output/vectors/temporal_hotspot_shift.gpkg",
           delete_dsn = TRUE, quiet = TRUE)
  cat("  Temporal hotspot shift map saved.\n")
}

cat("Module 3.1 complete.\n")


# =============================================================================
# MODULE 3.2 — CLIMATIC AND SURFACE CONDITION ANALYSIS
# =============================================================================
cat("\n===== MODULE 3.2: CLIMATE & SURFACE ANALYSIS =====\n")

# --- 3.2.1  PRISM Seasonal Aggregation (°C and mm) --------------------------
aggregate_by_season <- function(rast_stack, months_include) {
  mean(rast_stack[[months_include]], na.rm = TRUE)
}

temp_peak_mean   <- aggregate_by_season(prism_temp_clip,   PEAK_MONTHS)
precip_peak_mean <- aggregate_by_season(prism_precip_clip, PEAK_MONTHS)

# Long-term means for anomaly calculation (baseline = full record)
temp_lt_mean   <- mean(prism_temp_clip,   na.rm = TRUE)
precip_lt_mean <- mean(prism_precip_clip, na.rm = TRUE)

temp_anomaly   <- temp_peak_mean   - temp_lt_mean    # °C
precip_anomaly <- precip_peak_mean - precip_lt_mean  # mm

writeRaster(temp_peak_mean,   "output/rasters/temp_peak_mean_C.tif",     overwrite=TRUE)
writeRaster(precip_peak_mean, "output/rasters/precip_peak_mean_mm.tif",  overwrite=TRUE)
writeRaster(temp_anomaly,     "output/rasters/temp_anomaly_C.tif",       overwrite=TRUE)
writeRaster(precip_anomaly,   "output/rasters/precip_anomaly_mm.tif",    overwrite=TRUE)
cat("  PRISM seasonal means (°C, mm) and anomalies saved.\n")

# --- 3.2.2  Daymet Solar Radiation Peak Season Mean (W/m²) ------------------
daymet_peak_mean <- aggregate_by_season(daymet_res, PEAK_MONTHS)
writeRaster(daymet_peak_mean,
            "output/rasters/srad_peak_mean_Wm2.tif", overwrite = TRUE)
cat("  Daymet solar radiation peak mean (W/m²) saved.\n")

# --- 3.2.3  NDVI (unitless, from GEE GeoTIFF) -------------------------------
# GEE already computed the median composite — ndvi_res is ready to use
# Apply valid range masking
cat("  Processing NDVI ...\n")
ndvi_raster <- ndvi_res
ndvi_raster[ndvi_raster < -1] <- NA
ndvi_raster[ndvi_raster >  1] <- NA
writeRaster(ndvi_raster, "output/rasters/NDVI_peak_season.tif",
            overwrite = TRUE)

# NDVI zonal statistics: fire-affected vs. non-fire areas
# Buffer = 1600 m (2 × 800 m pixel size — captures fire pixel footprint)
fire_buf_v  <- vect(st_buffer(fire_points_clip, 1600))
fire_buf_u  <- aggregate(fire_buf_v)
nonfire_msk <- erase(wy_vect, fire_buf_u)

ndvi_fire_stats    <- global(mask(ndvi_raster, fire_buf_u),
                              c("mean","sd","min","max"), na.rm = TRUE)
ndvi_nonfire_stats <- global(mask(ndvi_raster, nonfire_msk),
                              c("mean","sd","min","max"), na.rm = TRUE)

write_csv(as.data.frame(ndvi_fire_stats),    "output/tables/NDVI_fire_stats.csv")
write_csv(as.data.frame(ndvi_nonfire_stats), "output/tables/NDVI_nonfire_stats.csv")
cat("  NDVI processed and zonal statistics saved.\n")

# --- 3.2.4  SPEI Drought Severity Classification (unitless) -----------------
cat("  Classifying SPEI drought severity ...\n")

classify_spei <- function(spei_r, scale_label) {
  # McKee et al. (1993) / Beguería & Vicente-Serrano (2017) thresholds
  rcl <- matrix(
    c(-Inf, -2.0,  5,   # extreme drought
      -2.0, -1.5,  4,   # severe drought
      -1.5, -1.0,  3,   # moderate drought
      -1.0, -0.5,  2,   # mild drought
      -0.5,  Inf,  1),  # near normal
    ncol = 3, byrow = TRUE
  )
  cls_r <- classify(spei_r, rcl, include.lowest = TRUE)
  levels(cls_r) <- data.frame(
    value = 1:5,
    label = c("near_normal", "mild_drought", "moderate_drought",
              "severe_drought", "extreme_drought")
  )
  writeRaster(cls_r,
              paste0("output/rasters/SPEI_", scale_label, "_classified.tif"),
              overwrite = TRUE)
  cat(sprintf("  SPEI %s classified and saved.\n", scale_label))
  cls_r
}

spei_1m_class  <- classify_spei(spei_1m_res,  "1m")
spei_3m_class  <- classify_spei(spei_3m_res,  "3m")
spei_6m_class  <- classify_spei(spei_6m_res,  "6m")
spei_12m_class <- classify_spei(spei_12m_res, "12m")

cat("Module 3.2 complete.\n")


# =============================================================================
# MODULE 3.3 — SPATIAL REGRESSION (GWR / GWLR)
# =============================================================================
cat("\n===== MODULE 3.3: GWR / GWLR =====\n")

# --- 3.3.1  Build Regression Input Table ------------------------------------
cat("  Building regression input table ...\n")

# Rasterize FRP: sum per 800 m cell (MW), FRP-valid only
fire_frpvalid  <- fire_points_clip %>% filter(frp_valid == TRUE)
frp_sum_raster <- rasterize(vect(fire_frpvalid), reference_grid,
                             field = "frp", fun = "sum", background = 0)
writeRaster(frp_sum_raster, "output/rasters/FRP_sum_per_cell_MW.tif",
            overwrite = TRUE)

# Binary fire presence (1/0)
fire_pres_raster <- ifel(frp_sum_raster > 0, 1, 0)

# Predictor stack — all 800 m, EPSG:32612, SI units
pred_stack <- c(
  temp_peak_mean,       # °C
  precip_peak_mean,     # mm
  daymet_peak_mean,     # W/m²
  ndvi_raster,          # unitless
  spei_3m_res[[1]],     # unitless — 3-month drought lag
  spei_12m_res[[1]]     # unitless — 12-month cumulative drought
)
names(pred_stack) <- c("temp_C", "precip_mm", "srad_Wm2",
                        "ndvi", "spei_3m", "spei_12m")

full_stack <- c(pred_stack, frp_sum_raster, fire_pres_raster)
names(full_stack)[7:8] <- c("frp_sum_MW", "fire_presence")

# One row per valid 800 m cell
regression_df <- as.data.frame(full_stack, xy = TRUE, na.rm = TRUE)
names(regression_df)[1:2] <- c("x_m", "y_m")   # metres — EPSG:32612
regression_df <- na.omit(regression_df)

cat(sprintf("  Regression input cells: %d\n", nrow(regression_df)))
write_csv(regression_df, "output/tables/regression_input.csv")

# --- 3.3.2  OLS Pre-check: Spatial Autocorrelation in Residuals -------------
cat("  OLS pre-check — Moran's I on residuals ...\n")

ols_model <- lm(
  frp_sum_MW ~ temp_C + precip_mm + srad_Wm2 + ndvi + spei_3m,
  data = regression_df
)
print(summary(ols_model))

coords_m  <- as.matrix(regression_df[, c("x_m", "y_m")])
nb_k8     <- knn2nb(knearneigh(coords_m, k = 8))
lw_k8     <- nb2listw(nb_k8, style = "W", zero.policy = TRUE)

moran_res <- moran.test(residuals(ols_model), lw_k8, zero.policy = TRUE)
cat(sprintf("  Moran's I = %.4f | p-value = %.4f\n",
            moran_res$statistic["Moran I statistic"],
            moran_res$p.value))

if (moran_res$p.value < 0.05) {
  cat("  Significant spatial autocorrelation — GWR/GWLR justified.\n")
} else {
  cat("  WARNING: Weak autocorrelation — consider whether GWR is appropriate.\n")
}

# --- 3.3.3  Geographically Weighted Regression (GWR) -----------------------
# Dependent: frp_sum_MW (continuous, MW)
cat("  Running GWR (adaptive bandwidth, AICc selection) ...\n")

regression_sp <- SpatialPointsDataFrame(
  coords      = coords_m,
  data        = regression_df,
  proj4string = CRS("+init=epsg:32612")
)

gwr_bw <- bw.gwr(
  formula  = frp_sum_MW ~ temp_C + precip_mm + srad_Wm2 + ndvi + spei_3m,
  data     = regression_sp,
  approach = "AICc",
  kernel   = "gaussian",
  adaptive = TRUE
)
cat(sprintf("  Optimal GWR bandwidth: %d neighbours\n", gwr_bw))

gwr_result <- gwr.basic(
  formula   = frp_sum_MW ~ temp_C + precip_mm + srad_Wm2 + ndvi + spei_3m,
  data      = regression_sp,
  bw        = gwr_bw,
  kernel    = "gaussian",
  adaptive  = TRUE,
  F123.test = FALSE
)
print(gwr_result)

gwr_sf     <- st_as_sf(gwr_result$SDF)
predictors <- c("temp_C", "precip_mm", "srad_Wm2", "ndvi", "spei_3m")

for (pred in predictors) {
  coef_r <- rasterize(vect(gwr_sf), reference_grid,
                      field = pred, fun = "mean")
  writeRaster(coef_r,
              paste0("output/rasters/GWR_coef_", pred, ".tif"),
              overwrite = TRUE)
}
r2_r <- rasterize(vect(gwr_sf), reference_grid,
                  field = "Local_R2", fun = "mean")
writeRaster(r2_r, "output/rasters/GWR_local_R2.tif", overwrite = TRUE)
cat("  GWR coefficient and local R² rasters saved.\n")

# --- 3.3.4  Geographically Weighted Logistic Regression (GWLR) -------------
# Dependent: fire_presence (binary 1/0)
cat("  Running GWLR (adaptive bandwidth, AICc selection) ...\n")

gwlr_bw <- bw.ggwr(
  formula  = fire_presence ~ temp_C + precip_mm + srad_Wm2 + ndvi + spei_3m,
  data     = regression_sp,
  family   = "binomial",
  approach = "AICc",
  kernel   = "gaussian",
  adaptive = TRUE
)
cat(sprintf("  Optimal GWLR bandwidth: %d neighbours\n", gwlr_bw))

gwlr_result <- ggwr.basic(
  formula  = fire_presence ~ temp_C + precip_mm + srad_Wm2 + ndvi + spei_3m,
  data     = regression_sp,
  family   = "binomial",
  bw       = gwlr_bw,
  kernel   = "gaussian",
  adaptive = TRUE
)
print(gwlr_result)

gwlr_sf <- st_as_sf(gwlr_result$SDF)
prob_r  <- rasterize(vect(gwlr_sf), reference_grid,
                     field = "yhat", fun = "mean")
writeRaster(prob_r, "output/rasters/GWLR_fire_probability.tif",
            overwrite = TRUE)

for (pred in predictors) {
  lo_r <- rasterize(vect(gwlr_sf), reference_grid,
                    field = pred, fun = "mean")
  writeRaster(lo_r,
              paste0("output/rasters/GWLR_logodds_", pred, ".tif"),
              overwrite = TRUE)
}
cat("  GWLR probability and log-odds rasters saved.\n")

# --- 3.3.5  Dominant Climate Driver Per Cell --------------------------------
cat("  Mapping dominant climate driver ...\n")

coef_stack   <- rast(paste0("output/rasters/GWR_coef_",
                             predictors, ".tif"))
names(coef_stack) <- predictors
driver_index <- which.max(abs(coef_stack))
levels(driver_index) <- data.frame(value = 1:5, label = predictors)
writeRaster(driver_index,
            "output/rasters/dominant_climate_driver.tif",
            overwrite = TRUE)
cat("  Dominant driver raster saved.\n")
cat("Module 3.3 complete.\n")


# =============================================================================
# STAGE 4 — MAP OUTPUT
# CRS: EPSG:32612 | Scale bars: km | Resolution: 300 dpi
# =============================================================================
cat("\n===== STAGE 4: MAP OUTPUT =====\n")
tmap_mode("plot")

# Reload final rasters
gwr_r2    <- rast("output/rasters/GWR_local_R2.tif")
gwlr_prob <- rast("output/rasters/GWLR_fire_probability.tif")
driver_r  <- rast("output/rasters/dominant_climate_driver.tif")
ndvi_r    <- rast("output/rasters/NDVI_peak_season.tif")
spei3_cls <- rast("output/rasters/SPEI_3m_classified.tif")
frp_last  <- frp_density[[as.character(max(years))]]

# Shared elements
scale_bar  <- tm_scale_bar(breaks    = c(0, 10, 20, 30),
                            text.size = 0.6,
                            unit      = "km",
                            position  = c("left", "bottom"))
compass    <- tm_compass(position = c("right", "top"), size = 1.5)
border_lyr <- tm_shape(study_area_prj) +
                tm_borders(lwd = 1.8, col = "black")

# --- Map 1: FRP-Weighted Hotspot Map ----------------------------------------
if (!is.null(hotspot_peak)) {
  map1 <- tm_shape(study_area_prj) + tm_fill(col = "grey92") +
    tm_shape(hotspot_peak %>%
               filter(hotspot_class != "not_significant")) +
    tm_dots(col     = "hotspot_class",
            palette = c("hotspot_99%" = "#67000D",
                        "hotspot_95%" = "#CB181D",
                        "hotspot_90%" = "#FC8D59",
                        "cold_spot"   = "#2166AC"),
            size    = 0.08,
            title   = "Gi* Class (FRP-weighted)") +
    border_lyr + scale_bar + compass +
    tm_layout(title          = "FRP-Weighted Wildfire Hotspots\nTeton County — Fire Season",
              legend.outside = TRUE)
  tmap_save(map1, "output/maps/Map1_FRP_Hotspot.png",
            dpi = 300, width = 10, height = 8)
  cat("  Map 1 saved.\n")
}

# --- Map 2: FRP Density Surface (MW/km²) ------------------------------------
if (!is.null(frp_last)) {
  map2 <- tm_shape(frp_last) +
    tm_raster(style   = "quantile", n = 7,
              palette = "YlOrRd",
              title   = "FRP Density (MW/km²)") +
    border_lyr + scale_bar + compass +
    tm_layout(title          = paste0("Fire Radiative Power Density\nTeton County — ",
                                      max(years)),
              legend.outside = TRUE)
  tmap_save(map2, "output/maps/Map2_FRP_Density.png",
            dpi = 300, width = 10, height = 8)
  cat("  Map 2 saved.\n")
}

# --- Map 3: NDVI Condition (unitless) ---------------------------------------
map3 <- tm_shape(ndvi_r) +
  tm_raster(breaks  = c(-1, 0.2, 0.5, 1),
            labels  = c("Stressed / Bare (<0.2)",
                        "Moderate (0.2–0.5)",
                        "Healthy (>0.5)"),
            palette = c("#8C510A", "#F6E8C3", "#1B7837"),
            title   = "NDVI (unitless)") +
  border_lyr + scale_bar + compass +
  tm_layout(title          = "NDVI Vegetation Condition\nTeton County — Jun–Sep 2012–2025",
            legend.outside = TRUE)
tmap_save(map3, "output/maps/Map3_NDVI.png",
          dpi = 300, width = 10, height = 8)
cat("  Map 3 saved.\n")

# --- Map 4: SPEI 3-Month Drought Severity -----------------------------------
map4 <- tm_shape(spei3_cls) +
  tm_raster(palette = c("near_normal"      = "#F7F7F7",
                         "mild_drought"     = "#FEE391",
                         "moderate_drought" = "#FEC44F",
                         "severe_drought"   = "#FE9929",
                         "extreme_drought"  = "#CC4C02"),
            title   = "SPEI Drought Class") +
  tm_shape(fire_points_clip) +
    tm_dots(col = "black", size = 0.03, alpha = 0.5) +
  border_lyr + scale_bar + compass +
  tm_layout(title          = "SPEI 3-Month Drought Severity\nTeton County",
            legend.outside = TRUE)
tmap_save(map4, "output/maps/Map4_SPEI_3m.png",
          dpi = 300, width = 10, height = 8)
cat("  Map 4 saved.\n")

# --- Map 5: GWR Local R² ----------------------------------------------------
map5 <- tm_shape(gwr_r2) +
  tm_raster(style   = "cont",
            palette = "Blues",
            title   = "Local R² (unitless)") +
  border_lyr + scale_bar + compass +
  tm_layout(title          = "GWR Local Model Fit (R²)\nTeton County",
            legend.outside = TRUE)
tmap_save(map5, "output/maps/Map5_GWR_R2.png",
          dpi = 300, width = 10, height = 8)
cat("  Map 5 saved.\n")

# --- Map 6: Dominant Climate Driver -----------------------------------------
map6 <- tm_shape(driver_r) +
  tm_raster(palette = c(temp_C    = "#D73027",
                         precip_mm = "#4575B4",
                         srad_Wm2  = "#F46D43",
                         ndvi      = "#1A9850",
                         spei_3m   = "#762A83"),
            title   = "Dominant Driver") +
  border_lyr + scale_bar + compass +
  tm_layout(title          = "Dominant Climate Driver of Fire Activity\nTeton County",
            legend.outside = TRUE)
tmap_save(map6, "output/maps/Map6_DominantDriver.png",
          dpi = 300, width = 10, height = 8)
cat("  Map 6 saved.\n")

# --- Map 7: GWLR Fire Probability (0–1) -------------------------------------
map7 <- tm_shape(gwlr_prob) +
  tm_raster(style   = "cont",
            palette = c("#FFFFCC", "#FEB24C", "#FC4E2A", "#800026"),
            title   = "P(Fire) — unitless") +
  border_lyr + scale_bar + compass +
  tm_layout(title          = "Predicted Fire Probability — GWLR\nTeton County",
            legend.outside = TRUE)
tmap_save(map7, "output/maps/Map7_GWLR_Probability.png",
          dpi = 300, width = 10, height = 8)
cat("  Map 7 saved.\n")

cat("\n===== FULL PIPELINE COMPLETE =====\n")
cat("Outputs:\n")
cat("  output/maps/     — 7 map PNGs at 300 dpi\n")
cat("  output/rasters/  — all GeoTIFF outputs\n")
cat("  output/vectors/  — hotspot GeoPackages\n")
cat("  output/tables/   — CSV summary tables\n")
