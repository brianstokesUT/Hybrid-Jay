# eBird processing + Maxent modeling for GRJA and BLJA
# Author: Brian R. Stokes
# Description: Clean pipeline to ingest eBird sampling/observation data, filter and thin
# occurrences, fit Maxent models on current climate, evaluate against background from
# non-detection sites, and project to a future climate stack.
#
# This script assumes you have already run `bioclim_dataprep.R`, which
# creates two climate stacks:
# cur : SpatRaster of current (1991–2020) variables (EPSG:4326)
# fut : SpatRaster of future (SSP245 2041–2060) variables (EPSG:4326)
# and defines the study-area extent `crop_window` (SpatExtent).

wd<-("~/PATH")
setwd(wd)

set.seed(223)
options(java.parameters = "-Xmx8000m")

suppressPackageStartupMessages({
  library(terra)
  library(raster)
  library(dismo)
  library(rJava)
  library(auk)
  library(dplyr)
  library(readr)
  library(lubridate)
  library(sf)
  library(ggplot2)
  library(GeoThinneR)
  library(ggpattern)
})

# Verify objects from previous script are present
stopifnot(exists("cur"), exists("fut"), exists("crop_window"))

# --- helpers ---
sanitize_names <- function(df) {
  nm <- names(df)
  nm <- trimws(nm)
  nm <- gsub("[^A-Za-z0-9_]+", "_", nm)  # kill stray chars/spaces
  names(df) <- nm
  df
}

prep_for_thinning <- function(df, crop_window,
                              lon_candidates = c("longitude","lon","LONGITUDE","longitude_x","lon_x"),
                              lat_candidates = c("latitude","lat","LATITUDE","latitude_x","lat_x")) {
  # drop geometry if sf
  if (inherits(df, "sf")) {
    coords <- sf::st_coordinates(df)
    df <- sf::st_drop_geometry(df)
    df$longitude <- coords[,1]
    df$latitude  <- coords[,2]
  }
  
  df <- sanitize_names(df)
  
  # find lon/lat columns by case-insensitive match
  lower <- tolower(names(df))
  lon_idx <- match(tolower(lon_candidates), lower)
  lat_idx <- match(tolower(lat_candidates), lower)
  lon_idx <- lon_idx[!is.na(lon_idx)][1]
  lat_idx <- lat_idx[!is.na(lat_idx)][1]
  if (is.na(lon_idx) || is.na(lat_idx)) {
    stop("Could not find lon/lat columns after sanitizing names.")
  }
  
  # normalize column names to lon/lat
  names(df)[lon_idx] <- "lon"
  names(df)[lat_idx] <- "lat"
  
  # coerce to numeric (auk sometimes yields character in weird cases)
  df$lon <- suppressWarnings(as.numeric(df$lon))
  df$lat <- suppressWarnings(as.numeric(df$lat))
  
  # toss rows where coercion produced NA
  df <- df[is.finite(df$lon) & is.finite(df$lat), , drop = FALSE]
  
  # crop to study window
  df <- dplyr::filter(
    df,
    lon >= terra::xmin(crop_window) & lon <= terra::xmax(crop_window) &
      lat >= terra::ymin(crop_window) & lat <= terra::ymax(crop_window)
  )
  
  df
}

# -----------------------
# 1) Prep climate layers for dismo
# -----------------------
cur_rs <- raster::stack(cur)
fut_rs <- raster::stack(fut)

# Coarse grid template for thinning (~factor of 10 relative to native res)
# NOTE: factor=10 is 10x cell size, not 10 km; adjust as needed.
raster_grid_coarse <- raster::aggregate(cur_rs[[1]], fact = 10)

# -----------------------
# 2) Helpers
# -----------------------
# Crop lon/lat by crop_window
.crop_pts <- function(df, x = "longitude", y = "latitude") {
  dplyr::filter(
    df,
    .data[[x]] >= terra::xmin(crop_window) & .data[[x]] <= terra::xmax(crop_window) &
      .data[[y]] >= terra::ymin(crop_window) & .data[[y]] <= terra::ymax(crop_window)
  )
}
# Ensure plain lon/lat columns exist (robust to sf or renamed columns after joins)
to_lonlat_df <- function(x, lon_col = NULL, lat_col = NULL) {
  if (inherits(x, "sf")) {
    coords <- sf::st_coordinates(x)
    x <- sf::st_drop_geometry(x)
    x$longitude <- coords[,1]
    x$latitude  <- coords[,2]
    return(x)
  }
  cand_lon <- c(lon_col, "longitude","lon","LONGITUDE","longitude.x","lon.x")
  cand_lat <- c(lat_col, "latitude","lat","LATITUDE","latitude.x","lat.x")
  lon <- cand_lon[cand_lon %in% names(x)][1]
  lat <- cand_lat[cand_lat %in% names(x)][1]
  if (is.na(lon) || is.na(lat)) stop("Could not find longitude/latitude columns in data.")
  names(x)[match(c(lon,lat), names(x))] <- c("longitude","latitude")
  x
}

# -----------------------
# 3) GRJA (Green Jay): eBird → thin → MaxEnt
# -----------------------
# Sampling Datasets - GRJA
smp_mx <- "ebd_MX_grnjay_smp_relJul-2024/ebd_MX_grnjay_smp_relJul-2024_sampling.txt"
smp_tx <- "ebd_US-TX_grnjay_smp_relJul-2024/ebd_US-TX_grnjay_smp_relJul-2024_sampling.txt"

checklists_mx <- auk::read_sampling(smp_mx, unique = TRUE)
checklists_tx <- auk::read_sampling(smp_tx, unique = TRUE)
checklists    <- dplyr::bind_rows(checklists_mx, checklists_tx)

# Keep checklists through 2023-05-31 (pre-hybrid capture)
checklists_jun2023 <- checklists %>%
  dplyr::filter(observation_date >= as.Date("1900-01-01"), observation_date <= as.Date("2023-05-31")) %>%
  dplyr::arrange(observation_date)

# Observations - GRJA (MX + US)
ebd_grja_mx <- "ebd_MX_grnjay_smp_relJul-2024/ebd_MX_grnjay_smp_relJul-2024.txt"
ebd_grja_us <- "ebd_US_grnjay_smp_relJul-2024/ebd_US_grnjay_smp_relJul-2024.txt"

obs_grja_mx <- auk::read_ebd(ebd_grja_mx)
obs_grja_us <- auk::read_ebd(ebd_grja_us)
obs_grja    <- dplyr::bind_rows(obs_grja_mx, obs_grja_us)

obs_grja <- obs_grja %>%
  dplyr::filter(all_species_reported, dplyr::between(lubridate::year(observation_date), 1900, 2024)) %>%
  dplyr::semi_join(checklists_jun2023, by = "checklist_id")

# Observations through 2023-05-31
obs_grja_jun2023 <- obs_grja %>%
  dplyr::filter(all_species_reported,
                observation_date >= as.Date("1900-01-01"),
                observation_date <= as.Date("2023-05-31")) %>%
  dplyr::arrange(observation_date)

# Filter out vagrant-heavy localities; retain recurring localities
grja_presence_jun2023 <- obs_grja_jun2023 %>%
  dplyr::filter(!locality_id %in% c("L8812047", "L13398708")) %>%                 # specific vagrancies
  dplyr::filter(!checklist_id %in% c("S53637846","G3934955","S53688913")) %>%     # associated lists
  dplyr::left_join(checklists_jun2023 %>% dplyr::distinct(checklist_id), by = "checklist_id") %>%
  dplyr::group_by(locality_id) %>%
  dplyr::mutate(
    matching_percent = sum(!is.na(checklist_id)) / dplyr::n(),
    years_recorded   = dplyr::n_distinct(lubridate::year(observation_date)[lubridate::year(observation_date) %in% c(2019:2023)])
  ) %>%
  dplyr::filter(matching_percent >= 0.5 | years_recorded >= 2) %>%
  dplyr::distinct(locality_id, .keep_all = TRUE) %>%
  dplyr::ungroup()

# Optional: sf view
grja_presence_jun2023_sf <- sf::st_as_sf(grja_presence_jun2023, coords = c("longitude","latitude"), crs = 4326)

# prepare for thinning (robust lon/lat handling)
grja_presence_for_thin <- prep_for_thinning(grja_presence_jun2023, crop_window)

thin_grja_presence_jun2023 <- GeoThinneR::thin_points(
  data       = grja_presence_for_thin,
  lon_col    = "lon",
  lat_col    = "lat",
  group_col  = "common_name",
  method     = "grid",
  raster_obj = raster_grid_coarse,
  trials     = 1,
  all_trials = TRUE
)
thin_grja_presence_jun2023<-largest(thin_grja_presence_jun2023)

# Filter to crop and format for dismo
filtered_grja_presence <- thin_grja_presence_jun2023
grja_p <- filtered_grja_presence %>%
  dplyr::select(lon, lat) %>%
  dplyr::rename(x = lon, y = lat) %>%
  as.data.frame()

# Non-detections for GRJA
nogrja_checklists <- checklists_jun2023 %>%
  dplyr::anti_join(obs_grja_jun2023, by = "locality_id")

# Require at least one checklist in each year 2019–2023
nogrja_checklists_min5 <- nogrja_checklists %>%
  dplyr::filter(lubridate::year(observation_date) %in% 2019:2023) %>%
  dplyr::group_by(locality_id) %>%
  dplyr::filter(dplyr::n_distinct(lubridate::year(observation_date)) == 5) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

nogrja_checklists_min5_sf <- sf::st_as_sf(nogrja_checklists_min5, coords = c("longitude","latitude"), crs = 4326)

filtered_nogrja_checklists_min5 <- .crop_pts(nogrja_checklists_min5, x = "longitude", y = "latitude")
nogrja_p <- filtered_nogrja_checklists_min5 %>%
  dplyr::select(longitude, latitude) %>%
  dplyr::rename(x = longitude, y = latitude) %>%
  as.data.frame()

# ---- MaxEnt (GRJA) ----
grja_maxent_model <- dismo::maxent(cur_rs, grja_p, removeDuplicates = TRUE, nbg = 10000, args = c('jackknife=true'))
print(grja_maxent_model)

e1 <- dismo::evaluate(grja_maxent_model, p = grja_p, a = nogrja_p, x = cur_rs)
plot(e1, 'ROC'); plot(e1, 'TPR'); boxplot(e1); density(e1); threshold(e1)

grja_predict_current <- raster::predict(cur_rs, grja_maxent_model)
grja_predict_current_df <- as.data.frame(grja_predict_current, xy = TRUE)
readr::write_csv(grja_predict_current_df, "grja_predict_current_df.csv")
saveRDS(grja_predict_current_df, "grja_predict_current_df.rds")

grja_p_sf   <- sf::st_as_sf(grja_p,   coords = c("x","y"), crs = 4326)
nogrja_p_sf <- sf::st_as_sf(nogrja_p, coords = c("x","y"), crs = 4326)

# Future projection (GRJA)
grja_predict_future <- raster::predict(fut_rs, grja_maxent_model)
grja_predict_future_df <- as.data.frame(grja_predict_future, xy = TRUE)
readr::write_csv(grja_predict_future_df, "grja_predict_future_df.csv")
saveRDS(grja_predict_future_df, "grja_predict_future_df.rds")

# -----------------------
# 4) BLJA (Blue Jay): eBird → thin → MaxEnt
# -----------------------
# Checklists additions (LA, OK)
smp_la <- "ebd_US-LA_blujay_smp_relJul-2024/ebd_US-LA_blujay_smp_relJul-2024_sampling.txt"
smp_ok <- "ebd_US-OK_blujay_smp_relJul-2024/ebd_US-OK_blujay_smp_relJul-2024_sampling.txt"
checklists_la <- auk::read_sampling(smp_la, unique = TRUE)
checklists_ok <- auk::read_sampling(smp_ok, unique = TRUE)

checklists_blja <- dplyr::bind_rows(checklists, checklists_ok, checklists_la)
checklists_blja_jun2023 <- checklists_blja %>%
  dplyr::filter(observation_date >= as.Date("1900-01-01"), observation_date <= as.Date("2023-05-31")) %>%
  dplyr::arrange(observation_date)

# Observations - BLJA (TX, LA, OK)
ebd_blja_la <- "ebd_US-LA_blujay_smp_relJul-2024/ebd_US-LA_blujay_smp_relJul-2024.txt"
ebd_blja_tx <- "ebd_US-TX_blujay_relJul-2024/ebd_US-TX_blujay_relJul-2024.txt"
ebd_blja_ok <- "ebd_US-OK_blujay_smp_relJul-2024/ebd_US-OK_blujay_smp_relJul-2024.txt"

obs_blja_la <- auk::read_ebd(ebd_blja_la)
obs_blja_tx <- auk::read_ebd(ebd_blja_tx)
obs_blja_ok <- auk::read_ebd(ebd_blja_ok)
obs_blja    <- dplyr::bind_rows(obs_blja_la, obs_blja_tx, obs_blja_ok)

obs_blja <- obs_blja %>%
  dplyr::filter(all_species_reported, dplyr::between(lubridate::year(observation_date), 1900, 2024)) %>%
  dplyr::semi_join(checklists_blja_jun2023, by = "checklist_id")

# Count of BLJA checklists rows (pre-hybrid window)
obs_blja_num <- dplyr::semi_join(obs_blja, checklists_blja_jun2023, by = "checklist_id")
cat("BLJA checklist rows (pre-hybrid window): ", nrow(obs_blja_num), "\n")

# Observations through 2023-05-31 (limited start for thinning perf)
obs_blja_jun2023 <- obs_blja %>%
  dplyr::filter(all_species_reported,
                observation_date >= as.Date("2018-01-01"),
                observation_date <= as.Date("2023-05-31")) %>%
  dplyr::arrange(observation_date)

# Filter vagrancies & retain recurring localities
blja_presence_jun2023 <- obs_blja_jun2023 %>%
  dplyr::left_join(checklists_blja_jun2023 %>% dplyr::distinct(checklist_id), by = "checklist_id") %>%
  dplyr::group_by(locality_id) %>%
  dplyr::filter(!locality_id %in% c("L355719","L1391443","L4204061")) %>%
  dplyr::mutate(
    matching_percent = sum(!is.na(checklist_id)) / dplyr::n(),
    years_recorded   = dplyr::n_distinct(lubridate::year(observation_date)[lubridate::year(observation_date) %in% c(2019:2023)])
  ) %>%
  dplyr::filter(matching_percent >= 0.5 | years_recorded >= 2) %>%
  dplyr::distinct(locality_id, .keep_all = TRUE) %>%
  dplyr::ungroup()

# Prepare lon/lat df and crop before thinning
blja_presence_for_thin <- to_lonlat_df(blja_presence_jun2023)
blja_presence_for_thin <- .crop_pts(blja_presence_for_thin, x = "longitude", y = "latitude")

thin_blja_presence_jun2023 <- GeoThinneR::thin_points(
  data          = blja_presence_for_thin,
  lon_col      = "longitude",
  lat_col       = "latitude",
  group_col     = "common_name",
  method        = "grid",
  raster_obj    = raster_grid_coarse,
  trials        = 1,
  all_trials    = TRUE
)
thin_blja_presence_jun2023<-largest(thin_blja_presence_jun2023)

# Presence (BLJA) for dismo
blja_p <- thin_blja_presence_jun2023 %>%
  dplyr::select(longitude, latitude) %>%
  dplyr::rename(x = longitude, y = latitude) %>%
  as.data.frame()

# Non-detections for BLJA
noblja_checklists <- checklists_jun2023 %>% dplyr::anti_join(obs_blja_jun2023, by = "locality_id")
noblja_checklists_min5 <- noblja_checklists %>%
  dplyr::filter(lubridate::year(observation_date) %in% 2019:2023) %>%
  dplyr::group_by(locality_id) %>%
  dplyr::filter(dplyr::n_distinct(lubridate::year(observation_date)) == 5) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

filtered_noblja_checklists_min5 <- .crop_pts(noblja_checklists_min5, x = "longitude", y = "latitude")
noblja_p <- filtered_noblja_checklists_min5 %>%
  dplyr::select(longitude, latitude) %>%
  dplyr::rename(x = longitude, y = latitude) %>%
  as.data.frame()

# ---- MaxEnt (BLJA) ----
blja_maxent_model <- dismo::maxent(cur_rs, blja_p, removeDuplicates = TRUE, nbg = 10000, args = c('jackknife=true'))
print(blja_maxent_model)

e2 <- dismo::evaluate(blja_maxent_model, p = blja_p, a = noblja_p, x = cur_rs)
plot(e2, 'ROC'); plot(e2, 'TPR'); boxplot(e2); density(e2); threshold(e2)

blja_predict_current <- raster::predict(cur_rs, blja_maxent_model)
blja_predict_current_df <- as.data.frame(blja_predict_current, xy = TRUE)
readr::write_csv(blja_predict_current_df, "blja_predict_current_df.csv")
saveRDS(blja_predict_current_df, "blja_predict_current_df.rds")

# Sanity check and write combined current
raster::compareRaster(blja_predict_current, grja_predict_current, extent=TRUE, rowcol=TRUE, crs=TRUE, res=TRUE, stopiffalse=FALSE)
combined_current_stack <- raster::stack(blja_predict_current, grja_predict_current)
names(combined_current_stack) <- c("blja", "grja")
writeRaster(combined_current_stack, filename="combined_current_stack.tif", format="GTiff", overwrite=TRUE)

# Spatial versions
blja_p_sf   <- sf::st_as_sf(blja_p,   coords = c("x","y"), crs = 4326)
noblja_p_sf <- sf::st_as_sf(noblja_p, coords = c("x","y"), crs = 4326)

# Future projection (BLJA)
blja_predict_future <- raster::predict(fut_rs, blja_maxent_model)
blja_predict_future_df <- as.data.frame(blja_predict_future, xy = TRUE)
readr::write_csv(blja_predict_future_df, "blja_predict_future_df.csv")
saveRDS(blja_predict_future_df, "blja_predict_future_df.rds")

raster::compareRaster(blja_predict_future, grja_predict_future, extent=TRUE, rowcol=TRUE, crs=TRUE, res=TRUE, stopiffalse=FALSE)
combined_future_stack <- raster::stack(blja_predict_future, grja_predict_future)
names(combined_future_stack) <- c("blja", "grja")
writeRaster(combined_future_stack, filename="combined_future_stack.tif", format="GTiff", overwrite=TRUE)

