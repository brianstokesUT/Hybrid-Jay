wd<-("~/PATH")
setwd(wd)

set.seed(223)
options(java.parameters = "-Xmx8000m")


library(auk)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(lubridate)
library(readr)
library(sf)
library("geosphere")
library("terra")
library("rJava")
library("dismo")

#### GRJA eBird Data ####

# OPTIONAL: Read in Env Rasters
# Normal_1991_2020_stack <- rast("Normal_1991_2020_stack.tif")
# future_stack            <- rast("future_stack.tif")

## ------------------------------------------------------------------
## Sampling datasets – GRJA
## ------------------------------------------------------------------

smp_mx <- "ebd_MX_grnjay_smp_relJul-2024/ebd_MX_grnjay_smp_relJul-2024_sampling.txt"
smp_tx <- "ebd_US-TX_grnjay_smp_relJul-2024/ebd_US-TX_grnjay_smp_relJul-2024_sampling.txt"
# smp_us <- "ebd_US_grnjay_smp_relJul-2024/ebd_US_grnjay_smp_relJul-2024_sampling.txt"

# Read in files, only keeping one checklist in the case of groups
checklists_mx <- read_sampling(smp_mx, unique = TRUE)
checklists_tx <- read_sampling(smp_tx, unique = TRUE)
# checklists_us <- read_sampling(smp_us)

# Combine checklists
checklists <- bind_rows(checklists_mx, checklists_tx)

# Keep only checklists before hybrid capture and sort by date
checklists_jun2023 <- checklists %>%
  filter(
    observation_date >= as.Date("1900-01-01"),
    observation_date <= as.Date("2023-05-31")
  ) %>%
  arrange(observation_date)

## ------------------------------------------------------------------
## Observations – GRJA
## ------------------------------------------------------------------

ebd_grja_mx <- "ebd_MX_grnjay_smp_relJul-2024/ebd_MX_grnjay_smp_relJul-2024.txt"
ebd_grja_us <- "ebd_US_grnjay_smp_relJul-2024/ebd_US_grnjay_smp_relJul-2024.txt"

obs_grja_mx <- read_ebd(ebd_grja_mx)
obs_grja_us <- read_ebd(ebd_grja_us)

# Combine obs dataframes
obs_grja <- bind_rows(obs_grja_mx, obs_grja_us)

# Filter observations to ensure dates are normal
obs_grja <- obs_grja %>%
  filter(
    all_species_reported,
    between(year(observation_date), 1900, 2024)
  )

# Remove observations without a matching checklist (also filters dates)
obs_grja <- semi_join(obs_grja, checklists_jun2023, by = "checklist_id")

# Keep only observations before hybrid capture and sort by date
obs_grja_jun2023 <- obs_grja %>%
  filter(
    all_species_reported,
    observation_date >= as.Date("1900-01-01"),
    observation_date <= as.Date("2023-05-31")
  ) %>%
  arrange(observation_date)

## ------------------------------------------------------------------
## GRJA presence data: filter vagrants and enforce sampling history
## ------------------------------------------------------------------

grja_presence_jun2023 <- obs_grja_jun2023 %>%
  # Manually toss out all checklists for two localities of highly observed vagrancies
  filter(!locality_id %in% c("L8812047", "L13398708")) %>%
  # Manually toss out other associated checklists (same general area)
  filter(!checklist_id %in% c("S53637846", "G3934955", "S53688913")) %>%
  # Join with checklists to get matching rows based on checklist_id
  left_join(checklists_jun2023 %>% distinct(checklist_id), by = "checklist_id") %>%
  group_by(locality_id) %>%
  # Percentage of matching checklist events for species across all time
  mutate(matching_percent = sum(!is.na(checklist_id)) / n()) %>%
  # Check if recordings occurred in at least 2 years between 2019–2023
  mutate(
    years_recorded = n_distinct(
      year(observation_date)[year(observation_date) %in% c(2019, 2020, 2021, 2022, 2023)]
    )
  ) %>%
  # Condition: matching_percent >= 50% OR at least 2 years recorded (2019–2023)
  filter(matching_percent >= 0.5 | years_recorded >= 2) %>%
  distinct(locality_id, .keep_all = TRUE) %>%
  ungroup()

# Convert to sf object
grja_presence_jun2023_sf <- grja_presence_jun2023 %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

## ------------------------------------------------------------------
## Thinning GRJA points (10 km grid)
## ------------------------------------------------------------------

raster_grid_10km <- aggregate(Normal_1991_2020_stack, fact = 10)

thin_grja_presence_jun2023 <- thin_points(
  grja_presence_jun2023,
  lon_col   = "longitude",
  lat_col   = "latitude",
  group_col = "common_name",
  method    = "grid",
  raster_obj = raster_grid_10km,
  trials    = 1,
  all_trials = TRUE
)

## ------------------------------------------------------------------
## GRJA non-detections (localities with no GRJA)
## ------------------------------------------------------------------

nogrja_checklists <- checklists_jun2023 %>%
  anti_join(obs_grja_jun2023, by = "locality_id")

# Throw out sampling locations with poor observation history
nogrja_checklists_min5 <- nogrja_checklists %>%
  filter(year(observation_date) %in% 2019:2023) %>%  # Keep only target years
  group_by(locality_id) %>%
  filter(n_distinct(year(observation_date)) == 5) %>%  # At least one checklist in each year
  slice(1) %>%                                         # One row per locality
  ungroup()

# sf for plotting
nogrja_checklists_min5_sf <- nogrja_checklists_min5 %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

#### GRJA Maxent Modeling ####

# Crop extent from env stack
crop_extent <- ext(Normal_1991_2020_stack)

# Get thinned data for trial 1 (changing trial will slightly affect results)
thin_grja_presence_jun2023_trial1 <- get_trial(thin_grja_presence_jun2023, trial = 1)

# Remove GRJA points outside environmental raster extent
filtered_grja_presence <- thin_grja_presence_jun2023_trial1 %>%
  filter(
    longitude >= xmin(crop_extent),
    longitude <= xmax(crop_extent),
    latitude  >= ymin(crop_extent),
    latitude  <= ymax(crop_extent)
  )

# Reformat presence dataframe for dismo specifications
grja_p <- filtered_grja_presence %>%
  dplyr::select(longitude, latitude) %>%
  dplyr::rename(x = longitude, y = latitude) %>%
  as.data.frame()

# Remove non-GRJA points outside environmental raster extent
filtered_nogrja_checklists_min5 <- nogrja_checklists_min5 %>%
  filter(
    longitude >= xmin(crop_extent),
    longitude <= xmax(crop_extent),
    latitude  >= ymin(crop_extent),
    latitude  <= ymax(crop_extent)
  )

# Reformat non-detection dataframe for dismo specifications
nogrja_p <- filtered_nogrja_checklists_min5 %>%
  ungroup() %>%
  dplyr::select(longitude, latitude) %>%
  dplyr::rename(x = longitude, y = latitude) %>%
  as.data.frame()

# Create maxent model based on "current" data using all bioclim layers
env <- raster::stack(Normal_1991_2020_stack)

grja_maxent_model <- maxent(
  env,
  grja_p,
  removeDuplicates = TRUE,
  nbg             = 10000,
  args            = c("jackknife=true")
)
grja_maxent_model

# Evaluate the model vs the non-GRJA points
e1 <- evaluate(
  grja_maxent_model,
  p = grja_p,
  a = nogrja_p,
  x = Normal_1991_2020_stack
)

plot(e1, "ROC")
plot(e1, "TPR")
boxplot(e1)
density(e1)
threshold(e1)

# Predict current species distribution
grja_predict_current <- predict(grja_maxent_model, env)

# Convert raster to df for ggplot and save
grja_predict_current_df <- as.data.frame(grja_predict_current, xy = TRUE)
write.csv(grja_predict_current_df, "grja_predict_current_df.csv", row.names = FALSE)
saveRDS(grja_predict_current_df, "grja_predict_current_df.rds")

# Training data as sf
grja_p_sf <- grja_p %>%
  st_as_sf(coords = c("x", "y"), crs = 4326)

# Validation data as sf
nogrja_p_sf <- nogrja_p %>%
  st_as_sf(coords = c("x", "y"), crs = 4326)

## ------------------------------------------------------------------
## GRJA Future Range Projection
## ------------------------------------------------------------------

fut_env <- raster::stack(future_stack)

# Coerce future names to match current
names(fut_env) <- names(env)

# Predict future species distribution
grja_predict_future <- predict(grja_maxent_model, fut_env)

# Convert raster to df and save
grja_predict_future_df <- as.data.frame(grja_predict_future, xy = TRUE)
write.csv(grja_predict_future_df, "grja_predict_future_df.csv", row.names = FALSE)
saveRDS(grja_predict_future_df, "grja_predict_future_df.rds")

#### BLJA eBird Data ####

## ------------------------------------------------------------------
## Sampling datasets – BLJA
## ------------------------------------------------------------------

# Update checklists – BLJA (checklists are the same regardless of species)
smp_la <- "ebd_US-LA_blujay_smp_relJul-2024/ebd_US-LA_blujay_smp_relJul-2024_sampling.txt"
smp_ok <- "ebd_US-OK_blujay_smp_relJul-2024/ebd_US-OK_blujay_smp_relJul-2024_sampling.txt"

checklists_la <- read_sampling(smp_la, unique = TRUE)
checklists_ok <- read_sampling(smp_ok, unique = TRUE)

# Add to previous MX/TX checklists
checklists_blja <- bind_rows(checklists, checklists_ok, checklists_la)

# Keep only checklists before hybrid capture and sort by date
checklists_blja_jun2023 <- checklists_blja %>%
  filter(
    observation_date >= as.Date("1900-01-01"),
    observation_date <= as.Date("2023-05-31")
  ) %>%
  arrange(observation_date)

## ------------------------------------------------------------------
## Observations – BLJA
## ------------------------------------------------------------------

# Only using TX, LA, OK to keep data manageable
ebd_blja_la <- "ebd_US-LA_blujay_smp_relJul-2024/ebd_US-LA_blujay_smp_relJul-2024.txt"
ebd_blja_tx <- "ebd_US-TX_blujay_relJul-2024/ebd_US-TX_blujay_relJul-2024.txt"
ebd_blja_ok <- "ebd_US-OK_blujay_smp_relJul-2024/ebd_US-OK_blujay_smp_relJul-2024.txt"

obs_blja_la <- read_ebd(ebd_blja_la)
obs_blja_tx <- read_ebd(ebd_blja_tx)
obs_blja_ok <- read_ebd(ebd_blja_ok)

# Combine obs dataframes
obs_blja <- bind_rows(obs_blja_la, obs_blja_tx, obs_blja_ok)

# Filter observations to ensure dates are normal
obs_blja <- obs_blja %>%
  filter(
    all_species_reported,
    between(year(observation_date), 1900, 2024)
  )

# Remove observations without a matching checklist
obs_blja <- semi_join(obs_blja, checklists_blja_jun2023, by = "checklist_id")

# Next two lines pull number used in manuscript for count of BLJA checklists in TX
obs_blja_num <- semi_join(obs_blja, checklists_blja_jun2023, by = "checklist_id")
nrow(obs_blja_num)

# Keep BLJA observations before hybrid capture (starting in 2018 here)
obs_blja_jun2023 <- obs_blja %>%
  filter(
    all_species_reported,
    observation_date >= as.Date("2018-01-01"),
    observation_date <= as.Date("2023-05-31")
  ) %>%
  arrange(observation_date)

## ------------------------------------------------------------------
## BLJA presence data: filter vagrants and enforce sampling history
## ------------------------------------------------------------------

blja_presence_jun2023 <- obs_blja_jun2023 %>%
  # Join with checklists to get matching rows based on checklist_id
  left_join(checklists_blja_jun2023 %>% distinct(checklist_id), by = "checklist_id") %>%
  group_by(locality_id) %>%
  # Remove known vagrancies
  filter(!locality_id %in% c("L355719", "L1391443", "L4204061")) %>%
  # Percentage of matching checklist events
  mutate(matching_percent = sum(!is.na(checklist_id)) / n()) %>%
  # Check if recordings occurred in at least 2 years between 2019–2023
  mutate(
    years_recorded = n_distinct(
      year(observation_date)[year(observation_date) %in% c(2019, 2020, 2021, 2022, 2023)]
    )
  ) %>%
  # Condition: matching_percent >= 50% OR at least 2 years recorded (2019–2023)
  filter(matching_percent >= 0.5 | years_recorded >= 2) %>%
  distinct(locality_id, .keep_all = TRUE) %>%
  ungroup()

# Crop BLJA presence to environmental raster extent (reduces memory load)
crop_blja_presence_jun2023 <- blja_presence_jun2023 %>%
  filter(
    longitude >= xmin(crop_extent),
    longitude <= xmax(crop_extent),
    latitude  >= ymin(crop_extent),
    latitude  <= ymax(crop_extent)
  )

## ------------------------------------------------------------------
## Thinning BLJA points (10 km grid)
## ------------------------------------------------------------------

raster_grid_10km <- aggregate(Normal_1991_2020_stack, fact = 10)

thin_blja_presence_jun2023 <- thin_points(
  crop_blja_presence_jun2023,
  lon_col   = "longitude",
  lat_col   = "latitude",
  group_col = "common_name",
  method    = "grid",
  raster_obj = raster_grid_10km,
  trials    = 1,
  all_trials = TRUE
)

# Get thinned data for trial 1
thin_blja_presence_jun2023_trial1 <- get_trial(thin_blja_presence_jun2023, trial = 1)

# Remove BLJA points outside environmental raster extent
filtered_blja_presence <- thin_blja_presence_jun2023_trial1 %>%
  filter(
    longitude >= xmin(crop_extent),
    longitude <= xmax(crop_extent),
    latitude  >= ymin(crop_extent),
    latitude  <= ymax(crop_extent)
  )

# Reformat presence dataframe for dismo specifications
blja_p <- filtered_blja_presence %>%
  dplyr::select(longitude, latitude) %>%
  dplyr::rename(x = longitude, y = latitude) %>%
  as.data.frame()

## ------------------------------------------------------------------
## BLJA non-detections
## ------------------------------------------------------------------

noblja_checklists <- checklists_jun2023 %>%
  anti_join(obs_blja_jun2023, by = "locality_id")

# Throw out sampling locations with poor observation history
noblja_checklists_min5 <- noblja_checklists %>%
  filter(year(observation_date) %in% 2019:2023) %>%
  group_by(locality_id) %>%
  filter(n_distinct(year(observation_date)) == 5) %>%  # At least one checklist in each year
  slice(1) %>%                                         # One row per locality
  ungroup()

# Remove non-BLJA points outside environmental raster extent
filtered_noblja_checklists_min5 <- noblja_checklists_min5 %>%
  filter(
    longitude >= xmin(crop_extent),
    longitude <= xmax(crop_extent),
    latitude  >= ymin(crop_extent),
    latitude  <= ymax(crop_extent)
  )

# Reformat non-detection dataframe for dismo specifications
noblja_p <- filtered_noblja_checklists_min5 %>%
  ungroup() %>%
  dplyr::select(longitude, latitude) %>%
  dplyr::rename(x = longitude, y = latitude) %>%
  as.data.frame()

#### BLJA Maxent Modeling ####

blja_maxent_model <- maxent(
  env,
  blja_p,
  removeDuplicates = TRUE,
  nbg             = 10000,
  args            = c("jackknife=true")
)
blja_maxent_model

# Evaluate the model vs the non-BLJA points
e2 <- evaluate(
  blja_maxent_model,
  p = blja_p,
  a = noblja_p,
  x = Normal_1991_2020_stack
)

plot(e2, "ROC")
plot(e2, "TPR")
boxplot(e2)
density(e2)
threshold(e2)

# Predict current BLJA distribution
blja_predict_current <- predict(blja_maxent_model, env)

# Convert raster to df and save
blja_predict_current_df <- as.data.frame(blja_predict_current, xy = TRUE)
write.csv(blja_predict_current_df, "blja_predict_current_df.csv", row.names = FALSE, digits = 15)
saveRDS(blja_predict_current_df, "blja_predict_current_df.rds")

# Compare GRJA and BLJA current predictions
compareRaster(
  blja_predict_current,
  grja_predict_current,
  extent      = TRUE,
  rowcol      = TRUE,
  crs         = TRUE,
  res         = TRUE,
  stopiffalse = FALSE
)

combined_current_stack <- stack(blja_predict_current, grja_predict_current)
names(combined_current_stack) <- c("blja", "grja")

writeRaster(
  combined_current_stack,
  filename = "combined_current_stack.tif",
  format   = "GTiff",
  overwrite = TRUE
)

# Training data as sf
blja_p_sf <- blja_p %>%
  st_as_sf(coords = c("x", "y"), crs = 4326)

# Validation data as sf
noblja_p_sf <- noblja_p %>%
  st_as_sf(coords = c("x", "y"), crs = 4326)

## ------------------------------------------------------------------
## BLJA Future Range Projection
## ------------------------------------------------------------------

blja_predict_future <- predict(blja_maxent_model, fut_env)

blja_predict_future_df <- as.data.frame(blja_predict_future, xy = TRUE)
write.csv(blja_predict_future_df, "blja_predict_future_df.csv", row.names = FALSE, digits = 15)
saveRDS(blja_predict_future_df, "blja_predict_future_df.rds")

# Compare GRJA and BLJA future predictions
compareRaster(
  blja_predict_future,
  grja_predict_future,
  extent      = TRUE,
  rowcol      = TRUE,
  crs         = TRUE,
  res         = TRUE,
  stopiffalse = FALSE
)

combined_future_stack <- stack(blja_predict_future, grja_predict_future)
names(combined_future_stack) <- c("blja", "grja")

writeRaster(
  combined_future_stack,
  filename = "combined_future_stack.tif",
  format   = "GTiff",
  overwrite = TRUE
)
