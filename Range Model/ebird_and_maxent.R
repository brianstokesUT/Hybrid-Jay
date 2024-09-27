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

# Sampling Datasets - GRJA 
smp_mx <- "ebd_MX_grnjay_smp_relJul-2024/ebd_MX_grnjay_smp_relJul-2024_sampling.txt"
smp_tx <- "ebd_US-TX_grnjay_smp_relJul-2024/ebd_US-TX_grnjay_smp_relJul-2024_sampling.txt"
#smp_us <- "ebd_US_grnjay_smp_relJul-2024/ebd_US_grnjay_smp_relJul-2024_sampling.txt"

# read in files, only keeping one checklist in the case of group samping events
checklists_mx <- read_sampling(smp_mx, unique = TRUE)
checklists_tx <- read_sampling(smp_tx, unique = TRUE)
#checklists_us <- read_sampling(smp_us)

# combine checklists
checklists <- bind_rows(checklists_mx, checklists_tx)

# output only checklists which occurred before we caught hybrid
checklists_jun2023 <- checklists %>%
  filter(observation_date >= as.Date("1900-01-01") & observation_date <= as.Date("2023-05-31"))
checklists_jun2023 <- checklists_jun2023 %>% arrange(observation_date)


# Observations Dataset - GRJA
ebd_grja_mx <- "ebd_MX_grnjay_smp_relJul-2024/ebd_MX_grnjay_smp_relJul-2024.txt"
ebd_grja_us <- "ebd_US_grnjay_smp_relJul-2024/ebd_US_grnjay_smp_relJul-2024.txt"

obs_grja_mx <- read_ebd(ebd_grja_mx)
obs_grja_us <- read_ebd(ebd_grja_us)

# combine obs dataframes
obs_grja <- bind_rows(obs_grja_mx, obs_grja_us)

# Filter Observations to ensure dates are normal
obs_grja <- obs_grja |> 
  filter(all_species_reported,
         between(year(observation_date), 1900, 2024))

# remove observations without a matching checklist (will also filter based on the previously applied date filter)
obs_grja <- semi_join(obs_grja, checklists_jun2023, by = "checklist_id")

# output only observations where green jay were observed before we caught hybrid 
obs_grja_jun2023 <- obs_grja |> 
  filter(all_species_reported,
         observation_date >= as.Date("1900-01-01") & observation_date <= as.Date("2023-05-31"))
obs_grja_jun2023 <- obs_grja_jun2023 %>% arrange(observation_date)

# Filter our observational data to remove outliers / vagranices but still keep recent expansion regions and under sampled points
grja_presence_jun2023 <- obs_grja_jun2023 %>%
  # Manually toss out two localities which were highly observed vagrancies - One in downtown Houston and one at a private home which had no other checklists ever.
  filter(!locality_id %in% c("L8812047", "L13398708")) %>%
  # A double check those localities were removed + an additional in downtown Houston
  filter(!checklist_id %in% c("S53637846", "G3934955", "S53688913")) %>%
  # Join with checklists_unique to get matching rows based on checklist_id
  left_join(checklists_jun2023 %>% distinct(checklist_id), by = "checklist_id") %>%
  group_by(locality_id) %>%
  # Calculate the percentage of matching checklist events
  mutate(matching_percent = sum(!is.na(checklist_id)) / n()) %>%
  # Check if recordings occurred in at least 3 years between 2020 and 2023
  mutate(years_recorded = n_distinct(year(observation_date)[year(observation_date) %in% c(2019, 2020, 2021, 2022, 2023)])) %>%
  # Apply the two conditions: matching percent >= 50% OR at least 3 years recorded between 2019-2023
  # This means that GRJA must be seen in at least half of all checklists - or there checklists within at least 1 year in the previous 4.5 (2019 - May 2023)
  # Ultimately this step is pretty unnecessary based on the next thinnign steps we follow - but is probably good practice for a species with a biased observational process like GRJA...
  filter(matching_percent >= 0.5 | years_recorded >= 1) %>%
  # Keep a single checklist per locality
  distinct(locality_id, .keep_all = TRUE) %>%
  ungroup()

# Convert to sf object
grja_presence_jun2023_sf <- grja_presence_jun2023 %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# Thin out the datapoints to correct to spatial bias in observational process (primarily an issue between TX and Mexico) 
# Previously tested without this step and causes strange maxent model behvior related to US/Mexico border observational differences
thin_raster <- terra::rast(xmin = -105, xmax = -93, ymin = 20, ymax = 35, res = 0.25) #creates raster "grid" to base thinning on - saves one point per grid block
thin_grja_presence_jun2023 <- thin_points(
  grja_presence_jun2023,
  long_col = "longitude",
  lat_col = "latitude",
  group_col = "common_name",
  method = "grid",
  raster_obj = thin_raster,
  trials = 1,
  all_trials = TRUE
)

# Make dataset of localities which have never had GRJA sightings
nogrja_checklists <- checklists_jun2023 %>%
  anti_join(obs_grja_jun2023, by = "locality_id")

#want to throw out sampling locations with poor observation history (must have minimum 5 different years of sampling events)
nogrja_checklists_min5 <- nogrja_checklists %>%
  group_by(locality_id) %>%
  filter(n_distinct(year(observation_date)) >= 5) %>%
  slice(1)  # Select a single row per locality_id (e.g., first row)

#make sf frame for ggplot
nogrja_checklists_min5_sf <- nogrja_checklists_min5 %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

#### GRJA Maxent Modeling ####

# Remove grja points not within our environmental "predictors" rasters (Normal_1991_2020_stack)
filtered_grja_presence <- thin_grja_presence_jun2023[[1]] %>%
  filter(longitude >= xmin(crop_extent) & longitude <= xmax(crop_extent) &
           latitude >= ymin(crop_extent) & latitude <= ymax(crop_extent))
# Reformat presence dataframe for dismo specifications
grja_p <- filtered_grja_presence %>%
  dplyr::select(longitude, latitude) %>%
  dplyr::rename(x = longitude, y = latitude) %>%
  as.data.frame()

# Remove nogrja points not within out rasters
filtered_nogrja_checklists_min5 <- nogrja_checklists_min5 %>%
  filter(longitude >= xmin(crop_extent) & longitude <= xmax(crop_extent) &
           latitude >= ymin(crop_extent) & latitude <= ymax(crop_extent))
# Reformat nogrja dataframe for dismo specifications
nogrja_p <- filtered_nogrja_checklists_min5 %>%
  ungroup() %>%  # Ensure no grouping variables are retained
  dplyr::select(longitude, latitude) %>%
  dplyr::rename(x = longitude, y = latitude) %>%
  as.data.frame()

# Create maxent model based on "current" data
grja_maxent_model<-maxent(Normal_1991_2020_stack, grja_p, removeDuplicates=TRUE, nbg = 3350)
grja_maxent_model

# Optionally, evaluate the model vs the nogrja points
#e1<-evaluate(grja_maxent_model, p=grja_p, a=nogrja_p, x=Normal_1991_2020_stack)
#plot(e1, 'ROC')
#plot(e1, 'TPR')
#boxplot(e1)
#density(e1)

# Predict current spp distribution
grja_predict_current<-predict(grja_maxent_model, Normal_1991_2020_stack)

# Convert raster to a format compatible with ggplot and save as csv
grja_predict_current_df <- as.data.frame(grja_predict_current, xy = TRUE)
write.csv(grja_predict_current_df, "grja_predict_current_df.csv", row.names = FALSE)

#make training data a spatial object
grja_p_sf <- grja_p %>%
  st_as_sf(coords = c("x", "y"), crs = 4326)
#make validation data a spatial object
nogrja_p_sf <- nogrja_p %>%
  st_as_sf(coords = c("x", "y"), crs = 4326)

# Future Range Projection
# predict current spp distribution
grja_predict_future<-predict(grja_maxent_model, future_stack)

# Convert raster to a format compatible with ggplot and save as csv
grja_predict_future_df <- as.data.frame(grja_predict_future, xy = TRUE)
write.csv(grja_predict_future_df, "grja_predict_future_df.csv", row.names = FALSE)


#### BLJA eBird Data ####
# update checklists - BLJA (checklists are the same regardless of species in the eBird dataset)
smp_la <- "ebd_US-LA_blujay_smp_relJul-2024/ebd_US-LA_blujay_smp_relJul-2024_sampling.txt"
checklists_la <- read_sampling(smp_la, unique = TRUE)
smp_ok <- "ebd_US-OK_blujay_smp_relJul-2024/ebd_US-OK_blujay_smp_relJul-2024_sampling.txt"
checklists_ok <- read_sampling(smp_ok, unique = TRUE)

# add these to previous checklists of MX/TX we made for GRJA
checklists_blja <- bind_rows(checklists, checklists_ok, checklists_la)

# output only checklists which occurred before we caught hybrid
checklists_blja_jun2023 <- checklists_blja %>%
  filter(observation_date >= as.Date("1900-01-01") & observation_date <= as.Date("2023-05-31"))
checklists_blja_jun2023 <- checklists_blja_jun2023 %>% arrange(observation_date)

# Observations Dataset - BLJA (just usinng Texas, Lousiana, Oklahoma so that the data remains a reasonable size)
ebd_blja_la <- "ebd_US-LA_blujay_smp_relJul-2024/ebd_US-LA_blujay_smp_relJul-2024.txt"
ebd_blja_tx <- "ebd_US-TX_blujay_relJul-2024/ebd_US-TX_blujay_relJul-2024.txt"
ebd_blja_ok <- "ebd_US-OK_blujay_smp_relJul-2024/ebd_US-OK_blujay_smp_relJul-2024.txt"

obs_blja_la <- read_ebd(ebd_blja_la)
obs_blja_tx <- read_ebd(ebd_blja_tx)
obs_blja_ok <- read_ebd(ebd_blja_ok)

# combine obs dataframes
obs_blja <- bind_rows(obs_blja_la, obs_blja_tx, obs_blja_ok)

# Filter Observations to ensure dates are normal
obs_blja <- obs_blja |> 
  filter(all_species_reported,
         between(year(observation_date), 1900, 2024))

# remove observations without a matching checklist (should also filter out dates)
obs_blja <- semi_join(obs_blja, checklists_blja_jun2023, by = "checklist_id")

# output only observations where blue jay were observed before we caught hybrid (changing starting date from GRJA bc I'm having vector size issues in the thinning step)
obs_blja_jun2023 <- obs_blja |> 
  filter(all_species_reported,
         observation_date >= as.Date("2018-01-01") & observation_date <= as.Date("2023-05-31"))
obs_blja_jun2023 <- obs_blja_jun2023 %>% arrange(observation_date)

# Filter our observational data to remove outliers / vagranices but still keep recent expansion regions and under sampled points
blja_presence_jun2023 <- obs_blja_jun2023 %>%
  left_join(checklists_blja_jun2023 %>% distinct(checklist_id), by = "checklist_id") %>%
  group_by(locality_id) %>%
  mutate(matching_percent = sum(!is.na(checklist_id)) / n()) %>%
  mutate(years_recorded = n_distinct(year(observation_date)[year(observation_date) %in% c(2019, 2020, 2021, 2022, 2023)])) %>%
  # doing 3 years instead - BLJA range edge has been slihtly more stable then GRJA
  # Similar to GRJA - makes minimal difference in model behavior for majority of range, but does effectivley remove the influence of extreme vagrancies
  filter(matching_percent >= 0.5 | years_recorded >= 3) %>%
  distinct(locality_id, .keep_all = TRUE) %>%
  ungroup()

# Crop these points by Lat and Long (Thinning step was too memory intensive for my machine)
crop_blja_presence_jun2023 <- blja_presence_jun2023 %>%
  filter(longitude >= xmin(crop_extent) & longitude <= xmax(crop_extent) &
           latitude >= ymin(crop_extent) & latitude <= ymax(crop_extent))

# Thin out the datapoints to match the sampling of GRJA dataset
thin_raster <- terra::rast(xmin = -105, xmax = -93, ymin = 20, ymax = 35, res = 0.25)
thin_blja_presence_jun2023 <- thin_points(
  crop_blja_presence_jun2023,
  long_col = "longitude",
  lat_col = "latitude",
  group_col = "common_name",
  method = "grid",
  raster_obj = thin_raster,
  trials = 1,
  all_trials = TRUE
)

# Reformat presence dataframe for dismo specifications
blja_p <- thin_blja_presence_jun2023[[1]] %>%
  dplyr::select(longitude, latitude) %>%
  dplyr::rename(x = longitude, y = latitude) %>%
  as.data.frame()


#### BLJA MAXENT ####
# Create maxent model based on "current" data
blja_maxent_model<-maxent(Normal_1991_2020_stack, blja_p, removeDuplicates=TRUE, nbg = 6000)
blja_maxent_model

# predict current spp distribution
blja_predict_current<-predict(blja_maxent_model, Normal_1991_2020_stack)

# Convert raster to a format compatible with ggplot
blja_predict_current_df <- as.data.frame(blja_predict_current, xy = TRUE)
write.csv(blja_predict_current_df, "blja_predict_current_df.csv", row.names = FALSE)

#make training data a spatial object
blja_p_sf <- blja_p %>%
  st_as_sf(coords = c("x", "y"), crs = 4326)
#make validation data a spatial object
noblja_p_sf <- noblja_p %>%
  st_as_sf(coords = c("x", "y"), crs = 4326)

# Future Range Projection
# predict current spp distribution
blja_predict_future<-predict(blja_maxent_model, future_stack)

# Convert raster to a format compatible with ggplot
blja_predict_future_df <- as.data.frame(blja_predict_future, xy = TRUE)
write.csv(blja_predict_future_df, "blja_predict_future_df.csv", row.names = FALSE)
