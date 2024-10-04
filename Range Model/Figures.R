wd<-("~/PATH")
setwd(wd)

set.seed(223)

library("dplyr")
library("ggplot2")
library("GeoThinneR")
library("tidyterra")
library("rnaturalearth")
library("rnaturalearthdata")
library("ggpattern")
library("ggnewscale")
library("ggpattern")
library("geosphere")
library("terra")
library("smoothr")
library("devtools")
devtools::install_github("ropensci/rnaturalearthhires")



# Create new grja tibble with one row per unique value of locality_id
obs_grja_jun2023_loc <- grja_presence_jun2023 %>%
  distinct(locality_id, .keep_all = TRUE)
# Create new blja tibble with one row per unique value of locality_id
obs_blja_jun2023_loc <- blja_presence_jun2023 %>%
  distinct(locality_id, .keep_all = TRUE)

# Create overlap (Both species in a single checklist) tibble for plotting - using full dataset
bothspp <- obs_grja %>%
  semi_join(obs_blja, by = "checklist_id")


# Remove extreme vagrancies (A few in from a single GRJA in downtown Houston and one BLJA in LRGV)
bothspp_filter <- bothspp %>%
  # Remove funky vagrancies based on locality_id
  filter(!locality_id %in% c("L8812047", "L13398708")) %>%
  # Double check removal with checklist_id
  filter(!checklist_id %in% c("S53637846", "G3934955", "S53688913", "S43173964")) %>%
  #filter up to end of 2023 and remove data before 1999
  #Note a single observation did occur in 1984 of dubious accuracy
  filter(observation_date >= as.Date("1900-01-01") & observation_date <= as.Date("2023-5-31"))

# Ouput of next line is used in manuscript
nrow(bothspp_filter)

# Thin to single point per locality for plotting
bothspp_loc <- bothspp_filter %>%
  distinct(locality_id, .keep_all = TRUE)

# used in manuscript
nrow(bothspp_loc)

# Create a new column with just the year extracted from observation_date
bothspp_filter <- bothspp_filter %>%
  #One occurance happened in 1984 which ruins the plot - make a note in the figure text
  filter(observation_date >= as.Date("1997-01-01") & observation_date <= as.Date("2023-5-31")) %>%
  mutate(year = as.numeric(format(as.Date(observation_date), "%Y")))
# Count the number of rows per year
yearly_counts <- bothspp_filter %>%
  group_by(year) %>%
  summarise(count = n())

# Count the number of unique locality_ids per year
unique_locality_counts <- bothspp_filter %>%
  group_by(year) %>%
  summarise(unique_localities = n_distinct(locality_id))

# Merging the two summaries into one data frame
combined_counts <- yearly_counts %>%
  left_join(unique_locality_counts, by = "year")

### FIGURE 2B ###
# Plotting the total number of rows per year and overlaying the number of unique locality_ids
fig2b<-ggplot(combined_counts, aes(x = year)) +
  # First bar for total number of observations
  geom_bar(aes(y = count, fill = "Total Observations"), stat = "identity") +
  # Overlay bar for the number of unique locality_ids, with transparency
  geom_bar(aes(y = unique_localities, fill = "Observations at New Locations"), stat = "identity") +
  ggtitle("Green Jay & Blue Jay Co-occurrences") +
  xlab("Year") +
  ylab("Number of Observations") +
  # Custom Fill Legend
  scale_fill_manual(name = "Legend",
                    values = c(
                      "Total Observations" = "grey65",
                      "Observations at New Locations" = "grey15"),
                    breaks = c("Total Observations", "Observations at New Locations")) +
  scale_x_continuous(breaks=seq(1997,2023,2)) +
  theme_classic() +
  theme(
    legend.title=element_blank(),
    aspect.ratio=1,
    plot.title = element_text(size=20, face = "bold"),
    axis.title = element_text(size=16),
    axis.text.x = element_text(size=12, angle=45, vjust=0.5),
    axis.text.y = element_text(size=12),
    legend.position = c(0.30,0.90),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black"),
    legend.text=element_text(size=12))
fig2b

#load dataframes made in "hybird_maxent.R"
blja_predict_current_df <- read.csv("blja_predict_current_df.csv")
blja_predict_current_df %>% drop_na()
blja_predict_current_df <- mutate_all(blja_predict_current_df, function(x) as.numeric(as.character(x)))

grja_predict_current_df <- read.csv("grja_predict_current_df.csv")
grja_predict_current_df %>% drop_na()
grja_predict_current_df <- mutate_all(grja_predict_current_df, function(x) as.numeric(as.character(x)))

blja_predict_future_df <- read.csv("blja_predict_future_df.csv")
blja_predict_future_df %>% drop_na()
blja_predict_future_df <- mutate_all(blja_predict_future_df, function(x) as.numeric(as.character(x)))

grja_predict_future_df <- read.csv("grja_predict_future_df.csv")
grja_predict_future_df %>% drop_na()
grja_predict_future_df <- mutate_all(grja_predict_future_df, function(x) as.numeric(as.character(x)))

# Filter both current dataframes for rows where the layer value is greater than 0.5
blja_curr_min <- blja_predict_current_df %>% filter(layer > 0.5)
grja_curr_min <- grja_predict_current_df %>% filter(layer > 0.5)
# Merge the two dataframes on the x and y columns
current_overlap_df <- inner_join(blja_curr_min, grja_curr_min, by = c("x", "y"), suffix = c("_blja", "_grja"))
# Select relevant columns for the new dataframe
current_overlap_df <- current_overlap_df %>%
  dplyr::select(x, y, layer_blja, layer_grja)
#convert to polygon for ggpattern fill 
current_overlap_spatraster<-rast(current_overlap_df)
crs(current_overlap_spatraster)<-"EPSG:4326"
current_overlap_polygon<-terra::as.polygons(current_overlap_spatraster)
plot(current_overlap_polygon)
current_overlap_sf <- st_as_sf(current_overlap_polygon)


#now output areas where the current ranges do not overlap
blja_no_overlap <- anti_join(blja_predict_current_df, current_overlap_df, by = c("x", "y"))
blja_no_overlap <- blja_no_overlap[blja_no_overlap$layer > 0.5, ]
grja_no_overlap <- anti_join(grja_predict_current_df, current_overlap_df, by = c("x", "y"))
grja_no_overlap <- grja_no_overlap[grja_no_overlap$layer > 0.5, ]

# Threshold the PREDICTION/FUTURE raster to keep only values > 0.25 (ie binary of areas where occupancy probability is greater than 0.25 or not)
# We use a lower threshold than in "current" bc this is predictive and poorly integrates spatiality
# Aggregate to help remove some small outlier patches and have smaller data
blja_predict_future_spatrast<-rast(blja_predict_future_df)
crs(blja_predict_future_spatrast)<-"EPSG:4326"
blja_future_range <- terra::aggregate(blja_predict_future_spatrast, fact=2, fun="max")
blja_future_range[blja_future_range <= 0.25] <- NA
blja_future_range[blja_future_range >= 0.25] <- 1

# Convert the thresholded raster to polygons
blja_future_range_polygon<-terra::as.polygons(blja_future_range)
plot(blja_future_range_polygon)
# use "smoothr" to remove jagged edges for plotting
blja_future_range_polygon_smooth<-smooth(blja_future_range_polygon, method = "ksmooth")
plot(blja_future_range_polygon_smooth)



# Perform same steps for GRJA
grja_predict_future_spatrast<-rast(grja_predict_future_df)
crs(grja_predict_future_spatrast)<-"EPSG:4326"
grja_future_range <- terra::aggregate(grja_predict_future_spatrast, fact=2, fun="max")
grja_future_range[grja_future_range <= 0.25] <- NA
grja_future_range[grja_future_range >= 0.25] <- 1

# Convert the thresholded raster to polygons
grja_future_range_polygon<-terra::as.polygons(grja_future_range)
plot(grja_future_range_polygon)
# use "smoothr" to remove jagged edges for plotting
grja_future_range_polygon_smooth<-smooth(grja_future_range_polygon, method = "ksmooth")
plot(grja_future_range_polygon_smooth)


# Download int'l boundaries using rnatural earth
regional_map_countries<-ne_countries(scale = 10, continent = "North America")
# State boundaries
US_states<-ne_states(country = "united states of america")
MX_states<-ne_states(country = "mexico")

### FIGURE 2A ###
fig2a<- ggplot(data = regional_map_countries) + geom_sf(fill = "white") + 
  geom_sf(data = US_states, color = "black", fill = "transparent", lwd = 0.05) + 
  geom_sf(data = MX_states, color = "black", fill = "transparent", lwd = 0.05) +
  geom_rect(aes(xmin = -101, xmax = -94, ymin = 25.5, ymax = 31), color = "red", lwd=1, fill = NA) +
  xlim(-105, -80) + ylim(15, 35) +
  theme_test() + 
  ggtitle("Region of Study") +
  xlab("Latitude") +
  ylab("Longitude") +
  theme(panel.background = element_rect(fill = "lightblue"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.title = element_text(size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        aspect.ratio = 1,
        plot.title = element_text(size=20, face = "bold"))

fig2a










blja_no_overlap_spatraster<-rast(blja_no_overlap)
crs(current_overlap_spatraster)<-"EPSG:4326"
blja_no_overlap_polygon<-terra::as.polygons(blja_no_overlap_spatraster)
plot(blja_no_overlap_polygon)
blja_no_overlap_sf <- st_as_sf(blja_no_overlap_polygon)
st_crs(blja_no_overlap_sf)<-"EPSG:4326"

grja_no_overlap_spatraster<-rast(grja_no_overlap)
crs(current_overlap_spatraster)<-"EPSG:4326"
grja_no_overlap_polygon<-terra::as.polygons(grja_no_overlap_spatraster)
plot(grja_no_overlap_polygon)
grja_no_overlap_sf <- st_as_sf(grja_no_overlap_polygon)
st_crs(grja_no_overlap_sf)<-"EPSG:4326"

### FIGURE 2C ###
# NOTE: I removed guides and was forced to add a working guide in post-processing via Adobe Illustrator. I battled valiantly against the mighty powers of ggplot2 in order to produce a working legend, but I like many grad students before me was defeated. Please reach out if you know how to produce an adequate legend that tracks the patterns used...
fig2c<-ggplot() +
  geom_sf(data = regional_map_countries, color = "black", fill = "white", lwd = 0.75, aes(linetype = "Country Borders"), show.legend = FALSE) + 
  
  # BLJA Current Range Layer (Raster)
  geom_sf_pattern(data = blja_no_overlap_sf, aes(fill = "BLJA Range - 2023", pattern = "BLJA Range - 2023", color = "BLJA Range - 2023", pattern_fill = "BLJA Range - 2023", pattern_density = "BLJA Range - 2023", alpha = "BLJA Range - 2023"), color = "black", inherit.aes = FALSE, show.legend = FALSE) +
  # GRJA Current Range Layer (Raster)
  geom_sf_pattern(data = grja_no_overlap_sf, aes(fill = "GRJA Range - 2023", pattern = "GRJA Range - 2023", color = "GRJA Range - 2023", pattern_fill = "GRJA Range - 2023", pattern_density = "GRJA Range - 2023", alpha = "GRJA Range - 2023"), color = "black", inherit.aes = FALSE, show.legend = FALSE) +
  # BOTH SPP Current Range Layer (Raster)
  geom_sf_pattern(data = current_overlap_sf, aes(fill = "Range Overlap - 2023", pattern = "Range Overlap - 2023", color = "Range Overlap - 2023", pattern_fill = "Range Overlap - 2023", pattern_density = "Range Overlap - 2023", pattern_alpha = "Range Overlap - 2023", alpha = "Range Overlap - 2023"), color = "black", inherit.aes = FALSE, show.legend = FALSE) +
  
  # Regional Map Layers (Countries, US States, MX States)
  geom_sf(data = US_states, color = "black", fill = "transparent", lwd = 0.25, aes(linetype = "US State Borders"), show.legend = FALSE) + 
  geom_sf(data = MX_states, color = "black", fill = "transparent", lwd = 0.25, aes(linetype = "MX State Borders"), show.legend = FALSE) +
  geom_sf(data = regional_map_countries, color = "black", fill = "transparent", lwd = 0.75, aes(linetype = "Country Borders"), show.legend = FALSE) + 
  
  # BLJA Future Range Prediction (SpatVector Polygon)
  geom_spatvector(data = blja_future_range_polygon, aes(color = "BLJA Range ~ 2050 (ssp245)", fill = "BLJA Range ~ 2050 (ssp245)", alpha = "BLJA Range ~ 2050 (ssp245)"), lwd = 0.5, show.legend = FALSE) +
  # GRJA Future Range Prediction (SpatVector Polygon)
  geom_spatvector(data = grja_future_range_polygon, aes(color = "GRJA Range ~ 2050 (ssp245)", fill = "GRJA Range ~ 2050 (ssp245)", alpha = "GRJA Range ~ 2050 (ssp245)"), lwd = 0.5, show.legend = FALSE) +
  
  # Coordinate Limits
  coord_sf(xlim = c(-101.5, -95.5), ylim = c(26, 31)) +
  
  # Custom Fill Legend
  scale_fill_manual(name = "Legend",
                    values = c(
                      "BLJA Range - 2023" = "white",
                      "GRJA Range - 2023" = "#81C784",
                      "Range Overlap - 2023" = "#81C784",
                      "GRJA Range ~ 2050 (ssp245)" = "darkgreen",
                      "BLJA Range ~ 2050 (ssp245)" = "blue"),
                    breaks = c("BLJA Range - 2023", "GRJA Range - 2023", "Range Overlap - 2023"),
                    labels = c("BLJA Range - 2023", "GRJA Range - 2023", "Range Overlap - 2023")) +
  #custom alpha
  scale_alpha_manual(name = "Legend",
                   values = c(
                     "BLJA Range - 2023" = 0.5,
                     "GRJA Range - 2023" = 0.5,
                     "Range Overlap - 2023" = 0.5,
                     "GRJA Range ~ 2050 (ssp245)" = 0.25,
                     "BLJA Range ~ 2050 (ssp245)" = 0.35)) +
  # Custom Pattern
  scale_pattern_manual(name = "Legend",
                       values = c(
                         "BLJA Range - 2023" = "stripe",
                         "GRJA Range - 2023" = "stripe",
                         "Range Overlap - 2023" = "stripe"),
                       breaks = c("BLJA Range - 2023", "GRJA Range - 2023", "Range Overlap - 2023"),
                       labels = c("BLJA Range - 2023", "GRJA Range - 2023", "Range Overlap - 2023")) + 
  
  # Custom Pattern FILL
  scale_pattern_fill_manual(name = "Legend",
                            values = c(
                              "BLJA Range - 2023" = "#64B5F6",
                              "GRJA Range - 2023" = "white",
                              "Range Overlap - 2023" = "#64B5F6"),
                            breaks = c("BLJA Range - 2023", "GRJA Range - 2023", "Range Overlap - 2023"),
                            labels = c("BLJA Range - 2023", "GRJA Range - 2023", "Range Overlap - 2023")) + 
  
  scale_pattern_alpha_manual(name = "Legend",
                             values = c(
                               "BLJA Range - 2023" = 0.5,
                               "GRJA Range - 2023" = 0.5,
                               "Range Overlap - 2023" = 0.5)) +
  
  # Custom Pattern DENSITY
  scale_pattern_density_manual(name = "Legend",
                               values = c(
                                 "BLJA Range - 2023" = 0.5,
                                 "GRJA Range - 2023" = 0.5,
                                 "Range Overlap - 2023" = 0.5)) + 
  
  # Custom Linetype Legend for Borders
  scale_linetype_manual(name = "Borders",
                        values = c(
                          "Country Borders" = "solid",
                          "US State Borders" = "dotted",
                          "MX State Borders" = "dotted"
                        )) +
  
  scale_color_manual(name = "futr",
                     values = c("GRJA Range ~ 2050 (ssp245)" = "darkgreen",
                                "BLJA Range ~ 2050 (ssp245)" = "blue"),
                     breaks = c("BLJA Range ~ 2050 (ssp245)", "GRJA Range ~ 2050 (ssp245)"),
                     labels = c("BLJA Range ~ 2050 (ssp245)", "GRJA Range ~ 2050 (ssp245)"))+
  
  # Titles
  ggtitle("Blue & Green Jay Range Overlap") +
  xlab("Latitude") +
  ylab("Longitude") +
  
  #Add theme
  theme(panel.background = element_rect(fill = "lightblue"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.title=element_blank(),
        aspect.ratio=1,
        plot.title = element_text(size=20, face = "bold"),
        axis.title = element_text(size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.position = "none")

fig2c






#### Figure 2d ####

fig2d<-ggplot() +
  geom_sf(data = regional_map_countries, color = "black", fill = "white", lwd = 0.75, aes(linetype = "Country Borders"), show.legend = FALSE) + 
  # Regional Map Layers (Countries, US States, MX States)
  geom_sf(data = US_states, color = "black", fill = "transparent", lwd = 0.25, aes(linetype = "US State Borders"), show.legend = FALSE) + 
  geom_sf(data = MX_states, color = "black", fill = "transparent", lwd = 0.25, aes(linetype = "MX State Borders"), show.legend = FALSE) +
  geom_sf(data = regional_map_countries, color = "black", fill = "transparent", lwd = 0.75, aes(linetype = "Country Borders"), show.legend = FALSE) + 
  
  # Coordinate Limits
  coord_sf(xlim = c(-101.5, -95.5), ylim = c(26, 31)) +
  geom_point(data=grja_presence_jun2023, aes(y=latitude,x=longitude, col = 'Green Jay only'), shape=1, size=1.5) +
  geom_point(data=blja_presence_jun2023, aes(y=latitude,x=longitude, col = 'Blue Jay only'), shape=1, size=1.5) +
  geom_point(data=bothspp_filter, aes(y=latitude,x=longitude, col = 'Both Species'), shape=19, size=2) +
  
  ggtitle("Blue & Green Jay Observations - eBird") +
  xlab("Latitude") +
  ylab("Longitude") +
  
  theme(panel.background = element_rect(fill = "lightblue"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        aspect.ratio=1,
        plot.title = element_text(size=20, face = "bold"),
        axis.title = element_text(size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.background = element_blank(),
        legend.position = c(0.17,0.90),
        legend.box.background = element_rect(colour = "black"),
        legend.key = element_rect(fill = "white")) +
  
  scale_color_manual(name='eBird Checklist Observations',
                     breaks=c('Blue Jay only', 'Green Jay only', 'Both Species'),
                     values=c('Blue Jay only'='#64B5F6', 'Green Jay only'='#81C784', 'Both Species'='#000000'),
                     aesthetics = c("color")) +

  guides(color=guide_legend(override.aes=list(fill="white")))


fig2d
  
  
  
### Optionally put figures together in a single grid
grid.arrange(fig2a, fig2b, fig2c, fig2d, nrow = 2)
