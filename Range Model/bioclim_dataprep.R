wd<-("~/PATH")
setwd(wd)

set.seed(223)
options(java.parameters = "-Xmx8000m")


library("auk")
library("dplyr")
library("ggplot2")
library("gridExtra")
library("lubridate")
library("readr")
library("sf")
library("geosphere")
library("terra")
library("rJava")
library("dismo")

#### Create Environmental Dataframe ####

## ------------------------------------------------------------------
## Current Env Data (Normal 1991–2020)
## Source: https://adaptwest.databasin.org/pages/adaptwest-climatena/
## ------------------------------------------------------------------

# List all .tif files in the directory
current_raster_files <- list.files(
  "Normal_1991_2020_bioclim",
  pattern    = "Normal_1991_2020_.*\\.tif$",
  full.names = TRUE
)

# Initialize an empty list to store the cropped and reprojected rasters
current_cropped_rasters <- list()

# Define the cropping extent in WGS 84 (lon/lat) – keeps within study area to save memory
crop_extent <- extent(-105, -93, 20, 35)

# Loop through the raster files: read, reproject, crop, store
for (i in seq_along(current_raster_files)) {
  current_r <- raster(current_raster_files[i])  # Read the raster file
  
  # Check the CRS of the raster (diagnostic)
  print(crs(current_r))
  
  # Reproject to WGS84 (EPSG:4326)
  current_r <- projectRaster(current_r, crs = 4326)
  
  # Crop the reprojected raster
  current_cropped_rasters[[i]] <- crop(current_r, crop_extent)
}

# Stack all the cropped rasters into a RasterStack
Normal_1991_2020_stack <- raster::stack(current_cropped_rasters)

# Get names for each layer of raster stack
raster_names <- gsub(
  pattern = "Normal_1991_2020_|\\.tif$",
  replacement = "",
  x = basename(current_raster_files)
)

# Set names for the RasterStack layers
names(Normal_1991_2020_stack) <- raster_names


## ------------------------------------------------------------------
## Future Env Data (ssp245, 2041–2060)
## ------------------------------------------------------------------

# List all .tif files in the directory
future_raster_files <- list.files(
  "ensemble_8GCMs_ssp245_2041_2060_bioclim",
  pattern    = "ensemble_8GCMs_ssp245_2041_2060_.*\\.tif$",
  full.names = TRUE
)

# Initialize an empty list to store the cropped and reprojected rasters
cropped_future_rasters <- list()

# Extract names from file names
future_raster_names <- gsub(
  pattern     = "ensemble_8GCMs_ssp245_2041_2060_|\\.tif$",
  replacement = "",
  x           = basename(future_raster_files)
)

# Loop through the raster files: read, reproject, crop, store
for (i in seq_along(future_raster_files)) {
  future_r <- raster(future_raster_files[i])  # Read the raster file
  
  # Check the CRS of the raster (diagnostic)
  print(crs(future_r))
  
  # Reproject to WGS84 (EPSG:4326)
  future_r <- projectRaster(future_r, crs = 4326)
  
  # Crop the reprojected raster
  cropped_future_rasters[[i]] <- crop(future_r, crop_extent)
}

# Need to remove _MAR file because it is missing from current dataset
cropped_future_rasters <- cropped_future_rasters[-16]
future_raster_names    <- future_raster_names[-16]

# Stack for Maxent formatting
future_stack <- raster::stack(cropped_future_rasters)

# Set names for the RasterStack layers
names(future_stack) <- future_raster_names


## ------------------------------------------------------------------
## MAT values at hybrid locality (for manuscript)
## ------------------------------------------------------------------

# Approximate coordinates (29.54, -98.31) – slightly shifted to protect privacy
cur_mat_layer <- Normal_1991_2020_stack[["MAT"]]
fut_mat_layer <- future_stack[["MAT"]]

# Coordinates must be: longitude (x), latitude (y)
hyb_coords <- matrix(c(-98.31, 29.54), ncol = 2)

# Extract values
raster::extract(cur_mat_layer, hyb_coords)
raster::extract(fut_mat_layer, hyb_coords)


## ------------------------------------------------------------------
## OPTIONAL: Write out rasters
## ------------------------------------------------------------------

# writeRaster(
#   Normal_1991_2020_stack,
#   filename = "Normal_1991_2020_stack.tif",
#   options  = "INTERLEAVE=BAND",
#   overwrite = TRUE
# )
#
# writeRaster(
#   future_stack,
#   filename = "future_stack.tif",
#   options  = "INTERLEAVE=BAND",
#   overwrite = TRUE
# )
