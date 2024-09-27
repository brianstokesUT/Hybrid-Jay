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

#### Create Environmental dataframe ####

# Current Env Data # https://adaptwest.databasin.org/pages/adaptwest-climatena/
# List all .tif files in the directory
current_raster_files <- list.files("Normal_1991_2020_bioclim", pattern = "Normal_1991_2020_.*\\.tif$", full.names = TRUE)

# Initialize an empty list to store the cropped and reprojected rasters
current_cropped_rasters <- list()
# Define the cropping extent in WGS 84 (lon/lat) - within our study area to save some memory
crop_extent <- extent(-105, -93, 20, 35)

# Loop through the raster files, read, reproject, crop, and store them in the list
for (i in seq_along(current_raster_files)) {
  current_r <- raster(current_raster_files[i])  # Read the raster file
  
  # Check the CRS of the raster
  print(crs(current_r))  # This will tell you the current projection system
  
  current_r <- projectRaster(current_r, crs = 4326)  # Reproject to 4326
  
  # Now crop the reprojected raster
  current_cropped_rasters[[i]] <- crop(current_r, crop_extent)
}
# Stack all the cropped rasters into a RasterStack
Normal_1991_2020_stack <- raster::stack(current_cropped_rasters)

#get names for each layer of raster stack
raster_names <- gsub("Normal_1991_2020_|\\.tif$", "", basename(current_raster_files))
# Set names for the RasterStack layers
names(Normal_1991_2020_stack) <- raster_names

# Future Env Data (ssp245 2041-2060)
# List all .tif files in the directory
future_raster_files <- list.files("ensemble_8GCMs_ssp245_2041_2060_bioclim", pattern = "ensemble_8GCMs_ssp245_2041_2060_.*\\.tif$", full.names = TRUE)

# Initialize an empty list to store the cropped and reprojected rasters
cropped_future_rasters <- list()

# Extract unique identifiers from file names (e.g., "Normal_1991_2020_1.tif" -> "1")
future_raster_names <- gsub("ensemble_8GCMs_ssp245_2041_2060_|\\.tif$", "", basename(future_raster_files))

# Loop through the raster files, read, reproject, crop, and store them in the list
for (i in seq_along(future_raster_files)) {
  future_r <- raster(future_raster_files[i])  # Read the raster file
  
  # Check the CRS of the raster
  print(crs(future_r))  # This will tell you the current projection system
  
  future_r <- projectRaster(future_r, crs = 4326)  # Reproject to 4326
  
  # Now crop the reprojected raster
  cropped_future_rasters[[i]] <- crop(future_r, crop_extent)
}

# need to remove _MAR file becasue its missing from current dataset
cropped_future_rasters<-cropped_future_rasters[-16]
future_raster_names<-future_raster_names[-16]

# stack for Maxent formatting
future_stack<-raster::stack(cropped_future_rasters)
# Set names for the RasterStack layers
names(future_stack) <- future_raster_names

# OPTIONAL: Write out Raster
writeRaster(Normal_1991_2020_stack, filename="Normal_1991_2020_stack.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
writeRaster(future_stack, filename="future_stack.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
