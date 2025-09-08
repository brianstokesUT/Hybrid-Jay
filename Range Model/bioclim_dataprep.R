# ClimateNA preprocessing (terra-only): current (1991–2020) and future (SSP245 2041–2060)
# Author: Brian R. Stokes
# Description: Read, reproject to EPSG:4326, crop, harmonize variable sets, and extract MAT at a point.

wd<-("~/PATH")
setwd(wd)

set.seed(223)
options(java.parameters = "-Xmx8000m")

if (!requireNamespace("terra", quietly = TRUE) ||
    as.character(utils::packageVersion("terra")) != 1.8-60) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  remotes::install_version("terra", version = "1.8-60", upgrade = "never")
}

suppressPackageStartupMessages({
  library(terra)  # now guaranteed version
})


# -----------------------
# 1) User-configurable paths & settings
# -----------------------
# If you prefer a working directory, uncomment:
# setwd("~/PATH")


# Input directories
cur_dir <- "Normal_1991_2020_bioclim"
fut_dir <- "ensemble_8GCMs_ssp245_2041_2060_bioclim"


# Filename patterns (kept flexible)
cur_pat <- "Normal_1991_2020_.*\\.tif$"
fut_pat <- "ensemble_8GCMs_ssp245_2041_2060_.*\\.tif$"


# Crop extent in lon/lat (EPSG:4326)
# Study-area window to reduce memory/IO
crop_window <- ext(-105, -93, 20, 35)


# Target CRS
target_crs <- "EPSG:4326"


# Point for extraction (lon, lat)
hyb_lonlat <- c(-98.31, 29.54) # ~San Antonio, TX (minor jitter/low res for privacy)


# Choose which variable to extract at hyb_lonlat (prefer "MAT", fallback to "bio1")
preferred_vars <- c("MAT", "bio1")


# Optional: write outputs?
write_outputs <- FALSE
cur_out_tif <- "Normal_1991_2020_stack.tif"
fut_out_tif <- "future_ssp245_2041_2060_stack.tif"


# -----------------------
# 2) Helpers
# -----------------------
# Read all .tif files in a directory by pattern, name layers from filenames (stripped prefix/suffix)
read_named_stack <- function(dir_path, file_pat, strip_prefix) {
  files <- list.files(dir_path, pattern = file_pat, full.names = TRUE)
  if (length(files) == 0) stop(sprintf("No files matched in: %s", dir_path))
  # Derive layer names from file basenames
  nm <- basename(files)
  nm <- sub(paste0("^", strip_prefix), "", nm)
  nm <- sub("\\.tif$", "", nm)
  r <- rast(files) # one layer per file
  if (nlyr(r) != length(nm)) {
    warning("Number of layers != number of file-derived names; assigning what we can.")
  }
  names(r) <- nm
  r
}


# Reproject + crop to a lon/lat window (EPSG:4326)
reproject_and_crop <- function(r, target_crs, window_ext) {
  # Project (bilinear for continuous climate surfaces)
  if (!same.crs(r, target_crs)) {
    r <- project(r, target_crs, method = "bilinear")
  }
  crop(r, window_ext)
}


# Pick the first available variable from a preference list
pick_var <- function(r, candidates) {
  found <- intersect(candidates, names(r))
  if (length(found) == 0) stop(sprintf("None of %s found in layer names.", paste(candidates, collapse = ", ")))
  found[[1]]
}


# -----------------------
# 3) Read, project, crop
# -----------------------
cur <- read_named_stack(
  dir_path = cur_dir,
  file_pat = cur_pat,
  strip_prefix= "Normal_1991_2020_"
)
cur <- reproject_and_crop(cur, target_crs, crop_window)


fut <- read_named_stack(
  dir_path = fut_dir,
  file_pat = fut_pat,
  strip_prefix= "ensemble_8GCMs_ssp245_2041_2060_"
)
fut <- reproject_and_crop(fut, target_crs, crop_window)


# -----------------------
# 4) Harmonize variable sets (robust to missing layers like "MAR" which is not present in future dataset)
# -----------------------
common_vars <- intersect(names(cur), names(fut))
if (length(common_vars) == 0) stop("No overlapping layer names between current and future stacks.")
cur <- cur[[common_vars]]
fut <- fut[[common_vars]]


# -----------------------
# 5) Extract MAT (or bio1) at the hybrid location & compute delta
# -----------------------
hyb_pt <- vect(data.frame(lon = hyb_lonlat[1], lat = hyb_lonlat[2]), geom = c("lon", "lat"), crs = target_crs)
mat_name <- pick_var(cur, preferred_vars)


cur_mat_val <- as.numeric(extract(cur[[mat_name]], hyb_pt)[,2])
fut_mat_val <- as.numeric(extract(fut[[mat_name]], hyb_pt)[,2])
delta_mat <- fut_mat_val - cur_mat_val


cat(sprintf("%s at hybrid site (%.4f, %.4f)\n", mat_name, hyb_lonlat[1], hyb_lonlat[2]))
cat(sprintf(" • Current: %s\n", format(cur_mat_val, digits = 6)))
cat(sprintf(" • Future : %s\n", format(fut_mat_val, digits = 6)))
cat(sprintf(" • ΔFuture-Current: %s\n", format(delta_mat, digits = 6)))


# -----------------------
# 6) Optional: write compressed GeoTIFFs
# -----------------------
if (isTRUE(write_outputs)) {
  writeRaster(cur, cur_out_tif, overwrite = TRUE, gdal = c("COMPRESS=LZW", "TILED=YES"))
  writeRaster(fut, fut_out_tif, overwrite = TRUE, gdal = c("COMPRESS=LZW", "TILED=YES"))
}


# End of script
