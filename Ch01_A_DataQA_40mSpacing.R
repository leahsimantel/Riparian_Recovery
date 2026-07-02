#########################################################
###    CH01 Script A: Quality Assurance and computation of covariates
###    11/08/2025 - - 40m spacing (fixed 11/8)
###    03/02/2026 - - Updated QA for 40m spacing and deleted another 33 pixels. 

# Load necessary libraries -------------------------
library(dplyr)
library(readr)
library(ggplot2)
library(sf)
library(tidyr)
library(stringr)
library(e1071)  # for skewness
library(scales)
library(moments)  # for skewness
library(geojsonsf)
library(writexl)
library(haven)
library(foreign)
# --------------------------------------------------

# Define the folder containing the CSV files
folder_path <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/NDVI_covariates_by_Pixel/CSVs_08202025"

# Get list of all CSV files in the folder
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# Read and combine all CSV files (let read_csv guess column types)
raw_data_fresh <- bind_rows(lapply(csv_files, read_csv))

################################################################################
################################################################################
##################       DATA CLEANUP / / QA     ###############################
################################################################################
################################################################################

### Explicitly set column types for downstream analysis: -----------------------

# If this shows up anywhere, drop it
raw_data <- raw_data_fresh %>% select(-any_of("system:index"))

raw_data <- raw_data %>%
  mutate(
    # IDs / categorical
    COMID        = format(as.character(COMID)),
    huc12        = as.character(huc12),   
    FCODE        = as.character(FCODE),
    pixel_ID     = as.character(pixel_ID),
    fire_name    = as.character(fire_name),
    sev_class    = as.character(sev_class),
    .geo         = as.character(.geo),
    
    # Core numerics
    fire_size_ha = as.numeric(fire_size_ha),
    hli          = as.numeric(hli),
    twi          = as.numeric(twi),
    fire_year    = as.integer(fire_year),  ## set fire year as integer
    num_of_burns = as.integer(num_of_burns),
    latitude     = as.numeric(latitude),
    longitude    = as.numeric(longitude),
    sev_num      = as.numeric(sev_num),
    sev_cannon_num = as.numeric(.data[["sev_cannon_num"]])
  ) %>%
  
  # Time-series numeric columns by regex patterns
  mutate(
    across(
      matches("^SWE_\\d{6}_TC$|^AET_\\d{6}_TC$|^CWD_\\d{6}_TC$"),
      ~ suppressWarnings(as.numeric(.))
    ),
    across(
      matches("^NDVI_\\d{4}$"),
      ~ suppressWarnings(as.numeric(.))  ## set NDVI values to numeric
    ),
    across(
      matches("^cwd_30yrAvg_TC_0[4-8]$"),
      ~ suppressWarnings(as.numeric(.))
    ),
    across(
      matches("^CWD_\\d{4}_TT$"),
      ~ suppressWarnings(as.numeric(.))
    ),
    across(
      matches("^(ppt|tmean|tmax)_\\d{6}_PR$"),
      ~ suppressWarnings(as.numeric(.))
    )
  )

# Print to verify
#head(raw_data)
#colnames(raw_data)

### 1. Merge 'raw_data' with 'streams_widened' data to attach Stream Order, Wetted width -------

# Read DBF
dbf_path <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/ArcGIS/BMWA_Streams_Widened_Ch01.dbf"

streams_widened_raw <- foreign::read.dbf(dbf_path, as.is = TRUE) %>%
  as_tibble()

# -------------------------------------------------------
# 2. Select needed columns
# -------------------------------------------------------
streams_widened <- streams_widened_raw %>%
  select(
    COMID,
    LENGTHK,
    TtDASKM,
    WETTEDW,
    WttdW_S,
    StrmOrd,
    flw_ctg
  )

# -------------------------------------------------------
# 3. Convert column types appropriately
# -------------------------------------------------------
streams_widened <- streams_widened %>%
  mutate(
    COMID         = as.character(COMID),
    LENGTHK       = as.numeric(LENGTHK),
    TtDASKM       = as.numeric(TtDASKM),
    wetted_width_m = as.numeric(WETTEDW),  # rename + convert
    WttdW_S       = as.character(WttdW_S),
    flw_ctg       = as.character(flw_ctg),
    stream_order       = as.integer(StrmOrd)
  ) %>%
  select(-WETTEDW, -StrmOrd)  # drop original after renaming

# -------------------------------------------------------
# 4. Merge using COMID
# -------------------------------------------------------
raw_data_merged <- raw_data %>%
  left_join(streams_widened,
            by = "COMID")


################################################################################
###### =============   Checking out column values   ===================  #######
################################################################################

####################  FIRE SEVERITY CLASS AND VALUE  ##########################

# Switch out 'NA' in sev_class to 'None' instead
raw_data <- raw_data_merged %>%
  mutate(sev_class = if_else(is.na(sev_class), "None", sev_class))

raw_data <- raw_data %>%
  mutate(
    sev_class = case_when(
      is.na(sev_num) ~ "None",
      sev_num == 0 ~ "Unburned/Very Low",
      sev_num > 0 & sev_num <= 1.25 ~ "Low",
      sev_num > 1.25 & sev_num <= 2.25 ~ "Moderate",
      sev_num > 2.25 ~ "High",
      TRUE ~ NA_character_  # fallback, just in case
    )
  )

# Replace placeholder -9999 with NA in sev_num
raw_data <- raw_data %>%
  mutate(sev_num = ifelse(sev_num == -9999, NA_real_, sev_num))

# Identify rows where sev_class is "None" but sev_num is not NA
missing_class_with_value <- raw_data %>%
  filter(sev_class == "None" & !is.na(sev_num))

# Count how many rows meet this condition
n_rows <- nrow(missing_class_with_value)
cat("Number of rows where sev_class is NA but sev_num is not:", n_rows, "\n") #0

# Define severity bins from sev_num
sev_num_bins <- raw_data %>%
  mutate(
    sev_num_class = case_when(
      is.na(sev_num)                  ~ "None",
      sev_num == 0                    ~ "Unburned/Very Low",
      sev_num > 0 & sev_num <= 1.25   ~ "Low",
      sev_num > 1.25 & sev_num <= 2.25~ "Moderate",
      sev_num > 2.25                  ~ "High",
      TRUE                            ~ "Unclassified"
    )
  ) %>%
  count(sev_num_class, name = "n_sev_num") %>%
  rename(sev_class = sev_num_class)

# Count sev_class directly (with NA replaced for join compatibility)
sev_class_counts <- raw_data %>%
  mutate(sev_class = if_else(is.na(sev_class), "None", sev_class)) %>%
  count(sev_class, name = "n_sev_class")

# Join and arrange
severity_comparison <- full_join(sev_class_counts, sev_num_bins, by = "sev_class") %>%
  arrange(desc(n_sev_class + n_sev_num))

# Print
print(severity_comparison, n = Inf)

# A tibble: 5 × 3
# sev_class            n_sev_class    n_sev_num
#   <chr>                   <int>     <int>
#  1 None                    36440     36440
# 2 Low                      8260      8260
# 3 Moderate                 7809      7809
# 4 High                     6058      6058
# 5 Unburned/Very Low        1609      1609

###################### ASSIGN VALUES TO CANNON FIRE PIXELS   ###################
# sev_cannon_num column was messed up! Need to upload TIF and re-assign

library(terra)
library(sf)
library(dplyr)
library(geojsonsf)

### clean slate: set all values in sev_cannon_num to 'NA'
raw_data$sev_cannon_num <- NA_real_

# Check unique values and their counts in sev_cannon_num
raw_data %>%
  count(sev_cannon_num, sort = TRUE) %>%
  print(n = Inf) 

sev_cannon_summary <- raw_data %>%
  count(sev_cannon_num, sort = TRUE)

#print(sev_cannon_summary, n = Inf)  ## NA: 60176

## === PART 1: ready the raw data 
# Convert GeoJSON in `.geo` column to sfc geometry
raw_data_sf <- raw_data %>%
  filter(!is.na(.geo)) %>%
  mutate(geometry = geojsonsf::geojson_sfc(.geo)) %>%
  st_sf()
#st_crs(raw_data_sf)  # currently in WGS84 (EPSG:4326) (will measure distances in degrees instead of m)

# Reproject raw_data to UTM Zone 12N
raw_data_utm <- st_transform(raw_data_sf, crs = 32612)


## === PART 2: ready the sev_Cannon fire severity layer
# Path to the raster file
cannon_raster_path <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/ArcGIS/BMWC_Fire_Data/fire.severity_LYB_freq/Fire Severity BMWC/MT4771611356020220807_CannonFire_Severity.TIF"

# Load raster
cannon_severity_rast <- rast(cannon_raster_path)

# Check the raster summary
print(cannon_severity_rast)
plot(cannon_severity_rast, main = "Cannon Fire Severity Raster")

# Reproject raster to match UTM Zone 12N (EPSG:32612)
cannon_severity_rast_utm <- terra::project(cannon_severity_rast, "EPSG:32612")

# Quick check: plot both raster and points
plot(cannon_severity_rast_utm, main = "Cannon Fire Severity (UTM)")
plot(st_geometry(raw_data_utm), add = TRUE, col = "blue", pch = 20, cex = 0.3)

## === PART 3: Extract & assign severity values for the Cannon fire
# Convert sf to terra vector
raw_data_vect <- terra::vect(raw_data_utm)

# Extract values from Cannon sev raster
extracted_vals <- terra::extract(cannon_severity_rast_utm, raw_data_vect)

# Identify valid cells with non-NA severity values
valid_cells <- !is.na(extracted_vals$MT4771611356020220807_CannonFire_Severity)

# Subset the spatial data to include only those valid points
valid_points <- raw_data_utm[valid_cells, ]

# Extract corresponding pixel_IDs and severity values
valid_pixel_ids <- valid_points$pixel_ID
valid_severity_vals <- extracted_vals$MT4771611356020220807_CannonFire_Severity[valid_cells]

# Build lookup table
severity_lookup <- data.frame(
  pixel_ID = valid_pixel_ids,
  sev_cannon_num = valid_severity_vals
)

# Assign severity values back into raw_data by pixel_ID
raw_data <- raw_data %>%
  left_join(severity_lookup, by = "pixel_ID", suffix = c("", ".new")) %>%
  mutate(
    sev_cannon_num = coalesce(sev_cannon_num.new, sev_cannon_num)
  ) %>%
  select(-sev_cannon_num.new)  ## remove dummy column from raw_data

# Use the same logic as before to classify numeric severity into classes
raw_data <- raw_data %>%
  mutate(
    sev_class = case_when(
      !is.na(sev_cannon_num) & sev_cannon_num == 0                    ~ "Unburned/Very Low",
      !is.na(sev_cannon_num) & sev_cannon_num > 0 & sev_cannon_num <= 1.25   ~ "Low",
      !is.na(sev_cannon_num) & sev_cannon_num > 1.25 & sev_cannon_num <= 2.25~ "Moderate",
      !is.na(sev_cannon_num) & sev_cannon_num > 2.25                  ~ "High",
      TRUE                                                            ~ sev_class  # retain existing value
    )
  )

## === PART 4: CHECK RESULTS

# Check unique values and their counts in sev_cannon_num
raw_data %>%
  count(sev_cannon_num, sort = TRUE) %>%
  print(n = Inf) 

sev_cannon_summary <- raw_data %>%
  count(sev_cannon_num, sort = TRUE)

#print(sev_cannon_summary, n = Inf)  ## 126 pixels in Cannon fire

## === PART 5: PLOT CLASSIFIED PIXELS

# Custom color palette
sev_colors <- c(
  "None"              = "#d9d9d9",  
  "Unburned/Very Low" = "#4d4d4d",  
  "Low"               = "#ffff00",  
  "Moderate"          = "#ffa500",  
  "High"              = "#ff0000"   
)

# Subset just the pixels with a cannon severity value
cannon_pixels <- raw_data %>%
  filter(!is.na(sev_cannon_num)) %>%
  filter(!is.na(.geo)) %>%
  mutate(geometry = geojsonsf::geojson_sfc(.geo)) %>%
  st_sf(crs = 4326) %>%
  st_transform(crs = 32612)  # UTM Zone 12N for spatial alignment

#plot(cannon_severity_rast_utm, main = "Cannon Fire Severity with Classified Pixels")
#plot(
#  st_geometry(cannon_pixels),
#  add = TRUE,
#  col = sev_colors[cannon_pixels$sev_class],
#  pch = 20,
#  cex = 0.7
#)

# Confirm number of pixels in cannon_pixels
n_cannon <- nrow(cannon_pixels)
cat("Number of Cannon fire pixels:", n_cannon, "\n") #126

## === PART 6: ASSIGN FIRE NAME, YEAR TO CANNON PIXELS

# Count unique values in num_of_burns for Cannon pixels
num_burns_summary <- cannon_pixels %>%
  count(num_of_burns, sort = TRUE)

#print(num_burns_summary)

# Update burn counts for cannon_pixels
cannon_pixels <- cannon_pixels %>%
  mutate(num_of_burns = case_when(
    num_of_burns == 0 ~ 1L,
    num_of_burns == 1 ~ 2L,
    TRUE              ~ num_of_burns  # keep other values as-is (just in case)
  ))

num_burns_summary <- cannon_pixels %>%
  count(num_of_burns, sort = TRUE) 
print(num_burns_summary)

# Assign fire info to Cannon pixels where this is their only known fire
cannon_pixels <- cannon_pixels %>%
  mutate(
    fire_year     = ifelse(num_of_burns == 1, 2023L, fire_year),
    fire_name     = ifelse(num_of_burns == 1, "CANNON_2023", fire_name),
    fire_size_ha  = ifelse(num_of_burns == 1, 863.54, fire_size_ha)
  )

# Confirm updated values
cannon_pixels %>%
  filter(num_of_burns == 1) %>%
  count(fire_name, fire_year, fire_size_ha)

# Drop any previous Cannon pixel entries from raw_data
raw_data <- raw_data %>%
  filter(!(pixel_ID %in% cannon_pixels$pixel_ID))  

# Combine updated cannon_pixels back in
raw_data <- bind_rows(raw_data, cannon_pixels)

# Confirm all Cannon pixels are present and updated
raw_data %>%
  filter(fire_name == "CANNON_2023") %>%
  count(fire_year, fire_name, fire_size_ha, num_of_burns)  # 84 pixels...


#########  QA AND DATA CLEANUP: DELETE DUPLICATES, FIRE DISCREPANCIES  #########

## === PART 1: CHECK FOR DUPLICATE PIXELS AND REMOVE THEM
# Check for duplicate values in pixel_ID column
num_duplicate_pixel_IDs <- raw_data %>%
  count(pixel_ID) %>%
  filter(n > 1) %>%
  nrow()  # Count the number of pixel_IDs that have duplicates

# Print the number of duplicates found:
#cat("Number of duplicate pixel_IDs:", num_duplicate_pixel_IDs, "\n")   #Number of duplicate pixel_IDs: 2270  

# Remove duplicate pixel_IDs, keeping the first occurrence
raw_data <- raw_data %>%
  distinct(pixel_ID, .keep_all = TRUE)

# Recount duplicates to confirm removal - - VERIFY
num_duplicate_pixel_IDs <- raw_data %>%
  count(pixel_ID) %>%
  filter(n > 1) %>%
  nrow()

cat("Number of duplicate pixel_IDs after removal:", num_duplicate_pixel_IDs, "\n")
#Number of duplicate pixel_IDs after removal: 0 

## === PART 2: FILTERING: ONCE BURNED / UNBURNED 
# Note - deleting all twice and thrice-burned from the dataset

# Delete all rows where num_of_burns is NOT 0 or 1:
rows_before <- nrow(raw_data)  # Count total rows before filtering

once_and_unb_data <- raw_data %>%
  filter(num_of_burns %in% c(0, 1))  # Keep only rows where num_of_burns is 0 or 1

# Count total rows after filtering
rows_after <- nrow(once_and_unb_data)  
num_deleted_rows <- rows_before - rows_after  # Calculate the number of deleted rows
cat("Number of rows deleted:", num_deleted_rows, "\n")   #  Number of rows deleted: 5065 

# Rename once_and_unb_data back to raw_data
raw_data <- once_and_unb_data

# Switch out 'NA' in sev_class to 'None' instead
raw_data <- raw_data %>%
  mutate(sev_class = if_else(is.na(sev_class), "None", sev_class))

# Summarize how many pixels are left in each severity class
sev_summary <- raw_data %>%
  count(sev_class, sort = TRUE)

print(sev_summary)
#  A tibble: 5 × 2
#  sev_class             n
# <chr>               <int>
# 1 None              34230
# 2 Low                7084
# 3 Moderate           5380
# 4 High               4418
# 5 Unburned/Very Low  1543

#######  QUICK CHECK ON SPATIAL ALIGNMENT: PIXELS  +  FIRE ATLAS POLYGONS ######
raw_data_sf <- raw_data %>%
  filter(!is.na(.geo)) %>%
  mutate(geometry = geojsonsf::geojson_sfc(.geo)) %>%
  st_sf(crs = 4326)   

# Make sure severity is a factor with levels in the same order as sev_colors
raw_data_sf$sev_class <- factor(raw_data_sf$sev_class, levels = names(sev_colors))

#ggplot(raw_data_sf) +
#  geom_sf(aes(color = sev_class), size = 0.6, alpha = 0.9) +
#  scale_color_manual(values = sev_colors, drop = FALSE) +
#  labs(
#    title = "Pixels by Fire Severity",
#   color = "Severity"
#  ) +
#  theme_minimal() 

### add the fire atlas to the map:
library(sf)
library(dplyr)
library(ggplot2)

atlas_dir  <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/ArcGIS/BMWC_Fire_Data/fire_atlas"
atlas_path <- file.path(atlas_dir, "fire_atlas.shp")
stopifnot(file.exists(atlas_path))

# 1) Read fresh (no CRS changes yet, no make_valid yet)
fire_atlas_raw <- st_read(atlas_path, quiet = TRUE)

# 2) Inspect coordinate magnitudes to decide if coords are actually in meters
#    If |x| or |y| are >> 200, it’s not degrees -> it's projected meters.
#coords_sample <- st_coordinates(st_geometry(fire_atlas_raw)[[1]])  # sample first feature
#summary(abs(coords_sample))  # look at ranges

# 3) If the numbers are in the thousands/millions, assign the most likely CRS:
#    For US fire perimeters, CONUS Albers is very common: EPSG:5070
fire_atlas_fixed <- st_set_crs(fire_atlas_raw, 5070)  # assign (do NOT transform)

# 4) Now repair geometry and ensure polygons
fire_atlas_fixed <- fire_atlas_fixed %>%
  st_make_valid() %>%
  st_collection_extract("POLYGON")

# 5) Plot perimeters only, hollow purple outlines, with extent = layer bbox
bb <- st_bbox(fire_atlas_fixed)

#ggplot() +
#  geom_sf(data = fire_atlas_fixed, fill = NA, color = "#3f007d", linewidth = 0.5) +
#  coord_sf(
#    xlim = c(bb["xmin"], bb["xmax"]),
#    ylim = c(bb["ymin"], bb["ymax"]),
#    expand = FALSE
#  ) +
#  labs(title = "Fire Atlas Perimeters (fixed CRS, hollow outlines)") +
#  theme_minimal(base_size = 12)

## Add the pixels to the map:
raw_data_sf_5070 <- st_transform(raw_data_sf, st_crs(fire_atlas_fixed))

bb <- st_bbox(fire_atlas_fixed)
#ggplot() +
#  geom_sf(data = fire_atlas_fixed, fill = NA, color = "#3f007d", linewidth = 0.5) +
#  geom_sf(data = raw_data_sf_5070, aes(color = sev_class), size = 0.6, alpha = 0.9) +
#  scale_color_manual(values = sev_colors, limits = names(sev_colors), drop = FALSE) +
#  coord_sf(
#    xlim = c(bb["xmin"], bb["xmax"]),
#    ylim = c(bb["ymin"], bb["ymax"]),
#    expand = FALSE
#  ) +
#  labs(
#    title = "Pixels by Fire Severity with Fire Perimeters (aligned CRS)",
#    color = "Severity"
#  ) +
#  theme_minimal(base_size = 12)

## transforming to raw_pixel CRS instead:
fire_atlas_wgs84 <- st_transform(fire_atlas_fixed, 4326)
bb_ll <- st_bbox(fire_atlas_wgs84)

ggplot() +
  geom_sf(data = fire_atlas_wgs84, fill = NA, color = "#3f007d", linewidth = 0.5) +
  geom_sf(data = raw_data_sf, aes(color = sev_class), size = 0.6, alpha = 0.9) +
  scale_color_manual(values = sev_colors, limits = names(sev_colors), drop = FALSE) +
  coord_sf(xlim = c(bb_ll["xmin"], bb_ll["xmax"]),
           ylim = c(bb_ll["ymin"], bb_ll["ymax"]),
           expand = FALSE) +
  labs(title = "Pixels + Fire Perimeters (both in EPSG:4326)", color = "Severity") +
  theme_minimal(base_size = 12)
### ===================================================================== ######

## === PART 3: LOCATE / ADDRESS FIRE DISCREPANCIES

# Identify pixels where sev is there, but no fire atlas info
sev_only <- raw_data %>%
  filter(
    (num_of_burns == 0 | fire_name == "None") &  # No fire info
      sev_class %in% c("Unburned/Very Low", "Low", "Moderate", "High")   # Must have severity classification
  )

print(nrow(sev_only))  # Number of invalid cases: 1336

## plot the sev_only pixels alongside fire perimeters:
sev_only_sf <- sev_only %>%
  filter(!is.na(.geo)) %>%
  mutate(geometry = geojsonsf::geojson_sfc(.geo)) %>%
  st_sf(crs = 4326)

bb_ll <- st_bbox(fire_atlas_wgs84)

#ggplot() +
#  geom_sf(data = fire_atlas_wgs84, fill = NA, color = "#3f007d", linewidth = 0.5) +
#  geom_sf(data = sev_only_sf, aes(color = sev_class), size = 0.6, alpha = 0.9) +
#  scale_color_manual(values = sev_colors, limits = names(sev_colors), drop = FALSE) +
#  coord_sf(xlim = c(bb_ll["xmin"], bb_ll["xmax"]),
#           ylim = c(bb_ll["ymin"], bb_ll["ymax"]),
#           expand = FALSE) +
#  labs(title = "Pixels w/ Sev but No Fire Info (EPSG:4326)", color = "Severity") +
#  theme_minimal(base_size = 12)
## These pixels appear to be on edges of fire perimeters and primarily low-severity.


##############  ASSIGNING FIRE INFO TO SEV_ONLY PIXELS:  #####

## Assigning fire info (conservatively) to sev_only pixels
library(sf)
library(dplyr)
library(stringr)
library(purrr)

# 1) Transform just these pixels to 5070 (meters) and give them IDs
sev_only_sf_5070 <- st_transform(sev_only_sf, st_crs(fire_atlas_fixed)) %>%
  mutate(.pid = dplyr::row_number())

# 2) Prepare atlas attributes
atlas_keep <- fire_atlas_fixed %>%
  mutate(
    Year = suppressWarnings(as.numeric(Year)),
    Year = ifelse(is.na(Year),
                  suppressWarnings(as.numeric(str_extract(as.character(Year), "(19|20)\\d{2}"))),
                  Year)
  ) %>%
  select(Year, FireNam, Size_ha, geometry)

# 3) Define conservative edge distance (meters)
edge_dist <- 15  # adjust as-needed

# 4) Candidate matches: polygons within 'edge_dist' of each pixel
# st_is_within_distance is fast and returns a list of indices
cand_idx <- st_is_within_distance(sev_only_sf_5070, atlas_keep, dist = edge_dist)

# Expand into a long table of candidates
edge_candidates <- purrr::map2_dfr(
  .x = seq_along(cand_idx),
  .y = cand_idx,
  .f = function(i, idxs) {
    if (length(idxs) == 0) return(NULL)
    d <- st_distance(sev_only_sf_5070[i, ], atlas_keep[idxs, ], by_element = FALSE)
    tibble(
      .pid       = sev_only_sf_5070$.pid[i],
      atlas_row  = idxs,
      Year       = atlas_keep$Year[idxs],
      FireNam    = atlas_keep$FireNam[idxs],
      Size_ha    = atlas_keep$Size_ha[idxs],
      dist_m     = as.numeric(d)
    )
  }
)

# Quick sanity checks
message("rows(edge_candidates) = ", nrow(edge_candidates))
message("distinct .pid with >=1 candidate = ", dplyr::n_distinct(edge_candidates$.pid), " of ", nrow(sev_only_sf_5070))
print(edge_candidates %>% count(.pid, name = "n_cands") %>% summary())

# keep only the single nearest perimeter per pixel before the “latest year” resolution:
nearest_one <- edge_candidates %>%
  group_by(.pid) %>%
  slice_min(dist_m, with_ties = FALSE) %>%
  ungroup()

assignments <- nearest_one %>%
  transmute(
    .pid,
    fire_year_edge   = as.integer(round(Year)),
    fire_name_edge   = FireNam,
    fire_size_ha_edge= suppressWarnings(as.numeric(Size_ha)),
    edge_dist_m      = dist_m,
    fire_assign_flag = "assigned_edge_nearest_latest"
  )

# Join onto sev_only (the edge set you’re fixing)
sev_only_sf_5070_assigned <- sev_only_sf_5070 %>%
  left_join(assignments, by = ".pid") %>%
  mutate(
    fire_assign_flag = ifelse(is.na(fire_year_edge), "no_nearby_perimeter", fire_assign_flag)
  )

# Quick QA
sev_only_sf_5070_assigned %>%
  st_drop_geometry() %>%
  count(fire_assign_flag, sort = TRUE) %>%
  print(n = Inf)

# ---- Remove unassigned edge pixels (outside buffer) from raw_data 
unassigned_ids <- sev_only_sf_5070_assigned %>%
  st_drop_geometry() %>%
  filter(fire_assign_flag == "no_nearby_perimeter") %>%
  pull(pixel_ID) %>%
  unique()

before_rows <- nrow(raw_data)
raw_data <- raw_data %>% filter(!(pixel_ID %in% unassigned_ids))
cat("Rows removed (outside buffer): ", before_rows - nrow(raw_data), "\n")


#### Apply sev_only fire assignments back into raw_data:

# 1) Pull the assigned ones only
edge_updates <- sev_only_sf_5070_assigned %>%
  st_drop_geometry() %>%
  filter(fire_assign_flag == "assigned_edge_nearest_latest") %>%
  select(pixel_ID, fire_year_edge, fire_name_edge, fire_size_ha_edge)

# 2) Record baseline for QA
before_unassigned <- raw_data %>%
  filter(num_of_burns %in% c(0,1),
         sev_class %in% c("Unburned/Very Low","Low","Moderate","High"),
         fire_name %in% c(NA, "None")) %>%
  nrow()

# 3) Put edited pixels back into raw_data 
# Merge conservatively (only fill where fire info is missing!)
raw_data <- raw_data %>%
  left_join(edge_updates, by = "pixel_ID") %>%  # join sev_only pixels to raw_data by pixel_ID
  mutate(
    fire_year    = ifelse(!is.na(fire_year_edge), as.integer(round(fire_year_edge)), fire_year),  # use 'round' in case yr not stored as an integer (will ignore something like 2003.0)
    
    fire_name    = ifelse(!is.na(fire_name_edge), fire_name_edge, fire_name),
    fire_size_ha = ifelse(!is.na(fire_size_ha_edge), as.numeric(fire_size_ha_edge), fire_size_ha)
  ) %>%
  select(-fire_year_edge, -fire_name_edge, -fire_size_ha_edge)  # remove temporary join columns

# 4) QA report:
after_unassigned <- raw_data %>%
  filter(num_of_burns %in% c(0,1),
         sev_class %in% c("Unburned/Very Low","Low","Moderate","High"),
         fire_name %in% c(NA, "None")) %>%
  nrow()

cat("Previously unassigned (sev_only):", before_unassigned, "\n")
cat("Remaining unassigned after edge pass:", after_unassigned, "\n")
cat("Newly assigned via edge pass:", before_unassigned - after_unassigned, "\n")


# Mark which pixels were newly assigned in this pass
newly_assigned_ids <- edge_updates$pixel_ID

## build spatial data frame with newly-assigned pixels/fires
raw_new_sf <- raw_data %>%
  filter(pixel_ID %in% newly_assigned_ids, !is.na(.geo)) %>%
  mutate(geometry = geojsonsf::geojson_sfc(.geo)) %>%
  st_sf(crs = 4326) %>%
  st_transform(st_crs(fire_atlas_fixed)) %>%
  mutate(fire_event = fire_name)

# Add a common "fire_event" label to the newly assigned pixels
raw_new_sf <- raw_new_sf %>%
  mutate(fire_event = fire_name)   

# Subset the atlas to only those newly assigned fires, and add the same label
fires_for_plot <- sort(unique(raw_new_sf$fire_event))
fire_atlas_new <- fire_atlas_fixed %>%
  filter(FireNam %in% fires_for_plot) %>%
  mutate(fire_event = FireNam)

# Make a discrete palette keyed to each fire_event
library(scales)
pal <- setNames(hue_pal()(length(fires_for_plot)), fires_for_plot)

# Plot: perimeters (hollow) + pixels, both colored by fire_event
bb <- sf::st_bbox(raw_new_sf)

#ggplot() +
#  geom_sf(data = fire_atlas_new,
#          aes(color = fire_event),
#          fill = NA, linewidth = 0.6) +
#  geom_sf(data = raw_new_sf,
#          aes(color = fire_event),
#          size = 0.9, alpha = 0.95) +
#  scale_color_manual(values = pal, name = "Assigned fire") +
#  coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
#           ylim = c(bb["ymin"], bb["ymax"]),
#           expand = FALSE) +
#  labs(title = "Newly Assigned Edge Pixels by Fire Event") +
#  theme_minimal(base_size = 12)


# How many pixels were assigned via edge pass?
n_ids <- length(unique(newly_assigned_ids))
#cat("Newly assigned pixel_IDs:", n_ids, "\n")  # Newly assigned pixel_IDs: 847 

# pull just the updated rows
assigned_rows <- raw_data %>%
  filter(pixel_ID %in% newly_assigned_ids) %>%
  select(pixel_ID, fire_name, fire_year, fire_size_ha)

# collapse atlas to unique fires (name ~ year/size) for joining
atlas_ref <- fire_atlas_fixed %>%
  st_drop_geometry() %>%
  transmute(
    FireNam = as.character(FireNam),
    Year    = suppressWarnings(as.integer(round(Year))),
    Size_ha = suppressWarnings(as.numeric(Size_ha))
  ) %>%
  distinct(FireNam, Year, Size_ha)

# JOIN AND COMPARE:
chk <- assigned_rows %>%
  left_join(atlas_ref, by = c("fire_name" = "FireNam")) %>%
  mutate(
    has_atlas_match = !is.na(Year) | !is.na(Size_ha),
    year_ok  = !is.na(fire_year)   & !is.na(Year)    & (fire_year == Year),
    # allow a little numeric wiggle room on size
    size_ok  = !is.na(fire_size_ha) & !is.na(Size_ha) & abs(fire_size_ha - Size_ha) <= 1e-6,
    size_ok_tol = !is.na(fire_size_ha) & !is.na(Size_ha) & abs(fire_size_ha - Size_ha) <= 0.5
  )

#cat("Rows in assigned_rows: ", nrow(assigned_rows), "\n")
#cat("Unique pixel_IDs:      ", n_distinct(assigned_rows$pixel_ID), "\n")
#cat("Have atlas match:      ", sum(chk$has_atlas_match), " / ", nrow(chk), "\n")
#cat("Year exact match:      ", sum(chk$year_ok, na.rm = TRUE),  "\n")
#cat("Size exact match:      ", sum(chk$size_ok, na.rm = TRUE),  "\n")
#cat("Size match (±0.5 ha):  ", sum(chk$size_ok_tol, na.rm = TRUE), "\n")

# show any mismatches (top 10)
#mismatches <- chk %>% filter(!year_ok | !size_ok_tol | !has_atlas_match)
#if (nrow(mismatches) > 0) {
#  cat("\n-- Examples of mismatches (up to 10) --\n")
#  print(
#    mismatches %>%
#      select(pixel_ID, fire_name, fire_year, fire_size_ha, Year, Size_ha) %>%
#      head(10)
#  )
#} else {
#  cat("\nAll newly assigned pixels match atlas year/size (within tolerance).\n")
#}

# QA!!
total_rows <- nrow(raw_data)

with_sev <- raw_data %>%
  filter(!is.na(sev_num))

# Version A (strict): treat "None" as missing fire_name
missing_strict <- with_sev %>%
  filter(is.na(fire_name) | fire_name == "None" |
           is.na(fire_year) |
           is.na(fire_size_ha))

# Version B (lenient): only NA treated as missing
missing_lenient <- with_sev %>%
  filter(is.na(fire_name) |
           is.na(fire_year) |
           is.na(fire_size_ha))

#cat("Total pixels remaining in raw_data: ", total_rows, "\n")  # 52166 
#cat("Pixels with non-NA sev_num:         ", nrow(with_sev), "\n\n") # 17860 

#cat("-- Strict check (fire_name cannot be NA or 'None') --\n")
#cat("Missing any of {fire_name, fire_year, fire_size_ha}: ",  # 0
#    nrow(missing_strict), "\n")

#cat("\n-- Lenient check (only NA considered missing) --\n")
#cat("Missing any of {fire_name, fire_year, fire_size_ha}: ",  # 0
#    nrow(missing_lenient), "\n")

## ===========================================================================

#######  QA FOR NUMBER OF BURNS LAYER:  #########

data_clean <- raw_data %>%
  mutate(
    num_of_burns_old = num_of_burns,
    num_of_burns     = ifelse(sev_num > 0 & !is.na(fire_name), 1L, num_of_burns)
  )

# Count how many actually changed
n_updated <- sum(data_clean$num_of_burns != data_clean$num_of_burns_old, na.rm = TRUE)

cat("Pixels updated (num_of_burns set to 1):", n_updated, "\n")  # 482 

### ==========================================  ###
# Check for missing COMID or huc12
missing_comid_or_huc12 <- data_clean %>%
  filter(is.na(COMID) | is.na(huc12))

# Report result
if (nrow(missing_comid_or_huc12) == 0) {
  cat("All pixels have valid COMID and huc12 values.\n")
} else {
  cat("Found", nrow(missing_comid_or_huc12), "pixels with missing COMID or huc12.\n")
  print(missing_comid_or_huc12)
}
#  "All pixels have valid COMID and huc12 values."

# Plot the cleaned up data with severity classification

# Convert to sf from .geo column
data_clean_sf <- data_clean %>%
  filter(!is.na(.geo)) %>%
  mutate(geometry = geojsonsf::geojson_sfc(.geo)) %>%
  st_as_sf(crs = 4326)  # WGS84 lat/long

# plot
#ggplot() +
#  geom_sf(data = data_clean_sf, aes(color = sev_class), size = 1) +
#  scale_color_manual(values = c(
#    "High" = "red",
#    "Moderate" = "orange",
#    "Low" = "yellow",
#    "Unburned/Very Low" = "gray",
#    "None" = "black"
#  ), na.value = "black") +
#  theme_minimal() +
#  ggtitle("Fire Severity Classification of Cleaned Pixels") +
#  theme(legend.title = element_text(size = 10))



##################################################
#################   PLOTTING  -  PRE-THINNING     #########################

n_pixels <- nrow(data_clean_sf)

## Plot using ggplot
ggplot(data = data_clean_sf) +
  geom_sf(aes(color = sev_class), size = 1.2) +
  scale_color_manual(values = sev_colors, na.value = "black") +
  labs(
    title = "Spatial Distribution of Pixels by Fire Severity (Pre-Thinning)",
    subtitle = paste("Sample size:", format(n_pixels, big.mark = ",")),
    color = "Severity Class"
  ) +
  theme_minimal()

################################################################################
#####     SUB-SAMPLING:  0.04 KM SPACING    (THINNING)  - - FIXED!!!!!!
################################################################################

library(nngeo)  # For nearest-neighbor filtering

# Transform to UTM for accurate distance calculations
utm_crs <- 32612  # UTM Zone 12N 
data_utm <- st_transform(data_clean_sf, crs = utm_crs)

# Get unique HUC12 watersheds (character names are fine)
huc_groups <- unique(data_utm$huc12)

# Initialize an empty list to store results
thinned_list <- list()

# Define spacing cutoff, in meters  <<-- CHANGED TO 40 m
min_dist <- 40

# --- FIXED thinning function ---
thin_points_sf <- function(sf_points, min_dist = min_dist) {
  stopifnot(inherits(sf_points, "sf"))
  if (!requireNamespace("units", quietly = TRUE)) {
    dist_arg <- min_dist
  } else {
    dist_arg <- tryCatch(units::set_units(min_dist, "m"), error = function(e) min_dist)
  }
  n <- nrow(sf_points)
  if (n < 2) return(sf_points)
  
  keep <- rep(TRUE, n)
  for (i in seq_len(n - 1)) {
    if (!keep[i]) next
    cand_idx <- (i + 1):n
    rel <- sf::st_is_within_distance(sf_points[i, ], sf_points[cand_idx, , drop = FALSE], dist = dist_arg)[[1]]
    if (length(rel)) {
      keep[cand_idx[rel]] <- FALSE
    }
  }
  sf_points[keep, ]
}

# Loop through each HUC12 and apply 40m thinning -----------
for (huc in huc_groups) {
  cat("Processing HUC12: Thinning", huc, "\n")  # Print progress
  
  # Subset data for this HUC12
  subset_data <- data_utm %>% dplyr::filter(huc12 == huc)
  
  # Skip thinning if there are too few points
  if (nrow(subset_data) < 2) {
    thinned_list[[huc]] <- subset_data
    next
  }
  
  # Perform thinning (40 meter spacing)  <<-- CHANGED
  thinned_subset <- thin_points_sf(subset_data, min_dist = min_dist)
  
  # Store results
  thinned_list[[huc]] <- thinned_subset
}

# Combine all results
data_thinned_sf <- do.call(rbind, thinned_list)

# Print summary
cat("Total pixels before thinning:", nrow(data_utm), "\n")
#cat("Total pixels after thinning:", nrow(data_thinned_sf), "\n")
cat("Total pixels removed:", nrow(data_utm) - nrow(data_thinned_sf), "\n")

# Plot the thinned pixels
ggplot() +
  geom_sf(data = data_thinned_sf, aes(color = sev_class), size = 1) +
  scale_color_manual(values = c(
    "High" = "red",
    "Moderate" = "orange",
    "Low" = "yellow",
    "Unburned/Very Low" = "gray",
    "None" = "black"
  ), na.value = "black") +
  theme_minimal() +
  ggtitle("Cleaned, Thinned Pixels (Minimum 0.04 km Apart)") +  # <<-- CHANGED
  theme(legend.title = element_text(size = 10))


#########################################################
#  QA check: 40m spacing -------------------------------   
#########################################################

library(RANN)   # fast NN on coordinates

# --- Set this to a single HUC name, or NULL for all 
#huc_filter <- "Bartlett Creek"
huc_filter <- NULL   # check the whole study area. 

# Safety: projected CRS in meters
stopifnot(st_crs(data_thinned_sf)$epsg == 32612)

# Apply optional HUC filter
if (is.null(huc_filter)) {
  data_thin_utm <- data_thinned_sf
  huc_label <- "ALL WATERSHEDS"
} else {
  data_thin_utm <- dplyr::filter(data_thinned_sf, huc12 == huc_filter)
  huc_label <- paste("Watershed:", huc_filter)
}

## Run QA on thinning step, and produce map & histogram of thinned dataset: ----
# If too few points, exit gracefully
if (nrow(data_thin_utm) < 2) {
  message("Not enough points to compute nearest neighbors for ", huc_label, ".")
} else {
  # Stable ID
  data_thin_utm <- data_thin_utm %>% mutate(pid = seq_len(n()))
  
  # --- Fast nearest neighbor via RANN on coordinates 
  coords <- sf::st_coordinates(data_thin_utm)         # meters (UTM)
  nn     <- RANN::nn2(coords, k = 2)                  # self + nearest other
  nn_idx <- nn$nn.idx[, 2]
  nn_m   <- nn$nn.dists[, 2]
  
  nn_tbl <- tibble::tibble(
    from_id = data_thin_utm$pid,
    to_id   = data_thin_utm$pid[nn_idx],
    nn_m    = as.numeric(nn_m)
  ) %>%
    dplyr::left_join(sf::st_drop_geometry(data_thin_utm) %>% dplyr::select(pid, huc12, sev_class),
                     by = c("from_id" = "pid"))
  
  # Violations (< 40 m) and geometry only for those   <<-- CHANGED
  viol_idx <- which(nn_tbl$nn_m < 40 & !is.na(nn_tbl$to_id))
  
  if (length(viol_idx) > 0) {
    build_line <- function(i) {
      sf::st_linestring(rbind(coords[i, ], coords[nn_idx[i], ]))
    }
    viol_lines_sfc <- sf::st_sfc(lapply(viol_idx, build_line), crs = sf::st_crs(data_thin_utm))
    
    nn_violate <- sf::st_as_sf(
      nn_tbl[viol_idx, ],
      geometry = viol_lines_sfc,
      crs = sf::st_crs(data_thin_utm)
    )
  } else {
    # empty sf with correct CRS if no violations
    nn_violate <- sf::st_as_sf(nn_tbl[integer(0), ], geometry = sf::st_sfc(crs = sf::st_crs(data_thin_utm)))
  }
  
  # Console summary
  cat(huc_label, "\n",
      "Thinned points:", nrow(data_thin_utm), "\n",
      "NN pairs < 40 m:", nrow(nn_violate), "\n",           
      if (nrow(nn_violate) > 0)
        paste0("Min offending distance (m): ", round(min(nn_violate$nn_m), 2), "\n")
      else "", sep = "")
  
  # 40 m buffers for a quick visual check  
  thin_buf40 <- sf::st_buffer(data_thin_utm, dist = 40)
  
  # Midpoints for labeling violating distances (only if any)
  if (nrow(nn_violate) > 0) {
    mid_pts <- nn_violate %>%
      sf::st_line_sample(sample = 0.5) %>%
      sf::st_cast("POINT") %>%
      sf::st_as_sf(crs = sf::st_crs(nn_violate)) %>%
      dplyr::mutate(nn_m = nn_violate$nn_m)
    
    # extract coordinates for text labels
    mid_xy <- sf::st_coordinates(mid_pts)
    mid_pts_df <- mid_pts %>%
      sf::st_drop_geometry() %>%
      dplyr::mutate(x = mid_xy[, 1], y = mid_xy[, 2])
  } else {
    mid_pts_df <- NULL
  }
  
  # Map: buffers (light blue), violations (red), points by severity
  p_map <- ggplot() +
    geom_sf(data = thin_buf40, fill = NA, color = "lightblue", linewidth = 0.3) +  # <<-- CHANGED
    geom_sf(data = nn_violate, color = "red", linewidth = 1) +
    geom_sf(data = data_thin_utm, aes(color = sev_class), size = 1) +
    { if (!is.null(mid_pts_df))
      geom_text(
        data = mid_pts_df,
        aes(x = x, y = y, label = round(nn_m, 1)),
        size = 2.8
      ) else NULL } +
    scale_color_manual(values = c(
      "High" = "red",
      "Moderate" = "orange",
      "Low" = "yellow",
      "Unburned/Very Low" = "gray",
      "None" = "black"
    ), na.value = "black") +
    coord_sf() +
    labs(
      title = paste("Nearest-Neighbor Spacing QA (violations only) —", huc_label),
      subtitle = paste0(
        "Points: ", nrow(data_thin_utm),
        " | Violations (<40 m): ", nrow(nn_violate)         # <<-- CHANGED
      ),
      color = "Severity Class"
    ) +
    theme_minimal()
  
  print(p_map)   # ensure the map renders
  
  # Histogram of NN distances for this watershed
  p_hist <- ggplot(nn_tbl, aes(x = nn_m)) +
    geom_histogram(binwidth = 10, fill = "gray70", color = "black") +
    labs(
      title = paste("Distribution of Nearest-Neighbor Distances —", huc_label),
      x = "Distance to Nearest Neighbor (m)", y = "Count"
    ) +
    theme_minimal()
  
  print(p_hist)  
}
# -----------------------------------------------------

### 3/2/26: RAN THINNING QA ON ALL WATERSHEDS; 33 POINTS IN VIOLATION. #########

## CONSOLE OUPUT:
#ALL WATERSHEDS
#Thinned points:27467
#NN pairs < 40 m:33
#Min offending distance (m): 30

# =============================================================================
# TARGETED FINAL FIX: remove only pixels violating 40m rule (global enforcement)
# =============================================================================

min_dist <- 40

# Ensure UTM meters
stopifnot(inherits(data_thinned_sf, "sf"))
if (sf::st_crs(data_thinned_sf)$epsg != 32612) {
  d_utm <- sf::st_transform(data_thinned_sf, 32612)
} else {
  d_utm <- data_thinned_sf
}

# Stable row id for bookkeeping
d_utm <- d_utm %>% dplyr::mutate(pid = dplyr::row_number())

repeat {
  coords <- sf::st_coordinates(d_utm)
  
  if (nrow(coords) < 2) break
  
  nn <- RANN::nn2(coords, k = 2)  # self + nearest neighbor
  to_id <- d_utm$pid[nn$nn.idx[, 2]]
  nn_m  <- as.numeric(nn$nn.dists[, 2])
  
  viol <- tibble::tibble(
    from_id = d_utm$pid,
    to_id   = to_id,
    nn_m    = nn_m
  ) %>%
    dplyr::filter(!is.na(to_id), nn_m < min_dist)
  
  if (nrow(viol) == 0) break
  
  # Make undirected edges (avoid A-B and B-A duplicates)
  edges <- viol %>%
    dplyr::transmute(
      a = pmin(from_id, to_id),
      b = pmax(from_id, to_id)
    ) %>%
    dplyr::distinct()
  
  # Greedy vertex cover-ish: drop the node with highest "violation degree"
  # (usually removes very few points)
  deg <- c(table(edges$a), table(edges$b))
  deg <- tapply(deg, names(deg), sum)  # combine counts for nodes appearing in both cols
  deg <- sort(deg, decreasing = TRUE)
  
  # Choose pid to drop: highest degree; tie-breaker drops larger pid (stable)
  max_deg <- deg[1]
  cands <- as.integer(names(deg[deg == max_deg]))
  drop_pid <- max(cands)
  
  d_utm <- d_utm %>% dplyr::filter(pid != drop_pid)
}

# Final QA (should be zero)
coords_final <- sf::st_coordinates(d_utm)
if (nrow(coords_final) >= 2) {
  nn_final <- RANN::nn2(coords_final, k = 2)
  nn_m_final <- as.numeric(nn_final$nn.dists[, 2])
  n_viol_final <- sum(nn_m_final < min_dist, na.rm = TRUE)
  
  cat("Final points:", nrow(d_utm), "\n")
  cat("Final NN pairs < ", min_dist, " m: ", n_viol_final, "\n", sep = "")
  cat("Min NN distance (m): ", round(min(nn_m_final, na.rm = TRUE), 2), "\n", sep = "")
  
  stopifnot(n_viol_final == 0)
}

# console output:
#Final points: 27451   ## 16 additional pixels removed!
#Final NN pairs < 40 m: 0
#Min NN distance (m): 42.43

# Drop bookkeeping id
data_thinned_global_utm <- d_utm %>% dplyr::select(-pid)
data_thinned_global_sf  <- sf::st_transform(data_thinned_global_utm, 4326)

### IMPORTANT TO RUN THIS: RE-NAME DATASET FOR DOWNSTREAM CONTINUITY ###########
## Replace the working dataset used downstream:
data_thinned_sf <- data_thinned_global_sf
################################################################################

################################################################################
##############    SPATIALLY-CONSTRAIN UNBURNED PIXELS   ########################

### Create sev_group and burn_status columns in data_final_sf
data_final_sf <- data_thinned_sf %>%     ## NOTE RENAME TO 'FINAL_SF'
  mutate(
    # Create sev_group column based on sev_class
    sev_group = case_when(
      sev_class %in% c("None", "Unburned/Very Low") ~ "Unburned",
      sev_class == "Low" ~ "Low",
      sev_class == "Moderate" ~ "Moderate",
      sev_class == "High" ~ "High",
      TRUE ~ NA_character_  # In case there's an unexpected value
    ),
    
    # Create burn_status column based on sev_group
    burn_status = case_when(
      sev_group == "Unburned" ~ "Unburned",
      sev_group %in% c("Low", "Moderate", "High") ~ "Burned",
      TRUE ~ NA_character_  # Handle missing or unexpected values
    )
  )

# Drop geometry column for a data frame
#data_final_df <- st_drop_geometry(data_final_sf)

## Starting dataset: data_final_sf (Post-thinning)  ~26K PIXELS

# Load required libraries
library(sf)
library(dplyr)
library(ggplot2)
library(purrr)
library(geojsonsf)
library(RColorBrewer)

# Set separate buffer distances
fire_perim_buffer_dist <- 500     # meters (fire perimeter context)
burned_pixel_buffer_dist <- 1000  # meters (proximity to burned pixels)

# Helper function for any missing fire geometries:
safe_intersects <- function(a, b) {
  if (inherits(a, "sfc")) {
    a <- st_as_sf(data.frame(geometry = a), sf_column_name = "geometry")
  }
  
  result <- tryCatch(
    st_intersects(a, b, sparse = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(result)) {
    return(rep(FALSE, nrow(a)))
  } else {
    return(rowSums(result) > 0)
  }
}

## Upload shapefiles for fire perimeters and huc12 boundaries

# Reproject fire atlas polygons → EPSG:32612 (UTM Zone 12N)
fire_perims <- st_transform(fire_atlas_fixed, 32612)

# Load HUC12 boundaries
huc12_polys <- st_read("C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/ArcGIS/BMWC_Fire_Data/HUC12_BMWA/HUC12_BobMarshall.shp")
# Reproject HUC12 polygons to UTM Zone 12N (meters)
huc12_polys <- st_transform(huc12_polys, 32612)

# Reproject data_final_sf to UTM 12N.
data_final_sf <- st_transform(data_final_sf, 32612)

# Check CRS of all three
st_crs(fire_perims)         # EPSG:32612 (UTM 12N), meters
st_crs(huc12_polys)         # EPSG:32612 (UTM 12N), meters
st_crs(data_final_sf)       # EPSG:32612 (UTM 12N), meters

# Calculate 'pools' for UNB pixel harvest
burned_pixels <- data_final_sf %>% filter(burn_status == "Burned")
unburned_pool <- data_final_sf %>% filter(burn_status == "Unburned")

# Initialize result list and log storage
final_results <- list()
skipped_log <- data.frame(huc12 = character(), fire_name = character(), reason = character(), stringsAsFactors = FALSE)

###### Assigning fire_year_control to reduced Unburned pixels and ref_year to all #####
# Loop over all HUC12s in dataset
for (huc in unique(data_final_sf$huc12)) {
  burned_huc <- burned_pixels %>% filter(huc12 == huc)
  unburned_huc <- unburned_pool %>% filter(huc12 == huc)
  fires_in_huc <- unique(burned_huc$fire_name)
  
  for (fire in fires_in_huc) {
    message(paste("Processing:", huc, "-", fire))  # Progress message
    
    burned_fire <- burned_huc %>% filter(fire_name == fire)
    
    fire_year_lookup <- burned_fire %>%
      distinct(fire_year) %>%
      pull(fire_year)
    
    if (length(fire_year_lookup) != 1 || is.na(fire_year_lookup) || fire_year_lookup == 0) {
      warning(paste("Could not resolve fire_year for:", fire, "- skipping."))
      skipped_log <- rbind(skipped_log, data.frame(huc12 = huc, fire_name = fire, reason = "Missing or invalid fire_year"))
      next
    }
    
    fire_poly <- fire_perims %>% filter(FireNam == fire)
    
    if (nrow(fire_poly) == 0 || any(is.na(st_is_valid(fire_poly))) || !all(st_is_valid(fire_poly))) {
      warning(paste("No valid fire perimeter found for:", fire, "— skipping."))
      skipped_log <- rbind(skipped_log, data.frame(huc12 = huc, fire_name = fire, reason = "Missing or invalid fire perimeter"))
      next
    }
    
    fire_poly <- st_union(fire_poly) %>% st_make_valid() %>% st_as_sf()
    fire_buffer <- st_buffer(fire_poly, dist = fire_perim_buffer_dist)
    
    # Define buffer around burned pixels for proximity constraint
    burned_fire_buffer <- st_buffer(st_union(burned_fire), dist = burned_pixel_buffer_dist)
    
    max_n <- burned_fire %>%
      filter(sev_group %in% c("Low", "Moderate", "High")) %>%
      group_by(sev_group) %>%
      summarise(n = n(), .groups = "drop") %>%
      pull(n) %>%
      max(na.rm = TRUE)
    
    # Identify unburned pixels that fall inside any fire perimeter (for all fires)
    unburned_within_any_fire <- unburned_huc %>%
      filter(safe_intersects(geometry, fire_perims))
    
    # Exclude unburned pixels that fall within other fire polygons, unless within this fire
    unburned_exclusive <- unburned_huc %>%
      filter(
        !pixel_ID %in% unburned_within_any_fire$pixel_ID |
          safe_intersects(geometry, fire_poly)
      )
    
    # Priority 1: unburned pixels that intersect *this* fire polygon AND are near burned pixels
    controls_in_poly <- unburned_exclusive %>%
      filter(
        safe_intersects(geometry, fire_poly),
        safe_intersects(geometry, burned_fire_buffer)
      )
    
    # Priority 2: unburned pixels in fire buffer AND near burned pixels
    controls_in_buffer <- unburned_exclusive %>%
      filter(
        safe_intersects(geometry, fire_buffer),
        safe_intersects(geometry, burned_fire_buffer)
      ) %>%
      filter(!pixel_ID %in% controls_in_poly$pixel_ID)
    
    # Combine and sample
    available_controls <- bind_rows(controls_in_poly, controls_in_buffer)
    if (nrow(available_controls) == 0) {
      warning(paste("No available unburned controls for:", fire, "- skipping."))
      skipped_log <- rbind(skipped_log, data.frame(huc12 = huc, fire_name = fire, reason = "No available controls"))
      next
    }
    
    sampled_unburned <- available_controls %>%
      slice_sample(n = min(nrow(available_controls), max_n)) %>%
      mutate(
        fire_name_control = fire,
        fire_year_control = fire_year_lookup
      )
    
    burned_fire <- burned_fire %>%
      mutate(
        fire_name_control = fire,
        fire_year_control = NA_integer_
      )
    
    final_results[[paste(huc, fire, sep = "_")]] <- bind_rows(burned_fire, sampled_unburned)
  }
}

#Warning messages: false spacing
#1: No available unburned controls for: RAILLEY_MOUNTAIN.2007 - skipping. 
#2: No available unburned controls for: NR.1995.89 - skipping. 
#3: No valid fire perimeter found for: CANNON_2023 — skipping. 
#4: No available unburned controls for: RAILLEY_MOUNTAIN.2007 - skipping. 
#5: No available unburned controls for: BIRK.2001 - skipping.

# Combine and calculate ref_year
data_sf_reduced <- bind_rows(final_results) %>%   ## RE-NAME: data_sf_reduced
  mutate(
    ref_year = case_when(
      burn_status == "Unburned" ~ fire_year_control,
      TRUE ~ fire_year
    )
  )
################################################################################


##########    DATA SUMMARIES:   ###########
###########################################

# Total number of pixels:
total_pixels <- nrow(data_sf_reduced)
cat("Total pixels:", total_pixels, "\n")
# Pixel count by severity group
data_sf_reduced %>%
  count(sev_group, name = "pixel_count") %>%
  print()

# Define severity level ordering
sev_order <- c("Unburned", "Low", "Moderate", "High")

# Build pxl count summary table for each watershed and fire event
huc12_summary <- data_sf_reduced %>%
  st_drop_geometry() %>%
  group_by(huc12, fire_name_control, fire_year_control, sev_group) %>%
  summarise(pixel_count = n(), .groups = "drop") %>%
  mutate(sev_group = factor(sev_group, levels = sev_order)) %>%
  arrange(huc12, fire_name_control, sev_group)

# Write to Excel
#write_xlsx(
#  huc12_summary,
#  path = "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/huc12_sevgroup_PixelCounts_Fixed40mSpacing.xlsx"
#)

# Optional: write skipped log to CSV 
# write.csv(skipped_log, "skipped_fire_log.csv", row.names = FALSE)

#######################################
########## PLOTTING LOOP ##############
#######################################

fire_colors <- c("red", "orange", "gold", "tomato", "darkorange")

# Ensure HUC names are trimmed and character type
huc12_polys$name <- trimws(as.character(huc12_polys$name))
data_sf_reduced$huc12 <- trimws(as.character(data_sf_reduced$huc12))

# Initialize storage
plots <- list()
successful_hucs <- c()
skipped_hucs <- c()

for (huc in unique(data_sf_reduced$huc12)) {
  message("Processing:", huc)
  
  try({
    # Validate that HUC name exists in huc12_polys
    if (!(huc %in% huc12_polys$name)) {
      warning(paste("HUC", huc, "not found in huc12_polys$name. Skipping."))
      skipped_hucs <- c(skipped_hucs, huc)
      next
    }
    
    # Subset relevant data and geometry
    huc_data <- data_sf_reduced %>% filter(huc12 == huc)
    huc_name <- huc
    huc_geom <- huc12_polys %>%
      filter(name == huc_name) %>%
      st_union() %>%
      st_make_valid()  # Removed st_collection_extract
    
    # Construct sf boundary object
    huc_boundary <- st_sf(
      geometry = st_sfc(huc_geom),
      crs = st_crs(huc12_polys),
      sf_column_name = "geometry"
    )
    
    # Validate geometry
    if (nrow(huc_boundary) == 0 || any(is.na(st_is_valid(huc_boundary))) || !all(st_is_valid(huc_boundary))) {
      warning(paste("Skipping", huc_name, "- invalid boundary."))
      skipped_hucs <- c(skipped_hucs, huc_name)
      next
    }
    
    # Extract relevant fire perimeters
    fires_this_huc <- unique(huc_data$fire_name_control)
    fire_layer <- fire_perims %>%
      filter(FireNam %in% fires_this_huc) %>%
      st_crop(huc_boundary)
    
    fire_fill_map <- setNames(fire_colors[1:length(fires_this_huc)], fires_this_huc)
    
    # Burned/unburned pixel separation
    unburned_data <- huc_data %>%
      filter(sev_group == "Unburned") %>%
      mutate(control_label = paste0("Unburned: ", fire_name_control))
    
    burned_data <- huc_data %>% filter(sev_group != "Unburned")
    
    # Color map for unburned control groups
    local_control_fires <- unique(unburned_data$control_label)
    local_control_colors <- setNames(
      RColorBrewer::brewer.pal(n = max(3, length(local_control_fires)), name = "Set2")[1:length(local_control_fires)],
      local_control_fires
    )
    
    # Summary table for annotation
    summary_table <- huc_data %>%
      st_drop_geometry() %>%
      group_by(fire_name_control, sev_group) %>%
      summarise(n = n(), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = sev_group, values_from = n, values_fill = 0)
    
    for (sev in c("Unburned", "Low", "Moderate", "High")) {
      if (!(sev %in% names(summary_table))) summary_table[[sev]] <- 0
    }
    
    summary_text <- summary_table %>%
      mutate(
        summary_text = paste0(fire_name_control, ": ",
                              "Unb=", Unburned, ", Low=", Low,
                              ", Mod=", Moderate, ", High=", High)
      ) %>%
      pull(summary_text) %>%
      paste(collapse = "\n")
    
    # Final plot
    plots[[huc_name]] <- ggplot() +
      geom_sf(data = fire_layer, aes(fill = FireNam), color = "black", alpha = 0.15, size = 0.3) +
      geom_sf(data = huc_boundary, color = "navy", fill = NA, size = 2.0, show.legend = FALSE) +
      geom_sf(data = burned_data, aes(color = sev_group), size = 2.2) +
      geom_sf(data = unburned_data, aes(color = control_label), size = 2.2) +
      scale_color_manual(
        name = "Pixel Classification",
        values = c(sev_colors, local_control_colors),
        breaks = c(names(sev_colors), names(local_control_colors))
      ) +
      scale_fill_manual(values = fire_fill_map, name = "Fire Perimeter") +
      annotate("text", x = Inf, y = -Inf, label = summary_text,
               hjust = 1.1, vjust = -0.2, size = 4, fontface = "italic") +
      coord_sf(
        xlim = st_bbox(huc_boundary)[c("xmin", "xmax")],
        ylim = st_bbox(huc_boundary)[c("ymin", "ymax")],
        expand = FALSE
      ) +
      theme_minimal() +
      labs(
        title = paste("Pixel Assignments:", huc_name),
        subtitle = "Prioritized Unburned Controls by Fire Event"
      )
    
    successful_hucs <- c(successful_hucs, huc_name)
  }, silent = FALSE)
}

# ===========================================
# POST-LOOP: Sort, Rename, and Save Maps
# ===========================================

# QUALITY CHECK:  Print all plots to screen using original names
#for (plot_name in names(plots)) {
#  print(plots[[plot_name]])
#}

# Define output directory
plots_output_dir <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/NDVI_BoxPlots_OnceBurned_08252025/Pre-Fire_and_ΔNDVI_Boxplots_Maps_byHUC12_OB+UNB_ReducedControlPixels_40m_spacing/MAPS_HUC12_PixelsBySevGroup_BufferReduction"

# Set max filename length 
max_name_length <- 32

# Optional log for truncated names
truncation_log <- list()

# Save all plots
for (i in seq_along(plots)) {
  
  plot_obj <- plots[[i]]
  full_name <- names(plots)[i]
  
  # Extract and abbreviate HUC12 portion only for filename
  huc12_clean <- gsub("^\\d+\\.\\s*", "", full_name)
  huc12_clean <- gsub("Middle Fork Flathead River", "MFFR", huc12_clean)
  huc12_clean <- gsub("South Fork Flathead River", "SFFR", huc12_clean)
  huc12_clean <- gsub("South Fork White River", "SFWR", huc12_clean)
  huc12_clean <- gsub("South Fork River", "SFR", huc12_clean)
  huc12_clean <- gsub("Lower", "Lwr", huc12_clean)
  huc12_clean <- gsub("Spotted Bear", "SpB", huc12_clean)
  huc12_clean <- gsub("\\s+", "", huc12_clean)
  
  # Truncate if necessary
  if (nchar(huc12_clean) > max_name_length) {
    truncation_log[[full_name]] <- huc12_clean
    huc12_clean <- substr(huc12_clean, 1, max_name_length)
  }
  
  # Define file path
  filename <- paste0("map_", huc12_clean, ".jpg")
  filepath <- file.path(plots_output_dir, filename)
  
  # Save plots to folder
  ggsave(
    filename = filepath,
    plot = plot_obj,
    width = 12,
    height = 10.2,
    dpi = 500
  )
}

################################################################################
######  ====  QA data check on TerraClimate variables:  ======== ###############

############################
###  TerraClimate CWD:
############################
# 1) Identify CWD columns precisely: "CWD_YYYYMM_TC"
cwd_cols <- grep("^CWD_\\d{6}_TC$", names(raw_data), value = TRUE)

# 2) If any of those columns are character, coerce "null"/"NA"/"NaN"/"" to NA, then numeric (fast lapply)
char_mask <- vapply(raw_data[cwd_cols], is.character, logical(1))
if (any(char_mask)) {
  raw_data[cwd_cols[char_mask]] <- lapply(raw_data[cwd_cols[char_mask]], function(x) {
    x2 <- tolower(trimws(x))
    x2[x2 %in% c("null","na","nan","")] <- NA_character_
    suppressWarnings(as.numeric(x2))
  })
}

# 3) Matrix view for very fast checks
cwd_mat <- as.matrix(raw_data[, cwd_cols, drop = FALSE])

# 4) Quick global flags
any_missing <- anyNA(cwd_mat)                         # NA or NaN
any_nan     <- any(is.nan(cwd_mat))
any_inf     <- any(is.infinite(cwd_mat))

cat("Any missing (NA/NaN)?", any_missing, "\n",
    "Any NaN?", any_nan, "\n",
    "Any Inf?", any_inf, "\n", sep = "")

# 5) Column-level summaries (counts)
col_na   <- colSums(is.na(cwd_mat))
col_nan  <- colSums(is.nan(cwd_mat))                  # subset of NA
col_inf  <- colSums(is.infinite(cwd_mat))
col_neg  <- colSums(cwd_mat < 0, na.rm = TRUE)        # CWD should be >= 0; negatives are suspicious

col_summary <- tibble(
  column  = cwd_cols,
  n_NA    = as.integer(col_na),
  n_NaN   = as.integer(col_nan),
  n_Inf   = as.integer(col_inf),
  n_neg   = as.integer(col_neg)
) %>%
  mutate(n_bad = n_NA + n_Inf + n_neg) %>%
  arrange(desc(n_bad))

print(head(col_summary, 20))  # worst 20 columns

# 6) Row-level summaries (per pixel)
row_na_cnt  <- rowSums(is.na(cwd_mat))
row_inf_cnt <- rowSums(is.infinite(cwd_mat))
row_neg_cnt <- rowSums(cwd_mat < 0, na.rm = TRUE)

row_summary <- tibble(
  pixel_ID      = raw_data$pixel_ID,
  n_NA          = as.integer(row_na_cnt),
  n_Inf         = as.integer(row_inf_cnt),
  n_neg         = as.integer(row_neg_cnt),
  any_problem   = (n_NA + n_Inf + n_neg) > 0
)

#cat("Pixels with any problem:", sum(row_summary$any_problem), "of", nrow(row_summary), "\n") # 0


############################
###  TerraClimate SWE:
############################

# 1) Identify SWE columns precisely: "SWE_YYYYMM_TC" and keep months 01–05
swe_cols_all <- grep("^SWE_\\d{6}_TC$", names(raw_data), value = TRUE)
swe_mm       <- substr(swe_cols_all, 9, 10)
swe_cols     <- swe_cols_all[swe_mm %in% sprintf("%02d", 1:5)]  # 01–05

# 2) If any of those columns are character, coerce "null"/"NA"/"NaN"/"" to NA, then numeric (fast lapply)
char_mask <- vapply(raw_data[swe_cols], is.character, logical(1))
if (any(char_mask)) {
  raw_data[swe_cols[char_mask]] <- lapply(raw_data[swe_cols[char_mask]], function(x) {
    x2 <- tolower(trimws(x))
    x2[x2 %in% c("null","na","nan","")] <- NA_character_
    suppressWarnings(as.numeric(x2))
  })
}

# 3) Matrix view for very fast checks
swe_mat <- as.matrix(raw_data[, swe_cols, drop = FALSE])

# 4) Quick global flags
any_missing <- anyNA(swe_mat)                         # NA or NaN
any_nan     <- any(is.nan(swe_mat))
any_inf     <- any(is.infinite(swe_mat))

cat("SWE — Any missing (NA/NaN)? ", any_missing, "\n",
    "SWE — Any NaN? ",              any_nan,     "\n",
    "SWE — Any Inf? ",              any_inf,     "\n", sep = "")

# 5) Column-level summaries (counts)
col_na   <- colSums(is.na(swe_mat))
col_nan  <- colSums(is.nan(swe_mat))                  # subset of NA
col_inf  <- colSums(is.infinite(swe_mat))
col_neg  <- colSums(swe_mat < 0, na.rm = TRUE)        # SWE should be >= 0; negatives are suspicious

col_summary <- tibble(
  column  = swe_cols,
  n_NA    = as.integer(col_na),
  n_NaN   = as.integer(col_nan),
  n_Inf   = as.integer(col_inf),
  n_neg   = as.integer(col_neg)
) %>%
  mutate(n_bad = n_NA + n_Inf + n_neg) %>%
  arrange(desc(n_bad))

#print(head(col_summary, 20))  # worst 20 columns

# 6) Row-level summaries (per pixel)
row_na_cnt  <- rowSums(is.na(swe_mat))
row_inf_cnt <- rowSums(is.infinite(swe_mat))
row_neg_cnt <- rowSums(swe_mat < 0, na.rm = TRUE)

row_summary <- tibble(
  pixel_ID      = raw_data$pixel_ID,
  n_NA          = as.integer(row_na_cnt),
  n_Inf         = as.integer(row_inf_cnt),
  n_neg         = as.integer(row_neg_cnt),
  any_problem   = (n_NA + n_Inf + n_neg) > 0
)

#cat("SWE — Pixels with any problem: ",
#    sum(row_summary$any_problem), " of ", nrow(row_summary), "\n", sep = "")  # 0


# 1) Identify columns exactly: "cwd_30yrAvg_TC_04" ... "_08"
cwd30_cols_expected <- paste0("cwd_30yrAvg_TC_", sprintf("%02d", 4:8))
cwd30_cols_present  <- intersect(cwd30_cols_expected, names(raw_data))
cwd30_cols_missing  <- setdiff(cwd30_cols_expected, cwd30_cols_present)

if (length(cwd30_cols_missing) > 0) {
  warning("Missing expected columns: ", paste(cwd30_cols_missing, collapse = ", "))
}
stopifnot(length(cwd30_cols_present) > 0)

# 2) If any of those columns are character, coerce "null"/"NA"/"NaN"/"" to NA, then numeric (fast lapply)
char_mask <- vapply(raw_data[cwd30_cols_present], is.character, logical(1))
if (any(char_mask)) {
  raw_data[cwd30_cols_present[char_mask]] <- lapply(raw_data[cwd30_cols_present[char_mask]], function(x) {
    x2 <- tolower(trimws(x))
    x2[x2 %in% c("null","na","nan","")] <- NA_character_
    suppressWarnings(as.numeric(x2))
  })
}

# 3) Matrix view for very fast checks
cwd30_mat <- as.matrix(raw_data[, cwd30_cols_present, drop = FALSE])

# 4) Quick global flags
any_missing <- anyNA(cwd30_mat)                         # NA or NaN
any_nan     <- any(is.nan(cwd30_mat))
any_inf     <- any(is.infinite(cwd30_mat))

cat("CWD 30yrAvg — Any missing (NA/NaN)? ", any_missing, "\n",
    "CWD 30yrAvg — Any NaN? "           , any_nan,     "\n",
    "CWD 30yrAvg — Any Inf? "           , any_inf,     "\n", sep = "")

# 5) Column-level summaries (counts)
col_na   <- colSums(is.na(cwd30_mat))
col_nan  <- colSums(is.nan(cwd30_mat))                  # subset of NA
col_inf  <- colSums(is.infinite(cwd30_mat))
col_neg  <- colSums(cwd30_mat < 0, na.rm = TRUE)        # CWD should be >= 0; negatives are suspicious

col_summary <- tibble(
  column  = cwd30_cols_present,
  n_NA    = as.integer(col_na),
  n_NaN   = as.integer(col_nan),
  n_Inf   = as.integer(col_inf),
  n_neg   = as.integer(col_neg)
) %>%
  mutate(n_bad = n_NA + n_Inf + n_neg) %>%
  arrange(desc(n_bad))

# print(col_summary)

# 6) Row-level summaries (per pixel)
row_na_cnt  <- rowSums(is.na(cwd30_mat))
row_inf_cnt <- rowSums(is.infinite(cwd30_mat))
row_neg_cnt <- rowSums(cwd30_mat < 0, na.rm = TRUE)

row_summary <- tibble(
  pixel_ID      = raw_data$pixel_ID,
  n_NA          = as.integer(row_na_cnt),
  n_Inf         = as.integer(row_inf_cnt),
  n_neg         = as.integer(row_neg_cnt),
  any_problem   = (n_NA + n_Inf + n_neg) > 0
)

cat("CWD 30yrAvg — Pixels with any problem: ",
    sum(row_summary$any_problem), " of ", nrow(row_summary), "\n", sep = "")




#############################################################################
## Compute cwd_3yr_zscore (monthly), cwd_5yr_zscore (monthly), 
##         veg climate index, ndvi_prefire_3yr_ stats, ndvi_prefire_5yr_avg 
#############################################################################

# 11.22.2025

library(dplyr)
library(purrr)
library(scales)  # for rescale()
library(sf)

# Make sure necessary columns exist
required_cols <- c("ref_year", "cwd_30yrAvg_TC_04", "cwd_30yrAvg_TC_05", "cwd_30yrAvg_TC_06", "cwd_30yrAvg_TC_07", "cwd_30yrAvg_TC_08", "hli", "twi", "huc12")
missing_cols <- setdiff(required_cols, colnames(data_sf_reduced))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# List of unique HUC12 groups to process
huc_groups <- unique(data_sf_reduced$huc12)

# Define the 5 monthly 30-year average CWD columns (Apr–Aug)
cwd_30yr_avg_cols <- paste0("cwd_30yrAvg_TC_", sprintf("%02d", 4:8))
aug_cols_86_15 <- paste0("CWD_", 1986:2015, "08_TC")

# Get latest CWD year-month (YYYYMM) dynamically from column names
latest_cwd_yearmo <- max(as.numeric(sub("CWD_(\\d{6})_TC", "\\1", 
                                        grep("^CWD_\\d{6}_TC$", colnames(data_sf_reduced), value = TRUE))))

# Initialize list to collect results
processed_list <- list()

# Loop through each HUC group
for (i in seq_along(huc_groups)) {
  huc <- huc_groups[i]
  cat("Processing:", huc, "(", i, "of", length(huc_groups), ")\n")
  
  # Subset to HUC12
  subset_df <- data_sf_reduced %>%
    filter(huc12 == huc)
  
  # --------------------------------------------------
  # 1) CWD metrics (rowwise) and 30-yr SD
  # --------------------------------------------------
  subset_df <- subset_df %>%
    rowwise() %>%
    mutate(
      ## ---- CWD Metrics ----
      # 3-year postfire August mean: ref_year+1, +2, +3
      cwd_3yr_postfire_avg_08 = {
        if (is.na(ref_year)) {
          NA_real_
        } else {
          yrs  <- (ref_year + 1):(ref_year + 3)
          cols <- paste0("CWD_", yrs, "08_TC")  ## August only
          vals <- sapply(cols, function(col) {
            if (col %in% names(pick(everything()))) {
              suppressWarnings(as.numeric(pick(everything())[[col]]))
            } else {
              NA_real_
            }
          })
          if (sum(!is.na(vals)) == 3) mean(vals, na.rm = TRUE) else NA_real_
        }
      },
      # 5-year postfire August mean: ref_year+1 .. ref_year+5
      cwd_5yr_postfire_avg_08 = {
        if (is.na(ref_year)) {
          NA_real_
        } else {
          yrs  <- (ref_year + 1):(ref_year + 5)
          cols <- paste0("CWD_", yrs, "08_TC")  ## August CWD only
          vals <- sapply(cols, function(col) {
            if (col %in% names(pick(everything()))) {
              suppressWarnings(as.numeric(pick(everything())[[col]]))
            } else {
              NA_real_
            }
          })
          if (sum(!is.na(vals)) == 5) mean(vals, na.rm = TRUE) else NA_real_
        }
      },
      # 5-year PREFIRE August mean: ref_year-5 .. ref_year-1
      cwd_5yr_prefire_avg_08 = {
        if (is.na(ref_year)) {
          NA_real_
        } else {
          yrs  <- (ref_year - 5):(ref_year - 1)
          cols <- paste0("CWD_", yrs, "08_TC")
          vals <- sapply(cols, function(col) {
            if (col %in% names(pick(everything()))) {
              suppressWarnings(as.numeric(pick(everything())[[col]]))
            } else {
              NA_real_
            }
          })
          if (sum(!is.na(vals)) == 5) mean(vals, na.rm = TRUE) else NA_real_
        }
      },
      # August-only 30-yr SD (1985–2015) from individual August columns
      cwd_30yr_sd_08 = {
        vals <- c_across(any_of(aug_cols_86_15))
        if (sum(!is.na(vals)) >= 2) sd(vals, na.rm = TRUE) else NA_real_
      }
    ) %>%
    ungroup()
  
  # --------------------------------------------------
  # 2) August CWD z-scores and veg climate index
  # --------------------------------------------------
  # Latest & earliest AUGUST years available in this HUC subset (from CWD_YYYY08_TC columns)
  latest_cwd_aug_year <- {
    aug_cols <- grep("^CWD_\\d{6}_TC$", colnames(subset_df), value = TRUE)
    if (length(aug_cols) == 0) NA_integer_ else {
      yyyymm   <- as.integer(sub("CWD_(\\d{6})_TC", "\\1", aug_cols))
      aug_only <- yyyymm[yyyymm %% 100 == 8L]
      if (length(aug_only) == 0) NA_integer_ else max(aug_only %/% 100)
    }
  }
  earliest_cwd_aug_year <- {
    aug_cols <- grep("^CWD_\\d{6}_TC$", colnames(subset_df), value = TRUE)
    if (length(aug_cols) == 0) NA_integer_ else {
      yyyymm   <- as.integer(sub("CWD_(\\d{6})_TC", "\\1", aug_cols))
      aug_only <- yyyymm[yyyymm %% 100 == 8L]
      if (length(aug_only) == 0) NA_integer_ else min(aug_only %/% 100)
    }
  }
  
  subset_df <- subset_df %>%
    mutate(
      cwd_3yr_zscore_08 = dplyr::case_when(
        is.na(ref_year) | is.na(latest_cwd_aug_year) ~ NA_real_,
        (latest_cwd_aug_year - ref_year) < 3 ~ NA_real_,
        is.na(cwd_3yr_postfire_avg_08) | is.na(cwd_30yrAvg_TC_08) | is.na(cwd_30yr_sd_08) ~ NA_real_,
        cwd_30yr_sd_08 == 0 ~ NA_real_,
        TRUE ~ (cwd_3yr_postfire_avg_08 - cwd_30yrAvg_TC_08) / cwd_30yr_sd_08
      ),
      cwd_5yr_zscore_08 = dplyr::case_when(
        is.na(ref_year) | is.na(latest_cwd_aug_year) ~ NA_real_,
        (latest_cwd_aug_year - ref_year) < 5 ~ NA_real_,
        is.na(cwd_5yr_postfire_avg_08) | is.na(cwd_30yrAvg_TC_08) | is.na(cwd_30yr_sd_08) ~ NA_real_,
        cwd_30yr_sd_08 == 0 ~ NA_real_,
        TRUE ~ (cwd_5yr_postfire_avg_08 - cwd_30yrAvg_TC_08) / cwd_30yr_sd_08
      ),
      # Prefire 5-yr August z-score
      cwd_5yr_prefire_zscore_08 = dplyr::case_when(
        is.na(ref_year) | is.na(earliest_cwd_aug_year) ~ NA_real_,
        (ref_year - earliest_cwd_aug_year) < 5 ~ NA_real_,  # need 5 prefire Augusts available
        is.na(cwd_5yr_prefire_avg_08) | is.na(cwd_30yrAvg_TC_08) | is.na(cwd_30yr_sd_08) ~ NA_real_,
        cwd_30yr_sd_08 == 0 ~ NA_real_,
        TRUE ~ (cwd_5yr_prefire_avg_08 - cwd_30yrAvg_TC_08) / cwd_30yr_sd_08
      ),
      ##  Compute veg climate index
      hli_scaled          = scales::rescale(hli, na.rm = TRUE),
      cwd_scaled_08       = scales::rescale(cwd_30yrAvg_TC_08, na.rm = TRUE),
      veg_climate_index_08 = (hli_scaled + cwd_scaled_08) / 2
    )
  
  # --------------------------------------------------
  # 3) NDVI prefire 3-year summaries (mean, min, range, median)
  #    Years used: ref_year - 3, ref_year - 2, ref_year - 1
  # --------------------------------------------------
  subset_df <- subset_df %>%
    rowwise() %>%
    mutate(
      ndvi_prefire_3yr_avg = {
        if (is.na(ref_year)) {
          NA_real_
        } else {
          # Uses NDVI from years: ref_year - 3, ref_year - 2, ref_year - 1
          yrs  <- (ref_year - 3):(ref_year - 1)
          cols <- paste0("NDVI_", yrs)
          vals <- sapply(cols, function(col) {
            if (col %in% names(pick(everything()))) {
              suppressWarnings(as.numeric(pick(everything())[[col]]))
            } else {
              NA_real_
            }
          })
          if (sum(!is.na(vals)) == 3) mean(vals, na.rm = TRUE) else NA_real_
        }
      },
      ndvi_prefire_3yr_min = {
        if (is.na(ref_year)) {
          NA_real_
        } else {
          yrs  <- (ref_year - 3):(ref_year - 1)
          cols <- paste0("NDVI_", yrs)
          vals <- sapply(cols, function(col) {
            if (col %in% names(pick(everything()))) {
              suppressWarnings(as.numeric(pick(everything())[[col]]))
            } else {
              NA_real_
            }
          })
          if (sum(!is.na(vals)) == 3) min(vals, na.rm = TRUE) else NA_real_
        }
      },
      ndvi_prefire_3yr_range = {
        if (is.na(ref_year)) {
          NA_real_
        } else {
          yrs  <- (ref_year - 3):(ref_year - 1)
          cols <- paste0("NDVI_", yrs)
          vals <- sapply(cols, function(col) {
            if (col %in% names(pick(everything()))) {
              suppressWarnings(as.numeric(pick(everything())[[col]]))
            } else {
              NA_real_
            }
          })
          if (sum(!is.na(vals)) == 3) {
            diff(range(vals, na.rm = TRUE))
          } else {
            NA_real_
          }
        }
      },
      ndvi_prefire_3yr_med = {
        if (is.na(ref_year)) {
          NA_real_
        } else {
          yrs  <- (ref_year - 3):(ref_year - 1)
          cols <- paste0("NDVI_", yrs)
          vals <- sapply(cols, function(col) {
            if (col %in% names(pick(everything()))) {
              suppressWarnings(as.numeric(pick(everything())[[col]]))
            } else {
              NA_real_
            }
          })
          if (sum(!is.na(vals)) == 3) stats::median(vals, na.rm = TRUE) else NA_real_
        }
      }
    ) %>%
    ungroup()
  
  # Store result
  processed_list[[i]] <- subset_df
}

# Combine all HUCs into final output
data_sf_reduced <- bind_rows(processed_list)

# Final check
cat("Final dataset dimensions:", nrow(data_sf_reduced), "rows ×", ncol(data_sf_reduced), "columns.\n")

################################################################################
######     QA CHECK   ##############
####################################

# Check for missing values in newly computed columns
cat("Missing values in key columns:\n")
cols_to_check <- c("ref_year", "ndvi_prefire_3yr_avg", "ndvi_prefire_3yr_min", 
                   "ndvi_prefire_3yr_med", "ndvi_prefire_3yr_range", "cwd_3yr_zscore_08", 
                   "cwd_5yr_zscore_08", "hli", "twi", "sev_group", "veg_climate_index_08")
na_summary <- sapply(data_sf_reduced[cols_to_check], function(x) sum(is.na(x)))
print(na_summary)

# Invalid or missing geometry
missing_geometry <- sum(is.na(st_is_valid(data_sf_reduced)))
cat("\nPixels with invalid or missing geometry:", missing_geometry, "\n")

# Identify missing type for each row
data_missing <- data_sf_reduced %>%
  mutate(
    missing_type = case_when(
      !st_is_valid(geometry) ~ "Invalid geometry",
      is.na(ref_year) ~ "ref_year",
      is.na(cwd_3yr_zscore_08) ~ "cwd_3yr_zscore_08",
      is.na(cwd_5yr_zscore_08) ~ "cwd_5yr_zscore_08",
      is.na(ndvi_prefire_3yr_avg) ~ "ndvi-prefire-avg",
      is.na(ndvi_prefire_3yr_med) ~ "ndvi-prefire-med",
      is.na(ndvi_prefire_3yr_range) ~ "ndvi_prefire_range",
      is.na(ndvi_prefire_3yr_min) ~ "ndvi-prefire-minimum",
      is.na(hli) ~ "hli",
      is.na(twi) ~ "twi",
      is.na(sev_group) ~ "sev_group",
      is.na(veg_climate_index_08) ~ "veg_climate_index_08",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(missing_type))

# Spatial plot of pixels with missing data, colored by missing type
ggplot() +
  geom_sf(data = data_sf_reduced, color = "grey85", size = 0.2) +  # background pixels
  geom_sf(data = data_missing, aes(color = missing_type), size = 2) +  # colored by issue
  scale_color_viridis_d(option = "C", end = 0.9, name = "Missing or Invalid") +
  labs(title = "Pixels with Missing or Invalid Covariates") +
  theme_minimal() +
  coord_sf()

## There are several pixels missing pre-fire and post-fire data because
##    they're from fires at the beginning/end of data collection period (eg 1985, 2020)

# Re-define key columns to check (exclude cwd_5yr_zscore)
cols_to_check <- c("ref_year", "cwd_3yr_zscore_08", "hli", "twi", "sev_group", 
                   "veg_climate_index_08", "ndvi_prefire_3yr_avg", "ndvi_prefire_3yr_range",
                   "ndvi_prefire_3yr_min", "ndvi_prefire_3yr_med")

# Remove pixels with any NA in key columns OR with invalid geometry
data_sf_clean <- data_sf_reduced %>%
  filter(if_all(all_of(cols_to_check), ~ !is.na(.)) & st_is_valid(geometry))

################################################################################
#########   ============   COMPUTING PRISM COVARIATES  ============  ###########

### ================ Precipitation Covariates  =================================

######## == Total summer precip (JJA): sum of June+July+Aug for each pixel ###########
# 1) Geometry-free copy
data_no_geo <- sf::st_drop_geometry(data_sf_clean)

# 2) Keep only years where all three monthly columns exist
years <- 1984:2024
# helper: safely get a column or an NA vector if it doesn't exist
get_col <- function(df, nm) if (nm %in% names(df)) df[[nm]] else rep(NA_real_, nrow(df))

# 3) Build list of JJA totals for each year
# Toggle:
# na_if_missing_month = TRUE  -> JJA is NA if any month is missing
# na_if_missing_month = FALSE -> JJA sums available months (partial sum)
make_JJA_totals <- function(df, years, na_if_missing_month = TRUE) {
  out <- lapply(years, function(y) {
    c6 <- paste0("ppt_", y, "06_PR")
    c7 <- paste0("ppt_", y, "07_PR")
    c8 <- paste0("ppt_", y, "08_PR")
    v6 <- get_col(df, c6); v7 <- get_col(df, c7); v8 <- get_col(df, c8)
    
    if (na_if_missing_month) {
      total <- v6 + v7 + v8
      total[is.na(v6) | is.na(v7) | is.na(v8)] <- NA_real_
      total
    } else {
      rowSums(cbind(v6, v7, v8), na.rm = TRUE)
    }
  })
  names(out) <- paste0("ppttot_JJA_", years, "_PR")
  out
}

new_vals_list <- make_JJA_totals(data_no_geo, years, na_if_missing_month = TRUE)

# Add JJA columns to both data frames
data_no_geo   <- dplyr::bind_cols(data_no_geo,   as.data.frame(new_vals_list))
data_sf_clean <- dplyr::bind_cols(data_sf_clean, as.data.frame(new_vals_list))

# quick check on original data frame
data_sf_clean |>
  dplyr::select(starts_with("ppttot_JJA_")) |>
  dplyr::glimpse()

###### Tot precip QA: test 10 random pixels (1984–2024) ######

library(dplyr)
library(sf)

set.seed(58)

# 1) Geometry-free copy
df <- data_sf_clean %>% sf::st_drop_geometry()

# 2) Try to find a pixel ID column (case-insensitive); falls back to none if not found
nm  <- names(df); low <- tolower(nm)
id_candidates <- c("pixel_id", "pixelid", "pixel_id_num", "pixel", "pixelkey", "pixel_key", "pix_id", "id")
id_col <- nm[match(id_candidates, low, nomatch = 0)[1]]
# If not found, id_col will be NA -> we'll just skip it in the select(any_of())

# 3) Build the full list of JJA total columns for 1984–2024
jja_cols <- paste0("ppttot_JJA_", 1984:2024, "_PR")

# 4) Select columns (pixel ID if present, plus all 41 JJA columns).
#    Also include a couple of context columns you may want to see (optional).
sel_cols <- c(id_col, "ref_year", "sev_group", "ppttot_JJA_burnYr", jja_cols)

# 5) Sample 10 rows and print
sample10 <- df %>%
  select(any_of(sel_cols)) %>%
  slice_sample(n = min(10L, nrow(.)))

glimpse(sample10)
print(sample10, n = nrow(sample10), width = Inf)

##### Plot avg summer precip over time w/ avg NDVI — Gorge Creek #####

library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)

# --- Geometry-free copy (so we can pivot fast) ---
df <- data_sf_clean %>% sf::st_drop_geometry()

# --- Find a HUC12 "name" column and filter to Gorge Creek ---
nm  <- names(df); low <- tolower(nm)
huc_name_col <- nm[match(c("huc12","huc12_name","huc_name","huc12label","huc12_label"), low, nomatch = 0)[1]]
stopifnot(!is.na(huc_name_col))

df_gc <- df %>% filter(.data[[huc_name_col]] == "Gorge Creek")
stopifnot(nrow(df_gc) > 0)

# --- Identify columns ---
ndvi_cols <- grep("^NDVI_\\d{4}$", names(df_gc), value = TRUE)
ppt_cols  <- grep("^ppttot_JJA_\\d{4}_PR$", names(df_gc), value = TRUE)

stopifnot(length(ndvi_cols) > 0, length(ppt_cols) > 0)

# --- Pivot NDVI (all NDVI_YYYY) ---
ndvi_long <- df_gc %>%
  select(all_of(ndvi_cols)) %>%
  pivot_longer(
    everything(),
    names_to = "year",
    names_pattern = "NDVI_(\\d{4})",
    values_to = "ndvi"
  ) %>%
  mutate(year = as.integer(year)) %>%
  filter(year >= 1984, year <= 2024)

# --- Pivot JJA precip (all ppttot_JJA_YYYY_PR) ---
ppt_long <- df_gc %>%
  select(all_of(ppt_cols)) %>%
  pivot_longer(
    everything(),
    names_to = "year",
    names_pattern = "ppttot_JJA_(\\d{4})_PR",
    values_to = "ppt"
  ) %>%
  mutate(year = as.integer(year)) %>%
  filter(year >= 1984, year <= 2024)

# --- Compute yearly stats across pixels (mean ± 1 SD) ---
ndvi_stats <- ndvi_long %>%
  group_by(year) %>%
  summarize(
    mean = mean(ndvi, na.rm = TRUE),
    sd   = sd(ndvi,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(metric = "NDVI")

ppt_stats <- ppt_long %>%
  group_by(year) %>%
  summarize(
    mean = mean(ppt, na.rm = TRUE),
    sd   = sd(ppt,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(metric = "Total Summer Precip (JJA)")

stats <- bind_rows(
  ndvi_stats %>% mutate(lower = mean - sd, upper = mean + sd),
  ppt_stats  %>% mutate(lower = mean - sd, upper = mean + sd)
)

# --- Normalize each metric to 0–1 so they share an axis ---
stats_norm <- stats %>%
  group_by(metric) %>%
  mutate(
    minv = min(mean, na.rm = TRUE),
    maxv = max(mean, na.rm = TRUE),
    rng  = ifelse(maxv > minv, maxv - minv, NA_real_),
    mean_norm  = (mean  - minv) / rng,
    lower_norm = (lower - minv) / rng,
    upper_norm = (upper - minv) / rng
  ) %>%
  ungroup()

# --- Plot both series on same normalized axis with ±1 SD ribbons ---
ggplot(stats_norm, aes(x = year, y = mean_norm, color = metric, fill = metric)) +
  geom_ribbon(aes(ymin = lower_norm, ymax = upper_norm), alpha = 0.20, color = NA) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = seq(1984, 2024, by = 4)) +
  labs(
    title    = "Gorge Creek: NDVI and Total Summer Precip (JJA), 1984–2024",
    subtitle = "Means across pixels in HUC12; shaded band = ±1 SD",
    x = "Year",
    y = "Normalized units (0–1)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.title = element_blank())
################################################################################

######## == Precip 30-yr means and seasonal z-scores, for each year ##########
library(dplyr)
library(sf)

# 1) Geometry-free copy
data_no_geo <- sf::st_drop_geometry(data_sf_clean)

years_all <- 1984:2024
years_30  <- 1986:2015           # 30-yr baseline
mm_codes  <- c("06","07","08")   # JJA months

# Helper: numeric matrix from selected columns
mk_mat <- function(d, cols) {
  if (length(cols) == 0) return(matrix(NA_real_, nrow(d), 0))
  as.matrix(as.data.frame(lapply(d[, cols, drop = FALSE], as.numeric)))
}

# Helper: row-wise SD with NA handling
row_sds <- function(m) {
  mu <- rowMeans(m, na.rm = TRUE)
  sqrt(rowMeans((m - mu)^2, na.rm = TRUE))
}

# Ensure per-year JJA totals exist (ppttot_JJA_YYYY_PR); create/overwrite if missing
jja_tot_names <- paste0("ppttot_JJA_", years_all, "_PR")
have_jja_tot  <- jja_tot_names %in% names(data_no_geo)

if (!all(have_jja_tot)) {
  jja_tot_vals <- lapply(years_all, function(y) {
    cols_y <- paste0("ppt_", y, mm_codes, "_PR")
    cols_y <- intersect(cols_y, names(data_no_geo))
    if (length(cols_y) == 0) return(rep(NA_real_, nrow(data_no_geo)))
    # Sum J+J+A per pixel for that year (set na.rm=TRUE so partial months still contribute)
    rowSums(mk_mat(data_no_geo, cols_y), na.rm = TRUE)
  })
  names(jja_tot_vals) <- jja_tot_names
  
  data_no_geo   <- data_no_geo   %>% mutate( !!!jja_tot_vals )
  data_sf_clean <- data_sf_clean %>% mutate( !!!jja_tot_vals )
}

# Matrices of seasonal totals (all years, and baseline years)
jja_cols_all <- intersect(paste0("ppttot_JJA_", years_all, "_PR"), names(data_no_geo))
jja_cols_30  <- intersect(paste0("ppttot_JJA_", years_30,  "_PR"), names(data_no_geo))

M_all <- mk_mat(data_no_geo, jja_cols_all)
M_30  <- mk_mat(data_no_geo, jja_cols_30)

# 30-yr JJA mean & SD per pixel (1986–2015)
mu_30 <- rowMeans(M_30, na.rm = TRUE)
sd_30 <- row_sds(M_30)
sd_30_safe <- ifelse(sd_30 > 0, sd_30, NA_real_)  # avoid divide-by-zero

# Attach to both data frames
data_no_geo$pptavg_30yr_JJA_PR <- mu_30
data_no_geo$pptsd_30yr_JJA_PR  <- sd_30
data_sf_clean$pptavg_30yr_JJA_PR <- mu_30
data_sf_clean$pptsd_30yr_JJA_PR  <- sd_30

# Per-year JJA z-scores for all available years: [(seasonal total - 30yrmean) / 30yrsd]
if (ncol(M_all) > 0) {
  Z_all   <- (M_all - mu_30) / sd_30_safe
  z_names <- sub("^ppttot", "pptz", jja_cols_all)  # pptz_JJA_YYYY_PR
  colnames(Z_all) <- z_names
  Zdf <- as.data.frame(Z_all)
  
  # Write back (create or overwrite)
  for (nm in z_names) {
    data_no_geo[[nm]]   <- Zdf[[nm]]
    data_sf_clean[[nm]] <- Zdf[[nm]]
  }
}

# Quick check
data_sf_clean %>%
  dplyr::select(
    starts_with("ppttot_JJA_"),   # seasonal totals per year
    pptavg_30yr_JJA_PR, pptsd_30yr_JJA_PR,
    starts_with("pptz_JJA_")      # seasonal z per year
  ) %>%
  glimpse()



###### Precip Seasonal Z-score QA: Test print 30 random pixels ######
set.seed(50)

# Geometry-free copy
df <- data_sf_clean %>% sf::st_drop_geometry()

years_all <- 1984:2024

# Column names for per-year seasonal totals and z-scores
tot_cols_all <- paste0("ppttot_JJA_", years_all, "_PR")
z_cols_all   <- paste0("pptz_JJA_",    years_all, "_PR")

# Keep only those that actually exist in the data
tot_cols <- intersect(tot_cols_all, names(df))
z_cols   <- intersect(z_cols_all,   names(df))

# (Optional) interleave totals and z by year for easier side-by-side reading
yrs_both   <- years_all[
  paste0("ppttot_JJA_", years_all, "_PR") %in% tot_cols &
    paste0("pptz_JJA_",    years_all, "_PR") %in% z_cols
]
pair_cols  <- as.vector(rbind(paste0("ppttot_JJA_", yrs_both, "_PR"),
                              paste0("pptz_JJA_",    yrs_both, "_PR")))
# Add any leftover totals/z columns that may not have a partner
leftovers  <- setdiff(c(tot_cols, z_cols), pair_cols)

# Columns to display
sel_cols <- c(
  "pixel_ID",                    # pixel identifier
  "huc12",                       # optional: basin/grouping context
  "pptavg_30yr_JJA_PR",          # 30-yr mean (1986–2015) for this pixel
  "pptsd_30yr_JJA_PR",           # 30-yr sd   (1986–2015) for this pixel
  pair_cols, leftovers           # per-year totals & z-scores, interleaved where possible
)

# Sample 10 rows
n_take <- min(10L, nrow(df))
test10 <- df %>%
  select(any_of(sel_cols)) %>%
  slice_sample(n = n_take)

# Print
glimpse(test10)
print(test10, n = n_take, width = Inf)
################################################################################


### ================ Temperature Covariates  ===================================

###### ======= Tmax ==== 
data_no_geo <- sf::st_drop_geometry(data_sf_clean)
years_all <- 1984:2024
years_30  <- 1986:2015
mm_jja    <- c("06","07","08")  # June, July, August

# Helper: safe matrix builder
mk_mat <- function(d, cols) {
  if (length(cols) == 0) return(matrix(NA_real_, nrow(d), 0))
  as.matrix(as.data.frame(lapply(d[, cols, drop = FALSE], as.numeric)))
}

# 1) Per-year seasonal JJA MEAN tmax for every year (1984–2024), strict:
#    - If ANY of the three months are NA for a pixel -> NA
#    - If ANY of the three monthly columns are missing for that year -> NA for all pixels
tmaxmean_names <- paste0("tmaxmean_JJA_", years_all, "_PR")
tmaxmean_vals  <- lapply(years_all, function(y) {
  cols <- paste0("tmax_", y, mm_jja, "_PR")  # tmax_YYYY06_PR, _07, _08
  # require that all three monthly columns exist for this year
  if (!all(cols %in% names(data_no_geo))) {
    return(rep(NA_real_, nrow(data_no_geo)))
  }
  mat <- mk_mat(data_no_geo, cols)          # nrow x 3
  # strict mean: NA if any month is NA
  rowMeans(mat, na.rm = FALSE)
})

# Add into both data frames
data_no_geo   <- data_no_geo   %>% dplyr::mutate( !!!setNames(tmaxmean_vals, tmaxmean_names) )
data_sf_clean <- data_sf_clean %>% dplyr::mutate( !!!setNames(tmaxmean_vals, tmaxmean_names) )

# 2) Build 30-yr baseline (1986–2015) from these strict seasonal means
base_cols <- intersect(paste0("tmaxmean_JJA_", years_30, "_PR"), names(data_no_geo))
base_mat  <- mk_mat(data_no_geo, base_cols)

tmaxavg_30yr <- if (ncol(base_mat) > 0) rowMeans(base_mat, na.rm = TRUE) else rep(NA_real_, nrow(data_no_geo))
tmaxsd_30yr  <- if (ncol(base_mat) > 0) apply(base_mat, 1, sd, na.rm = TRUE) else rep(NA_real_, nrow(data_no_geo))

# Add 30-year metrics
data_no_geo   <- data_no_geo   %>% dplyr::mutate(tmaxavg_30yr_JJA_PR = tmaxavg_30yr,
                                                 tmaxsd_30yr_JJA_PR  = tmaxsd_30yr)
data_sf_clean <- data_sf_clean %>% dplyr::mutate(tmaxavg_30yr_JJA_PR = tmaxavg_30yr,
                                                 tmaxsd_30yr_JJA_PR  = tmaxsd_30yr)

# 3) Per-year seasonal z-scores for 1984–2024
tmaxz_names <- paste0("tmaxz_JJA_", years_all, "_PR")
tmaxz_vals  <- lapply(years_all, function(y) {
  v <- data_no_geo[[paste0("tmaxmean_JJA_", y, "_PR")]]
  z <- (v - tmaxavg_30yr) / tmaxsd_30yr
  z[!is.finite(z)] <- NA_real_  # guard against sd=0 or NA baseline
  z
})

data_no_geo   <- data_no_geo   %>% dplyr::mutate( !!!setNames(tmaxz_vals, tmaxz_names) )
data_sf_clean <- data_sf_clean %>% dplyr::mutate( !!!setNames(tmaxz_vals, tmaxz_names) )

# quick check on original data frame
data_sf_clean |>
  dplyr::select(
    dplyr::starts_with("tmaxmean_JJA_"),
    dplyr::starts_with("tmaxz_JJA_"),
    tmaxavg_30yr_JJA_PR,
    tmaxsd_30yr_JJA_PR
  ) |>
  dplyr::glimpse()


###### Tmax QA: Test 10 random pixels ######
set.seed(58)
df <- data_sf_clean %>% sf::st_drop_geometry()

tmaxmean_cols <- intersect(paste0("tmaxmean_JJA_", years_all, "_PR"), names(df))
tmaxz_cols    <- intersect(paste0("tmaxz_JJA_",    years_all, "_PR"), names(df))
clim_cols     <- c("tmaxavg_30yr_JJA_PR", "tmaxsd_30yr_JJA_PR")
clim_cols     <- clim_cols[clim_cols %in% names(df)]

sel_cols <- c("pixel_ID", clim_cols, tmaxmean_cols, tmaxz_cols)

n_take <- min(10L, nrow(df))
test10_tmax <- df %>%
  dplyr::select(dplyr::any_of(sel_cols)) %>%
  dplyr::slice_sample(n = n_take)

glimpse(test10_tmax)
print(test10_tmax, n = n_take, width = Inf)


##### Plot NDVI with new Tmax metrics — Gorge Creek #####
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)

# Geometry-free copy
df <- data_sf_clean %>% sf::st_drop_geometry()

huc_target <- "Gorge Creek"
df_gc <- df %>% filter(huc12 == huc_target)
stopifnot(nrow(df_gc) > 0)

years <- 1984:2024

# --- NDVI (annual) ---
ndvi_cols <- grep("^NDVI_\\d{4}$", names(df_gc), value = TRUE)
ndvi_long <- df_gc %>%
  select(all_of(ndvi_cols)) %>%
  pivot_longer(
    everything(),
    names_to      = "year",
    names_pattern = "NDVI_(\\d{4})",
    values_to     = "ndvi"
  ) %>%
  mutate(year = as.integer(year)) %>%
  filter(year >= 1984, year <= 2024)

ndvi_stats <- ndvi_long %>%
  group_by(year) %>%
  summarize(
    mean = mean(ndvi, na.rm = TRUE),
    sd   = sd(ndvi,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(metric = "NDVI")

# --- JJA Tmax (seasonal mean) ---
tmaxmean_cols <- grep("^tmaxmean_JJA_\\d{4}_PR$", names(df_gc), value = TRUE)
tmaxmean_long <- df_gc %>%
  select(all_of(tmaxmean_cols)) %>%
  pivot_longer(
    everything(),
    names_to      = "year",
    names_pattern = "tmaxmean_JJA_(\\d{4})_PR",
    values_to     = "tmaxmean_jja"
  ) %>%
  mutate(year = as.integer(year)) %>%
  filter(year >= 1984, year <= 2024)

tmaxmean_stats <- tmaxmean_long %>%
  group_by(year) %>%
  summarize(
    mean = mean(tmaxmean_jja, na.rm = TRUE),
    sd   = sd(tmaxmean_jja,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(metric = "JJA Tmax (seasonal mean)")

# --- JJA Tmax z-scores (native z-scale) ---
tmaxz_cols <- grep("^tmaxz_JJA_\\d{4}_PR$", names(df_gc), value = TRUE)
tmaxz_long <- df_gc %>%
  select(all_of(tmaxz_cols)) %>%
  pivot_longer(
    everything(),
    names_to      = "year",
    names_pattern = "tmaxz_JJA_(\\d{4})_PR",
    values_to     = "tmaxz_jja"
  ) %>%
  mutate(year = as.integer(year)) %>%
  filter(year >= 1984, year <= 2024)

tmaxz_stats <- tmaxz_long %>%
  group_by(year) %>%
  summarize(
    mean = mean(tmaxz_jja, na.rm = TRUE),
    sd   = sd(tmaxz_jja,   na.rm = TRUE),
    .groups = "drop"
  )

# --- Plot A: NDVI + JJA Tmax (seasonal mean) on a normalized axis ---
stats_A <- bind_rows(
  ndvi_stats,
  tmaxmean_stats
) %>%
  mutate(lower = mean - sd, upper = mean + sd)

stats_A_norm <- stats_A %>%
  group_by(metric) %>%
  mutate(
    minv = min(mean, na.rm = TRUE),
    maxv = max(mean, na.rm = TRUE),
    rng  = ifelse(maxv > minv, maxv - minv, NA_real_),
    mean_norm  = (mean  - minv) / rng,
    lower_norm = (lower - minv) / rng,
    upper_norm = (upper - minv) / rng
  ) %>%
  ungroup()

pA <- ggplot(stats_A_norm, aes(x = year, y = mean_norm, color = metric, fill = metric)) +
  geom_ribbon(aes(ymin = lower_norm, ymax = upper_norm), alpha = 0.20, color = NA) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = seq(1984, 2024, by = 4)) +
  labs(
    title    = paste0(huc_target, ": NDVI & JJA Tmax (seasonal mean), 1984–2024"),
    subtitle = "Means across all pixels within HUC12; ribbons are ±1 SD; each metric normalized to 0–1",
    x = "Year", y = "Normalized units (0–1)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.title = element_blank())

print(pA)

# --- Plot B: JJA Tmax z-scores (native z-scale) ---
pB <- ggplot(tmaxz_stats, aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.20) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = seq(1984, 2024, by = 4)) +
  labs(
    title    = paste0(huc_target, ": JJA Tmax z-score (vs 1986–2015), 1984–2024"),
    subtitle = "Means across all pixels within HUC12; ribbon is ±1 SD",
    x = "Year", y = "Z-score"
  ) +
  theme_minimal(base_size = 12)

print(pB)
################################################################################

###### ======= Tmean ==== 
data_no_geo <- sf::st_drop_geometry(data_sf_clean)
years_all <- 1984:2024
years_30  <- 1986:2015
mm_jja    <- c("06","07","08")  # June, July, August

# Helper: safe matrix builder
mk_mat <- function(d, cols) {
  if (length(cols) == 0) return(matrix(NA_real_, nrow(d), 0))
  as.matrix(as.data.frame(lapply(d[, cols, drop = FALSE], as.numeric)))
}

# 1) Per-year seasonal JJA MEAN tmean for every year (1984–2024), strict:
#    - If ANY of the three months are NA for a pixel -> NA
#    - If ANY of the three monthly columns are missing for that year -> NA for all pixels
tmeanavg_names <- paste0("tmeanavg_JJA_", years_all, "_PR")
tmeanavg_vals  <- lapply(years_all, function(y) {
  cols <- paste0("tmean_", y, mm_jja, "_PR")  # tmean_YYYY06_PR, _07, _08
  # require that all three monthly columns exist for this year
  if (!all(cols %in% names(data_no_geo))) {
    return(rep(NA_real_, nrow(data_no_geo)))
  }
  mat <- mk_mat(data_no_geo, cols)            # nrow x 3
  # strict mean: NA if any month is NA
  rowMeans(mat, na.rm = FALSE)
})

# Add into both data frames
data_no_geo   <- data_no_geo   %>% dplyr::mutate( !!!setNames(tmeanavg_vals, tmeanavg_names) )
data_sf_clean <- data_sf_clean %>% dplyr::mutate( !!!setNames(tmeanavg_vals, tmeanavg_names) )

# 2) Build 30-yr baseline (1986–2015) from these strict seasonal means
base_cols <- intersect(paste0("tmeanavg_JJA_", years_30, "_PR"), names(data_no_geo))
base_mat  <- mk_mat(data_no_geo, base_cols)

tmeanavg_30yr <- if (ncol(base_mat) > 0) rowMeans(base_mat, na.rm = TRUE) else rep(NA_real_, nrow(data_no_geo))
tmeansd_30yr  <- if (ncol(base_mat) > 0) apply(base_mat, 1, sd, na.rm = TRUE) else rep(NA_real_, nrow(data_no_geo))

# Add 30-year metrics
data_no_geo   <- data_no_geo   %>% dplyr::mutate(tmeanavg_30yr_JJA_PR = tmeanavg_30yr,
                                                 tmeansd_30yr_JJA_PR  = tmeansd_30yr)
data_sf_clean <- data_sf_clean %>% dplyr::mutate(tmeanavg_30yr_JJA_PR = tmeanavg_30yr,
                                                 tmeansd_30yr_JJA_PR  = tmeansd_30yr)

# 3) Per-year seasonal z-scores for 1984–2024
tmeanz_names <- paste0("tmeanz_JJA_", years_all, "_PR")
tmeanz_vals  <- lapply(years_all, function(y) {
  v <- data_no_geo[[paste0("tmeanavg_JJA_", y, "_PR")]]
  z <- (v - tmeanavg_30yr) / tmeansd_30yr
  z[!is.finite(z)] <- NA_real_  # guard against sd=0 or NA baseline
  z
})

data_no_geo   <- data_no_geo   %>% dplyr::mutate( !!!setNames(tmeanz_vals, tmeanz_names) )
data_sf_clean <- data_sf_clean %>% dplyr::mutate( !!!setNames(tmeanz_vals, tmeanz_names) )

# quick check on original data frame
data_sf_clean |>
  dplyr::select(
    dplyr::starts_with("tmeanavg_JJA_"),
    dplyr::starts_with("tmeanz_JJA_"),
    tmeanavg_30yr_JJA_PR,
    tmeansd_30yr_JJA_PR
  ) |>
  dplyr::glimpse()


###### Tmean QA: Test 10 random pixels ######
set.seed(58)
df <- data_sf_clean %>% sf::st_drop_geometry()

tmeanavg_cols <- intersect(paste0("tmeanavg_JJA_", years_all, "_PR"), names(df))
tmeanz_cols    <- intersect(paste0("tmeanz_JJA_",    years_all, "_PR"), names(df))
clim_cols      <- c("tmeanavg_30yr_JJA_PR", "tmeansd_30yr_JJA_PR")
clim_cols      <- clim_cols[clim_cols %in% names(df)]

sel_cols <- c("pixel_ID", clim_cols, tmeanavg_cols, tmeanz_cols)

n_take <- min(10L, nrow(df))
test10_tmean <- df %>%
  dplyr::select(dplyr::any_of(sel_cols)) %>%
  dplyr::slice_sample(n = n_take)

glimpse(test10_tmean)
print(test10_tmean, n = n_take, width = Inf)


##### Plot NDVI with new Tmean metrics — Gorge Creek #####
# Geometry-free copy
df <- data_sf_clean %>% sf::st_drop_geometry()

huc_target <- "Gorge Creek"
df_gc <- df %>% dplyr::filter(huc12 == huc_target)
stopifnot(nrow(df_gc) > 0)

years <- 1984:2024

# --- NDVI (annual) ---
ndvi_cols <- grep("^NDVI_\\d{4}$", names(df_gc), value = TRUE)
ndvi_long <- df_gc %>%
  dplyr::select(dplyr::all_of(ndvi_cols)) %>%
  tidyr::pivot_longer(
    everything(),
    names_to      = "year",
    names_pattern = "NDVI_(\\d{4})",
    values_to     = "ndvi"
  ) %>%
  dplyr::mutate(year = as.integer(year)) %>%
  dplyr::filter(year >= 1984, year <= 2024)

ndvi_stats <- ndvi_long %>%
  dplyr::group_by(year) %>%
  dplyr::summarize(
    mean = mean(ndvi, na.rm = TRUE),
    sd   = sd(ndvi,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(metric = "NDVI")

# --- JJA Tmean (seasonal mean) ---
tmeanavg_cols <- grep("^tmeanavg_JJA_\\d{4}_PR$", names(df_gc), value = TRUE)
tmeanavg_long <- df_gc %>%
  dplyr::select(dplyr::all_of(tmeanavg_cols)) %>%
  tidyr::pivot_longer(
    everything(),
    names_to      = "year",
    names_pattern = "tmeanavg_JJA_(\\d{4})_PR",
    values_to     = "tmeanavg_jja"
  ) %>%
  dplyr::mutate(year = as.integer(year)) %>%
  dplyr::filter(year >= 1984, year <= 2024)

tmeanavg_stats <- tmeanavg_long %>%
  dplyr::group_by(year) %>%
  dplyr::summarize(
    mean = mean(tmeanavg_jja, na.rm = TRUE),
    sd   = sd(tmeanavg_jja,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(metric = "JJA Tmean (seasonal mean)")

# --- JJA Tmean z-scores (native z-scale) ---
tmeanz_cols <- grep("^tmeanz_JJA_\\d{4}_PR$", names(df_gc), value = TRUE)
tmeanz_long <- df_gc %>%
  dplyr::select(dplyr::all_of(tmeanz_cols)) %>%
  tidyr::pivot_longer(
    everything(),
    names_to      = "year",
    names_pattern = "tmeanz_JJA_(\\d{4})_PR",
    values_to     = "tmeanz_jja"
  ) %>%
  dplyr::mutate(year = as.integer(year)) %>%
  dplyr::filter(year >= 1984, year <= 2024)

tmeanz_stats <- tmeanz_long %>%
  dplyr::group_by(year) %>%
  dplyr::summarize(
    mean = mean(tmeanz_jja, na.rm = TRUE),
    sd   = sd(tmeanz_jja,   na.rm = TRUE),
    .groups = "drop"
  )

# --- Plot A: NDVI + JJA Tmean (seasonal mean) on a normalized axis ---
stats_A <- dplyr::bind_rows(
  ndvi_stats,
  tmeanavg_stats
) %>%
  dplyr::mutate(lower = mean - sd, upper = mean + sd)

stats_A_norm <- stats_A %>%
  dplyr::group_by(metric) %>%
  dplyr::mutate(
    minv = min(mean, na.rm = TRUE),
    maxv = max(mean, na.rm = TRUE),
    rng  = ifelse(maxv > minv, maxv - minv, NA_real_),
    mean_norm  = (mean  - minv) / rng,
    lower_norm = (lower - minv) / rng,
    upper_norm = (upper - minv) / rng
  ) %>%
  dplyr::ungroup()

pA <- ggplot(stats_A_norm, aes(x = year, y = mean_norm, color = metric, fill = metric)) +
  geom_ribbon(aes(ymin = lower_norm, ymax = upper_norm), alpha = 0.20, color = NA) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = seq(1984, 2024, by = 4)) +
  labs(
    title    = paste0(huc_target, ": NDVI & JJA Tmean (seasonal mean), 1984–2024"),
    subtitle = "Means across all pixels within HUC12; ribbons are ±1 SD; each metric normalized to 0–1",
    x = "Year", y = "Normalized units (0–1)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.title = element_blank())

print(pA)

# --- Plot B: JJA Tmean z-scores (native z-scale) ---
pB <- ggplot(tmeanz_stats, aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.20) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = seq(1984, 2024, by = 4)) +
  labs(
    title    = paste0(huc_target, ": JJA Tmean z-score (vs 1986–2015), 1984–2024"),
    subtitle = "Means across all pixels within HUC12; ribbon is ±1 SD",
    x = "Year", y = "Z-score"
  ) +
  theme_minimal(base_size = 12)

print(pB)
################################################################################

################################################################################
#########   ============   COMPUTING TERRA COVARIATES  ============  ###########

###### ======= SWE ===== 
# Geometry-free copy for fast numeric ops
data_no_geo <- data_sf_clean %>% sf::st_drop_geometry()

years_all <- 1984:2024
mm_MAM    <- c("03","04","05")  # Mar, Apr, May

# Helper: numeric matrix from selected columns
mk_mat <- function(d, cols) {
  if (length(cols) == 0) return(matrix(NA_real_, nrow(d), 0))
  as.matrix(as.data.frame(lapply(d[, cols, drop = FALSE], as.numeric)))
}

# Build per-year peak SWE across M/A/M (NA if all three months missing for a row)
swepeak_names <- paste0("swe_peak_MAM_", years_all)
swepeak_vals  <- lapply(years_all, function(y) {
  cols_have <- intersect(paste0("SWE_", y, mm_MAM, "_TC"), names(data_no_geo))
  if (length(cols_have) == 0) return(rep(NA_real_, nrow(data_no_geo)))
  mat <- mk_mat(data_no_geo, cols_have)  # nrow x (1..3)
  out <- apply(mat, 1, function(r) if (all(is.na(r))) NA_real_ else max(r, na.rm = TRUE))
  out
})

# Add to both data frames (sf + geometry-free, optional)
data_no_geo   <- data_no_geo   %>% mutate( !!!setNames(swepeak_vals, swepeak_names) )
data_sf_clean <- data_sf_clean %>% mutate( !!!setNames(swepeak_vals, swepeak_names) )

# Quick check
data_sf_clean |>
  dplyr::select(starts_with("swe_peak_MAM_")) |>
  dplyr::glimpse()

###### SWE QA: Test 10 random pixels #######
set.seed(5)

# -------------------- Ensure per-year PEAK SWE columns exist 
df <- data_sf_clean %>% sf::st_drop_geometry()

years <- 1984:2024
mm    <- c("03","04","05")  # Mar/Apr/May

mk_mat <- function(d, cols) {
  if (length(cols) == 0) return(matrix(NA_real_, nrow(d), 0))
  as.matrix(as.data.frame(lapply(d[, cols, drop = FALSE], as.numeric)))
}

peak_names <- paste0("swe_peak_MAM_", years)
peak_vals  <- lapply(years, function(y) {
  cols      <- paste0("SWE_", y, mm, "_TC")
  cols_have <- intersect(cols, names(df))
  if (length(cols_have) == 0) return(rep(NA_real_, nrow(df)))
  M <- mk_mat(df, cols_have)
  # elementwise max over available months; NA if all three months are NA
  pk <- do.call(pmax, c(as.data.frame(M), na.rm = TRUE))
  pk[rowSums(!is.na(M)) == 0] <- NA_real_
  pk
})

# Add peak columns to the main sf data
data_sf_clean <- data_sf_clean %>%
  dplyr::mutate( !!!setNames(peak_vals, peak_names) )

# -------------------- QA print: 10 pixels x all years 
df <- data_sf_clean %>% sf::st_drop_geometry()

# Monthly SWE long -> wide (SWE_Mar/SWE_Apr/SWE_May per pixel-year)
swe_months_wide <- df %>%
  dplyr::select(pixel_ID, dplyr::matches("^SWE_\\d{4}(03|04|05)_TC$")) %>%
  tidyr::pivot_longer(
    -pixel_ID,
    names_to = c("year", "mm"),
    names_pattern = "SWE_(\\d{4})(\\d{2})_TC",
    values_to = "SWE_val"
  ) %>%
  dplyr::mutate(
    year = as.integer(year),
    mm   = dplyr::recode(mm, "03" = "SWE_Mar", "04" = "SWE_Apr", "05" = "SWE_May")
  ) %>%
  tidyr::pivot_wider(names_from = mm, values_from = SWE_val)

# Peak SWE long
swe_peak_long <- df %>%
  dplyr::select(pixel_ID, dplyr::matches("^swe_peak_MAM_\\d{4}$")) %>%
  tidyr::pivot_longer(
    -pixel_ID,
    names_to = "year",
    names_pattern = "swe_peak_MAM_(\\d{4})",
    values_to = "swe_peak_MAM"
  ) %>%
  dplyr::mutate(year = as.integer(year))

# Join monthly + peak into one table
swe_all <- swe_months_wide %>%
  dplyr::left_join(swe_peak_long, by = c("pixel_ID", "year")) %>%
  dplyr::arrange(pixel_ID, year)

# Sample 10 pixel IDs and print all years for those pixels
pix10 <- sample(unique(swe_all$pixel_ID), size = min(10L, dplyr::n_distinct(swe_all$pixel_ID)))
test10_swe <- swe_all %>% dplyr::filter(pixel_ID %in% pix10)

glimpse(test10)
#print(test10, n = nrow(test10), width = Inf)

##### Plot NDVI with April SWE & Peak SWE — Gorge Creek #####
# Geometry-free copy
df <- data_sf_clean %>% sf::st_drop_geometry()

# --- Filter to the requested HUC12 ---
huc_target <- "Gorge Creek"
df_gc <- df %>% filter(huc12 == huc_target)
stopifnot(nrow(df_gc) > 0)

years <- 1984:2024

# --- Pivot NDVI ---
ndvi_cols <- grep("^NDVI_\\d{4}$", names(df_gc), value = TRUE)
ndvi_long <- df_gc %>%
  select(all_of(ndvi_cols)) %>%
  pivot_longer(
    everything(),
    names_to      = "year",
    names_pattern = "NDVI_(\\d{4})",
    values_to     = "ndvi"
  ) %>%
  mutate(year = as.integer(year)) %>%
  filter(year >= 1984, year <= 2024)

ndvi_stats <- ndvi_long %>%
  group_by(year) %>%
  summarize(
    mean = mean(ndvi, na.rm = TRUE),
    sd   = sd(ndvi,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(metric = "NDVI")

# --- Pull April SWE (SWE_YYYY04_TC) ------
apr_cols <- grep("^SWE_\\d{4}04_TC$", names(df_gc), value = TRUE)
apr_long <- df_gc %>%
  select(all_of(apr_cols)) %>%
  pivot_longer(
    everything(),
    names_to      = "year",
    names_pattern = "SWE_(\\d{4})04_TC",
    values_to     = "swe_apr"
  ) %>%
  mutate(year = as.integer(year)) %>%
  filter(year >= 1984, year <= 2024)

apr_stats <- apr_long %>%
  group_by(year) %>%
  summarize(
    mean = mean(swe_apr, na.rm = TRUE),
    sd   = sd(swe_apr,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(metric = "SWE April")

# --- Compute Peak SWE (MAM): use precomputed swe_peak_MAM_YYYY if available; else, derive ----
peak_cols <- grep("^swe_peak_MAM_\\d{4}$", names(df_gc), value = TRUE)

if (length(peak_cols) > 0) {
  peak_long <- df_gc %>%
    select(all_of(peak_cols)) %>%
    pivot_longer(
      everything(),
      names_to      = "year",
      names_pattern = "swe_peak_MAM_(\\d{4})",
      values_to     = "swe_peak"
    ) %>%
    mutate(year = as.integer(year)) %>%
    filter(year >= 1984, year <= 2024)
} else {
  # Fallback: compute peak from monthly SWE (Mar/Apr/May) if columns exist
  peak_long <- lapply(years, function(y) {
    c03 <- paste0("SWE_", y, "03_TC")
    c04 <- paste0("SWE_", y, "04_TC")
    c05 <- paste0("SWE_", y, "05_TC")
    v03 <- if (c03 %in% names(df_gc)) as.numeric(df_gc[[c03]]) else NA_real_
    v04 <- if (c04 %in% names(df_gc)) as.numeric(df_gc[[c04]]) else NA_real_
    v05 <- if (c05 %in% names(df_gc)) as.numeric(df_gc[[c05]]) else NA_real_
    pk  <- pmax(v03, v04, v05, na.rm = TRUE)
    all_na <- is.na(v03) & is.na(v04) & is.na(v05)
    pk[all_na] <- NA_real_
    tibble(year = y, swe_peak = pk)
  }) %>% bind_rows()
}

peak_stats <- peak_long %>%
  group_by(year) %>%
  summarize(
    mean = mean(swe_peak, na.rm = TRUE),
    sd   = sd(swe_peak,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(metric = "SWE Peak (MAM)")

# --- Combine & add ±1 SD bounds 
stats <- bind_rows(
  ndvi_stats,
  apr_stats,
  peak_stats
) %>%
  mutate(lower = mean - sd, upper = mean + sd)

# --- Normalize (0–1) per metric for comparability 
stats_norm <- stats %>%
  group_by(metric) %>%
  mutate(
    minv = min(mean, na.rm = TRUE),
    maxv = max(mean, na.rm = TRUE),
    rng  = ifelse(maxv > minv, maxv - minv, NA_real_),
    mean_norm  = (mean  - minv) / rng,
    lower_norm = (lower - minv) / rng,
    upper_norm = (upper - minv) / rng
  ) %>%
  ungroup()

# --- Plot: NDVI + April SWE + Peak SWE (normalized) 
ggplot(stats_norm, aes(x = year, y = mean_norm, color = metric, fill = metric)) +
  geom_ribbon(aes(ymin = lower_norm, ymax = upper_norm), alpha = 0.20, color = NA) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = seq(1984, 2024, by = 4)) +
  labs(
    title    = paste0(huc_target, ": NDVI, April SWE, and Peak SWE (MAM), 1984–2024"),
    subtitle = "Means across all pixels within HUC12; ribbons are ±1 SD; each metric normalized to 0–1",
    x = "Year", y = "Normalized units (0–1)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.title = element_blank())

# --- Compute April SWE z-scores  ########
suppressPackageStartupMessages({
  library(dplyr); library(sf)
})

# --- Geometry-free copy
data_no_geo <- sf::st_drop_geometry(data_sf_clean)

# --- Years
years_all <- 1984:2024
years_30  <- 1986:2015  # 30-yr baseline (inclusive)

# --- Helpers
mk_mat <- function(d, cols) {
  if (length(cols) == 0) return(matrix(NA_real_, nrow(d), 0))
  as.matrix(as.data.frame(lapply(d[, cols, drop = FALSE], as.numeric)))
}
row_sds <- function(m) {
  mu <- rowMeans(m, na.rm = TRUE)
  sqrt(rowMeans((m - mu)^2, na.rm = TRUE))
}

# --- Collect available April SWE columns
apr_cols_all <- intersect(paste0("SWE_", years_all, "04_TC"), names(data_no_geo))
apr_cols_30  <- intersect(paste0("SWE_", years_30, "04_TC"), names(data_no_geo))

# --- Build matrices
M_all <- mk_mat(data_no_geo, apr_cols_all)   # April SWE for all available years
M_30  <- mk_mat(data_no_geo, apr_cols_30)    # April SWE for 30-yr baseline

# --- 30-yr per-pixel baseline stats
mu_30 <- if (ncol(M_30) > 0) rowMeans(M_30, na.rm = TRUE) else rep(NA_real_, nrow(data_no_geo))
sd_30 <- if (ncol(M_30) > 0) row_sds(M_30)   else rep(NA_real_, nrow(data_no_geo))
sd_30_safe <- ifelse(sd_30 > 0, sd_30, NA_real_)  # avoid divide-by-zero

# --- Attach baseline mean/sd to both frames
data_no_geo$sweavg_30yr_Apr_TC <- mu_30
data_no_geo$swesd_30yr_Apr_TC  <- sd_30
data_sf_clean$sweavg_30yr_Apr_TC <- mu_30
data_sf_clean$swesd_30yr_Apr_TC  <- sd_30

# --- Compute per-year April SWE z-scores: (April SWE - 30yr mean) / 30yr sd
if (ncol(M_all) > 0) {
  Z_all <- (M_all - mu_30) / sd_30_safe
  
  # Name z-score columns as swez_Apr_YYYY_TC
  z_names <- gsub("^SWE_(\\d{4})04_TC$", "swez_Apr_\\1_TC", apr_cols_all)
  colnames(Z_all) <- z_names
  Zdf <- as.data.frame(Z_all)
  
  # Write back (create or overwrite) to both frames
  for (nm in z_names) {
    data_no_geo[[nm]]   <- Zdf[[nm]]
    data_sf_clean[[nm]] <- Zdf[[nm]]
  }
}

# --- Quick check
data_sf_clean %>%
  dplyr::select(
    dplyr::matches("^SWE_\\d{4}04_TC$"),   # April SWE per year
    sweavg_30yr_Apr_TC, swesd_30yr_Apr_TC,
    dplyr::matches("^swez_Apr_\\d{4}_TC$")
  ) %>%
  glimpse()

################################################################################

###### ======= AET ===== 
# Geometry-free copy
df <- data_sf_clean %>% sf::st_drop_geometry()

years_all <- 1984:2024
mm_amjja  <- c("04","05","06","07","08")  # Apr–Aug
mm_jja    <- c("06","07","08")            # Jun–Aug

# Helper: numeric matrix from selected columns
mk_mat <- function(d, cols) {
  if (length(cols) == 0) return(matrix(NA_real_, nrow(d), 0))
  as.matrix(as.data.frame(lapply(d[, cols, drop = FALSE], as.numeric)))
}

# Helper: pull AET for a given month ("04".."08") and a vector of target years (per row)
pull_aet_for_year <- function(df, mm, year_vec) {
  cols_have <- intersect(paste0("AET_", years_all, mm, "_TC"), names(df))
  M    <- mk_mat(df, cols_have)
  want <- paste0("AET_", year_vec, mm, "_TC")
  j    <- match(want, cols_have)
  i    <- seq_len(nrow(df))
  out  <- rep(NA_real_, nrow(df))
  ok   <- !is.na(j) & ncol(M) > 0
  if (any(ok)) out[ok] <- M[cbind(i[ok], j[ok])]
  out
}

# Strict seasonal mean for an arbitrary vector of years and a month set.
# Returns NA if ANY of the months are NA or missing for that year.
strict_mean_for_year <- function(year_vec, months_vec) {
  mats <- lapply(months_vec, function(mm) pull_aet_for_year(df, mm, year_vec))
  mat  <- do.call(cbind, mats)
  rowMeans(mat, na.rm = FALSE)  # strict: any NA -> NA
}

# ======== Computing 5post metrics ==

# 1) AMJJA 5-post mean (ref_year+1 .. ref_year+5), strict across months & years
post_amjja_means <- sapply(1:5, function(k) strict_mean_for_year(df$ref_year + k, mm_amjja))
if (is.null(dim(post_amjja_means))) post_amjja_means <- matrix(post_amjja_means, ncol = 1)
aetavg_AMJJA_5post <- rowMeans(post_amjja_means, na.rm = FALSE)

# 2) August 5-post mean (ref_year+1 .. ref_year+5), strict across years
aug_5post <- sapply(1:5, function(k) pull_aet_for_year(df, "08", df$ref_year + k))
if (is.null(dim(aug_5post))) aug_5post <- matrix(aug_5post, ncol = 1)
aetavg_Aug_5post <- rowMeans(aug_5post, na.rm = FALSE)

# ======== Compute per-year JJA means (1984–2024) ==
jja_names <- paste0("aetavg_JJA_", years_all)
jja_vals  <- lapply(years_all, function(y) {
  cols <- paste0("AET_", y, mm_jja, "_TC")      # AET_YYYY06/07/08_TC
  # strict: if any monthly column for this year is missing -> NA for all pixels
  if (!all(cols %in% names(df))) return(rep(NA_real_, nrow(df)))
  rowMeans(mk_mat(df, cols), na.rm = FALSE)     # strict: any NA -> NA
})

# ---------- Attach to original sf data 
data_sf_clean <- data_sf_clean %>%
  mutate(
    aetavg_AMJJA_5post = aetavg_AMJJA_5post,
    aetavg_Aug_5post   = aetavg_Aug_5post,
    !!!setNames(jja_vals, jja_names)
  )

# ------------- Quick check 
data_sf_clean |>
  dplyr::select(
    aetavg_AMJJA_5post, aetavg_Aug_5post,
    dplyr::starts_with("aetavg_JJA_")
  ) |>
  dplyr::glimpse()
################################################################################

###############################
#### Final study area map: ####
# Pixel counts by sev_group
sev_counts <- data_sf_clean %>%
  count(sev_group) %>%
  mutate(label = paste0(sev_group, " (n = ", format(n, big.mark = ","), ")"))
sev_labels <- setNames(sev_counts$label, sev_counts$sev_group)

# Plot with sev_group counts and labels in legend
ggplot(data_sf_clean) +
  geom_sf(aes(color = sev_group), size = 0.2) +
  scale_color_manual(
    values = sev_colors,
    name = "Severity Group",
    labels = sev_labels  # show counts in the legend
  ) +
  labs(
    title = "Final - Pixels by Severity Group",
    subtitle = paste("N =", format(nrow(data_sf_clean), big.mark = ","), "pixels")
  ) +
  theme_minimal() +
  coord_sf()
################################################################################

# Write attribute table (no geometry) to CSV
# File: C:\Users\leahs\OneDrive\Documents\U_of_M\Masters_Project\Ch01_data_drop\data_sf_clean_08272025.csv
out_dir  <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop"

# Drop geometry; make factors explicit strings for safer CSV I/O
df_out <- data_sf_clean %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(dplyr::across(where(is.factor), as.character))

out_csv <- file.path(out_dir, "data_sf_clean_03022025_40mSPACING.csv")

# Write the CSV file
readr::write_csv(df_out, out_csv)

cat("Wrote", nrow(df_out), "rows ×", ncol(df_out), "cols to:\n", out_csv, "\n")
