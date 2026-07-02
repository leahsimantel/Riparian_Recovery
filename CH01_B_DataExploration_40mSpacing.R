###############################################################
###    CH01 Script B: Data Exploration, conversion to data_long
###    11/08/2025    40m spacing of pixels.

# Load necessary libraries --------------------------
library(dplyr)
library(readr)
library(ggplot2)
library(sf)
library(tidyr)
library(stringr)
library(e1071)   # for skewness
library(scales)
library(moments) # for skewness
library(geojsonsf)
library(writexl)
library(purrr)
library(jsonlite)
library(slider)
library(haven)
# ---------------------------------------------------

csv_path <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/data_sf_clean_03022025_40mSPACING.csv"

df <- readr::read_csv(csv_path, show_col_types = FALSE)  # tibble
dim(df); names(df)[1:10]  

## Rename file for downstream code:
data_sf_reduced <- df

# Check for any duplicate pixels 
## ---- Duplicate pixel_ID check (pre-pivot) -----------------------------------
stopifnot("pixel_ID" %in% names(data_sf_reduced))

dup_tbl <- data_sf_reduced %>%
  dplyr::count(pixel_ID, name = "n") %>%
  dplyr::filter(n > 1) %>%
  dplyr::arrange(dplyr::desc(n))

cat(
  "Total rows:", nrow(data_sf_reduced), "\n",
  "Unique pixel_IDs:", dplyr::n_distinct(data_sf_reduced$pixel_ID), "\n",
  "Duplicate pixel_IDs (pixel_IDs with n>1):", nrow(dup_tbl), "\n"
)

if (nrow(dup_tbl) > 0) {
  cat("\nTop duplicate pixel_IDs:\n")
  print(utils::head(dup_tbl, 20))
  
  cat("\nExample duplicate rows (first 10 rows):\n")
  data_sf_reduced %>%
    dplyr::semi_join(dup_tbl, by = "pixel_ID") %>%
    dplyr::arrange(pixel_ID) %>%
    dplyr::slice_head(n = 10) %>%
    print()
}

# 1) Count duplicates
dup_counts <- data_sf_reduced %>%
  count(pixel_ID, name = "n") %>%
  filter(n > 1)

cat("Duplicate pixel_IDs:", nrow(dup_counts), "\n")  #70, doubled (so, 140 rows)
## note... these were twice-burned pixels that slipped through original QA because
##         the shapefile did not quite cover them, and they were falsely assigned a 
##         0 or 1 to num_times_burned

## Delete the 140 rows of twice-burned 
# --- 1) Identify duplicated pixel_IDs (appear > 1 time)
dup_tbl <- data_sf_reduced %>%
  count(pixel_ID, name = "n") %>%
  filter(n > 1)

dup_ids <- dup_tbl$pixel_ID

cat("Total rows (before):", nrow(data_sf_reduced), "\n")  ## 11707
cat("Distinct pixel_IDs (before):", n_distinct(data_sf_reduced$pixel_ID), "\n")
cat("Duplicated pixel_IDs found:", length(dup_ids), "\n")  ## 151

# --- 2) Drop ALL rows for these duplicated pixel_IDs
data_sf_reduced <- data_sf_reduced %>%
  filter(!pixel_ID %in% dup_ids)

# --- 3) Verify
#cat("\nAfter removal:\n")
cat("Total rows (after):", nrow(data_sf_reduced), "\n")  ## 11405
cat("Distinct pixel_IDs (after):", n_distinct(data_sf_reduced$pixel_ID), "\n") ## 11405

### ===== Assign 'sev_num == 0' for all sev_group == Unburned pixels ------

data_sf_reduced <- data_sf_reduced %>%
  mutate(
    sev_num = case_when(
      is.na(sev_num) & sev_group == "Unburned" ~ 0,
      TRUE ~ sev_num
    )
  )

### Fix 2 mis-matched control pixels (correctly assign ref_year) ------------
# The two mis-assigned pixels:
gabe_bad_pixel_id <- "47.4320463_-113.2945576"      # should map to 2009
ls_bad_pixel_id   <- "47.6447795_-113.4380834"      # should map to 2003

data_sf_reduced <- data_sf_reduced %>%
  mutate(
    # Fix 1: Gabe Creek pixel
    fire_name = case_when(
      pixel_ID == gabe_bad_pixel_id ~ "None",
      pixel_ID == ls_bad_pixel_id   ~ "None",
      TRUE                          ~ fire_name
    ),
    
    sev_group = case_when(
      pixel_ID %in% c(gabe_bad_pixel_id, ls_bad_pixel_id) ~ "Unburned",
      TRUE ~ sev_group
    ),
    
    ref_year = case_when(
      # assign corrected ignition-year to the control pixel
      pixel_ID == gabe_bad_pixel_id ~ 2009,
      pixel_ID == ls_bad_pixel_id   ~ 2003,
      TRUE ~ ref_year
    )
  )

### Optional verification:
cat("Corrected Gabe Creek 2009 pixel:\n")
print(
  data_sf_reduced %>%
    filter(pixel_ID == gabe_bad_pixel_id) %>%
    distinct(pixel_ID, fire_name, sev_group, ref_year)
)

cat("Corrected Little Salmon 2003 pixel:\n")
print(
  data_sf_reduced %>%
    filter(pixel_ID == ls_bad_pixel_id) %>%
    distinct(pixel_ID, fire_name, sev_group, ref_year)
)

# Remove any lingering pixels that were assigned to the Cannon fire ------------
data_sf_reduced <- data_sf_reduced %>%
  dplyr::filter(fire_name != "CANNON_2023")

# Remove any pixels that don't have a 5-yr CWD record post-fire  ---------------
sum(is.na(data_sf_reduced$cwd_5yr_postfire_avg_08))  # 84

data_sf_reduced <- data_sf_reduced %>%
  filter(!is.na(cwd_5yr_postfire_avg_08))

sum(is.na(data_sf_reduced$cwd_5yr_postfire_avg_08))  # 0

# ==============================================================
# Uploading New Columns, 4/1/26: Elevation & TopoTerra AET, CWD, Slope, TOTDA
# ==============================================================
# Define the folder containing the new CSV files with AET, CWD, Elevation
folder_path_new <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/NDVI_covariates_by_Pixel/CSVs_04012026"

# Get list of all CSV files in the folder
csv_files_new <- list.files(
  path = folder_path_new,
  pattern = "\\.csv$",
  full.names = TRUE
)

# Read and combine all CSV files
raw_covars_new_fresh <- bind_rows(
  lapply(csv_files_new, read_csv, show_col_types = FALSE)
)

# Drop system:index if present
raw_covars_new <- raw_covars_new_fresh %>%
  select(-any_of("system:index"))

# Explicitly set column types and round AET/CWD
raw_covars_new <- raw_covars_new %>%
  mutate(
    # IDs
    pixel_ID = as.character(pixel_ID),
    .geo     = as.character(.geo),
    
    # Core numerics
    aet_30yrAvg_TT = round(as.numeric(aet_30yrAvg_TT), 2),
    cwd_30yrAvg_TT = round(as.numeric(cwd_30yrAvg_TT), 2),
    elevation_m    = as.numeric(elevation_m),
    latitude       = as.numeric(latitude),
    longitude      = as.numeric(longitude),
    slope_deg      = as.numeric(slope_deg),
    TotDASqKm      = as.numeric(TotDASqKm)
  )

# Verify structure
str(raw_covars_new)
summary(raw_covars_new)

#### MERGE DATASETS
covars_to_join <- raw_covars_new %>%
  select(pixel_ID, aet_30yrAvg_TT, cwd_30yrAvg_TT, elevation_m, slope_deg, TotDASqKm)

data_sf_reduced <- data_sf_reduced %>%
  mutate(pixel_ID = as.character(pixel_ID))
data_sf_reduced <- data_sf_reduced %>%
  left_join(covars_to_join, by = "pixel_ID")

####  QA - confirm no extra rows were added to data_sf_reduced
n_before <- nrow(data_sf_reduced)
# (re-run join in a temporary object if needed)
data_sf_reduced_new <- data_sf_reduced %>%
  left_join(covars_to_join, by = "pixel_ID")

n_after <- nrow(data_sf_reduced_new)
n_before == n_after

data_sf_reduced %>%
  summarise(
    missing_aet = sum(is.na(aet_30yrAvg_TT)),
    missing_cwd = sum(is.na(cwd_30yrAvg_TT)),
    missing_elev = sum(is.na(elevation_m)),
    missing_slope = sum(is.na(slope_deg)),
    missing_TotDA = sum(is.na(TotDASqKm))
  )

## QA... mapping out pixels with missing NEW variables
# Pixels in data_long_complete that failed to match
missing_pixels <- data_sf_reduced %>%
  filter(
    is.na(aet_30yrAvg_TT) |
      is.na(cwd_30yrAvg_TT) |
      is.na(elevation_m) |
      is.na(slope_deg) |
      is.na(TotDASqKm)
  ) %>%
  distinct(pixel_ID, latitude, longitude)

nrow(missing_pixels)  # 6
library(sf)

missing_sf <- missing_pixels %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# ---- Static map layers 
huc12 <- st_read("C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/ArcGIS/BMWC_Fire_Data/HUC12_BMWA/HUC12_BobMarshall.shp")
streams <- st_read("C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/ArcGIS/BMWC_Fire_Data/hydro/BMWA_Flowline.shp")
fires <- st_read("C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/ArcGIS/BMWC_Fire_Data/fire_atlas/fire_atlas_Jan2025/fire_atlas.shp")

#ggplot() +
# HUC12 boundaries (background context)
#  geom_sf(data = huc12, fill = NA, color = "black", size = 0.4) +
#  geom_sf(data = missing_sf, color = "red", size = 1) +
#  theme_minimal() +
#  labs(title = "Pixels Missing AET/CWD/Elevation",
#       subtitle = "Red = unmatched pixels (n = 6)")

### REMOVE PIXELS WITH MISSING TOPOTERRA / ELEV /TOTDASQKM DATA
pixels_to_remove <- missing_pixels %>%
  distinct(pixel_ID) %>%
  pull(pixel_ID)

length(pixels_to_remove)  # 6 pixels.

data_sf_reduced <- data_sf_reduced %>%
  filter(!pixel_ID %in% pixels_to_remove)
## =============================================================================
#------------------------------------------------------------------------------

#########   FIRE DATA SUMMARIES:   ###########   Using data_sf_reduced -----

# Summary of pixel count per fire, year, and severity class
fire_summary <- data_sf_reduced %>%
  filter(sev_class %in% c("Low", "Moderate", "High")) %>%
  group_by(fire_name, fire_year, sev_class) %>%
  summarise(pixel_count = n(), .groups = "drop")

print(n_distinct(fire_summary$fire_name)) ##  51 fire events.

# Total pixels per fire to calculate percentages
fire_totals <- fire_summary %>%
  group_by(fire_name, fire_year) %>%
  summarise(total_pixels = sum(pixel_count), .groups = "drop")

# Get unique fire sizes
fire_sizes <- data_sf_reduced %>%
  st_drop_geometry() %>%
  distinct(fire_name, fire_year, fire_size_ha)

# Create one-row-per-fire table with severity percentages in separate columns
fire_summary_with_pct <- fire_summary %>%
  left_join(fire_totals, by = c("fire_name", "fire_year")) %>%
  mutate(percent = 100 * pixel_count / total_pixels) %>%
  select(fire_name, fire_year, sev_class, percent) %>%
  tidyr::pivot_wider(
    names_from  = sev_class,
    values_from = percent
  ) %>%
  rename(
    Percent_Low  = Low,
    Percent_Mod  = Moderate,
    Percent_High = High
  ) %>%
  left_join(fire_totals, by = c("fire_name", "fire_year")) %>%
  left_join(fire_sizes,  by = c("fire_name", "fire_year")) %>%
  mutate(
    Percent_Low  = round(Percent_Low, 1),
    Percent_Mod  = round(Percent_Mod, 1),
    Percent_High = round(Percent_High, 1)
  ) %>%
  select(
    fire_name,
    fire_year,
    fire_size_ha,
    Percent_Low,
    Percent_Mod,
    Percent_High
  ) %>%
  arrange(fire_year)

print(fire_summary_with_pct)

# Reshape one-row-per-fire summary table to long format for plotting 
fire_plot_df <- fire_summary_with_pct %>%
  tidyr::pivot_longer(
    cols = c(Percent_Low, Percent_Mod, Percent_High),
    names_to = "sev_class",
    values_to = "percent"
  ) %>%
  dplyr::mutate(
    sev_class = dplyr::recode(
      sev_class,
      "Percent_Low"  = "Low",
      "Percent_Mod"  = "Moderate",
      "Percent_High" = "High"
    ),
    area_ha = fire_size_ha * (percent / 100)
  )

fire_plot_df <- fire_plot_df %>%
  dplyr::arrange(fire_year, fire_name) %>%
  dplyr::mutate(
    fire_name = factor(fire_name, levels = unique(fire_name))
  )

# Plot bar chart with uniform scale (fire events)
ggplot(fire_plot_df, aes(x = fire_name, y = percent, fill = sev_class)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(
    values = c("Low" = "yellow", "Moderate" = "orange", "High" = "red"),
    name = "Severity Class"
  ) +
  theme_minimal() +
  labs(
    title = "Distribution of Burn Severity by Fire (Chronological Order)",
    x = "Fire Name (Ordered by Year)",
    y = "Percent of Pixels",
    caption = "Stacked bar shows % of Low, Moderate, and High severity pixels per fire"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

## Plot bar chart with scale based on fire size
ggplot(fire_plot_df, aes(x = fire_name, y = area_ha, fill = sev_class)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(
    values = c("Low" = "yellow", "Moderate" = "orange", "High" = "red"),
    name = "Severity Class"
  ) +
  theme_minimal() +
  labs(
    title = "Fire Size by Severity Class (Stacked by Area)",
    x = "Fire (Ordered by Year)",
    y = "Area Burned (ha)",
    caption = "Bar height = fire size; colors show severity class distribution"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )
################################################################################

###### Insert, 8/5/25: what does raw NDVI variation look like over time? -----

# Filter to Upper Danaher Creek and unburned pixels
#udc_unburned <- data_sf_reduced %>%
#  filter(huc12 == "Upper Danaher Creek",
#         sev_class %in% c("None", "Unburned/Very Low"))
udc_unburned <- data_sf_reduced %>%
  filter(sev_class %in% c("None", "Unburned/Very Low"))


# Pivot NDVI columns to obtain YYYY
ndvi_long <- udc_unburned %>%
  st_drop_geometry() %>%  # remove geometry column
  select(pixel_ID, matches("^NDVI_\\d{4}$")) %>%
  pivot_longer(
    cols = -pixel_ID,
    names_to = "year",
    names_prefix = "NDVI_",
    values_to = "NDVI"
  ) %>%
  mutate(year = as.integer(year))

# Summarize mean and SD NDVI per year
ndvi_summary <- ndvi_long %>%
  group_by(year) %>%
  summarize(
    mean_ndvi = mean(NDVI, na.rm = TRUE),
    sd_ndvi = sd(NDVI, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate the number of unique pixels used
n_pixels <- n_distinct(ndvi_long$pixel_ID)

# Create plot with subtitle
ggplot(ndvi_summary, aes(x = year, y = mean_ndvi)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_ribbon(aes(ymin = mean_ndvi - sd_ndvi, ymax = mean_ndvi + sd_ndvi),
              alpha = 0.2, fill = "forestgreen") +
  labs(
    title = "NDVI Over Time (Unburned Pixels)",
    subtitle = paste("Sample size:", n_pixels, "pixels"),
    x = "Year",
    y = "Mean NDVI ± SD"
  ) +
  theme_minimal()

# Pivot CWD columns (Note- Using Avg August CWD values only) to obtain YYYY
cwd_long <- udc_unburned %>%
  st_drop_geometry() %>%  # remove geometry column
  select(pixel_ID, matches("^CWD_\\d{4}08_TC$")) %>%
  pivot_longer(
    cols = -pixel_ID,
    names_to = "year",
    names_pattern = "^CWD_(\\d{4})08_TC$",
    values_to = "CWD_aug"
  ) %>%
  mutate(year = as.integer(year))

# Merge NDVI and CWD long tables
ndvi_cov_combined <- ndvi_long %>%
  left_join(cwd_long, by = c("pixel_ID", "year"))

# Summarize both over time
ndvi_cov_summary <- ndvi_cov_combined %>%
  group_by(year) %>%
  summarize(
    mean_ndvi = mean(NDVI, na.rm = TRUE),
    sd_ndvi = sd(NDVI, na.rm = TRUE),
    mean_cwd = mean(CWD_aug, na.rm = TRUE),
    sd_cwd = sd(CWD_aug, na.rm = TRUE),
    .groups = "drop"
  )

ndvi_cov_summary <- ndvi_cov_summary %>%
  mutate(
    ndvi_scaled = (mean_ndvi - min(mean_ndvi)) / (max(mean_ndvi) - min(mean_ndvi)),
    cwd_scaled = (mean_cwd - min(mean_cwd)) / (max(mean_cwd) - min(mean_cwd))
  )

# Plot time series of average NDVI and average August CWD across years 
ggplot(ndvi_cov_summary, aes(x = year)) +
  geom_line(aes(y = ndvi_scaled, color = "NDVI"), linewidth = 1) +
  geom_line(aes(y = cwd_scaled, color = "CWD"), linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("NDVI" = "darkgreen", "CWD" = "orange")) +
  labs(
    title = "Normalized NDVI and August CWD Over Time",
    subtitle = paste("Unburned Pixels (n =", n_pixels, "pixels)"),
    x = "Year",
    y = "Scaled Value (0–1)",
    color = "Variable"
  ) +
  theme_minimal()

## Plotting using 2 y axes:
library(ggplot2)
library(scales)

ggplot(ndvi_cov_summary, aes(x = year)) +
  geom_line(aes(y = mean_ndvi, color = "NDVI"), linewidth = 1) +
  geom_line(aes(y = mean_cwd / 1000, color = "CWD"), linewidth = 1, linetype = "dashed") +  # scale down CWD
  scale_color_manual(values = c("NDVI" = "darkgreen", "CWD" = "orange")) +
  scale_y_continuous(
    name = "NDVI",
    sec.axis = sec_axis(~ . * 1000, name = "CWD (approx)")
  ) +
  labs(
    title = "NDVI and August CWD Over Time (Unscaled)",
    subtitle = paste("Unburned Pixels (n =", n_pixels, "pixels)"),
    x = "Year", color = "Variable"
  ) +
  theme_minimal()

# Plot time series of average NDVI and average tmeanz_JJA across years 
# Pivot TEMP columns to obtain YYYY
tmean_long <- udc_unburned %>%
  st_drop_geometry() %>%
  select(pixel_ID, matches("^tmeanz_JJA_\\d{4}_PR$")) %>%
  pivot_longer(
    cols = -pixel_ID,
    names_to = "year",
    names_pattern = "^tmeanz_JJA_(\\d{4})_PR$",
    values_to = "TMEAN_JJA_z"
  ) %>%
  mutate(year = as.integer(year))

# Merge NDVI and TEMP
ndvi_temp_combined <- ndvi_long %>%
  left_join(tmean_long, by = c("pixel_ID", "year"))

# Summarize over time
ndvi_temp_summary <- ndvi_temp_combined %>%
  group_by(year) %>%
  summarize(
    mean_ndvi = mean(NDVI, na.rm = TRUE),
    sd_ndvi   = sd(NDVI, na.rm = TRUE),
    mean_tz   = mean(TMEAN_JJA_z, na.rm = TRUE),
    sd_tz     = sd(TMEAN_JJA_z, na.rm = TRUE),
    .groups = "drop"
  )

# Scaled (0–1) overlay
ndvi_temp_summary <- ndvi_temp_summary %>%
  mutate(
    ndvi_scaled = (mean_ndvi - min(mean_ndvi, na.rm = TRUE)) /
      (max(mean_ndvi, na.rm = TRUE) - min(mean_ndvi, na.rm = TRUE)),
    tmean_scaled = (mean_tz - min(mean_tz, na.rm = TRUE)) /
      (max(mean_tz, na.rm = TRUE) - min(mean_tz, na.rm = TRUE))
  )

ggplot(ndvi_temp_summary, aes(x = year)) +
  geom_line(aes(y = ndvi_scaled, color = "NDVI"), linewidth = 1) +
  geom_line(aes(y = tmean_scaled, color = "Temp (JJA, z)"), linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("NDVI" = "darkgreen", "Temp (JJA, z)" = "#FFA500")) +
  labs(
    title = "Normalized NDVI and Standardized Summer Temp Over Time",
    subtitle = paste("Unburned Pixels (n =", n_pixels, "pixels)"),
    x = "Year",
    y = "Scaled Value (0–1)",
    color = "Variable"
  ) +
  theme_minimal()

# Dual-axis (unscaled) with a linear transform that maps Temp onto the NDVI axis
ndvi_range <- diff(range(ndvi_temp_summary$mean_ndvi, na.rm = TRUE))
tmean_range <- diff(range(ndvi_temp_summary$mean_tz,   na.rm = TRUE))
scale_t <- if (is.finite(tmean_range) && tmean_range > 0) ndvi_range / tmean_range else 1
offset_t <- min(ndvi_temp_summary$mean_ndvi, na.rm = TRUE) - scale_t * min(ndvi_temp_summary$mean_tz, na.rm = TRUE)

ggplot(ndvi_temp_summary, aes(x = year)) +
  geom_line(aes(y = mean_ndvi, color = "NDVI"), linewidth = 1) +
  geom_line(aes(y = scale_t * mean_tz + offset_t, color = "Temp (JJA, z)"),
            linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("NDVI" = "darkgreen", "Temp (JJA, z)" = "#FFA500")) +
  scale_y_continuous(
    name = "NDVI",
    sec.axis = sec_axis(~ (. - offset_t) / scale_t, name = "JJA Temperature (z)")
  ) +
  labs(
    title = "NDVI and JJA Temperature (z) Over Time (Unscaled Axes)",
    subtitle = paste("Unburned Pixels (n =", n_pixels, "pixels)"),
    x = "Year", color = "Variable"
  ) +
  theme_minimal()
# -------------------------------------------
# ============================
# NDVI + PRECIPITATION (pptz_JJA_%Y_PR)
# ============================

# Pivot PPT columns to obtain YYYY
ppt_long <- udc_unburned %>%
  st_drop_geometry() %>%
  select(pixel_ID, matches("^pptz_JJA_\\d{4}_PR$")) %>%
  pivot_longer(
    cols = -pixel_ID,
    names_to = "year",
    names_pattern = "^pptz_JJA_(\\d{4})_PR$",
    values_to = "PPT_JJA_z"
  ) %>%
  mutate(year = as.integer(year))

# Merge NDVI and PPT
ndvi_ppt_combined <- ndvi_long %>%
  left_join(ppt_long, by = c("pixel_ID", "year"))

# Summarize over time
ndvi_ppt_summary <- ndvi_ppt_combined %>%
  group_by(year) %>%
  summarize(
    mean_ndvi = mean(NDVI, na.rm = TRUE),
    sd_ndvi   = sd(NDVI, na.rm = TRUE),
    mean_pz   = mean(PPT_JJA_z, na.rm = TRUE),
    sd_pz     = sd(PPT_JJA_z, na.rm = TRUE),
    .groups = "drop"
  )

# Scaled (0–1) overlay
ndvi_ppt_summary <- ndvi_ppt_summary %>%
  mutate(
    ndvi_scaled = (mean_ndvi - min(mean_ndvi, na.rm = TRUE)) /
      (max(mean_ndvi, na.rm = TRUE) - min(mean_ndvi, na.rm = TRUE)),
    ppt_scaled  = (mean_pz - min(mean_pz, na.rm = TRUE)) /
      (max(mean_pz, na.rm = TRUE) - min(mean_pz, na.rm = TRUE))
  )

ggplot(ndvi_ppt_summary, aes(x = year)) +
  geom_line(aes(y = ndvi_scaled, color = "NDVI"), linewidth = 1) +
  geom_line(aes(y = ppt_scaled,  color = "Precip (JJA, z)"), linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("NDVI" = "darkgreen", "Precip (JJA, z)" = "#003366")) +
  labs(
    title = "Normalized NDVI and JJA Precipitation (z) Over Time",
    subtitle = paste("Unburned Pixels (n =", n_pixels, "pixels)"),
    x = "Year",
    y = "Scaled Value (0–1)",
    color = "Variable"
  ) +
  theme_minimal()

# Dual-axis (unscaled) with a linear transform that maps Precip onto the NDVI axis
ndvi_range2 <- diff(range(ndvi_ppt_summary$mean_ndvi, na.rm = TRUE))
ppt_range   <- diff(range(ndvi_ppt_summary$mean_pz,   na.rm = TRUE))
scale_p <- if (is.finite(ppt_range) && ppt_range > 0) ndvi_range2 / ppt_range else 1
offset_p <- min(ndvi_ppt_summary$mean_ndvi, na.rm = TRUE) - scale_p * min(ndvi_ppt_summary$mean_pz, na.rm = TRUE)

ggplot(ndvi_ppt_summary, aes(x = year)) +
  geom_line(aes(y = mean_ndvi, color = "NDVI"), linewidth = 1) +
  geom_line(aes(y = scale_p * mean_pz + offset_p, color = "Precip (JJA, z)"),
            linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("NDVI" = "darkgreen", "Precip (JJA, z)" = "#003366")) +
  scale_y_continuous(
    name = "NDVI",
    sec.axis = sec_axis(~ (. - offset_p) / scale_p, name = "JJA Precipitation (z)")
  ) +
  labs(
    title = "NDVI and JJA Precipitation (z) Over Time (Unscaled Axes)",
    subtitle = paste("Unburned Pixels (n =", n_pixels, "pixels)"),
    x = "Year", color = "Variable"
  ) +
  theme_minimal()

################################################################################

############    PIVOT DATA & CALCULATE DELTA NDVI: create data_long_reduced  ###
data_long_reduced <- data_sf_reduced %>%
  filter(!is.na(ndvi_prefire_3yr_avg), !is.na(sev_class)) %>%
  pivot_longer(
    cols = matches("^NDVI_\\d{4}$"),  
    names_to = "year",
    values_to = "ndvi_postfire"
  ) %>%
  mutate(
    year = as.integer(gsub("NDVI_", "", year)),
    years_since_fire = year - ref_year,
    delta_ndvi_avg = ndvi_postfire - ndvi_prefire_3yr_avg,
    delta_ndvi_min = ndvi_postfire - ndvi_prefire_3yr_min,
    delta_ndvi_med = ndvi_postfire - ndvi_prefire_3yr_med
  ) %>%
  # FILTER: keep 3 years before, fire year, and all postfire years 
filter(year >= ref_year - 3)
################################################################################

##### --------- Pivoting climate data and creating model_df:  ------------------

# Helper: pull a single climate column by sprintf() format; return NA vector if missing
pull_col <- function(df, fmt, y) {
  nm <- sprintf(fmt, y)
  if (nm %in% names(df)) as.numeric(df[[nm]]) else rep(NA_real_, nrow(df))
}

# Years actually present in pivoted table
years_needed <- sort(unique(data_long_reduced$year))

# Preallocate results list and a progress bar
res_list <- vector("list", length(years_needed))
pb <- txtProgressBar(min = 0, max = length(years_needed), style = 3)

for (i in seq_along(years_needed)) {
  y <- years_needed[i]
  
  # NDVI rows for this year only — enforce one row per (pixel_ID, year)
  ndvi_y <- data_long_reduced %>%
    dplyr::filter(year == !!y) %>%
    dplyr::arrange(dplyr::desc(ref_year)) %>%                 # keep most recent fire if duplicates
    dplyr::distinct(pixel_ID, year, .keep_all = TRUE)
  
  if (nrow(ndvi_y) > 0L) {
    # Climate chunk for this year (may have dup pixel_IDs in data_sf_reduced)
    clim_y_raw <- tibble::tibble(
      pixel_ID     = data_sf_reduced$pixel_ID,
      year         = y,
      ppttot_JJA   = pull_col(data_sf_reduced, "ppttot_JJA_%d_PR", y),
      pptz_JJA     = pull_col(data_sf_reduced, "pptz_JJA_%d_PR",   y),
      tmaxmean_JJA = pull_col(data_sf_reduced, "tmaxmean_JJA_%d_PR", y),
      tmaxz_JJA    = pull_col(data_sf_reduced, "tmaxz_JJA_%d_PR",    y),
      tmeanavg_JJA = pull_col(data_sf_reduced, "tmeanavg_JJA_%d_PR", y),
      tmeanz_JJA   = pull_col(data_sf_reduced, "tmeanz_JJA_%d_PR",   y),
      swe_peak_MAM = pull_col(data_sf_reduced, "swe_peak_MAM_%d",    y),
      SWE_Apr      = pull_col(data_sf_reduced, "SWE_%d04_TC",        y),
      swez_Apr     = pull_col(data_sf_reduced, "swez_Apr_%d_TC",     y),
      aetavg_JJA   = pull_col(data_sf_reduced, "aetavg_JJA_%d",      y)
    )
    
    # Columns to summarise (exclude keys explicitly to avoid tidyselect errors)
    data_cols <- setdiff(names(clim_y_raw), c("pixel_ID", "year"))
    
    # Collapse to one row per (pixel_ID, year) on the climate side
    first_non_na <- function(v) {
      idx <- which(!is.na(v))
      if (length(idx) == 0) NA_real_ else v[idx[1]]
    }
    
    clim_y <- clim_y_raw %>%
      dplyr::group_by(pixel_ID, year) %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(data_cols), first_non_na), .groups = "drop")
    
    # Safe many-to-one join for this year only
    res_list[[i]] <- ndvi_y %>%
      dplyr::left_join(clim_y, by = c("pixel_ID", "year"), relationship = "many-to-one")
  } else {
    res_list[[i]] <- ndvi_y
  }
  
  # Update counter; occasional GC to keep memory tidy
  setTxtProgressBar(pb, i)
  if (i %% 5 == 0) gc(FALSE)
}
close(pb)

# Bind all years
model_df <- dplyr::bind_rows(res_list)

# --- Quick QA ---
model_df %>%
  select(pixel_ID, year, years_since_fire, delta_ndvi_min, sev_group,
         ppttot_JJA, pptz_JJA,
         tmaxmean_JJA, tmaxz_JJA,
         tmeanavg_JJA, tmeanz_JJA,
         swe_peak_MAM, SWE_Apr, swez_Apr,
         aetavg_JJA) %>%
  glimpse()

# Counts:
# Count unique pixel_IDs
n_pixel_ids <- model_df %>%
  distinct(pixel_ID) %>%
  nrow()

# Count unique geographic pixels (latitude + longitude)
n_geo_pixels <- model_df %>%
  distinct(latitude, longitude) %>%
  nrow()

cat("Unique pixel_IDs:", n_pixel_ids, "\n")  # 11329
cat("Unique geographic pixels:", n_geo_pixels, "\n")  # 11329
################################################################################

######  Filtering out pixels with low MINIMUM pre-fire NDVI   ##############

# Function: robustly convert a .geo column (character or list) to sfc
geojson_to_sfc <- function(x) {
  if (inherits(x, "sfc")) {
    return(x)
  }
  if (is.list(x)) {
    x <- purrr::map_chr(x, function(el) {
      if (is.null(el)) {
        return(NA_character_)
      } else if (is.character(el)) {
        return(el)
      } else {
        return(jsonlite::toJSON(el, auto_unbox = TRUE))
      }
    })
  } else if (!is.character(x)) {
    x <- as.character(x)
  }
  sf::st_as_sfc(x, GeoJSON = TRUE)
}

# 1) Build a per-pixel geometry table from model_df's .geo
pixel_geom_sf <- model_df %>%
  dplyr::select(pixel_ID, .geo) %>%
  dplyr::distinct(pixel_ID, .keep_all = TRUE) %>%  # access each pixel only once
  dplyr::mutate(geometry = geojson_to_sfc(.geo)) %>%
  sf::st_as_sf(crs = 4326) %>%
  sf::st_make_valid()


# 2) Identify pixels with very low prefire NDVI (ndvi_prefire_3yr_min <= 0.2)
low_ids <- model_df %>%
  dplyr::filter(!is.na(ndvi_prefire_3yr_min), ndvi_prefire_3yr_min <= 0.2) %>%
  dplyr::distinct(pixel_ID) %>%
  dplyr::pull(pixel_ID)

# 3) Quick counts
n_total <- dplyr::n_distinct(model_df$pixel_ID)
n_low   <- length(low_ids)
pct_low <- round(100 * n_low / n_total, 2)
message("Low-prefire pixels (<= 0.2): ", n_low, " / ", n_total, " (", pct_low, "%)")

# 4) Map those pixels spatially
low_prefire_sf <- pixel_geom_sf %>%
  dplyr::filter(pixel_ID %in% low_ids)

# --- Read in HUC12 & streams shapefiles
huc12_sf <- sf::st_read("C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/ArcGIS/BMWC_Fire_Data/HUC12_BMWA/HUC12_BobMarshall.shp")
streams_sf <- sf::st_read("C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/ArcGIS/BMWC_Fire_Data/NHDPlus_V2_Streamlines/BMWA_streams_widened.shp")

# --- Make sure map layers match CRS of pixel geometries
huc12_bg <- huc12_sf %>% sf::st_transform(sf::st_crs(pixel_geom_sf))
streams_bg <- streams_sf %>% sf::st_transform(sf::st_crs(pixel_geom_sf))

#install.packages("extrafont")
#library(extrafont)
#font_import(prompt = FALSE)   # scans your system fonts (takes a few minutes)
#loadfonts(device = "win")     # or device = "pdf" for PDF output

# --- Map with HUC12 backdrop + pixels on top
gg_low_prefire <- ggplot() +
  geom_sf(data = huc12_bg, fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_sf(data = streams_bg, size = 1, color = "darkblue") +
  geom_sf(data = low_prefire_sf, size = 1.4, alpha = 0.7, color = "darkred") +
  coord_sf(expand = FALSE) +
  labs(
    title = "BMWA Pixels Excluded from Analysis",
    subtitle = paste0("N = ", n_low, " (", pct_low, "% of pixels)"),
    caption = "Filter Based on Low Pre-fire NDVI (≤ 0.2)"
  ) +
  theme_minimal(base_size = 14, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12),
    plot.caption = element_text(size = 10)
  )

print(gg_low_prefire)

# 5) Delete those pixels from data frame
# --- Remove low-prefire pixels (all years) from model_df ---
model_df <- model_df %>%
  dplyr::filter(!pixel_ID %in% low_ids)

# Quick check
n_total_after <- dplyr::n_distinct(model_df$pixel_ID)
message("Remaining pixels after removal: ", n_total_after, " / ", n_total,
        " (", round(100 * n_total_after / n_total, 2), "% kept)")
################################################################################

#####   Merge wetted width data to pixels, and plot results   ##################
# Read the simplified stream network CSV into R 
width_data <- read.csv(
  "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/stream_network_poly_simple.csv",
  stringsAsFactors = FALSE
)

# Keep only the columns you want to join from width_data (and make them unique by COMID)
width_data_join <- width_data %>%
  dplyr::select(COMID, WettedWidth_m, WettedWidth_Source, StreamOrde, StreamLeve) %>%
  dplyr::distinct(COMID, .keep_all = TRUE)

# Left join onto model_df (preserves all rows in model_df; repeats values for repeated COMIDs)
model_df <- model_df %>%
  dplyr::left_join(width_data_join, by = "COMID")

#####   Define 2.5m-interval Wetted Width classes --------
# One wetted width value per pixel
pixel_width_df <- model_df %>%
  dplyr::select(pixel_ID, WettedWidth_m) %>%
  dplyr::distinct(pixel_ID, .keep_all = TRUE)

# Define breaks (aligned at 0, 2.5, 5, ...)
width_breaks <- seq(
  from = 0,
  to   = ceiling(max(pixel_width_df$WettedWidth_m, na.rm = TRUE) / 2.5) * 2.5,
  by   = 2.5
)

# Create width class factor with readable labels
pixel_width_df <- pixel_width_df %>%
  mutate(
    WettedWidth_class = cut(
      WettedWidth_m,
      breaks = width_breaks,
      right  = FALSE,
      include.lowest = TRUE,
      labels = paste0(
        "[",
        head(width_breaks, -1), "–",
        tail(width_breaks, -1), " m)"
      )
    )
  )

# Correction for legend labels:
levels(pixel_width_df$WettedWidth_class) <-
  gsub("\\)", "]", levels(pixel_width_df$WettedWidth_class))

# Discrete color palette for map and histogram
width_levels <- levels(pixel_width_df$WettedWidth_class)

library(viridisLite)

# Levels for your discrete width classes
width_levels <- levels(pixel_width_df$WettedWidth_class)
n_bins <- length(width_levels)

# Nonlinear spacing along the gradient:
# power < 1 spreads colors out more near the low end (first bins)
power <- 0.3
t_vals <- seq(0, 1, length.out = n_bins) ^ power

# Build palette (still a smooth gradient)
width_palette <- setNames(
  viridisLite::viridis(n_bins, option = "C")[floor(t_vals * (n_bins - 1)) + 1],
  width_levels
)

# (Optional) quick check of first few colors
head(width_palette, 6)

# Histogram ----------------------------------------
gg_width_hist_labeled <- ggplot(
  pixel_width_df,
  aes(x = WettedWidth_m, fill = WettedWidth_class)
) +
  geom_histogram(
    binwidth = 2.5,
    boundary = 0,
    color = "black",
    na.rm = TRUE
  ) +
  stat_bin(
    data = pixel_width_df,        
    aes(
      x = WettedWidth_m,
      label = paste0(
        "[",
        round(after_stat(xmin), 1), "–",
        round(after_stat(xmax), 1), " ]"
      )
    ),
    binwidth = 2.5,
    boundary = 0,
    geom = "text",
    inherit.aes = FALSE,           
    vjust = -0.4,
    size = 3,
    na.rm = TRUE
  ) +
  scale_fill_manual(
    values = width_palette,
    na.value = "gray60",
    name = "Wetted width class (m)"
  ) +
  labs(
    title = "Distribution of Wetted Stream Widths",
    subtitle = "2.5m size classes",
    x = "Wetted width",
    y = "Pixel count"
  ) +
  theme_minimal(base_size = 14, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "none"    # change to "bottom" to add legend.
  )

print(gg_width_hist_labeled)


### Map pixels, color-coded by wetted width --------------
# Join width classes to geometry (one row per pixel)
pixel_width_sf <- pixel_width_df %>%
  left_join(
    pixel_geom_sf %>% dplyr::select(pixel_ID, geometry),
    by = "pixel_ID"
  ) %>%
  sf::st_as_sf(crs = sf::st_crs(pixel_geom_sf)) %>%
  sf::st_make_valid()

# Counts for subtitle
n_pix_total <- nrow(pixel_width_sf)
n_pix_na    <- sum(is.na(pixel_width_sf$WettedWidth_class))
pct_na      <- round(100 * n_pix_na / n_pix_total, 2)

gg_all_pixels_width <- ggplot() +
  geom_sf(data = huc12_bg, fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_sf(
    data = pixel_width_sf,
    aes(color = WettedWidth_class),
    size = 1.1,
    alpha = 0.75
  ) +
  scale_color_manual(
    values = width_palette,
    na.value = "gray60",
    name = "Wetted width class"
  ) +
  coord_sf(expand = FALSE) +
  labs(
    title = "BMWA Pixels Colored by Wetted Width Class",
    subtitle = paste0(
      "Pixels = ", n_pix_total,
      "; NA width = ", n_pix_na, " (", pct_na, "%)"
    ),
    caption = "Width classes identical to histogram; NA shown in gray."
  ) +
  theme_minimal(base_size = 14, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12),
    plot.caption = element_text(size = 10),
    legend.position = "bottom"
  )

print(gg_all_pixels_width)
################################################################################

### What is the variation in pre-fire NDVI?? (Histograms) 
## Subtitle on each: 
sub_text <- "Solid line is the mean; dashed line is the median."
## ===============================
## Histogram 1: Prefire 3-year mean NDVI
## ===============================
p_avg <- ggplot(model_df, aes(x = ndvi_prefire_3yr_avg)) +
  geom_histogram(bins = 40, color = "black", fill = "gray70") +
  geom_vline(aes(xintercept = mean(ndvi_prefire_3yr_avg, na.rm = TRUE)),
             color = "black", linewidth = 0.9) +
  geom_vline(aes(xintercept = median(ndvi_prefire_3yr_avg, na.rm = TRUE)),
             color = "black", linetype = "dashed", linewidth = 0.9) +
  theme_minimal(base_family = "Times New Roman", base_size = 12) +
  labs(
    title    = "Distribution of 3-Year Prefire Mean NDVI (per pixel)",
    subtitle = sub_text,
    x        = "Prefire NDVI (3-year mean)",
    y        = "Pixel count"
  )

print(p_avg)

## ===============================
## Histogram 2: Prefire 3-year median NDVI
## ===============================
p_med <- ggplot(model_df, aes(x = ndvi_prefire_3yr_med)) +
  geom_histogram(bins = 40, color = "black", fill = "gray70") +
  geom_vline(aes(xintercept = mean(ndvi_prefire_3yr_med, na.rm = TRUE)),
             color = "black", linewidth = 0.9) +
  geom_vline(aes(xintercept = median(ndvi_prefire_3yr_med, na.rm = TRUE)),
             color = "black", linetype = "dashed", linewidth = 0.9) +
  theme_minimal(base_family = "Times New Roman", base_size = 12) +
  labs(
    title    = "Distribution of 3-Year Prefire Median NDVI (per pixel)",
    subtitle = sub_text,
    x        = "Prefire NDVI (3-year median)",
    y        = "Pixel count"
  )

print(p_med)

## ===============================
## Histogram 3: Prefire 3-year NDVI range
## ===============================
p_range <- ggplot(model_df, aes(x = ndvi_prefire_3yr_range)) +
  geom_histogram(bins = 40, color = "black", fill = "gray70") +
  geom_vline(aes(xintercept = mean(ndvi_prefire_3yr_range, na.rm = TRUE)),
             color = "black", linewidth = 0.9) +
  geom_vline(aes(xintercept = median(ndvi_prefire_3yr_range, na.rm = TRUE)),
             color = "black", linetype = "dashed", linewidth = 0.9) +
  theme_minimal(base_family = "Times New Roman", base_size = 12) +
  labs(
    title    = "Distribution of 3-Year Prefire NDVI Range (per pixel)",
    subtitle = sub_text,
    x        = "Prefire NDVI range (max - min over 3 years)",
    y        = "Pixel count"
  )

print(p_range)

## ===============================
## Histogram 4: Prefire 3-year NDVI minimum
## ===============================
p_min <- ggplot(model_df, aes(x = ndvi_prefire_3yr_min)) +
  geom_histogram(bins = 40, color = "black", fill = "gray70") +
  geom_vline(aes(xintercept = mean(ndvi_prefire_3yr_min, na.rm = TRUE)),
             color = "black", linewidth = 0.9) +
  geom_vline(aes(xintercept = median(ndvi_prefire_3yr_min, na.rm = TRUE)),
             color = "black", linetype = "dashed", linewidth = 0.9) +
  theme_minimal(base_family = "Times New Roman", base_size = 12) +
  labs(
    title    = "Distribution of 3-Year Prefire NDVI Minimum (per pixel)",
    subtitle = sub_text,
    x        = "Prefire NDVI min (over 3 years)",
    y        = "Pixel count"
  )

print(p_min)

################################################################################


################################################################################
## Identify high-range pixels (prefire NDVI range >= 0.40) and make scatterplots
################################################################################

# Threshold for "high" prefire NDVI range
range_thr <- 0.40

# Collapse to one unique record per pixel_ID
pixel_prefire_df <- model_df %>%
  dplyr::filter(!is.na(ndvi_prefire_3yr_range)) %>%
  dplyr::arrange(pixel_ID, year) %>%   # ensure deterministic selection
  dplyr::group_by(pixel_ID) %>%
  dplyr::slice(1) %>%  # take the first row for each pixel_ID
  dplyr::ungroup() %>%
  dplyr::select(
    pixel_ID,
    sev_num,
    ndvi_prefire_3yr_min,
    ndvi_prefire_3yr_avg,
    ndvi_prefire_3yr_range,
    cwd_5yr_postfire_avg_08,
    tmeanz_JJA,
    twi
  ) %>%
  dplyr::mutate(
    high_range = ndvi_prefire_3yr_range >= range_thr
)

# Quick counts
n_pix_total <- nrow(pixel_prefire_df)
n_pix_high  <- sum(pixel_prefire_df$high_range, na.rm = TRUE)
message("Pixels with prefire NDVI range ≥ ", range_thr, ": ",
        n_pix_high, " / ", n_pix_total, " (",
        round(100 * n_pix_high / n_pix_total, 2), "% )")

# Common theme
scatter_theme <- theme_minimal(base_family = "Times New Roman", base_size = 12)

# Color mapping
range_cols <- c(`FALSE` = "gray70", `TRUE` = "darkred")

## Scatterplot 3: y = ndvi_prefire_3yr_range, x = sev_num
## ---------------------------------------------------------------------------
p_sc3 <- ggplot(
  pixel_prefire_df,
  aes(x = sev_num, y = ndvi_prefire_3yr_range, color = high_range)
) +
  geom_point(alpha = 0.7) +
  scale_color_manual(
    values = range_cols,
    name   = paste0("Prefire NDVI range ≥ ", range_thr)
  ) +
  labs(
    title = "Prefire NDVI Range vs Fire Severity",
    x     = "Fire severity (sev_num)",
    y     = "Prefire NDVI range (max − min over 3 years)"
  ) +
  scatter_theme

print(p_sc3)

## ---------------------------------------------------------------------------
## Scatterplot 4: y = ndvi_prefire_3yr_range, x = cwd_5yr_postfire_avg_08
## ---------------------------------------------------------------------------
p_sc4 <- pixel_prefire_df %>%
  dplyr::filter(!is.na(cwd_5yr_postfire_avg_08)) %>%
  ggplot(
    aes(x = cwd_5yr_postfire_avg_08, y = ndvi_prefire_3yr_range, color = high_range)
  ) +
  geom_point(alpha = 0.7) +
  scale_color_manual(
    values = range_cols,
    name   = paste0("Prefire NDVI range ≥ ", range_thr)
  ) +
  labs(
    title = "Prefire NDVI Range vs Postfire CWD (5-year mean)",
    x     = "CWD (5-year postfire average, August; standardized)",
    y     = "Prefire NDVI range (max − min over 3 years)"
  ) +
  scatter_theme

print(p_sc4)

## ---------------------------------------------------------------------------
## Scatterplot 5: y = ndvi_prefire_3yr_range, x = tmeanz_JJA
## ---------------------------------------------------------------------------
p_sc5 <- pixel_prefire_df %>%
  dplyr::filter(!is.na(tmeanz_JJA)) %>%
  ggplot(
    aes(x = tmeanz_JJA, y = ndvi_prefire_3yr_range, color = high_range)
  ) +
  geom_point(alpha = 0.7) +
  scale_color_manual(
    values = range_cols,
    name   = paste0("Prefire NDVI range ≥ ", range_thr)
  ) +
  labs(
    title = "Prefire NDVI Range vs Summer Temperature (JJA)",
    x     = "Mean summer temperature anomaly (tmeanz_JJA)",
    y     = "Prefire NDVI range (max − min over 3 years)"
  ) +
  scatter_theme

print(p_sc5)

## ---------------------------------------------------------------------------
## Scatterplot 6: y = ndvi_prefire_3yr_range, x = min
## ---------------------------------------------------------------------------
p_sc6 <- pixel_prefire_df %>%
  dplyr::filter(!is.na(ndvi_prefire_3yr_min)) %>%
  ggplot(
    aes(x = ndvi_prefire_3yr_min, y = ndvi_prefire_3yr_range, color = high_range)
  ) +
  geom_point(alpha = 0.7) +
  scale_color_manual(
    values = range_cols,
    name   = paste0("Prefire NDVI range ≥ ", range_thr)
  ) +
  labs(
    title = "Prefire NDVI Range vs Prefire min",
    x     = "Min",
    y     = "Prefire NDVI range (max − min over 3 years)"
  ) +
  scatter_theme

print(p_sc6)
# ----------------------------------------------------------------------------
## Map of high-range prefire NDVI pixels (range ≥ 0.40)
# ----------------------------------------------------------------------------
# Subset to high-range pixels (using pixel_prefire_df already created above)
high_range_ids <- pixel_prefire_df %>%
  dplyr::filter(high_range) %>%
  dplyr::pull(pixel_ID)

# Join to geometry
high_range_sf <- pixel_geom_sf %>%
  dplyr::filter(pixel_ID %in% high_range_ids) %>%
  dplyr::left_join(
    pixel_prefire_df %>%
      dplyr::select(pixel_ID, ndvi_prefire_3yr_range),
    by = "pixel_ID"
  )

# Quick count for subtitle
n_high_map <- nrow(high_range_sf)

# Map: HUC12 backdrop + streams + high-range pixels colored by prefire NDVI range ----
gg_high_range <- ggplot() +
  geom_sf(data = huc12_bg,
          fill = "grey95",
          color = "grey70",
          linewidth = 0.2) +
  #geom_sf(data = streams_bg,   # uncomment to add stream layer. 
  #        color = "darkblue",
  #        linewidth = 0.5) +
  geom_sf(data = high_range_sf,
          aes(color = ndvi_prefire_3yr_range),
          size = 1.4,
          alpha = 0.8) +
  coord_sf(expand = FALSE) +
  scale_color_viridis_c(
    name   = "Prefire NDVI range\n(max − min over 3 years)",
    option = "plasma"
  ) +
  labs(
    title    = "BMWA Pixels with High Prefire NDVI Variability",
    subtitle = paste0("Prefire NDVI range ≥ ", range_thr,
                      " (N = ", n_high_map, " pixels)"),
    caption  = "Based on 3-year prefire NDVI range"
  ) +
  theme_minimal(base_size = 14, base_family = "Times New Roman") +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12),
    plot.caption  = element_text(size = 10),
    legend.position = "right"
  )

print(gg_high_range)
### 
# ----------------------------------------------------------------------------
## Histogram 3 (Repeat): Prefire 3-year NDVI range ≥ 0.40
# Reduce this (instead of pixel-yrs, use ind. pixels)
# ----------------------------------------------------------------------------
p_range_thr_shaded_labeled <- ggplot(pixel_prefire_df, aes(x = ndvi_prefire_3yr_range)) +
  
  # Shaded region for mapped pixels (range ≥ 0.40)
  annotate(
    "rect",
    xmin = range_thr, xmax = Inf,
    ymin = -Inf, ymax = Inf,
    fill  = "darkred",
    alpha = 0.12
  ) +
  
  # Histogram (counts unique pixels)
  geom_histogram(
    bins  = 40,
    color = "black",
    fill  = "gray70"
  ) +
  
  # Bar-wise N labels (inside each bar; counts are per pixel)
  geom_text(
    stat = "bin",
    bins = 40,
    aes(label = after_stat(count)),
    vjust = 1.2,
    size  = 2.8,
    color = "black"
  ) +
  
  # Mean (per pixel)
  geom_vline(
    aes(xintercept = mean(ndvi_prefire_3yr_range, na.rm = TRUE)),
    color = "black", linewidth = 0.9
  ) +
  
  # Median (per pixel)
  geom_vline(
    aes(xintercept = median(ndvi_prefire_3yr_range, na.rm = TRUE)),
    color = "black", linetype = "dashed", linewidth = 0.9
  ) +
  
  # Explicit threshold line
  geom_vline(
    xintercept = range_thr,
    linetype   = "dotted",
    linewidth  = 1.1,
    color      = "darkred"
  ) +
  
  # Label for shaded region, anchored near x-axis
  annotate(
    "text",
    x     = range_thr,
    y     = 0,
    label = "Mapped pixels (≥ 0.40)",
    vjust = -0.6,
    angle = 90,
    size  = 3.2,
    color = "darkred"
  ) +
  
  theme_minimal(base_family = "Times New Roman", base_size = 12) +
  labs(
    title    = "Distribution of 3-Year Prefire NDVI Range (per pixel)",
    subtitle = sub_text,
    x        = "Prefire NDVI range (max − min over 3 years)",
    y        = "Pixel count"
  )

print(p_range_thr_shaded_labeled)
################################################################################

################################################################################
#######  =======     PULLING PIXELS FOR NDVI IMAGING OVER TIME (FIGURE S1)   ===
################################################################################
library(tidyr)
library(readr)

response_var <- "delta_ndvi_min"

targets <- tibble::tibble(
  target_label = c("dNDVI_neg020", "dNDVI_neg010", "dNDVI_neg005"),
  target_delta = c(-0.20, -0.10, -0.05)
)

bartlett_candidates <- model_df %>%
  filter(
    huc12 == "Bartlett Creek",
    years_since_fire >= 0,
    !is.na(.data[[response_var]]),
    !sev_group %in% c("Unburned", "Control")
  ) %>%
  select(
    pixel_ID, huc12, fire_name, fire_year, ref_year,
    year, years_since_fire,
    sev_num, sev_group, sev_class,
    ndvi_postfire,
    ndvi_prefire_3yr_min,
    delta_ndvi_min,
    delta_ndvi_avg,
    delta_ndvi_med
  )

pixel_target_matches <- bartlett_candidates %>%
  tidyr::crossing(targets) %>%
  mutate(abs_error = abs(.data[[response_var]] - target_delta)) %>%
  group_by(pixel_ID, target_label) %>%
  slice_min(abs_error, n = 1, with_ties = FALSE) %>%
  ungroup()

pixel_arc_candidates <- pixel_target_matches %>%
  select(
    pixel_ID, huc12, fire_name, fire_year, ref_year,
    sev_num, sev_group, sev_class,
    target_label, year, years_since_fire,
    ndvi_postfire,
    ndvi_prefire_3yr_min,
    observed_delta = all_of(response_var),
    abs_error
  ) %>%
  pivot_wider(
    id_cols = c(
      pixel_ID, huc12, fire_name, fire_year, ref_year,
      sev_num, sev_group, sev_class,
      ndvi_prefire_3yr_min
    ),
    names_from = target_label,
    values_from = c(
      year,
      years_since_fire,
      ndvi_postfire,
      observed_delta,
      abs_error
    )
  ) %>%
  mutate(
    ordered_recovery =
      years_since_fire_dNDVI_neg020 <= years_since_fire_dNDVI_neg010 &
      years_since_fire_dNDVI_neg010 <= years_since_fire_dNDVI_neg005,
    
    mean_abs_error = rowMeans(across(starts_with("abs_error_")), na.rm = TRUE),
    
    max_abs_error = pmax(
      abs_error_dNDVI_neg020,
      abs_error_dNDVI_neg010,
      abs_error_dNDVI_neg005,
      na.rm = TRUE
    )
  ) %>%
  arrange(desc(ordered_recovery), mean_abs_error, max_abs_error)

# QA: how many candidates follow the ideal chronological recovery order?
pixel_arc_candidates %>%
  count(ordered_recovery)

# Candidate summary: one row per candidate pixel
top_20_candidate_summary <- pixel_arc_candidates %>%
  filter(ordered_recovery) %>%
  slice_head(n = 20) %>%
  mutate(rank = row_number()) %>%
  select(
    rank,
    pixel_ID,
    huc12,
    fire_name,
    fire_year,
    ref_year,
    sev_num,
    sev_group,
    sev_class,
    ndvi_prefire_3yr_min,
    
    year_dNDVI_neg020,
    years_since_fire_dNDVI_neg020,
    ndvi_postfire_dNDVI_neg020,
    observed_delta_dNDVI_neg020,
    abs_error_dNDVI_neg020,
    
    year_dNDVI_neg010,
    years_since_fire_dNDVI_neg010,
    ndvi_postfire_dNDVI_neg010,
    observed_delta_dNDVI_neg010,
    abs_error_dNDVI_neg010,
    
    year_dNDVI_neg005,
    years_since_fire_dNDVI_neg005,
    ndvi_postfire_dNDVI_neg005,
    observed_delta_dNDVI_neg005,
    abs_error_dNDVI_neg005,
    
    ordered_recovery,
    mean_abs_error,
    max_abs_error
  )

top_20_pixels <- top_20_candidate_summary %>%
  pull(pixel_ID)

rank_lookup <- top_20_candidate_summary %>%
  select(pixel_ID, rank)

# Long table: all annual NDVI and deltaNDVI values for selected pixels
top_20_delta_ndvi_timeseries <- model_df %>%
  filter(
    pixel_ID %in% top_20_pixels,
    huc12 == "Bartlett Creek"
  ) %>%
  left_join(rank_lookup, by = "pixel_ID") %>%
  select(
    rank,
    pixel_ID,
    huc12,
    fire_name,
    fire_year,
    ref_year,
    year,
    years_since_fire,
    sev_num,
    sev_group,
    sev_class,
    ndvi_postfire,
    ndvi_prefire_3yr_min,
    delta_ndvi_min,
    delta_ndvi_avg,
    delta_ndvi_med
  ) %>%
  arrange(rank, year)

# Wide table: raw NDVI values shown as NDVI_YYYY columns
top_20_ndvi_wide <- top_20_delta_ndvi_timeseries %>%
  select(
    pixel_ID,
    rank,
    huc12,
    fire_name,
    fire_year,
    ref_year,
    sev_group,
    sev_num,
    year,
    ndvi_postfire,
    ndvi_prefire_3yr_min
  ) %>%
  mutate(year_col = paste0("NDVI_", year)) %>%
  select(-year) %>%
  pivot_wider(
    names_from = year_col,
    values_from = ndvi_postfire
  )

# View outputs
print(top_20_candidate_summary)
View(top_20_candidate_summary)
View(top_20_delta_ndvi_timeseries)
View(top_20_ndvi_wide)
################################################################################

################################################################################
####### Stats & Histograms: all potential predictors and NDVI ------------------
### All sev_groups, all pixel-years: -----------------
# Columns to summarize
vars_to_summarize <- c(
  #"ppttot_JJA", 
  "pptz_JJA",
  #"tmaxmean_JJA", 
  #"tmaxz_JJA",
  #"tmeanavg_JJA", 
  "tmeanz_JJA",
  #"swe_peak_MAM", 
  #"SWE_Apr", 
  "swez_Apr",
  #"aetavg_JJA",
  "twi", 
  #"veg_climate_index_08", 
  "cwd_5yr_zscore_08",
  #"ndvi_prefire_3yr_min", 
  #"ndvi_postfire", 
  "delta_ndvi_min"
)

# Keep only those that exist
present_vars <- intersect(vars_to_summarize, names(model_df))
missing_vars <- setdiff(vars_to_summarize, present_vars)
if (length(missing_vars)) {
  message("Missing in model_df (skipped): ", paste(missing_vars, collapse = ", "))
}

# Summarize -> tidy table (min, max, mean, sd) without name collisions
summary_stats <- model_df %>%
  select(all_of(present_vars)) %>%
  summarise(across(
    everything(),
    list(
      min  = ~min(.x, na.rm = TRUE),
      max  = ~max(.x, na.rm = TRUE),
      mean = ~mean(.x, na.rm = TRUE),
      sd   = ~sd(.x, na.rm = TRUE)
    ),
    .names = "{.col}__{.fn}"   # <- use double underscore to avoid clashes
  )) %>%
  pivot_longer(
    everything(),
    names_to = c("variable", "stat"),
    names_sep = "__",          # <- split on the double underscore only
    values_to = "value"
  ) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  arrange(variable)

print(summary_stats, n = Inf)

## Plot histograms:
vars_to_plot <- c(
  "ppttot_JJA", "pptz_JJA",
  "tmaxmean_JJA", "tmaxz_JJA",
  "tmeanavg_JJA", "tmeanz_JJA",
  #"swe_peak_MAM", 
  "SWE_Apr", "swez_Apr",
  #"aetavg_JJA",
  "twi", "veg_climate_index_08", "cwd_5yr_zscore_08",
  "ndvi_prefire_3yr_min", "ndvi_postfire", "delta_ndvi_min"
)

# Keep only those that are present
present_vars <- intersect(vars_to_plot, names(model_df))
missing_vars <- setdiff(vars_to_plot, present_vars)
if (length(missing_vars)) {
  message("Missing in model_df (skipped): ", paste(missing_vars, collapse = ", "))
}

# Long format for plotting
plot_df <- model_df %>%
  select(all_of(present_vars)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value")

# Per-variable stats for overlays
stat_df <- plot_df %>%
  group_by(variable) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    m_minus_sd = mean - sd,
    m_plus_sd  = mean + sd
  )

# Plot: histograms with mean and ±1 SD lines
ggplot(plot_df, aes(x = value)) +
  geom_histogram(bins = 30, na.rm = TRUE, fill = "grey80", color = "black") +
  # mean line
  geom_vline(data = stat_df, aes(xintercept = mean), linewidth = 0.6) +
  # ±1 SD lines
  geom_vline(data = stat_df, aes(xintercept = m_minus_sd), linetype = "dashed", linewidth = 0.5) +
  geom_vline(data = stat_df, aes(xintercept = m_plus_sd),  linetype = "dashed", linewidth = 0.5) +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  labs(
    title = "Distributions of covariates and response variable (All sev_groups)",
    subtitle = "Vertical lines: mean (solid) and mean ± 1 SD (dashed)",
    x = "Value", y = "Count"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())
### UNBURNED ONLY, all pixel-years: --------
### ---- Filter to only Unburned pixels 
unburned_df <- model_df %>%
  filter(sev_group == "Unburned")

### ---- Columns to summarize 
vars_to_summarize <- c(
  "ppttot_JJA", "pptz_JJA",
  "tmaxmean_JJA", "tmaxz_JJA",
  "tmeanavg_JJA", "tmeanz_JJA",
  #"swe_peak_MAM", 
  "SWE_Apr", "swez_Apr",
  "aetavg_JJA",
  "twi", "veg_climate_index_08", "cwd_5yr_zscore_08",
  "ndvi_prefire_3yr_avg", "ndvi_postfire", "delta_ndvi_min"
)

# Keep only those that exist
present_vars <- intersect(vars_to_summarize, names(unburned_df))
missing_vars <- setdiff(vars_to_summarize, present_vars)
if (length(missing_vars)) {
  message("Missing in unburned_df (skipped): ", paste(missing_vars, collapse = ", "))
}

### ---- Summary stats 
summary_stats <- unburned_df %>%
  select(all_of(present_vars)) %>%
  summarise(across(
    everything(),
    list(
      min  = ~min(.x, na.rm = TRUE),
      max  = ~max(.x, na.rm = TRUE),
      mean = ~mean(.x, na.rm = TRUE),
      sd   = ~sd(.x, na.rm = TRUE)
    ),
    .names = "{.col}__{.fn}"
  )) %>%
  pivot_longer(
    everything(),
    names_to = c("variable", "stat"),
    names_sep = "__",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  arrange(variable)

print(summary_stats, n = Inf)

### ---- Histograms 
vars_to_plot <- vars_to_summarize

# Keep only those present
present_vars <- intersect(vars_to_plot, names(unburned_df))
missing_vars <- setdiff(vars_to_plot, present_vars)
if (length(missing_vars)) {
  message("Missing in unburned_df (skipped): ", paste(missing_vars, collapse = ", "))
}

# Long format
plot_df <- unburned_df %>%
  select(all_of(present_vars)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value")

# Stats for overlays
stat_df <- plot_df %>%
  group_by(variable) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    m_minus_sd = mean - sd,
    m_plus_sd  = mean + sd
  )

# Plot
ggplot(plot_df, aes(x = value)) +
  geom_histogram(bins = 30, na.rm = TRUE, fill = "grey80", color = "black") +
  geom_vline(data = stat_df, aes(xintercept = mean), linewidth = 0.6) +
  geom_vline(data = stat_df, aes(xintercept = m_minus_sd), linetype = "dashed", linewidth = 0.5) +
  geom_vline(data = stat_df, aes(xintercept = m_plus_sd),  linetype = "dashed", linewidth = 0.5) +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  labs(
    title = "Distributions of covariates and response variable (Unburned pixels only)",
    subtitle = "Vertical lines: mean (solid) and mean ± 1 SD (dashed)",
    x = "Value", y = "Count"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())
################################################################################

###  Rename dataset to data_long_complete (post-pivoting) -------
data_long_complete <- model_df  
rm(model_df)

total_rows <- dplyr::n_distinct(data_long_complete$pixel_ID)
print(total_rows) ## 10185
# --------------------------------------------------------------

# Test a subset of Pixels to check climate variables: --------------------------
# pick 5 random pixel IDs present in the table
pix_take <- data_long_complete %>%
  distinct(pixel_ID) %>%
  pull(pixel_ID) %>%
  { sample(., size = min(5L, length(.))) }

# columns to keep (long, year-aligned versions)
keep_cols <- c(
  "pixel_ID", "year", "fire_name", "sev_group", "ref_year", "huc12",
  "delta_ndvi_min",
  "ppttot_JJA", "pptz_JJA",
  "tmaxmean_JJA", "tmaxz_JJA",
  "tmeanavg_JJA", "tmeanz_JJA",
  "swe_peak_MAM", 
  "SWE_Apr", 
  "swez_Apr"
)

# build the subset
test_subset <- data_long_complete %>%
  filter(pixel_ID %in% pix_take) %>%
  arrange(pixel_ID, year) %>%
  select(any_of(keep_cols))
################################################################################


### ===== QA check: data_long_complete ======= ######
##### PART 1. CHECK FOR INVALID VALUES  #############
# ---- Set target YSF 
target_ysf <- 5

# ---- List of predictors to check 
predictors <- c(
  "ndvi_postfire",
  "ndvi_prefire_3yr_avg",
  "ndvi_prefire_3yr_min",
  "ndvi_prefire_3yr_range",
  "ndvi_prefire_3yr_med",
  "delta_ndvi_min",
  "delta_ndvi_min",
  "delta_ndvi_med",
  "ppttot_JJA", 
  "pptz_JJA",
  "tmeanz_JJA", 
  "tmaxz_JJA", 
  "tmeanavg_JJA", 
  "tmaxmean_JJA",
  #"swe_peak_MAM", 
  "SWE_Apr",
  "swez_Apr",
  "aetavg_JJA",
  "veg_climate_index_08", 
  "cwd_5yr_zscore_08",
  "twi", 
  "sev_num"#, 
  #"fire_size_ha"
)

# ---- Filter to target YSF 
df_check <- data_long_complete %>%
  filter(years_since_fire == target_ysf) %>%
  st_drop_geometry()

# ---- QA summary per predictor 
qa_tbl <- map_dfr(predictors, function(var) {
  x <- df_check[[var]]
  
  tibble(
    variable    = var,
    type        = class(x)[1],
    n_total     = length(x),
    n_valid     = sum(is.finite(suppressWarnings(as.numeric(x)))),
    n_NA        = sum(is.na(x)),
    n_nonfinite = sum(!is.na(x) & !is.finite(suppressWarnings(as.numeric(x))))
  )
})

# ---- Show table 
qa_tbl %>%
  rename(
    "Variable"         = variable,
    "Type"             = type,
    "Total Rows"       = n_total,
    "Valid (finite)"   = n_valid,
    "NA Values"        = n_NA,
    "Non-finite/Other" = n_nonfinite
  ) %>%
  knitr::kable(caption = "QA Summary of Predictor Variables @ targetYSF")
## ---------------------------------------------------
##### PART 2. MAP SPATIAL RELATIONSHIPS  #############

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(sf)
})

# ---- Predictors to color by 
predictors <- c(
  "ndvi_postfire",
  "ndvi_prefire_3yr_avg",
  "delta_ndvi_min",
  "huc12",
  #"ppttot_JJA", 
  "pptz_JJA",
  "tmeanz_JJA", 
  "tmaxz_JJA", 
  #"tmeanavg_JJA", 
  #"tmaxmean_JJA",
  #"swe_peak_MAM", 
  #"SWE_Apr",
  "swez_Apr",
  "veg_climate_index_08", 
  "cwd_5yr_zscore_08",
  "twi", 
  "sev_num",
  "fire_size_ha"
)

# ============================
# ==== USER PARAMETERS =====
# ============================
target_ysf       <- 5
predictor        <- "fire_size_ha"   # must be one of `predictors`
sev_filter       <- "All"     # "All" or any of: "Unburned","Low","Moderate","High"
point_size       <- 0.8
sample_frac_pre  <- NULL      # e.g., NULL or decimal (0.1, 0.9) to pre-sample BEFORE GeoJSON parse
sample_frac_plot <- 1.0       # e.g., reduce to 0.5 to downsample (just for the map)

# ===============================
# ====== HELPER FUNCTIONS  ======
# ============================
# --- Function: convert .geo -> sf (general) 
as_sf_from_geo <- function(df, crs_guess = 4326) {
  if (inherits(df, "sf") && "geometry" %in% names(df)) return(df)
  stopifnot(".geo" %in% names(df))
  geom <- sf::st_as_sfc(df$.geo, GeoJSON = TRUE)
  out  <- df
  out$geometry <- geom
  out  <- sf::st_as_sf(out, sf_column_name = "geometry")
  if (is.na(sf::st_crs(out))) sf::st_crs(out) <- sf::st_crs(crs_guess)
  dplyr::filter(out, !sf::st_is_empty(geometry))
}

# --- Function: convert .geo -> sf for a subset 
as_sf_from_geo_subset <- function(df_subset, crs_guess = 4326) {
  if (inherits(df_subset, "sf") && "geometry" %in% names(df_subset)) return(df_subset)
  stopifnot(".geo" %in% names(df_subset))
  geom <- sf::st_as_sfc(df_subset$.geo, GeoJSON = TRUE)
  out  <- df_subset
  out$geometry <- geom
  out  <- sf::st_as_sf(out, sf_column_name = "geometry")
  if (is.na(sf::st_crs(out))) sf::st_crs(out) <- sf::st_crs(crs_guess)
  dplyr::filter(out, !sf::st_is_empty(geometry))
}

# --- Function: build reusable XY table for a single YSF 
# Returns plain data.frame with x, y, sev_group, years_since_fire, and requested columns
prep_points_for_ysf <- function(df_long_complete,
                                target_YSF,
                                keep_cols,
                                crs_guess = 4326,
                                to_crs = 32612,
                                sample_frac_pre = NULL) {
  stopifnot("years_since_fire" %in% names(df_long_complete))
  
  # Filter to YSF BEFORE any geometry work
  df_sub <- df_long_complete[df_long_complete$years_since_fire == target_YSF, , drop = FALSE]
  if (nrow(df_sub) == 0) stop(sprintf("No rows with years_since_fire == %s.", target_YSF))
  
  # Optional pre-sample to reduce GeoJSON parsing cost
  if (!is.null(sample_frac_pre)) {
    stopifnot(sample_frac_pre > 0, sample_frac_pre <= 1)
    set.seed(1L)
    df_sub <- dplyr::sample_frac(df_sub, sample_frac_pre)
  }
  
  # Keep only needed columns (+ .geo for geometry creation)
  cols_needed <- unique(c("sev_group", "years_since_fire", keep_cols, ".geo"))
  cols_needed <- intersect(cols_needed, names(df_sub))
  df_sub <- df_sub[, cols_needed, drop = FALSE]
  
  # Convert this subset to sf
  sf_sub  <- as_sf_from_geo_subset(df_sub, crs_guess = crs_guess)
  sf_proj <- suppressMessages(sf::st_transform(sf_sub, to_crs))
  
  # Centroids -> XY
  cents <- sf::st_centroid(sf_proj)
  xy    <- sf::st_coordinates(cents)
  
  # Return plain data.frame (fast to re-use/recolor)
  out <- cbind(sf::st_drop_geometry(sf_proj), x = xy[,1], y = xy[,2])
  as.data.frame(out)
}

# --- Plotting function: recolor pixels by any numeric predictor 
plot_prepped_points <- function(pts_ysf,
                                var,
                                sev_filter = "All",
                                point_size = 0.5,
                                sample_frac_plot = NULL) {
  if (!all(c("x","y") %in% names(pts_ysf))) {
    stop("Input must come from prep_points_for_ysf (missing x/y).")
  }
  if (!var %in% names(pts_ysf)) {
    stop(sprintf("Column '%s' not found in prepped data.", var))
  }
  
  df <- pts_ysf
  
  # Severity subset (case-insensitive)
  if (!identical(sev_filter, "All")) {
    if (!"sev_group" %in% names(df)) stop("'sev_group' not in prepped data.")
    keep <- tolower(df$sev_group) %in% tolower(sev_filter)
    df <- df[keep, , drop = FALSE]
  }
  
  if (nrow(df) == 0) {
    stop(sprintf("No rows to plot after sev_filter = %s",
                 if (identical(sev_filter, "All")) "All" else paste(sev_filter, collapse = ", ")))
  }
  
  # Optional light sampling for draw-time speed
  if (!is.null(sample_frac_plot)) {
    stopifnot(sample_frac_plot > 0, sample_frac_plot <= 1)
    set.seed(1L)
    df <- dplyr::sample_frac(df, sample_frac_plot)
  }
  
  subtitle_txt <- if (identical(sev_filter, "All")) {
    "All severity groups"
  } else {
    paste0("sev_group: ", paste(sev_filter, collapse = ", "))
  }
  
  ggplot(df, aes(x = x, y = y, color = .data[[var]])) +
    geom_point(size = point_size, stroke = 0, shape = 16, alpha = 0.9, na.rm = TRUE) +
    scale_color_viridis_c(na.value = "transparent") +
    guides(color = guide_colorbar(title = var)) +
    labs(
      title = sprintf("YSF = %s • %s", unique(df$years_since_fire), var),
      subtitle = subtitle_txt,
      x = NULL, y = NULL
    ) +
    coord_equal() +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank())
}
# ===============================
# ===============================
# ==== SET PLOT PARAMETERS  =====
# ===============================
# PREP ONCE PER YSF (keep sev_group + all predictors for fast recoloring)
pts_ysf <- prep_points_for_ysf(
  df_long_complete = data_long_complete,
  target_YSF       = target_ysf,
  keep_cols        = c("sev_group", predictors),
  sample_frac_pre  = sample_frac_pre
)

# PLOT (recolor instantly by changing `predictor` or `sev_filter`)
p <- plot_prepped_points(
  pts_ysf          = pts_ysf,
  var              = predictor,          
  sev_filter       = sev_filter,  # "All" or c("Unburned") or c("High","Moderate")
  point_size       = point_size,
  sample_frac_plot = sample_frac_plot
)
print(p)

# ---- List of variables to copy/paste for maps: ------
#  "ndvi_postfire",
#  "ndvi_prefire_3yr_avg",
#  "delta_ndvi_min",
#  "huc12",
#  "ppttot_JJA", 
#  "pptz_JJA",
#  "tmeanz_JJA", 
#  "tmaxz_JJA", 
#  "tmeanavg_JJA", 
#  "tmaxmean_JJA",
#  "swe_peak_MAM", 
#  "SWE_Apr",
#  "aetavg_JJA",
#  "veg_climate_index_08", 
#  "cwd_5yr_zscore_08",
#  "twi", 
#  "sev_num"
# ,"fire_size_ha"
## ----Plot different variable / sev_group / YSF combos: ----------------------
print(plot_prepped_points(pts_ysf, 
                          var = "delta_ndvi_min", 
                          sev_filter = "All"))

print(plot_prepped_points(pts_ysf, 
                          var = "cwd_5yr_zscore_08", 
                          sev_filter = c("Unburned")))

print(plot_prepped_points(pts_ysf, 
                          var = "tmeanz_JJA", 
                          sev_filter = c("Unburned")))

# Optional: add huc12 layer to the map and re-print.
huc12_32612 <- sf::st_transform(huc12_sf, 32612)  # match points' CRS
p <- (plot_prepped_points(pts_ysf, 
                          var = "fire_size_ha", 
                          sev_filter = "High"))
p <- p +
  geom_sf(data = huc12_32612, fill = NA, color = "grey30", linewidth = 0.2, inherit.aes = FALSE) +
  labs(
    title = "Fire Size, High Severity",
    subtitle = paste0("YSF = ", unique(pts_ysf$years_since_fire), " • sev_filter: ", paste(sev_filter, collapse = ", "))
  )

print(p)


# ---- FILTER TO SPECIFIC HUC12: -----------------------------------------------

# Names to keep (these must match the polygon 'name' field)
keep_huc12_names <- c(
  "Bartlett Creek",
  "Lower Gordon Creek"
)

# Polygon name column (from names(huc12_sf))
poly_name_col <- "name"

# 1) Filter polygons to selected names (reproject to match points' CRS)
huc12_sel <- sf::st_transform(huc12_sf, 32612) |>
  dplyr::filter(.data[[poly_name_col]] %in% keep_huc12_names)

if (nrow(huc12_sel) == 0) {
  stop("No polygons matched the provided names in column '", poly_name_col, "'.")
}

# 2) Filter points by those same watershed names (points use huc12 column)
pts_sel <- dplyr::filter(pts_ysf, huc12 %in% keep_huc12_names)

# 3) Plot
p_sel <- plot_prepped_points(
  pts_ysf          = pts_sel,
  var              = "sev_num",
  sev_filter       = "All",
  point_size       = 2.0,
  sample_frac_plot = sample_frac_plot
) +
  geom_sf(data = huc12_sel, fill = NA, color = "grey30", linewidth = 0.2, inherit.aes = FALSE) +
  labs(
    title = "NDVI Recovery (Selected HUC12 Watersheds)",
    subtitle = paste0(
      "YSF = ", paste(sort(unique(pts_sel$years_since_fire)), collapse = ", "),
      " • Watersheds: ", paste(keep_huc12_names, collapse = ", ")
    )
  )

print(p_sel)

# ============================
# ==== EXPORT FILTERED POINTS
# ============================

suppressPackageStartupMessages({ library(dplyr); library(sf) })

# 1) Turn filtered XY into sf and add lon/lat in WGS84
pts_sel_sf <- sf::st_as_sf(pts_sel, coords = c("x","y"), crs = 32612, remove = FALSE) %>%
  sf::st_transform(4326)

coords <- sf::st_coordinates(pts_sel_sf)
pts_sel_sf$longitude <- coords[,1]
pts_sel_sf$latitude  <- coords[,2]

# 2) Build the export table with requested fields
export_df <- pts_sel_sf %>%
  sf::st_drop_geometry() %>%
  transmute(
    latitude,
    longitude,
    huc12,                      # watershed name in points
    #ref_year,
    sev_group,
    years_since_fire,
    delta_ndvi_min,
    ndvi_prefire_3yr_avg
  )

# 3) Write CSV (in the current working directory)
out_csv <- file.path(getwd(), paste0("filtered_points_", format(Sys.Date(), "%Y%m%d"), ".csv"))
utils::write.csv(export_df, out_csv, row.names = FALSE)

message("Wrote: ", out_csv)

## --- PLOT/EXTRACT PIXELS FOR ENTIRE STUDY AREA: ---------------------------

# ==== PLOT (no HUC12 filter)
# Reproject polygons to match points' CRS (assumes pts_ysf is in EPSG:32612)
huc12_utm <- sf::st_transform(huc12_sf, 32612)

# Use ALL points (no filtering)
pts_sel <- pts_ysf

p_sel <- plot_prepped_points(
  pts_ysf          = pts_sel,
  var              = "sev_num",
  sev_filter       = "All",
  point_size       = 2.0,
  sample_frac_plot = sample_frac_plot
) +
  geom_sf(data = huc12_utm, fill = NA, color = "grey30", linewidth = 0.2, inherit.aes = FALSE) +
  labs(
    title = "NDVI Recovery (All HUC12 Watersheds)",
    subtitle = paste0("YSF = ", paste(sort(unique(pts_sel$years_since_fire)), collapse = ", "))
  )

print(p_sel)

# ============================
# ==== EXPORT ALL PIXELS FOR STUDY AREA (shpfile, CSV)
# ============================

suppressPackageStartupMessages({ library(dplyr); library(sf) })

# 1) Ensure sf with correct CRS; if x/y exist, construct from them
pts_sel_sf <- sf::st_as_sf(pts_ysf, coords = c("x","y"), crs = 32612, remove = FALSE) %>%
  sf::st_transform(4326)

# 2) Add lon/lat in WGS84
coords <- sf::st_coordinates(pts_sel_sf)
pts_sel_sf$longitude <- coords[,1]
pts_sel_sf$latitude  <- coords[,2]

# 3) Print out column names first
cat("Columns in exported shapefile:\n")
print(names(pts_sel_sf))

# 4) WRITE SHAPEFILE
out_dir  <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/ArcGIS/BMWC_Fire_Data/CH01_file_drop"
out_name <- "all_pixels_filtered_by_preFire_MIN.shp"

sf::st_write(
  pts_sel_sf,
  dsn = file.path(out_dir, out_name),
  delete_dsn = TRUE
)

cat("Shapefile written to:\n", file.path(out_dir, out_name), "\n")

# 5) "METADATA" OPTIONS FOR A SHAPEFILE
meta_txt <- c(
  "Dataset: pts_ysf_all_points",
  paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "CRS: WGS84",
  "Description: Exported pixels for entire BMWA, filtered by pre-fire minimum NDVI",
  "",
  "Fields:",
  paste0(" - ", names(pts_sel_sf), collapse = "\n")
)

writeLines(meta_txt, con = file.path(out_dir, "all_pixels_filtered_by_preFire_MIN_METADATA.txt"))

## EXPORT AS CSV FILE: 
# Build export table 
export_df <- pts_sel_sf %>%
  sf::st_drop_geometry() %>%
  transmute(
    latitude,
    longitude,
    huc12,
    sev_num,                 
    sev_group,
    years_since_fire,
    delta_ndvi_min,
    ndvi_prefire_3yr_avg
  )

# Write CSV
out_csv <- file.path(getwd(), paste0("all_points_", format(Sys.Date(), "%Y%m%d"), ".csv"))
utils::write.csv(export_df, out_csv, row.names = FALSE)
message("Wrote: ", out_csv)
# ------------------------------------------------------------------------------

### ===== Plotting covariates over time, 1984-2024: ----------------------------
# Make a copy 
df <- data_long_complete

# ---- 1) Identify columns 
ppt_cols   <- grep("^ppttot_JJA_\\d{4}", names(df), value = TRUE)
tmean_cols <- grep("^tmeanavg_JJA_\\d{4}", names(df), value = TRUE)

# SWE columns look like SWE_YYYYMM_TC; keep only April (MM == 04)
swe_cols   <- grep("^SWE_\\d{4}04_", names(df), value = TRUE)

if (length(swe_cols) == 0) {
  stop("No April SWE columns found. Check that names look like 'SWE_YYYY04_*'.")
}

# ---- 2) Wide -> long and summarize yearly means 
ppt_year <- df %>%
  select(all_of(ppt_cols)) %>%
  pivot_longer(everything(), names_to = "name", values_to = "ppttot_JJA") %>%
  mutate(Year = as.integer(str_extract(name, "\\d{4}"))) %>%
  group_by(Year) %>%
  summarize(ppttot_JJA = mean(ppttot_JJA, na.rm = TRUE), .groups = "drop")

tmean_year <- df %>%
  select(all_of(tmean_cols)) %>%
  pivot_longer(everything(), names_to = "name", values_to = "tmean_JJA") %>%
  mutate(Year = as.integer(str_extract(name, "\\d{4}"))) %>%
  group_by(Year) %>%
  summarize(tmean_JJA = mean(tmean_JJA, na.rm = TRUE), .groups = "drop")

swe_year <- df %>%
  select(all_of(swe_cols)) %>%
  pivot_longer(everything(), names_to = "name", values_to = "swe_Apr") %>%
  # SWE_YYYYMM_* : extract the first four digits as year
  mutate(Year = as.integer(str_extract(name, "^\\D*(\\d{4})"))) %>%
  group_by(Year) %>%
  summarize(swe_Apr = mean(swe_Apr, na.rm = TRUE), .groups = "drop")

# ---- 3) Merge to one data frame and complete year range 
years_full <- tibble(Year = 1984:2024)

merged <- years_full %>%
  left_join(ppt_year,  by = "Year") %>%
  left_join(tmean_year, by = "Year") %>%
  left_join(swe_year,   by = "Year")

# ---- 4) Normalize 0–1 (min–max) per variable 
normalize <- function(x) {
  if (all(is.na(x))) return(x)
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(x * 0)  # constant series -> all zeros
  (x - rng[1]) / (rng[2] - rng[1])
}

plot_df <- merged %>%
  mutate(
    ppttot_JJA_norm = normalize(ppttot_JJA),
    tmean_JJA_norm  = normalize(tmean_JJA),
    swe_Apr_norm    = normalize(swe_Apr)
  ) %>%
  select(Year, ppttot_JJA_norm, tmean_JJA_norm, swe_Apr_norm) %>%
  pivot_longer(-Year, names_to = "Variable", values_to = "Value") %>%
  mutate(Variable = recode(Variable,
                           ppttot_JJA_norm = "ppttot_JJA",
                           tmean_JJA_norm  = "tmean_JJA",
                           swe_Apr_norm    = "swe_Apr"))

# ---- 5) Plot 
ggplot(plot_df, aes(x = Year, y = Value, color = Variable)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = seq(1984, 2024, 4), limits = c(1984, 2024)) +
  labs(x = "Year", y = "Normalized (0–1)",
       title = "Yearly Mean (Normalized) — ppttot_JJA, tmean_JJA, SWE (April), 1984–2024") +
  theme_minimal() +
  theme(legend.title = element_blank())
################################################################################

# ======================================================================
# === Further examining spatial relationships using variable thresholds: 
# ======================================================================
## Load libraries -----
suppressPackageStartupMessages({
  library(dplyr); library(sf); library(ggplot2); library(rlang)
})

# ==== USER TOGGLES ==========
# ============================
target_ysf           <- c(5)    # one YSF (e.g., 15) or multiple (e.g., "c(10, 15, 20)")
var_for_color_scale  <- "delta_ndvi_min"          # variable to COLOR by on plot
sev_filter           <- "All"    # "All" or any of: "Unburned","Low","Moderate","High"
point_size           <- 0.8
sample_frac_pre      <- NULL                # e.g., 0.25 to pre-sample BEFORE GeoJSON parse
sample_frac_plot     <- NULL                # e.g., 0.5 to downsample for plotting only

# ---- Dynamic, optional thresholds 
resp_var    <- "delta_ndvi_min"   # response variable
resp_op     <- ">"            # one of: ">", ">=", "<", "<=", "==", "!="
resp_thresh <- NULL              # set NULL to disable threshold

pred_var    <- "tmaxz_JJA"     # variable to filter on (can equal var_for_color_scale)
pred_op     <- ">"           # one of: ">", ">=", "<", "<=", "==", "!=" Sets as "swez_Apr <= 0"
pred_thresh <- NULL              # set NULL to disable

# Toggle: HUC12 base map (set to NULL to skip) ------
huc12_path  <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/ArcGIS/BMWC_Fire_Data/HUC12_BMWA/HUC12_BobMarshall.shp"
# huc12_path <- NULL  #(uncomment this to remove HUC12s from map)

# ---- Function: eval generic threshold with operator string ----
apply_threshold <- function(df, var, op, thr) {
  if (is.null(thr) || is.null(var) || is.null(op)) return(df)
  if (!var %in% names(df)) {
    warning(sprintf("apply_threshold: variable '%s' not found in data; skipping this filter.", var))
    return(df)
  }
  op <- match.arg(op, c(">", ">=", "<", "<=", "==", "!="))
  expr <- switch(
    op,
    ">"  = expr(.data[[var]] >  !!thr),
    ">=" = expr(.data[[var]] >= !!thr),
    "<"  = expr(.data[[var]] <  !!thr),
    "<=" = expr(.data[[var]] <= !!thr),
    "==" = expr(.data[[var]] == !!thr),
    "!=" = expr(.data[[var]] != !!thr)
  )
  dplyr::filter(df, !!expr)
}

# ---- .geo -> sf for filtered subset ----
as_sf_from_geo_subset <- function(df_subset, crs_guess = 4326) {
  if (inherits(df_subset, "sf") && "geometry" %in% names(df_subset)) return(df_subset)
  stopifnot(".geo" %in% names(df_subset))
  geom <- sf::st_as_sfc(df_subset$.geo, GeoJSON = TRUE)
  out  <- df_subset
  out$geometry <- geom
  out  <- sf::st_as_sf(out, sf_column_name = "geometry")
  if (is.na(sf::st_crs(out))) sf::st_crs(out) <- sf::st_crs(crs_guess)
  dplyr::filter(out, !sf::st_is_empty(geometry))
}

# ---- PREP FUNCTION: one combined XY table (robust to missing vars; resolves toggles inside) ----
prep_points_multi_ysf <- function(df,
                                  ysf_vec          = NULL,
                                  keep_cols        = NULL,
                                  crs_guess        = 4326,
                                  to_crs           = 32612,
                                  sample_frac_pre  = NULL,
                                  resp_var         = NULL,
                                  resp_op          = NULL,
                                  resp_thresh      = NULL,
                                  pred_var         = NULL,
                                  pred_op          = NULL,
                                  pred_thresh      = NULL) {
  stopifnot("years_since_fire" %in% names(df))
  
  # Can assign toggles here, if not provided above
  if (is.null(ysf_vec))         ysf_vec         <- get0("target_ysf",        inherits = TRUE)
  if (is.null(sample_frac_pre)) sample_frac_pre <- get0("sample_frac_pre",   inherits = TRUE)
  if (is.null(resp_var))        resp_var        <- get0("resp_var",          inherits = TRUE)
  if (is.null(resp_op))         resp_op         <- get0("resp_op",           inherits = TRUE)
  if (is.null(resp_thresh))     resp_thresh     <- get0("resp_thresh",       inherits = TRUE)
  if (is.null(pred_var))        pred_var        <- get0("pred_var",          inherits = TRUE)
  if (is.null(pred_op))         pred_op         <- get0("pred_op",           inherits = TRUE)
  if (is.null(pred_thresh))     pred_thresh     <- get0("pred_thresh",       inherits = TRUE)
  if (is.null(keep_cols)) {
    v_for_color <- get0("var_for_color_scale", inherits = TRUE)
    keep_cols <- unique(c("sev_group", v_for_color, resp_var, pred_var))
  }
  
  if (is.null(ysf_vec)) stop("`target_ysf` is not set.")
  ysf_vec <- unique(ysf_vec)
  
  # Determine which requested columns actually exist
  requested_vars <- unique(c(resp_var, pred_var, keep_cols))
  vars_present   <- intersect(requested_vars, names(df))
  missing_vars   <- setdiff(requested_vars, vars_present)
  if (length(missing_vars)) {
    warning(sprintf(
      "prep_points_multi_ysf: dropping %d missing column(s): %s",
      length(missing_vars), paste(missing_vars, collapse = ", ")
    ))
  }
  
  # Keep only needed columns (+ .geo for geometry)
  cols_needed <- unique(c("sev_group", "years_since_fire", "year",
                          vars_present, ".geo"))
  cols_needed <- intersect(cols_needed, names(df))
  
  df_sub <- df %>%
    dplyr::filter(years_since_fire %in% ysf_vec) %>%
    dplyr::select(dplyr::all_of(cols_needed))
  
  # Apply thresholds BEFORE geometry work (only if the variable is present)
  if (!is.null(resp_thresh) && !is.null(resp_var) && resp_var %in% names(df_sub)) {
    df_sub <- apply_threshold(df_sub, resp_var, resp_op, resp_thresh)
  } else if (!is.null(resp_thresh) && !is.null(resp_var)) {
    warning(sprintf("Response variable '%s' not present after selection; skipping its threshold.", resp_var))
  }
  
  if (!is.null(pred_thresh) && !is.null(pred_var) && pred_var %in% names(df_sub)) {
    df_sub <- apply_threshold(df_sub, pred_var, pred_op, pred_thresh)
  } else if (!is.null(pred_thresh) && !is.null(pred_var)) {
    warning(sprintf("Predictor variable '%s' not present after selection; skipping its threshold.", pred_var))
  }
  
  if (!nrow(df_sub)) stop("No rows after YSF/threshold filtering.")
  
  # Optional pre-sampling (speed)
  if (!is.null(sample_frac_pre)) {
    stopifnot(sample_frac_pre > 0, sample_frac_pre <= 1)
    set.seed(1L)
    df_sub <- dplyr::sample_frac(df_sub, sample_frac_pre)
  }
  
  # GeoJSON -> sf -> project -> centroids -> XY
  sf_sub  <- as_sf_from_geo_subset(df_sub, crs_guess = crs_guess)
  sf_proj <- suppressMessages(sf::st_transform(sf_sub, to_crs))
  cents   <- sf::st_centroid(sf_proj)
  xy      <- sf::st_coordinates(cents)
  
  out <- cbind(sf::st_drop_geometry(sf_proj), x = xy[,1], y = xy[,2])
  as.data.frame(out)
}

# ---- PLOTTING FUNCTION: color by selected variable; HUC12 handled ONLY here ----
plot_points_multi <- function(pts_xy,
                              var_for_color_scale,   # REQUIRED (no lazy default)
                              sev_filter       = "All",
                              point_size       = 0.7,
                              sample_frac_plot = NULL,
                              huc12_path       = NULL,
                              title_suffix     = NULL) {
  df <- pts_xy
  if (!var_for_color_scale %in% names(df)) {
    stop(sprintf("Column '%s' not found in pts_xy.", var_for_color_scale))
  }
  
  # Severity filter
  if (!identical(sev_filter, "All")) {
    if (!"sev_group" %in% names(df)) stop("'sev_group' not available.")
    keep <- tolower(df$sev_group) %in% tolower(sev_filter)
    df   <- df[keep, , drop = FALSE]
  }
  if (!nrow(df)) stop("No rows to plot after sev_filter.")
  
  # Optional plotting downsample
  if (!is.null(sample_frac_plot)) {
    stopifnot(sample_frac_plot > 0, sample_frac_plot <= 1)
    set.seed(1L)
    df <- dplyr::sample_frac(df, sample_frac_plot)
  }
  
  # HUC12 backdrop (only resolved here)
  huc12 <- NULL
  if (!is.null(huc12_path)) {
    htmp <- try(sf::st_read(huc12_path, quiet = TRUE), silent = TRUE)
    if (!inherits(htmp, "try-error")) {
      huc12 <- suppressMessages(sf::st_transform(htmp, 32612))
    }
  }
  
  p_base <- ggplot()
  if (!is.null(huc12)) {
    p_base <- p_base + geom_sf(data = huc12, fill = NA, color = "grey60", linewidth = 0.4)
  }
  
  p <- p_base +
    geom_point(
      data = df,
      aes(x = x, y = y, color = .data[[var_for_color_scale]]),
      size = point_size, alpha = 0.9, shape = 16, na.rm = TRUE
    ) +
    scale_color_viridis_c(na.value = "transparent") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank()) +
    labs(
      title = sprintf(
        "YSF: %s • %s%s",
        paste(sort(unique(df$years_since_fire)), collapse = ", "),
        var_for_color_scale,
        if (!is.null(title_suffix)) paste0(" • ", title_suffix) else ""
      ),
      subtitle = if (identical(sev_filter, "All"))
        "All severity groups"
      else
        paste0("sev_group: ", paste(sev_filter, collapse = ", ")),
      x = NULL, y = NULL, color = var_for_color_scale
    )
  
  if (!is.null(huc12)) {
    p <- p + coord_sf(datum = NA)   # required if HUC12 sf layer is present
  } else {
    p <- p + coord_equal(expand = FALSE)
  }
  
  if (length(unique(df$years_since_fire)) > 1) {
    p <- p + facet_wrap(~ years_since_fire, ncol = 2, scales = "fixed")
  }
  
  p
}

# ========================= RUN  =========================

# Build plotting points based on user toggles above
pts_xy <- prep_points_multi_ysf(
  df        = data_long_complete,
  ysf_vec   = target_ysf,
  keep_cols = c("sev_group", var_for_color_scale, resp_var, pred_var),
  resp_var  = resp_var,  resp_op  = resp_op,  resp_thresh  = resp_thresh,
  pred_var  = pred_var,  pred_op  = pred_op,  pred_thresh  = pred_thresh,
  sample_frac_pre = sample_frac_pre
)

# Title suffix documenting active filters (always include var names; add thresholds if set)
mk_filter_label <- function(var, op, thr) {
  if (is.null(var)) return(NA_character_)
  if (is.null(thr)) return(sprintf("%s (no threshold)", var))
  sprintf("%s %s %s", var, op, thr)
}

labels <- c(
  mk_filter_label(resp_var, resp_op, resp_thresh),
  mk_filter_label(pred_var, pred_op, pred_thresh)
)

# keep only real strings; drop NA/empty
labels <- labels[!is.na(labels) & nzchar(labels)]

title_suffix <- if (length(labels)) paste(labels, collapse = " • ") else NULL

# Plot — pass everything explicitly (no global lookups at call time)
print(
  plot_points_multi(
    pts_xy              = pts_xy,
    var_for_color_scale = var_for_color_scale,
    sev_filter          = sev_filter,
    point_size          = point_size,
    sample_frac_plot    = sample_frac_plot,
    huc12_path          = huc12_path,
    title_suffix        = title_suffix
  )
)

# ======================================================================
# === OPTIONAL EXPORT: selected pixels (same filters as plot) to SHP ===
# ======================================================================

# ---- USER OUTPUT SETTINGS ----
out_dir  <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/Shapefiles"
out_name <- "5YSF_Pixels_deltaNDVI_010826"  

# ---- Rebuild the filtered subset (based on 'User Toggles' above) ----
# This mirrors the above plotting pipeline: YSF filter + optional thresholds + optional pre-sampling,
# then GeoJSON -> sf, then project.
df_sub_export <- data_long_complete %>%
  dplyr::filter(years_since_fire %in% target_ysf) %>%
  dplyr::select(dplyr::any_of(c(
    "pixel_ID", "latitude", "longitude", 
    "sev_group", "sev_num", "fire_name", "ref_year", "year",
    "years_since_fire", 
    "delta_ndvi_min", 
    "twi", 
    "cwd_5yr_postfire_avg_08",
    "pptz_JJA",
    "tmaxz_JJA", "tmeanz_JJA",
    "veg_climate_index_08",
    "swez_Apr",
    var_for_color_scale, resp_var, pred_var,
    ".geo"
  )))

# Apply the SAME thresholds (only if thresholds are enabled)
if (!is.null(resp_thresh) && !is.null(resp_var) && resp_var %in% names(df_sub_export)) {
  df_sub_export <- apply_threshold(df_sub_export, resp_var, resp_op, resp_thresh)
}
if (!is.null(pred_thresh) && !is.null(pred_var) && pred_var %in% names(df_sub_export)) {
  df_sub_export <- apply_threshold(df_sub_export, pred_var, pred_op, pred_thresh)
}

if (!nrow(df_sub_export)) stop("No rows after applying YSF + threshold filters.")

# Apply the SAME pre-sampling (only if enabled)
if (!is.null(sample_frac_pre)) {
  stopifnot(sample_frac_pre > 0, sample_frac_pre <= 1)
  set.seed(1L)
  df_sub_export <- dplyr::sample_frac(df_sub_export, sample_frac_pre)
}

# ---- GeoJSON -> sf ----
# Use existing function so CRS guessing behavior is consistent
pts_sf <- as_sf_from_geo_subset(df_sub_export, crs_guess = 4326)

# Project to UTM 12N (matches your mapping)
pts_sf <- suppressMessages(sf::st_transform(pts_sf, 32612))

# Optional: create geometry explicitly as centroids (POINTS instead of pixels)
# pts_sf <- sf::st_centroid(pts_sf)

# ---- Clean shapefile field names + types (SHP limits) ----
names(pts_sf) <- gsub("[^A-Za-z0-9_]", "_", names(pts_sf))
names(pts_sf) <- substr(names(pts_sf), 1, 12)  # SHP field name limit

# Convert POSIXct (if any) to character to avoid driver issues
time_cols <- vapply(pts_sf, inherits, logical(1), what = "POSIXct")
if (any(time_cols)) pts_sf[time_cols] <- lapply(pts_sf[time_cols], as.character)

# ---- Write shapefile ----
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_path <- file.path(out_dir, paste0(out_name, ".shp"))

# delete existing layer if present
if (file.exists(out_path)) {
  sf::st_delete(out_path, quiet = TRUE)
}

## UNCOMMENT THESE LINES TO WRITE SHAPEFILE!! ##
#sf::st_write(pts_sf, out_path, driver = "ESRI Shapefile", append = FALSE)
#message("Wrote shapefile to: ", out_path)

################################################################################


#####   Summary of sample sizes by YSF and sev_group   #########################

# Step 1: Total N per YSF
total_n <- data_long_complete %>%
  filter(!is.na(delta_ndvi_min)) %>%
  group_by(years_since_fire) %>%
  summarise(Total_N = n(), .groups = "drop")

# Step 2: N per severity group per YSF 
sev_counts <- data_long_complete %>%
  filter(!is.na(delta_ndvi_min)) %>%
  st_drop_geometry() %>%    ###### DROP THE GEOMETRY
  group_by(years_since_fire, sev_group) %>%
  summarise(N = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = sev_group,
    values_from = N,
    values_fill = 0
  ) %>%
  rename(
    Unburned_N = `Unburned`,
    Low_N = `Low`,
    Moderate_N = `Moderate`,
    High_N = `High`
  )

# Step 3: Join total and severity counts
ysf_summary <- total_n %>%
  left_join(sev_counts, by = "years_since_fire") %>%
  arrange(years_since_fire)

################################################################################


#############################################################################
###    DELTA NDVI BOXPLOTS: NDVI OVER TIME       ############################

# Box Plot Font controls:
title_size <- 20
subtitle_size <- 15
axis_title_size <- 13
axis_text_size <- 12
caption_size <- 12
annotation_size <- 4
base_font <- "Times New Roman"

################################################################################
####### ============  HIGH SEVERITY BOXPLOT ====================================

high_sev <- data_long_complete %>%
  filter(sev_group == "High", years_since_fire <= 30)

label_counts <- high_sev %>%
  group_by(years_since_fire) %>%
  summarise(n = n(), .groups = "drop")

plot_high <- ggplot(high_sev, aes(x = factor(years_since_fire), y = delta_ndvi_min)) +
  geom_boxplot(fill = "firebrick", alpha = 0.7, outlier.size = 0.5) +
  geom_text(data = label_counts, aes(label = n, y = 0.98),
            color = "black", size = annotation_size, family = base_font) +
  coord_cartesian(ylim = c(-1, 1)) +
  theme_minimal(base_family = base_font, base_size = 12) +
  labs(
    title = "ΔNDVI Over Time Since Fire (High Severity)",
    subtitle = "40-meter minimum pixel spacing",
    x = "Years Since Fire",
    y = "ΔNDVI (Postfire - Prefire)",
    caption = "Numbers above boxes indicate sample size."
  ) +
  theme(
    plot.title    = element_text(size = title_size, face = "bold", family = base_font),
    plot.subtitle = element_text(size = subtitle_size, family = base_font),
    axis.title.x  = element_text(size = axis_title_size, face = "bold", family = base_font, margin = margin(t = 12)),
    axis.title.y  = element_text(size = axis_title_size, face = "bold", family = base_font),
    axis.text.x   = element_text(size = axis_text_size, family = base_font),
    axis.text.y   = element_text(size = axis_text_size, family = base_font),
    plot.caption  = element_text(size = caption_size, family = base_font)
  )

print(plot_high)

################################################################################
####### ====== ====  MODERATE SEVERITY BOXPLOT ================================

moderate_sev <- data_long_complete %>%
  filter(sev_group == "Moderate", years_since_fire <= 30)

label_counts_mod <- moderate_sev %>%
  group_by(years_since_fire) %>%
  summarise(n = n(), .groups = "drop")

plot_mod <- ggplot(moderate_sev, aes(x = factor(years_since_fire), y = delta_ndvi_min)) +
  geom_boxplot(fill = "orange", alpha = 0.7, outlier.size = 0.5) +
  geom_text(data = label_counts_mod, aes(label = n, y = 0.98),
            color = "black", size = annotation_size, family = base_font) +
  coord_cartesian(ylim = c(-1, 1)) +
  theme_minimal(base_family = base_font, base_size = 12) +
  labs(
    title = "ΔNDVI Over Time Since Fire (Moderate Severity)",
    subtitle = "40-meter minimum pixel spacing",
    x = "Years Since Fire",
    y = "ΔNDVI (Postfire - Prefire)",
    caption = "Numbers above boxes indicate sample size."
  ) +
  theme(
    plot.title    = element_text(size = title_size, face = "bold", family = base_font),
    plot.subtitle = element_text(size = subtitle_size, family = base_font),
    axis.title.x  = element_text(size = axis_title_size, face = "bold", family = base_font, margin = margin(t = 12)),
    axis.title.y  = element_text(size = axis_title_size, face = "bold", family = base_font),
    axis.text.x   = element_text(size = axis_text_size, family = base_font),
    axis.text.y   = element_text(size = axis_text_size, family = base_font),
    plot.caption  = element_text(size = caption_size, family = base_font)
  )

print(plot_mod)

################################################################################
####### =========  LOW SEVERITY BOXPLOT ========================================

low_sev <- data_long_complete %>%
  filter(sev_group == "Low", years_since_fire <= 30)

label_counts_low <- low_sev %>%
  group_by(years_since_fire) %>%
  summarise(n = n(), .groups = "drop")

plot_low <- ggplot(low_sev, aes(x = factor(years_since_fire), y = delta_ndvi_min)) +
  geom_boxplot(fill = "yellow", alpha = 0.7, outlier.size = 0.5) +
  geom_text(data = label_counts_low, aes(label = n, y = 0.98),
            color = "black", size = annotation_size, family = base_font) +
  coord_cartesian(ylim = c(-1, 1)) +
  theme_minimal(base_family = base_font, base_size = 12) +
  labs(
    title = "ΔNDVI Over Time Since Fire (Low Severity)",
    subtitle = "40-meter minimum pixel spacing",
    x = "Years Since Fire",
    y = "ΔNDVI (Postfire - Prefire)",
    caption = "Numbers above boxes indicate sample size."
  ) +
  theme(
    plot.title    = element_text(size = title_size, face = "bold", family = base_font),
    plot.subtitle = element_text(size = subtitle_size, family = base_font),
    axis.title.x  = element_text(size = axis_title_size, face = "bold", family = base_font, margin = margin(t = 12)),
    axis.title.y  = element_text(size = axis_title_size, face = "bold", family = base_font),
    axis.text.x   = element_text(size = axis_text_size, family = base_font),
    axis.text.y   = element_text(size = axis_text_size, family = base_font),
    plot.caption  = element_text(size = caption_size, family = base_font)
  )

print(plot_low)

################################################################################
##########    UNBURNED BOXPLOT ================================================

unburned <- data_long_complete %>%
  filter(sev_class %in% c("None", "Unburned/Very Low"),
         years_since_fire >= 1,
         years_since_fire <= 20)

label_counts_unburned <- unburned %>%
  group_by(years_since_fire) %>%
  summarise(n = n(), .groups = "drop")

plot_unburned <- ggplot(unburned, aes(x = factor(years_since_fire), y = delta_ndvi_min)) +
  geom_boxplot(fill = "gray60", alpha = 0.7, outlier.size = 0.5) +
  #geom_text(data = label_counts_unburned, aes(label = n, y = 0.98),
  #          color = "black", size = annotation_size, family = base_font) +
  coord_cartesian(ylim = c(-1, 1)) +
  theme_minimal(base_family = base_font, base_size = 12) +
  labs(
    title = "ΔNDVI Over Time (Unburned Pixels)",
    #subtitle = "40-meter minimum pixel spacing",
    x = "Years Since Matched Fire Event",
    y = "ΔNDVI (Postfire - Prefire)" #,
    #caption = "Numbers above boxes indicate sample size."
  ) +
  theme(
    plot.title    = element_text(size = title_size, face = "bold", family = base_font),
    plot.subtitle = element_text(size = subtitle_size, family = base_font),
    axis.title.x  = element_text(size = axis_title_size, face = "bold", family = base_font, margin = margin(t = 12)),
    axis.title.y  = element_text(size = axis_title_size, face = "bold", family = base_font),
    axis.text.x   = element_text(size = axis_text_size, family = base_font),
    axis.text.y   = element_text(size = axis_text_size, family = base_font),
    plot.caption  = element_text(size = caption_size, family = base_font)
  )

print(plot_unburned)
################################################################################
##### Combining all 4 into one panel: ===========================
library(patchwork)

(plot_high | plot_mod) / (plot_low | plot_unburned)
################################################################################
########  Plot ΔNDVI 30 years post-fire (all sev_groups on one figure)  =========
data_long_30 <- data_long_complete %>%
  mutate(
    years_since_fire = as.integer(year - ref_year),
    sev_group = factor(sev_group, levels = c("Unburned", "Low", "Moderate", "High"))
  ) %>%
  filter(years_since_fire %in% 1:30)

# Get sample sizes
sample_sizes_30 <- data_long_30 %>%
  group_by(years_since_fire, sev_group) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = paste0("n = ", n))

# Plotting 30 years post-fire:
# Determine dynamic y-position just above the max observed value
label_y_position <- max(data_long_30$delta_ndvi_min, na.rm = TRUE) + 0.05
# Create plot
ggplot(data_long_30, aes(x = factor(years_since_fire), y = delta_ndvi_min, fill = sev_group)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8, position = position_dodge(width = 0.8)) +
  geom_text(data = sample_sizes_30,
            aes(x = factor(years_since_fire), y = label_y_position, label = label, group = sev_group),
            position = position_dodge(width = 0.8),
            angle = 90,           # <<< rotate vertically
            hjust = -0.1,         # <<< adjust horizontal justification
            vjust = 0.5,
            size = 2.5,
            inherit.aes = FALSE) +
  scale_fill_manual(
    values = c("Unburned" = "gray50", "Low" = "yellow", "Moderate" = "orange", "High" = "firebrick"),
    name = "Severity Class"
  ) +
  theme_minimal() +
  labs(
    title = "ΔNDVI for 30 Years Post-Fire, by Severity Class",
    subtitle = "Dataset: data_long_complete, Nov '25",
    x = "Years Since Fire",
    y = "ΔNDVI",
    caption = "Sample sizes above each box represent pixel-year observations"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
################################################################################


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
################################################################################
#             -------------------------------------------                      #
########          RESEARCH QUESTION 1: INITIAL FIRE EFFECTS     ################
#             -------------------------------------------                      #
################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


################################################################################
# DeltaNDVI Box Plots, ANOVA & Tukey for 1, 2, 3 YSF - full dataset & subsets
## Note: Delta NDVI calculated based on pre-fire Minimum NDVI

# Libraries ---------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(purrr)
library(rstatix)  # tidy ANOVA & posthoc
library(tibble)
library(scales)   # for comma() in captions
library(multcompView)   # for compact letter display on plots
library(ggplot2)

# ------------------------------------------------------------------------------
set.seed(42)
pct = 0.2  # USE THIS TO SET PERCENTAGE FOR ANY RANDOMIZED SUBSETS (BELOW)
# ------------------------------------------------------------------------------
# HELPER FUNCTIONS (APPLY TO ANY DATASET WITH years_since_fire, sev_group, delta_ndvi_min)
# ------------------------------------------------------------------------------

# 1) Safe ANOVA table for a single data chunk (one YSF)
safe_anova_tbl <- function(df) {
  df2 <- df %>%
    dplyr::filter(!is.na(delta_ndvi_min)) %>%
    droplevels()
  
  if (dplyr::n_distinct(df2$sev_group) < 2L) {
    return(tibble(
      Effect = NA_character_,
      DFn    = NA_real_,
      DFd    = NA_real_,
      F      = NA_real_,
      p      = NA_real_,
      ges    = NA_real_
    ))
  }
  
  out <- tryCatch(
    rstatix::anova_test(data = df2, delta_ndvi_min ~ sev_group),
    error = function(e) NULL
  )
  
  if (is.null(out)) {
    return(tibble(
      Effect = NA_character_,
      DFn    = NA_real_,
      DFd    = NA_real_,
      F      = NA_real_,
      p      = NA_real_,
      ges    = NA_real_
    ))
  }
  
  out %>%
    rstatix::get_anova_table() %>%
    as_tibble()
}

# 2) Safe Tukey table for a single data chunk (one YSF)
safe_tukey_tbl <- function(df) {
  df2 <- df %>%
    dplyr::filter(!is.na(delta_ndvi_min)) %>%
    droplevels()
  
  # Need ≥ 2 groups and at least 2 obs in each group
  if (dplyr::n_distinct(df2$sev_group) < 2L || min(table(df2$sev_group)) < 2L) {
    return(tibble(
      term         = NA_character_,
      group1       = NA_character_,
      group2       = NA_character_,
      estimate     = NA_real_,
      conf.low     = NA_real_,
      conf.high    = NA_real_,
      p.adj        = NA_real_,
      p.adj.signif = NA_character_
    ))
  }
  
  fit <- tryCatch(
    aov(delta_ndvi_min ~ sev_group, data = df2),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(tibble(
      term         = NA_character_,
      group1       = NA_character_,
      group2       = NA_character_,
      estimate     = NA_real_,
      conf.low     = NA_real_,
      conf.high    = NA_real_,
      p.adj        = NA_real_,
      p.adj.signif = NA_character_
    ))
  }
  
  out <- tryCatch(
    rstatix::tukey_hsd(fit),
    error = function(e) NULL
  )
  if (is.null(out)) {
    return(tibble(
      term         = NA_character_,
      group1       = NA_character_,
      group2       = NA_character_,
      estimate     = NA_real_,
      conf.low     = NA_real_,
      conf.high    = NA_real_,
      p.adj        = NA_real_,
      p.adj.signif = NA_character_
    ))
  }
  
  as_tibble(out)
}

# 3) Per-YSF Ns and compact string, for any dataset
make_n_by_year <- function(df) {
  df %>%
    dplyr::count(years_since_fire, sev_group, name = "n_sampled") %>%
    tidyr::pivot_wider(names_from = sev_group, values_from = n_sampled) %>%
    dplyr::mutate(
      N_total = rowSums(dplyr::across(c(Unburned, Low, Moderate, High)), na.rm = TRUE),
      N_by_group = paste0(
        "Unburned=", dplyr::coalesce(Unburned, 0), "; ",
        "Low=",      dplyr::coalesce(Low, 0),      "; ",
        "Moderate=", dplyr::coalesce(Moderate, 0), "; ",
        "High=",     dplyr::coalesce(High, 0)
      )
    )
}

# 4) Generic ANOVA summary for any dataset (by YSF) - NO sample size columns
build_anova_summary <- function(df,
                                n_by_year_tbl = NULL, # kept for compatibility; not used
                                alpha = 0.01,
                                N_total_name = "N_total",
                                N_by_group_name = "N_by_group") {
  df %>%
    dplyr::group_by(years_since_fire) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      anova_tbl = purrr::map(data, ~ safe_anova_tbl(.x))
    ) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(anova_tbl) %>%
    dplyr::filter(Effect == "sev_group" | is.na(Effect)) %>%
    dplyr::mutate(
      `ANOVA F(3, DFd)` = dplyr::if_else(
        is.na(DFd),
        NA_character_,
        paste0("F(3, ", DFd, ") = ", round(F, 1))
      ),
      `p-value` = dplyr::case_when(
        is.na(p)  ~ NA_character_,
        p < 0.001 ~ "< 0.001",
        TRUE      ~ format(p, digits = 3, scientific = TRUE)
      ),
      `Sig (p < alpha)` = dplyr::case_when(
        is.na(p)  ~ NA_character_,
        p < alpha ~ "Yes",
        TRUE      ~ "No"
      ),
      `Generalized eta^2 (ges)` = round(ges, 3)
    ) %>%
    dplyr::select(
      years_since_fire,
      `ANOVA F(3, DFd)`,
      `p-value`,
      `Sig (p < alpha)`,
      `Generalized eta^2 (ges)`
    ) %>%
    dplyr::arrange(years_since_fire)
}


# 5) Generic Tukey summary for any dataset (by YSF)
build_tukey_summary <- function(df, n_by_year_tbl) {
  df %>%
    dplyr::group_by(years_since_fire) %>%
    tidyr::nest() %>%
    dplyr::mutate(tukey_tbl = purrr::map(data, safe_tukey_tbl)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(tukey_tbl) %>%
    dplyr::left_join(n_by_year_tbl, by = "years_since_fire") %>%
    dplyr::arrange(years_since_fire, group1, group2)
}

# 6) Stratified percentage subset (by YSF × sev_group) using "pct"
make_pct_subset <- function(df,
                            pct            = pct,
                            min_n_per_cell = 2,
                            subset_label   = NULL) 
{
  if (is.null(subset_label)) {
    subset_label <- paste0(
      scales::percent(pct, accuracy = 1),
      " per YSF × severity"
    )
  }
}

# 7) Compact Letter Display (CLD) table per YSF (safe)
build_cld_tbl <- function(df,
                          ysf_col       = "years_since_fire",
                          group_col     = "sev_group",
                          response_col  = "delta_ndvi_min",
                          group_levels  = c("Unburned", "Low", "Moderate", "High"),
                          alpha         = 0.01) {
  
  stopifnot(all(c(ysf_col, group_col, response_col) %in% names(df)))
  
  df2 <- df %>%
    dplyr::filter(!is.na(.data[[response_col]])) %>%
    dplyr::mutate(
      !!group_col := factor(.data[[group_col]], levels = group_levels)
    ) %>%
    droplevels()
  
  df2 %>%
    dplyr::group_by(.data[[ysf_col]]) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      # Tukey per YSF (safe)
      tukey_tbl = purrr::map(data, ~ safe_tukey_tbl(.x)),
      # Letters per YSF (use alpha threshold)
      letters = purrr::map(tukey_tbl, ~ {
        tt <- .x %>% dplyr::filter(!is.na(p.adj))
        
        # If Tukey couldn't be computed / no comparisons, fall back to all "a"
        if (nrow(tt) == 0) {
          levs <- levels(df2[[group_col]])
          out  <- rep("a", length(levs))
          names(out) <- levs
          return(out)
        }
        
        pvec <- tt$p.adj
        names(pvec) <- paste(tt$group1, tt$group2, sep = "-")
        
        multcompView::multcompLetters(pvec, threshold = alpha)$Letters
      })
    ) %>%
    dplyr::select(all_of(ysf_col), letters) %>%
    tidyr::unnest_wider(letters, names_sep = "_") %>%
    tidyr::pivot_longer(
      cols      = starts_with("letters_"),
      names_to  = group_col,
      values_to = "cld"
    ) %>%
    dplyr::mutate(
      !!group_col := gsub("^letters_", "", .data[[group_col]]),
      !!group_col := factor(.data[[group_col]], levels = group_levels)
    )
}

# 8) Alphabetize CLD letters so they start at "a" for the first group in group_levels
canonicalize_cld_letters <- function(L, group_levels) {
  # L is a named character vector from multcompLetters()$Letters
  # names(L) should be group names
  
  # Ensure all groups exist in L (fill missing with "a" as fallback)
  L_full <- setNames(rep("a", length(group_levels)), group_levels)
  L_full[names(L)] <- as.character(L)
  
  # Collect the distinct letter symbols in the order they appear
  # when scanning groups left-to-right in group_levels.
  symbols <- character(0)
  for (g in group_levels) {
    chars <- strsplit(L_full[[g]], "", fixed = TRUE)[[1]]
    for (ch in chars) {
      if (!ch %in% symbols) symbols <- c(symbols, ch)
    }
  }
  
  # Map old symbols -> new symbols: a, b, c, ...
  new_symbols <- letters[seq_along(symbols)]
  map <- setNames(new_symbols, symbols)
  
  # Apply mapping to each group's string without accidental cascading
  L_new <- vapply(L_full, function(s) {
    chars <- strsplit(s, "", fixed = TRUE)[[1]]
    paste0(map[chars], collapse = "")
  }, FUN.VALUE = character(1))
  
  # Return as named vector in the same group_levels order
  setNames(L_new[group_levels], group_levels)
}

# ------------------------------------------------------------------------------
# BUILD SNAPSHOT (YSF = X) FROM FULL DATASET
# ------------------------------------------------------------------------------
snapshot <- data_long_complete %>%
  dplyr::mutate(
    years_since_fire = year - ref_year,
    sev_group = factor(sev_group, levels = c("Unburned", "Low", "Moderate", "High"))
  ) %>%
  dplyr::filter(years_since_fire %in% c(5, 10, 15, 20)) %>%
  dplyr::filter(!is.na(delta_ndvi_min)) %>%
  droplevels()

# ------------------------------------------------------------------------------
# FULL DATASET: ANOVA, Tukey, N TABLE, BOXPLOT
# ------------------------------------------------------------------------------
# N table for full dataset
full_n_by_year <- make_n_by_year(snapshot)

# ANOVA summary for full dataset
anova_results_full <- build_anova_summary(
  df            = snapshot,
  n_by_year_tbl = full_n_by_year
)

# Tukey summary for full dataset
tukey_results_full <- build_tukey_summary(
  df            = snapshot,
  n_by_year_tbl = full_n_by_year
)

# Mean + 95% CI per YSF x severity (FULL)
mean_ci_full <- snapshot %>%
  dplyr::group_by(years_since_fire, sev_group) %>%
  dplyr::summarise(
    n     = dplyr::n(),
    mean  = mean(delta_ndvi_min, na.rm = TRUE),
    sd    = sd(delta_ndvi_min, na.rm = TRUE),
    se    = sd / sqrt(n),
    tcrit = qt(0.975, df = n - 1),
    ci_low  = mean - tcrit * se,
    ci_high = mean + tcrit * se,
    .groups = "drop"
  )

# Mean deltaNDVI summary table — FULL DATASET
mean_deltaNDVI_by_sev_time_full <- snapshot %>%
  dplyr::group_by(years_since_fire, sev_group) %>%
  dplyr::summarise(
    mean = mean(delta_ndvi_min, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(mean = round(mean, 2)) %>%
  tidyr::pivot_wider(
    names_from  = sev_group,
    values_from = mean
  ) %>%
  dplyr::arrange(years_since_fire)

print(mean_deltaNDVI_by_sev_time_full)


# Compact Letter Display (CLD) for FULL
cld_full <- build_cld_tbl(snapshot)

# Y positions for CLD letters (FULL) — match 250 formatting
cld_pos_full <- snapshot %>%
  dplyr::group_by(years_since_fire, sev_group) %>%
  dplyr::summarise(
    y_pos = max(delta_ndvi_min, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(cld_full, by = c("years_since_fire", "sev_group")) %>%
  dplyr::mutate(y_pos = y_pos + 0.03)

# Full dataset boxplot (WITH MEAN+CI + CLD)
sample_sizes_full <- snapshot %>%
  dplyr::group_by(years_since_fire, sev_group) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(label = paste0("n = ", n))

N_min_full <- min(sample_sizes_full$n, na.rm = TRUE)
N_max_full <- max(sample_sizes_full$n, na.rm = TRUE)

#N_caption_full <- paste0(
#  "Full dataset. Boxplots grouped by severity class for each time point. N ranges from ",
#  scales::comma(N_min_full), "–", scales::comma(N_max_full), "."
#)

show_N_full <- FALSE

# ---- Dynamic title from whatever YSFs are actually in snapshot
ysf_vals_full <- snapshot %>%
  dplyr::distinct(years_since_fire) %>%
  dplyr::arrange(years_since_fire) %>%
  dplyr::pull(years_since_fire)

title_text_full <- paste0(
  "ΔNDVI at ",
  paste(ysf_vals_full, collapse = ", "),
  " Years Post-Fire"
)

# Full dataset boxplot (MATCH 250 formatting: mean-centered "boxplot" + CI + CLD)
p_full <- ggplot(
  snapshot,
  aes(x = factor(years_since_fire), y = delta_ndvi_min, fill = sev_group)
) +
  stat_summary(
    fun.data = function(y) {
      y <- y[is.finite(y)]
      if (length(y) == 0) {
        return(data.frame(ymin = NA, lower = NA, middle = NA, upper = NA, ymax = NA))
      }
      
      q1  <- unname(quantile(y, 0.25, na.rm = TRUE))
      q3  <- unname(quantile(y, 0.75, na.rm = TRUE))
      iqr <- q3 - q1
      
      low_fence  <- q1 - 1.5 * iqr
      high_fence <- q3 + 1.5 * iqr
      
      whisk_low  <- min(y[y >= low_fence],  na.rm = TRUE)
      whisk_high <- max(y[y <= high_fence], na.rm = TRUE)
      
      data.frame(
        ymin   = whisk_low,
        lower  = q1,
        middle = mean(y, na.rm = TRUE),  # center line = MEAN (matches your 250 section)
        upper  = q3,
        ymax   = whisk_high
      )
    },
    geom      = "boxplot",
    alpha     = 0.85,
    linewidth = 0.6,
    width     = 0.75,
    position  = position_dodge(width = 0.85)
  ) +
  # 95% CI error bars for the mean (same as 250 section)
  geom_errorbar(
    data = mean_ci_full,
    aes(
      x = factor(years_since_fire),
      y = mean,
      ymin = ci_low,
      ymax = ci_high,
      group = sev_group
    ),
    position = position_dodge(width = 0.85),
    width = 0.18,
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  # CLD letters (Tukey HSD within each YSF)
  geom_text(
    data = cld_pos_full,
    aes(
      x = factor(years_since_fire),
      y = y_pos,
      label = cld,
      group = sev_group
    ),
    position = position_dodge(width = 0.85),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold",
    family = "Times New Roman"
  ) +
  scale_fill_manual(
    values = c(
      "Unburned" = "gray50",
      "Low"      = "yellow",
      "Moderate" = "orange",
      "High"     = "firebrick"
    ),
    name = "Severity Class"
  ) +
  theme_minimal(base_family = "Times New Roman", base_size = 12) +
  labs(
    title    = title_text_full,
    subtitle = "Full Dataset",
    x        = "Years Since Fire",
    y        = "ΔNDVI"
    #caption = N_caption_full
  ) +
  theme(
    plot.title    = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    axis.title.x  = element_text(size = 13, face = "bold"),
    axis.title.y  = element_text(size = 13, face = "bold"),
    axis.text.x   = element_text(size = 11, face = "bold"),
    axis.text.y   = element_text(size = 11, face = "bold"),
    legend.title  = element_text(size = 13, face = "bold"),
    legend.text   = element_text(size = 12, face = "bold")
    #plot.caption  = element_text(size = 12)
  )

print(p_full)
# ------------------------------------------------------------------------------
# 250-PER-SEV_GROUP SUBSET: ANOVA, Tukey, N TABLE
# (CONSENSUS / MONTE CARLO SEED SELECTION; STRATIFIED BY FIRE EVENT)
# ------------------------------------------------------------------------------

# BOOTSTRAP / MONTE CARLO CONSENSUS SAMPLING (500 runs)
# Goal: Use HUC-balanced sampling + Tukey/CLD, then choose the MOST COMMON
# CLD (compact-letter-display) pattern across "B" runs; re-run using "best" seed

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(rstatix)
  library(multcompView)
  library(future)
  library(furrr)
})

# ---- USER SETTINGS 
B <- 500
target_n_per_cell <- 250
alpha <- 0.01   # <-- set desired significance threshold here
ysf_levels <- sort(unique(snapshot$years_since_fire))
sev_levels <- c("Unburned","Low","Moderate","High")

# FUNCTIONS:
# 1) Fire event - balanced sub-sampling for each [YSFxSeverity] group.
sample_event_balanced_one_group <- function(df_group, n_target = 250, event_col = "fire_name") {
  stopifnot(event_col %in% names(df_group))
  
  df_group <- df_group %>% dplyr::filter(!is.na(.data[[event_col]]))
  if (!nrow(df_group)) return(df_group[0, , drop = FALSE])
  
  df_group <- df_group %>%
    dplyr::arrange(.data[[event_col]], pixel_ID)
  
  # If not enough rows, take all
  if (nrow(df_group) <= n_target) return(df_group)
  
  # Distribute target across available events (approximately equal allocation)
  events <- sort(unique(df_group[[event_col]]))
  E      <- length(events)
  if (E == 0L) return(dplyr::slice_sample(df_group, n = n_target))
  
  base_n <- floor(n_target / E)
  rem    <- n_target - base_n * E
  
  # Randomly choose which events get +1 to account for remainder
  extra_events <- if (rem > 0) sample(events, size = rem, replace = FALSE) else character(0)
  
  # Sample within each event allocation (capped by availability)
  out <- df_group %>%
    dplyr::group_by(.data[[event_col]]) %>%
    dplyr::group_modify(~{
      ev_val <- .y[[1]]
      take_n <- base_n + as.integer(ev_val %in% extra_events)
      take_n <- min(take_n, nrow(.x))
      if (take_n <= 0) return(.x[0, , drop = FALSE])
      dplyr::slice_sample(.x, n = take_n)
    }) %>%
    dplyr::ungroup()
  
  # Top off to reach n_target using remaining unsampled rows
  if (nrow(out) < n_target) {
    need <- n_target - nrow(out)
    
    out_ids <- NULL
    if ("pixel_ID" %in% names(out)) out_ids <- out$pixel_ID
    
    remaining <- df_group
    if (!is.null(out_ids) && "pixel_ID" %in% names(df_group)) {
      remaining <- df_group %>% dplyr::filter(!pixel_ID %in% out_ids)
    }
    
    if (nrow(remaining) > 0) {
      out <- dplyr::bind_rows(
        out,
        dplyr::slice_sample(remaining, n = min(need, nrow(remaining)))
      )
    }
  }
  
  # Safety: return exactly n_target if possible
  if (nrow(out) > n_target) {
    out <- dplyr::slice_sample(out, n = n_target)
  }
  
  out
}

# Pre-split snapshot ONCE (trim columns) - run in parallel
keep_cols <- c("years_since_fire", "sev_group", "delta_ndvi_min", "fire_name", "pixel_ID")
keep_cols <- intersect(keep_cols, names(snapshot))

snapshot_small <- snapshot %>%
  dplyr::select(dplyr::all_of(keep_cols)) %>%
  dplyr::mutate(sev_group = factor(sev_group, levels = sev_levels)) %>%
  droplevels()

key <- interaction(snapshot_small$years_since_fire, snapshot_small$sev_group, drop = TRUE)
groups <- split(snapshot_small, key)

# make group iteration order deterministic
groups <- groups[order(names(groups))]

rm(snapshot_small)
gc()

# 3) CONSENSUS SEARCH (PER-YSF): choose modal CLD independently for each YSF
cld_signature_one_ysf <- function(L, sev_levels = c("Unburned","Low","Moderate","High")) {
  L_full <- setNames(rep("a", length(sev_levels)), sev_levels)
  L_full[names(L)] <- as.character(L)
  paste0(
    paste0(substr(sev_levels, 1, 1), "=", unname(L_full[sev_levels])),
    collapse = "|"
  )
}

run_once_signature_one_ysf <- function(seed,
                                       groups_one_ysf,
                                       target_n_per_cell = 250,
                                       sev_levels = c("Unburned","Low","Moderate","High"),
                                       alpha = 0.01) {
  set.seed(seed)
  
  sampled_list <- lapply(groups_one_ysf, function(g) {
    sample_event_balanced_one_group(g, n_target = target_n_per_cell, event_col = "fire_name")
  })
  samp <- dplyr::bind_rows(sampled_list)
  samp$sev_group <- factor(samp$sev_group, levels = sev_levels)
  samp <- droplevels(samp)
  
  # Safe fallback
  if (nrow(samp) == 0L ||
      length(unique(samp$sev_group)) < 2L ||
      min(table(samp$sev_group)) < 2L) {
    L <- setNames(rep("a", length(sev_levels)), sev_levels)
    return(tibble(seed = seed, signature = cld_signature_one_ysf(L, sev_levels)))
  }
  
  fit <- tryCatch(aov(delta_ndvi_min ~ sev_group, data = samp), error = function(e) NULL)
  if (is.null(fit)) {
    L <- setNames(rep("a", length(sev_levels)), sev_levels)
    return(tibble(seed = seed, signature = cld_signature_one_ysf(L, sev_levels)))
  }
  
  tk <- tryCatch(rstatix::tukey_hsd(fit), error = function(e) NULL)
  if (is.null(tk) || !nrow(tk)) {
    L <- setNames(rep("a", length(sev_levels)), sev_levels)
    return(tibble(seed = seed, signature = cld_signature_one_ysf(L, sev_levels)))
  }
  
  tk <- as.data.frame(tk)
  tk <- tk[!is.na(tk$p.adj), , drop = FALSE]
  if (!nrow(tk)) {
    L <- setNames(rep("a", length(sev_levels)), sev_levels)
    return(tibble(seed = seed, signature = cld_signature_one_ysf(L, sev_levels)))
  }
  
  pvec <- tk$p.adj
  names(pvec) <- paste(tk$group1, tk$group2, sep = "-")
  
  L_raw <- multcompView::multcompLetters(pvec, threshold = alpha)$Letters
  L <- canonicalize_cld_letters(L_raw, group_levels = sev_levels)
  
  # sort within-string letters ("ba" -> "ab")
  L <- vapply(
    L,
    function(s) paste0(sort(strsplit(s, "", fixed = TRUE)[[1]]), collapse = ""),
    FUN.VALUE = character(1)
  )
  
  tibble(seed = seed, signature = cld_signature_one_ysf(L, sev_levels))
}

# ---- Make sure snapshot_small exists here (you can reuse yours; keeping explicit is safest)
snapshot_small <- snapshot %>%
  dplyr::select(dplyr::all_of(keep_cols)) %>%
  dplyr::mutate(sev_group = factor(sev_group, levels = sev_levels)) %>%
  droplevels()

# Shared candidate seeds across YSFs (optional but consistent)
set.seed(1)
seeds_to_try <- sample.int(1e7, B)

# Parallel plan
workers_use <- max(1, parallel::detectCores() - 1)
options(future.globals.maxSize = 2 * 1024^3)
future::plan(future::multisession, workers = workers_use)

consensus_by_ysf <- lapply(ysf_levels, function(ysf) {
  
  snap_y <- snapshot_small %>% dplyr::filter(years_since_fire == ysf)
  
  # Split into severity groups (deterministic order)
  groups_y <- split(snap_y, snap_y$sev_group, drop = TRUE)
  groups_y <- groups_y[sev_levels[sev_levels %in% names(groups_y)]]
  
  sig_tbl_y <- furrr::future_map_dfr(
    seeds_to_try,
    ~ run_once_signature_one_ysf(
      seed = .x,
      groups_one_ysf = groups_y,
      target_n_per_cell = target_n_per_cell,
      sev_levels = sev_levels,
      alpha = alpha
    ),
    .options  = furrr::furrr_options(seed = TRUE),
    .progress = TRUE
  )
  
  mode_tbl_y <- sig_tbl_y %>%
    dplyr::count(signature, name = "n_runs") %>%
    dplyr::arrange(desc(n_runs)) %>%
    dplyr::mutate(prop_runs = n_runs / sum(n_runs))
  
  best_sig_y <- mode_tbl_y$signature[1]
  
  candidate_seeds_y <- sig_tbl_y %>%
    dplyr::filter(signature == best_sig_y) %>%
    dplyr::pull(seed)
  
  best_seed_y <- candidate_seeds_y[
    which(vapply(candidate_seeds_y, function(s) {
      run_once_signature_one_ysf(
        seed = s,
        groups_one_ysf = groups_y,
        target_n_per_cell = target_n_per_cell,
        sev_levels = sev_levels,
        alpha = alpha
      )$signature == best_sig_y
    }, logical(1)))[1]
  ]
  
  if (is.na(best_seed_y) || length(best_seed_y) == 0) {
    stop(paste0("No candidate seed reproduced modal CLD for YSF ", ysf, "."), call. = FALSE)
  }
  
  list(
    ysf = ysf,
    groups_y = groups_y,
    best_signature = best_sig_y,
    best_seed = best_seed_y,
    mode_tbl_y = mode_tbl_y
  )
})

future::plan(future::sequential)

best_seeds_tbl <- dplyr::bind_rows(lapply(consensus_by_ysf, function(x) {
  tibble(
    years_since_fire = x$ysf,
    best_seed = x$best_seed,
    best_signature = x$best_signature
  )
}))
print(best_seeds_tbl)

# FINAL (LOCKED) SAMPLE: build sampled_250 by YSF using each YSF's best_seed
sampled_250 <- dplyr::bind_rows(lapply(consensus_by_ysf, function(x) {
  set.seed(x$best_seed)
  
  sampled_list <- lapply(x$groups_y, function(g) {
    sample_event_balanced_one_group(g, n_target = target_n_per_cell, event_col = "fire_name")
  })
  
  dplyr::bind_rows(sampled_list)
})) %>%
  dplyr::mutate(sev_group = factor(sev_group, levels = sev_levels)) %>%
  droplevels()

# Quick N check: Ns should be <= 250 where groups are small, otherwise = 250
#sampled_250 %>%
#  dplyr::count(years_since_fire, sev_group, name = "n") %>%
#  print(n = Inf)

# ANALYSES (ANOVA / TUKEY / MEAN+CI / CLD / CLD POSITIONS)
n_by_year_250 <- make_n_by_year(sampled_250)

# ANOVA results table
anova_results_250 <- build_anova_summary(
  df            = sampled_250,
  alpha         = alpha
)

# Tukey's HSD results table
tukey_results_250 <- build_tukey_summary(
  df            = sampled_250,
  n_by_year_tbl = n_by_year_250
) %>%
  dplyr::mutate(
    sig_at_alpha = dplyr::if_else(is.na(p.adj), NA, p.adj < alpha),
    p_adj_display = dplyr::case_when(
      is.na(p.adj)  ~ NA_character_,
      p.adj < alpha ~ paste0("< ", format(alpha, scientific = FALSE, trim = TRUE), "*"),
      TRUE          ~ paste0("> ", format(alpha, scientific = FALSE, trim = TRUE))
    )
  ) %>%
  dplyr::select(
    years_since_fire,
    group1,
    group2,
    mean_diff_deltaNDVI = estimate,
    conf.low,
    conf.high,
    p.adj,
    p_adj_display,
    sig_at_alpha
  ) %>%
  dplyr::mutate(
    dplyr::across(c(mean_diff_deltaNDVI, conf.low, conf.high), ~ round(.x, 2))
  ) %>%
  dplyr::arrange(years_since_fire, group1, group2)

# Mean deltaNDVI summary table
mean_deltaNDVI_by_sev_time_wide <- sampled_250 %>%
  dplyr::group_by(years_since_fire, sev_group) %>%
  dplyr::summarise(
    mean = mean(delta_ndvi_min, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(mean = round(mean, 2)) %>%
  tidyr::pivot_wider(
    names_from  = sev_group,
    values_from = mean
  ) %>%
  dplyr::arrange(years_since_fire)


# FUNCTION: Build compact letter display (CLD) from Tukey table results
build_cld_from_tukey <- function(tukey_tbl,
                                 ysf_col      = "years_since_fire",
                                 group1_col   = "group1",
                                 group2_col   = "group2",
                                 p_col        = "p.adj",
                                 group_levels = c("Unburned", "Low", "Moderate", "High"),
                                 alpha        = alpha) {  #alpha signif. threshold defined above
  
  stopifnot(all(c(ysf_col, group1_col, group2_col, p_col) %in% names(tukey_tbl)))
  
  tukey_tbl %>%
    dplyr::filter(!is.na(.data[[p_col]])) %>%
    dplyr::group_by(.data[[ysf_col]]) %>%
    dplyr::group_modify(~{
      tt <- .x
      
      # Build a full p-value vector for multcompLetters
      pvec <- tt[[p_col]]
      names(pvec) <- paste(tt[[group1_col]], tt[[group2_col]], sep = "-")
      
      # If somehow empty, return all "a"
      if (length(pvec) == 0) {
        return(tibble::tibble(
          sev_group = factor(group_levels, levels = group_levels),
          cld = "a"
        ))
      }
      
      L_raw <- multcompView::multcompLetters(pvec, threshold = alpha)$Letters
      
      # Canonicalize to force Unburned -> "a" (left-to-right)
      L <- canonicalize_cld_letters(L_raw, group_levels = group_levels)
      
      # Optional: sort within-string letters ("ba" -> "ab")
      L <- vapply(
        L,
        function(s) paste0(sort(strsplit(s, "", fixed = TRUE)[[1]]), collapse = ""),
        FUN.VALUE = character(1)
      )
      
      out <- tibble::tibble(
        sev_group = factor(group_levels, levels = group_levels),
        cld = unname(L[group_levels])
      )
      
      out$cld[is.na(out$cld)] <- "a"
      out
    }) %>%
    dplyr::ungroup() %>%
    dplyr::rename(!!ysf_col := .data[[ysf_col]])
}

# Use the SAME tukey_results_250 object already built:
cld_250 <- build_cld_from_tukey(
  tukey_tbl     = tukey_results_250,
  alpha         = alpha,
  group_levels  = sev_levels
)

# Mean + 95% CI per YSF x severity
mean_ci_tbl <- sampled_250 %>%
  dplyr::group_by(years_since_fire, sev_group) %>%
  dplyr::summarise(
    n     = dplyr::n(),
    mean  = mean(delta_ndvi_min, na.rm = TRUE),
    sd    = sd(delta_ndvi_min, na.rm = TRUE),
    se    = sd / sqrt(n),
    tcrit = qt(0.975, df = n - 1),
    ci_low  = mean - tcrit * se,
    ci_high = mean + tcrit * se,
    .groups = "drop"
  )

cld_pos_250 <- sampled_250 %>%
  dplyr::group_by(years_since_fire, sev_group) %>%
  dplyr::summarise(
    y_pos = max(delta_ndvi_min, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(cld_250, by = c("years_since_fire", "sev_group")) %>%
  dplyr::mutate(y_pos = y_pos + 0.03)

# BOXPLOT 
# Dynamic title text
ysf_vals <- sampled_250 %>%
  dplyr::distinct(years_since_fire) %>%
  dplyr::arrange(years_since_fire) %>%
  dplyr::pull(years_since_fire)

title_text <- paste0("ΔNDVI at ", paste(ysf_vals, collapse = ", "), " Years Post-Fire")

### 250 SUBSAMPLE BOX PLOT:
p_250 <- ggplot(
  sampled_250,
  aes(x = factor(years_since_fire), y = delta_ndvi_min, fill = sev_group)
) +
  geom_hline(
    yintercept = 0,
    color = "darkgrey",
    linewidth = 0.8
  ) +
  stat_summary(
    fun.data = function(y) {
      y <- y[is.finite(y)]
      if (length(y) == 0) {
        return(data.frame(ymin = NA, lower = NA, middle = NA, upper = NA, ymax = NA))
      }
      
      q1  <- unname(quantile(y, 0.25, na.rm = TRUE))
      q3  <- unname(quantile(y, 0.75, na.rm = TRUE))
      iqr <- q3 - q1
      
      low_fence  <- q1 - 1.5 * iqr
      high_fence <- q3 + 1.5 * iqr
      
      whisk_low  <- min(y[y >= low_fence],  na.rm = TRUE)
      whisk_high <- max(y[y <= high_fence], na.rm = TRUE)
      
      data.frame(
        ymin   = whisk_low,
        lower  = q1,
        middle = mean(y, na.rm = TRUE),
        upper  = q3,
        ymax   = whisk_high
      )
    },
    geom      = "boxplot",
    alpha     = 0.85,
    linewidth = 0.6,
    width     = 0.75,
    position  = position_dodge(width = 0.85)
  ) +
  geom_errorbar(
    data = mean_ci_tbl,
    aes(
      x = factor(years_since_fire),
      y = mean,
      ymin = ci_low,
      ymax = ci_high,
      group = sev_group
    ),
    position = position_dodge(width = 0.85),
    width = 0.18,
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  geom_text(
    data = cld_pos_250,
    aes(
      x = factor(years_since_fire),
      y = y_pos,
      label = cld,
      group = sev_group
    ),
    position = position_dodge(width = 0.85),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold",
    family = "Times New Roman"
  ) +
  scale_fill_manual(
    values = c(
      "Unburned" = "gray50",
      "Low"      = "yellow",
      "Moderate" = "orange",
      "High"     = "firebrick"
    ),
    name = "Severity Class"
  ) +
  theme_minimal(base_family = "Times New Roman", base_size = 12) +
  labs(
    title    = title_text,
    subtitle = "250 pixels per severity class",
    x        = "Years Since Fire",
    y        = "ΔNDVI"
  ) +
  theme(
    plot.title    = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title.x  = element_text(size = 13, face = "bold"),
    axis.title.y  = element_text(size = 13, face = "bold"),
    axis.text.x   = element_text(size = 11, face = "bold"),
    axis.text.y   = element_text(size = 11, face = "bold"),
    legend.title  = element_text(size = 13, face = "bold"),
    legend.text   = element_text(size = 12, face = "bold"),
    
    panel.grid = element_blank(),
    
    axis.line.x = element_line(color = "black", linewidth = 0.6),
    axis.line.y = element_line(color = "black", linewidth = 0.6),
    axis.ticks  = element_line(color = "black", linewidth = 0.5)
  )
  

print(p_250)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# BUILDING UP "250" BOXPLOTS FOR SFS PRESENTATION 2026: --------------------------

sev_cols <- c(
  "Unburned" = "gray50",
  "Low"      = "yellow",
  "Moderate" = "orange",
  "High"     = "firebrick"
)

make_p250_build <- function(
    sev_to_show = character(0),
    show_hline  = FALSE,
    show_boxes  = FALSE,
    show_ci     = FALSE,
    show_cld    = FALSE,
    show_legend = TRUE
) {
  
  sampled_plot <- sampled_250 %>%
    dplyr::mutate(
      sev_group = factor(sev_group, levels = sev_levels),
      years_since_fire = factor(years_since_fire, levels = sort(unique(sampled_250$years_since_fire)))
    )
  
  mean_ci_plot <- mean_ci_tbl %>%
    dplyr::mutate(
      sev_group = factor(sev_group, levels = sev_levels),
      years_since_fire = factor(years_since_fire, levels = sort(unique(sampled_250$years_since_fire)))
    )
  
  cld_plot <- cld_pos_250 %>%
    dplyr::mutate(
      sev_group = factor(sev_group, levels = sev_levels),
      years_since_fire = factor(years_since_fire, levels = sort(unique(sampled_250$years_since_fire)))
    )
  
  p <- ggplot(
    sampled_plot,
    aes(
      x = years_since_fire,
      y = delta_ndvi_min
    )
  ) +
    scale_x_discrete(drop = FALSE) +
    labs(
      title    = title_text,
      subtitle = "250 pixels per severity class",
      x        = "Years Since Fire",
      y        = "ΔNDVI"
    ) +
    theme_minimal(base_family = "Aptos Display", base_size = 12) +
    theme(
      plot.title    = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title.x  = element_text(size = 13, face = "bold"),
      axis.title.y  = element_text(size = 13, face = "bold"),
      axis.text.x   = element_text(size = 11, face = "bold"),
      axis.text.y   = element_text(size = 11, face = "bold"),
      legend.title  = element_text(size = 13, face = "bold"),
      legend.text   = element_text(size = 12, face = "bold"),
      panel.grid    = element_blank(),
      axis.line.x   = element_line(color = "black", linewidth = 0.6),
      axis.line.y   = element_line(color = "black", linewidth = 0.6),
      axis.ticks    = element_line(color = "black", linewidth = 0.5),
      legend.position = "right"
    )
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  if (show_hline) {
    p <- p +
      geom_hline(
        yintercept = 0,
        color = "darkgrey",
        linewidth = 0.8
      )
  }
  
  if (show_boxes) {
    p <- p +
      stat_summary(
        data = sampled_plot %>%
          dplyr::mutate(
            alpha_group  = ifelse(sev_group %in% sev_to_show, 0.85, 0),
            box_line_col = ifelse(sev_group %in% sev_to_show, "black", NA)
          ),
        aes(
          fill   = sev_group,
          alpha  = alpha_group,
          colour = box_line_col,
          group  = interaction(years_since_fire, sev_group)
        ),
        fun.data = function(y) {
          y <- y[is.finite(y)]
          
          if (length(y) == 0) {
            return(data.frame(
              ymin = NA,
              lower = NA,
              middle = NA,
              upper = NA,
              ymax = NA
            ))
          }
          
          q1  <- unname(quantile(y, 0.25, na.rm = TRUE))
          q3  <- unname(quantile(y, 0.75, na.rm = TRUE))
          iqr <- q3 - q1
          
          low_fence  <- q1 - 1.5 * iqr
          high_fence <- q3 + 1.5 * iqr
          
          whisk_low  <- min(y[y >= low_fence],  na.rm = TRUE)
          whisk_high <- max(y[y <= high_fence], na.rm = TRUE)
          
          data.frame(
            ymin   = whisk_low,
            lower  = q1,
            middle = mean(y, na.rm = TRUE),
            upper  = q3,
            ymax   = whisk_high
          )
        },
        geom      = "boxplot",
        linewidth = 0.6,
        width     = 0.75,
        position  = position_dodge(width = 0.85)
      ) +
      scale_alpha_identity() +
      scale_colour_identity()
  }
  
  if (show_ci) {
    p <- p +
      geom_errorbar(
        data = mean_ci_plot %>%
          dplyr::mutate(
            alpha_group = ifelse(sev_group %in% sev_to_show, 1, 0)
          ),
        aes(
          x = years_since_fire,
          y = mean,
          ymin = ci_low,
          ymax = ci_high,
          alpha = alpha_group,
          group = interaction(years_since_fire, sev_group)
        ),
        position = position_dodge(width = 0.85),
        width = 0.18,
        inherit.aes = FALSE,
        linewidth = 0.6
      ) +
      scale_alpha_identity()
  }
  
  if (show_cld) {
    p <- p +
      geom_text(
        data = cld_plot %>%
          dplyr::mutate(
            cld_show = ifelse(sev_group %in% sev_to_show, cld, "")
          ),
        aes(
          x = years_since_fire,
          y = y_pos,
          label = cld_show,
          group = interaction(years_since_fire, sev_group)
        ),
        position = position_dodge(width = 0.85),
        inherit.aes = FALSE,
        size = 5,
        fontface = "bold",
        family = "Aptos Display"
      )
  }
  
  if (show_boxes || show_ci || show_cld) {
    p <- p +
      scale_fill_manual(
        values = sev_cols,
        limits = sev_levels,
        drop = FALSE,
        name = "Severity Class"
      )
  }
  
  p
}


# CREATE 6 BUILD-UP PLOTS
p_250_step1_axes <- make_p250_build(
  show_legend = FALSE
)

p_250_step2_hline <- make_p250_build(
  show_hline  = TRUE,
  show_legend = FALSE
)

p_250_step3_unburned <- make_p250_build(
  sev_to_show = c("Unburned"),
  show_hline  = TRUE,
  show_boxes  = TRUE,
  show_ci     = TRUE,
  show_legend = TRUE
)

p_250_step4_unburned_low <- make_p250_build(
  sev_to_show = c("Unburned", "Low"),
  show_hline  = TRUE,
  show_boxes  = TRUE,
  show_ci     = TRUE,
  show_cld    = TRUE,
  show_legend = TRUE
)

p_250_step5_all_boxes <- make_p250_build(
  sev_to_show = c("Unburned", "Low", "Moderate", "High"),
  show_hline  = TRUE,
  show_boxes  = TRUE,
  show_ci     = TRUE,
  show_legend = TRUE
)

p_250_step6_final_cld <- make_p250_build(
  sev_to_show = c("Unburned", "Low", "Moderate", "High"),
  show_hline  = TRUE,
  show_boxes  = TRUE,
  show_ci     = TRUE,
  show_cld    = TRUE,
  show_legend = TRUE
)

# Print one at a time
p_250_step1_axes
p_250_step2_hline
p_250_step3_unburned
p_250_step4_unburned_low
p_250_step5_all_boxes
p_250_step6_final_cld

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# RANDOM % OF EACH SEV_GRP SUBSET: ANOVA, Tukey, N TABLE, BOXPLOT
# ------------------------------------------------------------------------------
min_n_per_cell <- 2  # minimum number of pixels per YSF × sev_group cell

# Stratified random sample by years_since_fire × sev_group
sampled_pct <- snapshot %>%
  dplyr::group_by(years_since_fire, sev_group) %>%
  dplyr::group_modify(~ dplyr::slice_sample(
    .x,
    n = max(min_n_per_cell, ceiling(pct * nrow(.x)))
  )) %>%
  dplyr::ungroup()

# N table for this percentage-based subset
n_by_year_pct <- make_n_by_year(sampled_pct)

# ANOVA summary for percentage-based subset
anova_results_pct <- build_anova_summary(
  df             = sampled_pct,
  n_by_year_tbl   = n_by_year_pct,
  N_total_name    = "N_total_pct",
  N_by_group_name = "N_by_group_pct"
)

# Tukey summary for percentage-based subset
tukey_results_pct <- build_tukey_summary(
  df            = sampled_pct,
  n_by_year_tbl = n_by_year_pct
)

# Mean + 95% CI per YSF x severity (PCT)
mean_ci_pct <- sampled_pct %>%
  dplyr::group_by(years_since_fire, sev_group) %>%
  dplyr::summarise(
    n     = dplyr::n(),
    mean  = mean(delta_ndvi_min, na.rm = TRUE),
    sd    = sd(delta_ndvi_min, na.rm = TRUE),
    se    = sd / sqrt(n),
    tcrit = qt(0.975, df = n - 1),
    ci_low  = mean - tcrit * se,
    ci_high = mean + tcrit * se,
    .groups = "drop"
  )

# CLD letters (PCT)
cld_pct <- build_cld_tbl(sampled_pct)

# Y positions for CLD letters (PCT)
cld_pos_pct <- sampled_pct %>%
  dplyr::group_by(years_since_fire, sev_group) %>%
  dplyr::summarise(
    y_pos = quantile(delta_ndvi_min, probs = 0.98, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(cld_pct, by = c("years_since_fire", "sev_group")) %>%
  dplyr::mutate(
    y_pos = y_pos + 0.03
  )

# Boxplot for percentage-based subset
sample_sizes_pct <- sampled_pct %>%
  dplyr::group_by(years_since_fire, sev_group) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(label = paste0("n = ", n))

N_min_pct <- min(sample_sizes_pct$n, na.rm = TRUE)
N_max_pct <- max(sample_sizes_pct$n, na.rm = TRUE)

pct_label <- scales::percent(pct, accuracy = 1)

N_caption_pct <- paste0(
  pct_label, " stratified sample (per YSF × severity). N ranges from ",
  scales::comma(N_min_pct), "–", scales::comma(N_max_pct), "."
)

show_N_pct <- FALSE

p_pct <- ggplot(
  sampled_pct,
  aes(x = factor(years_since_fire), y = delta_ndvi_min, fill = sev_group)
) +
  geom_boxplot(
    outlier.size = 0.6,
    alpha        = 0.85,
    linewidth    = 0.6,
    position     = position_dodge(width = 0.85)
  ) +
  # 95% CI error bars for the mean
  geom_errorbar(
    data = mean_ci_pct,
    aes(
      x = factor(years_since_fire),
      y = mean,
      ymin = ci_low,
      ymax = ci_high,
      group = sev_group
    ),
    position = position_dodge(width = 0.85),
    width = 0.18,
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  # mean points
  geom_point(
    data = mean_ci_pct,
    aes(
      x = factor(years_since_fire),
      y = mean,
      group = sev_group
    ),
    position = position_dodge(width = 0.85),
    inherit.aes = FALSE,
    size = 2.2,
    shape = 21
  ) +
  # CLD letters
  geom_text(
    data = cld_pos_pct,
    aes(
      x = factor(years_since_fire),
      y = y_pos,
      label = cld,
      group = sev_group
    ),
    position = position_dodge(width = 0.85),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  {
    if (show_N_pct) {
      geom_text(
        data        = sample_sizes_pct,
        aes(x = factor(years_since_fire), y = 0.75, label = label, group = sev_group),
        position    = position_dodge(width = 0.85),
        inherit.aes = FALSE,
        vjust       = -0.5,
        size        = 3.8
      )
    }
  } +
  scale_fill_manual(
    values = c(
      "Unburned" = "gray50",
      "Low"      = "yellow",
      "Moderate" = "orange",
      "High"     = "firebrick"
    ),
    name = "Severity Class"
  ) +
  theme_minimal(base_family = "Times New Roman", base_size = 12) +
  labs(
    title   = paste0("ΔNDVI by Severity Across Time (", pct_label, " Subset)"),
    x       = "Years Since Fire",
    y       = "ΔNDVI (Postfire - Prefire)",
    caption = N_caption_pct
  ) +
  theme(
    plot.title    = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14),
    axis.title.x  = element_text(size = 13, face = "bold"),
    axis.title.y  = element_text(size = 13, face = "bold"),
    axis.text.x   = element_text(size = 11, face = "bold", angle = 0,
                                 vjust = 0.5, hjust = 0.5),
    axis.text.y   = element_text(size = 11, face = "bold"),
    legend.title  = element_text(size = 13, face = "bold"),
    legend.text   = element_text(size = 12, face = "bold"),
    plot.caption  = element_text(size = 12)
  )

print(p_pct)
# ------------------------------------------------------------------------------




################################################################################
# SAME-PIXEL COHORT OVER TIME (selected at 5YSF; shown at 1, 3, 5 YSF)
# - Select up to 250 pixel_IDs per severity class from the 5YSF rows
# - Keep ONLY those same pixel_IDs at 1YSF and 3YSF (must have all 3 timepoints)
################################################################################

# --- ASSUMES you already have:
# data_long_complete with columns: pixel_ID, year, ref_year, sev_group, delta_ndvi_min
# and your helper fns: safe_anova_tbl(), safe_tukey_tbl(), make_n_by_year(),
# build_anova_summary(), build_tukey_summary()
# If you already defined build_cld_tbl() earlier, replace it with the updated version below.

# -------------------------------------------------------------------------
# Updated CLD helper
# -------------------------------------------------------------------------
build_cld_tbl <- function(df,
                          ysf_col       = "years_since_fire",
                          group_col     = "sev_group",
                          response_col  = "delta_ndvi_min",
                          group_levels  = c("Unburned", "Low", "Moderate", "High")) {
  
  stopifnot(all(c(ysf_col, group_col, response_col) %in% names(df)))
  
  df2 <- df %>%
    dplyr::filter(!is.na(.data[[response_col]])) %>%
    dplyr::mutate(
      !!group_col := factor(.data[[group_col]], levels = group_levels)
    ) %>%
    droplevels()
  
  df2 %>%
    dplyr::group_by(.data[[ysf_col]]) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      tukey_tbl = purrr::map(data, ~ safe_tukey_tbl(.x)),
      letters   = purrr::map(tukey_tbl, ~ {
        tt <- .x %>% dplyr::filter(!is.na(p.adj))
        
        # If Tukey couldn't be computed / no comparisons, fall back to all "a"
        if (nrow(tt) == 0) {
          levs <- levels(df2[[group_col]])
          out  <- rep("a", length(levs))
          names(out) <- levs
          return(out)
        }
        
        pvec <- tt$p.adj
        names(pvec) <- paste(tt$group1, tt$group2, sep = "-")
        multcompView::multcompLetters(pvec, threshold = alpha)$Letters
      })
    ) %>%
    dplyr::select(all_of(ysf_col), letters) %>%
    tidyr::unnest_wider(letters, names_sep = "_") %>%
    tidyr::pivot_longer(
      cols      = starts_with("letters_"),
      names_to  = group_col,
      values_to = "cld"
    ) %>%
    dplyr::mutate(
      !!group_col := gsub("^letters_", "", .data[[group_col]]),
      !!group_col := factor(.data[[group_col]], levels = group_levels)
    )
}

# ------------------------------------------------------------------------------
# 1) Build 1/3/5 snapshot from full dataset
# ------------------------------------------------------------------------------
#snapshot_all <- data_long_complete %>%
  dplyr::mutate(
    years_since_fire = year - ref_year,
    sev_group = factor(sev_group, levels = c("Unburned", "Low", "Moderate", "High"))
  ) %>%
  dplyr::filter(years_since_fire %in% c(1, 3, 5)) %>%
  dplyr::filter(!is.na(delta_ndvi_min)) %>%
  droplevels()

# ------------------------------------------------------------------------------
# 2) Restrict to pixel_IDs that have all 3 timepoints (1,3,5)
# ------------------------------------------------------------------------------
eligible_pixels <- snapshot_all %>%
  dplyr::distinct(pixel_ID, years_since_fire) %>%
  dplyr::count(pixel_ID, name = "n_timepoints") %>%
  dplyr::filter(n_timepoints == 3L) %>%
  dplyr::pull(pixel_ID)

snapshot_complete <- snapshot_all %>%
  dplyr::filter(pixel_ID %in% eligible_pixels) %>%
  droplevels()

# ------------------------------------------------------------------------------
# 3) Choose the cohort ONCE from the 5YSF rows: up to 250 pixel_IDs per severity
#    (This is the key change: selection occurs only at YSF=5.)
# ------------------------------------------------------------------------------
target_n_per_class <- 250

cohort_ids_5ysf <- snapshot_complete %>%
  dplyr::filter(years_since_fire == 5) %>%
  dplyr::distinct(pixel_ID, sev_group) %>%
  dplyr::group_by(sev_group) %>%
  dplyr::group_modify(~ dplyr::slice_sample(.x, n = min(target_n_per_class, nrow(.x)))) %>%
  dplyr::ungroup() %>%
  dplyr::pull(pixel_ID)

# Keep this same cohort at 1, 3, and 5 YSF
sampled_same_cohort <- snapshot_complete %>%
  dplyr::filter(pixel_ID %in% cohort_ids_5ysf) %>%
  droplevels()

# ------------------------------------------------------------------------------
# 4) ANOVA / Tukey / Ns (same structure as above workflow)
# ------------------------------------------------------------------------------
n_by_year_same <- make_n_by_year(sampled_same_cohort)

anova_results_same <- build_anova_summary(
  df            = sampled_same_cohort,
  n_by_year_tbl = n_by_year_same
)

tukey_results_same <- build_tukey_summary(
  df            = sampled_same_cohort,
  n_by_year_tbl = n_by_year_same
)

# Mean + 95% CI per YSF x severity (same as before)
mean_ci_same <- sampled_same_cohort %>%
  dplyr::group_by(years_since_fire, sev_group) %>%
  dplyr::summarise(
    n     = dplyr::n(),
    mean  = mean(delta_ndvi_min, na.rm = TRUE),
    sd    = sd(delta_ndvi_min, na.rm = TRUE),
    se    = sd / sqrt(n),
    tcrit = qt(0.975, df = n - 1),
    ci_low  = mean - tcrit * se,
    ci_high = mean + tcrit * se,
    .groups = "drop"
  )

# CLD letters per YSF (Tukey within each YSF)
cld_same <- build_cld_tbl(sampled_same_cohort)

# Y positions for CLD labels
cld_pos_same <- sampled_same_cohort %>%
  dplyr::group_by(years_since_fire, sev_group) %>%
  dplyr::summarise(
    y_pos = max(delta_ndvi_min, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(cld_same, by = c("years_since_fire", "sev_group")) %>%
  dplyr::mutate(
    y_pos = y_pos + 0.03  # tweak if needed
  )

# ------------------------------------------------------------------------------
# 5) Boxplot: same cohort shown across 1, 3, 5
# ------------------------------------------------------------------------------
p_same_cohort <- ggplot(
  sampled_same_cohort,
  aes(x = factor(years_since_fire), y = delta_ndvi_min, fill = sev_group)
) +
  geom_boxplot(
    outlier.size = 0.6,
    alpha        = 0.85,
    linewidth    = 0.6,
    position     = position_dodge(width = 0.85)
  ) +
  # 95% CI error bars for the mean
  geom_errorbar(
    data = mean_ci_same,
    aes(
      x = factor(years_since_fire),
      y = mean,
      ymin = ci_low,
      ymax = ci_high,
      group = sev_group
    ),
    position = position_dodge(width = 0.85),
    width = 0.18,
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  # mean points
  geom_point(
    data = mean_ci_same,
    aes(
      x = factor(years_since_fire),
      y = mean,
      group = sev_group
    ),
    position = position_dodge(width = 0.85),
    inherit.aes = FALSE,
    size = 2.2,
    shape = 21
  ) +
  # CLD letters
  geom_text(
    data = cld_pos_same,
    aes(
      x = factor(years_since_fire),
      y = y_pos,
      label = cld,
      group = sev_group
    ),
    position = position_dodge(width = 0.85),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Unburned" = "gray50",
      "Low"      = "yellow",
      "Moderate" = "orange",
      "High"     = "firebrick"
    ),
    name = "Severity Class"
  ) +
  theme_minimal(base_family = "Times New Roman", base_size = 12) +
  labs(
    title = "ΔNDVI at 1, 3, and 5 Years Post Fire by Severity",
    subtitle = "Same pixel cohort over time (250 pixels per severity class, selected from 5YSF dataset)",
    x = "Years Since Fire",
    y = "ΔNDVI (Postfire - Prefire)"
  ) +
  theme(
    plot.title    = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    axis.title.x  = element_text(size = 13, face = "bold"),
    axis.title.y  = element_text(size = 13, face = "bold"),
    axis.text.x   = element_text(size = 11, face = "bold"),
    axis.text.y   = element_text(size = 11, face = "bold"),
    legend.title  = element_text(size = 13, face = "bold"),
    legend.text   = element_text(size = 12, face = "bold")
  )

print(p_same_cohort)

################################################################################




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
################################################################################
#          -------------------------------------------                         #
####     RESEARCH QUESTION 2: TEMPORAL RECOVERY PATTERNS     ###################
#          -------------------------------------------                         #
################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# ==========================================================
# FILTERING PIXELS BY THEIR AVAILABLE TEMPORAL RECORD FOR EACH YSF: 
# ==========================================================

## ================================================================
## Part 1A: MASTER SWITCHES + TEMPORAL RECORD FILTERING  -----------------
# Preconditions   
stopifnot(all(c("pixel_ID",
                "year",
                "ref_year",
                "ndvi_postfire",
                "ndvi_prefire_3yr_min",
                "sev_group",
                "fire_name",
                "fire_size_ha") %in% names(data_long_complete)))

# --- Parameters (MASTER SWITCHES) 
cutoff_year        <- 2024
ysf_targets        <- c(5, 10, 15, 20) # YSF targets (5, 10, 15, 20)
sev_levels         <- c("Unburned","Low","Moderate","High")
recovery_threshold <- 1.0   # recovery threshold to meet 
recovery_metric <- "min" # Which delta NDVI metric is used for the recovery ratio?
                          #  Options: 3-yr pre-fire "avg", "med", "min"

metric_label <- dplyr::case_when(
  recovery_metric == "avg"    ~ "Mean NDVI",
  recovery_metric == "med"    ~ "Median NDVI",
  recovery_metric == "min"    ~ "Minimum NDVI",
  TRUE                        ~ "Pre-Fire NDVI"
)

# --- Per-pixel cutoff year (how many yrs of NDVI available up to 2024?) 
has_fire_name <- "fire_name" %in% names(data_long_complete)

pixel_meta <- data_long_complete %>%
  dplyr::filter(year <= cutoff_year) %>%
  dplyr::arrange(pixel_ID, year) %>%   
  dplyr::group_by(pixel_ID) %>%
  dplyr::summarise(
    ref_year        = dplyr::first(ref_year),
    sev_group       = dplyr::first(sev_group),
    fire_name       = if (has_fire_name) dplyr::first(.data$fire_name) else NA_character_,
    last_year_avail = max(year, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    max_ysf_avail = pmax(0L, last_year_avail - ref_year)
  )

## ================================================================
## Part 1B: Pixel recovery flags, 1-30ysf + Calculate years_to_recovery ------------
max_flag_ysf <- 30L  # 30YSF cutoff for table. 

# 1) Long form: one row per pixel-YSF with monotonic recovery_flag
pixel_recovery_long <- data_long_complete %>%
  dplyr::filter(year <= cutoff_year, years_since_fire >= 1, years_since_fire <= max_flag_ysf) %>%
  dplyr::mutate(
    ndvi_ratio = dplyr::case_when(
      recovery_metric == "avg" ~ ndvi_postfire / ndvi_prefire_3yr_avg,
      recovery_metric == "med" ~ ndvi_postfire / ndvi_prefire_3yr_med,
      recovery_metric == "min" ~ ndvi_postfire / ndvi_prefire_3yr_min,
      TRUE ~ NA_real_
    ),
    # raw per-year flag (not yet monotonic)
    recovery_flag_raw = dplyr::if_else(
      is.finite(ndvi_ratio) & (ndvi_ratio >= recovery_threshold),
      1L, 0L
    )
  ) %>%
  dplyr::select(pixel_ID, years_since_fire, recovery_flag_raw) %>%
  dplyr::distinct(pixel_ID, years_since_fire, .keep_all = TRUE) %>%
  dplyr::group_by(pixel_ID) %>%
  dplyr::mutate(
    # first YSF where raw flag hits 1; NA if never hits 1 in observed record
    years_to_recovery = {
      y <- years_since_fire[recovery_flag_raw == 1L]
      if (length(y) == 0L) NA_integer_ else as.integer(min(y))
    },
    # monotonic (absorbing) flag: once recovered, always recovered (where record exists)
    recovery_flag = dplyr::if_else(
      !is.na(years_to_recovery) & years_since_fire >= years_to_recovery,
      1L, 0L
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(pixel_ID, years_since_fire, recovery_flag, years_to_recovery)

# 2) Wide form: columns rec_ysf1 ... rec_ysf30 with NA where no record exists
pixel_recovery_wide <- pixel_recovery_long %>%
  dplyr::select(pixel_ID, years_since_fire, recovery_flag) %>%
  tidyr::pivot_wider(
    names_from   = years_since_fire,
    values_from  = recovery_flag,
    names_prefix = "rec_ysf"
  ) %>%
  dplyr::mutate(
    dplyr::across(dplyr::starts_with("rec_ysf"), as.integer)
  )

# Ensure every expected YSF column exists (rec_ysf1..rec_ysf30); add missing as NA
needed_cols  <- paste0("rec_ysf", 1:max_flag_ysf)
missing_cols <- setdiff(needed_cols, names(pixel_recovery_wide))
if (length(missing_cols) > 0) {
  pixel_recovery_wide[missing_cols] <- NA_integer_
}
pixel_recovery_wide <- pixel_recovery_wide %>%
  dplyr::select(pixel_ID, dplyr::all_of(needed_cols))

# 3) years_to_recovery: first YSF where recovery_flag == 1; NA if never recovered
pixel_ytr <- pixel_recovery_long %>%
  dplyr::group_by(pixel_ID) %>%
  dplyr::summarise(
    years_to_recovery = dplyr::if_else(
      any(recovery_flag == 1L, na.rm = TRUE),
      as.integer(min(years_since_fire[recovery_flag == 1L], na.rm = TRUE)),
      NA_integer_
    ),
    .groups = "drop"
  )

# 4) Combine flags + ytr into one pixel-level table
pixel_recovery_flags <- pixel_recovery_wide %>%
  dplyr::left_join(pixel_ytr, by = "pixel_ID")

# 5) Add rec_ysf columns onto data_long_complete (every row for a pixel gets same flags/ytr)
# --- remove previously-added recovery columns so reruns don't create .x/.y suffixes
drop_cols <- c("years_to_recovery", paste0("rec_ysf", 1:max_flag_ysf))
data_long_complete <- data_long_complete %>%
  dplyr::select(-dplyr::any_of(drop_cols))
data_long_complete <- data_long_complete %>%
  dplyr::left_join(pixel_recovery_flags, by = "pixel_ID")

# Checks
stopifnot("years_to_recovery" %in% names(data_long_complete))
stopifnot(all(paste0("rec_ysf", 1:30) %in% names(data_long_complete)))
# =================================================================
## ================================================================
## Part 1C: Eligible pixels by YSF cohort + Summaries =========================

# Function: eligible cohort for a given target YSF
# (pixels with >= target_ysf years of NDVI data post-fire)
eligible_for <- function(target_ysf) {
  pixel_meta %>%
    dplyr::filter(max_ysf_avail >= target_ysf, !is.na(sev_group)) %>%
    dplyr::select(pixel_ID, sev_group, fire_name, ref_year)
}
# YSF Summary Tables: 1 concise, 1 detailed
#   - N_pixels_YSF: # of eligible pixels
#   - N_fire_events_YSF: # of burned fire events (excludes "None")
#   - N_ignition_years_YSF: # of distinct ignition years among those fires

####### Concise summary by YSF #######
YSF_summary <- lapply(ysf_targets, function(tgt) {
  elig <- eligible_for(tgt)
  
  # All eligible pixels (burned + unburned)
  N_pixels_YSF <- dplyr::n_distinct(elig$pixel_ID)
  
  # Counting fires and ignition years (for fire events, exclude "None")
  fires_tbl <- elig %>%
    dplyr::filter(!is.na(fire_name), fire_name != "None") %>%
    dplyr::distinct(fire_name, ref_year)
  
  N_fire_events_YSF    <- nrow(fires_tbl)
  N_ignition_years_YSF <- dplyr::n_distinct(fires_tbl$ref_year)
  
  tibble::tibble(
    YSF                   = tgt,
    N_pixels_YSF          = N_pixels_YSF,
    N_fire_events_YSF     = N_fire_events_YSF,
    N_ignition_years_YSF  = N_ignition_years_YSF
  )
}) %>%
  dplyr::bind_rows()

print(YSF_summary)

# Named vectors for referencing them elsewhere:
N_pixels_by_YSF <- stats::setNames(YSF_summary$N_pixels_YSF,
                                   paste0("N_YSF", YSF_summary$YSF))
N_fires_by_YSF  <- stats::setNames(YSF_summary$N_fire_events_YSF,
                                   paste0("number_of_fires_YSF", YSF_summary$YSF))

####### Detailed summary: fire-by-fire within each YSF #######
YSF_fire_detail <- lapply(ysf_targets, function(tgt) {
  # Eligible pixels for this YSF (burned + unburned)
  elig <- eligible_for(tgt)
  
  # If no eligible pixels at all, return an empty row for clarity
  if (nrow(elig) == 0L) {
    return(
      tibble::tibble(
        YSF               = tgt,
        fire_name         = NA_character_,
        ref_year          = NA_integer_,
        N_pixels_YSF_fire = 0L
      )
    )
  }
  
  # Group by fire_name (including "None") and ref_year
  # If you'd prefer ref_year = NA for unburned, uncomment the mutate below
  # elig <- elig %>%
  #   dplyr::mutate(
  #     ref_year = dplyr::if_else(fire_name == "None" | is.na(fire_name),
  #                               NA_integer_, ref_year)
  #   )
  
  elig %>%
    dplyr::group_by(YSF = tgt, fire_name, ref_year) %>%
    dplyr::summarise(
      N_pixels_YSF_fire = dplyr::n_distinct(pixel_ID),
      .groups = "drop"
    )
}) %>%
  dplyr::bind_rows() %>%
  dplyr::arrange(YSF, ref_year, fire_name)

print(YSF_fire_detail, n = Inf)
## ================================================================
## Part 1D: Maps — Contributing Fires & Ignition Yrs at Each Target YSF ---------
# Join attributes needed for mapping to existing pixel geometries
pixel_map_df <- pixel_geom_sf %>%
  dplyr::left_join(
    data_long_complete %>%
      dplyr::distinct(pixel_ID, sev_group, fire_name, ref_year),
    by = "pixel_ID"
  ) %>%
  dplyr::filter(!is.na(sev_group))   # burned + controls

# Ensure backdrops match CRS of pixel geometries (do this ONCE)
huc12_bg   <- huc12_sf   %>% sf::st_transform(sf::st_crs(pixel_geom_sf))
streams_bg <- streams_sf %>% sf::st_transform(sf::st_crs(pixel_geom_sf))
study_bbox <- sf::st_bbox(huc12_bg)

# Font
library(showtext)
font_add("Times New Roman", "C:/Windows/Fonts/times.ttf")
showtext_auto()

# Build one map per target YSF
maps_by_target <- lapply(ysf_targets, function(tgt) {
  
  # 1) Eligible cohort for this YSF (burned + unburned)
  base_elig <- eligible_for(tgt)   # pixel_ID, sev_group, fire_name, ref_year
  
  # If there are no eligible pixels at this YSF, draw only background
  if (nrow(base_elig) == 0L) {
    message("No eligible pixels found for YSF = ", tgt)
    return(
      ggplot() +
        geom_sf(data = huc12_bg, fill = "grey97", color = "grey70", linewidth = 0.2) +
        geom_sf(data = streams_bg, color = "steelblue4", linewidth = 0.6, alpha = 0.7) +
        coord_sf(
          xlim = c(study_bbox["xmin"], study_bbox["xmax"]),
          ylim = c(study_bbox["ymin"], study_bbox["ymax"]),
          expand = FALSE
        ) +
        labs(
          title = paste0("Contributing Fires & Pixels — Target YSF ≥ ", tgt),
          subtitle = "No eligible pixels for this YSF group",
          x = NULL, y = NULL
        ) +
        theme_minimal(base_size = 12)
    )
  }
  
  # ---- 2) Compute summary metrics INSIDE the map loop ----
  N_pixels_YSF <- nrow(base_elig)
  
  burned_only <- base_elig %>% 
    dplyr::filter(!is.na(fire_name), fire_name != "None")
  
  N_fire_events_YSF <- burned_only %>% 
    dplyr::distinct(fire_name) %>% 
    nrow()
  
  N_ignition_years_YSF <- burned_only %>%
    dplyr::distinct(ref_year) %>%
    nrow()
  
  # ---- 3) Attach geometry + ignition year ----
  elig <- base_elig %>%
    dplyr::mutate(
      ignition_year = factor(ref_year, levels = sort(unique(ref_year)))
    ) %>%
    dplyr::left_join(
      pixel_geom_sf %>% dplyr::select(pixel_ID, geometry),
      by = "pixel_ID"
    ) %>%
    sf::st_as_sf()
  
  ignition_levels <- levels(elig$ignition_year)
  n_ignition_years <- length(ignition_levels)
  
  # ---- 4) Palette ----
  pal <- viridis::viridis(n_ignition_years, option = "D", direction = 1)
  names(pal) <- ignition_levels
  
  # ---- 5) Plot ----
  ggplot() +
    geom_sf(data = huc12_bg, fill = "grey97", color = "grey70", linewidth = 0.2) +
    geom_sf(data = streams_bg, color = "steelblue4", linewidth = 0.6, alpha = 0.7) +
    
    # Pixel aesthetics
    geom_sf(    
      data  = elig,
      aes(fill = ignition_year),
      size  = 3.5,
      alpha = 0.75,
      shape = 21,
      stroke = 0
    ) +
    
    scale_fill_manual(
      name         = "Ignition Year",
      values       = pal,
      drop         = FALSE,
      na.translate = FALSE
    ) +
    
    coord_sf(
      xlim   = c(study_bbox["xmin"], study_bbox["xmax"]),
      ylim   = c(study_bbox["ymin"], study_bbox["ymax"]),
      expand = FALSE
    ) +
    
    labs(
      title = paste0("Contributing Fires & Pixels — Target YSF ≥ ", tgt),
      subtitle = paste0(
        "N(pixels) = ", scales::comma(N_pixels_YSF),
        " • N(fire events) = ", N_fire_events_YSF,
        " • N(ignition years) = ", N_ignition_years_YSF
      ),
      x = NULL, y = NULL
    ) +
    theme_minimal(base_family = "Times New Roman", base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      
      legend.position  = "right",
      legend.title     = element_text(face = "bold", size = 15, family = "Times New Roman"),
      legend.text      = element_text(size = 13, family = "Times New Roman"),
      plot.title       = element_text(face = "bold", size = 18, family = "Times New Roman"),
      plot.subtitle    = element_text(size = 14, family = "Times New Roman"),
      axis.text        = element_text(family = "Times New Roman"),
      axis.title       = element_text(family = "Times New Roman")
    )
})

# Print all target-specific maps
for (m in maps_by_target) print(m)
################################################################################

## ==============================================================
## Part 2A: Filter/Join. pixel-level years_to_recovery by YSF cohort
## ==============================================================
# One record per pixel with the already-computed recovery time from Part 1B
pixel_ytr_master <- data_long_complete %>%
  dplyr::distinct(pixel_ID, sev_group, fire_name, ref_year, years_to_recovery, fire_size_ha) %>%
  dplyr::rename(yrs_to_recovery = years_to_recovery)  # <-- standardize name for downstream code

# Build cohort-specific pixel tables (eligibility is the only thing that changes by YSF target)
pixel_recovery_by_YSF <- lapply(ysf_targets, function(target_ysf) {
  eligible_for(target_ysf) %>%
    dplyr::left_join(pixel_ytr_master, by = c("pixel_ID","sev_group","fire_name","ref_year")) %>%
    dplyr::mutate(YSF = target_ysf) %>%
    dplyr::relocate(YSF, pixel_ID, sev_group, fire_name, ref_year, yrs_to_recovery)
}) %>%
  dplyr::bind_rows()

stopifnot("yrs_to_recovery" %in% names(pixel_recovery_by_YSF))

# ===============================================================
## ==============================================================
## Part 2B: Cumulative recovery by YSF Cohort
## "Given a recovery observation window of X years, 
##     what fraction of pixels of each severity have recovered by each year?"
## ============================================================

# Helper: build cumulative recovery table for one target YSF
build_recovery_summary_for_target <- function(target_ysf, ysf_min = 1L) {
  
  eligible <- pixel_recovery_by_YSF %>%
    dplyr::filter(YSF == target_ysf, sev_group %in% sev_levels)
  
  ysf_seq <- seq(ysf_min, target_ysf, by = 1L)
  
  build_sev_table <- function(sev_label, short_label) {
    
    sub <- eligible %>% dplyr::filter(sev_group == sev_label)
    
    N_total <- dplyr::n_distinct(sub$pixel_ID)
    
    n_rec_vec <- vapply(
      ysf_seq,
      function(k) {
        sum(!is.na(sub$yrs_to_recovery) & sub$yrs_to_recovery <= k)
      },
      integer(1)
    )
    
    tibble::tibble(
      YSF = ysf_seq,
      !!paste0("N_", short_label) := N_total,
      !!paste0("N_Recovered_", short_label) := n_rec_vec,
      !!paste0("Prop_Recovered_", short_label) := n_rec_vec / N_total
    )
  }
  
  tab_unb  <- build_sev_table("Unburned",  "Unb")
  tab_low  <- build_sev_table("Low",       "Low")
  tab_mod  <- build_sev_table("Moderate",  "Mod")
  tab_high <- build_sev_table("High",      "High")
  
  tab_unb %>%
    dplyr::left_join(tab_low,  by = "YSF") %>%
    dplyr::left_join(tab_mod,  by = "YSF") %>%
    dplyr::left_join(tab_high, by = "YSF") %>%
    dplyr::select(
      YSF,
      Prop_Recovered_Unb,  N_Unb,  N_Recovered_Unb,
      Prop_Recovered_Low,  N_Low,  N_Recovered_Low,
      Prop_Recovered_Mod,  N_Mod,  N_Recovered_Mod,
      Prop_Recovered_High, N_High, N_Recovered_High
    ) %>%
    dplyr::mutate(YSF_target = target_ysf, .before = 1)
}

recovery_summaries_by_target <- lapply(ysf_targets, build_recovery_summary_for_target)

# Name the list elements for convenience (e.g., "YSF5", "YSF10", ...)
names(recovery_summaries_by_target) <- paste0("YSF", ysf_targets)

# Example: inspect the table for YSF = X
if ("YSF10" %in% names(recovery_summaries_by_target)) {
  cat("\n=====================================================\n")
  cat("Cumulative Recovery by Target YSF — YSF-eligible cohort\n")
  cat("Columns are per severity group: proportion, N eligible, N recovered\n")
  cat("=====================================================\n\n")
  
  print(recovery_summaries_by_target[["YSF10"]], n = Inf)
}
## ==============================================================
# ===============================================================
## Part 2C: Bar Graphs — Cumulative recovery by YSF and severity
# ===============================================================
sev_cols <- c(
  Unburned = "gray70",
  Low      = "yellow",
  Moderate = "orange",
  High     = "firebrick"
)

# Safety: require metric_label to exist (fail fast if not)
if (!exists("metric_label")) {
  stop(
    "metric_label is not defined. ",
    "Define recovery_metric and metric_label in Part 1 before making plots.",
    call. = FALSE
  )
}

# Build one plot per target YSF using recovery_summaries_by_target
recovery_plots_by_target <- lapply(names(recovery_summaries_by_target), function(nm) {
  
  summary_tbl <- recovery_summaries_by_target[[nm]]
  target_ysf  <- unique(summary_tbl$YSF_target)
  
  # Look up cohort metadata from YSF_summary (from Part 1)
  cohort_row <- YSF_summary %>%
    dplyr::filter(YSF == target_ysf)
  
  total_N        <- if (nrow(cohort_row)) cohort_row$N_pixels_YSF[1] else NA_integer_
  total_fire     <- if (nrow(cohort_row)) cohort_row$N_fire_events_YSF[1] else NA_integer_
  total_ign_year <- if (nrow(cohort_row)) cohort_row$N_ignition_years_YSF[1] else NA_integer_
  
  # Subtitle formatting (metric-agnostic but correctly labeled)
  subtitle_txt <- paste0(
    "Using Pre-fire ", metric_label,
    " • Pixels w/ Record ≥ ", target_ysf, " Years",
    " • N = ", if (!is.na(total_N)) scales::comma(total_N) else "NA",
    if (!is.na(total_fire))     paste0(" • Number of Fires: ", total_fire) else "",
    if (!is.na(total_ign_year)) paste0(" • Ignition Years: ", total_ign_year) else ""
  )
  
  # -------- LONG FORMAT --------
  summary_long <- summary_tbl %>%
    dplyr::select(
      YSF,
      Prop_Recovered_Unb,
      Prop_Recovered_Low,
      Prop_Recovered_Mod,
      Prop_Recovered_High
    ) %>%
    tidyr::pivot_longer(
      cols      = starts_with("Prop_Recovered_"),
      names_to  = "sev_group",
      values_to = "prop_rec"
    ) %>%
    dplyr::mutate(
      sev_group = dplyr::recode(
        sev_group,
        "Prop_Recovered_Unb"  = "Unburned",
        "Prop_Recovered_Low"  = "Low",
        "Prop_Recovered_Mod"  = "Moderate",
        "Prop_Recovered_High" = "High"
      ),
      sev_group = factor(sev_group, levels = sev_levels)
    )
  
  # -------- PLOT --------
  ggplot(summary_long, aes(x = YSF, y = prop_rec, fill = sev_group)) +
    geom_col(
      position = position_dodge2(
        preserve = "single",
        padding = 0.1
      ),
      color = "black",
      width = 0.85
    ) +
    scale_fill_manual(values = sev_cols, name = "Severity") +
    scale_y_continuous(breaks = seq(0.5, 1, by = 0.1)) +
    coord_cartesian(ylim = c(0.5, 1)) +
    scale_x_continuous(breaks = min(summary_long$YSF):max(summary_long$YSF)) +
    labs(
      title    = paste0("Proportion Recovered by ", target_ysf, "YSF"),
      subtitle = subtitle_txt,
      x = "Years Since Fire",
      y = "Proportion of Pixels Recovered"
    ) +
    theme_minimal(base_family = "Times New Roman") +
    theme(
      text            = element_text(family = "Times New Roman"),
      axis.text.x     = element_text(angle = 45, hjust = 1, family = "Times New Roman"),
      axis.text.y     = element_text(family = "Times New Roman"),
      axis.title.x    = element_text(face = "bold", family = "Times New Roman"),
      axis.title.y    = element_text(face = "bold", family = "Times New Roman"),
      legend.position = "right",
      legend.title    = element_text(face = "bold", family = "Times New Roman"),
      legend.text     = element_text(family = "Times New Roman"),
      plot.title      = element_text(face = "bold", family = "Times New Roman"),
      plot.subtitle   = element_text(family = "Times New Roman")
    )
})

# Name the list elements (e.g., YSF5_plot, YSF10_plot, ...)
names(recovery_plots_by_target) <- paste0(names(recovery_summaries_by_target), "_plot")

# Print all target-specific plots
for (nm in names(recovery_plots_by_target)) {
  print(recovery_plots_by_target[[nm]])
}

# 2C OPTIONAL: Rebuild 1 cohort plot; show bars only at chosen YSF years =======

# (a) Choose one of the cohorts you already created (must match names(recovery_summaries_by_target)) 
choose_cohort <- names(recovery_summaries_by_target)[which.max(ysf_targets)] 
# (b) Choose which YSF years to show bars at (must be within the cohort's available YSF range) 
show_years <- c(5, 10, 15, 20)

# Pull the summary table for the selected cohort and filter to desired years
summary_tbl_sel <- recovery_summaries_by_target[[choose_cohort]]
target_ysf_sel  <- unique(summary_tbl_sel$YSF_target)

summary_long_sel <- summary_tbl_sel %>%
  dplyr::filter(YSF %in% show_years) %>%
  dplyr::select(
    YSF,
    Prop_Recovered_Unb,
    Prop_Recovered_Low,
    Prop_Recovered_Mod,
    Prop_Recovered_High
  ) %>%
  tidyr::pivot_longer(
    cols      = starts_with("Prop_Recovered_"),
    names_to  = "sev_group",
    values_to = "prop_rec"
  ) %>%
  dplyr::mutate(
    sev_group = dplyr::recode(
      sev_group,
      "Prop_Recovered_Unb"  = "Unburned",
      "Prop_Recovered_Low"  = "Low",
      "Prop_Recovered_Mod"  = "Moderate",
      "Prop_Recovered_High" = "High"
    ),
    sev_group = factor(sev_group, levels = sev_levels)
  )

# Recompute subtitle from YSF_summary (same logic as above)
cohort_row_sel <- YSF_summary %>%
  dplyr::filter(YSF == target_ysf_sel)

total_N_sel        <- if (nrow(cohort_row_sel)) cohort_row_sel$N_pixels_YSF[1] else NA_integer_
total_fire_sel     <- if (nrow(cohort_row_sel)) cohort_row_sel$N_fire_events_YSF[1] else NA_integer_
total_ign_year_sel <- if (nrow(cohort_row_sel)) cohort_row_sel$N_ignition_years_YSF[1] else NA_integer_

# Subtitle text:
#subtitle_txt_sel <- paste0(
#  "Using Pre-fire ", metric_label,
#  " • Pixels w/ Record ≥ ", target_ysf_sel, " Years",
#  " • N = ", if (!is.na(total_N_sel)) scales::comma(total_N_sel) else "NA",
#  if (!is.na(total_fire_sel))     paste0(" • Number of Fires: ", total_fire_sel) else "",
#  if (!is.na(total_ign_year_sel)) paste0(" • Ignition Years: ", total_ign_year_sel) else ""
#)

# Plot only the selected years
p_selected_years <- ggplot(summary_long_sel, aes(x = factor(YSF), y = prop_rec, fill = sev_group)
  ) +
  geom_col(
    position = position_dodge2(preserve = "single", padding = 0.1),
    color = "black",
    width = 0.85
  ) +
  scale_fill_manual(values = sev_cols, name = "Severity") +
  scale_y_continuous(breaks = seq(0.5, 1, by = 0.1)) +
  coord_cartesian(ylim = c(0.5, 1)) +
  scale_x_discrete(drop = FALSE) +
  labs(
    title    = paste0("Proportion Recovered by ", target_ysf_sel, "YSF"),
    #subtitle = subtitle_txt_sel,
    x = "Years Since Fire",
    y = "Proportion of Pixels Recovered"
  ) +
  theme_minimal(base_family = "Times New Roman") +
  theme(
    text            = element_text(family = "Times New Roman", size = 14),
    axis.text.x     = element_text(angle = 45, hjust = 1, family = "Times New Roman"),
    axis.text.y     = element_text(family = "Times New Roman"),
    axis.title.x    = element_text(face = "bold", family = "Times New Roman"),
    axis.title.y    = element_text(face = "bold", family = "Times New Roman"),
    legend.position = "right",
    legend.title    = element_text(face = "bold", size = 13, family = "Times New Roman"),
    legend.text     = element_text(family = "Times New Roman", size = 12),
    plot.title      = element_text(face = "bold", family = "Times New Roman")
    #plot.subtitle   = element_text(family = "Times New Roman")
  )

print(p_selected_years)

# ===============================================================
# ===============================================================
## Part 2D: Line Graphs — Cumulative recovery by YSF and severity =======

# Choose which severity curves to show
# Comment out any severity groups you do NOT want plotted
sev_to_plot <- c(
  #"Unburned",
  #"Low",
  #"Moderate",
  #"High"
)

sev_cols <- c(
  Unburned = "gray40",
  Low      = "goldenrod2",
  Moderate = "darkorange",
  High     = "firebrick"
)

# Safety: require metric_label to exist 
if (!exists("metric_label")) {
  stop(
    "metric_label is not defined. ",
    "Define recovery_metric and metric_label in Part 1 before making plots.",
    call. = FALSE
  )
}

# Range of years to show
ysf_min <- 1
ysf_max <- 20

# Build one plot per target YSF using recovery_summaries_by_target
recovery_plots_by_target <- lapply(names(recovery_summaries_by_target), function(nm) {
  
  summary_tbl <- recovery_summaries_by_target[[nm]]
  target_ysf  <- unique(summary_tbl$YSF_target)
  
  # Look up cohort metadata from YSF_summary (from Part 1)
  cohort_row <- YSF_summary %>%
    dplyr::filter(YSF == target_ysf)
  
  total_N        <- if (nrow(cohort_row)) cohort_row$N_pixels_YSF[1] else NA_integer_
  total_fire     <- if (nrow(cohort_row)) cohort_row$N_fire_events_YSF[1] else NA_integer_
  total_ign_year <- if (nrow(cohort_row)) cohort_row$N_ignition_years_YSF[1] else NA_integer_
  
  subtitle_txt <- paste0(
    " Pixels with Record ≥ ", target_ysf, " Years",
    " • n = ", if (!is.na(total_N)) scales::comma(total_N) else "NA"
    #,
    #if (!is.na(total_fire))     paste0(" • Number of Fires: ", total_fire) else "",
    #if (!is.na(total_ign_year)) paste0(" • Ignition Years: ", total_ign_year) else ""
  )
  
  # -------- LONG FORMAT + filter to min/max years 
  summary_long <- summary_tbl %>%
    dplyr::filter(YSF >= ysf_min, YSF <= ysf_max) %>%
    dplyr::select(
      YSF,
      Prop_Recovered_Unb,
      Prop_Recovered_Low,
      Prop_Recovered_Mod,
      Prop_Recovered_High
    ) %>%
    tidyr::pivot_longer(
      cols      = starts_with("Prop_Recovered_"),
      names_to  = "sev_group",
      values_to = "prop_rec"
    ) %>%
    dplyr::mutate(
      sev_group = dplyr::recode(
        sev_group,
        "Prop_Recovered_Unb"  = "Unburned",
        "Prop_Recovered_Low"  = "Low",
        "Prop_Recovered_Mod"  = "Moderate",
        "Prop_Recovered_High" = "High"
      ),
      sev_group = factor(sev_group, levels = sev_to_plot)
    ) %>%
    dplyr::filter(sev_group %in% sev_to_plot
    )
  
  # -------- PLOT (smooth curves) --------
  ggplot(summary_long, aes(x = YSF, y = prop_rec, color = sev_group, group = sev_group)) +
    geom_hline(yintercept = 0.9, linetype = "solid", color = "darkgray", linewidth = 0.8) +
    geom_line(linewidth = 2.0) +
    geom_point(size = 3.2) +
    scale_color_manual(values = sev_cols[sev_to_plot], name = "Severity") +
    scale_y_continuous(breaks = seq(0.2, 1, by = 0.1)) +
    coord_cartesian(ylim = c(0.2, 1)) +
    scale_x_continuous(breaks = ysf_min:ysf_max) +
    labs(
      title    = paste0("Recovery Curves (", ysf_min, "–", ysf_max, " YSF)"),
      subtitle = subtitle_txt,
      x = "Years Since Fire",
      y = "Proportion of Pixels Recovered"
    ) +
    theme_minimal(base_family = "Aptos") +
    theme(
      text            = element_text(family = "Aptos"),
      axis.text.x     = element_text(angle = 45, hjust = 1, family = "Aptos", size = 13),
      axis.text.y     = element_text(angle = 45, hjust = 1, family = "Aptos", size = 13),
      axis.title.x    = element_text(face = "bold", family = "Aptos", size = 16),
      axis.title.y    = element_text(face = "bold", family = "Aptos", size = 16),
      axis.line       = element_line(color = "black", linewidth = 0.6),
      axis.ticks      = element_line(color = "black", linewidth = 0.5),
      legend.position = "right",
      legend.title    = element_text(face = "bold", family = "Aptos", size = 16),
      legend.text     = element_text(family = "Aptos", size = 14),
      plot.title      = element_text(face = "bold", family = "Aptos", size = 18),
      plot.subtitle   = element_text(family = "Aptos", size = 14)
    )
})

names(recovery_plots_by_target) <- paste0(names(recovery_summaries_by_target), "_smooth")

for (nm in names(recovery_plots_by_target)) {
  print(recovery_plots_by_target[[nm]])
}

# ===============================================================
# ===============================================================
### Part 2E: Export unrecovered pixels at 20 YSF as shapefile
# ===============================================================
# ---- Output folder ----
out_dir <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/Shapefiles"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_shp <- file.path(out_dir, "unrecovered_pixels_20YSF.shp")

# ---- Safety checks ----
stopifnot(exists("data_long_complete"))
stopifnot(exists("pixel_geom_sf"))
stopifnot(exists("eligible_for"))

stopifnot("pixel_ID" %in% names(data_long_complete))
stopifnot("rec_ysf20" %in% names(data_long_complete))


# 1) Isolate the 20 YSF eligible cohort

cohort_20ysf <- eligible_for(20) %>%
  dplyr::distinct(pixel_ID, sev_group, fire_name, ref_year)

cat("\n=====================================================\n")
cat("20 YSF eligible cohort\n")
cat("=====================================================\n")
cat("N eligible pixels:", dplyr::n_distinct(cohort_20ysf$pixel_ID), "\n")
cat("Severity counts:\n")
print(cohort_20ysf %>% dplyr::count(sev_group, sort = TRUE))


# 2) Filter to unrecovered pixels based on rec_ysf20
# Pixel-level attributes from data_long_complete.
# This keeps one row per pixel and retains as many non-year-specific columns
# as possible without duplicating pixels.
pixel_attrs_20ysf <- data_long_complete %>%
  dplyr::filter(pixel_ID %in% cohort_20ysf$pixel_ID) %>%
  dplyr::arrange(pixel_ID, year) %>%
  dplyr::group_by(pixel_ID) %>%
  dplyr::summarise(
    dplyr::across(
      .cols = -dplyr::any_of(c("year", "years_since_fire")),
      .fns  = dplyr::first
    ),
    last_year_avail = max(year, na.rm = TRUE),
    max_ysf_avail   = max(years_since_fire, na.rm = TRUE),
    .groups = "drop"
  )

unrecovered_20ysf_attrs <- pixel_attrs_20ysf %>%
  dplyr::filter(rec_ysf20 == 0L)

cat("\n=====================================================\n")
cat("Unrecovered pixels at 20 YSF\n")
cat("=====================================================\n")
cat("N unrecovered pixels:", dplyr::n_distinct(unrecovered_20ysf_attrs$pixel_ID), "\n")
cat("Recovery flag counts:\n")
print(unrecovered_20ysf_attrs %>% dplyr::count(rec_ysf20))

cat("\nRandom example unrecovered pixels:\n")
set.seed(546)

n_examples <- min(5L, nrow(unrecovered_20ysf_attrs))

print(
  unrecovered_20ysf_attrs %>%
    dplyr::select(
      dplyr::any_of(c(
        "pixel_ID", "sev_group", "fire_name", "ref_year",
        "years_to_recovery", "rec_ysf20",
        "fire_size_ha", "last_year_avail", "max_ysf_avail"
      ))
    ) %>%
    dplyr::slice_sample(n = n_examples)
)


# 3) Attach geometry and export shapefile
unrecovered_20ysf_sf <- pixel_geom_sf %>%
  dplyr::select(pixel_ID, geometry) %>%
  dplyr::inner_join(unrecovered_20ysf_attrs, by = "pixel_ID") %>%
  sf::st_as_sf()

cat("\n=====================================================\n")
cat("Columns being extracted to shapefile\n")
cat("=====================================================\n")
print(names(unrecovered_20ysf_sf))

cat("\nGeometry type:\n")
print(sf::st_geometry_type(unrecovered_20ysf_sf, by_geometry = FALSE))

cat("\nCRS:\n")
print(sf::st_crs(unrecovered_20ysf_sf))

cat("\nPreview of shapefile attribute table:\n")
print(
  unrecovered_20ysf_sf %>%
    sf::st_drop_geometry() %>%
    dplyr::slice_head(n = 5)
)

# Drop high-volume climate / time-series columns before shapefile export
drop_patterns <- paste(
  c(
    "^AET_",
    "aetavg_",
    "^ppt_",
    "^pptt_",
    "^ppttot_",
    "^pptavg_",
    "^pptsd_",
    "^pptz_",
    "^SWE_",
    "^swe_peak_",
    "^sweavg_",
    "^swesd_",
    "^swez_",
    "^CWD_",
    "^cwd_",
    "^tmax_",
    "^tmaxmean_",
    "^tmaxavg_",
    "^tmaxsd_",
    "^tmaxz_",
    "^tmean_",
    "^tmeanavg_",
    "^tmeansd_",
    "^tmeanz_"
  ),
  collapse = "|"
)

cols_to_drop <- names(unrecovered_20ysf_sf)[
  grepl(drop_patterns, names(unrecovered_20ysf_sf))
]

cat("\n=====================================================\n")
cat("Dropping high-volume climate columns before export\n")
cat("=====================================================\n")
cat("N columns before drop:", ncol(unrecovered_20ysf_sf), "\n")
cat("N columns to drop:", length(cols_to_drop), "\n")
cat("Drop regex:", drop_patterns, "\n")

unrecovered_20ysf_sf <- unrecovered_20ysf_sf %>%
  dplyr::select(-dplyr::all_of(cols_to_drop))
## Filter out the 6 low-sev and 1 unburned pixel. 
unrecovered_20ysf_sf <- unrecovered_20ysf_sf %>%
  dplyr::filter(sev_group %in% c("Moderate", "High"))

cat("N columns after drop:", ncol(unrecovered_20ysf_sf), "\n")
cat("Remaining columns:\n")
print(names(unrecovered_20ysf_sf))

## Print some example pixels prior to export
print(
  unrecovered_20ysf_sf %>%
    dplyr::slice_sample(n = n_examples)
)

unrecovered_20ysf_sf <- unrecovered_20ysf_sf %>%
  dplyr::select(
    -ndvi_prefire_3yr_avg,
    -ndvi_prefire_3yr_range,
    -ndvi_prefire_3yr_med
  )

out_shp <- file.path(
  out_dir,
  "unrecovered_pixels_20YSF_v3.shp"
)

sf::st_write(
  unrecovered_20ysf_sf,
  out_shp,
  delete_layer = TRUE
)







################################################################################
# ===============================================================
# COHORT ISOLATION + ONE-ROW-PER-PIXEL + RECOVERED-BY-LAST-YEAR FLAG
# Uses above pixel_recovery_by_YSF, and 1-year recovery rule
# Run 1A, 1B, 2A, 2B beforehand
# ===============================================================

stopifnot(exists("pixel_recovery_by_YSF"))
stopifnot(all(c("YSF","pixel_ID","sev_group","fire_name","ref_year","yrs_to_recovery", "fire_size_ha") %in% names(pixel_recovery_by_YSF)))

cohorts_to_build <- c(10L, 15L, 20L)

build_clean_cohort <- function(target_ysf) {
  
  cohort_raw <- pixel_recovery_by_YSF %>%
    dplyr::filter(YSF == target_ysf) %>%
    dplyr::select(YSF, pixel_ID, sev_group, fire_name, ref_year, yrs_to_recovery)
  
  # Ensure each pixel appears only once in this cohort
  cohort_one_row <- cohort_raw %>%
    dplyr::group_by(pixel_ID) %>%
    dplyr::arrange(YSF) %>%  # not strictly necessary, but harmless
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  # Sanity check: confirm uniqueness
  n_total <- nrow(cohort_one_row)
  n_unique <- dplyr::n_distinct(cohort_one_row$pixel_ID)
  if (n_total != n_unique) {
    warning("Cohort YSF", target_ysf, ": still has duplicate pixel_IDs after deduping. Check upstream logic.")
  }
  
  # Binary: recovered by the last year of this observation window
  cohort_one_row <- cohort_one_row %>%
    dplyr::mutate(
      last_year_observed = as.integer(target_ysf),
      Recovered_by_last_year = dplyr::if_else(
        !is.na(yrs_to_recovery) & yrs_to_recovery <= target_ysf,
        1L, 0L
      )
    )
  
  # Confirm the flag exists and is binary
  stopifnot("Recovered_by_last_year" %in% names(cohort_one_row))
  stopifnot(all(cohort_one_row$Recovered_by_last_year %in% c(0L, 1L)))
  
  cohort_one_row
}

# Build each cohort as its own object
cohort_10 <- build_clean_cohort(10L)
cohort_15 <- build_clean_cohort(15L)
cohort_20 <- build_clean_cohort(20L)

# Optional: quick summaries so you can sanity-check in console
summarize_cohort <- function(df) {
  df %>%
    dplyr::group_by(YSF, sev_group) %>%
    dplyr::summarise(
      N_pixels = dplyr::n(),
      N_recovered_by_last = sum(Recovered_by_last_year == 1L, na.rm = TRUE),
      Prop_recovered_by_last = mean(Recovered_by_last_year == 1L, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(YSF, sev_group)
}

summary_10 <- summarize_cohort(cohort_10)
summary_15 <- summarize_cohort(cohort_15)
summary_20 <- summarize_cohort(cohort_20)

summary_10
summary_15
summary_20
# ===============================================================
#  Mapping Recovered vs Unrecovered Pixels ----------------------

stopifnot(exists("pixel_geom_sf"))
stopifnot("pixel_ID" %in% names(pixel_geom_sf))
stopifnot(exists("cohort_10"), exists("cohort_15"), exists("cohort_20"))

stopifnot(exists("huc12_bg"))
# Otherwise, define here from existing objects:
# huc12_bg   <- huc12_sf
# streams_bg <- streams_sf

rec_cols <- c(
  "Recovered"   = "palegreen3",
  "Unrecovered" = "red"
)

## Read in fire polygons:
fire_atlas_raw <- sf::st_read(
  "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/ArcGIS/BMWC_Fire_Data/fire_atlas/fire_atlas.shp",
  quiet = TRUE
)

# Assign the TRUE CRS (no reprojection yet)
fire_atlas_raw <- sf::st_set_crs(fire_atlas_raw, 5070)  # NAD83 / Conus Albers

# Now transform to match pixel geometry CRS (WGS84 lon/lat)
fire_atlas_bg <- sf::st_transform(fire_atlas_raw, sf::st_crs(pixel_geom_sf))

# Sanity check: bbox should now be ~(-113, 47)
sf::st_bbox(fire_atlas_bg)

fire_atlas_bg <- fire_atlas_bg %>%
  sf::st_make_valid() %>%
  sf::st_cast("MULTIPOLYGON", warn = FALSE)
fire_atlas_bg <- sf::st_crop(fire_atlas_bg, sf::st_bbox(huc12_bg))


## Mapping function ---------------------------
build_one_cohort_map <- function(cohort_df,
                                 pixel_geom_sf,
                                 huc12_bg,
                                 streams_bg = NULL,
                                 fire_atlas_bg = NULL,
                                 unrecovered_only = FALSE,
                                 add_streams = FALSE,
                                 point_size = 1.4,
                                 # ---- FIRE POLYGON SYMBOLOGY SWITCHES 
                                 fire_outline_color = "darkorange3",
                                 fire_outline_width = 0.2,
                                 fire_fill_color    = "orange1",      # NA = no fill
                                 fire_alpha         = 0.25) {     # 0–1
  
  stopifnot(all(c("YSF", "pixel_ID", "Recovered_by_last_year") %in% names(cohort_df)))
  stopifnot(inherits(pixel_geom_sf, "sf"))
  stopifnot(inherits(huc12_bg, "sf"))
  
  cohort_ysf <- unique(cohort_df$YSF)
  if (length(cohort_ysf) != 1L) stop("cohort_df must contain exactly one YSF value.")
  
  # Join recovery status to pixel geometry
  px_map <- pixel_geom_sf %>%
    dplyr::mutate(pixel_ID = as.character(pixel_ID)) %>%
    dplyr::left_join(
      cohort_df %>%
        dplyr::mutate(pixel_ID = as.character(pixel_ID)) %>%
        dplyr::select(pixel_ID, Recovered_by_last_year),
      by = "pixel_ID"
    ) %>%
    dplyr::mutate(
      recovery_status = dplyr::case_when(
        Recovered_by_last_year == 1L ~ "Recovered",
        Recovered_by_last_year == 0L ~ "Unrecovered",
        TRUE                         ~ NA_character_
      ),
      recovery_status = factor(recovery_status, levels = c("Recovered", "Unrecovered"))
    ) %>%
    dplyr::filter(!is.na(recovery_status))
  
  if (unrecovered_only) {
    px_map <- px_map %>% dplyr::filter(recovery_status == "Unrecovered")
  }
  
  n_total <- nrow(px_map)
  n_rec   <- sum(px_map$recovery_status == "Recovered", na.rm = TRUE)
  n_unrec <- sum(px_map$recovery_status == "Unrecovered", na.rm = TRUE)
  
  p <- ggplot() +
    geom_sf(
      data = huc12_bg,
      fill = "grey95",
      color = "grey70",
      linewidth = 0.2
    )
  
  # Fire boundaries (only for unrecovered-only maps)
  if (unrecovered_only && !is.null(fire_atlas_bg)) {
    p <- p +
      geom_sf(
        data = fire_atlas_bg,
        fill = fire_fill_color,
        color = fire_outline_color,
        linewidth = fire_outline_width,
        alpha = fire_alpha
      )
  }
  
  # Optional streams
  if (add_streams && !is.null(streams_bg)) {
    p <- p +
      geom_sf(
        data = streams_bg,
        color = "darkblue",
        linewidth = 0.5
      )
  }
  
  if (unrecovered_only) {
    p <- p +
      geom_sf(
        data = px_map,
        color = "red",
        size = point_size,
        alpha = 0.85
      ) +
      labs(
        title    = paste0("Unrecovered Pixels, ", cohort_ysf, " YSF"),
        subtitle = paste0("N = ", n_total), hjust = 0.5,
        caption  = "Fire polygons shown in orange"
      )
  } else {
    p <- p +
      geom_sf(
        data = px_map,
        aes(color = recovery_status),
        size = point_size,
        alpha = 0.8
      ) +
      scale_color_manual(
        values = c("Recovered" = "palegreen3", "Unrecovered" = "red"),
        name   = NULL
      ) +
      labs(
        title    = paste0("Recovery Status: ", cohort_ysf, " YSF"),
        subtitle = paste0("Recovered vs Unrecovered by ", cohort_ysf, " years since fire (N = ", n_total, ")"),
        caption  = "'Recovered' based on 1-year recovery rule"
      )
  }
  
  p +
    coord_sf(expand = FALSE) +
    theme_minimal(base_size = 14, base_family = "Times New Roman") +
    theme(
      plot.title      = element_text(face = "bold"),
      plot.subtitle   = element_text(size = 12),
      plot.caption    = element_text(size = 10),
      legend.position = "right",
      legend.title    = element_blank()
    )
}

# --- Print all maps --------------------------------
print_all_cohort_maps <- function(cohort_list,
                                  pixel_geom_sf,
                                  huc12_bg,
                                  fire_atlas_bg,
                                  streams_bg = NULL,
                                  add_streams = FALSE) {
  
  stopifnot(is.list(cohort_list))
  
  for (nm in names(cohort_list)) {
    
    cohort_df <- cohort_list[[nm]]
    cohort_ysf <- unique(cohort_df$YSF)
    
    message("\n==============================")
    message("Mapping cohort ≥ ", cohort_ysf, " YSF")
    message("==============================")
    
    # ---- Full recovery map
    p_all <- build_one_cohort_map(
      cohort_df       = cohort_df,
      pixel_geom_sf   = pixel_geom_sf,
      huc12_bg        = huc12_bg,
      streams_bg      = streams_bg,
      add_streams     = add_streams,
      unrecovered_only = FALSE
    )
    
    print(p_all)
    
    # ---- Unrecovered-only + fire boundaries
    p_unrec <- build_one_cohort_map(
      cohort_df        = cohort_df,
      pixel_geom_sf    = pixel_geom_sf,
      huc12_bg         = huc12_bg,
      fire_atlas_bg    = fire_atlas_bg,
      streams_bg       = streams_bg,
      add_streams      = add_streams,
      unrecovered_only = TRUE
    )
    
    print(p_unrec)
  }
}

cohort_list <- list(
  YSF10 = cohort_10,
  YSF15 = cohort_15,
  YSF20 = cohort_20
)

print_all_cohort_maps(
  cohort_list     = cohort_list,
  pixel_geom_sf   = pixel_geom_sf,
  huc12_bg        = huc12_bg,
  fire_atlas_bg   = fire_atlas_bg,
  streams_bg      = streams_bg,   # optional
  add_streams     = FALSE         # flip TRUE if you want streams
)
# ---------------------------------------------------------------
################################################################################
# ===============================================================
################################################################################
# ===============================================================
#  Unrecovered Pixels Colored by Categorical Fire Severity
# ===============================================================

stopifnot(exists("pixel_geom_sf"))
stopifnot(exists("huc12_bg"))
stopifnot(exists("fire_atlas_bg"))
stopifnot(exists("cohort_10"), exists("cohort_15"), exists("cohort_20"))

sev_colors <- c(
  "Unburned" = "grey60",
  "Low"      = "yellow",
  "Moderate" = "orange2",
  "High"     = "darkred"
)

## Mapping function:
build_one_cohort_unrec_sev_map <- function(cohort_df,
                                           pixel_geom_sf,
                                           huc12_bg,
                                           streams_bg = NULL,
                                           fire_atlas_bg = NULL,
                                           add_streams = FALSE,
                                           point_size = 2,
                                           fire_outline_color = "orange4",
                                           fire_outline_width = 0.2,
                                           fire_fill_color    = "orange",
                                           fire_alpha         = 0.25) {
  
  stopifnot(all(c("YSF","pixel_ID","Recovered_by_last_year","sev_group") %in% names(cohort_df)))
  stopifnot(inherits(pixel_geom_sf, "sf"))
  stopifnot(inherits(huc12_bg, "sf"))
  
  cohort_ysf <- unique(cohort_df$YSF)
  if (length(cohort_ysf) != 1L) stop("cohort_df must contain exactly one YSF value.")
  
  px_map <- pixel_geom_sf %>%
    dplyr::mutate(pixel_ID = as.character(pixel_ID)) %>%
    dplyr::left_join(
      cohort_df %>%
        dplyr::mutate(pixel_ID = as.character(pixel_ID)) %>%
        dplyr::select(pixel_ID, Recovered_by_last_year, sev_group),
      by = "pixel_ID"
    ) %>%
    dplyr::filter(Recovered_by_last_year == 0L) %>%
    dplyr::mutate(
      sev_group = factor(sev_group, levels = c("Unburned","Low","Moderate","High"))
    ) %>%
    dplyr::filter(!is.na(sev_group))
  
  # ---- Counts + percents for subtitle
  n_total <- nrow(px_map)
  n_mod   <- sum(px_map$sev_group == "Moderate", na.rm = TRUE)
  n_high  <- sum(px_map$sev_group == "High",     na.rm = TRUE)
  
  pct_mod  <- if (n_total > 0) 100 * n_mod  / n_total else NA_real_
  pct_high <- if (n_total > 0) 100 * n_high / n_total else NA_real_
  
  # ---- Base map
  p <- ggplot() +
    geom_sf(
      data = huc12_bg,
      fill = "grey95",
      color = "grey70",
      linewidth = 0.2
    )
  
  # ---- Fire polygons
  if (!is.null(fire_atlas_bg)) {
    p <- p +
      geom_sf(
        data = fire_atlas_bg,
        fill = fire_fill_color,
        color = fire_outline_color,
        linewidth = fire_outline_width,
        alpha = fire_alpha
      )
  }
  
  # ---- Optional streams
  if (add_streams && !is.null(streams_bg)) {
    p <- p +
      geom_sf(
        data = streams_bg,
        color = "darkblue",
        linewidth = 0.5
      )
  }
  
  # ---- Pixels (plot order: Unburned -> Low -> Moderate -> High)
  px_map <- px_map %>% dplyr::arrange(sev_group)
  
  p <- p +
    geom_sf(
      data = px_map,
      aes(color = sev_group),
      size = point_size,
      alpha = 0.85
    ) +
    scale_color_manual(
      values = sev_colors,
      name = "Fire Severity",
      drop = FALSE
    ) +
    labs(
      title = paste0("Unrecovered Pixels by Fire Severity, ", cohort_ysf, " YSF"),
      subtitle = paste0(
        "N = ", format(n_total, big.mark = ","),
        " | Moderate = ", format(n_mod, big.mark = ","), " (", sprintf("%.1f%%", pct_mod), ")",
        " | High = ",     format(n_high, big.mark = ","), " (", sprintf("%.1f%%", pct_high), ")"
      ),
      caption = "Fire polygons shown in orange"
    ) +
    coord_sf(expand = FALSE) +
    theme_minimal(base_size = 14, base_family = "Times New Roman") +
    theme(
      plot.title      = element_text(face = "bold"),
      plot.subtitle   = element_text(size = 12),
      plot.caption    = element_text(size = 10),
      legend.position = "right"
    )
  
  p
}

print_all_unrec_sev_maps <- function(cohort_list,
                                     pixel_geom_sf,
                                     huc12_bg,
                                     fire_atlas_bg,
                                     streams_bg = NULL,
                                     add_streams = FALSE) {
  
  stopifnot(is.list(cohort_list))
  
  for (nm in names(cohort_list)) {
    
    cohort_df <- cohort_list[[nm]]
    cohort_ysf <- unique(cohort_df$YSF)
    
    message("\n==============================")
    message("Mapping UNRECOVERED pixels by severity: ", cohort_ysf, " YSF")
    message("==============================")
    
    p <- build_one_cohort_unrec_sev_map(
      cohort_df     = cohort_df,
      pixel_geom_sf = pixel_geom_sf,
      huc12_bg      = huc12_bg,
      fire_atlas_bg = fire_atlas_bg,
      streams_bg    = streams_bg,
      add_streams   = add_streams
    )
    
    print(p)
  }
}

cohort_list <- list(
  YSF10 = cohort_10,
  YSF15 = cohort_15,
  YSF20 = cohort_20
)

print_all_unrec_sev_maps(
  cohort_list    = cohort_list,
  pixel_geom_sf  = pixel_geom_sf,
  huc12_bg       = huc12_bg,
  fire_atlas_bg  = fire_atlas_bg,
  streams_bg     = streams_bg,
  add_streams    = FALSE
)

################################################################################
################################################################################
# ===============================================================
#  Explore whether UNRECOVERED HIGH-severity pixels occur in larger burns
#  Output:
#   (A) Summary tables comparing fire_size_ha (unrecovered only; includes High-only comparison)
#   (B) Maps of UNRECOVERED HIGH-severity pixels colored by fire size (continuous, log scale)
# ===============================================================

stopifnot(exists("pixel_recovery_by_YSF"))
stopifnot(all(c("YSF","pixel_ID","sev_group","fire_name","ref_year","yrs_to_recovery","fire_size_ha") %in%
                names(pixel_recovery_by_YSF)))

stopifnot(exists("pixel_geom_sf"))
stopifnot("pixel_ID" %in% names(pixel_geom_sf))

stopifnot(exists("huc12_bg"))
stopifnot(inherits(huc12_bg, "sf"))

# -----------------------------
# 1) Rebuild cohorts INCLUDING fire_size_ha
# -----------------------------
cohorts_to_build <- c(10L, 15L, 20L)

build_clean_cohort <- function(target_ysf) {
  
  cohort_raw <- pixel_recovery_by_YSF %>%
    dplyr::filter(YSF == target_ysf) %>%
    dplyr::select(YSF, pixel_ID, sev_group, fire_name, ref_year, yrs_to_recovery, fire_size_ha)
  
  cohort_one_row <- cohort_raw %>%
    dplyr::group_by(pixel_ID) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  # Sanity check: confirm uniqueness
  n_total  <- nrow(cohort_one_row)
  n_unique <- dplyr::n_distinct(cohort_one_row$pixel_ID)
  if (n_total != n_unique) {
    warning("Cohort YSF ", target_ysf, ": duplicate pixel_IDs after deduping. Check upstream logic.")
  }
  
  cohort_one_row <- cohort_one_row %>%
    dplyr::mutate(
      last_year_observed = as.integer(target_ysf),
      Recovered_by_last_year = dplyr::if_else(
        !is.na(yrs_to_recovery) & yrs_to_recovery <= target_ysf,
        1L, 0L
      )
    )
  
  stopifnot("Recovered_by_last_year" %in% names(cohort_one_row))
  stopifnot(all(cohort_one_row$Recovered_by_last_year %in% c(0L, 1L)))
  
  cohort_one_row
}

cohort_10 <- build_clean_cohort(10L)
cohort_15 <- build_clean_cohort(15L)
cohort_20 <- build_clean_cohort(20L)

cohort_list <- list(YSF10 = cohort_10, YSF15 = cohort_15, YSF20 = cohort_20)

# -----------------------------
# 2) Summary: UNRECOVERED pixels only, fire_size_ha by severity
# -----------------------------
unrec_fire_size_summary <- dplyr::bind_rows(cohort_list, .id = "cohort_name") %>%
  dplyr::filter(Recovered_by_last_year == 0L) %>%
  dplyr::mutate(
    sev_group = factor(sev_group, levels = c("Unburned","Low","Moderate","High"))
  ) %>%
  dplyr::group_by(YSF, sev_group) %>%
  dplyr::summarise(
    N = dplyr::n(),
    fire_size_ha_median = stats::median(fire_size_ha, na.rm = TRUE),
    fire_size_ha_mean   = mean(fire_size_ha, na.rm = TRUE),
    fire_size_ha_p25    = stats::quantile(fire_size_ha, 0.25, na.rm = TRUE),
    fire_size_ha_p75    = stats::quantile(fire_size_ha, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(YSF, sev_group)

unrec_fire_size_summary

# Optional: high vs other (unrecovered only)
unrec_high_vs_other <- dplyr::bind_rows(cohort_list, .id = "cohort_name") %>%
  dplyr::filter(Recovered_by_last_year == 0L) %>%
  dplyr::mutate(
    sev_group = factor(sev_group, levels = c("Unburned","Low","Moderate","High")),
    high_vs_other = dplyr::if_else(sev_group == "High", "High", "Not High")
  ) %>%
  dplyr::group_by(YSF, high_vs_other) %>%
  dplyr::summarise(
    N = dplyr::n(),
    fire_size_ha_median = stats::median(fire_size_ha, na.rm = TRUE),
    fire_size_ha_mean   = mean(fire_size_ha, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(YSF, dplyr::desc(high_vs_other))

unrec_high_vs_other

# Optional: boxplot (unrecovered only; log y)
ggplot(
  dplyr::bind_rows(cohort_list, .id = "cohort_name") %>%
    dplyr::filter(Recovered_by_last_year == 0L) %>%
    dplyr::mutate(sev_group = factor(sev_group, levels = c("Unburned","Low","Moderate","High"))),
  aes(x = sev_group, y = fire_size_ha)
) +
  geom_boxplot(outlier.alpha = 0.25) +
  scale_y_continuous(trans = "log10") +
  facet_wrap(~ YSF, nrow = 1) +
  labs(
    title = "Unrecovered pixels: fire size (ha) by severity (log10 scale)",
    x = "Fire severity",
    y = "Fire size (ha, log10)"
  ) +
  theme_minimal(base_size = 13, base_family = "Times New Roman")

# -----------------------------
# 3) Map: UNRECOVERED HIGH-severity pixels colored by fire_size_ha
# -----------------------------
build_one_cohort_unrec_high_firesize_map <- function(cohort_df,
                                                     pixel_geom_sf,
                                                     huc12_bg,
                                                     streams_bg = NULL,
                                                     fire_atlas_bg = NULL,
                                                     add_streams = FALSE,
                                                     point_size = 2,
                                                     fire_outline_color = "orange4",
                                                     fire_outline_width = 0.2,
                                                     fire_fill_color    = "orange",
                                                     fire_alpha         = 0.20,
                                                     log_color = TRUE) {
  
  stopifnot(all(c("YSF","pixel_ID","Recovered_by_last_year","fire_size_ha","sev_group") %in% names(cohort_df)))
  stopifnot(inherits(pixel_geom_sf, "sf"))
  stopifnot(inherits(huc12_bg, "sf"))
  
  cohort_ysf <- unique(cohort_df$YSF)
  if (length(cohort_ysf) != 1L) stop("cohort_df must contain exactly one YSF value.")
  
  # Join fire size to geometry, keep UNRECOVERED + HIGH severity only
  px_map <- pixel_geom_sf %>%
    dplyr::mutate(pixel_ID = as.character(pixel_ID)) %>%
    dplyr::left_join(
      cohort_df %>%
        dplyr::mutate(pixel_ID = as.character(pixel_ID)) %>%
        dplyr::select(pixel_ID, Recovered_by_last_year, fire_size_ha, sev_group),
      by = "pixel_ID"
    ) %>%
    dplyr::filter(
      Recovered_by_last_year == 0L,
      sev_group == "High"
    ) %>%
    dplyr::filter(!is.na(fire_size_ha))
  
  n_total <- nrow(px_map)
  
  p <- ggplot() +
    geom_sf(
      data = huc12_bg,
      fill = "grey95",
      color = "grey70",
      linewidth = 0.2
    )
  
  # Optional fire polygons for context
  if (!is.null(fire_atlas_bg)) {
    p <- p +
      geom_sf(
        data = fire_atlas_bg,
        fill = fire_fill_color,
        color = fire_outline_color,
        linewidth = fire_outline_width,
        alpha = fire_alpha
      )
  }
  
  if (add_streams && !is.null(streams_bg)) {
    p <- p +
      geom_sf(
        data = streams_bg,
        color = "darkblue",
        linewidth = 0.5
      )
  }
  
  p <- p +
    geom_sf(
      data = px_map,
      aes(color = fire_size_ha),
      size = point_size,
      alpha = 0.85
    ) +
    labs(
      title = paste0("Unrecovered HIGH-severity pixels colored by fire size, ", cohort_ysf, " YSF"),
      subtitle = paste0(
        "N (unrecovered high severity) = ", format(n_total, big.mark = ","),
        if (log_color) " | color scale is log10(ha)" else ""
      ),
      color = "Fire size (ha)",
      caption = "Fire polygons shown in orange (for context)"
    ) +
    coord_sf(expand = FALSE) +
    theme_minimal(base_size = 14, base_family = "Times New Roman") +
    theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(size = 12),
      plot.caption  = element_text(size = 10),
      legend.position = "right"
    )
  
  # Prefer viridis if available (best for continuous), otherwise fall back to default gradient
  if (requireNamespace("viridis", quietly = TRUE)) {
    if (log_color) {
      p <- p + viridis::scale_color_viridis(trans = "log10", option = "C")
    } else {
      p <- p + viridis::scale_color_viridis(option = "C")
    }
  } else {
    if (log_color) {
      p <- p + scale_color_gradient(trans = "log10")
    } else {
      p <- p + scale_color_gradient()
    }
  }
  
  p
}

print_all_unrec_high_firesize_maps <- function(cohort_list,
                                               pixel_geom_sf,
                                               huc12_bg,
                                               fire_atlas_bg = NULL,
                                               streams_bg = NULL,
                                               add_streams = FALSE,
                                               log_color = TRUE) {
  
  stopifnot(is.list(cohort_list))
  
  for (nm in names(cohort_list)) {
    
    cohort_df <- cohort_list[[nm]]
    cohort_ysf <- unique(cohort_df$YSF)
    
    message("\n==============================")
    message("Mapping UNRECOVERED HIGH pixels by fire size: ", cohort_ysf, " YSF")
    message("==============================")
    
    p <- build_one_cohort_unrec_high_firesize_map(
      cohort_df     = cohort_df,
      pixel_geom_sf = pixel_geom_sf,
      huc12_bg      = huc12_bg,
      fire_atlas_bg = fire_atlas_bg,
      streams_bg    = streams_bg,
      add_streams   = add_streams,
      log_color     = log_color
    )
    
    print(p)
  }
}

# Run maps (log scale strongly recommended; fire sizes tend to be highly skewed)
print_all_unrec_high_firesize_maps(
  cohort_list    = cohort_list,
  pixel_geom_sf  = pixel_geom_sf,
  huc12_bg       = huc12_bg,
  fire_atlas_bg  = fire_atlas_bg,
  streams_bg     = streams_bg,
  add_streams    = FALSE,
  log_color      = TRUE
)
# ===============================================================
# 4) Histograms: Unrecovered Moderate + High severity pixels by fire size
# ===============================================================
hist_dat <- dplyr::bind_rows(cohort_list, .id = "cohort_name") %>%
  dplyr::filter(
    YSF %in% c(10, 15),
    Recovered_by_last_year == 0L,
    sev_group %in% c("Moderate", "High"),
    !is.na(fire_size_ha),
    fire_size_ha > 0
  ) %>%
  dplyr::mutate(
    sev_group = factor(sev_group, levels = c("Moderate", "High")),
    YSF = factor(YSF, levels = c(10, 15))
  ) %>%
  droplevels()

# Quick checks
hist_dat %>% dplyr::count(YSF, sev_group)
table(hist_dat$YSF, useNA = "ifany")

# ---- Version 1: raw fire size on x-axis --------
ggplot(hist_dat, aes(x = fire_size_ha, fill = sev_group)) +
  geom_histogram(
    bins = 30,
    color = "black",
    linewidth = 0.2
  ) +
  facet_grid(sev_group ~ YSF, scales = "free_y", drop = TRUE) +
  scale_fill_manual(values = c(
    "Moderate" = "orange",
    "High"     = "darkred"
  )) +
  labs(
    title = "Unrecovered Moderate and High-severity Pixels by Fire Size",
    #subtitle = "X-axis = fire size (ha)",
    x = "Fire size (ha)",
    y = "Number of pixels"
  ) +
  theme_minimal(base_size = 13, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

# ---- Version 2: recommended log10 x-axis for skewed fire sizes ---------
ggplot(hist_dat, aes(x = fire_size_ha, fill = sev_group)) +
  geom_histogram(
    bins = 30,
    color = "black",
    linewidth = 0.2
  ) +
  scale_x_log10() +
  facet_grid(sev_group ~ YSF, scales = "free_y", drop = TRUE) +
  scale_fill_manual(values = c(
    "Moderate" = "orange",
    "High"     = "darkred"
  )) +
  labs(
    title = "Unrecovered Moderate and High-severity Pixels by Fire Size",
    subtitle = "**Log10 x-axis",
    x = "Fire size (ha, log10 scale)",
    y = "Number of pixels"
  ) +
  theme_minimal(base_size = 13, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

## =============================================================
# BACKDROP HISTOGRAM DATA
# all pixels in each cohort, regardless of recovery/severity
# ===============================================================
hist_bg <- dplyr::bind_rows(cohort_list, .id = "cohort_name") %>%
  dplyr::filter(
    YSF %in% c(10, 15),
    !is.na(fire_size_ha),
    fire_size_ha > 0
  ) %>%
  dplyr::mutate(
    YSF = factor(YSF, levels = c(10, 15))
  ) %>%
  droplevels()

# ===============================================================
# FOREGROUND HISTOGRAM DATA
# unrecovered Moderate + High pixels only
# ===============================================================
hist_fg <- dplyr::bind_rows(cohort_list, .id = "cohort_name") %>%
  dplyr::filter(
    YSF %in% c(10, 15),
    Recovered_by_last_year == 0L,
    sev_group %in% c("Moderate", "High"),
    !is.na(fire_size_ha),
    fire_size_ha > 0
  ) %>%
  dplyr::mutate(
    sev_group = factor(sev_group, levels = c("Moderate", "High")),
    YSF = factor(YSF, levels = c(10, 15))
  ) %>%
  droplevels()

# === Regular scale histogram: =====
ggplot() +
  geom_histogram(
    data = hist_bg,
    aes(x = fire_size_ha),
    bins = 30,
    fill = "grey80",
    color = "black",
    linewidth = 0.2,
    alpha = 0.7
  ) +
  geom_histogram(
    data = hist_fg,
    aes(x = fire_size_ha, fill = sev_group),
    bins = 30,
    color = "black",
    linewidth = 0.2,
    alpha = 0.75,
    position = "identity"
  ) +
  facet_wrap(~ YSF, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c(
    "Moderate" = "orange",
    "High"     = "darkred"
  )) +
  labs(
    title = "Unrecovered Moderate and High-severity Pixels",
    subtitle = "Gray bars show all pixels in each cohort; colored bars show unrecovered Moderate and High pixels",
    x = "Fire size (ha)",
    y = "Number of pixels",
    fill = "Severity"
  ) +
  theme_minimal(base_size = 13, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.minor = element_blank()
  )

# === Log10 scale histogram: ====
ggplot() +
  geom_histogram(
    data = hist_bg,
    aes(x = fire_size_ha),
    bins = 30,
    fill = "grey80",
    color = "black",
    linewidth = 0.2,
    alpha = 0.7
  ) +
  geom_histogram(
    data = hist_fg,
    aes(x = fire_size_ha, fill = sev_group),
    bins = 30,
    color = "black",
    linewidth = 0.2,
    alpha = 0.75,
    position = "identity"
  ) +
  scale_x_log10() +
  facet_wrap(~ YSF, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c(
    "Moderate" = "orange",
    "High"     = "darkred"
  )) +
  labs(
    title = "Unrecovered Moderate and High-severity Pixels",
    subtitle = "Gray bars show all pixels; colored bars show unrecovered Moderate and High pixels",
    x = "Fire size (ha, log10 scale)",
    y = "Number of pixels",
    fill = "Severity"
  ) +
  theme_minimal(base_size = 13, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.minor = element_blank()
  )
# === What % of unrecovered pixels in largest fires? =====
fire_thresh <- 10000  # hectares
cohort_list <- list(
  YSF10 = cohort_10,
  YSF15 = cohort_15
)

unrec_mh <- dplyr::bind_rows(cohort_list) %>%
  dplyr::filter(
    Recovered_by_last_year == 0L,
    sev_group %in% c("Moderate", "High"),
    !is.na(fire_size_ha),
    fire_size_ha > 0
  )
pct_large_overall <- mean(unrec_mh$fire_size_ha >= fire_thresh) * 100
pct_large_overall

# By severity:
pct_by_sev <- unrec_mh %>%
  dplyr::group_by(sev_group) %>%
  dplyr::summarise(
    N = dplyr::n(),
    N_large = sum(fire_size_ha >= fire_thresh),
    pct_large = 100 * N_large / N,
    .groups = "drop"
  )
pct_by_sev

# By YSF:
pct_by_ysf <- unrec_mh %>%
  dplyr::group_by(YSF, sev_group) %>%
  dplyr::summarise(
    N = dplyr::n(),
    N_large = sum(fire_size_ha >= fire_thresh),
    pct_large = 100 * N_large / N,
    .groups = "drop"
  )
pct_by_ysf

# === Print out info for largest fires ======
large_fires <- dplyr::bind_rows(cohort_list) %>%
  dplyr::filter(
    fire_name != "None",
    !is.na(fire_size_ha),
    fire_size_ha >= 10000
  ) %>%
  dplyr::group_by(fire_name) %>%
  dplyr::summarise(
    ref_year = unique(ref_year)[1],
    fire_size_ha = max(fire_size_ha, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(ref_year)

large_fires
# ============================================
# Map unrecovered Moderate + High pixels in large fires only
# Large fires defined as fire_size_ha >= 10,000 ha
# ===============================================================

stopifnot(exists("pixel_geom_sf"))
stopifnot(exists("huc12_bg"))
stopifnot(exists("fire_atlas_bg"))
stopifnot(exists("cohort_10"), exists("cohort_15"))

sev_colors_large <- c(
  "Moderate" = "orange",
  "High"     = "darkred"
)

build_one_largefire_unrec_map <- function(cohort_df,
                                          pixel_geom_sf,
                                          huc12_bg,
                                          streams_bg = NULL,
                                          fire_atlas_bg = NULL,
                                          add_streams = FALSE,
                                          point_size = 2,
                                          fire_thresh = 10000,
                                          fire_outline_color = "orange4",
                                          fire_outline_width = 0.2,
                                          fire_fill_color    = "orange",
                                          fire_alpha         = 0.20) {
  
  stopifnot(all(c("YSF", "pixel_ID", "Recovered_by_last_year", "sev_group",
                  "fire_size_ha", "fire_name") %in% names(cohort_df)))
  stopifnot(inherits(pixel_geom_sf, "sf"))
  stopifnot(inherits(huc12_bg, "sf"))
  
  cohort_ysf <- unique(cohort_df$YSF)
  if (length(cohort_ysf) != 1L) stop("cohort_df must contain exactly one YSF value.")
  
  # Join filtered cohort data to geometry
  px_map <- pixel_geom_sf %>%
    dplyr::mutate(pixel_ID = as.character(pixel_ID)) %>%
    dplyr::left_join(
      cohort_df %>%
        dplyr::mutate(pixel_ID = as.character(pixel_ID)) %>%
        dplyr::filter(
          Recovered_by_last_year == 0L,
          sev_group %in% c("Moderate", "High"),
          !is.na(fire_size_ha),
          fire_size_ha >= fire_thresh,
          !is.na(fire_name),
          fire_name != "None"
        ) %>%
        dplyr::select(pixel_ID, sev_group, fire_name, fire_size_ha),
      by = "pixel_ID"
    ) %>%
    dplyr::filter(!is.na(sev_group)) %>%
    dplyr::mutate(
      sev_group = factor(sev_group, levels = c("Moderate", "High"))
    )
  
  # Counts for subtitle
  n_total <- nrow(px_map)
  n_mod   <- sum(px_map$sev_group == "Moderate", na.rm = TRUE)
  n_high  <- sum(px_map$sev_group == "High", na.rm = TRUE)
  n_fires <- px_map %>%
    sf::st_drop_geometry() %>%
    dplyr::distinct(fire_name) %>%
    nrow()
  
  # Plot order so high plots on top
  px_map <- px_map %>% dplyr::arrange(sev_group)
  
  p <- ggplot() +
    geom_sf(
      data = huc12_bg,
      fill = "grey95",
      color = "grey70",
      linewidth = 0.2
    )
  
  # Fire polygons for context
  if (!is.null(fire_atlas_bg)) {
    p <- p +
      geom_sf(
        data = fire_atlas_bg,
        fill = fire_fill_color,
        color = fire_outline_color,
        linewidth = fire_outline_width,
        alpha = fire_alpha
      )
  }
  
  # Optional streams
  if (add_streams && !is.null(streams_bg)) {
    p <- p +
      geom_sf(
        data = streams_bg,
        color = "darkblue",
        linewidth = 0.5
      )
  }
  
  p +
    geom_sf(
      data = px_map,
      aes(color = sev_group),
      size = point_size,
      alpha = 0.9
    ) +
    scale_color_manual(
      values = sev_colors_large,
      name = "Fire Severity",
      drop = FALSE
    ) +
    labs(
      title = paste0("Unrecovered Pixels in Large Fires, ", cohort_ysf, " YSF"),
      subtitle = paste0(
        "Large fires defined as \u2265 ", format(fire_thresh, big.mark = ","), " ha",
        " | N pixels = ", format(n_total, big.mark = ","),
        " | Moderate = ", format(n_mod, big.mark = ","),
        " | High = ", format(n_high, big.mark = ","),
        " | N fires = ", n_fires
      ),
      caption = "Fire polygons shown in orange"
    ) +
    coord_sf(expand = FALSE) +
    theme_minimal(base_size = 14, base_family = "Times New Roman") +
    theme(
      plot.title      = element_text(face = "bold"),
      plot.subtitle   = element_text(size = 12),
      plot.caption    = element_text(size = 10),
      legend.position = "right"
    )
}

print_largefire_unrec_maps <- function(cohort_list,
                                       pixel_geom_sf,
                                       huc12_bg,
                                       fire_atlas_bg = NULL,
                                       streams_bg = NULL,
                                       add_streams = FALSE,
                                       fire_thresh = 10000) {
  
  stopifnot(is.list(cohort_list))
  
  for (nm in names(cohort_list)) {
    cohort_df  <- cohort_list[[nm]]
    cohort_ysf <- unique(cohort_df$YSF)
    
    message("\n==============================")
    message("Mapping unrecovered Moderate/High pixels in large fires: ", cohort_ysf, " YSF")
    message("==============================")
    
    p <- build_one_largefire_unrec_map(
      cohort_df     = cohort_df,
      pixel_geom_sf = pixel_geom_sf,
      huc12_bg      = huc12_bg,
      streams_bg    = streams_bg,
      fire_atlas_bg = fire_atlas_bg,
      add_streams   = add_streams,
      fire_thresh   = fire_thresh
    )
    
    print(p)
  }
}

# Only 10 and 15 YSF
cohort_list_large <- list(
  YSF10 = cohort_10,
  YSF15 = cohort_15
)

print_largefire_unrec_maps(
  cohort_list    = cohort_list_large,
  pixel_geom_sf  = pixel_geom_sf,
  huc12_bg       = huc12_bg,
  fire_atlas_bg  = fire_atlas_bg,
  streams_bg     = streams_bg,
  add_streams    = FALSE,
  fire_thresh    = 10000
)

# ===============================================================
# Map unrecovered Moderate + High pixels, one map per LARGE fire
# ===============================================================
stopifnot(exists("pixel_geom_sf"))
stopifnot(exists("huc12_bg"))
stopifnot(exists("fire_atlas_bg"))
stopifnot(exists("cohort_10"), exists("cohort_15"), exists("cohort_20"))

sev_colors_large <- c(
  "Moderate" = "orange",
  "High"     = "darkred"
)

# transform CRS
safe_transform_to <- function(x, target) {
  if (is.na(sf::st_crs(x))) return(x)
  if (sf::st_crs(x) == sf::st_crs(target)) return(x)
  sf::st_transform(x, sf::st_crs(target))
}

# Function: crop layer to bbox
safe_crop_sf <- function(x, bbox_obj) {
  tryCatch(
    suppressWarnings(sf::st_crop(x, sf::st_bbox(bbox_obj))),
    error = function(e) x
  )
}
# ---------------------------------------------------------------
# FUNCTION: get fire perimeter from fire atlas
# Uses FireNam + Fire_Year from fire_atlas_bg
# ---------------------------------------------------------------
get_fire_perimeter <- function(fire_atlas_bg, fire_name_i, ref_year_i) {
  
  out <- fire_atlas_bg %>%
    dplyr::filter(
      FireNam == fire_name_i,
      Fire_Year == ref_year_i
    )
  
  if (nrow(out) == 0) {
    out <- fire_atlas_bg %>%
      dplyr::filter(FireNam == fire_name_i)
  }
  
  if (nrow(out) == 0) return(NULL)
  
  out %>%
    dplyr::summarise(.groups = "drop")
}

# ---------------------------------------------------------------
# FUNCTION: Build one map for one large fire
# ---------------------------------------------------------------
build_one_largefire_perfire_map <- function(cohort_df,
                                            fire_name_i,
                                            ref_year_i,
                                            pixel_geom_sf,
                                            huc12_bg,
                                            streams_bg = NULL,
                                            fire_atlas_bg = NULL,
                                            add_streams = TRUE,
                                            point_size = 3,  # pixel size
                                            fire_thresh = 10000,
                                            fire_outline_color = "orange4",
                                            fire_outline_width = 0.35,
                                            fire_fill_color    = "orange",
                                            fire_alpha         = 0.20,
                                            bbox_buffer        = 1500) {
  
  stopifnot(all(c("YSF", "pixel_ID", "Recovered_by_last_year", "sev_group",
                  "fire_size_ha", "fire_name", "ref_year") %in% names(cohort_df)))
  stopifnot(inherits(pixel_geom_sf, "sf"))
  stopifnot(inherits(huc12_bg, "sf"))
  
  cohort_ysf <- unique(cohort_df$YSF)
  if (length(cohort_ysf) != 1L) stop("cohort_df must contain exactly one YSF value.")
  
  # align CRS
  huc12_bg <- safe_transform_to(huc12_bg, pixel_geom_sf)
  if (!is.null(streams_bg)) streams_bg <- safe_transform_to(streams_bg, pixel_geom_sf)
  if (!is.null(fire_atlas_bg)) fire_atlas_bg <- safe_transform_to(fire_atlas_bg, pixel_geom_sf)
  
  # unrecovered Moderate/High pixels for this specific large fire
  px_sub <- cohort_df %>%
    dplyr::mutate(pixel_ID = as.character(pixel_ID)) %>%
    dplyr::filter(
      Recovered_by_last_year == 0L,
      sev_group %in% c("Moderate", "High"),
      !is.na(fire_size_ha),
      fire_size_ha >= fire_thresh,
      !is.na(fire_name),
      fire_name != "None",
      !is.na(ref_year),
      fire_name == fire_name_i,
      ref_year == ref_year_i
    ) %>%
    dplyr::select(pixel_ID, sev_group, fire_name, ref_year, fire_size_ha)
  
  if (nrow(px_sub) == 0) return(NULL)
  
  # join to geometry
  px_map <- pixel_geom_sf %>%
    dplyr::mutate(pixel_ID = as.character(pixel_ID)) %>%
    dplyr::inner_join(px_sub, by = "pixel_ID") %>%
    dplyr::mutate(
      sev_group = factor(sev_group, levels = c("Moderate", "High"))
    ) %>%
    dplyr::arrange(sev_group)
  
  # get fire polygon
  fire_poly <- NULL
  if (!is.null(fire_atlas_bg)) {
    fire_poly <- get_fire_perimeter(
      fire_atlas_bg = fire_atlas_bg,
      fire_name_i   = fire_name_i,
      ref_year_i    = ref_year_i
    )
  }
  
  # if no matching polygon found, fall back to pixel extent
  if (is.null(fire_poly) || nrow(fire_poly) == 0) {
    bbox_target <- sf::st_buffer(
      sf::st_as_sfc(sf::st_bbox(px_map)),
      dist = bbox_buffer
    )
    fire_crop <- NULL
  } else {
    # keep full perimeter as much as possible within study area
    study_area <- huc12_bg %>%
      sf::st_union() %>%
      sf::st_as_sf()
    
    fire_poly_in_study <- tryCatch(
      suppressWarnings(sf::st_intersection(fire_poly, study_area)),
      error = function(e) NULL
    )
    
    if (is.null(fire_poly_in_study) || nrow(fire_poly_in_study) == 0) {
      fire_poly_in_study <- fire_poly
    }
    
    bbox_target <- tryCatch(
      sf::st_buffer(fire_poly_in_study, dist = bbox_buffer),
      error = function(e) fire_poly_in_study
    )
    
    fire_crop <- safe_crop_sf(fire_poly_in_study, bbox_target)
  }
  
  huc12_crop <- safe_crop_sf(huc12_bg, bbox_target)
  
  streams_crop <- NULL
  if (add_streams && !is.null(streams_bg)) {
    streams_crop <- safe_crop_sf(streams_bg, bbox_target)
  }
  
  p <- ggplot() +
    geom_sf(
      data = huc12_crop,
      fill = "grey95",
      color = "grey70",
      linewidth = 0.2
    )
  
  if (!is.null(fire_crop)) {
    p <- p +
      geom_sf(
        data = fire_crop,
        fill = fire_fill_color,
        color = fire_outline_color,
        linewidth = fire_outline_width,
        alpha = fire_alpha
      )
  }
  
  if (add_streams && !is.null(streams_crop)) {
    p <- p +
      geom_sf(
        data = streams_crop,
        color = "blue4",
        linewidth = 0.2
      )
  }
  
  p +
    geom_sf(
      data = px_map,
      aes(color = sev_group),
      size = point_size,
      alpha = 1  # transparency
    ) +
    scale_color_manual(
      values = sev_colors_large,
      name = "Fire Severity",
      drop = FALSE
    ) +
    labs(
      title = paste0("Unrecovered Pixels in Large Fires, ", cohort_ysf, " YSF"),
      subtitle = paste0(fire_name_i, " | Ignition year: ", ref_year_i),
      caption = "Pixels shown are unrecovered Moderate- and High-severity pixels only"
    ) +
    coord_sf(expand = FALSE) +
    theme_minimal(base_size = 14, base_family = "Times New Roman") +
    theme(
      plot.title      = element_text(face = "bold"),
      plot.subtitle   = element_text(size = 12),
      plot.caption    = element_text(size = 10),
      legend.position = "right",
      panel.grid      = element_blank()
    )
}

# ---------------------------------------------------------------
# FUNCTION: Print one map per LARGE fire
# ---------------------------------------------------------------
print_largefire_perfire_maps <- function(cohort_df,
                                         pixel_geom_sf,
                                         huc12_bg,
                                         fire_atlas_bg = NULL,
                                         streams_bg = NULL,
                                         add_streams = FALSE,
                                         fire_thresh = 10000) {
  
  stopifnot(is.data.frame(cohort_df))
  
  fires_to_map <- cohort_df %>%
    dplyr::filter(
      Recovered_by_last_year == 0L,
      sev_group %in% c("Moderate", "High"),
      !is.na(fire_size_ha),
      fire_size_ha >= fire_thresh,
      !is.na(fire_name),
      fire_name != "None",
      !is.na(ref_year)
    ) %>%
    dplyr::distinct(fire_name, ref_year) %>%
    dplyr::arrange(ref_year, fire_name)
  
  if (nrow(fires_to_map) == 0) {
    message("No unrecovered Moderate/High pixels found in large fires for this cohort.")
    return(invisible(NULL))
  }
  
  for (i in seq_len(nrow(fires_to_map))) {
    fire_name_i <- fires_to_map$fire_name[i]
    ref_year_i  <- fires_to_map$ref_year[i]
    
    message("Printing map: ", fire_name_i, " | ", ref_year_i,
            " | YSF ", unique(cohort_df$YSF))
    
    p <- build_one_largefire_perfire_map(
      cohort_df      = cohort_df,
      fire_name_i    = fire_name_i,
      ref_year_i     = ref_year_i,
      pixel_geom_sf  = pixel_geom_sf,
      huc12_bg       = huc12_bg,
      streams_bg     = streams_bg,
      fire_atlas_bg  = fire_atlas_bg,
      add_streams    = add_streams,
      fire_thresh    = fire_thresh
    )
    
    if (!is.null(p)) print(p)
  }
  
  invisible(NULL)
}

# ---------------------------------------------------------------
# Run for 10, 15, 20 YSF
# ---------------------------------------------------------------
cohort_list_large_perfire <- list(
  YSF10 = cohort_10,
  YSF15 = cohort_15,
  YSF20 = cohort_20
)

for (nm in names(cohort_list_large_perfire)) {
  cohort_df <- cohort_list_large_perfire[[nm]]
  cohort_ysf <- unique(cohort_df$YSF)
  
  message("\n==============================")
  message("Mapping large fires one-by-one: ", cohort_ysf, " YSF")
  message("==============================")
  
  print_largefire_perfire_maps(
    cohort_df      = cohort_df,
    pixel_geom_sf  = pixel_geom_sf,
    huc12_bg       = huc12_bg,
    fire_atlas_bg  = fire_atlas_bg,
    streams_bg     = if (exists("streams_bg")) streams_bg else NULL,
    add_streams    = TRUE,
    fire_thresh    = 10000
  )
}
# ===============================================================

# ===============================================================
# PROPORTION UNRECOVERED BY IND. FIRE (Moderate & High Pixels only)
# ===============================================================

prop_unrec_by_fire <- dplyr::bind_rows(cohort_list, .id = "cohort_name") %>%
  dplyr::filter(sev_group %in% c("Moderate", "High")) %>%
  dplyr::group_by(YSF, fire_name, fire_size_ha, sev_group) %>%
  dplyr::summarise(
    N_total = dplyr::n(),
    N_unrec = sum(Recovered_by_last_year == 0L, na.rm = TRUE),
    prop_unrec = N_unrec / N_total,
    .groups = "drop"
  )

prop_table <- dplyr::bind_rows(cohort_list, .id = "cohort_name") %>%
  dplyr::filter(sev_group %in% c("Moderate", "High")) %>%
  dplyr::group_by(YSF, fire_name, fire_size_ha, sev_group) %>%
  dplyr::summarise(
    N_total = dplyr::n(),
    N_unrecovered = sum(Recovered_by_last_year == 0L, na.rm = TRUE),
    prop_unrecovered = N_unrecovered / N_total,
    .groups = "drop"
  ) %>%
  dplyr::arrange(YSF, fire_size_ha, fire_name, sev_group)

prop_table

# Optional: filter out tiny sample sizes
prop_plot_10 <- prop_table %>%
  dplyr::filter(
    YSF == 10,
    N_total >= 10
  )

ggplot(
  prop_plot_10,
  aes(
    x = reorder(fire_name, fire_size_ha),
    y = prop_unrecovered,
    fill = sev_group
  )
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  
  # ---- ADD LABELS 
geom_text(
  aes(label = N_total),
  position = position_dodge(width = 0.8),
  vjust = -0.3,
  size = 3
) +
  
  scale_fill_manual(values = c(
    "Moderate" = "orange",
    "High"     = "darkred"
  )) +
  labs(
    title = "Proportion Unrecovered by Fire (10 YSF)",
    x = "Fire (ordered by size)",
    y = "Proportion unrecovered",
    fill = "Severity"
  ) +
  ylim(0, 1.05) +  # ensures space for labels above bars
  theme_minimal(base_size = 13, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.grid.minor = element_blank()
  )

### Bar graph for 15 YSF
prop_plot_15 <- prop_table %>%
  dplyr::filter(
    YSF == 15,
    N_total >= 10
  )

ggplot(
  prop_plot_15,
  aes(
    x = reorder(fire_name, fire_size_ha),
    y = prop_unrecovered,
    fill = sev_group
  )
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  
  # ---- ADD LABELS 
geom_text(
  aes(label = N_total),
  position = position_dodge(width = 0.8),
  vjust = -0.3,
  size = 3
) +
  
  scale_fill_manual(values = c(
    "Moderate" = "orange",
    "High"     = "darkred"
  )) +
  labs(
    title = "Proportion Unrecovered by Fire (15 YSF)",
    x = "Fire (ordered by size)",
    y = "Proportion unrecovered",
    fill = "Severity"
  ) +
  ylim(0, 1.05) +
  theme_minimal(base_size = 13, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.grid.minor = element_blank()
  )

################################################################################


################################################################################
# ===============================================================
#  Explore whether UNRECOVERED Moderate + High pixels share
#  common FCODEs, slope, elevation, upstream drainage areas, etc.
#
#  Output:
#   (A) Cohorts rebuilt with new variables
#   (B) FCODE summary tables for unrecovered Moderate/High pixels
#   (C) FCODE bar plots (counts + proportions)
#   (D) Elevation summaries
#   (E) Elevation histograms with gray backdrop = all pixels in cohort,
#       colored bars = unrecovered Moderate/High pixels
# ===============================================================

stopifnot(exists("pixel_recovery_by_YSF"))
stopifnot(exists("data_long_complete"))

# Make sure these fields exist in data_long_complete
stopifnot("pixel_ID"     %in% names(data_long_complete))
stopifnot("elevation_m"  %in% names(data_long_complete))
stopifnot("slope_deg"  %in% names(data_long_complete))
stopifnot("TotDASqKm"  %in% names(data_long_complete))

# If your FCODE field is named slightly differently, change this here:
fcode_col <- "FCODE"
stopifnot(fcode_col %in% names(data_long_complete))

# -----------------------------
# 1) Build one-row-per-pixel lookup table for static attributes
# -----------------------------
pixel_attr_lookup <- data_long_complete %>%
  dplyr::mutate(pixel_ID = as.character(pixel_ID)) %>%
  dplyr::select(
    pixel_ID,
    elevation_m,
    swez_Apr,
    pptz_JJA,
    tmeanz_JJA,
    hli,
    TotDASqKm,
    slope_deg,
    dplyr::all_of(fcode_col)
  ) %>%
  dplyr::distinct(pixel_ID, .keep_all = TRUE) %>%
  dplyr::rename(FCODE = dplyr::all_of(fcode_col))

# Quick QA
pixel_attr_lookup %>%
  dplyr::summarise(
    n_pixels      = dplyr::n(),
    missing_elev  = sum(is.na(elevation_m)),
    missing_swe   = sum(is.na(swez_Apr)),
    missing_ppt   = sum(is.na(pptz_JJA)),
    missing_temp  = sum(is.na(tmeanz_JJA)),
    missing_hli   = sum(is.na(hli)),
    missing_TotDA = sum(is.na(TotDASqKm)),
    missing_fcode = sum(is.na(FCODE)),
    missing_slope = sum(is.na(slope_deg))
  )
# -----------------------------
# 2) Rebuild YSF cohorts 
# -----------------------------
build_clean_cohort_attr <- function(target_ysf) {
  
  cohort_raw <- pixel_recovery_by_YSF %>%
    dplyr::filter(YSF == target_ysf) %>%
    dplyr::select(YSF, pixel_ID, sev_group, fire_name, ref_year, yrs_to_recovery, fire_size_ha)
  
  cohort_one_row <- cohort_raw %>%
    dplyr::group_by(pixel_ID) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(pixel_attr_lookup, by = "pixel_ID") %>%
    dplyr::mutate(
      last_year_observed = as.integer(target_ysf),
      Recovered_by_last_year = dplyr::if_else(
        !is.na(yrs_to_recovery) & yrs_to_recovery <= target_ysf,
        1L, 0L
      )
    )
  
  stopifnot("Recovered_by_last_year" %in% names(cohort_one_row))
  stopifnot(all(cohort_one_row$Recovered_by_last_year %in% c(0L, 1L)))
  
  cohort_one_row
}

cohort_10_attr <- build_clean_cohort_attr(10L)
cohort_15_attr <- build_clean_cohort_attr(15L)
cohort_20_attr <- build_clean_cohort_attr(20L)

cohort_list_attr <- list(
  YSF10 = cohort_10_attr,
  YSF15 = cohort_15_attr,
  YSF20 = cohort_20_attr
)

unrec_mh_attr <- dplyr::bind_rows(cohort_list_attr, .id = "cohort_name") %>%
  dplyr::filter(
    Recovered_by_last_year == 0L,
    sev_group %in% c("Moderate", "High")
  ) %>%
  dplyr::mutate(
    sev_group = factor(sev_group, levels = c("Moderate", "High")),
    YSF = factor(YSF, levels = c(10, 15, 20))
  )

# -----------------------------
# 3) FCODE summaries for UNRECOVERED Moderate + High pixels
# -----------------------------
unrec_mh_attr <- dplyr::bind_rows(cohort_list_attr, .id = "cohort_name") %>%
  dplyr::filter(
    Recovered_by_last_year == 0L,
    sev_group %in% c("Moderate", "High")
  ) %>%
  dplyr::mutate(
    sev_group = factor(sev_group, levels = c("Moderate", "High")),
    YSF = factor(YSF, levels = c(10, 15, 20))
  )

# Count FCODEs within each YSF x severity
fcode_summary_unrec_mh <- unrec_mh_attr %>%
  dplyr::filter(!is.na(FCODE)) %>%
  dplyr::count(YSF, sev_group, FCODE, name = "N_pixels") %>%
  dplyr::group_by(YSF, sev_group) %>%
  dplyr::mutate(prop_within_group = N_pixels / sum(N_pixels)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(YSF, sev_group, dplyr::desc(N_pixels))

fcode_summary_unrec_mh

# Optional: top FCODEs only
top_fcodes <- fcode_summary_unrec_mh %>%
  dplyr::group_by(FCODE) %>%
  dplyr::summarise(total_N = sum(N_pixels), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(total_N)) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::pull(FCODE)

fcode_summary_top <- fcode_summary_unrec_mh %>%
  dplyr::filter(FCODE %in% top_fcodes)

fcode_summary_top

# -----------------------------
# 4) FCODE bar plots
# -----------------------------
fcode_summary_top_labeled <- fcode_summary_top %>%
  dplyr::mutate(
    FCODE_label = dplyr::case_when(
      FCODE == 46003 ~ "46003: Intermittent",
      FCODE == 46006 ~ "46006: Perennial",
      TRUE ~ as.character(FCODE)
    ),
    # keep ordering behavior the same
    FCODE_label = factor(FCODE_label)
  )

ggplot(
  fcode_summary_top,
  aes(x = reorder(FCODE, N_pixels), y = N_pixels, fill = sev_group)
) +
  geom_col(position = "dodge") +
  facet_wrap(~ YSF, scales = "free_y") +
  coord_flip() +
  scale_fill_manual(values = c(
    "Moderate" = "orange",
    "High"     = "darkred"
  )) +
  labs(
    title = "Top FCODEs among unrecovered Moderate and High-severity pixels",
    x = "FCODE",
    y = "Number of pixels",
    fill = "Severity"
  ) +
  theme_minimal(base_size = 13, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.minor = element_blank()
  )

ggplot(
  fcode_summary_top_labeled,
  aes(x = reorder(FCODE_label, prop_within_group),
      y = prop_within_group,
      fill = sev_group)
) +
  geom_col(position = "dodge") +
  facet_wrap(~ YSF, scales = "free_y") +
  coord_flip() +
  scale_fill_manual(values = c(
    "Moderate" = "orange",
    "High"     = "darkred"
  )) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Relative FCODE composition of unrecovered Moderate and High-severity pixels",
    x = "FCODE",
    y = "Proportion within group",
    fill = "Severity"
  ) +
  theme_minimal(base_size = 13, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.minor = element_blank()
  )

# -----------------------------
# 5) Elevation summaries
# Compare unrecovered M/H pixels against the full cohort
# -----------------------------
elev_summary_bg <- dplyr::bind_rows(cohort_list_attr, .id = "cohort_name") %>%
  dplyr::filter(!is.na(elevation_m)) %>%
  dplyr::group_by(YSF) %>%
  dplyr::summarise(
    N_all         = dplyr::n(),
    elev_mean_all = mean(elevation_m, na.rm = TRUE),
    elev_med_all  = median(elevation_m, na.rm = TRUE),
    elev_p25_all  = quantile(elevation_m, 0.25, na.rm = TRUE),
    elev_p75_all  = quantile(elevation_m, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

elev_summary_fg <- unrec_mh_attr %>%
  dplyr::filter(!is.na(elevation_m)) %>%
  dplyr::group_by(YSF, sev_group) %>%
  dplyr::summarise(
    N_unrec       = dplyr::n(),
    elev_mean     = mean(elevation_m, na.rm = TRUE),
    elev_med      = median(elevation_m, na.rm = TRUE),
    elev_p25      = quantile(elevation_m, 0.25, na.rm = TRUE),
    elev_p75      = quantile(elevation_m, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

elev_summary_bg
elev_summary_fg

# Optional: direct "high vs rest of cohort" / "moderate vs rest of cohort" check
elev_compare_table <- dplyr::bind_rows(cohort_list_attr, .id = "cohort_name") %>%
  dplyr::filter(!is.na(elevation_m)) %>%
  dplyr::mutate(
    group_compare = dplyr::case_when(
      Recovered_by_last_year == 0L & sev_group == "Moderate" ~ "Unrecovered Moderate",
      Recovered_by_last_year == 0L & sev_group == "High"     ~ "Unrecovered High",
      TRUE                                                   ~ "All other pixels"
    )
  ) %>%
  dplyr::group_by(YSF, group_compare) %>%
  dplyr::summarise(
    N          = dplyr::n(),
    elev_mean  = mean(elevation_m, na.rm = TRUE),
    elev_med   = median(elevation_m, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(YSF, group_compare)

elev_compare_table

# -----------------------------
# 6) Elevation histograms with gray backdrop
# Gray = all pixels in cohort
# Color = unrecovered Moderate/High only
# -----------------------------
elev_bg <- dplyr::bind_rows(cohort_list_attr, .id = "cohort_name") %>%
  dplyr::filter(
    YSF %in% c(10, 15, 20),
    !is.na(elevation_m)
  ) %>%
  dplyr::mutate(
    YSF = factor(YSF, levels = c(10, 15, 20))
  )

elev_fg <- unrec_mh_attr %>%
  dplyr::filter(
    YSF %in% c(10, 15, 20),
    !is.na(elevation_m)
  ) %>%
  dplyr::mutate(
    sev_group = factor(sev_group, levels = c("Moderate", "High")),
    YSF = factor(YSF, levels = c(10, 15, 20))
  )

# Overlaid histogram, same idea as fire-size version
ggplot() +
  geom_histogram(
    data = elev_bg,
    aes(x = elevation_m),
    bins = 30,
    fill = "grey80",
    color = "black",
    linewidth = 0.2,
    alpha = 0.7
  ) +
  geom_histogram(
    data = elev_fg,
    aes(x = elevation_m, fill = sev_group),
    bins = 30,
    color = "black",
    linewidth = 0.2,
    alpha = 0.75,
    position = "identity"
  ) +
  facet_wrap(~ YSF, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c(
    "Moderate" = "orange",
    "High"     = "darkred"
  )) +
  labs(
    title = "Elevation of unrecovered Moderate and High-severity pixels",
    subtitle = "Gray bars show all pixels in each cohort; colored bars show unrecovered Moderate and High pixels",
    x = "Elevation (m)",
    y = "Number of pixels",
    fill = "Severity"
  ) +
  theme_minimal(base_size = 13, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.minor = element_blank()
  )

# -----------------------------
# 7) Faceted elevation histograms by severity
# Helpful if Moderate and High overlap too much
# -----------------------------
# Restrict to Mod/High severity
elev_bg_mh <- elev_bg %>%
  dplyr::filter(sev_group %in% c("Moderate", "High"))

elev_fg_mh <- elev_fg %>%
  dplyr::filter(sev_group %in% c("Moderate", "High"))

## Plot histograms in single panel
ggplot() +
  geom_histogram(
    data = elev_bg_mh,
    aes(x = elevation_m),
    bins = 30,
    fill = "grey80",
    color = "black",
    linewidth = 0.2,
    alpha = 0.7
  ) +
  geom_histogram(
    data = elev_fg_mh,
    aes(x = elevation_m, fill = sev_group),
    bins = 30,
    color = "black",
    linewidth = 0.2,
    alpha = 0.8,
    position = "identity"
  ) +
  facet_grid(sev_group ~ YSF, scales = "free_y", drop = TRUE) +
  scale_fill_manual(values = c(
    "Moderate" = "orange",
    "High"     = "darkred"
  )) +
  labs(
    title = "Elevation of unrecovered Moderate and High-severity pixels",
    subtitle = "Gray = all Moderate/High pixels in cohort; color = unrecovered subset",
    x = "Elevation (m)",
    y = "Number of pixels"
  ) +
  theme_minimal(base_size = 13, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.minor = element_blank()
  )

# -----------------------------
# 8) Boxplots: Elevation, Recovery Status - M/H pixels
# -----------------------------
# N labels per box
n_labels <- elev_box_dat %>%
  dplyr::group_by(YSF, recovery_group) %>%
  dplyr::summarise(
    N = dplyr::n(),
    y_pos = max(elevation_m, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    label = paste0("n = ", N),
    y_pos = y_pos + 50   # adjust if needed depending on elevation range
  )


elev_box_dat <- dplyr::bind_rows(cohort_list_attr, .id = "cohort_name") %>%
  dplyr::filter(
    !is.na(elevation_m),
    YSF %in% c(10, 15, 20),
    sev_group %in% c("Moderate", "High")
  ) %>%
  dplyr::mutate(
    recovery_group = dplyr::case_when(
      sev_group == "Moderate" & Recovered_by_last_year == 1L ~ "Recovered Moderate",
      sev_group == "Moderate" & Recovered_by_last_year == 0L ~ "Unrecovered Moderate",
      sev_group == "High"     & Recovered_by_last_year == 1L ~ "Recovered High",
      sev_group == "High"     & Recovered_by_last_year == 0L ~ "Unrecovered High"
    ),
    recovery_group = factor(
      recovery_group,
      levels = c(
        "Recovered Moderate",
        "Unrecovered Moderate",
        "Recovered High",
        "Unrecovered High"
      )
    ),
    YSF = factor(YSF, levels = c(10, 15, 20))
  )

ggplot(elev_box_dat, aes(x = recovery_group, y = elevation_m, fill = recovery_group)) +
  geom_boxplot(outlier.alpha = 0.2) +
  facet_wrap(~ YSF, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c(
    "Recovered Moderate"   = "yellow2",
    "Unrecovered Moderate" = "darkorange",
    "Recovered High"       = "salmon",
    "Unrecovered High"     = "red4"
  )) +
  geom_text(
    data = n_labels,
    aes(
      x = recovery_group,
      y = y_pos,
      label = label
    ),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold",
    family = "Times New Roman"
  ) +
  labs(
    title = "Elevation of recovered vs unrecovered Moderate- and High-severity pixels",
    x = NULL,
    y = "Elevation (m)"
  ) +
  theme_minimal(base_size = 13, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1),
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.minor = element_blank()
  )

# -----------------------------
# 9) Scatterplots: Fire Size x Prop Unrecovered - M/H pixels
# -----------------------------
# Optional: filter out very small groups
prop_plot <- prop_table %>%
  dplyr::filter(
    YSF %in% c(10, 15),
    sev_group %in% c("Moderate", "High"),
    N_total >= 20,
    !is.na(fire_size_ha),
    fire_size_ha > 0,
    !is.na(prop_unrecovered)
  ) %>%
  dplyr::mutate(
    YSF = factor(YSF, levels = c(10, 15)),
    sev_group = factor(sev_group, levels = c("Moderate", "High")),
    log_fire_size = log10(fire_size_ha)
  )

# Compute panel-specific R² values
r2_df <- prop_plot %>%
  dplyr::group_by(YSF, sev_group) %>%
  dplyr::summarise(
    r2 = summary(lm(prop_unrecovered ~ log_fire_size, data = dplyr::cur_data()))$r.squared,
    x_pos = min(fire_size_ha, na.rm = TRUE),
    y_pos = 0.95,
    label = paste0("R² = ", sprintf("%.2f", r2)),
    .groups = "drop"
  )

ggplot(
  prop_plot,
  aes(x = fire_size_ha, y = prop_unrecovered)
) +
  geom_point(
    aes(color = sev_group),
    size = 2.8,
    alpha = 0.8
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ log10(x),
    se = FALSE,
    color = "darkgray",
    linewidth = 0.8
  ) +
  geom_text(
    data = r2_df,
    aes(x = x_pos, y = y_pos, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.5,
    color = "black"
  ) +
  scale_color_manual(values = c(
    "Moderate" = "orange",
    "High"     = "darkred"
  )) +
  scale_x_log10() +
  facet_grid(sev_group ~ YSF, scales = "fixed", drop = TRUE) +
  labs(
    title = "Fire Size vs. Proportion Unrecovered",
    x = "Fire size (ha, log10 scale)",
    y = "Proportion unrecovered",
    color = "Severity"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 13, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )
# -----------------------------
# 10) Recovery Binary Status Histograms for Climate Variables
# -----------------------------
# Snow Water Equivalent ==============
swe_bg <- dplyr::bind_rows(cohort_list_attr, .id = "cohort_name") %>%
  dplyr::filter(
    YSF %in% c(10, 15, 20),
    !is.na(swez_Apr)
  ) %>%
  dplyr::mutate(
    YSF = factor(YSF, levels = c(10, 15, 20))
  )

swe_fg <- unrec_mh_attr %>%
  dplyr::filter(
    YSF %in% c(10, 15, 20),
    !is.na(swez_Apr)
  ) %>%
  dplyr::mutate(
    sev_group = factor(sev_group, levels = c("Moderate", "High")),
    YSF = factor(YSF, levels = c(10, 15, 20))
  )

# Restrict background to Moderate/High only for cleaner comparison
swe_bg_mh <- swe_bg %>%
  dplyr::filter(sev_group %in% c("Moderate", "High"))

swe_fg_mh <- swe_fg %>%
  dplyr::filter(sev_group %in% c("Moderate", "High"))

ggplot() +
  geom_histogram(
    data = swe_bg_mh,
    aes(x = swez_Apr),
    bins = 30,
    fill = "grey80",
    color = "black",
    linewidth = 0.2,
    alpha = 0.7
  ) +
  geom_histogram(
    data = swe_fg_mh,
    aes(x = swez_Apr, fill = sev_group),
    bins = 30,
    color = "black",
    linewidth = 0.2,
    alpha = 0.8,
    position = "identity"
  ) +
  facet_grid(sev_group ~ YSF, scales = "free_y", drop = TRUE) +
  scale_fill_manual(values = c(
    "Moderate" = "orange",
    "High"     = "darkred"
  )) +
  labs(
    title = "April SWE of unrecovered Moderate and High-severity pixels",
    subtitle = "Gray = all Moderate/High pixels in cohort; color = unrecovered subset",
    x = "April SWE (swez_Apr)",
    y = "Number of pixels"
  ) +
  theme_minimal(base_size = 13, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.minor = element_blank()
  )

# Summer precipitation ===========
ppt_bg <- dplyr::bind_rows(cohort_list_attr, .id = "cohort_name") %>%
  dplyr::filter(
    YSF %in% c(10, 15, 20),
    !is.na(pptz_JJA)
  ) %>%
  dplyr::mutate(
    YSF = factor(YSF, levels = c(10, 15, 20))
  )

ppt_fg <- unrec_mh_attr %>%
  dplyr::filter(
    YSF %in% c(10, 15, 20),
    !is.na(pptz_JJA)
  ) %>%
  dplyr::mutate(
    sev_group = factor(sev_group, levels = c("Moderate", "High")),
    YSF = factor(YSF, levels = c(10, 15, 20))
  )

ppt_bg_mh <- ppt_bg %>%
  dplyr::filter(sev_group %in% c("Moderate", "High"))

ppt_fg_mh <- ppt_fg %>%
  dplyr::filter(sev_group %in% c("Moderate", "High"))

ggplot() +
  geom_histogram(
    data = ppt_bg_mh,
    aes(x = pptz_JJA),
    bins = 30,
    fill = "grey80",
    color = "black",
    linewidth = 0.2,
    alpha = 0.7
  ) +
  geom_histogram(
    data = ppt_fg_mh,
    aes(x = pptz_JJA, fill = sev_group),
    bins = 30,
    color = "black",
    linewidth = 0.2,
    alpha = 0.8,
    position = "identity"
  ) +
  facet_grid(sev_group ~ YSF, scales = "free_y", drop = TRUE) +
  scale_fill_manual(values = c(
    "Moderate" = "orange",
    "High"     = "darkred"
  )) +
  labs(
    title = "Summer precipitation of unrecovered Moderate and High-severity pixels",
    subtitle = "Gray = all Moderate/High pixels in cohort; color = unrecovered subset",
    x = "Summer precipitation (pptz_JJA)",
    y = "Number of pixels"
  ) +
  theme_minimal(base_size = 13, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.minor = element_blank()
  )

# Summer temperature ========================
temp_bg <- dplyr::bind_rows(cohort_list_attr, .id = "cohort_name") %>%
  dplyr::filter(
    YSF %in% c(10, 15, 20),
    !is.na(tmeanz_JJA)
  ) %>%
  dplyr::mutate(
    YSF = factor(YSF, levels = c(10, 15, 20))
  )

temp_fg <- unrec_mh_attr %>%
  dplyr::filter(
    YSF %in% c(10, 15, 20),
    !is.na(tmeanz_JJA)
  ) %>%
  dplyr::mutate(
    sev_group = factor(sev_group, levels = c("Moderate", "High")),
    YSF = factor(YSF, levels = c(10, 15, 20))
  )

temp_bg_mh <- temp_bg %>%
  dplyr::filter(sev_group %in% c("Moderate", "High"))

temp_fg_mh <- temp_fg %>%
  dplyr::filter(sev_group %in% c("Moderate", "High"))

ggplot() +
  geom_histogram(
    data = temp_bg_mh,
    aes(x = tmeanz_JJA),
    bins = 30,
    fill = "grey80",
    color = "black",
    linewidth = 0.2,
    alpha = 0.7
  ) +
  geom_histogram(
    data = temp_fg_mh,
    aes(x = tmeanz_JJA, fill = sev_group),
    bins = 30,
    color = "black",
    linewidth = 0.2,
    alpha = 0.8,
    position = "identity"
  ) +
  facet_grid(sev_group ~ YSF, scales = "free_y", drop = TRUE) +
  scale_fill_manual(values = c(
    "Moderate" = "orange",
    "High"     = "darkred"
  )) +
  labs(
    title = "Summer temperature of unrecovered Moderate and High-severity pixels",
    subtitle = "Gray = all Moderate/High pixels in cohort; color = unrecovered subset",
    x = "Summer temperature (tmeanz_JJA)",
    y = "Number of pixels"
  ) +
  theme_minimal(base_size = 13, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.minor = element_blank()
  )
# HLI =============================
hli_bg <- dplyr::bind_rows(cohort_list_attr, .id = "cohort_name") %>%
  dplyr::filter(
    YSF %in% c(10, 15, 20),
    !is.na(hli)
  ) %>%
  dplyr::mutate(
    YSF = factor(YSF, levels = c(10, 15, 20))
  )

hli_fg <- unrec_mh_attr %>%
  dplyr::filter(
    YSF %in% c(10, 15, 20),
    !is.na(hli)
  ) %>%
  dplyr::mutate(
    sev_group = factor(sev_group, levels = c("Moderate", "High")),
    YSF = factor(YSF, levels = c(10, 15, 20))
  )

# Restrict background to Moderate/High only for cleaner comparison
hli_bg_mh <- hli_bg %>%
  dplyr::filter(sev_group %in% c("Moderate", "High"))

hli_fg_mh <- hli_fg %>%
  dplyr::filter(sev_group %in% c("Moderate", "High"))

ggplot() +
  geom_histogram(
    data = hli_bg_mh,
    aes(x = hli),
    bins = 30,
    fill = "grey80",
    color = "black",
    linewidth = 0.2,
    alpha = 0.7
  ) +
  geom_histogram(
    data = hli_fg_mh,
    aes(x = hli, fill = sev_group),
    bins = 30,
    color = "black",
    linewidth = 0.2,
    alpha = 0.8,
    position = "identity"
  ) +
  facet_grid(sev_group ~ YSF, scales = "free_y", drop = TRUE) +
  scale_fill_manual(values = c(
    "Moderate" = "orange",
    "High"     = "darkred"
  )) +
  labs(
    title = "Heat Load Index of unrecovered Moderate and High-severity pixels",
    subtitle = "Gray = all Moderate/High pixels in cohort; color = unrecovered subset",
    x = "Heat Load Index (hli)",
    y = "Number of pixels"
  ) +
  theme_minimal(base_size = 13, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.minor = element_blank()
  )
# Drainage Area TotDASqKm =============================
TotDA_bg <- dplyr::bind_rows(cohort_list_attr, .id = "cohort_name") %>%
  dplyr::filter(
    YSF %in% c(10, 15),
    !is.na(TotDASqKm)
  ) %>%
  dplyr::mutate(
    YSF = factor(YSF, levels = c(10, 15))
  )

TotDA_fg <- unrec_mh_attr %>%
  dplyr::filter(
    YSF %in% c(10, 15),
    !is.na(TotDASqKm)
  ) %>%
  dplyr::mutate(
    sev_group = factor(sev_group, levels = c("Moderate", "High")),
    YSF = factor(YSF, levels = c(10, 15))
  )

## Filter for moderate/high severity only.
TotDA_bg_mh <- TotDA_bg %>%
  dplyr::filter(sev_group %in% c("Moderate", "High"))

TotDA_fg_mh <- TotDA_fg %>%
  dplyr::filter(sev_group %in% c("Moderate", "High"))

## Plot histogram
ggplot() +
  geom_histogram(
    data = TotDA_bg_mh,
    aes(x = TotDASqKm),
    bins = 30,
    fill = "grey80",
    color = "black",
    linewidth = 0.2,
    alpha = 0.7
  ) +
  geom_histogram(
    data = TotDA_fg_mh,
    aes(x = TotDASqKm, fill = sev_group),
    bins = 30,
    color = "black",
    linewidth = 0.2,
    alpha = 0.8,
    position = "identity"
  ) +
  facet_grid(sev_group ~ YSF, scales = "free_y", drop = TRUE) +
  scale_fill_manual(values = c(
    "Moderate" = "orange",
    "High"     = "darkred"
  )) +
  labs(
    title = "Drainage area of unrecovered Moderate and High-severity pixels",
    subtitle = "Gray = all Moderate/High pixels in cohort; color = unrecovered subset",
    x = "Total drainage area (km²)",
    y = "Number of pixels"
  ) +
  theme_minimal(base_size = 13, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.minor = element_blank()
  )

# Log scale:
ggplot() +
  geom_histogram(
    data = TotDA_bg_mh %>% dplyr::filter(TotDASqKm > 0),
    aes(x = TotDASqKm),
    bins = 30,
    fill = "grey80",
    color = "black",
    linewidth = 0.2,
    alpha = 0.7
  ) +
  geom_histogram(
    data = TotDA_fg_mh %>% dplyr::filter(TotDASqKm > 0),
    aes(x = TotDASqKm, fill = sev_group),
    bins = 30,
    color = "black",
    linewidth = 0.2,
    alpha = 0.8,
    position = "identity"
  ) +
  scale_x_log10() +
  facet_grid(sev_group ~ YSF, scales = "free_y", drop = TRUE) +
  scale_fill_manual(values = c(
    "Moderate" = "orange",
    "High"     = "darkred"
  )) +
  labs(
    title = "Drainage area of unrecovered Moderate and High-severity pixels",
    subtitle = "Gray = all Moderate/High pixels in cohort; color = unrecovered subset",
    x = "Total drainage area (km², log10 scale)",
    y = "Number of pixels"
  ) +
  theme_minimal(base_size = 13, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.minor = element_blank()
  )


# Slope =============================
slope_bg <- dplyr::bind_rows(cohort_list_attr, .id = "cohort_name") %>%
  dplyr::filter(
    YSF %in% c(10, 15),
    !is.na(slope_deg)
  ) %>%
  dplyr::mutate(
    YSF = factor(YSF, levels = c(10, 15))
  )

slope_fg <- unrec_mh_attr %>%
  dplyr::filter(
    YSF %in% c(10, 15),
    !is.na(slope_deg)
  ) %>%
  dplyr::mutate(
    sev_group = factor(sev_group, levels = c("Moderate", "High")),
    YSF = factor(YSF, levels = c(10, 15))
  )

# Restrict background to Moderate/High only for cleaner comparison
slope_bg_mh <- slope_bg %>%
  dplyr::filter(sev_group %in% c("Moderate", "High"))

slope_fg_mh <- slope_fg %>%
  dplyr::filter(sev_group %in% c("Moderate", "High"))

ggplot() +
  geom_histogram(
    data = slope_bg_mh,
    aes(x = slope_deg),
    bins = 30,
    fill = "grey80",
    color = "black",
    linewidth = 0.2,
    alpha = 0.7
  ) +
  geom_histogram(
    data = slope_fg_mh,
    aes(x = slope_deg, fill = sev_group),
    bins = 30,
    color = "black",
    linewidth = 0.2,
    alpha = 0.8,
    position = "identity"
  ) +
  facet_grid(sev_group ~ YSF, scales = "free_y", drop = TRUE) +
  scale_fill_manual(values = c(
    "Moderate" = "orange",
    "High"     = "darkred"
  )) +
  labs(
    title = "Slope of unrecovered Moderate and High-severity pixels",
    subtitle = "Gray = all Moderate/High pixels in cohort; color = unrecovered subset",
    x = "Slope (degrees)",
    y = "Number of pixels"
  ) +
  theme_minimal(base_size = 13, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.minor = element_blank()
  )

data_full <- pixel_recovery_by_YSF %>%
  dplyr::left_join(
    dplyr::bind_rows(cohort_list_attr, .id = "cohort_name") %>%
      dplyr::select(pixel_ID, swez_Apr, pptz_JJA, tmeanz_JJA, TotDASqKm, slope_deg, hli),
    by = "pixel_ID"
  )
data_rec <- data_full %>%
  dplyr::filter(!is.na(yrs_to_recovery))
################################################################################










################################################################################


# =========================================================================
# ============  MAPS: PIXELS COLORED BY RF COVARIATES  ====================
# =========================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(ggplot2)
  library(viridis)
})

# ---------------------------
# USER CONTROLS 
# ---------------------------
map_ctl <- list(
  # --- choose ONE or MANY targets (will return one map per target)
  target_ysf   = c(5,10,15,20),         # e.g., c(5,10,20) or c(20)
  
  # --- filters (set to "All" to skip)
  sev_group    = "All",         # "All" OR c("Unburned","Low","Moderate","High")
  ignition_year = "All",  #c(2003),     # "All",    OR    c(2003, 2017) etc  
  
  include_unburned = TRUE,      # if TRUE, keep fire_name == "None" pixels (ignition_year filter won't apply to those)
  burned_only      = FALSE,     # if TRUE, drop unburned regardless of include_unburned
  
  # --- color variable (must exist in data_long_complete)
  # from your list:
  # "ndvi_postfire","delta_ndvi_min","huc12","swez_Apr", pptz_JJA","tmeanz_JJA","tmaxz_JJA",
  # "veg_climate_index_08","cwd_5yr_zscore_08","twi","sev_num","fire_size_ha"
  color_var   = "swez_Apr",     # set to any variable above
  
  # --- plotting style
  point_size  = 3.5,
  point_alpha = 0.75,
  shape       = 21,             # 21 uses fill; 16 uses color (see below)
  stroke      = 0,              # outline width (0 removes outlines)
  
  # --- basemap style
  huc_fill    = "grey97",
  huc_color   = "grey70",
  huc_lwd     = 0.2,
  stream_color = "steelblue4",
  stream_lwd   = 0.6,
  stream_alpha = 0.7,
  
  # --- theme
  base_family = "Times New Roman",
  base_size   = 12,
  legend_pos  = "right"
)

# ---------------------------
# HUC12/Stream backdrops and clean variable labels
# ---------------------------
huc12_bg   <- huc12_sf   %>% sf::st_transform(sf::st_crs(pixel_geom_sf))
streams_bg <- streams_sf %>% sf::st_transform(sf::st_crs(pixel_geom_sf))
study_bbox <- sf::st_bbox(huc12_bg)

# Define variable labels for maps
var_labels <- c(
  "sev_num"              = "Fire Severity",
  "tmeanz_JJA"           = "Summer Mean Temperature",
  "twi"                  = "Topographic Wetness Index",
  "veg_climate_index_08" = "Vegetation–Climate Index",
  "cwd_5yr_zscore_08"    = "5-Year Post-Fire CWD",
  "swez_Apr"             = "Snow Water Equivalent",
  "pptz_JJA"             = "Summer Precipitation",
  "delta_ndvi_min"       = "ΔNDVI"
)

# ---------------------------
# BUILD: pixel attributes for mapping
# (one row per pixel; joins any variables you want to color/filter by)
# ---------------------------
# Make a pixel-level table with the union of variables we might use
map_vars_needed <- unique(c(
  "pixel_ID", "sev_group", "fire_name", "ref_year",
  map_ctl$color_var,
  # include common map vars in case you toggle later
  "ndvi_postfire","delta_ndvi_min","huc12","pptz_JJA","tmeanz_JJA","tmaxz_JJA",
  "veg_climate_index_08","cwd_5yr_zscore_08","twi","sev_num","fire_size_ha"
))

pixel_attr <- data_long_complete %>%
  dplyr::distinct(pixel_ID, .keep_all = TRUE) %>%
  dplyr::select(dplyr::any_of(map_vars_needed))

pixel_map_sf <- pixel_geom_sf %>%
  dplyr::left_join(pixel_attr, by = "pixel_ID") %>%
  dplyr::filter(!is.na(sev_group)) %>%
  sf::st_as_sf()

# ---------------------------
# Function: apply user filters to eligible cohort for a target_ysf
# ---------------------------
filter_elig_for_map <- function(target_ysf, ctl) {
  
  # eligible pixels by record length (from your Part 1C function)
  base_elig <- eligible_for(target_ysf)  # pixel_ID, sev_group, fire_name, ref_year
  
  if (nrow(base_elig) == 0L) return(base_elig)
  
  # burned/unburned handling
  if (isTRUE(ctl$burned_only)) {
    base_elig <- base_elig %>%
      dplyr::filter(!is.na(fire_name), fire_name != "None")
  } else if (!isTRUE(ctl$include_unburned)) {
    base_elig <- base_elig %>%
      dplyr::filter(!is.na(fire_name), fire_name != "None")
  }
  
  # sev_group filter
  if (!identical(ctl$sev_group, "All")) {
    base_elig <- base_elig %>%
      dplyr::filter(sev_group %in% ctl$sev_group)
  }
  
  # ignition year filter (applies to burned pixels; unburned retained only if include_unburned=TRUE & burned_only=FALSE)
  if (!identical(ctl$ignition_year, "All")) {
    yrs <- as.integer(ctl$ignition_year)
    
    if (isTRUE(ctl$burned_only) || !isTRUE(ctl$include_unburned)) {
      base_elig <- base_elig %>% dplyr::filter(ref_year %in% yrs)
    } else {
      base_elig <- base_elig %>%
        dplyr::filter((!is.na(fire_name) & fire_name != "None" & ref_year %in% yrs) |
                        (is.na(fire_name) | fire_name == "None"))
    }
  }
  
  base_elig
}

# ---------------------------
# Function: build map for target_ysf(s)
# ---------------------------

# Implement nice variable labels in function below
nice_var_label <- function(var, var_labels, tgt = NULL) {
  if (identical(var, "ndvi_postfire") && !is.null(tgt)) {
    return(paste0("NDVI at ", tgt, " YSF"))
  }
  if (!is.null(var_labels) && var %in% names(var_labels)) {
    return(var_labels[[var]])
  }
  var
}

# Mapping function:
make_map_for_target <- function(tgt, ctl) {
  
  # Caption text: ignition year(s) used
  ign_txt <- if (identical(ctl$ignition_year, "All")) {
    "Ignition year(s): All"
  } else {
    paste0("Ignition year(s): ", paste(ctl$ignition_year, collapse = ", "))
  }
  
  elig_ids <- filter_elig_for_map(tgt, ctl) %>% dplyr::pull(pixel_ID)
  
  # background-only if nothing passes filters
  if (length(elig_ids) == 0L) {
    return(
      ggplot() +
        geom_sf(data = huc12_bg, fill = ctl$huc_fill, color = ctl$huc_color, linewidth = ctl$huc_lwd) +
        geom_sf(data = streams_bg, color = ctl$stream_color, linewidth = ctl$stream_lwd, alpha = ctl$stream_alpha) +
        coord_sf(
          xlim = c(study_bbox["xmin"], study_bbox["xmax"]),
          ylim = c(study_bbox["ymin"], study_bbox["ymax"]),
          expand = FALSE
        ) +
        labs(
          title = paste0("Map — Target YSF ≥ ", tgt),
          subtitle = "No pixels after filters",
          caption = ign_txt,
          x = NULL, y = NULL
        ) +
        theme_minimal(base_family = ctl$base_family, base_size = ctl$base_size) +
        theme(legend.position = "none")
    )
  }
  
  # subset sf to eligible pixel_IDs
  mdat <- pixel_map_sf %>% dplyr::filter(pixel_ID %in% elig_ids)
  
  # summary counts for subtitle (burned-only stats)
  burned_sub <- mdat %>% dplyr::filter(!is.na(fire_name), fire_name != "None")
  N_pixels   <- nrow(mdat)
  N_fires    <- burned_sub %>% dplyr::distinct(fire_name) %>% nrow()
  N_years    <- burned_sub %>% dplyr::distinct(ref_year) %>% nrow()
  
  # color variable + pretty label
  var <- ctl$color_var
  var_pretty <- nice_var_label(var, var_labels, tgt = tgt)
  v <- mdat[[var]]
  is_cont <- is.numeric(v)
  
  # Use fill aesthetic for shape 21; use color for solid shapes like 16
  use_fill <- ctl$shape %in% c(21, 22, 23, 24, 25)
  
  p <- ggplot() +
    geom_sf(data = huc12_bg, fill = ctl$huc_fill, color = ctl$huc_color, linewidth = ctl$huc_lwd) +
    geom_sf(data = streams_bg, color = ctl$stream_color, linewidth = ctl$stream_lwd, alpha = ctl$stream_alpha)
  
  if (use_fill) {
    p <- p + geom_sf(
      data   = mdat,
      aes(fill = .data[[var]]),
      size   = ctl$point_size,
      alpha  = ctl$point_alpha,
      shape  = ctl$shape,
      stroke = ctl$stroke
    )
  } else {
    p <- p + geom_sf(
      data   = mdat,
      aes(color = .data[[var]]),
      size   = ctl$point_size,
      alpha  = ctl$point_alpha,
      shape  = ctl$shape
    )
  }
  
  # scales (use nice label for legend title)
  if (use_fill) {
    if (is_cont) {
      p <- p + scale_fill_viridis_c(name = var_pretty, option = "D", na.value = "transparent")
    } else {
      p <- p + scale_fill_viridis_d(name = var_pretty, option = "D", na.translate = FALSE)
    }
  } else {
    if (is_cont) {
      p <- p + scale_color_viridis_c(name = var_pretty, option = "D", na.value = "transparent")
    } else {
      p <- p + scale_color_viridis_d(name = var_pretty, option = "D", na.translate = FALSE)
    }
  }
  
  # coords + labels + theme (use pretty label in subtitle too)
  p +
    coord_sf(
      xlim = c(study_bbox["xmin"], study_bbox["xmax"]),
      ylim = c(study_bbox["ymin"], study_bbox["ymax"]),
      expand = FALSE
    ) +
    labs(
      title = paste0("BMWA Pixels — Target YSF ≥ ", tgt),
      subtitle = paste0(
        "Color = ", var_pretty,
        " • N(pixels) = ", scales::comma(N_pixels),
        " • N(fire events) = ", N_fires,
        " • N(ignition years) = ", N_years
      ),
      caption = ign_txt,
      x = NULL, y = NULL
    ) +
    theme_minimal(base_family = ctl$base_family, base_size = ctl$base_size) +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      legend.position  = ctl$legend_pos,
      plot.title       = element_text(face = "bold", size = 18, family = ctl$base_family),
      plot.subtitle    = element_text(size = 14, family = ctl$base_family),
      plot.caption = element_text(
        size   = 12,
        face   = "italic",
        family = ctl$base_family,
        color  = "grey30",
        hjust  = 0.8,
        angle  = 0,
        margin = ggplot2::margin(t = 10, b = 8)
      )
    )
}

# ---------------------------
# PRINT MAPS (one for each target_ysf)
# ---------------------------
maps_by_target <- lapply(map_ctl$target_ysf, function(tgt) make_map_for_target(tgt, map_ctl))
for (m in maps_by_target) print(m)

# Optional: access one map directly (e.g., first)
# maps_by_target[[1]]
################################################################################


################################################################################
### ===== HISTOGRAMS + PAIRED MAPS: YEARS TO RECOVERY BY TARGET_YSF COHORT =====
#### (NOTE: Run Parts 1A–1C for RQ2 first.) ####################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(showtext)
})

# ---------------------------
# HISTOGRAM USER CONTROLS 
# ---------------------------
hist_ctl <- list(
  ysf_targets = c(5, 10, 15, 20),     # YSF targets (5, 10, 15, 20)             
  sev_group   = c("High"),               # match your map cohort, or set c("High","Moderate"), etc.
  burned_only = TRUE,                    # TRUE = exclude "None"
  include_unburned = FALSE,              # usually FALSE if burned_only = TRUE
  
  n_bins_legend = 6,                     # stepped intervals
  base_family   = "Times New Roman",
  base_size     = 12
)

# Ensure Times New Roman font
font_add("Times New Roman", regular = "C:/Windows/Fonts/times.ttf")
showtext_auto()

# ---------------------------
# Function: compute breaks + palette mapping for ONE cohort table
# ---------------------------
make_ytr_breaks <- function(ytr_vec, n_bins) {
  ytr_vec <- ytr_vec[is.finite(ytr_vec)]
  if (!length(ytr_vec)) return(c(NA_real_, NA_real_))
  
  ytr_min <- min(ytr_vec, na.rm = TRUE)
  ytr_max <- max(ytr_vec, na.rm = TRUE)
  
  # integer-ish handling (your same logic)
  is_int <- all(abs(ytr_vec - round(ytr_vec)) < 1e-8, na.rm = TRUE)
  
  if (isTRUE(is_int)) {
    if ((ytr_max - ytr_min + 1) <= (n_bins + 1)) {
      breaks <- seq(ytr_min, ytr_max, by = 1)
    } else {
      breaks <- unique(round(seq(ytr_min, ytr_max, length.out = n_bins + 1)))
      breaks[1] <- ytr_min
      breaks[length(breaks)] <- ytr_max
    }
  } else {
    breaks <- seq(ytr_min, ytr_max, length.out = n_bins + 1)
  }
  
  breaks <- sort(unique(breaks))
  if (length(breaks) < 2) breaks <- c(ytr_min, ytr_max)
  breaks
}
# ---------------------------
# Function: build ONE histogram for a given target_ysf
# ---------------------------
build_hist_for_target <- function(tgt) {
  
  # 1) Pixel-level recovery info
  pix_attr <- data_long_complete %>%
    distinct(pixel_ID, sev_group, fire_name, ref_year, years_to_recovery)
  
  elig <- eligible_for(tgt) %>%
    left_join(pix_attr, by = c("pixel_ID", "sev_group", "fire_name", "ref_year")) %>%
    filter(sev_group %in% hist_ctl$sev_group)
  
  if (isTRUE(hist_ctl$burned_only)) {
    elig <- elig %>% filter(!is.na(fire_name), fire_name != "None")
  }
  
  # ---------------------------
  # 2) ADMINISTRATIVE CENSORING AT tgt
  # ---------------------------
  elig <- elig %>%
    mutate(
      years_to_recovery = as.integer(years_to_recovery),
      
      # censor recovery time at tgt
      ytr_censored = case_when(
        is.na(years_to_recovery)        ~ tgt + 1L,
        years_to_recovery > tgt         ~ tgt + 1L,
        TRUE                            ~ years_to_recovery
      )
    )
  
  # ---------------------------
  # 3) Breaks for recovered pixels ONLY (≤ tgt)
  # ---------------------------
  ytr_recovered <- elig$ytr_censored[elig$ytr_censored <= tgt]
  
  breaks_ytr <- make_ytr_breaks(ytr_recovered, hist_ctl$n_bins_legend)
  
  # Force last break to be exactly tgt
  breaks_ytr[length(breaks_ytr)] <- tgt
  
  # ---------------------------
  # 4) Bin assignment
  # ---------------------------
  elig <- elig %>%
    mutate(
      ytr_bin = case_when(
        ytr_censored > tgt ~ paste0(tgt + 1, "+"),
        TRUE ~ as.character(cut(
          ytr_censored,
          breaks = breaks_ytr,
          include.lowest = TRUE,
          right = TRUE
        ))
      )
    )
  
  # ---------------------------
  # 5) Palette (map-consistent)
  # ---------------------------
  bin_levels <- c(
    levels(factor(cut(ytr_recovered, breaks = breaks_ytr, include.lowest = TRUE))),
    paste0(tgt + 1, "+")
  )
  
  n_int <- length(bin_levels) - 1
  
  ramp_cols <- grDevices::colorRampPalette(
    c("darkgreen", "yellow4", "red3")
  )(n_int)
  
  pal <- setNames(
    c(ramp_cols, "darkred"),
    bin_levels
  )
  
  # ---------------------------
  # 6) Histogram table
  # ---------------------------
  hist_tbl <- elig %>%
    count(ytr_bin, name = "N") %>%
    mutate(ytr_bin = factor(ytr_bin, levels = bin_levels))
  
  # ---------------------------
  # 7) Plot
  # ---------------------------
  ggplot(hist_tbl, aes(x = ytr_bin, y = N, fill = ytr_bin)) +
    geom_col(color = "black", linewidth = 0.2, width = 0.9) +
    scale_fill_manual(values = pal, name = "Years to recovery") +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = paste0("Years to Recovery (YSF ≥ ", tgt, ")"),
      subtitle = paste0(
        "Recovery evaluated through ", tgt, " years post-fire • ",
        "Severity: ", paste(hist_ctl$sev_group, collapse = ", ")
      ),
      x = "Years to recovery",
      y = "Number of pixels"
    ) +
    theme_minimal(base_family = hist_ctl$base_family, base_size = hist_ctl$base_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
}
# ---------------------------
# Build and Print Histogram(s)
# ---------------------------
hist_by_target <- lapply(hist_ctl$ysf_targets, build_hist_for_target)

names(hist_by_target) <- paste0("YSF", hist_ctl$ysf_targets, "_hist")
for (p in hist_by_target) {
  print(p)
}
# ==============================================================================
################################################################################
### ===== MAPS: YEARS TO RECOVERY BINS (MATCH HISTOGRAMS) BY TARGET_YSF =========
#### (NOTE: Run Parts 1A–1C for RQ2 first.) ####################################

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(ggplot2)
  library(scales)
  library(showtext)
})

# ---------------------------
# MAPPING USER CONTROLS (re-uses histogram settings!)
# ---------------------------
map_ctl <- list(
  ysf_targets = hist_ctl$ysf_targets,     # reuse histogram targets
  sev_group   = hist_ctl$sev_group,       # reuse histogram sev_group
  burned_only = hist_ctl$burned_only,
  include_unburned = hist_ctl$include_unburned,
  
  n_bins_legend = hist_ctl$n_bins_legend, # reuse histogram bin count
  base_family   = hist_ctl$base_family,
  base_size     = hist_ctl$base_size,
  legend_pos    = "right",
  
  # point styling
  point_size  = 2.2,
  point_alpha = 0.85,
  
  # basemap styling
  huc_fill    = "grey97",
  huc_color   = "grey70",
  huc_lwd     = 0.2,
  stream_color = "steelblue4",
  stream_lwd   = 0.6,
  stream_alpha = 0.7
)

# Font (Windows)
font_add("Times New Roman", regular = "C:/Windows/Fonts/times.ttf")
showtext_auto()

# ---------------------------
# Backdrops (same CRS as pixels)
# ---------------------------
huc12_bg   <- huc12_sf   %>% st_transform(st_crs(pixel_geom_sf))
streams_bg <- streams_sf %>% st_transform(st_crs(pixel_geom_sf))
study_bbox <- st_bbox(huc12_bg)

# ---------------------------
# Reuse histogram break helper (must exist in env)
# make_ytr_breaks(ytr_vec, n_bins)
# ---------------------------
stopifnot(exists("make_ytr_breaks"))
stopifnot(exists("eligible_for"))

# ---------------------------
# Function: build ONE map for a given target_ysf (bins match histogram)
# ---------------------------
build_map_for_target <- function(tgt) {
  
  # 1) Pixel-level recovery info (one row per pixel)
  pix_attr <- data_long_complete %>%
    distinct(pixel_ID, sev_group, fire_name, ref_year, years_to_recovery)
  
  # 2) Eligible cohort + severity + burned/unburned handling
  elig <- eligible_for(tgt) %>%
    left_join(pix_attr, by = c("pixel_ID", "sev_group", "fire_name", "ref_year")) %>%
    filter(sev_group %in% map_ctl$sev_group)
  
  if (isTRUE(map_ctl$burned_only)) {
    elig <- elig %>% filter(!is.na(fire_name), fire_name != "None")
  } else if (!isTRUE(map_ctl$include_unburned)) {
    elig <- elig %>% filter(!is.na(fire_name), fire_name != "None")
  }
  
  # 3) ADMINISTRATIVE CENSORING AT tgt (exactly like histogram)
  elig <- elig %>%
    mutate(
      years_to_recovery = as.integer(years_to_recovery),
      ytr_censored = case_when(
        is.na(years_to_recovery) ~ tgt + 1L,
        years_to_recovery > tgt  ~ tgt + 1L,
        TRUE                     ~ years_to_recovery
      )
    )
  
  # 4) Breaks computed ONLY from recovered pixels (≤ tgt), same as histogram
  ytr_recovered <- elig$ytr_censored[elig$ytr_censored <= tgt]
  
  if (length(ytr_recovered) == 0L) {
    # No recovered pixels at all; still map the cohort (everything will be tgt+1+)
    # Create a trivial breaks vector just to keep downstream code happy
    breaks_ytr <- c(1, tgt)
  } else {
    breaks_ytr <- make_ytr_breaks(ytr_recovered, map_ctl$n_bins_legend)
    breaks_ytr[length(breaks_ytr)] <- tgt  # force last break to be exactly tgt
    breaks_ytr <- sort(unique(breaks_ytr))
    if (length(breaks_ytr) < 2) breaks_ytr <- c(min(ytr_recovered), tgt)
  }
  
  # 5) Bin assignment (same labels as histogram)
  elig <- elig %>%
    mutate(
      ytr_bin = case_when(
        ytr_censored > tgt ~ paste0(tgt + 1, "+"),
        TRUE ~ as.character(cut(
          ytr_censored,
          breaks = breaks_ytr,
          include.lowest = TRUE,
          right = TRUE
        ))
      )
    )
  
  # 6) Bin levels + palette (map-consistent; "tgt+1+" gets darkest red)
  if (length(ytr_recovered) > 0L) {
    bin_levels_core <- levels(factor(
      cut(ytr_recovered, breaks = breaks_ytr, include.lowest = TRUE, right = TRUE)
    ))
  } else {
    # no recovered bins exist
    bin_levels_core <- character(0)
  }
  bin_levels <- c(bin_levels_core, paste0(tgt + 1, "+"))
  
  n_int <- length(bin_levels)
  if (n_int <= 1L) {
    # Only the "tgt+1+" bin exists
    pal <- setNames("red4", bin_levels)
  } else {
    ramp_cols <- grDevices::colorRampPalette(c("darkgreen", "yellow4", "red4"))(n_int)
    pal <- setNames(ramp_cols, bin_levels)
    pal[ paste0(tgt + 1, "+") ] <- "red4"  # ensure darkest red for censored bin
  }
  
  # 7) Attach geometry
  map_sf <- pixel_geom_sf %>%
    select(pixel_ID, geometry) %>%
    inner_join(
      elig %>% select(pixel_ID, ytr_bin, ytr_censored, years_to_recovery, fire_name, ref_year, sev_group),
      by = "pixel_ID"
    ) %>%
    mutate(ytr_bin = factor(ytr_bin, levels = bin_levels)) %>%
    st_as_sf()
  
  N_map <- n_distinct(map_sf$pixel_ID)
  
  # 8) Plot
  ggplot() +
    geom_sf(data = huc12_bg, fill = map_ctl$huc_fill, color = map_ctl$huc_color, linewidth = map_ctl$huc_lwd) +
    geom_sf(data = streams_bg, color = map_ctl$stream_color, linewidth = map_ctl$stream_lwd, alpha = map_ctl$stream_alpha) +
    geom_sf(
      data = map_sf,
      aes(color = ytr_bin),
      shape = 16,
      size  = map_ctl$point_size,
      alpha = map_ctl$point_alpha
    ) +
    scale_color_manual(
      name   = "Years to Recovery",
      values = pal,
      drop   = FALSE
    ) +
    coord_sf(
      xlim = c(study_bbox["xmin"], study_bbox["xmax"]),
      ylim = c(study_bbox["ymin"], study_bbox["ymax"]),
      expand = FALSE
    ) +
    labs(
      title = paste0("Years to Recovery (YSF ≥ ", tgt, ")"),
      subtitle = paste0(
        "Severity: ", paste(map_ctl$sev_group, collapse = ", "),
        " • N = ", scales::comma(N_map)
      ),
      caption = paste0("Recovery evaluated through ", tgt, " years post-fire; ", tgt + 1, "+ = not recovered within window")
    ) +
    theme_minimal(base_family = map_ctl$base_family, base_size = map_ctl$base_size) +
    theme(
      legend.position = map_ctl$legend_pos,
      plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle   = element_text(size = 12, hjust = 0.5),
      plot.caption    = element_text(size = 10, face = "italic", margin = ggplot2::margin(t = 8))
    )
}

# ---------------------------
# Build + print one map per target_ysf (parallel to hist_by_target)
# ---------------------------
maps_by_target <- lapply(map_ctl$ysf_targets, build_map_for_target)
names(maps_by_target) <- paste0("YSF", map_ctl$ysf_targets, "_map")

for (m in maps_by_target) print(m)
# ==============================================================================
################################################################################
### ===== COMBINED PANELS: MAP + HISTOGRAM (MATCHED BINS) BY TARGET_YSF =========
#### (NOTE: Run Parts 1A–1C for RQ2 first.) ####################################
suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(ggplot2)
  library(scales)
  library(showtext)
  library(patchwork)
})

# ---------------------------
# USER CONTROLS (shared)
# ---------------------------
panel_ctl <- list(
  ysf_targets = hist_ctl$ysf_targets,
  sev_group   = hist_ctl$sev_group,
  burned_only = hist_ctl$burned_only,
  include_unburned = hist_ctl$include_unburned,
  
  n_bins_legend = hist_ctl$n_bins_legend,
  
  base_family = hist_ctl$base_family,
  base_size   = hist_ctl$base_size,
  
  # map styling
  point_size  = 2.2,
  point_alpha = 0.85,
  huc_fill    = "grey97",
  huc_color   = "grey70",
  huc_lwd     = 0.2,
  stream_color = "steelblue4",
  stream_lwd   = 0.6,
  stream_alpha = 0.7,
  legend_pos  = "right",
  
  # layout
  map_width  = 1.35,
  hist_width = 1.00
)

stopifnot(exists("eligible_for"))
stopifnot(exists("make_ytr_breaks"))

# Font (Windows)
font_add("Times New Roman", regular = "C:/Windows/Fonts/times.ttf")
showtext_auto()

# Backdrops (same CRS as pixels)
huc12_bg   <- huc12_sf   %>% st_transform(st_crs(pixel_geom_sf))
streams_bg <- streams_sf %>% st_transform(st_crs(pixel_geom_sf))
study_bbox <- st_bbox(huc12_bg)

# ---------------------------
# Function: build cohort table + bins + palette ONCE per target
# ---------------------------
build_cohort_binning <- function(tgt) {
  
  pix_attr <- data_long_complete %>%
    distinct(pixel_ID, sev_group, fire_name, ref_year, years_to_recovery)
  
  elig <- eligible_for(tgt) %>%
    left_join(pix_attr, by = c("pixel_ID", "sev_group", "fire_name", "ref_year")) %>%
    filter(sev_group %in% panel_ctl$sev_group)
  
  if (isTRUE(panel_ctl$burned_only)) {
    elig <- elig %>% filter(!is.na(fire_name), fire_name != "None")
  } else if (!isTRUE(panel_ctl$include_unburned)) {
    elig <- elig %>% filter(!is.na(fire_name), fire_name != "None")
  }
  
  elig <- elig %>%
    mutate(
      years_to_recovery = as.integer(years_to_recovery),
      ytr_censored = case_when(
        is.na(years_to_recovery) ~ tgt + 1L,
        years_to_recovery > tgt  ~ tgt + 1L,
        TRUE                     ~ years_to_recovery
      )
    )
  
  ytr_recovered <- elig$ytr_censored[elig$ytr_censored <= tgt]
  
  if (length(ytr_recovered) == 0L) {
    breaks_ytr <- c(1, tgt)
    bin_levels_core <- character(0)
  } else {
    breaks_ytr <- make_ytr_breaks(ytr_recovered, panel_ctl$n_bins_legend)
    breaks_ytr[length(breaks_ytr)] <- tgt
    breaks_ytr <- sort(unique(breaks_ytr))
    if (length(breaks_ytr) < 2) breaks_ytr <- c(min(ytr_recovered), tgt)
    
    bin_levels_core <- levels(factor(
      cut(ytr_recovered, breaks = breaks_ytr, include.lowest = TRUE, right = TRUE)
    ))
  }
  
  elig <- elig %>%
    mutate(
      ytr_bin = case_when(
        ytr_censored > tgt ~ paste0(tgt + 1, "+"),
        TRUE ~ as.character(cut(
          ytr_censored,
          breaks = breaks_ytr,
          include.lowest = TRUE,
          right = TRUE
        ))
      )
    )
  
  bin_levels <- c(bin_levels_core, paste0(tgt + 1, "+"))
  
  n_bins_total <- length(bin_levels)
  if (n_bins_total <= 1L) {
    pal <- setNames("red4", bin_levels)
  } else {
    ramp_cols <- grDevices::colorRampPalette(c("darkgreen", "yellow4", "red4"))(n_bins_total)
    pal <- setNames(ramp_cols, bin_levels)
    pal[paste0(tgt + 1, "+")] <- "red4"
  }
  
  list(
    elig = elig,
    tgt  = tgt,
    breaks_ytr = breaks_ytr,
    bin_levels = bin_levels,
    pal = pal
  )
}

# ---------------------------
# Functions - build map, histogram
# ---------------------------
# Build map
build_map_from_binning <- function(bin_obj) {
  
  elig <- bin_obj$elig
  tgt  <- bin_obj$tgt
  
  map_sf <- pixel_geom_sf %>%
    select(pixel_ID, geometry) %>%
    inner_join(
      elig %>% select(pixel_ID, ytr_bin, fire_name, ref_year, sev_group),
      by = "pixel_ID"
    ) %>%
    mutate(ytr_bin = factor(ytr_bin, levels = bin_obj$bin_levels)) %>%
    st_as_sf()
  
  N_map <- n_distinct(map_sf$pixel_ID)
  
  ggplot() +
    geom_sf(data = huc12_bg, fill = panel_ctl$huc_fill, color = panel_ctl$huc_color, linewidth = panel_ctl$huc_lwd) +
    geom_sf(data = streams_bg, color = panel_ctl$stream_color, linewidth = panel_ctl$stream_lwd, alpha = panel_ctl$stream_alpha) +
    geom_sf(
      data = map_sf,
      aes(color = ytr_bin),
      shape = 16,
      size  = panel_ctl$point_size,
      alpha = panel_ctl$point_alpha
    ) +
    scale_color_manual(
      name   = "Years to Recovery",
      values = bin_obj$pal,
      drop   = FALSE
    ) +
    coord_sf(
      xlim = c(study_bbox["xmin"], study_bbox["xmax"]),
      ylim = c(study_bbox["ymin"], study_bbox["ymax"]),
      expand = FALSE
    ) +
    labs(
      title = paste0("Map (YSF ≥ ", tgt, ")"),
      subtitle = paste0(
        "Severity: ", paste(panel_ctl$sev_group, collapse = ", "),
        " • N = ", scales::comma(N_map)
      ),
      caption = paste0("Recovery evaluated through ", tgt, " years; ", tgt + 1, "+ = not recovered within window")
    ) +
    theme_minimal(base_family = panel_ctl$base_family, base_size = panel_ctl$base_size) +
    theme(
      legend.position = panel_ctl$legend_pos,
      plot.title      = element_text(face = "bold", hjust = 0),
      plot.caption    = element_text(size = 9, face = "italic", margin = ggplot2::margin(t = 6))
    )
}

# Function: build histogram
build_hist_from_binning <- function(bin_obj) {
  
  elig <- bin_obj$elig
  tgt  <- bin_obj$tgt
  
  hist_tbl <- elig %>%
    mutate(ytr_bin = factor(ytr_bin, levels = bin_obj$bin_levels)) %>%
    count(ytr_bin, name = "N")
  
  ggplot(hist_tbl, aes(x = ytr_bin, y = N, fill = ytr_bin)) +
    geom_col(color = "black", linewidth = 0.2, width = 0.9) +
    scale_fill_manual(values = bin_obj$pal, drop = FALSE, name = "Years to recovery") +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = paste0("Histogram (YSF ≥ ", tgt, ")"),
      subtitle = paste0("Recovery evaluated through ", tgt, " years post-fire"),
      x = "Years to Recovery",
      y = "Number of pixels"
    ) +
    theme_minimal(base_family = panel_ctl$base_family, base_size = panel_ctl$base_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0)
    )
}

# ---------------------------
# Build one panel per target_ysf: map | histogram
# ---------------------------
panel_plots <- lapply(panel_ctl$ysf_targets, function(tgt) {
  bin_obj <- build_cohort_binning(tgt)
  
  p_map  <- build_map_from_binning(bin_obj)
  p_hist <- build_hist_from_binning(bin_obj)
  
  (p_map | p_hist) +
    plot_layout(widths = c(panel_ctl$map_width, panel_ctl$hist_width)) +
    plot_annotation(
      title = paste0("Years to Recovery (YSF ≥ ", tgt, ")"),
      theme = theme(
        text = element_text(family = panel_ctl$base_family),
        plot.title = element_text(face = "bold", size = 14)
      )
    )
})

names(panel_plots) <- paste0("YSF", panel_ctl$ysf_targets, "_panel")

# Print all panels
for (p in panel_plots) print(p)

# ---------------------------
# Re-build maps, adding Severity
# ---------------------------

# Read in severity TIF
library(terra)
sev_path <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/ArcGIS/BMWC_Fire_Data/fire.severity_LYB_freq/Fire Severity BMWC/fire.severity.mosaic.1985.2020.tif"
sev_r <- rast(sev_path)
# Reproject
sev_r <- project(sev_r, st_crs(pixel_geom_sf)$wkt)

# Crop to study area & mask
sev_r <- crop(sev_r, vect(huc12_bg))
sev_r <- mask(sev_r, vect(huc12_bg)) 

# Class raster for severity color bins
sev_class_r <- classify(
  sev_r,
  rcl = matrix(
    c(
      #-Inf, 0,    NA,  # set anything <0 to NA (optional safety)
      0,    0,    0,   # exactly 0 stays 0
      0,    1.25, 1,   # >0 up to 1.25
      1.25, 2.25, 2,   # >1.25 up to 2.25
      2.25, Inf,  3    # >2.25
    ),
    ncol = 3,
    byrow = TRUE
  ),
  include.lowest = TRUE,
  right = TRUE
)

# Convert to data frame for ggplot
sev_df <- as.data.frame(sev_class_r, xy = TRUE, na.rm = TRUE)
names(sev_df)[3] <- "sev_class"

sev_df$sev_class <- factor(
  sev_df$sev_class,
  levels = c(0, 1, 2, 3),
  labels = c(
    "0",
    ">0–1.25",
    ">1.25–2.25",
    ">2.25"
  )
)

# Assign severity class colors
sev_scale <- scale_fill_manual(
  name   = "Fire severity",
  values = c(
    "0"         = "grey60",
    ">0–1.25"   = "khaki2",
    ">1.25–2.25"= "orange2",
    ">2.25"     = "darkred"
  ),
  drop = FALSE,
  na.value = NA
)
# ----------------------------------------
## NEW: Function to build maps, with severity option available 
# ----------------------------------------
build_map_from_binning <- function(bin_obj, 
                                   sev_df = NULL, 
                                   sev_scale = NULL, 
                                   sev_alpha = 0.45) {
  
  elig <- bin_obj$elig
  tgt  <- bin_obj$tgt
  
  map_sf <- pixel_geom_sf %>%
    dplyr::select(pixel_ID, geometry) %>%
    dplyr::inner_join(
      elig %>% dplyr::select(pixel_ID, ytr_bin, fire_name, ref_year, sev_group),
      by = "pixel_ID"
    ) %>%
    dplyr::mutate(ytr_bin = factor(ytr_bin, levels = bin_obj$bin_levels)) %>%
    sf::st_as_sf()
  
  N_map <- dplyr::n_distinct(map_sf$pixel_ID)
  
  p <- ggplot()
  
  # ---- severity raster background (optional)
  if (!is.null(sev_df)) {
    p <- p +
      geom_raster(
        data = sev_df,
        aes(x = x, y = y, fill = sev_class),
        alpha = sev_alpha
      )
    
    # add fill scale if provided
    if (!is.null(sev_scale)) {
      p <- p + sev_scale
    }
  }
  
  # ---- vector basemap + pixels
  p +
    geom_sf(data = huc12_bg, fill = NA, color = panel_ctl$huc_color, linewidth = panel_ctl$huc_lwd) +
    geom_sf(data = streams_bg, color = panel_ctl$stream_color, linewidth = panel_ctl$stream_lwd, alpha = panel_ctl$stream_alpha) +
    geom_sf(
      data = map_sf,
      aes(color = ytr_bin),
      shape = 16,
      size  = panel_ctl$point_size,
      alpha = panel_ctl$point_alpha
    ) +
    scale_color_manual(
      name   = "Years to Recovery",
      values = bin_obj$pal,
      drop   = FALSE
    ) +
    coord_sf(
      xlim = c(study_bbox["xmin"], study_bbox["xmax"]),
      ylim = c(study_bbox["ymin"], study_bbox["ymax"]),
      expand = FALSE
    ) +
    labs(
      title = paste0("Map (YSF ≥ ", tgt, ")"),
      subtitle = paste0(
        "Severity: ", paste(panel_ctl$sev_group, collapse = ", "),
        " • N = ", scales::comma(N_map)
      ),
      caption = paste0("Recovery evaluated through ", tgt, " years; ", tgt + 1, "+ = not recovered within window")
    ) +
    theme_minimal(base_family = panel_ctl$base_family, base_size = panel_ctl$base_size) +
    theme(
      legend.position = panel_ctl$legend_pos,
      plot.title      = element_text(face = "bold", hjust = 0),
      plot.caption    = element_text(size = 9, face = "italic", margin = ggplot2::margin(t = 6))
    )
}

# ----------------------------------------
# Rebuild maps for all target YSF (with upland severity) 
# ----------------------------------------
panel_plots_sev <- lapply(panel_ctl$ysf_targets, function(tgt) {
  bin_obj <- build_cohort_binning(tgt)
  
  p_map  <- build_map_from_binning(bin_obj, sev_df = sev_df, sev_scale = sev_scale, sev_alpha = 0.45)
  p_hist <- build_hist_from_binning(bin_obj)
  
  (p_map | p_hist) +
    plot_layout(widths = c(panel_ctl$map_width, panel_ctl$hist_width)) +
    plot_annotation(
      title = paste0("Years to Recovery (YSF ≥ ", tgt, ")"),
      theme = theme(
        text = element_text(family = panel_ctl$base_family),
        plot.title = element_text(face = "bold", size = 14)
      )
    )
})

names(panel_plots_sev) <- paste0("YSF", panel_ctl$ysf_targets, "_panel_sev")
for (p in panel_plots_sev) print(p)
# ------------------------------------------------------------------------------

################################################################################
### ===== SAVE INTERACTIVE LEAFLET MAPS (HTML): Severity raster + HUC12 + Streams
###       + Pixels colored by YTR bin (matched to histogram bins)
suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(terra)
  library(leaflet)
  library(htmlwidgets)
})

# ---------------------------
# USER CONTROLS (HTML export)
# ---------------------------
html_ctl <- list(
  ysf_targets = panel_ctl$ysf_targets,
  out_dir     = "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/ArcGIS/BMWC_Fire_Data/CH01_file_drop",
  html_prefix = "YTR_Severity",
  selfcontained = TRUE,   # set FALSE if file size is too large
  pixel_radius  = 3,
  pixel_opacity = 0.85,
  sev_opacity   = 0.45,
  
  # severity colors (DISCRETE)
  sev_colors = c(
    "0"         = "grey60",
    ">0–1.25"   = "khaki2",
    ">1.25–2.25"= "orange2",
    ">2.25"     = "darkred"
  )
)

# ---------------------------
# Ensure backdrops are WGS84 for leaflet (EPSG:4326)
# ---------------------------
huc12_ll   <- st_transform(huc12_bg, 4326)
streams_ll <- st_transform(streams_bg, 4326)

# ---------------------------
# Read + project + crop/mask severity raster, then CLASS to 4 bins
# ---------------------------
sev_path <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/ArcGIS/BMWC_Fire_Data/fire.severity_LYB_freq/Fire Severity BMWC/fire.severity.mosaic.1985.2020.tif"

sev_r <- rast(sev_path)

# Project raster to WGS84 for leaflet
sev_r_ll <- project(sev_r, "EPSG:4326")

# Crop/mask to study area (use huc12 in WGS84)
sev_r_ll <- crop(sev_r_ll, vect(huc12_ll))
sev_r_ll <- mask(sev_r_ll, vect(huc12_ll))

# ---- classify to discrete classes: 0, (0–1.25], (1.25–2.25], (2.25+)
# IMPORTANT: your original rcl had overlapping boundaries; this is clean + non-overlapping:
# class 0: exactly 0
# class 1: >0 to 1.25
# class 2: >1.25 to 2.25
# class 3: >2.25
sev_class_r <- classify(
  sev_r_ll,
  rcl = matrix(
    c(
      #-Inf, 0,     NA,  # optional: anything <0 -> NA
      0,    0,     0,   # exactly 0
      0,    1.25,  1,   # >0 to 1.25 (we'll treat 0 separately via exact rule above)
      1.25, 2.25,  2,   # >1.25 to 2.25
      2.25, Inf,   3    # >2.25
    ),
    ncol = 3,
    byrow = TRUE
  ),
  include.lowest = TRUE,
  right = TRUE
)

# Convert 0/1/2/3 -> labeled factor via raster levels 
sev_levels <- data.frame(
  value = c(0, 1, 2, 3),
  label = c("0", ">0–1.25", ">1.25–2.25", ">2.25")
)
levels(sev_class_r) <- list(sev_levels)

sev_pal_fun <- colorFactor(
  palette = html_ctl$sev_colors,
  domain  = sev_levels$label,
  ordered = TRUE
)

# addRasterImage() wants a color function that accepts raster values.
# Because we set raster levels, we can map numeric -> label -> color:
sev_color_fun <- function(values) {
  labs <- sev_levels$label[match(values, sev_levels$value)]
  sev_pal_fun(labs)
}

# ---------------------------
# Function: build ONE YSF leaflet map and save HTML
# (Requires build_cohort_binning(tgt) from code above)
# ---------------------------
save_leaflet_for_target <- function(tgt) {
  
  # ---- cohort binning object (your existing function)
  bin_obj <- build_cohort_binning(tgt)
  
  # ---- build sf for pixels (include attributes for popup)
  map_sf <- pixel_geom_sf %>%
    dplyr::select(pixel_ID, geometry) %>%
    dplyr::inner_join(
      bin_obj$elig %>% dplyr::select(pixel_ID, ytr_bin, years_to_recovery, ytr_censored,
                                     fire_name, ref_year, sev_group),
      by = "pixel_ID"
    ) %>%
    dplyr::mutate(
      ytr_bin = factor(ytr_bin, levels = bin_obj$bin_levels)
    ) %>%
    sf::st_as_sf() %>%
    sf::st_transform(4326)
  
  N_map <- dplyr::n_distinct(map_sf$pixel_ID)
  
  # ---- YTR palette (already computed per target; names are bin labels)
  pal_vec <- bin_obj$pal
  pal_fun <- colorFactor(palette = pal_vec, domain = names(pal_vec), ordered = TRUE)
  
  # ---- leaflet
  m <- leaflet(options = leafletOptions(preferCanvas = TRUE)) %>%
    addProviderTiles("CartoDB.Positron", group = "Basemap") %>%
    
    # severity raster (background)
    addRasterImage(
      sev_class_r,
      colors  = sev_color_fun,
      opacity = html_ctl$sev_opacity,
      group   = "Fire severity (binned)",
      project = FALSE   # already in EPSG:4326
    ) %>%
    
    # HUC + streams
    addPolygons(
      data = huc12_ll,
      fill = FALSE, color = "grey70", weight = 1,
      group = "HUC12 boundary"
    ) %>%
    addPolylines(
      data = streams_ll,
      color = "steelblue4", weight = 1, opacity = 0.7,
      group = "Streams"
    ) %>%
    
    # pixels
    addCircleMarkers(
      data = map_sf,
      radius      = html_ctl$pixel_radius,
      stroke      = FALSE,
      fillOpacity = html_ctl$pixel_opacity,
      color       = ~pal_fun(ytr_bin),
      group       = "Pixels (years to recovery)",
      popup = ~paste0(
        "<b>pixel_ID:</b> ", pixel_ID,
        "<br><b>YSF target:</b> ", tgt,
        "<br><b>Years-to-recovery bin:</b> ", as.character(ytr_bin),
        "<br><b>yrs_to_recovery (raw):</b> ", ifelse(is.na(years_to_recovery), "NA", years_to_recovery),
        "<br><b>yrs_to_recovery (censored):</b> ", ytr_censored,
        "<br><b>sev_group:</b> ", sev_group,
        "<br><b>fire_name:</b> ", fire_name,
        "<br><b>ref_year:</b> ", ref_year
      )
    ) %>%
    
    # legends
    addLegend(
      position = "topright",
      title    = paste0("Years to recovery (YSF ≥ ", tgt, ")<br>N = ", scales::comma(N_map)),
      colors   = unname(pal_vec),
      labels   = names(pal_vec),
      opacity  = 1,
      group    = "Pixels (years to recovery)"
    ) %>%
    addLegend(
      position = "bottomright",
      title    = "Fire severity (binned)",
      colors   = unname(html_ctl$sev_colors),
      labels   = names(html_ctl$sev_colors),
      opacity  = 1,
      group    = "Fire severity (binned)"
    ) %>%
    
    # layer controls (toggle raster/streams/pixels)
    addLayersControl(
      overlayGroups = c("Fire severity (binned)", "HUC12 boundary", "Streams", "Pixels (years to recovery)"),
      options = layersControlOptions(collapsed = FALSE)
    )
  
  out_file <- file.path(
    html_ctl$out_dir,
    paste0(html_ctl$html_prefix, "_YSF", tgt, "_", paste(panel_ctl$sev_group, collapse = "-"), ".html")
  )
  
  htmlwidgets::saveWidget(m, out_file, selfcontained = html_ctl$selfcontained)
  message("Saved: ", out_file)
  
  invisible(m)
}

# ---------------------------
# 3) Save one HTML per target_ysf
# ---------------------------
#dir.create(html_ctl$out_dir, recursive = TRUE, showWarnings = FALSE)
#leaflet_maps <- lapply(html_ctl$ysf_targets, save_leaflet_for_target)
#names(leaflet_maps) <- paste0("YSF", html_ctl$ysf_targets, "_leaflet")

################################################################################





################################################################################









# ==============================================================================
# ==============================================================================
## FINAL DATA FRAME EDITS PRIOR TO EXPORT #######
aet_monthly_cols <- grep("^AET_\\d{6}_TC$", names(data_long_complete), value = TRUE)
cols_to_drop <- unique(c(aet_monthly_cols))
# Drop monthly AET columns (no longer needed)
data_long_complete <- data_long_complete %>%
  dplyr::select(-dplyr::any_of(cols_to_drop))

# Rename the data frame
data_long_rfmodel <- data_long_complete

write_csv(data_long_rfmodel,
          "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/data_long_rfmodel_03032026_40mSpacing.csv")

################################################################################