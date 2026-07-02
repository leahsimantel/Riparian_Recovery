################################################################################
###    CH01 Script C: Env Blocking and RF Modeling 
###    K-means clustering
###    40-meter spacing 
###    12/18/2025


# ==============================================================
# Libraries
# ==============================================================
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
  library(randomForest)
  library(cluster)      
  library(RColorBrewer)
  library(corrplot)
  library(car)
  library(forcats)
  library(lme4)
  library(purrr)
  library(patchwork)
  library(readr)
  library(stringr)
})
# -------------------------------------------------------------

# ==============================================================
# Upload CSV and rename it
# ==============================================================
csv_path <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/data_long_rfmodel_12182025_40mSpacing.csv"

## false spacing (~60m) 10k pixels: data_long_rfmodel_09102025.csv 
## fixed 100m spacing:  data_long_rfmodel_11052025_100mSpacingFR.csv  
# 40-m spacing:         data_long_rfmodel_12182025_40mSpacing.csv

df <- readr::read_csv(csv_path, show_col_types = FALSE)  # tibble
#dim(df); names(df)[1:10]  

## Rename file for downstream code:
data_long_complete <- df
# --------------------------------------------------------------

# ==============================================================
# Adding New Columns, 2/19/26: Elevation & TopoTerra AET, CWD, Slope, DA
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

data_long_complete <- data_long_complete %>%
  mutate(pixel_ID = as.character(pixel_ID))
data_long_complete <- data_long_complete %>%
  left_join(covars_to_join, by = "pixel_ID")

####  QA - confirm no extra rows were added to data_long_complete
#n_before <- nrow(data_long_complete)
# (re-run join in a temporary object if needed)
#data_long_complete_new <- data_long_complete %>%
#  left_join(covars_to_join, by = "pixel_ID")

#n_after <- nrow(data_long_complete_new)
#n_before == n_after

data_long_complete %>%
  summarise(
    missing_aet = sum(is.na(aet_30yrAvg_TT)),
    missing_cwd = sum(is.na(cwd_30yrAvg_TT)),
    missing_elev = sum(is.na(elevation_m)),
    missing_slope = sum(is.na(slope_deg)),
    missing_TotDA = sum(is.na(TotDASqKm))
  )

## QA... mapping out pixels with missing NEW variables
# Pixels in data_long_complete that failed to match
missing_pixels <- data_long_complete %>%
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

data_long_complete <- data_long_complete %>%
  filter(!pixel_ID %in% pixels_to_remove)
# --------------------------------------------------------------

# ==============================================================
# Data Prep
# ==============================================================

# ---- Sample size summaries 
total_n <- data_long_complete %>%
  filter(!is.na(delta_ndvi_min)) %>%
  group_by(years_since_fire) %>%
  summarise(Total_N = n(), .groups = "drop")

sev_counts <- data_long_complete %>%
  filter(!is.na(delta_ndvi_min)) %>%
  st_drop_geometry() %>%
  group_by(years_since_fire, sev_group) %>%
  summarise(N = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = sev_group,
    values_from = N,
    values_fill = 0
  ) %>%
  rename(
    Unburned_N = `Unburned`,
    Low_N      = `Low`,
    Moderate_N = `Moderate`,
    High_N     = `High`
  )

ysf_summary <- total_n %>%
  left_join(sev_counts, by = "years_since_fire") %>%
  arrange(years_since_fire)

# ---- Ensure an sf object 
if (!inherits(data_long_complete, "sf")) {
  data_long_complete[["longitude"]] <- suppressWarnings(as.numeric(data_long_complete[["longitude"]]))
  data_long_complete[["latitude"]]  <- suppressWarnings(as.numeric(data_long_complete[["latitude"]]))
  
  data_long_reduced_sf <- sf::st_as_sf(
    data_long_complete,
    coords  = c("longitude", "latitude"),
    crs     = 4326,
    remove  = FALSE
  )
} else {
  data_long_reduced_sf <- data_long_complete
}
# --------------------------------------------------------------


##### ====== SCATTERPLOTS: deltaNDVI against Predictor Variables ====== 
## Plotting Prep: set filters and variables 
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork); library(rlang)
})

target_ysf_vals <- c(10, 15)                       # YSFs to plot
sev_keep        <- c("Low", "Moderate", "High")  # which sev_group to include

sev_colors <- c(
  "Unburned" = "grey60",
  "Low"      = "yellow",
  "Moderate" = "orange",
  "High"     = "firebrick"
)

# Variables to iterate over
vars_to_plot <- c(
  "pptz_JJA",
  #"ppttot_JJA",
  #"tmeanz_JJA", 
  "tmaxz_JJA", 
  #"tmeanavg_JJA", "tmaxmean_JJA",
  #"swe_peak_MAM", 
  #"SWE_Apr", 
  "swez_Apr",
  "cwd_5yr_zscore_08",
  "twi",
  "sev_num",
  "fire_size_ha",
  "elevation_m"
)

# x-axis labels 
var_labels <- c(
  pptz_JJA            = "JJA Precip Anomaly (z-score)",
  #ppttot_JJA          = "JJA Precip",
  tmeanz_JJA          = "JJA Mean Temp Anomaly (z-score)",
  tmaxz_JJA           = "JJA Max Temp Anomaly (z-score)",
  #tmeanavg_JJA        = "JJA Mean Temp (°C)",
  #tmaxmean_JJA        = "JJA Max Temp (°C)",
  swez_Apr            = "Apr SWE Anomaly (z-score)",
  #SWE_Apr             = "April SWE (mm)",
  cwd_5yr_zscore_08   = "Post-fire 5-yr CWD (Aug; z-score)",
  twi                 = "Topographic Wetness Index",
  sev_num             = "Fire Severity (numeric)",
  fire_size_ha        = "Fire Size (hectares)",
  elevation_m         = "Elevation (m)"
)

# ---- Base data prep (once): filter severity, coerce ΔNDVI, keep needed 
dlr_base <- data_long_complete %>%
  filter(sev_group %in% sev_keep) %>%
  mutate(
    delta_ndvi_min = suppressWarnings(as.numeric(delta_ndvi_min)),
    sev_group  = factor(sev_group, levels = names(sev_colors))
  ) %>%
  filter(!is.na(delta_ndvi_min))

# ---- Function to build one scatter for a given variable and YSF 
plot_one_ysf <- function(df, xvar, yvar = "delta_ndvi_min", ysf) {
  df_y <- df %>% filter(years_since_fire == ysf)
  
  # coerce the xvar on the fly
  x_sym <- sym(xvar)
  df_y <- df_y %>%
    mutate(!!x_sym := suppressWarnings(as.numeric(.data[[xvar]]))) %>%
    filter(!is.na(.data[[xvar]]))
  
  n_pts <- nrow(df_y)
  r_val <- suppressWarnings(cor(df_y[[xvar]], df_y[[yvar]], use = "complete.obs"))
  
  ggplot(df_y, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_point(aes(color = sev_group), alpha = 0.25, size = 1) +
    scale_color_manual(values = sev_colors, drop = TRUE) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    labs(
      title    = paste0("YSF = ", ysf),
      subtitle = paste0("N = ", n_pts, " | Pearson r = ",
                        ifelse(is.finite(r_val), sprintf("%.3f", r_val), "NA")),
      x = unname(ifelse(xvar %in% names(var_labels), var_labels[[xvar]], xvar)),
      y = "ΔNDVI",
      color = "Severity"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title   = element_text(face = "bold"),
      plot.margin  = grid::unit(c(10, 20, 10, 10), "pt"),
      legend.position = "right"
    )
}

# ---- Wrapper: make a 2x2 panel for a single variable across YSFs 
panel_for_var <- function(xvar) {
  # skip gracefully if column missing
  if (!xvar %in% names(dlr_base)) {
    message("Skipping '", xvar, "' (column not found).")
    return(NULL)
  }
  
  plots <- lapply(target_ysf_vals, function(y) plot_one_ysf(dlr_base, xvar, "delta_ndvi_min", y))
  # remove NULLs in case any YSF had no data
  plots <- Filter(Negate(is.null), plots)
  
  wrap_plots(plots, ncol = 2) +
    plot_annotation(
      title = paste0("ΔNDVI vs ", unname(ifelse(xvar %in% names(var_labels), var_labels[[xvar]], xvar))),
      subtitle = paste("Severities included:", paste(sev_keep, collapse = ", "))
    )
}

# ---- Generate all panels 
panels <- lapply(vars_to_plot, panel_for_var)
names(panels) <- vars_to_plot

# ---- Print panels one-by-one (or index into `panels[['pptz_JJA']]`) 
for (v in names(panels)) {
  if (!is.null(panels[[v]])) {
    print(panels[[v]])
  }
}
# ===============================================================

# ---- Scatterplot: Fire size vs sev_num (burned pixels only) ----
df_fs <- data_long_complete %>%
  mutate(
    fire_size_ha = as.numeric(fire_size_ha),
    sev_num      = as.numeric(sev_num)
  ) %>%
  filter(
    !is.na(fire_size_ha),
    !is.na(sev_num),
    sev_num > 0,                 # drop unburned
    sev_group != "Unburned"
  ) %>%
  distinct(fire_name, .keep_all = TRUE)  # one row per fire

# linear model
lm_fs <- lm(sev_num ~ fire_size_ha, data = df_fs)
r2_fs <- summary(lm_fs)$r.squared

ggplot(df_fs, aes(x = fire_size_ha, y = sev_num)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  labs(
    x = "Fire Size (ha)",
    y = "Fire Severity (numeric)",
    title = "Fire Size vs Fire Severity (Burned Riparian Only)",
    subtitle = paste0("R² = ", sprintf("%.3f", r2_fs))
  ) +
  theme_minimal(base_size = 12)
######################################


##### SCATTERPLOTS: COLOR-CODING BY RECOVERY STATUS USING YSF COHORTS ##########
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork); library(rlang)
})

target_ysf_vals <- c(5, 10, 15, 20)
sev_keep        <- c("Unburned", "Low", "Moderate", "High")

# Recovery colors (requested)
rec_colors <- c(
  "Recovered"   = "palegreen3",
  "Unrecovered" = "red"
)

# ---- Base data prep (once)
dlr_base <- data_long_complete %>%
  filter(sev_group %in% sev_keep) %>%
  mutate(
    delta_ndvi_min = suppressWarnings(as.numeric(delta_ndvi_min)),
    sev_group      = factor(sev_group, levels = sev_keep)
  ) %>%
  filter(!is.na(delta_ndvi_min))

# ---- Helper: for a given YSF, add cohort + recovery status using rec_ysf{YSF}
add_recovery_for_ysf <- function(df, ysf) {
  rec_col <- paste0("rec_ysf", ysf)
  
  if (!rec_col %in% names(df)) {
    stop("Missing column: ", rec_col, call. = FALSE)
  }
  
  df %>%
    filter(years_since_fire == ysf) %>%
    # cohort filter: keep only pixels with an actual recovery label at this YSF
    filter(!is.na(.data[[rec_col]])) %>%
    mutate(
      recovery_status = dplyr::case_when(
        .data[[rec_col]] == 1 ~ "Recovered",
        .data[[rec_col]] == 0 ~ "Unrecovered",
        TRUE                  ~ NA_character_
      ),
      recovery_status = factor(recovery_status, levels = c("Recovered", "Unrecovered"))
    ) %>%
    filter(!is.na(recovery_status))
}

# ---- Helper to build one scatter for a given predictor and YSF (cohort-filtered)
plot_one_ysf <- function(df, xvar, yvar = "delta_ndvi_min", ysf) {
  
  # cohort-filter + status
  df_y <- add_recovery_for_ysf(df, ysf)
  
  # coerce xvar on the fly
  df_y <- df_y %>%
    mutate("{xvar}" := suppressWarnings(as.numeric(.data[[xvar]]))) %>%
    filter(!is.na(.data[[xvar]]))
  
  n_pts <- nrow(df_y)
  
  # Pearson r within this YSF cohort (overall)
  r_val <- suppressWarnings(cor(df_y[[xvar]], df_y[[yvar]], use = "complete.obs"))
  
  # also report recovered/unrecovered counts
  n_rec <- sum(df_y$recovery_status == "Recovered", na.rm = TRUE)
  n_unr <- sum(df_y$recovery_status == "Unrecovered", na.rm = TRUE)
  
  ggplot(df_y, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_point(aes(color = recovery_status), alpha = 0.25, size = 1) +
    scale_color_manual(values = rec_colors, drop = FALSE) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    labs(
      title    = paste0("YSF = ", ysf, " (cohort-filtered)"),
      subtitle = paste0(
        "N = ", n_pts,
        " (Rec=", n_rec, ", Unrec=", n_unr, ")",
        " | Pearson r = ", ifelse(is.finite(r_val), sprintf("%.3f", r_val), "NA")
      ),
      x = unname(ifelse(xvar %in% names(var_labels), var_labels[[xvar]], xvar)),
      y = "ΔNDVI",
      color = "Recovery status"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title   = element_text(face = "bold"),
      plot.margin  = grid::unit(c(10, 20, 10, 10), "pt"),
      legend.position = "right"
    )
}

# ---- Wrapper: make a 2x2 panel for a single variable across YSFs
panel_for_var <- function(xvar) {
  
  if (!xvar %in% names(dlr_base)) {
    message("Skipping '", xvar, "' (column not found).")
    return(NULL)
  }
  
  plots <- lapply(target_ysf_vals, function(y) {
    
    # if the rec_ysf{y} column doesn't exist, skip that YSF
    rec_col <- paste0("rec_ysf", y)
    if (!rec_col %in% names(dlr_base)) {
      message("Skipping YSF=", y, " for ", xvar, " (missing ", rec_col, ").")
      return(NULL)
    }
    
    plot_one_ysf(dlr_base, xvar, "delta_ndvi_min", y)
  })
  
  plots <- Filter(Negate(is.null), plots)
  
  wrap_plots(plots, ncol = 2) +
    plot_annotation(
      title    = paste0("ΔNDVI vs ", unname(ifelse(xvar %in% names(var_labels), var_labels[[xvar]], xvar))),
      subtitle = paste("Cohort-filtered: only pixels with non-NA rec_ysf{YSF}. Severities included:", paste(sev_keep, collapse = ", "))
    )
}

# ---- Generate all panels
panels <- lapply(vars_to_plot, panel_for_var)
names(panels) <- vars_to_plot

for (v in names(panels)) {
  if (!is.null(panels[[v]])) print(panels[[v]])
}
# ==========================================================

#####  SCATTERPLOTS: COLORING POINTS BY YEAR ####################################
# ---- Function: same scatter but colored by YEAR 
plot_one_ysf_by_year <- function(df, xvar, yvar = "delta_ndvi_min", ysf) {
  df_y <- df %>% filter(years_since_fire == ysf)
  
  x_sym <- rlang::sym(xvar)
  df_y <- df_y %>%
    mutate(
      !!x_sym := suppressWarnings(as.numeric(.data[[xvar]])),
      year    = suppressWarnings(as.numeric(.data[["year"]]))
    ) %>%
    filter(!is.na(.data[[xvar]]), !is.na(year))
  
  n_pts <- nrow(df_y)
  r_val <- suppressWarnings(cor(df_y[[xvar]], df_y[[yvar]], use = "complete.obs"))
  
  ggplot(df_y, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_point(aes(color = year), alpha = 0.25, size = 1) +
    # continuous year scale; swap to scale_color_viridis_c(option = "C") if you prefer
    scale_color_viridis_c() +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1, color = "black") +
    labs(
      title    = paste0("YSF = ", ysf),
      subtitle = paste0("N = ", n_pts, " | Pearson r = ",
                        ifelse(is.finite(r_val), sprintf("%.3f", r_val), "NA")),
      x = unname(ifelse(xvar %in% names(var_labels), var_labels[[xvar]], xvar)),
      y = "ΔNDVI",
      color = "Year"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title   = element_text(face = "bold"),
      plot.margin  = grid::unit(c(10, 20, 10, 10), "pt"),
      legend.position = "right"
    )
}

# Restrict to Unburned only (keep your original dlr_base for severity plots)
dlr_base_unburned <- dlr_base %>% dplyr::filter(sev_group == "Unburned")

# ---- Wrapper: 2x2 panel for a single variable across YSFs, colored by YEAR 
panel_for_var_by_year <- function(xvar) {
  if (!xvar %in% names(dlr_base_unburned)) {
    message("Skipping '", xvar, "' (column not found).")
    return(NULL)
  }
  plots <- lapply(target_ysf_vals, function(y) plot_one_ysf_by_year(dlr_base_unburned, xvar, "delta_ndvi_min", y))
  plots <- Filter(Negate(is.null), plots)
  wrap_plots(plots, ncol = 2) +
    plot_annotation(
      title = paste0("ΔNDVI vs ", unname(ifelse(xvar %in% names(var_labels), var_labels[[xvar]], xvar))),
      subtitle = "Colored by Year — Unburned pixels only"
    )
}

# ---- Generate & print the YEAR-colored panels (unburned only) 
panels_by_year <- lapply(vars_to_plot, panel_for_var_by_year)
names(panels_by_year) <- vars_to_plot

for (v in names(panels_by_year)) {
  if (!is.null(panels_by_year[[v]])) print(panels_by_year[[v]])
}
################################################################################


# ================================================================
### ----- Individual deltaNDVI-covariate Scatterplots ------- ####
### ----- Updated 12/18/25 to use the delta_ndvi_min  ------- ####
# ================================================================

# Set target YSF:
target_ysf <- 5

sev_colors <- c(
  "Unburned" = "grey60",
  "Low"      = "yellow",
  "Moderate" = "orange",
  "High"     = "firebrick"
)

### =====  pptz_JJA  (summer precipitation z-score) ========

# Filter to target YSF and prep variables
dlr_ysf <- data_long_complete %>%
  filter(years_since_fire == target_ysf) %>%
  mutate(
    pptz_JJA   = suppressWarnings(as.numeric(pptz_JJA)),
    delta_ndvi_min = suppressWarnings(as.numeric(delta_ndvi_min))
  ) %>%
  filter(!is.na(pptz_JJA), !is.na(delta_ndvi_min))

# Compute N and Pearson r for subtitle
n_pts <- nrow(dlr_ysf)
r_val <- suppressWarnings(cor(dlr_ysf$pptz_JJA, dlr_ysf$delta_ndvi_min, use = "complete.obs"))

# Plot
p <- ggplot(dlr_ysf, aes(x = pptz_JJA, y = delta_ndvi_min)) +
  geom_point(aes(color = sev_group), alpha = 0.25, size = 1) +
  scale_color_manual(values = sev_colors, drop = FALSE) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  labs(
    title    = paste0("ΔNDVI vs Summer Precip Anomaly — YSF = ", target_ysf),
    subtitle = paste0("N = ", n_pts, " | Pearson r = ",
                      ifelse(is.finite(r_val), sprintf("%.3f", r_val), "NA")),
    x = "JJA Precip Anomaly (z-score)",
    y = "ΔNDVI",
    color = "Severity"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.margin = grid::unit(c(10, 20, 10, 10), "pt")
  )

print(p)

# Replicate plot, but color points by tmeanz_JJA:  -----------------------------
# Filter to target YSF and prep variables
dlr_ysf <- data_long_complete %>%
  filter(years_since_fire == target_ysf) %>%
  mutate(
    pptz_JJA        = suppressWarnings(as.numeric(pptz_JJA)),
    delta_ndvi_min  = suppressWarnings(as.numeric(delta_ndvi_min)),
    tmeanz_JJA      = suppressWarnings(as.numeric(tmeanz_JJA))
  ) %>%
  filter(
    !is.na(pptz_JJA),
    !is.na(delta_ndvi_min),
    !is.na(tmeanz_JJA)
  )

# ---- Bin tmeanz_JJA (quantile-based) 
n_bins <- 5  # adjust as desired

dlr_ysf <- dlr_ysf %>%
  mutate(
    tmean_bin = cut(
      tmeanz_JJA,
      breaks = quantile(tmeanz_JJA, probs = seq(0, 1, length.out = n_bins + 1),
                        na.rm = TRUE),
      include.lowest = TRUE,
      ordered_result = TRUE
    )
  )

# Compute N and Pearson r for subtitle
n_pts <- nrow(dlr_ysf)
r_val <- suppressWarnings(
  cor(dlr_ysf$pptz_JJA, dlr_ysf$delta_ndvi_min, use = "complete.obs")
)

# Plot
p <- ggplot(dlr_ysf, aes(x = pptz_JJA, y = delta_ndvi_min)) +
  geom_point(aes(color = tmean_bin), alpha = 0.25, size = 1) +
  scale_color_viridis_d(
    option = "C",
    name = "JJA Mean Temp\n(z-score)"
  ) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1, color = "black") +
  labs(
    title    = paste0("ΔNDVI vs Summer Precip Anomaly — YSF = ", target_ysf),
    subtitle = paste0(
      "N = ", n_pts,
      " | Pearson r = ",
      ifelse(is.finite(r_val), sprintf("%.3f", r_val), "NA")
    ),
    x = "JJA Precip Anomaly (z-score)",
    y = "ΔNDVI"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title  = element_text(face = "bold"),
    plot.margin = grid::unit(c(10, 20, 10, 10), "pt")
  )

print(p)

# Replicate plot, but color points by swez_Apr:  -------------------------------
# Filter to target YSF and prep variables
dlr_ysf <- data_long_complete %>%
  filter(years_since_fire == target_ysf) %>%
  mutate(
    pptz_JJA        = suppressWarnings(as.numeric(pptz_JJA)),
    delta_ndvi_min  = suppressWarnings(as.numeric(delta_ndvi_min)),
    swez_Apr        = suppressWarnings(as.numeric(swez_Apr))
  ) %>%
  filter(
    !is.na(pptz_JJA),
    !is.na(delta_ndvi_min),
    !is.na(swez_Apr)
  )

# ---- Bin swez_Apr (quantile-based)
n_bins <- 5  # adjust as desired

dlr_ysf <- dlr_ysf %>%
  mutate(
    swe_bin = cut(
      swez_Apr,
      breaks = quantile(swez_Apr,
                        probs = seq(0, 1, length.out = n_bins + 1),
                        na.rm = TRUE),
      include.lowest = TRUE,
      ordered_result = TRUE
    )
  )

# Compute N and Pearson r for subtitle
n_pts <- nrow(dlr_ysf)
r_val <- suppressWarnings(
  cor(dlr_ysf$pptz_JJA, dlr_ysf$delta_ndvi_min, use = "complete.obs")
)

# Plot
p <- ggplot(dlr_ysf, aes(x = pptz_JJA, y = delta_ndvi_min)) +
  geom_point(aes(color = swe_bin), alpha = 0.25, size = 1) +
  scale_color_viridis_d(
    option = "C",
    name = "April SWE\n(z-score)"
  ) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1, color = "black") +
  labs(
    title    = paste0("ΔNDVI vs Summer Precip Anomaly — YSF = ", target_ysf),
    subtitle = paste0(
      "N = ", n_pts,
      " | Pearson r = ",
      ifelse(is.finite(r_val), sprintf("%.3f", r_val), "NA")
    ),
    x = "JJA Precip Anomaly (z-score)",
    y = "ΔNDVI"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title  = element_text(face = "bold"),
    plot.margin = grid::unit(c(10, 20, 10, 10), "pt")
  )

print(p)
#######################################################

### =====  tmeanz_JJA  (summer mean temperature z-score) ========

# Filter to target YSF and prep variables
dlr_ysf <- data_long_complete %>%
  filter(years_since_fire == target_ysf) %>%
  mutate(
    tmeanz_JJA = suppressWarnings(as.numeric(tmeanz_JJA)),
    delta_ndvi_min = suppressWarnings(as.numeric(delta_ndvi_min))
  ) %>%
  filter(!is.na(tmeanz_JJA), !is.na(delta_ndvi_min))

# Compute N and Pearson r for subtitle
n_pts <- nrow(dlr_ysf)
r_val <- suppressWarnings(cor(dlr_ysf$tmeanz_JJA, dlr_ysf$delta_ndvi_min, use = "complete.obs"))

# Plot
p <- ggplot(dlr_ysf, aes(x = tmeanz_JJA, y = delta_ndvi_min)) +
  geom_point(aes(color = sev_group), alpha = 0.25, size = 1) +
  scale_color_manual(values = sev_colors, drop = FALSE) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  labs(
    title    = paste0("ΔNDVI vs Summer Mean Temp Anomaly — YSF = ", target_ysf),
    subtitle = paste0("N = ", n_pts, " | Pearson r = ",
                      ifelse(is.finite(r_val), sprintf("%.3f", r_val), "NA")),
    x = "JJA Mean Temperature (z-score)",
    y = "ΔNDVI",
    color = "Severity"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.margin = grid::unit(c(10, 20, 10, 10), "pt")
  )

print(p)
#######################################################

### =====  tmaxz_JJA  (summer maximum temperature z-score) ========

# Filter to target YSF and prep variables
dlr_ysf <- data_long_complete %>%
  filter(years_since_fire == target_ysf) %>%
  mutate(
    tmaxz_JJA  = suppressWarnings(as.numeric(tmaxz_JJA)),
    delta_ndvi_min = suppressWarnings(as.numeric(delta_ndvi_min))
  ) %>%
  filter(!is.na(tmaxz_JJA), !is.na(delta_ndvi_min))

# Compute N and Pearson r for subtitle
n_pts <- nrow(dlr_ysf)
r_val <- suppressWarnings(cor(dlr_ysf$tmaxz_JJA, dlr_ysf$delta_ndvi_min, use = "complete.obs"))

# Plot
p <- ggplot(dlr_ysf, aes(x = tmaxz_JJA, y = delta_ndvi_min)) +
  geom_point(aes(color = sev_group), alpha = 0.25, size = 1) +
  scale_color_manual(values = sev_colors, drop = FALSE) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  labs(
    title    = paste0("ΔNDVI vs Summer Max Temp Anomaly — YSF = ", target_ysf),
    subtitle = paste0("N = ", n_pts, " | Pearson r = ",
                      ifelse(is.finite(r_val), sprintf("%.3f", r_val), "NA")),
    x = "JJA Maximum Temperature (z-score)",
    y = "ΔNDVI",
    color = "Severity"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.margin = grid::unit(c(10, 20, 10, 10), "pt")
  )

print(p)
#######################################################

### =====  swe_peak_MAM  (peak SWE during Mar–May) ========

# Filter to target YSF and prep variables
dlr_ysf <- data_long_complete %>%
  filter(years_since_fire == target_ysf) %>%
  mutate(
    swe_peak_MAM = suppressWarnings(as.numeric(swe_peak_MAM)),
    delta_ndvi_min   = suppressWarnings(as.numeric(delta_ndvi_min))
  ) %>%
  filter(!is.na(swe_peak_MAM), !is.na(delta_ndvi_min))

# Compute N and Pearson r for subtitle
n_pts <- nrow(dlr_ysf)
r_val <- suppressWarnings(cor(dlr_ysf$swe_peak_MAM, dlr_ysf$delta_ndvi_min, use = "complete.obs"))

# Plot
p <- ggplot(dlr_ysf, aes(x = swe_peak_MAM, y = delta_ndvi_min)) +
  geom_point(aes(color = sev_group), alpha = 0.25, size = 1) +
  scale_color_manual(values = sev_colors, drop = FALSE) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  labs(
    title    = paste0("ΔNDVI vs Peak SWE (MAM) — YSF = ", target_ysf),
    subtitle = paste0("N = ", n_pts, " | Pearson r = ",
                      ifelse(is.finite(r_val), sprintf("%.3f", r_val), "NA")),
    x = "Peak SWE (Mar–May)",
    y = "ΔNDVI",
    color = "Severity"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.margin = grid::unit(c(10, 20, 10, 10), "pt")
  )

print(p)
######################################################

### =====  SWE_Apr  (April SWE) ========

# Filter to target YSF and prep variables
dlr_ysf <- data_long_complete %>%
  filter(years_since_fire == target_ysf) %>%
  mutate(
    SWE_Apr    = suppressWarnings(as.numeric(SWE_Apr)),
    delta_ndvi_min = suppressWarnings(as.numeric(delta_ndvi_min))
  ) %>%
  filter(!is.na(SWE_Apr), !is.na(delta_ndvi_min))

# Compute N and Pearson r for subtitle
n_pts <- nrow(dlr_ysf)
r_val <- suppressWarnings(cor(dlr_ysf$SWE_Apr, dlr_ysf$delta_ndvi_min, use = "complete.obs"))

# Plot
p <- ggplot(dlr_ysf, aes(x = SWE_Apr, y = delta_ndvi_min)) +
  geom_point(aes(color = sev_group), alpha = 0.25, size = 1) +
  scale_color_manual(values = sev_colors, drop = FALSE) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  labs(
    title    = paste0("ΔNDVI vs April SWE — YSF = ", target_ysf),
    subtitle = paste0("N = ", n_pts, " | Pearson r = ",
                      ifelse(is.finite(r_val), sprintf("%.3f", r_val), "NA")),
    x = "April SWE",
    y = "ΔNDVI",
    color = "Severity"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.margin = grid::unit(c(10, 20, 10, 10), "pt")
  )

print(p)
#######################################################

### =====  veg_climate_index_08  (vegetation–climate index, August) ========
# Filter to target YSF and prep variables
dlr_ysf <- data_long_complete %>%
  filter(years_since_fire == target_ysf) %>%
  mutate(
    veg_climate_index_08 = suppressWarnings(as.numeric(veg_climate_index_08)),
    delta_ndvi_min           = suppressWarnings(as.numeric(delta_ndvi_min))
  ) %>%
  filter(!is.na(veg_climate_index_08), !is.na(delta_ndvi_min))

# Compute N and Pearson r for subtitle
n_pts <- nrow(dlr_ysf)
r_val <- suppressWarnings(cor(dlr_ysf$veg_climate_index_08, dlr_ysf$delta_ndvi_min, use = "complete.obs"))

# Plot
p <- ggplot(dlr_ysf, aes(x = veg_climate_index_08, y = delta_ndvi_min)) +
  geom_point(aes(color = sev_group), alpha = 0.25, size = 1) +
  scale_color_manual(values = sev_colors, drop = FALSE) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  labs(
    title    = paste0("ΔNDVI vs Veg–Climate Index — YSF = ", target_ysf),
    subtitle = paste0("N = ", n_pts, " | Pearson r = ",
                      ifelse(is.finite(r_val), sprintf("%.3f", r_val), "NA")),
    x = "Veg–Climate Index (Using August CWD)",
    y = "ΔNDVI",
    color = "Severity"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.margin = grid::unit(c(10, 20, 10, 10), "pt")
  )

print(p)
#######################################################

### =====  cwd_5yr_zscore_08  (5-yr postfire CWD anomaly, August) ========
# Filter to target YSF and prep variables
dlr_ysf <- data_long_complete %>%
  filter(years_since_fire == target_ysf) %>%
  mutate(
    cwd_5yr_zscore_08 = suppressWarnings(as.numeric(cwd_5yr_zscore_08)),
    delta_ndvi_min        = suppressWarnings(as.numeric(delta_ndvi_min))
  ) %>%
  filter(!is.na(cwd_5yr_zscore_08), !is.na(delta_ndvi_min))

# Compute N and Pearson r for subtitle
n_pts <- nrow(dlr_ysf)
r_val <- suppressWarnings(cor(dlr_ysf$cwd_5yr_zscore_08, dlr_ysf$delta_ndvi_min, use = "complete.obs"))

# Plot
p <- ggplot(dlr_ysf, aes(x = cwd_5yr_zscore_08, y = delta_ndvi_min)) +
  geom_point(aes(color = sev_group), alpha = 0.25, size = 1) +
  scale_color_manual(values = sev_colors, drop = FALSE) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  labs(
    title    = paste0("ΔNDVI vs 5-yr Postfire CWD Anomaly (08) — YSF = ", target_ysf),
    subtitle = paste0("N = ", n_pts, " | Pearson r = ",
                      ifelse(is.finite(r_val), sprintf("%.3f", r_val), "NA")),
    x = "CWD Anomaly (5-yr Postfire, August)",
    y = "ΔNDVI",
    color = "Severity"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.margin = grid::unit(c(10, 20, 10, 10), "pt")
  )

print(p)
#######################################################

### =====  twi  (Topographic Wetness Index) ========
# Filter to target YSF and prep variables
dlr_ysf <- data_long_complete %>%
  filter(years_since_fire == target_ysf) %>%
  mutate(
    twi        = suppressWarnings(as.numeric(twi)),
    delta_ndvi_min = suppressWarnings(as.numeric(delta_ndvi_min))
  ) %>%
  filter(!is.na(twi), !is.na(delta_ndvi_min))

# Compute N and Pearson r for subtitle
n_pts <- nrow(dlr_ysf)
r_val <- suppressWarnings(cor(dlr_ysf$twi, dlr_ysf$delta_ndvi_min, use = "complete.obs"))

# Plot
p <- ggplot(dlr_ysf, aes(x = twi, y = delta_ndvi_min)) +
  geom_point(aes(color = sev_group), alpha = 0.25, size = 1) +
  scale_color_manual(values = sev_colors, drop = FALSE) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  labs(
    title    = paste0("ΔNDVI vs Topographic Wetness Index — YSF = ", target_ysf),
    subtitle = paste0("N = ", n_pts, " | Pearson r = ",
                      ifelse(is.finite(r_val), sprintf("%.3f", r_val), "NA")),
    x = "TWI",
    y = "ΔNDVI",
    color = "Severity"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.margin = grid::unit(c(10, 20, 10, 10), "pt")
  )

print(p)
#######################################################

### =====  sev_num  (Fire severity, numeric) ========
# Filter to target YSF and prep variables
dlr_ysf <- data_long_complete %>%
  filter(years_since_fire == target_ysf) %>%
  mutate(
    sev_num    = suppressWarnings(as.numeric(sev_num)),
    delta_ndvi_min = suppressWarnings(as.numeric(delta_ndvi_min))
  ) %>%
  filter(!is.na(sev_num), !is.na(delta_ndvi_min))

# Compute N and Pearson r for subtitle
n_pts <- nrow(dlr_ysf)
r_val <- suppressWarnings(cor(dlr_ysf$sev_num, dlr_ysf$delta_ndvi_min, use = "complete.obs"))

# Plot
p <- ggplot(dlr_ysf, aes(x = sev_num, y = delta_ndvi_min)) +
  geom_point(aes(color = sev_group), alpha = 0.25, size = 1) +
  scale_color_manual(values = sev_colors, drop = FALSE) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  labs(
    title    = paste0("ΔNDVI vs Fire Severity (numeric) — YSF = ", target_ysf),
    subtitle = paste0("N = ", n_pts, " | Pearson r = ",
                      ifelse(is.finite(r_val), sprintf("%.3f", r_val), "NA")),
    x = "Severity (numeric)",
    y = "ΔNDVI",
    color = "Severity"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.margin = grid::unit(c(10, 20, 10, 10), "pt")
  )

print(p)
#######################################################

#### ====== deltaNDVI vs All Predictors  ===== 
# ---- load necessary libraries 
library(dplyr)
library(tidyr)
library(purrr)
library(forcats)
library(ggplot2)

# ---- predictors to include (same order as in family_map) 
predictors <- c(
  "ppttot_JJA", "pptz_JJA",                   # Precip
  "tmeanz_JJA", "tmaxz_JJA", "tmeanavg_JJA", "tmaxmean_JJA",  # Temp
  "swe_peak_MAM", "SWE_Apr",                  # SWE
  "aetavg_JJA", "veg_climate_index_08", "cwd_5yr_zscore_08",
  "twi", "sev_num" ,                           # Other
  # ,"fire_size_ha"
  "aet_30yrAvg_TT", "cwd_30yrAvg_TT",
  "elevation_m"
)

# ---- Define predictor families (counts must match length(predictors)) 
family_map <- tibble(
  predictor = predictors,
  family = c(
    rep("Precip", 2),
    rep("Temp",   4),
    rep("SWE",    2),
    rep("Other",  8)
  )
)

# ---- Filter to target YSF and coerce to numeric to mirror individual plots 
df_ysf <- data_long_complete %>%
  filter(years_since_fire == target_ysf) %>%
  st_drop_geometry() %>%
  mutate(
    delta_ndvi_min = suppressWarnings(as.numeric(delta_ndvi_min)),
    across(all_of(predictors), ~ suppressWarnings(as.numeric(.)))
  )

# For subtitle: number of pixel–year rows after YSF filtering (NOT pivoted)
n_pts <- nrow(df_ysf)

# ---- Build plot_long using the SAME per-predictor NA filter as individual plots 
# For each predictor, keep only rows where that predictor AND delta_ndvi_min are non-NA
plot_long <- map_dfr(predictors, function(var) {
  df_ysf %>%
    select(delta_ndvi_min, all_of(var)) %>%
    filter(!is.na(.data[[var]]), !is.na(delta_ndvi_min)) %>%
    transmute(predictor = var, value = .data[[var]], delta_ndvi_min)
})

# ---- Min–max scale within predictor for comparable x-axis (does not affect r) 
plot_long <- plot_long %>%
  group_by(predictor) %>%
  mutate(
    minv = min(value, na.rm = TRUE),
    maxv = max(value, na.rm = TRUE),
    denom = maxv - minv,
    value_01 = ifelse(denom > 0, (value - minv) / denom, NA_real_)
  ) %>%
  ungroup() %>%
  filter(!is.na(value_01))  # drop predictors that are constant at this YSF

# ---- Compute Pearson r ON THE SAME ROWS used for scatterplots 
corr_tbl <- plot_long %>%
  group_by(predictor) %>%
  summarise(
    r     = suppressWarnings(cor(value, delta_ndvi_min, use = "complete.obs")),
    abs_r = abs(r),
    n_used = n(),  # rows used for that predictor
    .groups = "drop"
  ) %>%
  left_join(family_map, by = "predictor")

# ---- Order legend: by family, then descending |r| within family 
family_levels <- c("SWE", "Precip", "Temp", "Other")
corr_tbl <- corr_tbl %>%
  mutate(
    family = factor(family, levels = family_levels),
    predictor = fct_reorder2(predictor, family, abs_r, .desc = TRUE)
  )

# ---- Legend labels: "var (r = 0.123)" 
label_map <- setNames(
  paste0(corr_tbl$predictor, " (r = ", ifelse(is.finite(corr_tbl$r), sprintf("%.3f", corr_tbl$r), "NA"), ")"),
  corr_tbl$predictor
)

# ---- Use the ordered predictor factor in plotting data 
plot_long <- plot_long %>%
  mutate(predictor = factor(predictor, levels = levels(corr_tbl$predictor)))

# ---- Plot 
ggplot(plot_long, aes(x = value_01, y = delta_ndvi_min,
                      color = predictor, linetype = predictor)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_color_discrete(labels = function(x) label_map[x]) +
  scale_linetype_discrete(labels = function(x) label_map[x]) +
  labs(
    title    = paste0("ΔNDVI vs Multiple Predictors — Trendlines (YSF = ", target_ysf, ")"),
    subtitle = paste0("N = ", n_pts,
                      " pixel-year rows | Predictors scaled 0–1"),
    x = "Predictor value (scaled)",
    y = "ΔNDVI",
    color   = "Pearson's r",
    linetype = "Pearson's r"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.margin = grid::unit(c(10, 20, 10, 10), "pt"),
    legend.key.height = unit(12, "pt")
  )

################################################################################



################################################################################
###    Check multicollinearity in final covariates     #########################
################################################################################

####### ===== ORIGINAL CODE =========== ########################################
# --- Pearson correlations: delta_ndvi_min vs climate covariates   #######
stopifnot("delta_ndvi_min" %in% names(data_long_complete))

# Candidate covariates (long-format names; one value per pixel-year combo (so, all data))
cand_vars <- c(
  "ppttot_JJA", 
  "pptz_JJA",
  "tmaxmean_JJA", 
  "tmaxz_JJA", 
  "tmeanavg_JJA", 
  "tmeanz_JJA",
  #"swe_peak_MAM", 
  "SWE_Apr", 
  "swez_Apr"
)

# Keep only candidates that actually exist in the data
cand_vars <- cand_vars[cand_vars %in% names(data_long_complete)]
if (length(cand_vars) == 0) stop("No candidate covariates found in the selected dataset.")

# Function to compute Pearson r, p-value, CI, and N (pairwise complete)
corr_one <- function(var) {
  d <- data_long_complete %>% 
    dplyr::select(delta_ndvi_min, !!rlang::sym(var)) %>% tidyr::drop_na()
  n <- nrow(d)
  if (n < 3) {
    tibble::tibble(variable = var, n = n, r = NA_real_, p = NA_real_, ci_low = NA_real_, ci_high = NA_real_)
  } else {
    ct <- suppressWarnings(cor.test(d$delta_ndvi_min, d[[var]], method = "pearson"))
    tibble::tibble(
      variable = var,
      n        = n,
      r        = unname(ct$estimate),
      p        = ct$p.value,
      ci_low   = unname(ct$conf.int[1]),
      ci_high  = unname(ct$conf.int[2])
    )
  }
}

# Run and sort by absolute correlation magnitude
pearson_results <- purrr::map_dfr(cand_vars, corr_one) %>%
  dplyr::arrange(dplyr::desc(abs(r)))

print(pearson_results, n = Inf)

# --- Prep for plotting (add family + significance; order by |r|) ---
plot_df <- pearson_results %>%
  dplyr::mutate(
    family = dplyr::case_when(
      grepl("^ppt", variable)              ~ "Precip",
      grepl("^tmax|^tmean", variable)      ~ "Temp",
      grepl("^swe|^SWE", variable)         ~ "SWE",
      grepl("^aet", variable)              ~ "AET",
      TRUE                                 ~ "Other"
    ),
    sig = dplyr::case_when(
      is.na(p)         ~ "p NA",
      p < 0.001        ~ "p < 0.001",
      p < 0.01         ~ "p < 0.01",
      p < 0.05         ~ "p < 0.05",
      TRUE             ~ "n.s."
    )
  ) %>%
  dplyr::arrange(dplyr::desc(abs(r))) %>%
  dplyr::mutate(variable = factor(variable, levels = rev(variable)))

# --- Lollipop / pointrange plot of r with 95% CI ---
ggplot(plot_df, aes(x = variable, y = r, color = family)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymin = ci_low, ymax = ci_high),
                  position = position_dodge(width = 0.4)) +
  coord_flip() +
  geom_text(
    aes(label = sprintf("r = %.3f", r), y = r),
    hjust = 1.05, size = 3, color = "black"
  ) +
  geom_point(
    aes(shape = sig),
    size = 2.5, fill = NA, stroke = 1
  ) +
  labs(
    title = "Pearson correlations with ΔNDVI (all pixel-years)",
    subtitle = "Point locations show r; bars show 95% CI",
    x = NULL, y = "Pearson r"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.title = element_blank())

### ---- Correlation heatmaps: ΔNDVI + climate covariates  ---- #####
vars <- c(
  #"delta_ndvi_min",
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
  "sev_num",
  "twi",
  #"veg_climate_index",
  #"cwd_5yr_postfire_avg_08",
  "cwd_5yr_zscore_08",
  #"aet_30yrAvg_TT",
  "cwd_30yrAvg_TT"
  #"elevation_m"
)

compute_corr_df <- function(d, vars) {
  vars <- vars[vars %in% names(d)]
  cmb <- expand.grid(x = vars, y = vars, stringsAsFactors = FALSE)
  cmb$r <- purrr::pmap_dbl(
    list(cmb$x, cmb$y),
    ~ suppressWarnings(cor(d[[..1]], d[[..2]], use = "pairwise.complete.obs", method = "pearson"))
  )
  cmb$n_pair <- purrr::pmap_int(
    list(cmb$x, cmb$y),
    ~ sum(stats::complete.cases(d[[..1]], d[[..2]]))
  )
  cmb$x <- factor(cmb$x, levels = vars)
  cmb$y <- factor(cmb$y, levels = vars)
  cmb
}

## Heatmap plotting function:
plot_corr_heatmap <- function(corr_df, title, subtitle, n_rows) {
  
  # Ensure consistent factor ordering
  vars <- unique(corr_df$x)
  
  corr_df <- corr_df %>%
    dplyr::mutate(
      x = factor(x, levels = vars),
      y = factor(y, levels = vars)
    ) %>%
    # keep only lower triangle (drop self-comparisons too)
    dplyr::filter(as.integer(y) > as.integer(x))
  
  ggplot(corr_df, aes(x = x, y = y, fill = r)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", r)),
              size = 3,
              family = "Times New Roman") +
    scale_fill_gradient2(
      limits   = c(-1, 1),
      low      = "blue",
      mid      = "white",
      high     = "red",
      midpoint = 0,
      oob      = scales::squish
    ) +
    coord_equal() +
    labs(
      title    = title,
      subtitle = paste0(subtitle, "\ n = ", format(n_rows, big.mark = ",")),
      x = NULL,
      y = NULL,
      fill = "r"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(),
      legend.title = element_text(),
      legend.text = element_text(),
      panel.grid = element_blank()
    )
}

# Use the same ordered variable list for all heatmaps
vars_order <- c(
  #"delta_ndvi_min",
  #"ppttot_JJA",
  "pptz_JJA",
  #"tmaxmean_JJA",
  #"tmaxz_JJA",
  #"tmeanavg_JJA",
  "tmeanz_JJA",
  #"swe_peak_MAM",
  #"SWE_Apr",
  "swez_Apr",
  "sev_num",
  "twi",
  #"veg_climate_index_08",          
  #"cwd_5yr_postfire_avg_08",  ## note this is NOT the z-score
  "cwd_5yr_zscore_08",
  #"aet_30yrAvg_TT",
  "cwd_30yrAvg_TT"
  #"elevation_m"
)

##### ==== Print heat maps:  ==================
# --- Heat map 1: ALL post-fire pixel–years 
d_all   <- data_long_complete
n_all   <- nrow(d_all)
corr_all <- compute_corr_df(d_all, vars_order)
p_all <- plot_corr_heatmap(
  corr_all,
  "Heat Map: Pearson's r Pairwise Correlations",  # Title
  "All Pixel-Years;",   # Subtible
  n_all
)
print(p_all)

# --- Heat map 2: YSF == 5
d_ysf5   <- data_long_complete %>% dplyr::filter(years_since_fire == 5)
n_ysf5   <- nrow(d_ysf5)
corr_ysf5 <- compute_corr_df(d_ysf5, vars_order)
p_ysf5 <- plot_corr_heatmap(
  corr_ysf5,
  "Heat Map: Pearson's r Pairwise Correlations",
  "YSF = 5;",
  n_ysf5
)
print(p_ysf5)

# --- Heat map 3: YSF == 10
d_ysf10   <- data_long_complete %>% dplyr::filter(years_since_fire == 10)
n_ysf10   <- nrow(d_ysf10)
corr_ysf10 <- compute_corr_df(d_ysf10, vars_order)
p_ysf10 <- plot_corr_heatmap(
  corr_ysf10,
  "Heat Map: Pearson's r Pairwise Correlations",
  "YSF = 10;",
  n_ysf10
)
print(p_ysf10)

# --- Heat map 4: YSF == 15
d_ysf15   <- data_long_complete %>% dplyr::filter(years_since_fire == 15)
n_ysf15   <- nrow(d_ysf15)
corr_ysf15 <- compute_corr_df(d_ysf15, vars_order)
p_ysf15 <- plot_corr_heatmap(
  corr_ysf15,
  "Heat Map: Pearson's r Pairwise Correlations",
  "YSF = 15;",
  n_ysf15
)
print(p_ysf15)

# --- Heat map 5: YSF == 20
d_ysf20   <- data_long_complete %>% dplyr::filter(years_since_fire == 20)
n_ysf20   <- nrow(d_ysf20)
corr_ysf20 <- compute_corr_df(d_ysf20, vars_order)
p_ysf20 <- plot_corr_heatmap(
  corr_ysf20,
  "Heat Map: Pearson's r Pairwise Correlations",
  "YSF = 20;",
  n_ysf20
)
print(p_ysf20)
################################################################################

####### ===== NEW CODE (FOR-LOOP) =========== ##################################
# --- PART 1 ---Extract and prepare all potential covariates ------
ysf_values <- c(5, 10, 15, 20)

covariates_df <- data_long_complete %>%
  st_drop_geometry() %>%
  select(
    # Response variable
    #delta_ndvi_min,
    
    # Originals - topography, VCI, fire severity
    twi, 
    #veg_climate_index_08, 
    sev_num, 
    #fire_size_ha,
    
    # CWD measures 
    #cwd_3yr_zscore_08,
    #cwd_5yr_postfire_avg_08,
    cwd_5yr_zscore_08, 
    #cwd_5yr_prefire_zscore_08,
    
    # June, July, and August precip measures
    #ppttot_JJA,  # total summer precip (JJA)
    pptz_JJA, # z-score of average precip across JJA
    
    # Avg and maximum temps: summer (june, july Aug)
    #tmaxmean_JJA, # Mean of daily maximum temperature across JJA (°C). 
    #tmaxz_JJA,    # Standardized anomaly (z-score) of JJA mean daily max temperature.
    #tmeanavg_JJA, # Mean of daily mean temperature across JJA (°C).
    tmeanz_JJA,   # Standardized anomaly (z-score) of JJA mean daily mean temperature.
    
    # Actual Evapotranspiration
    #aetavg_JJA,
    
    # Snow Water Equivalent
    #swe_peak_MAM, 
    #SWE_Apr,
    swez_Apr,
    
    # Elevation
    #elevation_m,
    
    # 30-year Normals
    #aet_30yrAvg_TT,
    cwd_30yrAvg_TT
  )

# Ensure sev_num is numeric (and not accidentally a factor/character)
if ("sev_num" %in% names(covariates_df)) {
  if (!is.numeric(covariates_df$sev_num)) {
    covariates_df$sev_num <- suppressWarnings(as.numeric(covariates_df$sev_num))
  }
}

# --- PART 2 ----- Compute correlation matrix with consistent variable ordering ----
var_order <- names(covariates_df)  # preserve the select() order
cor_matrix <- cor(covariates_df[, var_order], use = "pairwise.complete.obs")

# Count complete pixel-year rows used in correlation
n_obs <- nrow(na.omit(covariates_df))

### Original code using corrplot: 
# Visualize correlations in a heat map
#corrplot(cor_matrix,
#         method = "color",
#         type = "upper",       # upper triangle only
#         order = "original",   # preserve input order
#         addCoef.col = "black",
#         tl.col = "black",
#         tl.cex = 1,
#         number.cex = 0.9,
#         diag = FALSE)         # hide 1's on the diagonal

# Add title and subtitle (with row count)
#title(main = "Correlation Matrix of Candidate Covariates",
#      sub  = paste("Based on", n_obs, "pixel-years"))

# ---------------- Variable labels 
var_labels <- c(
  "sev_num"              = "Fire Severity Class",
  "twi"                  = "Topographic Wetness Index (TWI)",
  "pptz_JJA"             = "Summer Precipitation",
  "tmaxz_JJA"            = "Summer Maximum Temperature",
  "tmeanz_JJA"           = "Summer Mean Temperature",
  "swez_Apr"             = "Snow Water Equivalent, April",
  "cwd_5yr_zscore_08"    = "5-Year Post-Fire Aug CWD",
  "cwd_30yrAvg_TT"       = "CWD 30-Year Normal"
)

## Convert correlation matrix to long format
cor_long <- as.data.frame(cor_matrix) %>%
  tibble::rownames_to_column("var1") %>%
  tidyr::pivot_longer(
    -var1,
    names_to  = "var2",
    values_to = "Pearson_r"   # <- keep syntactic
  ) %>%
  dplyr::mutate(
    var1 = factor(var1, levels = var_order),
    var2 = factor(var2, levels = var_order),
    i    = as.integer(var1),
    j    = as.integer(var2),
    YSF  = ysf
  ) %>%
  dplyr::filter(i < j)

corr_table_all[[as.character(ysf)]] <- cor_long

## Generate heatmap
p <- ggplot(cor_long, aes(x = var2, y = var1, fill = Pearson_r)) +
  geom_tile() +
  geom_text(
    aes(label = sprintf("%.2f", Pearson_r)),
    size = 4,
    family = "Times New Roman"
  ) +
  scale_fill_gradient2(
    limits   = c(-1, 1),
    low      = "blue",
    mid      = "white",
    high     = "red",
    midpoint = 0,
    oob      = scales::squish
  ) +
  coord_equal() +
  scale_x_discrete(labels = var_labels) +
  scale_y_discrete(position = "right", labels = var_labels) +
  labs(
    title    = "Correlation Matrix of Candidate Covariates",
    subtitle = paste("YSF =", ysf, "| n =", n_obs),
    x = NULL,
    y = NULL,
    fill = "Pearson's r"   # <- legend title shows exactly what you want
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    plot.subtitle = element_text(hjust = 0.5, size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(hjust = 0, size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, margin = ggplot2::margin(b = 6)),
    panel.grid = element_blank()
  )

print(p)
## List highly correlated pairs ---

# Function to extract unique var pairs meeting an abs(corr) threshold
get_corr_pairs <- function(cor_mat, threshold = 0.3) {
  m <- cor_mat
  # keep only upper triangle (no duplicates) and drop diagonal
  m[lower.tri(m, diag = TRUE)] <- NA
  
  df <- as.data.frame(as.table(m))
  names(df) <- c("var1", "var2", "cor")
  df <- df[!is.na(df$cor), , drop = FALSE]
  df$abs_cor <- abs(df$cor)
  
  # filter and sort
  df <- df[df$abs_cor >= threshold, , drop = FALSE]
  df <- df[order(-df$abs_cor), , drop = FALSE]
  rownames(df) <- NULL
  as_tibble(df)
}

# Tables at two thresholds
corr_pairs_0_30 <- get_corr_pairs(cor_matrix, threshold = 0.30)
corr_pairs_0_70 <- get_corr_pairs(cor_matrix, threshold = 0.70)

# Inspect in console
#cat("\nPairs with |correlation| >= 0.30 (n =", nrow(corr_pairs_0_30), ")\n")
#print(corr_pairs_0_30)

#cat("\nPairs with |correlation| >= 0.70 (n =", nrow(corr_pairs_0_70), ")\n")
print(corr_pairs_0_70)

## --- PART 3 -- Check Variance Inflation Factors (VIF) -------
# Drop columns with all NA or zero variance
nzv_cols <- names(covariates_df)[
  vapply(covariates_df, function(x) {
    ux <- unique(x[!is.na(x)])
    length(ux) <= 1
  }, logical(1))
]
if (length(nzv_cols)) {
  message("Dropping all-NA / zero-variance columns: ", paste(nzv_cols, collapse = ", "))
  covariates_df <- covariates_df %>% select(-all_of(nzv_cols))
}

set.seed(123)
covariates_df$dummy_response <- rnorm(nrow(covariates_df))

# Build formula from all predictors 
predictor_vars <- setdiff(names(covariates_df), "dummy_response")
vif_formula <- as.formula(paste("dummy_response ~", paste(predictor_vars, collapse = " + ")))

# Fit linear model
dummy_fit <- lm(vif_formula, data = covariates_df, na.action = na.omit)

# If there are aliased (NA) coefficients, refit using only estimable terms
na_coefs <- names(coef(dummy_fit))[is.na(coef(dummy_fit))]
if (length(na_coefs)) {
  keep_terms <- setdiff(names(coef(dummy_fit))[!is.na(coef(dummy_fit))], "(Intercept)")
  message("Refitting without aliased terms: ", paste(na_coefs, collapse = ", "))
  dummy_fit <- lm(reformulate(keep_terms, response = "dummy_response"),
                  data = covariates_df, na.action = na.omit)
}

# Compute VIFs (all terms numeric now -> standard VIFs)
raw_vif <- car::vif(dummy_fit)

vif_table <- data.frame(
  Term = names(raw_vif),
  VIF  = as.numeric(raw_vif),
  row.names = NULL,
  check.names = FALSE
) %>% arrange(desc(VIF))

print(vif_table)
message("\nInterpretation: VIF < 5 is generally acceptable (10 is a lenient cutoff).")

# flag high-collinearity terms
threshold <- 5
high <- vif_table %>% filter(VIF >= threshold)
if (nrow(high)) {
  message("\nTerms with VIF ≥ ", threshold, ": ", paste(high$Term, collapse = ", "))
}

## --- PART 4 -- multicollinearity for-loop for each target YSF -------
corr_table_all <- list()
vif_table_all  <- list()

### Variable labels

var_labels <- c(
  "sev_num"           = "Fire Severity",
  "twi"               = "TWI",
  "pptz_JJA"          = "Precipitation",
  "tmeanz_JJA"        = "Temperature",
  "swez_Apr"          = "Snowpack",
  "cwd_5yr_zscore_08" = "5-Year CWD",
  "cwd_30yrAvg_TT"    = "CWD 30-Year Normal"
)

### Loop through YSF values
for (ysf in ysf_values) {
  
  message("Processing YSF = ", ysf)
  
  covariates_df <- data_long_complete %>%
    dplyr::filter(years_since_fire == ysf) %>%
    sf::st_drop_geometry() %>%
    dplyr::select(
      twi,
      sev_num,
      cwd_5yr_zscore_08,
      pptz_JJA,
      tmeanz_JJA,
      swez_Apr,
      cwd_30yrAvg_TT
    )
  
  if (!is.numeric(covariates_df$sev_num)) {
    covariates_df$sev_num <- suppressWarnings(as.numeric(covariates_df$sev_num))
  }
  
  ## Compute correlation matrix
  var_order <- names(covariates_df)
  
  cor_matrix <- cor(
    covariates_df[, var_order],
    use = "pairwise.complete.obs"
  )
  
  n_obs <- nrow(na.omit(covariates_df))
  
  ## Convert correlation matrix to long format
  cor_long <- as.data.frame(cor_matrix) %>%
    tibble::rownames_to_column("var1") %>%
    tidyr::pivot_longer(
      -var1,
      names_to = "var2",
      values_to = "Pearson's r"
    ) %>%
    dplyr::mutate(
      var1 = factor(var1, levels = var_order),
      var2 = factor(var2, levels = var_order),
      i = as.integer(var1),
      j = as.integer(var2),
      YSF = ysf
    ) %>%
    dplyr::filter(i < j)
  
  corr_table_all[[as.character(ysf)]] <- cor_long
  
  ## Generate heatmap
  print(
    ggplot(cor_long, aes(x = var2, y = var1, fill = `Pearson's r`)) +
      geom_tile() +
      geom_text(   ## section for text inside boxes.
        aes(label = sprintf("%.2f", `Pearson's r`)),
        size = 5,
        fontface = "bold",
        family = "Times New Roman"
      ) +
      scale_fill_gradient2(
        limits   = c(-1, 1),
        low      = "blue",
        mid      = "white",
        high     = "darkred",
        midpoint = 0,
        oob      = scales::squish
      ) +
      coord_equal() + 
      scale_x_discrete(labels = var_labels) +
      scale_y_discrete(position = "right", labels = var_labels) +
      labs(
        title = "Covariate Heat Map",
        subtitle = paste("YSF =", ysf, "| n = ", n_obs),
        x = NULL,
        y = NULL
      ) +
      theme_minimal() +
      theme(
        text = element_text(family = "Times New Roman"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        plot.subtitle = element_text(hjust = 0.5, size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(hjust = 0, size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, margin = ggplot2::margin(b = 12)),
        panel.grid = element_blank()
      )
  )

  ## Check Variance Inflation Factors
  nzv_cols <- names(covariates_df)[
    vapply(covariates_df, function(x) {
      ux <- unique(x[!is.na(x)])
      length(ux) <= 1
    }, logical(1))
  ]
  
  if (length(nzv_cols)) {
    covariates_df <- covariates_df %>% dplyr::select(-all_of(nzv_cols))
  }
  
  set.seed(123)
  covariates_df$dummy_response <- rnorm(nrow(covariates_df))
  
  predictor_vars <- setdiff(names(covariates_df), "dummy_response")
  
  vif_formula <- as.formula(
    paste("dummy_response ~", paste(predictor_vars, collapse = " + "))
  )
  
  dummy_fit <- lm(vif_formula, data = covariates_df, na.action = na.omit)
  
  raw_vif <- car::vif(dummy_fit)
  
  vif_table <- data.frame(
    Term = names(raw_vif),
    VIF  = as.numeric(raw_vif),
    YSF  = ysf,
    row.names = NULL,
    check.names = FALSE
  )
  
  vif_table_all[[as.character(ysf)]] <- vif_table
}

### ==== Combine all 4 heat maps into one panel =========
library(patchwork)

# Order for display
var_display_order <- c(
  "TWI",
  "Fire Severity",
  "5-Year CWD",
  "Precipitation",
  "Temperature",
  "Snowpack",
  "CWD 30-Year Normal"
)

# Prepare plotting table
pearson_plot_df <- pearson_table %>%
  dplyr::mutate(
    var1 = factor(var1, levels = var_display_order),
    var2 = factor(var2, levels = var_display_order),
    YSF  = factor(YSF, levels = c(5, 10, 15, 20))
  )

# Function to make one panel
make_corr_panel <- function(target_ysf) {
  
  panel_df <- pearson_plot_df %>%
    dplyr::filter(YSF == target_ysf)
  
  ggplot(panel_df, aes(x = var2, y = var1, fill = `Pearson's r`)) +
    geom_tile() +
    geom_text(
      aes(label = sprintf("%.2f", `Pearson's r`)),
      size = 5,
      fontface = "bold",
      family = "Times New Roman"
    ) +
    scale_fill_gradient2(
      limits   = c(-1, 1),
      low      = "blue",
      mid      = "white",
      high     = "darkred",
      midpoint = 0,
      oob      = scales::squish,
      name     = "Pearson's r"
    ) +
    coord_equal() +
    labs(
      title = paste("YSF =", target_ysf),
      x = NULL,
      y = NULL
    ) +
    theme_minimal() +
    theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
      axis.text.y = element_text(size = 13),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, margin = ggplot2::margin(b = 10)),
      panel.grid = element_blank()
    )
}

# Build the four panels
p5  <- make_corr_panel(5)
p10 <- make_corr_panel(10)
p15 <- make_corr_panel(15)
p20 <- make_corr_panel(20)

# Combine with panel labels a, b, c, d
combined_corr_panel <- (p5 + p10) / (p15 + p20) +
  plot_annotation(tag_levels = "a")

# Force tag styling
combined_corr_panel <- combined_corr_panel &
  theme(
    plot.tag = element_text(
      size = 20,           # adjust size here
      face = "bold",
      family = "Times New Roman"
    )
  )

print(combined_corr_panel)
# -------------------------------------------------------------------------

### Combine Pearson correlation results ######
pearson_table <- dplyr::bind_rows(corr_table_all) %>%
  dplyr::select(YSF, var1, var2, `Pearson's r`) %>%
  dplyr::mutate(
    var1 = var_labels[as.character(var1)],
    var2 = var_labels[as.character(var2)],
    `Pearson's r` = round(`Pearson's r`, 2)
  )
print(pearson_table)

### Combine VIF results #########
vif_table_final <- dplyr::bind_rows(vif_table_all) %>%
  dplyr::mutate(
    Term = var_labels[Term],
    VIF  = round(VIF, 2)
  ) %>%
  dplyr::select(YSF, Term, VIF) %>%
  dplyr::arrange(YSF, dplyr::desc(VIF))
print(vif_table_final)

message("\nInterpretation: VIF < 5 is generally acceptable (10 is a lenient cutoff).")
  
################################################################################


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
################################################################################
#          -------------------------------------------                         #
#########     RESEARCH QUESTION 3: ENVIRONMENTAL DRIVERS     ###################
#          -------------------------------------------                         #
################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


#==============================================================
  # MASTER SWITCHES (edit HERE ONLY.)
# ==============================================================
ysf_set <- c(5, 10, 15, 20)                 # <<< choose Years Since Fire to process
cwd_zscore_var <- "cwd_5yr_zscore_08"   # e.g., "cwd_5yr_zscore_08" or "cwd_3yr_zscore"

# Random Forest model / Cross-validation controls
k      <- 4
ntree  <- 500
seed   <- 460                       # global seed

# Fold-balancing weights (cluster -> fold assignment)
w_size     <- 0.35   # weight for balancing fold sample sizes
w_sev      <- 0.65   # weight for balancing sev_group proportions
max_dev    <- 0.05  # max allowed absolute deviation per class from global prop
big_pen    <- 700   # penalty once max_dev exceeded
n_restarts <- 500   # multi-start attempts for cluster->fold assignment


#          -------------------------------------------                         #
#########     ENVIRONMENTAL BLOCKING / CLUSTERING     ##########################
#          -------------------------------------------                         #
# ==============================================================
# CLUSTERING CONTROLS — silhouette-driven with adaptive scaling
# ==============================================================

# k-means cluster health guards
sil_min               <- 0.30  # mininum silhouette value for clusters (moderate structure)
sil_tol               <- 0.10
min_clusters_per_fold <- 6
min_cluster_size      <- 100
k_max_mult            <- 15  # Higher value = tests more k → slower processing
strict_enforce_kmeans <- FALSE # set to TRUE to stop processing if sil isn't met

# Enable adaptive scaling of k-means cluster size/range by sample size
use_dynamic_clustering <- TRUE

# Function: scale cluster parameters based on available sample size (n)
get_cluster_params <- function(n_pts, k_folds = k) {
  # Adaptive heuristics for smaller sample sets
  if (n_pts < 500) {
    min_cluster_size <- max(20, floor(n_pts / 20))
    min_clusters_per_fold <- 2
    k_max_mult <- 5
  } else if (n_pts < 1500) {
    min_cluster_size <- max(30, floor(n_pts / 40))
    min_clusters_per_fold <- 3
    k_max_mult <- 8
  } else if (n_pts < 3000) {
    min_cluster_size <- max(40, floor(n_pts / 60))
    min_clusters_per_fold <- 4
    k_max_mult <- 10
  } else {
    min_cluster_size <- max(50, floor(n_pts / 100))
    min_clusters_per_fold <- 5
    k_max_mult <- 15
  }
  
  message(sprintf(
    "Adaptive cluster settings: n=%d | min_cluster_size=%d | min_clusters_per_fold=%d | k_max_mult=%d",
    n_pts, min_cluster_size, min_clusters_per_fold, k_max_mult
  ))
  
  list(
    min_clusters_per_fold = min_clusters_per_fold,
    min_cluster_size = min_cluster_size,
    k_max_mult = k_max_mult
  )
}

# ==============================================================
# Functions
# ==============================================================

# Canonical ordering for`sev_group`
sev_levels <- c("Unburned","Low","Moderate","High")

# Count-align helper: ensure all sev_levels are present in a count named vector
align_counts <- function(tab, sev_levels) {
  out <- setNames(integer(length(sev_levels)), sev_levels)
  if (length(tab)) out[names(tab)] <- as.integer(tab)
  out[is.na(out)] <- 0L
  out
}

# New, 2/9/26: PCA variance % Table ----------------------
get_pca_variance_table <- function(pc_obj, ysf, max_pcs = 10) {
  
  # eigenvalues
  eig <- pc_obj$sdev^2
  
  var_explained <- eig / sum(eig)
  cum_var <- cumsum(var_explained)
  
  tibble::tibble(
    YSF = ysf,
    PC = paste0("PC", seq_along(eig)),
    var_explained = var_explained,
    pct_explained = 100 * var_explained,
    cum_var_explained = cum_var,
    cum_pct_explained = 100 * cum_var
  ) %>%
    dplyr::slice_head(n = max_pcs)
}
# ------------------ Dynamic k-means cluster selector ------------------
# Chooses k within safe bounds using a minimum silhouette score
select_kmeans_k <- function( 
    X,
    k_folds,
    min_clusters_per_fold = min_clusters_per_fold,
    min_cluster_size      = min_cluster_size,
    k_max_mult            = k_max_mult,
    seed                  = seed,
    sil_min               = sil_min,   
    sil_tol               = sil_tol,   
    strict_enforce        = strict_enforce_kmeans,    # if TRUE: stop when no k meets sil_min; if FALSE: warn & pick best
    
    # Stabilizers / performance knobs
    nstart                = 500,
    iter.max              = 4000,
    algorithm             = "Lloyd",
    
    # Cluster-size tolerance
    size_tolerance_frac   = 0.10,   # allow clusters down to (1 - frac) * min_cluster_size
    undersize_allow_frac  = 0.10,   # allow up to this fraction of clusters to be undersized
    
    # One-shot automatic relaxation if no candidates survive
    relax_on_fail         = TRUE,
    relax_factor          = 0.50,
    
    # Silhouette sampling (avoid full dist() on big n)
    sil_sample_cap        = 4500     # if n > cap, compute silhouette on a random subset
){
  set.seed(seed)
  n <- nrow(X)
  
  # Trivial case: fewer rows than folds
  if (n < k_folds) {
    km <- stats::kmeans(X, centers = n, nstart = 20, iter.max = 1000, algorithm = algorithm)
    return(list(k = n, km = km, silhouette = NA_real_, sizes = rep(1L, n), met_threshold = FALSE))
  }
  
  # ----- Bounds for k -----
  k_min <- max(k_folds * min_clusters_per_fold, k_folds)
  k_max <- min(
    n,
    max(k_min, floor(n / max(1, min_cluster_size))),
    k_folds * k_max_mult
  )
  if (k_max < k_min) k_max <- k_min
  ks <- seq(k_min, k_max)
  
  # ----- (Approx) silhouette setup -----
  use_sampled_sil <- n > sil_sample_cap
  sample_idx <- if (use_sampled_sil) sample.int(n, sil_sample_cap) else seq_len(n)
  D <- tryCatch(stats::dist(X[sample_idx, , drop = FALSE]), error = function(e) NULL)
  
  mean_sil <- function(labels) {
    if (is.null(D)) return(NA_real_)
    labs <- if (use_sampled_sil) labels[sample_idx] else labels
    if (length(unique(labs)) < 2L) return(NA_real_)
    suppressWarnings(mean(cluster::silhouette(labs, D)[, 3], na.rm = TRUE))
  }
  
  cand <- vector("list", length(ks))
  best <- NULL; best_score <- -Inf
  
  # Soften the hard minimum slightly with tolerance; allow a few undersized clusters
  hard_min <- max(1L, floor(min_cluster_size * (1 - size_tolerance_frac)))
  
  for (i in seq_along(ks)) {
    kk <- ks[i]
    km <- stats::kmeans(X, centers = kk, nstart = nstart, iter.max = iter.max, algorithm = algorithm)
    
    # sizes via tabulate to ensure length == kk (even if some empty)
    sizes <- as.integer(tabulate(km$cluster, nbins = kk))
    
    n_undersized <- sum(sizes < hard_min)
    if (n_undersized > undersize_allow_frac * kk) {
      cand[[i]] <- NULL
      next
    }
    
    sil <- mean_sil(km$cluster)
    if (is.na(sil)) {
      cand[[i]] <- NULL
      next
    }
    
    cand[[i]] <- list(k = kk, km = km, silhouette = sil, sizes = sizes)
    if (sil > best_score) {
      best_score <- sil
      best <- cand[[i]]
    }
  }
  
  cand <- Filter(Negate(is.null), cand)
  if (!length(cand)) {
    if (isTRUE(relax_on_fail)) {
      warning(sprintf(
        "k-means: no valid candidates (min_cluster_size=%d, k_min=%d). Relaxing guards and retrying once.",
        min_cluster_size, k_min
      ), call. = FALSE)
      return(select_kmeans_k(
        X = X, k_folds = k_folds,
        min_clusters_per_fold = max(2L, floor(min_clusters_per_fold * relax_factor)),
        min_cluster_size      = max(1L, floor(min_cluster_size * relax_factor)),
        k_max_mult            = max(5L, floor(k_max_mult * relax_factor)),
        seed = seed,
        sil_min = max(0, sil_min - 0.03),
        sil_tol = sil_tol,
        strict_enforce = FALSE,
        nstart = nstart, iter.max = iter.max, algorithm = algorithm,
        size_tolerance_frac   = min(0.20, size_tolerance_frac + 0.05),
        undersize_allow_frac  = min(0.25, undersize_allow_frac + 0.05),
        relax_on_fail = FALSE,  # only relax once
        relax_factor  = relax_factor,
        sil_sample_cap = sil_sample_cap
      ))
    } else {
      stop("k-means search produced no valid candidates (after relaxation disabled).")
    }
  }
  
  # ----- Enforce silhouette threshold -----
  sils <- vapply(cand, `[[`, numeric(1), "silhouette")
  ks_v <- vapply(cand, `[[`, integer(1), "k")
  ok   <- which(sils >= sil_min)
  
  if (length(ok)) {
    max_sil <- max(sils[ok])
    near_ok <- ok[ sils[ok] >= (max_sil - sil_tol) ]
    pick    <- near_ok[ which.max(ks_v[near_ok]) ]  # prefer larger k among near-best
    out     <- cand[[pick]]
    out$met_threshold <- TRUE
    return(out)
  } else {
    msg <- sprintf("No k met silhouette threshold (sil_min=%.2f). Best sil=%.3f at k=%d.",
                   sil_min, max(sils), ks_v[ which.max(sils) ])
    if (isTRUE(strict_enforce)) stop(msg)
    warning(msg, call. = FALSE)
    out <- cand[[ which.max(sils) ]]
    out$met_threshold <- FALSE
    return(out)
  }
}


# ------------------ Cluster -> Fold balancing ------------------
balance_folds <- function(clust_summ, global_prop, k = k, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  sev_levels_local <- sev_levels[sev_levels %in% names(global_prop)]
  if (length(sev_levels_local) == 0L) sev_levels_local <- names(global_prop)
  
  fold_sizes <- integer(k)
  fold_tabs  <- replicate(k, setNames(integer(length(sev_levels_local)), sev_levels_local), simplify = FALSE)
  
  score_state <- function(sizes, tabs) {
    size_score <- stats::var(sizes); if (is.na(size_score)) size_score <- 0
    l1s <- numeric(k); over_pen <- 0
    for (f in seq_len(k)) {
      vec <- tabs[[f]]; vec[is.na(vec)] <- 0L
      ntot <- sum(vec, na.rm = TRUE)
      if (ntot == 0) {
        l1s[f] <- 0
      } else {
        p <- as.numeric(vec) / ntot
        l1s[f] <- sum(abs(p - global_prop[sev_levels_local]))
        dev <- abs(p - global_prop[sev_levels_local])
        exceed <- pmax(dev - max_dev, 0)
        over_pen <- over_pen + sum(exceed)
      }
    }
    sev_score <- mean(l1s) + big_pen * over_pen
    w_size * size_score + w_sev * sev_score
  }
  
  cluster_to_fold <- integer(nrow(clust_summ))
  for (i in seq_len(nrow(clust_summ))) {
    n_i <- clust_summ$n[i]
    t_i <- align_counts(clust_summ$sev_tab[[i]], sev_levels_local)
    
    best_f <- NA_integer_; best_sc <- Inf
    for (f in seq_len(k)) {
      sizes_try <- fold_sizes; tabs_try <- fold_tabs
      sizes_try[f] <- sizes_try[f] + n_i
      tabs_try[[f]] <- tabs_try[[f]] + t_i
      sc <- score_state(sizes_try, tabs_try)
      if (sc < best_sc) { best_sc <- sc; best_f <- f }
    }
    cluster_to_fold[i] <- best_f
    fold_sizes[best_f] <- fold_sizes[best_f] + n_i
    fold_tabs[[best_f]] <- fold_tabs[[best_f]] + t_i
  }
  
  list(
    cluster_to_fold = cluster_to_fold,
    fold_sizes = fold_sizes,
    fold_tabs  = fold_tabs
  )
}

assign_with_restarts <- function(clust_summ, global_prop,
                                 k = k,
                                 n_restarts = n_restarts,
                                 seed = seed) {
  set.seed(seed)
  
  score_final <- function(fold_sizes, fold_tabs) {
    size_score <- stats::var(fold_sizes); if (is.na(size_score)) size_score <- 0
    l1s <- numeric(k); over_pen <- 0
    
    sev_levels_local <- sev_levels[sev_levels %in% names(global_prop)]
    if (length(sev_levels_local) == 0L) sev_levels_local <- names(global_prop)
    
    for (f in seq_len(k)) {
      vec <- fold_tabs[[f]]; vec[is.na(vec)] <- 0L
      ntot <- sum(vec, na.rm = TRUE)
      if (ntot == 0) { l1s[f] <- 0; next }
      p <- as.numeric(vec[sev_levels_local]) / ntot
      target <- as.numeric(global_prop[sev_levels_local])
      l1s[f] <- sum(abs(p - target))
      dev <- abs(p - target)
      exceed <- pmax(dev - max_dev, 0)
      over_pen <- over_pen + sum(exceed)
    }
    sev_score <- mean(l1s) + big_pen * over_pen
    w_size * size_score + w_sev * sev_score
  }
  
  best <- NULL; best_score <- Inf
  for (r in seq_len(n_restarts)) {
    ord <- sample(seq_len(nrow(clust_summ)))
    attempt <- balance_folds(
      clust_summ[ord, , drop = FALSE],
      global_prop,
      k = k,
      seed = NULL
    )
    sc <- score_final(attempt$fold_sizes, attempt$fold_tabs)
    if (sc < best_score) {
      best_score <- sc
      restore <- integer(length(ord)); restore[ord] <- attempt$cluster_to_fold
      best <- list(
        cluster_to_fold = restore,
        fold_sizes      = attempt$fold_sizes,
        fold_tabs       = attempt$fold_tabs,
        score           = sc,
        restart         = r,
        cluster_order   = ord
      )
    }
  }
  best
}


# ==============================================================
# YSF PIPELINE — ENV BLOCKING & FOLD/CLUSTER MAPPING 
# ==============================================================
run_env_block_one_ysf <- function(target_ysf,
                                  cwd_var      = cwd_zscore_var,
                                  k_folds      = k,
                                  seed_cluster = seed,
                                  show_panel   = TRUE,
                                  strict_enforce = strict_enforce_kmeans) {
  message("\n========================")
  message("Environmental Blocking — YSF = ", target_ysf)
  message("========================")
  
  # ---------------- Subset & prep points ----------------
  pts <- data_long_reduced_sf %>%
    dplyr::filter(years_since_fire == target_ysf,
                  !is.na(delta_ndvi_min),
                  !is.na(.data[[cwd_var]])) %>%
    sf::st_transform(crs = 32612)
  
  n_pts <- nrow(pts)
  if (n_pts < k_folds) {
    stop("Not enough points (", n_pts, ") for k_folds = ", k_folds,
         " after filtering YSF = ", target_ysf, ".")
  }
  
  # ---- Dynamically adjust cluster parameters ----
  if (isTRUE(use_dynamic_clustering)) {
    dyn <- get_cluster_params(n_pts, k_folds = k_folds)
    min_clusters_per_fold <- dyn$min_clusters_per_fold
    min_cluster_size <- dyn$min_cluster_size
    k_max_mult <- dyn$k_max_mult
  }
  
  # ----- Environmental columns used for clustering -----
  env_cols <- c(
    "twi",
    "pptz_JJA",
    #"tmaxz_JJA",  # max temp
    "tmeanz_JJA",    # mean temp
    "swez_Apr",
    cwd_var,
    #"elevation_m"  # elevation added 2/19.
    "cwd_30yrAvg_TT"  # added 2/19.
  )
  
  env_df <- sf::st_drop_geometry(pts)
  
  # Check presence and median-impute NA
  for (col in env_cols) {
    if (!col %in% names(env_df)) stop("Missing column: ", col)
    med <- suppressWarnings(stats::median(env_df[[col]], na.rm = TRUE)); if (is.na(med)) med <- 0
    nas <- is.na(env_df[[col]]); if (any(nas)) env_df[[col]][nas] <- med
  }
  
  # Drop zero-variance env columns
  X_raw <- as.matrix(env_df[, env_cols, drop = FALSE])
  sdv   <- apply(X_raw, 2, stats::sd, na.rm = TRUE)
  keep  <- is.finite(sdv) & (sdv > 0)
  if (!any(keep)) stop("All env columns have zero variance for YSF = ", target_ysf, ".")
  X <- scale(X_raw[, keep, drop = FALSE])
  
  # ---- PCA dimensionality reduction, pre-clustering ----
  pc <- prcomp(X, center = FALSE, scale. = FALSE)  # X already scaled
  # ---- Extract variance explained ----
  pca_var_tbl <- get_pca_variance_table(
    pc_obj = pc,
    ysf    = target_ysf,
    max_pcs = ncol(pc$x)
  )
  
  Xp <- pc$x
  
  if (is.null(dim(Xp))) Xp <- matrix(Xp, ncol = 1)
  Xp <- Xp[, 1:min(5, ncol(Xp)), drop = FALSE]     # up to 5 PCs
  
  # ---- k-means dynamic-k (thresholds set above) ----
  km_sel <- select_kmeans_k(
    Xp,
    k_folds = k_folds,
    min_clusters_per_fold = min_clusters_per_fold,
    min_cluster_size = min_cluster_size,
    k_max_mult = k_max_mult,
    seed = seed_cluster,
    sil_min = sil_min,       
    sil_tol = sil_tol,
    strict_enforce = strict_enforce
  )
  
  n_clusters <- km_sel$k
  pts <- pts %>% dplyr::mutate(cluster_id = km_sel$km$cluster)
  sizes <- as.integer(table(km_sel$km$cluster))
  
  message("Chosen method = k-means",
          " | No. of clusters = ", n_clusters,
          " | silhouette≈", sprintf("%.3f", km_sel$silhouette))
  message(sprintf("Cluster size summary: min=%d | max=%d | mean=%.1f | total n=%d",
                  min(sizes), max(sizes), mean(sizes), sum(sizes)))
  
  # ----- Severity summary by sev_group -----
  pts <- pts %>% dplyr::mutate(sev_group = factor(sev_group, levels = sev_levels))
  
  clust_summ <- pts |>
    sf::st_drop_geometry() |>
    dplyr::group_by(cluster_id) |>
    dplyr::summarise(
      n       = dplyr::n(),
      sev_tab = list(table(factor(sev_group, levels = sev_levels))),
      .groups = "drop"
    )
  
  global_tab  <- table(factor(sf::st_drop_geometry(pts)$sev_group, levels = sev_levels))
  global_prop <- as.numeric(global_tab) / sum(global_tab)
  names(global_prop) <- sev_levels
  
  # ----- Assign clusters -> folds (multi-start) -------
  res <- assign_with_restarts(
    clust_summ,
    global_prop,
    k           = k_folds,
    n_restarts  = n_restarts,
    seed        = seed
  )
  
  assign_df <- clust_summ %>%
    dplyr::mutate(fold = res$cluster_to_fold) %>%
    dplyr::select(cluster_id, fold)
  
  data_with_folds_env <- pts |>
    dplyr::left_join(assign_df, by = "cluster_id") |>
    dplyr::mutate(
      fold = factor(as.integer(fold), levels = 1:k_folds)
    )
  
  # ---------------- Fold-balance tables ----------------
  sev_bal_tbl <- data_with_folds_env |>
    sf::st_drop_geometry() |>
    dplyr::mutate(sev_group = factor(sev_group, levels = sev_levels)) |>
    dplyr::count(fold, sev_group, .drop = FALSE) |>
    dplyr::group_by(fold) |>
    dplyr::mutate(fold_total = sum(n),
                  pct_of_fold = n / fold_total) |>
    dplyr::ungroup() |>
    dplyr::arrange(fold, sev_group)
  
  sev_counts_wide <- sev_bal_tbl |>
    dplyr::select(fold, sev_group, n) |>
    tidyr::pivot_wider(names_from = sev_group, values_from = n, values_fill = 0)
  
  sev_props_wide <- sev_bal_tbl |>
    dplyr::select(fold, sev_group, pct_of_fold) |>
    tidyr::pivot_wider(names_from = sev_group, values_from = pct_of_fold, values_fill = 0,
                       names_glue = "{.value}_{sev_group}")
  
  fold_totals <- sev_bal_tbl |>
    dplyr::distinct(fold, fold_total)
  
  fold_balance_summary <- fold_totals |>
    dplyr::left_join(sev_counts_wide, by = "fold") |>
    dplyr::left_join(sev_props_wide, by = "fold") |>
    dplyr::arrange(fold)
  
  global_prop_tbl <- tibble::tibble(
    sev_group = sev_levels,
    global_prop = as.numeric(global_prop)
  )
  
  sev_dev_tbl <- sev_bal_tbl |>
    dplyr::left_join(global_prop_tbl, by = "sev_group") |>
    dplyr::mutate(abs_dev = abs(pct_of_fold - global_prop))
  
  message("Fold size + severity composition (wide):")
  print(fold_balance_summary)
  
  message("Per-fold severity proportions vs global (long):")
  print(sev_dev_tbl)
  
  # ---------------- Maps: clusters & folds --------------
  huc12_loc   <- sf::st_transform(huc12,  sf::st_crs(data_with_folds_env))
  streams_loc <- sf::st_transform(streams, sf::st_crs(data_with_folds_env))
  
  pal <- if (n_clusters <= 12) {
    RColorBrewer::brewer.pal(n_clusters, "Set3")
  } else if (n_clusters <= 20) {
    grDevices::hcl.colors(n_clusters, "Dark3")
  } else {
    grDevices::hcl.colors(n_clusters, "Spectral")
  }
  
  p_clusters <- ggplot() +
    geom_sf(data = huc12_loc, fill = NA, color = "grey30", linewidth = 0.3) +
    geom_sf(data = streams_loc, color = "grey60", linewidth = 0.2) +
    geom_sf(data = data_with_folds_env,
            aes(color = factor(cluster_id)),
            size = 0.7, alpha = 0.8) +
    scale_color_manual(values = pal, name = "Cluster ID") +
    labs(
      title = paste("Environmental Clusters — YSF =", target_ysf),
      subtitle = paste0("Clusters: ", n_clusters, " | ", k_folds, " Folds"),
      x = "Easting (m)", y = "Northing (m)")
  
  p_folds <- ggplot() +
    geom_sf(data = huc12_loc, fill = NA, color = "grey30", linewidth = 0.3) +
    geom_sf(data = streams_loc, color = "grey60", linewidth = 0.2) +
    geom_sf(data = data_with_folds_env,
            aes(color = fold),
            size = 0.9, alpha = 0.9) +
    scale_color_brewer(palette = "Set1", name = "Fold ID (Env-block)") +
    labs(
      title = paste("Env-Block Fold Assignment — YSF =", target_ysf),
      subtitle = paste0("k = ", k_folds,
                        " | Severity balance weights: w_size=", w_size,
                        ", w_sev=", w_sev),
      x = "Easting (m)", y = "Northing (m)")
  
  # ---------- Combine into one panel and draw ----------
  panel_grob <- gridExtra::arrangeGrob(p_clusters, p_folds, ncol = 2)
  if (isTRUE(show_panel)) {
    grid::grid.newpage()
    grid::grid.draw(panel_grob)
  }
  
  # Return
  list(
    YSF                  = target_ysf,
    n_clusters           = n_clusters,
    silhouette           = km_sel$silhouette,        
    cluster_sizes        = sizes,             
    data_with_folds_env  = data_with_folds_env,
    clust_summ           = clust_summ,
    global_prop          = global_prop,
    sev_balance_table    = sev_bal_tbl,
    fold_balance_summary = fold_balance_summary,
    sev_deviation_table  = sev_dev_tbl,
    p_clusters           = p_clusters,
    p_folds              = p_folds,
    panel_grob           = panel_grob,
    pca_variance_table = pca_var_tbl  # New 2/9. Principal components variance.
  )
}


# ==============================================================
# ENVIRONMENTAL BLOCKING FOR ALL TARGET YSFs
# ==============================================================
envblock_results <- lapply(
  ysf_set,
  function(y) run_env_block_one_ysf(
    target_ysf   = y,
    cwd_var      = cwd_zscore_var,
    k_folds      = k,
    seed_cluster = seed
  )
)

invisible(lapply(envblock_results, function(res) {
  print(res$p_clusters)
  print(res$p_folds)
}))

envblock_sev_balance_all <- dplyr::bind_rows(lapply(envblock_results, `[[`, "sev_balance_table"),
                                             .id = NULL) %>%
  dplyr::mutate(YSF = rep(ysf_set, times = sapply(envblock_results, function(x) nrow(x$sev_balance_table)))) %>%
  dplyr::relocate(YSF, .before = 1)

envblock_fold_summary_all <- dplyr::bind_rows(lapply(envblock_results, `[[`, "fold_balance_summary"),
                                              .id = NULL) %>%
  dplyr::mutate(YSF = rep(ysf_set, times = sapply(envblock_results, function(x) nrow(x$fold_balance_summary)))) %>%
  dplyr::relocate(YSF, .before = 1)

message("\n=== Combined: Per-fold counts & proportions (all YSFs) ===")
print(envblock_fold_summary_all)

# New 2/9 - PCA variance table.
pca_variance_all <- dplyr::bind_rows(
  lapply(envblock_results, `[[`, "pca_variance_table")
)
pca_summary <- pca_variance_all %>%
  dplyr::filter(PC %in% paste0("PC", 1:5)) %>%
  dplyr::group_by(YSF) %>%
  dplyr::summarise(
    pct_variance_PC1_5 = max(cum_pct_explained),
    .groups = "drop"
  )
print(pca_summary)
# ==============================================================
# SAVE ENVIRONMENTAL BLOCKING RESULTS (RDS)
# ==============================================================

out_dir <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/Ch_01_Figures/40m_spacing/ENV_BLOCKING"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Timestamped save (for provenance / rollback)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
#saveRDS(
#  envblock_results,
#  file = file.path(out_dir, paste0("envblock_results_DWD5CWD30_twi", stamp, ".rds"))  ## EDIT FILE NAME HERE
#)

# Rolling "latest" save
#saveRDS(
#  envblock_results,
#  file = file.path(out_dir, "envblock_results_latest.rds")
#)

#cat("Saved environmental blocking results to:\n", out_dir, "\n")
# ------------------------------------------------------------------------------




#          -------------------------------------------                         #
#########     LOAD ENVIRONMENTAL BLOCKING RESULTS (SKIP RECOMPUTE)   ###########
####          -------------------------------------------                   ####

# Load latest file. 
#rds_path <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/Ch_01_Figures/40m_spacing/ENV_BLOCKING/envblock_results_latest.rds"
#envblock_results <- readRDS(rds_path)

# Alternatively... load a specific file based on new covariates used for clustering. 

## NO VCI / ADDED ELEVATION TO CLUSTERING:
#rds_path_el <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/Ch_01_Figures/40m_spacing/ENV_BLOCKING/envblock_results_ELEVATION20260219_143034.rds" 
#envblock_results <- readRDS(rds_path_el)

## NO VCI / ADDED cwd_30yrAvg TO CLUSTERING:
#rds_path_cwd <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/Ch_01_Figures/40m_spacing/ENV_BLOCKING/envblock_results_CWD30YRNORMAL20260219_200950.rds" 
#envblock_results <- readRDS(rds_path_cwd)

## LOCO MODEL A: CWD5, TWI
#rds_path_A <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/Ch_01_Figures/40m_spacing/ENV_BLOCKING/envblock_results_CWD5_twi20260220_170447.rds" 
#envblock_results <- readRDS(rds_path_A)

## LOCO MODEL B: CWD5, NO TWI
#rds_path_B <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/Ch_01_Figures/40m_spacing/ENV_BLOCKING/envblock_results_CWD5_NOtwi20260221_232506.rds" 
#envblock_results <- readRDS(rds_path_B)

## LOCO MODEL C: CWD30, TWI
#rds_path_C <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/Ch_01_Figures/40m_spacing/ENV_BLOCKING/envblock_results_CWD30_twi20260223_163124.rds" 
#envblock_results <- readRDS(rds_path_C)

## LOCO MODEL D: CWD30, NO TWI
#rds_path_D <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/Ch_01_Figures/40m_spacing/ENV_BLOCKING/envblock_results_CWD30_NOtwi20260223_232416.rds" 
#envblock_results <- readRDS(rds_path_D)

## ==============================================

## LOCO MODEL E: CWD30, CWD5, TWI  (USE FOR 5 & 20 YSF)
rds_path_E <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/Ch_01_Figures/40m_spacing/ENV_BLOCKING/envblock_results_CWD5CWD30_twi20260307_095447.rds" 
envblock_results <- readRDS(rds_path_E)

## LOCO MODEL F: CWD30, CWD5, No TWI  (USE FOR 10 AND 15 YSF)
#rds_path_F <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/Ch_01_Figures/40m_spacing/ENV_BLOCKING/envblock_results_CWD5CWD30_NO.twi.rds" 
#envblock_results <- readRDS(rds_path_F)


####   AFTER UPLOADING ENV BLOCKING OBJECT:  ============
# Re-name by YSF (critical for downstream code!!)
envblock_by_ysf <- function(envblock_results) {
  ys <- vapply(envblock_results, function(x) x$YSF, numeric(1))
  names(envblock_results) <- as.character(ys)
  envblock_results
}
envblock_results_named <- envblock_by_ysf(envblock_results)
# ------------------------------------------------------------------------------


#          -------------------------------------------                         #
#########     RANDOM FOREST MODELING (NOTE THREE OPTIONS!)     #################
#          -------------------------------------------                         #

# =====================================================================
# 1. RANDOM FOREST (RF) — USING ENV-BLOCK FOLDS FROM PRIOR STEP: SEV_GROUP
# =====================================================================
# Create RF function -------------------------------------------
#rf_run_one_ysf_sevgroup <- function(envblock_result,
                           trees   = ntree,
                           seed_rf = 23,
                           cwd_var = cwd_zscore_var) {
  stopifnot(is.list(envblock_result),
            "data_with_folds_env" %in% names(envblock_result))
  
  target_ysf <- envblock_result$YSF
  data_with_folds_env <- envblock_result$data_with_folds_env
  k_folds <- length(levels(data_with_folds_env$fold))
  
  message("\n========================")
  message("RF CV with Env-Blocks - YSF = ", target_ysf)
  message("========================")
  
  rf_factor_cols  <- c("sev_group")
  rf_numeric_cols <- c("twi",
                       "pptz_JJA",
                       #"tmaxz_JJA",  # max temp
                       "tmeanz_JJA",  # mean temp
                       "swez_Apr",
                       cwd_var
                       #"elevation_m"
                       )
  
  # ---------- Pearson correlation screen (numeric covariates only) ----------
  # pairwise Pearson r across numeric predictors; keep |r| > 0.7
  num_df <- sf::st_drop_geometry(data_with_folds_env)[, rf_numeric_cols, drop = FALSE]
  # coerce to numeric just in case
  for (cc in names(num_df)) num_df[[cc]] <- suppressWarnings(as.numeric(num_df[[cc]]))
  cm <- tryCatch(stats::cor(num_df, use = "pairwise.complete.obs", method = "pearson"),
                 error = function(e) NULL)
  if (is.null(cm)) {
    corr_summary <- tibble::tibble(
      YSF = target_ysf, Var1 = character(), Var2 = character(),
      Pearson_r = numeric(), Abs_r = numeric()
    )
  } else {
    nm <- colnames(cm)
    idx <- which(upper.tri(cm), arr.ind = TRUE)
    corr_summary <- tibble::tibble(
      YSF = target_ysf,
      Var1 = nm[idx[, 1]],
      Var2 = nm[idx[, 2]],
      Pearson_r = cm[idx],
      Abs_r = abs(cm[idx])
    ) %>%
      dplyr::filter(Abs_r > 0.7) %>%
      dplyr::arrange(dplyr::desc(Abs_r))
  }
  # --------------------------------------------------------------------------
  
  fold_results    <- vector("list", length = k_folds)
  names(fold_results) <- paste0("Fold_", levels(data_with_folds_env$fold))
  importance_list <- vector("list", length = k_folds)
  names(importance_list) <- names(fold_results)
  
  for (kf in levels(data_with_folds_env$fold)) {
    train_data <- data_with_folds_env %>% dplyr::filter(fold != kf)
    test_data  <- data_with_folds_env %>% dplyr::filter(fold == kf)
    
    for (col in intersect(rf_factor_cols, names(train_data))) {
      tr <- factor(train_data[[col]]); te <- factor(test_data[[col]])
      all_lvls <- union(levels(tr), levels(te))
      train_data[[col]] <- factor(as.character(train_data[[col]]), levels = all_lvls)
      test_data[[col]]  <- factor(as.character(test_data[[col]]),  levels = all_lvls)
    }
    
    for (col in intersect(rf_numeric_cols, names(train_data))) {
      med <- suppressWarnings(stats::median(train_data[[col]], na.rm = TRUE)); if (is.na(med)) med <- 0
      train_data[[col]][is.na(train_data[[col]])] <- med
      test_data[[col]][is.na(test_data[[col]])]  <- med
    }
    
    set.seed(seed_rf + as.integer(as.character(kf)))
    rf_formula <- as.formula(paste(
      "delta_ndvi_min ~",
      paste(c("sev_group",
              "twi",
              "pptz_JJA",
              #"tmaxz_JJA",
              "tmeanz_JJA",
              "swez_Apr",
              "veg_climate_index_08",
              cwd_var),
            collapse = " + ")
    ))
    
    rf_model <- randomForest::randomForest(
      rf_formula,
      data = train_data,
      ntree = trees,
      importance = TRUE
    )
    
    predictions <- predict(rf_model, newdata = test_data)
    observed    <- test_data$delta_ndvi_min
    N           <- length(observed)
    mse         <- if (N) mean((predictions - observed)^2) else NA_real_
    rmse        <- if (is.na(mse)) NA_real_ else sqrt(mse)
    denom       <- sum((observed - mean(observed))^2)
    pseudo_r2   <- if (denom == 0) NA_real_ else 1 - (sum((observed - predictions)^2) / denom)
    
    fold_results[[paste0("Fold_", kf)]] <- list(
      N = N, mse = mse, rmse = rmse, pseudo_r2 = pseudo_r2,
      model = rf_model, predictions = predictions, observed = observed
    )
    
    importance_list[[paste0("Fold_", kf)]] <- randomForest::importance(rf_model, type = 1)
  }
  
  performance_df <- dplyr::bind_rows(lapply(names(fold_results), function(nm) {
    fr <- fold_results[[nm]]
    tibble::tibble(
      YSF = target_ysf,
      Fold = nm,
      N = fr$N,
      MSE = fr$mse,
      RMSE = fr$rmse,
      Pseudo_R2 = fr$pseudo_r2
    )
  }))
  
  mean_rmse  <- mean(performance_df$RMSE, na.rm = TRUE)
  mean_r2    <- mean(performance_df$Pseudo_R2, na.rm = TRUE)
  w          <- ifelse(is.na(performance_df$N), 0, performance_df$N)
  wmean_rmse <- stats::weighted.mean(performance_df$RMSE, w, na.rm = TRUE)
  wmean_r2   <- stats::weighted.mean(performance_df$Pseudo_R2, w, na.rm = TRUE)
  
  perf_row <- tibble::tibble(
    YSF = target_ysf,
    Mean_RMSE_unweighted     = mean_rmse,
    Mean_PseudoR2_unweighted = mean_r2,
    Mean_RMSE_weighted       = wmean_rmse,
    Mean_PseudoR2_weighted   = wmean_r2
  )
  
  importance_df <- dplyr::bind_rows(
    lapply(importance_list, function(mat) {
      if (is.null(mat)) return(NULL)
      df_out <- as.data.frame(mat)
      df_out$Variable <- rownames(df_out)
      df_out
    }),
    .id = "Fold"
  )
  
  if (!nrow(importance_df)) {
    vi_summary <- tibble::tibble(YSF = target_ysf,
                                 Variable = character(),
                                 Mean_IncMSE = numeric(),
                                 SD_IncMSE = numeric(),
                                 W_Mean_IncMSE = numeric())
  } else {
    fold_sizes <- performance_df %>% dplyr::select(Fold, N)
    importance_df <- importance_df %>% dplyr::left_join(fold_sizes, by = "Fold")
    vi_summary <- importance_df %>%
      dplyr::group_by(Variable) %>%
      dplyr::summarise(
        Mean_IncMSE   = mean(`%IncMSE`, na.rm = TRUE),
        SD_IncMSE     = stats::sd(`%IncMSE`, na.rm = TRUE),
        W_Mean_IncMSE = stats::weighted.mean(`%IncMSE`, w = ifelse(is.na(N), 0, N), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(W_Mean_IncMSE)) %>%
      dplyr::mutate(YSF = target_ysf, .before = 1)
  }
  
  list(
    YSF            = target_ysf,
    corr_summary   = corr_summary,   
    perf_row       = perf_row,
    perf_folds     = performance_df,
    vi_summary     = vi_summary,
    fold_results   = fold_results,
    importance     = importance_list
  )
}
# ---------------------------------------------------------
# ------------------ Run RF for all YSFs ------------------
envblock_by_ysf <- function(envblock_results) {
  ys <- vapply(envblock_results, function(x) x$YSF, numeric(1))
  names(envblock_results) <- as.character(ys)
  envblock_results
}
envblock_results_named <- envblock_by_ysf(envblock_results)

# Run RF:
#rf_results_envblock <- lapply(
  as.character(ysf_set),
  function(ysf_key) {
    rf_run_one_ysf_sevgroup(
      envblock_result = envblock_results_named[[ysf_key]],
      trees   = ntree,
      seed_rf = 168,
      cwd_var = cwd_zscore_var
    )
  }
)

# ------------------ Print Summary Tables ------------------
names(rf_results_envblock) <- as.character(ysf_set)

perf_rows_all  <- dplyr::bind_rows(lapply(rf_results_envblock, `[[`, "perf_row"))
perf_folds_all <- dplyr::bind_rows(lapply(rf_results_envblock, `[[`, "perf_folds"))
vi_summary_all <- dplyr::bind_rows(lapply(rf_results_envblock, `[[`, "vi_summary"))

# Gather correlation summaries (Pearson's r)
corrs_all <- dplyr::bind_rows(lapply(rf_results_envblock, `[[`, "corr_summary"))

cat("\n================ High Covariate Correlations (|r| > 0.7) (by YSF) ================\n")
if (nrow(corrs_all) == 0) {
  print(tibble::tibble(YSF = numeric(), Var1 = character(), Var2 = character(),
                       Pearson_r = numeric(), Abs_r = numeric()))
} else {
  print(corrs_all, n = Inf)
}

cat("\n================ Combined Performance (by YSF) ================\n")
print(perf_rows_all)

cat("\n================ Fold-level Performance (all YSF) ================\n")
print(perf_folds_all, n = Inf)

cat("\n================ Variable Importance Summary (all YSF) ================\n")
print(vi_summary_all, n = Inf)

# ---------------------------------------------------------
#  PLOT: Variable Importance across all YSFs
# ==============================================================

# Order variables by order they're described in the methods 1/19 (descending)
var_order_pub <- c(
  "sev_group",
  "pptz_JJA",
  "tmeanz_JJA",
  "swez_Apr",
  "cwd_5yr_zscore_08",
  "twi",
  "veg_climate_index_08"
)

#### Other ordering option: order by the earliest YSF rankings:
#first_ysf <- if (exists("ysf_set")) ysf_set[1] else sort(unique(vi_summary_all$YSF))[1]
#var_levels_first <- vi_summary_all %>%
#  dplyr::filter(YSF == first_ysf) %>%
#  dplyr::arrange(dplyr::desc(Mean_IncMSE)) %>%
#  dplyr::pull(Variable)

# Use descending YSF levels so 5 appears first in legend (then reversed)
ysf_levels_plot <- sort(unique(vi_summary_all$YSF), decreasing = TRUE)  # e.g., 20,15,10,5

vi_plot_df <- vi_summary_all %>%
  dplyr::mutate(
    Variable = factor(Variable, levels = rev(var_order_pub)),  # reverse for coord_flip
    YSF      = factor(YSF, levels = ysf_levels_plot)
  )

pd <- position_dodge(width = 0.6)

# Customize Variable Importance Plot ------------------------

# ---- Define YSF colors manually 
ysf_colors <- c(
  "1"  = "aquamarine",
  "5"  = "darkred",
  "10" = "darkgreen",
  "15" = "darkblue",
  "20" = "darkviolet"
)

# ---- Rename variables for nicer y-axis labels 
var_labels <- c(
  "sev_group"            = "Fire Severity Class",
  "twi"                  = "Topographic Wetness Index (TWI)",
  "pptz_JJA"             = "Summer Precipitation",
  "tmaxz_JJA"           = "Summer Maximum Temperature",
  "tmeanz_JJA"          = "Summer Mean Temperature",
  "swez_Apr"             = "Snow Water Equivalent, April",
  "cwd_5yr_zscore_08"    = "5-Year Post-Fire Aug CWD",
  "cwd_30yrAvg_TT"       = "CWD 30-Year Normal"
)

# ---- Generate plot -----------------------------------------------------------
ggplot(vi_plot_df,
       aes(x = Variable, y = Mean_IncMSE, color = YSF, group = YSF)) +
  geom_errorbar(aes(ymin = Mean_IncMSE - SD_IncMSE,
                    ymax = Mean_IncMSE + SD_IncMSE),
                width = 0.25,
                position = pd,
                linewidth = 1.2) +   # << bolder lines
  
  geom_point(size = 4.0,          # <<  point size
             stroke = 1.2,        # << thicker point outline
             position = pd) +
  coord_flip() +
  scale_color_manual(values = ysf_colors, name = "Years Since Fire") +
  scale_x_discrete(labels = var_labels) +  # rename variables on y-axis
  labs(
    title    = "Variable Importance Across Years Since Fire",
    #subtitle = paste("CWD variable:", cwd_zscore_var),
    x = "Predictor Variable",
    y = "Mean %IncMSE (± SD)"
  ) +
  guides(color = guide_legend(reverse = TRUE)) +
  theme_minimal() +
  theme(
    
    # ---- GLOBAL FONT
    text = element_text(family = "Times New Roman"),
    
    # ---- Titles & Subtitles 
    plot.title      = element_text(size = 25, face = "bold"),
    plot.subtitle   = element_text(size = 12, face = "italic"),
    
    # ---- Axis labels
    axis.title.x  = element_text(size = 16, margin = ggplot2::margin(t = 10), face = "bold"),
    axis.title.y  = element_text(size = 16, margin = ggplot2::margin(r = 10), face = "bold"),
    
    # ---- Axis tick labels
    axis.text.x   = element_text(size = 14, margin = ggplot2::margin(t = 8), face = "bold"),
    axis.text.y   = element_text(size = 14, margin = ggplot2::margin(r = 8), face = "bold"),  # variable names
    
    # ---- Legend 
    legend.title    = element_text(size = 14, face = "bold"),   # “Years Since Fire”
    legend.text     = element_text(size = 13),                  # YSF number size
    legend.key.size = unit(1.4, "lines"),                       # symbol size
    legend.key.width = unit(1.0, "lines"),                      # legend item width
    legend.spacing.y = unit(0.2, "cm"),                         # vertical spacing between legend items
    
    # ---- General panel tweaks 
    panel.grid.minor = element_blank()
  )
# ------------------------------------------------------------------------------
################################################################################

# ==============================================================
# 2. RANDOM FOREST (RF) — USING ENV-BLOCK FOLDS FROM PRIOR STEP: SEV_NUM
# =====================================================================
#   *sev_num version* (numeric severity instead of sev_group)
#   - SAME function names, SAME seeds, SAME user controls, ETC
# ==============================================================
# Create RF function -------------------------------------------
rf_run_one_ysf_sevnum <- function(envblock_result,
                           trees   = ntree,
                           seed_rf = 23,
                           cwd_var = cwd_zscore_var) {
  stopifnot(is.list(envblock_result),
            "data_with_folds_env" %in% names(envblock_result))
  
  target_ysf <- envblock_result$YSF
  data_with_folds_env <- envblock_result$data_with_folds_env
  k_folds <- length(levels(data_with_folds_env$fold))
  
  message("\n========================")
  message("RF CV with Env-Blocks (SEV_NUM) - YSF = ", target_ysf)
  message("========================")
  
  # --- CHANGE: sev_group -> sev_num (numeric)
  rf_factor_cols  <- character(0)  
  rf_numeric_cols <- c("sev_num",   # numeric severity
                       "twi",
                       "pptz_JJA",
                       "tmeanz_JJA",  
                       "swez_Apr",
                       cwd_var,                      
                       "cwd_30yrAvg_TT"
                       )
  
  # ---------- Pearson correlation screen (numeric covariates only) ----------
  # pairwise Pearson r across numeric predictors; keep |r| > 0.7
  num_df <- sf::st_drop_geometry(data_with_folds_env)[, rf_numeric_cols, drop = FALSE]
  # coerce to numeric just in case
  for (cc in names(num_df)) num_df[[cc]] <- suppressWarnings(as.numeric(num_df[[cc]]))
  cm <- tryCatch(stats::cor(num_df, use = "pairwise.complete.obs", method = "pearson"),
                 error = function(e) NULL)
  if (is.null(cm)) {
    corr_summary <- tibble::tibble(
      YSF = target_ysf, Var1 = character(), Var2 = character(),
      Pearson_r = numeric(), Abs_r = numeric()
    )
  } else {
    nm <- colnames(cm)
    idx <- which(upper.tri(cm), arr.ind = TRUE)
    corr_summary <- tibble::tibble(
      YSF = target_ysf,
      Var1 = nm[idx[, 1]],
      Var2 = nm[idx[, 2]],
      Pearson_r = cm[idx],
      Abs_r = abs(cm[idx])
    ) %>%
      dplyr::filter(Abs_r > 0.7) %>%
      dplyr::arrange(dplyr::desc(Abs_r))
  }
  # --------------------------------------------------------------------------
  
  fold_results    <- vector("list", length = k_folds)
  names(fold_results) <- paste0("Fold_", levels(data_with_folds_env$fold))
  importance_list <- vector("list", length = k_folds)
  names(importance_list) <- names(fold_results)
  
  for (kf in levels(data_with_folds_env$fold)) {
    train_data <- data_with_folds_env %>% dplyr::filter(fold != kf)
    test_data  <- data_with_folds_env %>% dplyr::filter(fold == kf)
    
    # factor handling (none now, but keep structure)
    for (col in intersect(rf_factor_cols, names(train_data))) {
      tr <- factor(train_data[[col]]); te <- factor(test_data[[col]])
      all_lvls <- union(levels(tr), levels(te))
      train_data[[col]] <- factor(as.character(train_data[[col]]), levels = all_lvls)
      test_data[[col]]  <- factor(as.character(test_data[[col]]),  levels = all_lvls)
    }
    
    # numeric median-impute (train median)
    for (col in intersect(rf_numeric_cols, names(train_data))) {
      train_data[[col]] <- suppressWarnings(as.numeric(train_data[[col]]))
      test_data[[col]]  <- suppressWarnings(as.numeric(test_data[[col]]))
      
      med <- suppressWarnings(stats::median(train_data[[col]], na.rm = TRUE)); if (is.na(med)) med <- 0
      train_data[[col]][is.na(train_data[[col]])] <- med
      test_data[[col]][is.na(test_data[[col]])]  <- med
    }
    
    set.seed(seed_rf + as.integer(as.character(kf)))
    rf_formula <- as.formula(paste(
      "delta_ndvi_min ~",
      paste(c("sev_num",
              "pptz_JJA",
              "tmeanz_JJA",
              "swez_Apr",
              #"veg_climate_index_08",
              "twi",
              "cwd_30yrAvg_TT",
              cwd_var),
            collapse = " + ")
    ))
    
    rf_model <- randomForest::randomForest(
      rf_formula,
      data = train_data,
      ntree = trees,
      importance = TRUE
    )
    
    predictions <- predict(rf_model, newdata = test_data)
    observed    <- test_data$delta_ndvi_min
    N           <- length(observed)
    mse         <- if (N) mean((predictions - observed)^2) else NA_real_
    rmse        <- if (is.na(mse)) NA_real_ else sqrt(mse)
    denom       <- sum((observed - mean(observed))^2)
    pseudo_r2   <- if (denom == 0) NA_real_ else 1 - (sum((observed - predictions)^2) / denom)
    
    fold_results[[paste0("Fold_", kf)]] <- list(
      N = N, mse = mse, rmse = rmse, pseudo_r2 = pseudo_r2,
      model = rf_model, predictions = predictions, observed = observed
    )
    
    importance_list[[paste0("Fold_", kf)]] <- randomForest::importance(rf_model, type = 1)
  }
  
  performance_df <- dplyr::bind_rows(lapply(names(fold_results), function(nm) {
    fr <- fold_results[[nm]]
    tibble::tibble(
      YSF = target_ysf,
      Fold = nm,
      N = fr$N,
      MSE = fr$mse,
      RMSE = fr$rmse,
      Pseudo_R2 = fr$pseudo_r2
    )
  }))
  
  mean_rmse  <- mean(performance_df$RMSE, na.rm = TRUE)
  mean_r2    <- mean(performance_df$Pseudo_R2, na.rm = TRUE)
  w          <- ifelse(is.na(performance_df$N), 0, performance_df$N)
  wmean_rmse <- stats::weighted.mean(performance_df$RMSE, w, na.rm = TRUE)
  wmean_r2   <- stats::weighted.mean(performance_df$Pseudo_R2, w, na.rm = TRUE)
  
  # ---- name columns:
  perf_row <- tibble::tibble(
    YSF = target_ysf,
    Mean_RMSE_unweighted     = mean_rmse,
    Mean_PseudoR2_unweighted = mean_r2,
    Mean_RMSE_weighted       = wmean_rmse,
    Mean_PseudoR2_weighted   = wmean_r2
  )
  
  importance_df <- dplyr::bind_rows(
    lapply(importance_list, function(mat) {
      if (is.null(mat)) return(NULL)
      df_out <- as.data.frame(mat)
      df_out$Variable <- rownames(df_out)
      df_out
    }),
    .id = "Fold"
  )
  
  if (!nrow(importance_df)) {
    vi_summary <- tibble::tibble(YSF = target_ysf,
                                 Variable = character(),
                                 Mean_IncMSE = numeric(),
                                 SD_IncMSE = numeric(),
                                 W_Mean_IncMSE = numeric())
  } else {
    fold_sizes <- performance_df %>% dplyr::select(Fold, N)
    importance_df <- importance_df %>% dplyr::left_join(fold_sizes, by = "Fold")
    vi_summary <- importance_df %>%
      dplyr::group_by(Variable) %>%
      dplyr::summarise(
        Mean_IncMSE   = mean(`%IncMSE`, na.rm = TRUE),
        SD_IncMSE     = stats::sd(`%IncMSE`, na.rm = TRUE),
        W_Mean_IncMSE = stats::weighted.mean(`%IncMSE`, w = ifelse(is.na(N), 0, N), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(W_Mean_IncMSE)) %>%
      dplyr::mutate(YSF = target_ysf, .before = 1)
  }
  
  list(
    YSF            = target_ysf,
    corr_summary   = corr_summary,   
    perf_row       = perf_row,
    perf_folds     = performance_df,
    vi_summary     = vi_summary,
    fold_results   = fold_results,
    importance     = importance_list
  )
}

# ---------------------------------------------------------
# ------------------ Run RF for all YSFs ------------------

envblock_by_ysf <- function(envblock_results) {
  ys <- vapply(envblock_results, function(x) x$YSF, numeric(1))
  names(envblock_results) <- as.character(ys)
  envblock_results
}
envblock_results_named <- envblock_by_ysf(envblock_results)

# Run RF:
rf_results_envblock <- lapply(
  as.character(ysf_set),
  function(ysf_key) {
    rf_run_one_ysf_sevnum(
      envblock_result = envblock_results_named[[ysf_key]],
      trees   = ntree,
      seed_rf = 168,
      cwd_var = cwd_zscore_var
    )
  }
)

# ------------------ Print Summary Tables ------------------
names(rf_results_envblock) <- as.character(ysf_set)

perf_rows_all  <- dplyr::bind_rows(lapply(rf_results_envblock, `[[`, "perf_row"))
perf_folds_all <- dplyr::bind_rows(lapply(rf_results_envblock, `[[`, "perf_folds"))

# ---- Variable labels 
var_labels <- c(
  "sev_num"              = "Fire Severity",
  "twi"                  = "Topographic Wetness Index (TWI)",
  "pptz_JJA"             = "Summer Precipitation",
  "tmaxz_JJA"            = "Summer Maximum Temperature",
  "tmeanz_JJA"           = "Summer Mean Temperature",
  "swez_Apr"             = "Snow Water Equivalent, April",
  #"veg_climate_index_08" = "Veg–Climate Index (Aug)",
  "cwd_5yr_zscore_08"    = "5-Year Post-Fire Aug CWD",
  "elevation_m"          = "Elevation", 
  "cwd_30yrAvg_TT"       = "CWD 30-Year Normal"
)

# Variable Importance Summary
vi_summary_all <- dplyr::bind_rows(lapply(rf_results_envblock, `[[`, "vi_summary"))
vi_summary_all <- vi_summary_all %>%
  dplyr::mutate(
    Variable_label = dplyr::recode(Variable, !!!var_labels)
  )


# Gather correlation summaries (Pearson's r)
corrs_all <- dplyr::bind_rows(lapply(rf_results_envblock, `[[`, "corr_summary"))

cat("\n================ High Covariate Correlations (|r| > 0.7) (by YSF) ================\n")
if (nrow(corrs_all) == 0) {
  print(tibble::tibble(YSF = numeric(), Var1 = character(), Var2 = character(),
                       Pearson_r = numeric(), Abs_r = numeric()))
} else {
  print(corrs_all, n = Inf)
}

cat("\n================ Combined Performance (by YSF) ================\n")
print(perf_rows_all)

cat("\n================ Fold-level Performance (all YSF) ================\n")
print(perf_folds_all, n = Inf)

cat("\n================ Variable Importance Summary (all YSF) ================\n")
print(vi_summary_all, n = Inf)

# ---------------------------------------------------------
#  PLOT: Variable Importance across all YSFs
#   (same plotting code as above; just update var order/labels so sev_num shows)
# ==============================================================

# Order variables by order they're described in the methods 1/19 (descending)
var_order_pub <- c(
  "sev_num",
  "pptz_JJA",
  "tmeanz_JJA",
  "swez_Apr",
  "cwd_5yr_zscore_08",
  "twi",
  "cwd_30yrAvg_TT",
  "veg_climate_index_08"
)

ysf_levels_plot <- sort(unique(vi_summary_all$YSF), decreasing = TRUE)  # e.g., 20,15,10,5

vi_plot_df <- vi_summary_all %>%
  dplyr::mutate(
    Variable = factor(Variable, levels = rev(var_order_pub)),  # reverse for coord_flip
    YSF      = factor(YSF, levels = ysf_levels_plot)
  )

pd <- position_dodge(width = 0.6)

# ---- Define YSF colors manually 
ysf_colors <- c(
  "1"  = "aquamarine",
  "5"  = "darkred",
  "10" = "darkgreen",
  "15" = "darkblue",
  "20" = "darkviolet"
)

# ---- Rename variables for nicer y-axis labels 
var_labels <- c(
  "sev_num"              = "Fire Severity",
  "twi"                  = "Topographic Wetness Index (TWI)",
  "pptz_JJA"             = "Summer Precipitation",
  "tmaxz_JJA"            = "Summer Maximum Temperature",
  "tmeanz_JJA"           = "Summer Mean Temperature",
  "swez_Apr"             = "Snow Water Equivalent, April",
  "veg_climate_index_08" = "Veg–Climate Index (Aug)",
  "cwd_5yr_zscore_08"    = "5-Year Post-Fire Aug CWD",
  "elevation_m"          = "Elevation",
  "cwd_30yrAvg_TT"       = "CWD 30-Year Normal"
)

ggplot(vi_plot_df,
       aes(x = Variable, y = Mean_IncMSE, color = YSF, group = YSF)) +
  geom_errorbar(aes(ymin = Mean_IncMSE - SD_IncMSE,
                    ymax = Mean_IncMSE + SD_IncMSE),
                width = 0.25,
                position = pd,
                linewidth = 1.2) +
  geom_point(size = 4.0,
             stroke = 1.2,
             position = pd) +
  coord_flip() +
  scale_color_manual(values = ysf_colors, name = "Years Since Fire") +
  scale_x_discrete(labels = var_labels) +
  labs(
    title    = "Variable Importance Across Years Since Fire",
    x = "Predictor Variable",
    y = "Mean %IncMSE (± SD)"
  ) +
  guides(color = guide_legend(reverse = TRUE)) +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title      = element_text(size = 25, face = "bold"),
    plot.subtitle   = element_text(size = 12, face = "italic"),
    axis.title.x  = element_text(size = 16, margin = ggplot2::margin(t = 10), face = "bold"),
    axis.title.y  = element_text(size = 16, margin = ggplot2::margin(r = 10), face = "bold"),
    axis.text.x   = element_text(size = 14, margin = ggplot2::margin(t = 8), face = "bold"),
    axis.text.y   = element_text(size = 14, margin = ggplot2::margin(r = 8), face = "bold"),
    legend.title    = element_text(size = 14, face = "bold"),
    legend.text     = element_text(size = 13),
    legend.key.size = unit(1.4, "lines"),
    legend.key.width = unit(1.0, "lines"),
    legend.spacing.y = unit(0.2, "cm"),
    panel.grid.minor = element_blank()
  )
################################################################################

# ==============================================================
# 3. RANDOM FOREST (RF) — SEV_NUM + LEAVE-ONE-OUT MODEL COMPARISON
# =====================================================================
#   *sev_num version* (numeric severity instead of sev_group)
# LOCO (Leave-One-Covariate-Out) predictor sets: ---------------------
rf_vars_full <- c(
  "sev_num",
  "twi",
  "pptz_JJA",
  "tmeanz_JJA",
  "swez_Apr",
  cwd_zscore_var,
  #"elevation_m",
  "cwd_30yrAvg_TT"
)

rf_var_sets <- c(
  list(full = rf_vars_full),
  lapply(rf_vars_full, function(v) setdiff(rf_vars_full, v))
)

names(rf_var_sets) <- c(
  "full",
  paste0("drop_", rf_vars_full)
)

# Create RF function -----------------------------------------------
rf_run_one_ysf_loco <- function(envblock_result,
                                rf_vars,
                                trees   = ntree,
                                seed_rf = 23) {
  
  stopifnot(is.list(envblock_result),
            "data_with_folds_env" %in% names(envblock_result))
  
  target_ysf <- envblock_result$YSF
  df <- envblock_result$data_with_folds_env
  k_folds <- length(levels(df$fold))
  
  message("  → YSF = ", target_ysf,
          " | predictors = ", paste(rf_vars, collapse = ", "))
  
  fold_results <- vector("list", length = k_folds)
  names(fold_results) <- paste0("Fold_", levels(df$fold))
  
  for (kf in levels(df$fold)) {
    
    message("    • Fold ", kf, "/", k_folds)
    
    train_data <- df %>% dplyr::filter(fold != kf)
    test_data  <- df %>% dplyr::filter(fold == kf)
    
    # ---- numeric median imputation (train median) ----
    for (col in rf_vars) {
      med <- suppressWarnings(stats::median(train_data[[col]], na.rm = TRUE))
      if (is.na(med)) med <- 0
      train_data[[col]][is.na(train_data[[col]])] <- med
      test_data[[col]][is.na(test_data[[col]])]   <- med
    }
    
    set.seed(seed_rf + as.integer(as.character(kf)))
    
    rf_formula <- as.formula(
      paste("delta_ndvi_min ~", paste(rf_vars, collapse = " + "))
    )
    
    rf_model <- randomForest::randomForest(
      rf_formula,
      data  = train_data,
      ntree = trees,
      importance = TRUE
    )
    
    # ---- Variable importance (%IncMSE) ----
    vi_tbl <- tryCatch({
      as.data.frame(randomForest::importance(rf_model, type = 1)) |>
        tibble::rownames_to_column("variable") |>
        dplyr::rename(pctIncMSE = `%IncMSE`)
    }, error = function(e) NULL)
    
    # ---- OOB metrics (safe extraction) ----
    oob_mse_val <- tryCatch(tail(rf_model$mse, 1), error = function(e) NA_real_)
    oob_r2_val  <- tryCatch(tail(rf_model$rsq, 1), error = function(e) NA_real_)
    
    preds <- predict(rf_model, newdata = test_data)
    obs   <- test_data$delta_ndvi_min
    N     <- length(obs)
    
    mse  <- if (N) mean((preds - obs)^2) else NA_real_
    rmse <- if (is.na(mse)) NA_real_ else sqrt(mse)
    
    denom <- sum((obs - mean(obs))^2)
    r2    <- if (denom == 0) NA_real_ else 1 - sum((obs - preds)^2) / denom
    
    fold_results[[paste0("Fold_", kf)]] <- list(
      model   = rf_model,
      preds   = preds,
      obs     = obs,
      N       = N,
      rmse    = rmse,
      r2      = r2,
      vi      = vi_tbl,
      oob_mse = oob_mse_val,
      oob_r2  = oob_r2_val
    )
  }
  
  # ---- Aggregate performance ----------------------------------------
  rmse_vals <- sapply(fold_results, `[[`, "rmse")
  r2_vals   <- sapply(fold_results, `[[`, "r2")
  Ns        <- sapply(fold_results, `[[`, "N")
  
  perf <- tibble::tibble(
    YSF = target_ysf,
    Mean_RMSE_unweighted     = mean(rmse_vals, na.rm = TRUE),
    Mean_PseudoR2_unweighted = mean(r2_vals,   na.rm = TRUE),
    Mean_RMSE_weighted       = stats::weighted.mean(rmse_vals, Ns, na.rm = TRUE),
    Mean_PseudoR2_weighted   = stats::weighted.mean(r2_vals,   Ns, na.rm = TRUE)
  )
  
  # ---- Aggregate VI across folds ----
  vi_all <- purrr::map_dfr(
    fold_results,
    "vi",
    .id = "Fold"
  )
  
  vi_summary <- vi_all |>
    dplyr::group_by(variable) |>
    dplyr::summarise(
      Mean_IncMSE = mean(pctIncMSE, na.rm = TRUE),
      SD_IncMSE   = sd(pctIncMSE, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(desc(Mean_IncMSE))
  
  list(
    YSF          = target_ysf,
    rf_vars      = rf_vars,
    fold_models  = fold_results,
    perf         = perf,
    vi_summary   = vi_summary,   # Aggregated VI (main table)
    vi_by_fold   = vi_all        # fold-level VI for diagnostics
  )
}

# ---- RUN LOCO RF MODELS ---------------------------------------------
envblock_by_ysf <- function(envblock_results) {
  ys <- vapply(envblock_results, function(x) x$YSF, numeric(1))
  names(envblock_results) <- as.character(ys)
  envblock_results
}
envblock_results_named <- envblock_by_ysf(envblock_results)

rf_loco_results <- purrr::imap(
  rf_var_sets,
  function(vars, model_name) {
    
    message("\n==============================")
    message("LOCO RF model: ", model_name)
    message("==============================")
    
    res <- lapply(
      as.character(ysf_set),
      function(ysf_key) {
        rf_run_one_ysf_loco(
          envblock_result = envblock_results_named[[ysf_key]],
          rf_vars = vars,
          trees   = ntree,
          seed_rf = 168
        )
      }
    )
    
    names(res) <- as.character(ysf_set)
    
    list(
      model_name = model_name,
      rf_vars    = vars,
      results    = res
    )
  }
)


# ---- PERFORMANCE COMPARISON (ΔRMSE, Δpseudo-R²) ----------------------

perf_loco_all <- purrr::map_dfr(
  rf_loco_results,
  function(m) {
    dplyr::bind_rows(lapply(m$results, `[[`, "perf")) %>%
      dplyr::mutate(
        model  = m$model_name,
        n_vars = length(m$rf_vars)
      )
  }
)

### Compare all models across each YSF 
perf_loco_summary <- perf_loco_all %>%
  
  # Join full-model weighted metrics (baseline per YSF)
  dplyr::left_join(
    perf_loco_all %>%
      dplyr::filter(model == "full") %>%
      dplyr::select(
        YSF,
        RMSE_full = Mean_RMSE_weighted,
        R2_full   = Mean_PseudoR2_weighted
      ),
    by = "YSF"
  ) %>%
  
  # Compute deltas relative to full model
  dplyr::mutate(
    delta_RMSE = Mean_RMSE_weighted     - RMSE_full,
    delta_R2   = Mean_PseudoR2_weighted - R2_full
  ) %>%
  
  # Rank models WITHIN each YSF
  dplyr::group_by(YSF) %>%
  dplyr::arrange(
    Mean_RMSE_weighted,                 # primary: lower RMSE is better
    dplyr::desc(Mean_PseudoR2_weighted),# secondary: higher R² is better
    .by_group = TRUE
  ) %>%
  dplyr::mutate(
    rank_within_YSF = dplyr::row_number()
  ) %>%
  dplyr::ungroup()

print(perf_loco_summary, n = Inf)


# EXTRACT VI FOR TOP-RANKED MODEL WITHIN EACH YSF 
# ---- 1. Identify the winning model per YSF ------
top_models_by_ysf <- perf_loco_summary %>%
  dplyr::filter(rank_within_YSF == 1) %>%
  dplyr::select(YSF, model)

print(top_models_by_ysf)

# ---- 2. Helper to pull VI from rf_loco_results ----
get_vi_for_model_ysf <- function(model_name, ysf_val) {
  
  # locate the model object
  model_obj <- purrr::keep(
    rf_loco_results,
    ~ .x$model_name == model_name
  )[[1]]
  
  # pull VI summary for that YSF
  vi_tbl <- model_obj$results[[as.character(ysf_val)]]$vi_summary
  
  if (is.null(vi_tbl)) return(NULL)
  
  vi_tbl %>%
    dplyr::mutate(
      YSF        = ysf_val,
      top_model  = model_name
    ) %>%
    dplyr::select(YSF, top_model, dplyr::everything())
}

# ---- 3. Build combined VI table ----
vi_top_models_all <- purrr::map_dfr(
  seq_len(nrow(top_models_by_ysf)),
  function(i) {
    get_vi_for_model_ysf(
      model_name = top_models_by_ysf$model[i],
      ysf_val    = top_models_by_ysf$YSF[i]
    )
  }
)

# ---- 4. Arrange nicely ----
vi_top_models_all <- vi_top_models_all %>%
  dplyr::arrange(YSF, dplyr::desc(Mean_IncMSE))

print(vi_top_models_all, n = Inf)

# -------------------------------------------------------
# SAVE LOCO RF RESULTS (timestamped + rolling "latest") ------------------------
out_dir <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/Ch01_RF_LOCO"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

loco_save <- list(
  meta = list(
    created      = Sys.time(),
    ysf_set      = ysf_set,
    predictors   = rf_vars_full,
    ntree        = ntree,
    seed_rf      = 168,
    severity_mode = "sev_num",
    description  = "LOCO Random Forest results with env-block CV"
  ),
  rf_loco_results   = rf_loco_results,
  perf_loco_all     = perf_loco_all,
  perf_loco_summary = perf_loco_summary
)

# ---- Timestamped archive, new covariates 2/19
#saveRDS(
#  loco_save,
#  file = file.path(out_dir, paste0("rf_loco_results_CWD5_no.TWI_", stamp, ".rds")),
#  compress = "xz"
#)
# ---- Rolling latest version 
#saveRDS(
#  loco_save,
#  file = file.path(out_dir, "rf_loco_results_latest.rds")
#)

cat("✅ LOCO RF results saved to:\n", out_dir, "\n")
# ====================================================



# ===============================================================
# Standalone LOCO RF VI Plot, using existing Excel table
# ===============================================================
### load libraries
# ==============================================================
# Make VI plots — Build-Up Versions
# ==============================================================

library(readxl)
library(dplyr)
library(ggplot2)
library(grid)

# ---- File path
folder_path <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/Ch_01_Figures/40m_spacing"
file_name   <- "VI_Table_Final.xlsx"
excel_path  <- file.path(folder_path, file_name)

# ---- Read data
vi_summary_all <- readxl::read_excel(excel_path)

# ---- Standardize column names if needed
vi_summary_all <- vi_summary_all %>%
  dplyr::rename(
    Variable    = dplyr::any_of(c("Variable", "variable")),
    Mean_IncMSE = dplyr::any_of(c("Mean_IncMSE", "mean_incMSE", "Mean_pctIncMSE")),
    SD_IncMSE   = dplyr::any_of(c("SD_IncMSE", "sd_incMSE", "SD_pctIncMSE"))
  )

# ---- Variable order
var_order_pub <- c(
  "sev_num",
  "pptz_JJA",
  "tmeanz_JJA",
  "swez_Apr",
  "cwd_5yr_zscore_08",
  "cwd_30yrAvg_TT",
  "twi"
)

# ---- YSFs to include in plot layout
ysf_all <- c(5, 10, 15, 20)

# ---- Colors
#ysf_colors <- c(
#  "1"  = "aquamarine",
#  "5"  = "darkred",
#  "10" = "darkgreen",
#  "15" = "darkblue",
#  "20" = "darkviolet"
#)

## colorblind safe:
ysf_colors <- c(
  "5"  = "#D55E00",  # orange
  "10" = "#009E73",  # bluish green
  "15" = "#0072B2",  # blue
  "20" = "#CC79A7"   # reddish purple
)

# ---- Axis labels
var_labels <- c(
  "sev_num"              = "Fire Severity",
  "twi"                  = "Soil Moisture Capacity",
  "pptz_JJA"             = "Summer Total Precipitation",
  "tmeanz_JJA"           = "Summer Mean Temperature",
  "swez_Apr"             = "April Snow Water Equivalent",
  "cwd_5yr_zscore_08"    = "Post-Fire Drought",
  "cwd_30yrAvg_TT"       = "Underlying Conditions"
)

# ---- Positioning
pd <- position_dodge(width = 0.6)

# Helper function: build VI plot while preserving layout
make_vi_build <- function(
    ysf_to_show = ysf_all,
    show_legend = TRUE
) {
  
  vi_plot_df <- vi_summary_all %>%
    dplyr::filter(YSF %in% ysf_all) %>%
    dplyr::mutate(
      Variable = factor(Variable, levels = rev(var_order_pub)),
      YSF      = factor(YSF, levels = sort(ysf_all, decreasing = TRUE)),
      
      alpha_group = ifelse(as.numeric(as.character(YSF)) %in% ysf_to_show, 1, 0),
      color_group = ifelse(as.numeric(as.character(YSF)) %in% ysf_to_show, as.character(YSF), NA)
    )
  
  p <- ggplot(
    vi_plot_df,
    aes(
      x = Variable,
      y = Mean_IncMSE,
      group = YSF
    )
  ) +
    geom_errorbar(
      aes(
        ymin  = Mean_IncMSE - SD_IncMSE,
        ymax  = Mean_IncMSE + SD_IncMSE,
        color = color_group,
        alpha = alpha_group
      ),
      width = 0.25,
      position = pd,
      linewidth = 1.2
    ) +
    geom_point(
      aes(
        color = color_group,
        alpha = alpha_group
      ),
      size = 4.0,
      stroke = 1.2,
      position = pd
    ) +
    coord_flip() +
    scale_color_manual(
      values = ysf_colors,
      limits = as.character(sort(ysf_all, decreasing = TRUE)),
      drop = FALSE,
      na.value = NA,
      name = "Years Since Fire"
    ) +
    scale_alpha_identity() +
    scale_x_discrete(
      labels = var_labels,
      drop = FALSE
    ) +
    labs(
      title = "Variable Importance Across Years Since Fire",
      x = "Predictor Variable",
      y = "Mean %IncMSE (± SD)"
    ) +
    guides(color = guide_legend(reverse = TRUE)) +
    theme_minimal() +
    theme(
      text = element_text(family = "Aptos"),
      plot.title = element_text(size = 25, face = "bold"),
      axis.title.x = element_text(size = 16, margin = ggplot2::margin(t = 10), face = "bold"),
      axis.title.y = element_text(size = 16, margin = ggplot2::margin(r = 10), face = "bold"),
      axis.text.x = element_text(size = 14, margin = ggplot2::margin(t = 8), face = "bold"),
      axis.text.y = element_text(size = 14, margin = ggplot2::margin(r = 8), face = "bold"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 13),
      legend.key.size = unit(1.4, "lines"),
      legend.key.width = unit(1.0, "lines"),
      legend.spacing.y = unit(0.2, "cm"),
      panel.grid.major.x = element_line(
        color = "gray85",
        linewidth = 0.5
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      
      axis.line.x = element_line(color = "black", linewidth = 0.7),
      axis.line.y = element_line(color = "black", linewidth = 0.7),
      
      axis.ticks.x = element_line(color = "black", linewidth = 0.6),
      axis.ticks.y = element_line(color = "black", linewidth = 0.6),
      legend.position = "right"
    )
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  p
}


# Build-up plots for presentation
p_vi_5 <- make_vi_build(
  ysf_to_show = c(5),
  show_legend = TRUE
)

p_vi_5_10 <- make_vi_build(
  ysf_to_show = c(5, 10),
  show_legend = TRUE
)

p_vi_15_20 <- make_vi_build(
  ysf_to_show = c(15, 20),
  show_legend = TRUE
)

p_vi_all <- make_vi_build(
  ysf_to_show = c(5, 10, 15, 20),
  show_legend = TRUE
)

# ---- Print plots
p_vi_5
p_vi_5_10
p_vi_15_20
p_vi_all
# ==============================================================




# LOAD LOCO RF RESULTS (skip re-running LOCO models)
# ==============================================================

## Edit this pathway to load the correct file. 
#rds_path_upload <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/Ch01_RF_LOCO/rf_loco_results_latest.rds"
rds_path_upload <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/Ch01_RF_LOCO/rf_loco_results_cwd5_TWI_20260220_205732.rds"
loco_save <- readRDS(rds_path_upload)

# Organize uploaded data - save variables
rf_loco_results   <- loco_save$rf_loco_results
perf_loco_all     <- loco_save$perf_loco_all
perf_loco_summary <- loco_save$perf_loco_summary

str(perf_loco_summary)

# Double check - need envblock_results_named
stopifnot(exists("envblock_results_named"))
# Reconstruct YSF set from loaded env-blocking results
ysf_set <- as.numeric(names(envblock_results_named))
stopifnot(length(ysf_set) > 0)

# ========================================================
# Single-variable PDPs for LOCO models -----------------------------------------
# ========================================================
# ---- controls ----
pdp_grid_n   <- 60
pdp_q_lo     <- 0.01
pdp_q_hi     <- 0.99
pdp_row_cap  <- 2500
pdp_seed     <- 999

cwd_col  <- "cwd_5yr_zscore_08"
sev_mode <- "sev_num"   # LOCO RFs ALWAYS use numeric severity

# ---- predictors to plot PDP for ----
base_vars_effect <- c(
  "tmeanz_JJA",
  #"twi",
  #"veg_climate_index_08",
  "pptz_JJA",
  "swez_Apr",
  "cwd_30yrAvg_TT",
  cwd_col
)

# Include severity PDP only when numeric
vars_effect <- c("sev_num", base_vars_effect)

# ---- FUNCTION: robust grid helper (numeric) ----
.master_grid_cont <- function(x, n, q_lo, q_hi) {
  xnum <- suppressWarnings(as.numeric(x))
  if (!length(xnum) || length(unique(xnum[is.finite(xnum)])) < 2) return(numeric(0))
  ql <- suppressWarnings(stats::quantile(xnum, probs = q_lo, na.rm = TRUE))
  qh <- suppressWarnings(stats::quantile(xnum, probs = q_hi, na.rm = TRUE))
  if (!is.finite(ql) || !is.finite(qh) || ql == qh) return(numeric(0))
  seq(ql, qh, length.out = n)
}
# ---- FUNCTION: get fold training data (mirror RF prep, dynamic factor vs numeric) ----
.get_train_data_for_fold <- function(ysf_key, fold_name,
                                     rf_numeric_cols,
                                     rf_factor_cols = character(0),
                                     row_cap = pdp_row_cap,
                                     seed = pdp_seed) {
  kf_val <- sub("^Fold_", "", fold_name)
  
  df_all <- envblock_results_named[[as.character(ysf_key)]]$data_with_folds_env
  df2    <- df_all %>% dplyr::mutate(.fold_chr = as.character(fold))
  
  train_raw <- df2 %>% dplyr::filter(.fold_chr != kf_val)
  if (!nrow(train_raw)) train_raw <- df2
  
  train <- train_raw %>%
    dplyr::select(-.fold_chr) %>%
    sf::st_drop_geometry()
  
  # factor handling
  for (col in intersect(rf_factor_cols, names(train))) {
    train[[col]] <- factor(as.character(train[[col]]))
  }
  
  # numeric imputation (median from full training fold)
  for (col in intersect(rf_numeric_cols, names(train))) {
    train[[col]] <- suppressWarnings(as.numeric(train[[col]]))
    med <- suppressWarnings(stats::median(train[[col]], na.rm = TRUE))
    if (is.na(med)) med <- 0
    train[[col]][is.na(train[[col]])] <- med
  }
  
  # speed cap (after imputation)
  if (is.finite(row_cap) && nrow(train) > row_cap) {
    set.seed(seed)
    train <- train[sample.int(nrow(train), row_cap), , drop = FALSE]
  }
  
  train
}
# ---- FUNCTION: brute-force PD at a grid value ----
.predict_grid_cont <- function(model, train_df, var, grid_x) {
  if (is.null(model) || !is.data.frame(train_df) || !length(grid_x)) return(tibble())
  purrr::map_dfr(grid_x, function(xx) {
    newdat <- train_df
    newdat[[var]] <- xx
    preds <- try(predict(model, newdata = newdat), silent = TRUE)
    if (inherits(preds, "try-error")) return(tibble())
    pnum <- as.numeric(preds)
    tibble::tibble(
      x         = xx,
      fold_mean = mean(pnum, na.rm = TRUE),
      fold_sd   = stats::sd(pnum, na.rm = TRUE)
    )
  })
}
# ---- FUNCTION: build CV PDPs across YSFs ----
build_cv_pdp <- function(ysf_set,
                         rf_results_envblock,
                         envblock_results_named,
                         cwd_col,
                         vars_effect,
                         pdp_grid_n,
                         sev_mode = "sev_num") {
  
  rf_factor_cols  <- character(0)
  rf_numeric_cols <- c(
    "sev_num",
    #"twi",
    "pptz_JJA",
    "tmeanz_JJA",
    "swez_Apr",
    #"veg_climate_index_08",
    "cwd_30yrAvg_TT"
    #cwd_col
  )
  
  set.seed(pdp_seed)
  
  purrr::map_dfr(as.character(ysf_set), function(ysf_key) {
    
    rf_obj <- rf_results_envblock[[as.character(ysf_key)]]
    if (is.null(rf_obj) || is.null(rf_obj$fold_models)) return(tibble())
    
    df_all <- envblock_results_named[[as.character(ysf_key)]]$data_with_folds_env %>%
      sf::st_drop_geometry()
    
    # grids per variable
    master_grids <- lapply(vars_effect, function(vv) {
      if (!vv %in% names(df_all)) return(NULL)
      .master_grid_cont(df_all[[vv]], pdp_grid_n, pdp_q_lo, pdp_q_hi)
    })
    names(master_grids) <- vars_effect
    
    fold_names <- names(rf_obj$fold_models)
    
    fold_curves <- purrr::map_dfr(fold_names, function(fnm) {
      
      rf_model <- rf_obj$fold_models[[fnm]]$model
      if (is.null(rf_model)) return(tibble())
      
      tr_df <- .get_train_data_for_fold(
        ysf_key,
        fnm,
        rf_numeric_cols,
        rf_factor_cols,
        pdp_row_cap,
        pdp_seed
      )
      
      vars_here <- intersect(vars_effect, names(tr_df))
      
      purrr::map_dfr(vars_here, function(vv) {
        grid_x <- master_grids[[vv]]
        if (is.null(grid_x) || !length(grid_x)) return(tibble())
        out <- .predict_grid_cont(rf_model, tr_df, vv, grid_x)
        if (!nrow(out)) return(tibble())
        out %>%
          dplyr::mutate(
            YSF      = as.numeric(ysf_key),
            variable = vv,
            Fold     = fnm
          )
      })
    })
    
    if (!nrow(fold_curves)) return(tibble())
    
    fold_curves %>%
      dplyr::group_by(YSF, variable, x) %>%
      dplyr::summarise(
        y_mean = mean(fold_mean, na.rm = TRUE),
        
        # between-fold SD
        sd_between = if (dplyr::n() > 1)
          stats::sd(fold_mean, na.rm = TRUE) else 0,
        
        mean_within_var = {
          v <- fold_sd^2
          v <- v[is.finite(v)]
          if (length(v)) mean(v) else 0
        },
        
        # total SD (keep for diagnostics)
        y_sd = sqrt(pmax(0, sd_between^2 + mean_within_var)),
        
        n_folds = dplyr::n(),
        
        # 95% Confidence Interval for ribbon
        se_mean = ifelse(n_folds > 1,
                         sd_between / sqrt(n_folds),
                         NA_real_),
        t_crit  = stats::qt(0.975, df = pmax(n_folds - 1, 1)),
        ci_halfwidth = t_crit * se_mean,
        ci_lower = y_mean - ci_halfwidth,
        ci_upper = y_mean + ci_halfwidth,
        
        .groups = "drop"
      ) %>%
      dplyr::arrange(variable, x)
  })
}
# ========================================================

## Build PDPs; create pdp_loco_cv ========================
message("\n==============================")
message("Building LOCO single-variable PDPs")
message("==============================")

pdp_loco_results <- purrr::imap(
  rf_loco_results,
  function(m, model_name) {
    
    message("\n--- PDPs for LOCO model: ", model_name, " ---")
    
    build_cv_pdp(
      ysf_set                = ysf_set,
      rf_results_envblock    = m$results,
      envblock_results_named = envblock_results_named,
      cwd_col                = cwd_col,
      vars_effect            = intersect(vars_effect, m$rf_vars),
      pdp_grid_n             = pdp_grid_n,
      sev_mode               = "sev_num"
    ) %>%
      dplyr::mutate(model = model_name)
  }
)
pdp_loco_cv <- dplyr::bind_rows(pdp_loco_results)
stopifnot(nrow(pdp_loco_cv) > 0)

# ---- factor prep 
model_levels <- names(rf_loco_results)

pdp_loco_cv <- pdp_loco_cv %>%
  dplyr::mutate(
    YSF      = factor(YSF, levels = sort(unique(YSF), decreasing = TRUE)),
    model    = factor(model, levels = model_levels),
    variable = factor(variable, levels = vars_effect)
  )

# ---- QA
pdp_loco_cv %>%
  dplyr::group_by(model, YSF, variable) %>%
  dplyr::summarise(
    n_x = dplyr::n(),
    min_x = min(x, na.rm = TRUE),
    max_x = max(x, na.rm = TRUE),
    min_sd = min(y_sd, na.rm = TRUE),
    max_sd = max(y_sd, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print(n = Inf)

# PLOTTING =========

## height 1100, w 900

vars_to_plot <- c(
  #"tmeanz_JJA"
  #"swez_Apr"
  #"pptz_JJA"
  #"cwd_5yr_zscore_08"
  #"sev_num"
  #"twi"
  #"cwd_30yrAvg_TT"
  #"elevation_m"
)

vars_to_plot <- intersect(vars_to_plot, unique(pdp_loco_cv$variable))
# ========================================================

# ---------------------------------------------------------
# LABELS
# ---------------------------------------------------------
var_labels <- c(
  "sev_num"              = "Fire Severity",
  "tmeanz_JJA"           = "Temperature",
  "twi"                  = "Topographic Wetness Index",
  "veg_climate_index_08" = "Vegetation–Climate Index",
  "cwd_5yr_zscore_08"    = "5-Year Post-Fire CWD",
  "swez_Apr"             = "Snow Water Equivalent",
  "pptz_JJA"             = "Precipitation",
  "elevation_m"          = "Elevation",
  "cwd_30yrAvg_TT"       = "CWD 30-Year Normal"
)

model_labels <- c(
  full = "Full model",
  drop_sev_num = "Drop: Fire severity",
  drop_twi = "Drop: TWI",
  drop_pptz_JJA = "Drop: Summer precip",
  drop_tmeanz_JJA = "Drop: Summer temp",
  drop_swez_Apr = "Drop: April SWE",
  drop_cwd_5yr_zscore_08 = "Drop: Postfire CWD"
)

# ---------------------------------------------------------
# Prepare PDP plotting data (REMOVE EMPTY PANELS)
# ---------------------------------------------------------
pdp_plot_df <- pdp_loco_cv %>%
  dplyr::filter(variable %in% vars_to_plot) %>%
  dplyr::mutate(
    variable = factor(variable, levels = vars_to_plot),
    model_label = dplyr::recode(
      as.character(model),
      !!!model_labels,
      .default = as.character(model)
    )
  ) %>%
  dplyr::group_by(model_label, variable) %>%
  dplyr::filter(
    dplyr::n_distinct(x) > 1,
    any(is.finite(y_mean))
  ) %>%
  dplyr::ungroup()

# =========================================================
# R² ANNOTATION TABLE (ROW ACROSS TOP)
# =========================================================
# panel ranges for proper annotation placement
panel_ranges <- pdp_plot_df %>%
  dplyr::group_by(model_label, variable) %>%
  dplyr::summarise(
    x_min = min(x, na.rm = TRUE),
    x_max = max(x, na.rm = TRUE),
    .groups = "drop"
  )

annot_df <- perf_loco_all %>%
  dplyr::mutate(
    model = as.character(model),
    model_label = dplyr::recode(
      model,
      !!!model_labels,
      .default = model
    ),
    YSF = factor(YSF, levels = levels(pdp_plot_df$YSF)),
    label = sprintf("R² = %.2f", Mean_PseudoR2_weighted)
  ) %>%
  tidyr::crossing(variable = vars_to_plot) %>%
  dplyr::semi_join(
    pdp_plot_df %>% dplyr::distinct(model_label, variable),
    by = c("model_label", "variable")
  ) %>%
  dplyr::left_join(panel_ranges,
                   by = c("model_label", "variable")) %>%
  dplyr::group_by(model_label, variable) %>%
  dplyr::arrange(desc(YSF), .by_group = TRUE) %>%
  dplyr::mutate(
    x_pos = x_min + (x_max - x_min) *
      seq(0.08, 0.92, length.out = dplyr::n()),
    y_pos = Inf
  ) %>%
  dplyr::ungroup()

# =========================================================
# FINAL PDP PLOT
# =========================================================
ggplot(
  pdp_plot_df,
  aes(x = x, y = y_mean, color = YSF, group = YSF)
) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  
  # ---- 95% CI ribbon 
geom_ribbon(
  aes(
    ymin = ci_lower,
    ymax = ci_upper,
    fill = YSF
  ),
  alpha = 0.18,
  color = NA
) +
  
  geom_line(linewidth = 1.1) +
  
  # ---- R² labels across top 
geom_text(
  data = annot_df,
  aes(
    x = x_pos,
    y = y_pos,
    label = label,
    color = YSF
  ),
  hjust = 0.5,
  vjust = 1.2,
  size = 3.2,
  inherit.aes = FALSE,
  show.legend = FALSE
) +
  
  facet_grid(
    model_label ~ variable,
    scales = "free_x",
    drop = TRUE,
    labeller = labeller(
      variable = as_labeller(var_labels),
      model_label = label_value
    )
  ) +
  
  labs(
    title = "LOCO Partial Dependence of Predictors on ΔNDVI",
    subtitle = "Model: CWD5+CWD30, TWI",
    x = NULL,
    y = "Model-Predicted ΔNDVI",
    caption = "Shaded ribbons show 95% confidence intervals across CV folds"
  ) +
  
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(size = 19, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(size = 13),
    strip.text.y = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    axis.title.y = element_text(size = 14, face = "bold"),
    panel.grid.minor = element_blank() 
  )
################################################################################


### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #####
### Blurb... looking at 5YSF precip. ===========================================
ppt_low  <- -0.5
ppt_high <-  0.25
target_ysf <- 5
pixel_ysf5_sf <- data_long_reduced_sf %>%
  dplyr::filter(years_since_fire == target_ysf) %>%
  dplyr::filter(!is.na(pptz_JJA)) %>%
  dplyr::arrange(pixel_ID, year) %>%   
  dplyr::group_by(pixel_ID) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

## Flag the range of pptz to look at -----------
pixel_ysf5_sf <- pixel_ysf5_sf %>%
  dplyr::mutate(
    ppt_flag = dplyr::case_when(
      pptz_JJA >= ppt_low & pptz_JJA <= ppt_high ~ "Suspicious range",
      TRUE ~ "Other"
    )
  )
table(pixel_ysf5_sf$ppt_flag)
summary(pixel_ysf5_sf$pptz_JJA)

# --- Read shapefiles & match CRS to pixels ------------
huc12_sf <- sf::st_read(
  "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/ArcGIS/BMWC_Fire_Data/HUC12_BMWA/HUC12_BobMarshall.shp",
  quiet = TRUE
)
streams_sf <- sf::st_read(
  "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/ArcGIS/BMWC_Fire_Data/NHDPlus_V2_Streamlines/BMWA_streams_widened.shp",
  quiet = TRUE
)
huc12_bg  <- huc12_sf  %>% sf::st_transform(sf::st_crs(pixel_ysf5_sf))
streams_bg <- streams_sf %>% sf::st_transform(sf::st_crs(pixel_ysf5_sf))

## Plot the map -------
gg_ppt_problem <- ggplot() +
  
  # HUC12 background
  geom_sf(
    data = huc12_bg,
    fill = "grey95",
    color = "grey70",
    linewidth = 0.2
  ) +
  
  # all pixels (faint)
  geom_sf(
    data = pixel_ysf5_sf,
    color = "grey60",
    size = 0.4,
    alpha = 0.25
  ) +
  
  # highlighted suspicious pixels
  geom_sf(
    data = dplyr::filter(pixel_ysf5_sf, ppt_flag == "Suspicious range"),
    color = "red",
    size = 0.9,
    alpha = 0.9
  ) +
  
  coord_sf() +
  labs(
    title = "YSF = 5 Pixels with Suspicious pptz Range",
    subtitle = paste0("pptz_JJA in [", ppt_low, ", ", ppt_high, "]"),
    color = NULL
  ) +
  theme_bw()

gg_ppt_problem
## ------------------------------------------------
### What about fire severity / year? ========
suspicious_pixels <- pixel_ysf5_sf %>%
  dplyr::filter(ppt_flag == "Suspicious range")

# ------------------------------------------------------------
# ref_year summary
# ------------------------------------------------------------
year_summary <- suspicious_pixels %>%
  sf::st_drop_geometry() %>%
  dplyr::count(ref_year, name = "n_pixels") %>%
  dplyr::mutate(
    prop = n_pixels / sum(n_pixels),
    pct  = prop * 100
  ) %>%
  dplyr::arrange(ref_year)

print(year_summary)

# Quick diagnostics
year_range <- range(suspicious_pixels$ref_year, na.rm = TRUE)
year_stats <- summary(suspicious_pixels$ref_year)

print(year_range)
print(year_stats)

# ------------------------------------------------------------
# severity group proportions
# ------------------------------------------------------------

sev_summary <- suspicious_pixels %>%
  dplyr::count(sev_group, name = "n_pixels") %>%
  dplyr::mutate(
    prop = n_pixels / sum(n_pixels),
    pct  = prop * 100
  ) %>%
  dplyr::arrange(desc(n_pixels))

print(sev_summary)
### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #####


#          -------------------------------------------                         #
#########      POST-RF: PARTIAL DEPENDENCE PLOTTING     ########################
#          -------------------------------------------                         #

# Single-Variable PDP's:
# - If RF uses sev_group: DOES NOT create a severity PDP (by design)
# - If RF uses sev_num:  DOES create a sev_num PDP 
# ==============================================================================

# Load libraries ----------------
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(sf)
  library(ggplot2)
})
# -------------------------------

# ---- controls ----
pdp_grid_n      <- 60
pdp_q_lo        <- 0.01
pdp_q_hi        <- 0.99
pdp_row_cap     <- 2500
pdp_seed        <- 999

# ---- resolve CWD column ----
cwd_zscore_var <- "cwd_5yr_zscore_08"
cwd_col <- if (exists("cwd_zscore_var")) as.character(cwd_zscore_var) else "cwd_5yr_zscore_08"

# ---- detect which RF version is currently in rf_results_envblock ----
# Assumes: already run ONLY ONE RF option, and stored it in rf_results_envblock.
.detect_rf_severity_mode <- function(rf_results_envblock) {
  # Grab first available fold model + its training data column names.
  ysf_keys <- names(rf_results_envblock)
  ysf_keys <- ysf_keys[!is.na(ysf_keys) & ysf_keys != ""]
  if (!length(ysf_keys)) return(NA_character_)
  
  rf_obj <- rf_results_envblock[[ysf_keys[1]]]
  if (is.null(rf_obj$fold_results) || !length(rf_obj$fold_results)) return(NA_character_)
  fnm <- names(rf_obj$fold_results)[1]
  mdl <- rf_obj$fold_results[[fnm]]$model
  
  # randomForest stores the training data variables in the call/formula,
  # but safest is to inspect model$call$formula or model$terms when present.
  vars <- character(0)
  
  # 1) terms (preferred)
  trm <- tryCatch(stats::terms(mdl), error = function(e) NULL)
  if (!is.null(trm)) {
    vars <- attr(trm, "term.labels")
  }
  
  # 2) fallback to call$formula string parsing
  if (!length(vars) && !is.null(mdl$call$formula)) {
    f <- as.character(mdl$call$formula)
    # f is like c("~", "delta_ndvi_min", "sev_group + twi + ...")
    if (length(f) >= 3) {
      rhs <- f[3]
      vars <- unlist(strsplit(gsub("\\s+", "", rhs), "\\+"))
    }
  }
  
  if ("sev_num" %in% vars) return("sev_num")
  if ("sev_group" %in% vars) return("sev_group")
  NA_character_
}

sev_mode <- .detect_rf_severity_mode(rf_results_envblock)

if (is.na(sev_mode)) {
  message("PDP: Could not detect RF severity mode (sev_group vs sev_num). ",
          "Defaulting to NO severity PDP.")
  sev_mode <- "unknown"
} else {
  message("PDP: Detected RF severity mode = ", sev_mode)
}

# ---- predictors to plot PDP for (dynamic) ----
base_vars_effect <- c(
  #"cwd_30yrAvg_TT",
  #"twi",
  #"veg_climate_index_08",
  "tmeanz_JJA",
  "pptz_JJA",
  "swez_Apr",
  cwd_col
)

# *Only* include severity PDP if the severity variable is numeric (sev_num)
vars_effect <- if (identical(sev_mode, "sev_num")) c("sev_num", base_vars_effect) else base_vars_effect

# ---- robust grid helper (numeric) ----
.master_grid_cont <- function(x, n = pdp_grid_n, q_lo = pdp_q_lo, q_hi = pdp_q_hi) {
  xnum <- suppressWarnings(as.numeric(x))
  if (!length(xnum) || length(unique(xnum[is.finite(xnum)])) < 2) return(numeric(0))
  ql <- suppressWarnings(stats::quantile(xnum, probs = q_lo, na.rm = TRUE))
  qh <- suppressWarnings(stats::quantile(xnum, probs = q_hi, na.rm = TRUE))
  if (!is.finite(ql) || !is.finite(qh) || ql == qh) return(numeric(0))
  seq(ql, qh, length.out = n)
}

# ---- get fold training data (mirror RF prep, dynamic factor vs numeric) ----
.get_train_data_for_fold <- function(ysf_key, fold_name,
                                     rf_numeric_cols,
                                     rf_factor_cols = character(0),
                                     row_cap = pdp_row_cap, seed = pdp_seed) {
  kf_val <- sub("^Fold_", "", fold_name)
  
  df_all <- envblock_results_named[[as.character(ysf_key)]]$data_with_folds_env
  df2 <- df_all %>% mutate(.fold_chr = as.character(fold))
  
  train_raw <- df2 %>% filter(.fold_chr != kf_val)
  if (!nrow(train_raw)) train_raw <- df2
  
  train <- train_raw %>%
    select(-.fold_chr) %>%
    sf::st_drop_geometry()
  
  # factor handling
  for (col in intersect(rf_factor_cols, names(train))) {
    train[[col]] <- factor(as.character(train[[col]]))
  }
  
  # numeric imputation (median)
  for (col in intersect(rf_numeric_cols, names(train))) {
    train[[col]] <- suppressWarnings(as.numeric(train[[col]]))
    med <- suppressWarnings(stats::median(train[[col]], na.rm = TRUE))
    if (is.na(med)) med <- 0
    train[[col]][is.na(train[[col]])] <- med
  }
  
  # speed cap
  if (is.finite(row_cap) && nrow(train) > row_cap) {
    set.seed(seed)
    train <- train[sample.int(nrow(train), row_cap), , drop = FALSE]
  }
  
  train
}

# ---- brute-force PD at a given grid x: return fold_mean and fold_sd ----
.predict_grid_cont <- function(model, train_df, var, grid_x) {
  if (is.null(model) || !is.data.frame(train_df) || !length(grid_x)) return(tibble())
  map_dfr(grid_x, function(xx) {
    newdat <- train_df
    newdat[[var]] <- xx
    preds <- try(predict(model, newdata = newdat), silent = TRUE)
    if (inherits(preds, "try-error")) return(tibble())
    pnum <- as.numeric(preds)
    tibble(
      x         = xx,
      fold_mean = mean(pnum, na.rm = TRUE),
      fold_sd   = stats::sd(pnum, na.rm = TRUE)
    )
  })
}

# ---- FUNCTION: build CV PDPs across YSFs ----
build_cv_pdp <- function(ysf_set,
                         rf_results_envblock,
                         envblock_results_named,
                         cwd_col,
                         vars_effect,
                         pdp_grid_n = pdp_grid_n,
                         sev_mode = c("sev_group", "sev_num", "unknown")) {
  
  sev_mode <- match.arg(sev_mode)
  
  if (sev_mode == "sev_num") {
    rf_factor_cols  <- character(0)
    rf_numeric_cols <- c(
      "sev_num",
      "twi",
      "pptz_JJA",
      #"tmaxz_JJA",
      "tmeanz_JJA",
      "swez_Apr",
      #"veg_climate_index_08",
      cwd_col
    )
  } else {
    # sev_group or unknown: mirror your sev_group RF structure
    rf_factor_cols  <- c("sev_group")
    rf_numeric_cols <- c(
      "twi",
      "pptz_JJA",
      #"tmaxz_JJA",
      "tmeanz_JJA",
      "swez_Apr",
      #"veg_climate_index_08",
      cwd_col
    )
  }
  
  set.seed(pdp_seed)
  
  map_dfr(as.character(ysf_set), function(ysf_key) {
    
    rf_obj <- rf_results_envblock[[as.character(ysf_key)]]
    if (is.null(rf_obj) || is.null(rf_obj$fold_results)) return(tibble())
    
    df_all <- envblock_results_named[[as.character(ysf_key)]]$data_with_folds_env %>%
      sf::st_drop_geometry()
    
    # build grids per variable (per YSF)
    master_grids <- list()
    for (vv in vars_effect) {
      if (!vv %in% names(df_all)) next
      master_grids[[vv]] <- .master_grid_cont(df_all[[vv]], n = pdp_grid_n)
    }
    
    fold_names <- names(rf_obj$fold_results)
    
    fold_curves <- map_dfr(fold_names, function(fnm) {
      rf_model <- rf_obj$fold_results[[fnm]]$model
      if (is.null(rf_model)) return(tibble())
      
      tr_df <- .get_train_data_for_fold(
        ysf_key   = ysf_key,
        fold_name = fnm,
        rf_numeric_cols = rf_numeric_cols,
        rf_factor_cols  = rf_factor_cols,
        row_cap = pdp_row_cap,
        seed    = pdp_seed
      )
      
      vars_here <- intersect(vars_effect, names(tr_df))
      
      map_dfr(vars_here, function(vv) {
        grid_x <- master_grids[[vv]]
        if (is.null(grid_x) || !length(grid_x)) return(tibble())
        out <- .predict_grid_cont(rf_model, tr_df, vv, grid_x)
        if (!nrow(out)) return(tibble())
        out %>%
          mutate(
            YSF      = as.numeric(ysf_key),
            variable = vv,
            Fold     = fnm
          )
      })
    })
    
    if (!nrow(fold_curves)) return(tibble())
    
    fold_curves %>%
      group_by(YSF, variable, x) %>%
      dplyr::summarise(
        y_mean = mean(fold_mean, na.rm = TRUE),
        
        # between-fold SD
        sd_between = if (dplyr::n() > 1)
          stats::sd(fold_mean, na.rm = TRUE) else 0,
        
        mean_within_var = {
          v <- fold_sd^2
          v <- v[is.finite(v)]
          if (length(v)) mean(v) else 0
        },
        
        # keep total SD if you still want it
        y_sd = sqrt(pmax(0, sd_between^2 + mean_within_var)),
        
        n_folds = dplyr::n(),
        
        # -----------------------------
        # 95% CONFIDENCE INTERVAL
        # -----------------------------
        t_crit = stats::qt(0.975, df = pmax(n_folds - 1, 1)),
        se_mean = sd_between / sqrt(pmax(n_folds, 1)),
        ci_halfwidth = t_crit * se_mean,
        ci_lower = y_mean - ci_halfwidth,
        ci_upper = y_mean + ci_halfwidth,
        
        .groups = "drop"
      )
  })
}

# =========================
# BUILD SINGLE-VARIABLE PDPs
# =========================
pdp_cv <- build_cv_pdp(
  ysf_set               = ysf_set,
  rf_results_envblock    = rf_results_envblock,
  envblock_results_named = envblock_results_named,
  cwd_col               = cwd_col,
  vars_effect           = vars_effect,
  pdp_grid_n            = pdp_grid_n,
  sev_mode              = if (sev_mode %in% c("sev_group", "sev_num")) sev_mode else "unknown"
)

# sanity checks
pdp_cv %>%
  group_by(YSF, variable) %>%
  summarise(
    n_x = dplyr::n(),
    min_x = min(x, na.rm = TRUE),
    max_x = max(x, na.rm = TRUE),
    min_sd = min(y_sd, na.rm = TRUE),
    max_sd = max(y_sd, na.rm = TRUE),
    .groups = "drop"
  ) %>% print(n = Inf)

# =========================
# PLOT PDPs
# =========================
var_labels <- c(
  "sev_num"              = "Fire Severity",
  "tmeanz_JJA"           = "Summer Mean Temperature",
  "twi"                  = "Topographic Wetness Index",
  #"veg_climate_index_08" = "Vegetation–Climate Index",
  "cwd_30yrAvg_TT"       = "CWD 30-Year Normal",
  "cwd_5yr_zscore_08"    = "5-Year Post-Fire CWD",
  "swez_Apr"             = "Snow Water Equivalent",
  "pptz_JJA"             = "Summer Precipitation"
)

if (!cwd_col %in% names(var_labels)) var_labels[[cwd_col]] <- cwd_col

pdp_cv <- pdp_cv %>%
  mutate(
    YSF      = factor(YSF, levels = sort(unique(YSF), decreasing = TRUE)),
    variable = factor(variable, levels = vars_effect)
  )

# Explicitly list the PDP's to print:
vars_to_plot <- c(
                  "sev_num", 
                  "tmeanz_JJA", 
                  "swez_Apr",
                  "pptz_JJA",
                  #"veg_climate_index_08",
                  "cwd_5yr_zscore_08"
                  #"cwd_30yrAvg_TT",
                  #"twi"
                  )   

# ---------------- Filter PDPs to selected variables ----------------
vars_to_plot <- intersect(vars_to_plot, unique(pdp_cv$variable))

pdp_plot_df <- pdp_cv %>%
  dplyr::filter(variable %in% vars_to_plot) %>%
  dplyr::mutate(variable = factor(variable, levels = vars_to_plot))

# ---------------- Plot PDPs ----------------
ggplot(pdp_plot_df, aes(x = x, y = y_mean, color = YSF, group = YSF)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_ribbon(
    aes(
      ymin = ci_lower,
      ymax = ci_upper,
      fill = YSF
    ),
    alpha = 0.18,
    color = NA
  ) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ variable,
             scales = "free_x",
             labeller = as_labeller(var_labels),
             ncol = 2) +   # set number of columns here
  labs(
    title = "Independent Effects on \u0394NDVI",
    x = NULL,
    y = "Model-Predicted \u0394NDVI"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(size = 19, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 15),
    strip.text = element_text(size = 14),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    axis.title.y = element_text(size = 15, face = "bold"),
    panel.grid.minor = element_blank()
  )
# ==============================================================================




###############################################################################
# 2-VARIABLE PDP'S (ALL PAIRWISE COMBOS) — fold-averaged 2D surfaces
# - Builds 2D PDPs for *all numeric predictor pairs* (including sev_num if present)
# - Skips sev_group automatically (factor) and skips any non-numeric predictors
# - Uses RED (low) -> WHITE (0) -> GREEN (high) fill, midpoint at 0
# - Returns:
#     1) pdp2_all: one long data frame with columns (var1,var2,YSF,x,y,yhat_mean,...)
#     2) pdp2_plots: named list of ggplot objects (one per pair)
###############################################################################

# Load libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(purrr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
})
# ------------------------------------------------------------------------------

# USER CONTROLS ================================================================
pdp2_grid_res <- 35      # grid points per axis (n x n total)
pdp2_row_cap  <- 3500    # optional speed cap per fold train set; set Inf to disable
pdp2_seed     <- 999

# If you sometimes swap cwd variable name, keep this:
cwd_col <- as.character(cwd_zscore_var)

# ---- define candidate predictor set (include sev_num if present) ----
rf_vars_candidate <- c(
  "sev_num",    # will be included if present + numeric
  #"twi",
  "pptz_JJA",
  "tmeanz_JJA",
  "swez_Apr",
  #"veg_climate_index_08",
  cwd_col
)

var_labels <- c(
  "sev_num"              = "Fire Severity",
  "tmeanz_JJA"           = "Summer Mean Temperature",
  #"twi"                  = "Topographic Wetness Index",
  #"veg_climate_index_08" = "Vegetation–Climate Index",
  "cwd_5yr_zscore_08"    = "5-Year Post-Fire CWD",
  "swez_Apr"             = "Snow Water Equivalent",
  "pptz_JJA"             = "Summer Precipitation"
)

# Function: get training data for a fold (from RF above) ------------------------
# IMPORTANT: train data here corresponds to "fold != kf"
get_train_df_for_fold <- function(ysf_key, fold_name, cwd_col, row_cap = Inf, seed = 999) {
  
  kf_val <- sub("^Fold_", "", fold_name)
  
  df_all <- envblock_results_named[[as.character(ysf_key)]]$data_with_folds_env %>%
    mutate(.fold_chr = as.character(fold)) %>%
    filter(!is.na(delta_ndvi_min)) %>%
    sf::st_drop_geometry()
  
  # training rows for this fold
  train <- df_all %>%
    filter(.fold_chr != kf_val) %>%
    select(-.fold_chr)
  
  if (!nrow(train)) train <- df_all %>% select(-.fold_chr)
  
  # keep ONLY response + candidate predictors (avoids extra cols)
  keep_cols <- intersect(c("delta_ndvi_min", rf_vars_candidate), names(train))
  train <- train[, keep_cols, drop = FALSE]
  
  # Coerce candidate predictors to numeric when possible (sev_num etc.)
  for (cc in intersect(rf_vars_candidate, names(train))) {
    train[[cc]] <- suppressWarnings(as.numeric(train[[cc]]))
  }
  
  # numeric median-impute
  for (cc in intersect(rf_vars_candidate, names(train))) {
    med <- suppressWarnings(stats::median(train[[cc]], na.rm = TRUE))
    if (is.na(med)) med <- 0
    train[[cc]][is.na(train[[cc]])] <- med
  }
  
  # speed cap (optional)
  if (is.finite(row_cap) && nrow(train) > row_cap) {
    set.seed(seed)
    train <- train[sample.int(nrow(train), row_cap), , drop = FALSE]
  }
  
  train
}

# Function: create master grid -------------------------------------------------
.master_grid_cont2 <- function(x, n = 35, q_lo = 0.01, q_hi = 0.99) {
  xnum <- suppressWarnings(as.numeric(x))
  xnum <- xnum[is.finite(xnum)]
  if (length(xnum) < 2 || length(unique(xnum)) < 2) return(numeric(0))
  
  ql <- suppressWarnings(stats::quantile(xnum, probs = q_lo, na.rm = TRUE))
  qh <- suppressWarnings(stats::quantile(xnum, probs = q_hi, na.rm = TRUE))
  if (!is.finite(ql) || !is.finite(qh) || ql == qh) return(numeric(0))
  
  seq(ql, qh, length.out = n)
}

# Function: Build fold-averaged 2D PDP surfaces across all YSF -----------------
build_cv_pdp2_fixedgrid <- function(ysf_set,
                                    rf_results_envblock,
                                    envblock_results_named,
                                    var1,
                                    var2,
                                    grid_n = 35,
                                    q_lo = 0.01,
                                    q_hi = 0.99,
                                    row_cap = Inf,
                                    seed = 999,
                                    cwd_col = cwd_col) {
  
  purrr::map_dfr(as.character(ysf_set), function(ysf_key) {
    
    rf_obj <- rf_results_envblock[[as.character(ysf_key)]]
    if (is.null(rf_obj) || is.null(rf_obj$fold_results)) return(tibble::tibble())
    
    # Build ONE grid per YSF from the full YSF dataset (not fold-specific)
    df_all <- envblock_results_named[[as.character(ysf_key)]]$data_with_folds_env %>%
      sf::st_drop_geometry()
    
    if (!(var1 %in% names(df_all)) || !(var2 %in% names(df_all))) return(tibble::tibble())
    
    gx <- .master_grid_cont2(df_all[[var1]], n = grid_n, q_lo = q_lo, q_hi = q_hi)
    gy <- .master_grid_cont2(df_all[[var2]], n = grid_n, q_lo = q_lo, q_hi = q_hi)
    if (!length(gx) || !length(gy)) return(tibble::tibble())
    
    grid_df <- tidyr::expand_grid(x = gx, y = gy)
    
    fold_names <- names(rf_obj$fold_results)
    
    fold_surfaces <- purrr::map_dfr(fold_names, function(fnm) {
      
      mdl <- rf_obj$fold_results[[fnm]]$model
      if (is.null(mdl)) return(tibble::tibble())
      
      train_df <- get_train_df_for_fold(
        ysf_key   = ysf_key,
        fold_name = fnm,
        cwd_col   = cwd_col,
        row_cap   = row_cap,
        seed      = seed
      )
      
      # needed columns exist?
      if (!(var1 %in% names(train_df)) || !(var2 %in% names(train_df))) return(tibble::tibble())
      
      # Evaluate PD by overwriting BOTH vars on the SAME grid
      pd_fold <- purrr::pmap_dfr(grid_df, function(x, y) {
        newdat <- train_df
        newdat[[var1]] <- x
        newdat[[var2]] <- y
        
        preds <- try(predict(mdl, newdata = newdat), silent = TRUE)
        if (inherits(preds, "try-error")) return(tibble::tibble())
        
        tibble::tibble(
          x = x,
          y = y,
          yhat = mean(as.numeric(preds), na.rm = TRUE)
        )
      })
      
      if (!nrow(pd_fold)) return(tibble::tibble())
      
      pd_fold %>%
        mutate(YSF = as.numeric(ysf_key), Fold = fnm)
    })
    
    if (!nrow(fold_surfaces)) return(tibble::tibble())
    
    fold_surfaces %>%
      group_by(YSF, x, y) %>%
      summarise(
        yhat_mean = mean(yhat, na.rm = TRUE),
        yhat_sd   = if (dplyr::n() > 1) stats::sd(yhat, na.rm = TRUE) else 0,
        n_folds   = dplyr::n(),
        .groups   = "drop"
      )
  })
}

# =============================================================================
# 1) Decide which variables to include (numeric only; includes sev_num if present)
# =============================================================================
# Use one YSF's full dataset to detect types robustly:
ysf0 <- as.character(ysf_set[1])
df0  <- envblock_results_named[[ysf0]]$data_with_folds_env %>% sf::st_drop_geometry()

# Keep only predictors that exist + are numeric-ish (after coercion)
vars_present <- intersect(rf_vars_candidate, names(df0))

is_numericish <- function(v) {
  x <- suppressWarnings(as.numeric(df0[[v]]))
  x <- x[is.finite(x)]
  length(x) >= 2 && length(unique(x)) >= 2
}

vars_2d <- vars_present[vapply(vars_present, is_numericish, logical(1))]

# (Safety) explicitly drop sev_group if it snuck in
vars_2d <- setdiff(vars_2d, "sev_group")

message("2D PDP: variables included = ", paste(vars_2d, collapse = ", "))

# All pairwise combos
pair_mat <- utils::combn(vars_2d, 2)

# =============================================================================
# 2) Build all 2-Var PDPs (long table) + ggplot objects list  (~10 hours...)
# =============================================================================
pdp2_all <- purrr::pmap_dfr(
  as.data.frame(t(pair_mat), stringsAsFactors = FALSE),
  function(V1, V2) {
    message("Building 2D PDP: ", V1, " x ", V2)
    
    out <- build_cv_pdp2_fixedgrid(
      ysf_set               = ysf_set,
      rf_results_envblock    = rf_results_envblock,
      envblock_results_named = envblock_results_named,
      var1                  = V1,
      var2                  = V2,
      grid_n                = pdp2_grid_res,
      q_lo                  = 0.01,
      q_hi                  = 0.99,
      row_cap               = pdp2_row_cap,
      seed                  = pdp2_seed,
      cwd_col               = cwd_col
    )
    
    if (!nrow(out)) return(tibble::tibble())
    
    out %>% mutate(var1 = V1, var2 = V2, .before = 1)
  }
)

# Round values to 2 sig digits
pdp2_all <- pdp2_all %>%
  dplyr::mutate(
    yhat_mean = signif(yhat_mean, 2),
    yhat_sd   = signif(yhat_sd, 2)
  )

# Save the output to disk!!!
out_dir <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/Ch_01_Figures/40m_spacing/Q3_FIGURES/ModelA_dropTWI_03052026"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")  # DateTime saved
saveRDS(
  pdp2_all,
  file = file.path(out_dir, paste0("pdp2_all_", stamp, ".rds"))
)

# also save a rolling 'latest' version
saveRDS(
  pdp2_all,
  file = file.path(out_dir, "pdp2_all_latest.rds")
)
cat("Saved pdp2_all to:\n", out_dir, "\n")

# =============================================================================
# 3) TO SKIP RE-MAKING 2-Var PDP'S: UPLOAD RDS FILE DIRECTLY!!                  
# =============================================================================
# path to rolling latest file:
rds_path <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/Ch_01_Figures/40m_spacing/Q3_FIGURES/2_Var_PDPs_Copy/pdp2_all_latest.rds"
pdp2_all <- readRDS(rds_path)

# Resolve CWD variable error:
cwd_zscore_var <- "cwd_5yr_zscore_08"
cwd_col <- as.character(cwd_zscore_var)

vars_2d <- c(
  "sev_num",
  "twi",
  "pptz_JJA",
  "tmeanz_JJA",
  "swez_Apr",
  "veg_climate_index_08",
  cwd_col
)

# quick sanity checks
str(pdp2_all)
dplyr::glimpse(pdp2_all)

# =============================================================================
# 4) Build PDP Legend                   
# =============================================================================

# labels: snap tiny values to 0, then 2 sig digits
.pdp_label <- function(x) {
  x[abs(x) < 1e-8] <- 0
  formatC(signif(x, 2), format = "fg", digits = 2)
}

# Uniform color scale settings (run after "pdp2_all" exists)
step_fill  <- 0.05  # color bin width
step_label <- 0.15  # label frequency

global_min <- min(pdp2_all$yhat_mean, na.rm = TRUE)
global_max <- max(pdp2_all$yhat_mean, na.rm = TRUE)

# data-driven rounding outward to nearest *fill* step (keeps bin edges aligned)
global_limits <- c(
  floor(global_min / step_fill) * step_fill,
  ceiling(global_max / step_fill) * step_fill
)

# Option A: add one extra bin on each side
# global_limits <- global_limits + c(-step_fill, step_fill)

# Option B: extend ONLY the lower limit by one fill-step
global_limits <- c(global_limits[1] - step_fill, global_limits[2])

# breaks for labels (every step_fill)
global_breaks <- seq(global_limits[1], global_limits[2], by = step_fill)
# label only at every step_label increment
.label_every <- function(brks, by = step_label, anchor = 0) {
  idx <- abs(((brks - anchor) / by) - round((brks - anchor) / by)) < 1e-8
  labs <- rep("", length(brks))
  labs[idx] <- .pdp_label(brks[idx])
  labs
}

# minor breaks define the color “steps” (every 0.05)
global_minor_breaks <- seq(global_limits[1], global_limits[2], by = step_fill)

fill_scale_global <- scale_fill_steps2(
  low = "#4B3621",
  mid = "beige",
  high = "darkgreen",
  midpoint = 0,
  limits = global_limits,
  breaks = global_breaks,   # bin edges every step_fill
  oob = scales::squish,
  name = "Predicted ΔNDVI",
  labels = .label_every(global_breaks, by = step_label, anchor = 0),
  guide = guide_coloursteps(
    title.position = "top",
    title.hjust = 0.5,
    title.theme = element_text(
      family = "Times New Roman",
      face = "bold",
      size = 12,
      margin = ggplot2::margin(b = 14)
    ),
    barheight = unit(4.5, "cm"),
    barwidth  = unit(0.6, "cm"),
    ticks = TRUE,
    label.position = "right",
    label.hjust = 0
  )
)

cat("Data range:", global_min, "to", global_max, "\n")
cat("Legend limits:", global_limits[1], "to", global_limits[2], "\n")
cat("Legend labeled breaks:", paste(global_breaks, collapse = ", "), "\n")
cat("Legend bin width:", step_fill, "\n")

# =============================================================================
# 5A) Build and Print Plots (Panels by Predictor Variable)               
# =============================================================================
# Create var_labels
var_labels <- c(
  "sev_num"              = "Fire Severity",
  "tmeanz_JJA"           = "Summer Mean Temperature",
  #"twi"                  = "Topographic Wetness Index",
  "cwd_5yr_zscore_08"    = "5-Year Post-Fire CWD",
  "swez_Apr"             = "Snow Water Equivalent",
  "pptz_JJA"             = "Summer Precipitation"
)

# Function: make one plot per variable pair ------------------------------------
.make_pdp2_plot <- function(df_pair, fill_scale = fill_scale_global, var_labels = var_labels) {
  stopifnot(nrow(df_pair) > 0)
  v1 <- df_pair$var1[1]
  v2 <- df_pair$var2[1]
  
  # pretty labels (fallback to raw name if missing)
  v1_lab <- if (!is.null(var_labels[[v1]])) var_labels[[v1]] else v1
  v2_lab <- if (!is.null(var_labels[[v2]])) var_labels[[v2]] else v2
  
  ggplot(df_pair, aes(x = x, y = y, fill = yhat_mean)) +
    geom_raster() +
    facet_wrap(
      ~ YSF,
      scales = "free",
      labeller = labeller(YSF = function(x) paste(x, "YSF"))
    ) +
    geom_contour(aes(z = yhat_mean), color = "black", alpha = 0.35) +
    fill_scale +
    labs(
      title = paste0("Partial Dependence: ", v1_lab, " × ", v2_lab),
      x = v1_lab,
      y = v2_lab
    ) +
    theme_minimal() +
    theme(
      text = element_text(family = "Times New Roman"),
      plot.title = element_text(size = 19, face = "bold", hjust = 0.5),
      strip.text = element_text(size = 15, face = "bold"),
      legend.title = element_text(size = 13, face = "bold", hjust = 0.5, vjust = 0.8),
      legend.text  = element_text(size = 11),
      axis.title.y = element_text(size = 15, face = "bold"),
      axis.title.x = element_text(size = 15, face = "bold"),
      axis.text.x = element_text(size = 11, face = "plain"),
      axis.text.y = element_text(size = 11, face = "plain"),
      panel.grid.minor = element_blank()
    )
}

# Build named list of plots (one per pair) -------------------------------------
pdp2_plots <- list()

if (nrow(pdp2_all)) {
  pair_keys <- pdp2_all %>%
    dplyr::distinct(var1, var2) %>%
    dplyr::mutate(key = paste0(var1, "__x__", var2))
  
  for (i in seq_len(nrow(pair_keys))) {
    v1 <- pair_keys$var1[i]
    v2 <- pair_keys$var2[i]
    key <- pair_keys$key[i]
    
    df_pair <- pdp2_all %>% dplyr::filter(var1 == v1, var2 == v2)
    if (!nrow(df_pair)) next
    
    pdp2_plots[[key]] <- tryCatch(
      .make_pdp2_plot(df_pair, fill_scale = fill_scale_global, var_labels = var_labels),
      error = function(e) {
        message("Plot failed for ", key, ": ", conditionMessage(e))
        NULL
      }
    )
  }
}

# Print all plots (one after another) OR access them by name: -------------------
# ---- Option A: Print everything (will page through in RStudio plot pane) 
for (nm in names(pdp2_plots)) print(pdp2_plots[[nm]])

# --- Option B: Print one specific pair 
# available plot names for copy/paste

#[1] "sev_num__x__twi"                           
#[2] "sev_num__x__pptz_JJA"                      
#[3] "sev_num__x__tmeanz_JJA"                    
#[4] "sev_num__x__swez_Apr"                      
#[5] "sev_num__x__veg_climate_index_08"          
#[6] "sev_num__x__cwd_5yr_zscore_08"             
#[7] "twi__x__pptz_JJA"                          
#[8] "twi__x__tmeanz_JJA"                        
#[9] "twi__x__swez_Apr"                          
#[10] "twi__x__veg_climate_index_08"              
#[11] "twi__x__cwd_5yr_zscore_08"                 
#[12] "pptz_JJA__x__tmeanz_JJA"                   
#[13] "pptz_JJA__x__swez_Apr"                     
#[14] "pptz_JJA__x__veg_climate_index_08"         
#[15] "pptz_JJA__x__cwd_5yr_zscore_08"            
#[16] "tmeanz_JJA__x__swez_Apr"                   
#[17] "tmeanz_JJA__x__veg_climate_index_08"       
#[18] "tmeanz_JJA__x__cwd_5yr_zscore_08"          
#[19] "swez_Apr__x__veg_climate_index_08"         
#[20] "swez_Apr__x__cwd_5yr_zscore_08"            
#[21] "veg_climate_index_08__x__cwd_5yr_zscore_08"

print(pdp2_plots[["sev_num__x__swez_Apr"]])
print(pdp2_plots[["sev_num__x__tmeanz_JJA"]])
print(pdp2_plots[["sev_num__x__pptz_JJA"]])
print(pdp2_plots[["sev_num__x__twi"]])

# =============================================================================
# 5B) 6-panel PDP figure from selected pair labels
# =============================================================================

# ---- CHOOSE SIX PLOTS (copy from Option B list above) 
pairs_6 <- c(
  "sev_num__x__swez_Apr",
  "sev_num__x__tmeanz_JJA",
  "sev_num__x__pptz_JJA",
  "sev_num__x__twi",
  "sev_num__x__veg_climate_index_08",
  "sev_num__x__cwd_5yr_zscore_08"
)

var_labels <- c(
  "sev_num"              = "Fire Severity",
  "tmeanz_JJA"           = "Summer Mean Temperature",
  "twi"                  = "Topographic Wetness Index",
  "veg_climate_index_08" = "Vegetation–Climate Index",
  "cwd_5yr_zscore_08"    = "5-Year Post-Fire CWD",
  "swez_Apr"             = "Snow Water Equivalent",
  "pptz_JJA"             = "Summer Precipitation"
)

# YSF groups to print 
ysf_groups <- sort(unique(pdp2_all$YSF))

# Legend formatting for panel ----
fill_scale_pdp_panel <- scale_fill_steps2(
  low = "#4B3621",
  mid = "beige",
  high = "darkgreen",
  midpoint = 0,
  limits = global_limits,
  breaks = global_breaks,
  labels = .label_every(global_breaks, by = step_label, anchor = 0),
  oob = scales::squish,
  name = "Predicted ΔNDVI",
  
  guide = guide_coloursteps(
    title.position = "top",
    title.hjust = 0.5,
    
    title.theme = element_text(
      family = "Times New Roman",
      face   = "bold",
      size   = 22,      # <- panel-specific legend title size
      margin = ggplot2::margin(b = 10)
    ),
    
    label.theme = element_text(
      family = "Times New Roman",
      size   = 15
    ),
    
    barheight = unit(4.5, "cm"),
    barwidth  = unit(0.6, "cm"),
    ticks = TRUE,
    label.position = "right"
  )
)

# Function: make single panel for a YSF -----
.make_pdp2_plot_single_ysf <- function(df_pair, ysf,
                                       fill_scale = fill_scale_pdp_panel,
                                       var_labels = var_labels,
                                       show_legend = FALSE,
                                       show_xlab   = TRUE) {
  
  df_pair <- df_pair %>% dplyr::filter(YSF == ysf)
  stopifnot(nrow(df_pair) > 0)
  
  v1 <- df_pair$var1[1]
  v2 <- df_pair$var2[1]
  
  v1_lab <- if (!is.null(var_labels[[v1]])) var_labels[[v1]] else v1
  v2_lab <- if (!is.null(var_labels[[v2]])) var_labels[[v2]] else v2
  
  ggplot(df_pair, aes(x = x, y = y, fill = yhat_mean)) +
    geom_raster() +
    geom_contour(aes(z = yhat_mean), color = "black", alpha = 0.35) +
    fill_scale +
    labs(
      # title = paste0(v1_lab, " × ", v2_lab), # List both vars in plot title
      title = paste0(v2_lab),   # removes "Fire Severity × ..."
      #subtitle = "subtitle placeholder - individual plots",            
      x = if (show_xlab) v1_lab else NULL,
      y = v2_lab
    ) +
    theme_minimal() +
    theme(
      text = element_text(family = "Times New Roman"),
      
      # titles & text
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
      
      axis.title.x = element_text(
        size = 15,
        face = "bold",
        margin = ggplot2::margin(t = 10)
        ),
      
      axis.title.y = element_text(
        size = 15, face = "bold",
        margin = ggplot2::margin(r = 10)
        ),
      
      axis.text    = element_text(size = 15),
      
      # legend
      legend.position = if (show_legend) "right" else "none",
      
      # margins + square panel
      plot.margin = ggplot2::margin(t = 4, r = 7, b = 1, l = 6, unit = "pt"),
      panel.grid.minor = element_blank(),
      aspect.ratio = 1
    )
}

# PANEL LOOP: --------------
pdp_panels_by_ysf <- list()

for (ysf in ysf_groups) {
  
  plots_6 <- vector("list", length(pairs_6))
  
  for (i in seq_along(pairs_6)) {
    key <- pairs_6[i]
    parts <- strsplit(key, "__x__", fixed = TRUE)[[1]]
    v1 <- parts[1]
    v2 <- parts[2]
    
    df_pair <- pdp2_all %>%
      dplyr::filter(var1 == v1, var2 == v2)
    
    # bottom row gets x-axis labels (plots 4–6 in 2×3 layout)
    show_xlab <- i > 3
    
    plots_6[[i]] <- .make_pdp2_plot_single_ysf(
      df_pair,
      ysf = ysf,
      fill_scale = fill_scale_pdp_panel,
      var_labels = var_labels,
      show_legend = (i == 1),   # legend only once
      show_xlab   = show_xlab
    )
  }
  
  # ---- extract shared legend ----
  shared_legend <- get_legend(plots_6[[1]])
  
  # ---- remove legend from all individual PDPs ----
  plots_6 <- lapply(plots_6, function(p) p + theme(legend.position = "none"))
  
  panel <- plot_grid(
    plotlist = plots_6,
    ncol = 3,
    align = "hv",
    axis = "tblr"
  )
  
  title <- ggdraw() +
    draw_label(
      paste0("Partial Dependence Plots (", ysf, " YSF)"),
      fontfamily = "Times New Roman",
      fontface = "bold",
      size = 30,
      hjust = 0.5
    )
  
  subtitle <- ggdraw() +
    draw_label(
      "Variable interactions with fire severity",
      fontfamily = "Times New Roman",
      fontface = "plain",
      size = 24,
      hjust = 0.5,
      y = 0.9  # vertical adjustment
    )
  
  pdp_panels_by_ysf[[paste0("YSF_", ysf)]] <-
    plot_grid(
      NULL,  # top padding
      plot_grid(
        NULL,  # left padding
        plot_grid(
          title,
          subtitle,
          plot_grid(panel, shared_legend, rel_widths = c(1, 0.14)),
          ncol = 1,
          rel_heights = c(0.10, 0.06, 1)
        ),
        NULL,  # right padding  
        ncol = 3,
        rel_widths = c(0.015, 1, 0.015)  # small right buffer
      ),
      NULL,  # bottom padding
      ncol = 1,
      rel_heights = c(0.02, 1, 0.02)
    )
  
}

# ---- PRINT ---------
print(pdp_panels_by_ysf[["YSF_5"]])  # When exporting, width = 1850
print(pdp_panels_by_ysf[["YSF_10"]])
print(pdp_panels_by_ysf[["YSF_15"]])
print(pdp_panels_by_ysf[["YSF_20"]])

# =============================================================================
# 5C) Build and Print Triangle Plots (Panels by YSF)
# =============================================================================

pdp_target_ysf <- 10

var_labels <- c(
  "sev_num"              = "Fire Severity",
  "tmeanz_JJA"           = "Summer Mean Temperature",
  "twi"                  = "Topographic Wetness Index",
  "veg_climate_index_08" = "Vegetation–Climate Index",
  "cwd_5yr_zscore_08"    = "5-Year Post-Fire CWD",
  "swez_Apr"             = "Snow Water Equivalent",
  "pptz_JJA"             = "Summer Precipitation"
)

# Use your variable order (vars_2d) if it exists; otherwise infer from pdp2_all
vars_order <- if (exists("vars_2d")) vars_2d else {
  sort(unique(c(pdp2_all$var1, pdp2_all$var2)))
}

# Build a triangular (half-matrix) facet plot for ONE YSF
pdp_tri_df <- pdp2_all %>%
  dplyr::filter(YSF == pdp_target_ysf) %>%
  dplyr::mutate(
    var1 = factor(var1, levels = vars_order),
    var2 = factor(var2, levels = vars_order)
  ) %>%
  # keep ONLY one triangle so each pair appears once:
  # (var2 "below" var1 in the matrix)
  dplyr::filter(as.integer(var2) > as.integer(var1))

# Optional: if you want the opposite triangle, swap the inequality:
# dplyr::filter(as.integer(var2) < as.integer(var1))

# Plot  ---------------------------
ggplot(pdp_tri_df, aes(x = x, y = y, fill = yhat_mean)) +
  geom_raster() +
  geom_contour(aes(z = yhat_mean), color = "black", alpha = 0.35) +
  fill_scale_global +
  facet_grid(
    rows = vars(var2),
    cols = vars(var1),
    scales = "free",
    drop = TRUE,
    switch = "both",
    labeller = ggplot2::labeller(
      var1 = ggplot2::as_labeller(var_labels),
      var2 = ggplot2::as_labeller(var_labels)
    )
  ) +
  labs(
    title = paste0("Interactive Effects of Environmental Predictors on ΔNDVI"),
    subtitle = paste0(pdp_target_ysf, " Years Post-Fire"),
    x = NULL,
    y = NULL
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 15, hjust = 0.5),
    
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x.bottom = element_text(size = 11, face = "bold"),
    strip.text.y.left   = element_text(size = 11, face = "bold"),
    
    panel.spacing = unit(0.4, "lines"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    
    legend.title = element_text(
      margin = ggplot2::margin(t = 0, r = 0, b = 6, l = 0), # space below legend title
                   size = 12, face = "bold", hjust = 0.5),
    legend.text  = element_text(size = 10),
    panel.grid.minor = element_blank()
  )
###############################################################################


































################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################











################################################################################
############ ---- Just for kicks... RF with NO spatial CV -----  ###############
################################################################################
# ==============================================================
# RANDOM (NON-SPATIAL) K-FOLD CV: PIXEL -> FOLD ASSIGNMENT + RF 
# ==============================================================

# ---- Helper: make random folds for one YSF (optionally stratified) ----
make_random_folds_one_ysf <- function(target_ysf,
                                      k_folds      = k,
                                      seed_folds   = seed,
                                      stratify     = TRUE,      # set FALSE for pure random
                                      group_col    = "sev_group") {
  stopifnot(exists("data_long_reduced_sf"))
  set.seed(seed_folds + as.integer(target_ysf))
  
  # Subset
  pts <- data_long_reduced_sf %>%
    dplyr::filter(years_since_fire == target_ysf,
                  !is.na(delta_ndvi_min)) %>%
    # drop geometry; folds are random so we don't need coordinates
    sf::st_drop_geometry()
  
  if (nrow(pts) < k_folds) {
    stop("Not enough rows for YSF = ", target_ysf, " with k_folds = ", k_folds, ".")
  }
  
  # Build fold vector
  if (isTRUE(stratify) && group_col %in% names(pts)) {
    # Stratify by sev_group to keep proportions similar across folds
    fold_vec <- integer(nrow(pts))
    grp <- factor(pts[[group_col]])
    for (g in levels(grp)) {
      idx <- which(grp == g)
      # balanced assignment per group
      fold_ids <- rep(1:k_folds, length.out = length(idx))
      fold_ids <- sample(fold_ids, length(idx), replace = FALSE)
      fold_vec[idx] <- fold_ids
    }
  } else {
    # Pure random (no stratification)
    fold_vec <- rep(1:k_folds, length.out = nrow(pts))
    fold_vec <- sample(fold_vec, nrow(pts), replace = FALSE)
  }
  
  pts$fold <- factor(fold_vec, levels = 1:k_folds)
  
  list(
    YSF                    = target_ysf,
    data_with_folds_random = pts
  )
}

# ---- Convenience: name the results by YSF (same pattern as above) ----
randomcv_by_ysf <- function(randomcv_results) {
  ys <- vapply(randomcv_results, function(x) x$YSF, numeric(1))
  names(randomcv_results) <- as.character(ys)
  randomcv_results
}

# ---- RF function that expects the random-folds object (no spatial CV) ----
#rf_run_one_ysf_random <- function(randomcv_result,
                                  trees   = ntree,
                                  seed_rf = 23,
                                  cwd_var = cwd_zscore_var) {
  stopifnot(is.list(randomcv_result),
            "data_with_folds_random" %in% names(randomcv_result))
  
  # Reuse your existing RF code path as closely as possible
  target_ysf <- randomcv_result$YSF
  data_with_folds <- randomcv_result$data_with_folds_random
  k_folds <- length(levels(data_with_folds$fold))
  
  message("\n========================")
  message("RF CV with RANDOM Folds - YSF = ", target_ysf)
  message("========================")
  
  rf_factor_cols  <- c("sev_group")
  rf_numeric_cols <- c("twi","pptz_JJA","tmeanz_JJA","swez_Apr","veg_climate_index_08", cwd_var)
  
  # ---------- Pearson correlation screen (numeric covariates only) ----------
  num_df <- data_with_folds[, rf_numeric_cols, drop = FALSE]
  for (cc in names(num_df)) num_df[[cc]] <- suppressWarnings(as.numeric(num_df[[cc]]))
  cm <- tryCatch(stats::cor(num_df, use = "pairwise.complete.obs", method = "pearson"),
                 error = function(e) NULL)
  if (is.null(cm)) {
    corr_summary <- tibble::tibble(
      YSF = target_ysf, Var1 = character(), Var2 = character(),
      Pearson_r = numeric(), Abs_r = numeric()
    )
  } else {
    nm  <- colnames(cm)
    idx <- which(upper.tri(cm), arr.ind = TRUE)
    corr_summary <- tibble::tibble(
      YSF = target_ysf,
      Var1 = nm[idx[, 1]],
      Var2 = nm[idx[, 2]],
      Pearson_r = cm[idx],
      Abs_r = abs(cm[idx])
    ) %>%
      dplyr::filter(Abs_r > 0.7) %>%
      dplyr::arrange(dplyr::desc(Abs_r))
  }
  # --------------------------------------------------------------------------
  
  fold_results    <- vector("list", length = k_folds)
  names(fold_results) <- paste0("Fold_", levels(data_with_folds$fold))
  importance_list <- vector("list", length = k_folds)
  names(importance_list) <- names(fold_results)
  
  for (kf in levels(data_with_folds$fold)) {
    train_data <- data_with_folds %>% dplyr::filter(fold != kf)
    test_data  <- data_with_folds %>% dplyr::filter(fold == kf)
    
    # Factor level harmonization
    for (col in intersect(rf_factor_cols, names(train_data))) {
      tr <- factor(train_data[[col]]); te <- factor(test_data[[col]])
      all_lvls <- union(levels(tr), levels(te))
      train_data[[col]] <- factor(as.character(train_data[[col]]), levels = all_lvls)
      test_data[[col]]  <- factor(as.character(test_data[[col]]),  levels = all_lvls)
    }
    # Median impute
    for (col in intersect(rf_numeric_cols, names(train_data))) {
      med <- suppressWarnings(stats::median(train_data[[col]], na.rm = TRUE)); if (is.na(med)) med <- 0
      train_data[[col]][is.na(train_data[[col]])] <- med
      test_data[[col]][is.na(test_data[[col]])]  <- med
    }
    
    set.seed(seed_rf + as.integer(as.character(kf)))
    rf_formula <- as.formula(paste(
      "delta_ndvi_min ~",
      paste(c("sev_group","twi","pptz_JJA","tmeanz_JJA","swez_Apr","veg_climate_index_08",cwd_var),
            collapse = " + ")
    ))
    
    rf_model <- randomForest::randomForest(
      rf_formula,
      data = train_data,
      ntree = trees,
      importance = TRUE
    )
    
    predictions <- predict(rf_model, newdata = test_data)
    observed    <- test_data$delta_ndvi_min
    N           <- length(observed)
    mse         <- if (N) mean((predictions - observed)^2) else NA_real_
    rmse        <- if (is.na(mse)) NA_real_ else sqrt(mse)
    denom       <- sum((observed - mean(observed))^2)
    pseudo_r2   <- if (denom == 0) NA_real_ else 1 - (sum((observed - predictions)^2) / denom)
    
    fold_results[[paste0("Fold_", kf)]] <- list(
      N = N, mse = mse, rmse = rmse, pseudo_r2 = pseudo_r2,
      model = rf_model, predictions = predictions, observed = observed
    )
    importance_list[[paste0("Fold_", kf)]] <- randomForest::importance(rf_model, type = 1)
  }
  
  performance_df <- dplyr::bind_rows(lapply(names(fold_results), function(nm) {
    fr <- fold_results[[nm]]
    tibble::tibble(
      YSF = target_ysf,
      Fold = nm,
      N = fr$N,
      MSE = fr$mse,
      RMSE = fr$rmse,
      Pseudo_R2 = fr$pseudo_r2
    )
  }))
  
  mean_rmse  <- mean(performance_df$RMSE, na.rm = TRUE)
  mean_r2    <- mean(performance_df$Pseudo_R2, na.rm = TRUE)
  w          <- ifelse(is.na(performance_df$N), 0, performance_df$N)
  wmean_rmse <- stats::weighted.mean(performance_df$RMSE, w, na.rm = TRUE)
  wmean_r2   <- stats::weighted.mean(performance_df$Pseudo_R2, w, na.rm = TRUE)
  
  perf_row <- tibble::tibble(
    YSF = target_ysf,
    Mean_RMSE_unweighted     = mean_rmse,
    Mean_PseudoR2_unweighted = mean_r2,
    Mean_RMSE_weighted       = wmean_rmse,
    Mean_PseudoR2_weighted   = wmean_r2
  )
  
  importance_df <- dplyr::bind_rows(
    lapply(importance_list, function(mat) {
      if (is.null(mat)) return(NULL)
      df_out <- as.data.frame(mat)
      df_out$Variable <- rownames(df_out)
      df_out
    }),
    .id = "Fold"
  )
  
  if (!nrow(importance_df)) {
    vi_summary <- tibble::tibble(YSF = target_ysf,
                                 Variable = character(),
                                 Mean_IncMSE = numeric(),
                                 SD_IncMSE = numeric(),
                                 W_Mean_IncMSE = numeric())
  } else {
    fold_sizes <- performance_df %>% dplyr::select(Fold, N)
    importance_df <- importance_df %>% dplyr::left_join(fold_sizes, by = "Fold")
    vi_summary <- importance_df %>%
      dplyr::group_by(Variable) %>%
      dplyr::summarise(
        Mean_IncMSE   = mean(`%IncMSE`, na.rm = TRUE),
        SD_IncMSE     = stats::sd(`%IncMSE`, na.rm = TRUE),
        W_Mean_IncMSE = stats::weighted.mean(`%IncMSE`, w = ifelse(is.na(N), 0, N), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(W_Mean_IncMSE)) %>%
      dplyr::mutate(YSF = target_ysf, .before = 1)
  }
  
  list(
    YSF            = target_ysf,
    corr_summary   = corr_summary,
    perf_row       = perf_row,
    perf_folds     = performance_df,
    vi_summary     = vi_summary,
    fold_results   = fold_results,
    importance     = importance_list
  )
}

# ------------------ Build random folds for all YSFs ------------------
#randomcv_results <- lapply(
#  ysf_set,
#  function(y) make_random_folds_one_ysf(
#    target_ysf   = y,
#    k_folds      = k,
#    seed_folds   = seed,     # keep comparable to env-block runs
#    stratify     = TRUE      # set FALSE to remove sev_group stratification
#  )
#)
#randomcv_results_named <- randomcv_by_ysf(randomcv_results)

# ------------------ Run RF with RANDOM folds for all YSFs ------------------
#rf_random_results <- lapply(
#  as.character(ysf_set),
#  function(ysf_key) {
#    rf_run_one_ysf_random(
#      randomcv_result = randomcv_results_named[[ysf_key]],
#      trees   = ntree,
#      seed_rf = 23,
#      cwd_var = cwd_zscore_var
#    )
#  }
#)

# ------------------ Collect outputs (random CV) ------------------
#perf_rows_all_random  <- dplyr::bind_rows(lapply(rf_random_results, `[[`, "perf_row"))
#perf_folds_all_random <- dplyr::bind_rows(lapply(rf_random_results, `[[`, "perf_folds"))
#vi_summary_all_random <- dplyr::bind_rows(lapply(rf_random_results, `[[`, "vi_summary"))
#corrs_all_random      <- dplyr::bind_rows(lapply(rf_random_results, `[[`, "corr_summary"))

cat("\n================ RANDOM CV — Combined Performance (by YSF) ================\n")
print(perf_rows_all_random)

cat("\n================ RANDOM CV — Fold-level Performance (all YSF) ================\n")
print(perf_folds_all_random, n = Inf)

cat("\n================ RANDOM CV — Variable Importance Summary (all YSF) ================\n")
print(vi_summary_all_random, n = Inf)

cat("\n================ RANDOM CV — High Covariate Correlations (|r| > 0.7) (by YSF) ================\n")
if (nrow(corrs_all_random) == 0) {
  print(tibble::tibble(YSF = numeric(), Var1 = character(), Var2 = character(),
                       Pearson_r = numeric(), Abs_r = numeric()))
} else {
  print(corrs_all_random, n = Inf)
}

# ------------------ (Optional) quick side-by-side summary ------------------
#if (exists("perf_rows_all")) {
  comp_perf <- dplyr::full_join(
    perf_rows_all  %>% dplyr::mutate(CV = "Env-block"),
    perf_rows_all_random %>% dplyr::mutate(CV = "Random"),
    by = c("YSF","Mean_RMSE_unweighted","Mean_PseudoR2_unweighted",
           "Mean_RMSE_weighted","Mean_PseudoR2_weighted","CV")
  ) %>%
    dplyr::arrange(YSF, CV)
  cat("\n================ COMPARISON — Env-block vs Random (by YSF) ================\n")
  print(comp_perf, n = Inf)
}

# ===== Variable Importance Plot — RANDOM CV =====

stopifnot(exists("vi_summary_all_random"))

# Order variables by the earliest YSF in your ysf_set (same convention)
first_ysf_random <- ysf_set[1]

var_levels_first_random <- vi_summary_all_random %>%
  dplyr::filter(YSF == first_ysf_random) %>%
  dplyr::arrange(dplyr::desc(Mean_IncMSE)) %>%
  dplyr::pull(Variable)

# YSF levels so legend shows 5 first; then reversed via guides(reverse=TRUE)
ysf_levels_plot_random <- sort(unique(vi_summary_all_random$YSF), decreasing = TRUE)

vi_plot_df_random <- vi_summary_all_random %>%
  dplyr::mutate(
    Variable = factor(Variable, levels = rev(var_levels_first_random)),
    YSF      = factor(YSF, levels = ysf_levels_plot_random)
  )

pd <- position_dodge(width = 0.6)

# Reuse your colors/labels (subset to levels present, just in case)
ysf_levels_present <- levels(vi_plot_df_random$YSF)
ysf_colors_random  <- ysf_colors[names(ysf_colors) %in% ysf_levels_present]

var_labels_random <- var_labels[names(var_labels) %in% levels(vi_plot_df_random$Variable)]
# fall back to names if some labels missing
if (length(var_labels_random) == 0) var_labels_random <- levels(vi_plot_df_random$Variable)

ggplot(vi_plot_df_random,
       aes(x = Variable, y = Mean_IncMSE, color = YSF, group = YSF)) +
  geom_errorbar(aes(ymin = Mean_IncMSE - SD_IncMSE,
                    ymax = Mean_IncMSE + SD_IncMSE),
                width = 0.2, position = pd, linewidth = 0.7) +
  geom_point(size = 3, position = pd) +
  coord_flip() +
  scale_color_manual(values = ysf_colors_random, name = "Years Since Fire") +
  scale_x_discrete(labels = var_labels_random) +
  labs(
    title = "Variable Importance — Random (Non-Spatial) CV",
    x = "Predictor Variable",
    y = "Mean %IncMSE (± SD)"
  ) +
  guides(color = guide_legend(reverse = TRUE)) +
  theme_minimal(base_family = "Times New Roman") +
  theme(
    text             = element_text(family = "Times New Roman"),
    plot.title       = element_text(
      size = 25, face = "bold"),
    plot.subtitle    = element_text(
      size = 12, face = "italic"),
    axis.title.x     = element_text(
      size = 16, margin = ggplot2::margin(t = 10), face = "bold"),
    axis.title.y = element_text(
      size = 15,face = "bold", margin = ggplot2::margin(r = 10)),
    axis.text.x      = element_text(size = 14, face = "bold", margin = ggplot2::margin(t = 8)),
    axis.text.y      = element_text(size = 14, face = "bold", margin = ggplot2::margin(l = 15)),
    legend.title     = element_text(size = 15, face = "bold"),
    legend.text      = element_text(size = 14),
    legend.key.size  = unit(1.4, "lines"),
    legend.key.width = unit(1.0, "lines"),
    legend.spacing.y = unit(0.2, "cm"),
    panel.grid.minor = element_blank()
  )
# ==============================================================================

