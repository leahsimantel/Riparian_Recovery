### Widening stream shapefile

# Load required libraries
library(sf)       # For spatial operations
library(dplyr)    # For data manipulation
library(tidyr)    # For handling missing data
library(ggplot2)  # For visualization

# Read in CSV files
main_table <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/NHDPlus_Flowlines_EntireBMWA_NoLakes.csv"  
wetted_width <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/EntireBMWA_WettedWidth.csv"  

# Read the CSV files into data frames
main_table <- read.csv(main_table, stringsAsFactors = FALSE)
wetted_width <- read.csv(wetted_width, stringsAsFactors = FALSE)

# Count the number of NA values for TotDASqKM and WettedWidth.
num_na_totda <- sum(is.na(main_table$TotDASqKM))
print(num_na_totda)  # 0

num_na_ww <- sum(is.na(wetted_width$WETTEDWIDTH))
print(num_na_ww)  # 247

# Thus, need to extrapolate wetted width for 247 segments.

# Keep only the 'COMID', ' and 'TotDASqKM' columns in 'main_table'
#main_table <- main_table %>%
#  select(COMID, TotDASqKM)

# Delete two observations that are irrelevant
# COMID 22967406 and 22968560
main_table <- main_table %>%
  filter(!COMID %in% c(22967406, 22968560))
length(main_table$COMID)  #981

# Merge the wetted_width data into main_table based on the COMID column
merged <- merge(main_table, wetted_width, by = "COMID", all.x = TRUE)


# Load the shapefile 
stream_network <- st_read("C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/ArcGIS/BMWC_Fire_Data/NHDPlus_V2_Streamlines/NHDFlowline_Network_clipped_NoLakes.shp")
# plot(stream_network)

# Keep only the relevant columns: COMID, WETTEDWIDTH, TotDASqKM, and geometry-related fields
stream_network <- stream_network %>%
  select(COMID, geometry)

# Merge the attributes from the 'merged' data frame with the spatial data
stream_network <- stream_network %>%
  left_join(merged, by = "COMID")  #join using COMID


# View the updated dataframe
print(st_geometry_type(stream_network))  # Check geometry type
print(names(stream_network))  # Check remaining column names


################################################################################
###############    Handle missing WETTEDWIDTH Values    ########################
################################################################################

# Identify COMIDs in main_table that are missing from wetted_width
missing_comids <- setdiff(main_table$COMID, wetted_width$COMID)

cat("COMIDs missing from wetted_width:\n", paste(missing_comids, collapse = ", "), "\n")
# 2 missing values, meaning we need to extrapolate for 249 widths

# Predict missing WETTEDWIDTH based on a relationship with TotDASqKM (e.g., linear regression)
# Ensure TotDASqKM and WETTEDWIDTH are numeric
stream_network <- stream_network %>%
  mutate(
    WETTEDWIDTH = as.numeric(WETTEDWIDTH),
    TotDASqKM = as.numeric(TotDASqKM)
  )

# Build a linear model to predict WETTEDWIDTH from TotDASqKM
lm_model <- lm(WETTEDWIDTH ~ TotDASqKM, data = stream_network, na.action = na.exclude)

# Add a column to identify whether WETTEDWIDTH is observed or extrapolated
stream_network <- stream_network %>%
  mutate(
    WettedWidth_Source = ifelse(is.na(WETTEDWIDTH), "Extrapolated", "Observed")
  )

# Fill missing WETTEDWIDTH values with predicted values from the model
stream_network <- stream_network %>%
  mutate(
    WETTEDWIDTH = ifelse(
      is.na(WETTEDWIDTH),
      predict(lm_model, newdata = stream_network),
      WETTEDWIDTH
    )
  )

# Check the new column
table(stream_network$WettedWidth_Source)  # Count observed vs extrapolated

##############################################################################
#############       Creating the Polygons       ##############################
##############################################################################

# Add a column for half the wetted width
stream_network <- stream_network %>%
  mutate(half_width = WETTEDWIDTH / 2)
   
View(stream_network$half_width)  # Check!

# Drop Z and M components from the geometry
stream_network <- stream_network %>%
  mutate(geometry = st_zm(geometry, drop = TRUE))

# Get a summary of all unique FCODEs and their counts
fcode_summary <- stream_network %>%
  group_by(FCODE) %>%
  summarise(
    count = n(),
    stream_types = paste(unique(FTYPE), collapse = ", ") # Optional: to get corresponding FTYPEs
  ) %>%
  arrange(desc(count))  # Order by count, descending

print(fcode_summary)

# Add a column based on FCODE values
# If FCODE is 46003, categorize as 'Intermittent'; otherwise, map based on the NHD User Guide Fcode List
stream_network <- stream_network %>%
  mutate(flow_category = ifelse(
    FCODE == 46003, 
    "Intermittent", 
    case_when(
      FCODE == 46006 ~ "Perennial",  # Example mapping for FCODE 46006
      FCODE == 46007 ~ "Artificial",  # Example mapping for FCODE 46007
      # Add more FCODE mappings as needed based on the NHD User Guide
      TRUE ~ "Unknown"  # Default for unmapped FCODE values
    )
  ))

View(stream_network$flow_category)  # Check!

# Apply st_buffer using the half_width column
# Use row-wise operation for buffering with varying distances
stream_network_polygons <- st_sf(
  stream_network %>%
    rowwise() %>%
    mutate(geometry = st_buffer(geometry, dist = half_width)) %>%
    ungroup()
)

# Validate geometries after buffering
stream_network_polygons <- st_make_valid(stream_network_polygons)

# Reproject the spatial data to WGS84 (EPSG:4326)
stream_network_polygons <- st_transform(stream_network_polygons, crs = 4326)

###############################################################################
##########        Map the Results Using ggplot      ###########################
###############################################################################

# Plot the stream network polygons with color based on WETTEDWIDTH
ggplot(data = stream_network_polygons) +
  geom_sf(aes(fill = WETTEDWIDTH), color = NA) +  # Fill by WETTEDWIDTH, no border
  scale_fill_gradient(
    name = "Wetted Width",       # Legend title
    low = "lightgreen",           # Color for small widths
    high = "darkblue",           # Color for large widths
    na.value = "black"           # Color for missing values (if any)
  ) +
  labs(
    title = "Stream Network with Wetted Width Buffers",
    subtitle = "Wider polygons represent wider stream segments",
    caption = "Segments used for pixel harvest"
  ) +
  theme_minimal() +  # Apply a minimal theme
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title = element_blank()  # Remove axis titles
  )


##########################################################################
######### Create simplified CSV for output   #############################
##########################################################################
# 1) Create FCODE_cat 
stream_network_polygons <- stream_network_polygons %>%
  mutate(
    FCODE_cat = case_when(
      FCODE == 46006 ~ "Perennial",
      FCODE == 46007 ~ "Ephemeral",
      FCODE == 46003 ~ "Intermittent",
      TRUE           ~ "Unknown"
    )
  )

# Rename WETTEDWIDTH to WettedWidth_m
stream_network_polygons <- stream_network_polygons %>%
  dplyr::rename(WettedWidth_m = WETTEDWIDTH)

# 2) Create a simplified version with requested columns (in requested order)
stream_network_poly_simple <- stream_network_polygons %>%
  st_drop_geometry() %>%
  select(
    COMID,
    Shape..,
    GNIS_NAME,
    GNIS_ID,
    FCODE,
    FCODE_cat,
    StreamOrde,
    StreamLeve,
    LENGTHKM,
    WettedWidth_m,
    WettedWidth_Source,
    AreaSqKM,
    TotDASqKM
  )

# Quick check
glimpse(stream_network_poly_simple)

# Write simplified dataframe to CSV
write.csv(
  stream_network_poly_simple,
  file = "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/Ch01_data_drop/stream_network_poly_simple.csv",
  row.names = FALSE
)

##########################################################################
######### Writing the Shapefile   ########################################
##########################################################################

# Specify the output directory and filename
output_path <- "C:/Users/leahs/OneDrive/Documents/U_of_M/Masters_Project/ArcGIS/BMWC_Fire_Data/NHDPlus_V2_Streamlines/BMWA_streams_widened.shp"

# Write the shapefile
st_write(
  stream_network_polygons,
  dsn = output_path,
  delete_dsn = TRUE  # Overwrite if the file already exists
)

# Confirmation message
cat("Shapefile successfully written to:\n", output_path, "\n")


################################################################################
#######    tmap approach - interactive   #######################################
###############################################################################


# Allows you to interact with map features and make sure they make sense, prior 
# to writing the shapefile and uploading it in Arc. 


# Load tmap library
install.packages("tmap")
library(tmap)

# Set tmap mode to interactive
tmap_mode("view")

# Create an interactive map
tm_shape(stream_network_polygons) +
  tm_fill(
    col = "WETTEDWIDTH",                # Fill color by WETTEDWIDTH
    palette = "Blues",                  # Color palette
    title = "Wetted Width"              # Legend title
  ) +
  tm_layout(
    title = "Stream Network with Wetted Width Buffers",
    legend.outside = TRUE               # Place legend outside the map
  )
# Force tmap to open in a browser


# Create the tmap visualization
stream_map <- tm_shape(stream_network_polygons) +
  tm_fill(
    col = "WETTEDWIDTH",
    palette = "Blues",
    title = "Wetted Width"
  ) +
  tm_layout(title = "Stream Network with Wetted Width Buffers")

# Save the map to an HTML file
tmap_save(stream_map, "stream_network_map.html")  # Save to the specified path

# Open the saved map in your default browser for interactive viewing.
browseURL("stream_network_map.html")





