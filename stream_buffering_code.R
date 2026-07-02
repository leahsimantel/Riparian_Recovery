## Stream buffering.

# load sf library since we're using spatial data
library(sf) 

# Ensure that you have loaded the shapefile as an sf object using st_read 
# before applying the st_buffer function! USE 'st_read' command.
streams_subset_GEE <- st_read("C:/Users/leahs/OneDrive/Documents/U_of_M/WILD591_GEE/shapefiles/flowline_subset_dissolved/Flowline_PN17_GEE_subset_dissolved.shp")

# Create a 30-meter buffer around each stream line
streams_30_buffer <- st_buffer(streams_subset_GEE, dist = 30)

# check projection of other layers... fire polygons and stream buffer layer are both in WGS 84 projection!

fires_subset_GEE <- st_read("C:/Users/leahs/OneDrive/Documents/U_of_M/WILD591_GEE/shapefiles/Fires_2003_GEE_subset/Holbrook_Bartlett_2003Fires_Subset.shp")

# Save the new shapefile
st_write(streams_30_buffer, "C:/Users/leahs/OneDrive/Documents/U_of_M/WILD591_GEE/shapefiles/flowline_subset_dissolved/streams_30_buffer.shp")
