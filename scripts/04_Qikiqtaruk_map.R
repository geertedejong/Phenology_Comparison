### 4. Creating a map of the phenocam locations on Qikiqtaruk ###
### Geerte de Jong, geerte.de.jong@bioenv.gu.se ###
### Phenocam Project with Geerte de Jong, Joe Boyle, Maude Grenier & Elise Gallois ###
### Date: 5th of September 2024 ###

#### LOAD PACKAGES  #####
# Install necessary packages if not already installed
# install.packages("sf")
# install.packages("ggplot2")
install.packages("ggrepel")

# Load libraries
library(sf)       # For reading shapefiles
library(ggplot2)  # For plotting
library(ggrepel)  # For better label positioning

# The shapefile should contain the .shp, .shx, and .dbf files
shapefile_path <- "data/Ecological_classification_Herschel_Island.shp"  
shape_data <- st_read(shapefile_path)

# 2. Create a data frame of coordinates
coords <- data.frame(
  name = c("P1", "P2", "P3", "P4", "P5", "P6"),  # Names of the points
  lat = c(69.575599, 69.575698, 69.577812, 69.575149, 69.567507, 69.574829),  # Latitudes
  lon = c(-138.906, -138.90525, -138.91276, -138.89453, -138.86711, -138.86305) # Longitudes
)

# Convert coordinates to sf points (assuming WGS84 CRS, EPSG:4326)
coords_sf <- st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)

# 3. Check CRS compatibility and transform 
if (st_crs(shape_data) != st_crs(coords_sf)) {
  shape_data <- st_transform(shape_data, crs = st_crs(coords_sf))
}

# 3. Plot the shapefile and overlay the points
ggplot() +
  geom_sf(data = shape_data, fill = "lightgrey", color = "lightgrey") +  # Shapefile
  geom_sf(data = coords_sf,fill="red", size = 3) +      
  geom_text_repel(data = coords, aes(x = lon, y = lat, label = name), 
                  size = 3, color = "black", 
                  nudge_y = 0.002, box.padding = 0.3, max.overlaps = 10) +  # Adjust label positioning
  scale_color_viridis_c() +                                        
  labs(title = "Phenocam locations on Qikiqtaruk, Herschel Island")+
  theme_minimal()

