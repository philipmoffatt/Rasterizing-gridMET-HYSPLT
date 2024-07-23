# Load the raster
raster_file <- 'processed_data/rasters/r_2021-01-25_total_rain.tif'
raster_data <- rast(raster_file)

# Convert raster to data frame
rasters_df <- as.data.frame(raster_data, xy = TRUE)
colnames(rasters_df) <- c("lon", "lat", "total_rain")  # Ensure proper column names

# Get the map background
bbox <- c(left = -123.5, bottom = 42, right = -122.6, top = 42.6)
map_background <- get_map(bbox, zoom = 10, maptype = "stamen_terrain_background")

# Create the raster plot
main_plot <- ggmap(map_background, alpha = 0.7) +
  geom_raster(data = rasters_df, aes(x = lon, y = lat, fill = total_rain), alpha = 0.5) +
  borders("state", size = 1) +
  geom_text(data = locs, 
            aes(x = longitude, y = latitude, label = station)) +
  coord_fixed(xlim = c(-123.5, -122.6),
              ylim = c(42, 42.6)) +
  labs(
    title = 'Trajectory Local Rainfall',
    x = "Longitude",
    y = "Latitude",
    fill = 'Rainfall (mm)'
  ) +
  theme(
    plot.background = element_rect(fill = "white"),  # Set plot background color
    panel.background = element_rect(fill = "white"),  # Set panel background color
    strip.background = element_rect(fill = "white"),  # Set strip (facet label) background color
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    strip.text = element_text(size = 12),  # Adjust strip text size
    panel.spacing = unit(1, "lines")  # Increase spacing between facets
  )

# Create the inset plot of the USA outline
usa <- map_data("usa")
inset_plot <- ggplot() +
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), fill = "lightgrey", color = "black") +
  geom_point(data = locs %>% filter(station==5), aes(longitude, latitude), col = 'red')+
  theme_void()

# Combine the raster plot with the inset plot
final_plot <- ggdraw() +
  draw_plot(main_plot) +
  draw_plot(inset_plot, x = 0.1, y = 0.725, width = 0.2, height = 0.2)  # Adjust position and size as needed

# Display the raster plot
print(final_plot)