# study area locations. Medford Airport is point 5.
locs <- data.frame(
  station = c('1', '2', '3', '4', '5', '6', '7', '8', '9'),
  latitude = c(42.235032, 42.235032, 42.23503, 42.375032, 42.375032, 42.375032, 42.515032, 42.515032, 42.515032),
  longitude = c(-123.017016, -122.877016, -122.7370160, -123.017016, -122.877016, -122.737016, -123.017016, -122.877016, -122.737016)
)

# examine the locations chosen on a map
leaflet(data = locs) %>%
  addTiles() %>%
  addMarkers(
    label = ~station,
    clusterOptions = markerClusterOptions()
  )

# Convert locs to a simple features data frame with a CRS
locs_sf <- st_as_sf(locs, coords = c("longitude", "latitude"), crs = 4326)

# plot the locations 
ggplot() +
  geom_sf(data = locs_sf, aes(color = station)) +
  theme_minimal() +
  labs(title = "Station Locations",
       x = "Longitude",
       y = "Latitude") +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 12))


start_date <- '2021-01-25' # beginning date for trajectory data
end_date <- '2021-01-26'   # end date for the trajectory data

# Create a sequence of dates
dates <- seq.Date(as.Date(start_date), as.Date(end_date), by = "day") 

# Create a new dataframe with repeated rows for each date
locations_df <- locs %>%
  slice(rep(1:n(), each = length(dates))) %>%
  mutate(date = rep(dates, times = nrow(locs))) %>% 
  arrange(date, station)

# write_csv(locations_df, 'raw_data/locations.csv')

