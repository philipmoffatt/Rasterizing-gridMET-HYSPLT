---

---

# Rasterizing-HYSPLIT

The Hybrid Single-Particle Lagrangian Integrated Trajectory (HYSPLIT)
model (Stein et al., 2015) is used for relating air mass trajectories
with the stable isotopic compositions of precipitation (δ2H and
δ18O). This workflow manipulates climate data and creates rasterized
atmospheric indices from HYSPLIT trajectories. The script begins with
selecting a study domain around Medford, OR, and works through a
methodology for running HYSPLIT, summarizing trajectory data into
indices, and finally creating rasters. The rasters estimate the air mass
conditions impacting the selected study area.

1.  Install and load the necessary software
2.  Run HYSPLIT and obtain trajectory files
3.  Process trajectory files and create rasters from HYSPLIT indices
4.  Conclusion

### 1. Install and load software

Software requirements:\
R-Studio 2023.12.0.369:
<https://posit.co/products/open-source/rstudio/>\
R version 4.3.2 (2023-10-31 ucrt): <https://cran.rstudio.com/>

```{r packages, echo = FALSE}

# Libraries ####

# Function to check if a package is installed, and install it if necessary
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

# List of packages
packages <- c(
  "lattice", "AOI", "sf", "raster", "rasterVis", "prism", "data.table", 
  "multcompView", "cowplot", "mgcv", "geosphere", "ggcorrplot", "anytime", 
  "mblm", "asbio", "Kendall", "trend", "qqplotr", "ggpubr", "psych", 
  "GGally", "gt", "vip", "pdp", "openair", "mapproj", "openairmaps", 
  "leaflet.extras", "randomForest", "leaflet", "webshot", "cluster", 
  "rnoaa", "ratmos", "climateR", "climate", "tidyverse", "ggrepel", 
  "rnaturalearth", "patchwork", "elevatr", "ggmap", "grid", 
  "rnaturalearthhires", "ggspatial", "terra", "FAMEFMR", "splitr", 
  "maps", "plyr", "tictoc", "furrr"
)

# Install missing packages
install_if_missing(packages)

# Additional libraries that may need to be installed from GitHub
if (!require("devtools", character.only = TRUE)) {
  install.packages("devtools")
  library(devtools)
}

# Install 'splitr' from GitHub if not installed
if (!require("splitr", character.only = TRUE)) {
  devtools::install_github("rich-iannone/splitr")
  library(splitr)
}

# Install 'climateR' from GitHub if not installed
if (!require("climateR", character.only = TRUE)) {
  devtools::install_github("mikejohnson51/climateR")
  library(climateR)
}

# Install 'climate' from GitHub if not installed
if (!require("climate", character.only = TRUE)) {
  devtools::install_github("ropensci/climate")
  library(climate)
}

```

#### Load functions

```{r functions, echo = FALSE}
# Functions ####


#' @title Run HYSPLIT
#' @description Runs the HYSPLIT model
#' @param start_date the begining of the time period to model back trajectories over (character)
#' @param end_date end of the period (character)
#' @param df a dataframe consisting of one row and four columns:  latitude (num), longitude (num), station (num), date(Date yyyy-mm-dd)
#' @return trajectory dataframe
#' 
run_trajectory_model <- function(start_date, end_date, 
                                 df) {
  # Debugging statements
  print("Starting run_trajectory_model")
  print(paste("Start date:", start_date))
  print(paste("End date:", end_date))
  print("Locations DataFrame:")
  print(df)
  
  # Iterate over each location in the locations dataframe
  lat <- df$latitude
  lon <- df$longitude
  
  # Get the StationNumber associated with the current lat and lon
  station <- df$station
  
  # Create trajectory model for the current location
  trajectory_model <- create_trajectory_model() %>%
    # builds a matrix of points, removed for example purposes
    # add_grid(
    #   lat = lat,
    #   lon = lon,
    #   range = c(0.28, 0.28),
    #   division = c(0.14, 0.14))  %>%
    # load selected parameters
    add_trajectory_params(
      lat = lat,
      lon = lon,
      height = 1000, # m agl
      duration = 120, # hours
      days = seq(as.Date(start_date), as.Date(end_date), by = "day"),
      daily_hours = c(0, 6, 12, 18),
      direction = "backward",
      met_type = "reanalysis",
      met_dir = NULL,
      model_height = 10000) %>% # trajectory upper limit
    splitr::run_model()
  
  # Get trajectory output as a data frame for the current location
  trajectory_df <- trajectory_model %>% get_output_tbl()
  
  # Add a column to the dataframe with the corresponding StationNumber
  trajectory_df$station <- station
  
  # Handle trajectories that pass the prime meridian
  trajectory_df$lon[which(trajectory_df$lon > 0)] <- trajectory_df$lon[which(trajectory_df$lon > 0)] - (180*2)
  
  # Return the list of named trajectory data frames for all locations
  return(trajectory_df)
}

#' @title North America Coastline
#' @description provides the spatial data for the coastline necessary for delineating between the whole trajectory indices and those indices that reflect airmass characteristics overland. 
#' 
get_northAmerica <- function(){
  usa = ne_countries(country = "United States of America") 
  canada = ne_countries(country = "Canada") 
  mexico = ne_countries(country = "Mexico") 
  
  canada <- canada %>% st_as_sf
  usa <- usa %>% st_as_sf
  mexico <- mexico %>% st_as_sf
  
  north_america <- st_union(usa, canada) %>% st_union(., mexico) %>% 
    mutate(LC = factor("land")) %>% dplyr::select(LC)
}

#' @title Walk Along Trajectory
#' @description Estimates for the number of phase changes along a trajectory are calculated along the trajectory. New columns 'freeze_cycles' and sat_cycles' are created. 
#' @param traj_in a dataframe trajectory file already modified to include land boundaries from `get_northAmerica`
#' @return dataframe with HYSPLIT indices
#' 
walk_traj_function= function(traj_in = traj_in){
  traj_in$LC_group =1
  temp = 1
  freeze_cycles = 0
  sat_cycles = 0
  groupslc <- for(i in 2:nrow(traj_in)){
    if(traj_in$landcover[i]!=traj_in$landcover[i-1]){
      temp = temp +1
    }
    
    if(traj_in$air_temp[i]<273 & traj_in$air_temp[i-1]>273){ # count of times the air temperature rose from below freezing to above freezing
      freeze_cycles = freeze_cycles+1 
    }
    
    if(traj_in$rh[i]>95 & traj_in$rh[i-1]<95){ # count of times relative humidity changed from above 95% to below 95%
      sat_cycles = sat_cycles+1 
    }
    
    traj_in$LC_group[i] = temp  
  }
  traj_in = traj_in %>% mutate(
    freeze_cycles = freeze_cycles, 
    sat_cycles = sat_cycles)
  
  
  return(traj_in)
}

#' @title Sum Land
#' @description The functions takes a trajectory dataframe, filters to trajectory data occuring overland, and processes each column into daily indices.
#' @param test_sf an dataframe trajectory file
#' @return dataframe with HYSPLIT indices for trajectory paths overland
#' 
sum_land_fn = function(test_sf){
  test_sum_land <- test_sf %>% 
    filter(landcover =="land") %>% 
    dplyr::summarise(
      location = first(location),
      date = first(date), 
      run = first(run),
      temp_range_l = max(air_temp)-min(air_temp),
      temp_mean_l = mean(air_temp),
      elevation_range_l =  max(terr_msl)-min(terr_msl),
      height_range_l = max(height)-min(height),
      total_distance_l = sum(dist_test) %>% replace_na(0),
      total_time_l = length(dist_test) %>% replace_na(0),
      total_rain_l = sum(rainfall) %>% replace_na(0),
      lat_shore_l = last(lat),
      long_shore_l = last(lon),
      freeze_time_l = sum(air_temp<273) %>% replace_na(0),
      sat_time_l = sum(rh>95) %>% replace_na(0),
    ) 
  return(test_sum_land)
}


#' @title Sum All
#' @description The functions takes a trajectory dataframe and processes each column into daily indices.
#' @param test_sf an dataframe trajectory file
#' @return dataframe with HYSPLIT indices for whole trajectory paths
#' 
sum_all_fn = function(test_sf){
  test_sum <- test_sf %>% dplyr::summarise(
    lat_120hr = last(lat),
    long_120hr = last(lon),
    run = first(run),
    location = first(location),
    date = first(date),
    lat = first(lat_i),
    lon = first(lon_i), 
    receptor = first(receptor), 
    temp_range = max(air_temp)-min(air_temp),
    temp_mean = mean(air_temp),
    elevation_range = max(terr_msl)-min(terr_msl),
    height_range = max(height)-min(height),
    total_distance = sum(dist_test),
    total_rain = sum(rainfall),
    d1hr_rain = first(rainfall),
    freeze_time = sum(air_temp<273),
    sat_time = sum(rh>95),
    sat_cycles = max(sat_cycles), 
    freeze_cycles = max(freeze_cycles),
  )
  return(test_sum)
} 

#' @title Backtrajectory Summary
#' @description Given a trajectory file (rds), the helper functions `get_northAmerica`, `walk_traj_function`, `sum_land_fn`, and `sum_all_fn` are applied to the trajectory data to produce a datatable of HYSPLIT indices. 
#' @param H a raw, rds trajectory file stored in a folder named by date (2021-01-20).
#' @param north_America is a sf dataframe with columns 'LC' (landcover) and 'geometry'. 
#' @return datatable with HYSPLIT indices
#'
backtrj_summarize_k <- function(H = backtrj_dt, 
                                north_america = north_america){
  
  test <- read_rds(H) %>% 
    data.table() #$backtrj_files
  
  # Convert 'traj_dt_i' to Pacific time
  test$traj_dt_i_pacific <- with_tz(test$traj_dt_i, tz = "America/Los_Angeles")
  
  # Extract the number after the last underscore
  test$location <- H %>%
    str_extract("_\\d+$") %>%  # Extract the underscore followed by digits at the end of the string
    str_remove("_") %>%        # Remove the underscore
    as.numeric()               # Convert to numeric
  
  test$date = H %>% 
    dirname() %>% 
    basename() # $site_dates 
  
  # filter to the pacific time trajectory model runs
  test <- test %>% 
    filter(format(traj_dt_i_pacific, "%Y-%m-%d") == date)
  
  # get the distances of each step 
  l <- list( a = test$lon[1:nrow(test)-1],
             b = test$lat[1:nrow(test)-1],
             c = test$lon[2:nrow(test)],
             d = test$lat[2:nrow(test)])
  
  test[,dist_test := purrr::pmap(l, 
                                 function (a,b,c,d)  distVincentyEllipsoid(c(a,b),
                                                                           c(c,d))) %>%
         unlist(.) %>% c(0,.)]
  
  test[hour_along==0]$dist_test = 0  
  
  test_sf <- st_as_sf(test, 
                      coords = c("lon", "lat"), 
                      crs = st_crs(4326),remove = F)
  
  test_sf <- test_sf %>% 
    mutate(landcover = ifelse(st_intersects(test_sf ,
                                            north_america %>% 
                                              st_make_valid(),
                                            sparse = F),
                              "land", "water") %>% 
             factor)
  
  test_it_all = test_sf %>% 
    st_drop_geometry() %>% 
    split(test_sf$run)
  
  # running each trajectory run into walk function for sat/freeze cycles
  test_sf_list = purrr::map(test_it_all, 
                            walk_traj_function)
  # summarizing each trajectory with sum all function
  test_sum_all = purrr::map_df(test_sf_list, 
                               sum_all_fn)
  # summarizing each trajectory for land based calcs
  test_sum_land = purrr::map_df(test_sf_list, 
                                sum_land_fn)
  
  test_sum_all2 <- left_join(test_sum_all, test_sum_land)
  
  test_sum_all3 <- test_sum_all2 %>% 
    dplyr::select(date, location, 
                  lat, lon, lat_120hr, long_120hr, lat_shore_l, long_shore_l,
                  elevation_range, elevation_range_l,
                  height_range, height_range_l,
                  total_distance, total_distance_l,
                  total_time_l,
                  temp_mean,temp_mean_l, temp_range, temp_range_l,
                  total_rain, total_rain_l, d1hr_rain,
                  sat_time, sat_time_l, sat_cycles,
                  freeze_time, freeze_time_l, freeze_cycles,
    ) %>% # drop crap variables
    group_by(date) %>% 
    summarize_all(mean, na.rm = TRUE) %>% 
    ungroup()
  
}

#' @title Summarize Trajectory to Raster
#' @description Operating on a parent directory folder with subdirectories named by date (2021-01-20), the function uses `backtraj_summarize_k` to aggregate trajectory data files within the date folder. The datatable produced by `backtraj_summarize_k` for each trajectory file for a given date is rasterized and stored in a user defined output parent directory. The result is a spatRaster object (.tif) that contains all trajectory locations and has layers for each daily HYSPLIT indice. The output files are name r_'date'.tif. 
#' @param data_in_path parent directory where date folders containing trajectory files for each location are stored.
#' @param out_path output parent directory where daily spatRaster objects for the study domain are stored as .tif files. 
#' @return daily spatRaster objects as .tif files
#'
summarize_to_raster= function(data_in_path, out_path){
  
  data_in_file = list.files(data_in_path, full.names = T)
  # Define the function to be applied to each element in the list in parallel
  
  print(paste("the first file in the data path ",data_in_file[[1]]))
  
  # Regular expression to match a date in the format YYYY-MM-DD
  date_pattern <- "\\d{4}-\\d{2}-\\d{2}"
  
  # Extract the date
  temp_date <- str_extract(data_in_file[[1]], date_pattern)
  
  print(paste("the folder being processed is", temp_date))
  
  tic()
  traj_back = purrr::map(data_in_file,
                         function(H) backtrj_summarize_k(H = H, north_america = north_america))
  toc()
  print("traj_back has run")
  
  #### 4) make raster 
  
  traj_back_df_filled = traj_back %>% rbind.fill()
  # Determine extent and resolution
  extent_x <- range(traj_back_df_filled$lon)
  extent_y <- range(traj_back_df_filled$lat)
  res_x <- 0.1
  res_y <- 0.1
  
  # Create a blank raster with the defined extent and resolution
  r <- raster(xmn = extent_x[1], xmx = extent_x[2],
              ymn = extent_y[1], ymx = extent_y[2],
              res = c(res_x, res_y), crs = crs("epsg:4326"))
  
  # Define the columns to rasterize (adjust as necessary)
  value_column <- names(traj_back_df_filled)[5:28]  # Columns to rasterize
  
  # Rasterize each value column
  for (value_col in value_column) {
    print(paste("Rasterizing column:", value_col))
    
    # Create a raster for the current value column
    r_value <- rasterize(traj_back_df_filled[, c("lon", "lat")],
                         r, field = traj_back_df_filled[[value_col]],
                         fun = mean, na.rm = TRUE)
    
    # Save the raster
    temp_out_path <- file.path(out_path, paste0("r_", temp_date, "_", value_col, ".tif"))
    writeRaster(r_value, temp_out_path, overwrite = TRUE)
  }
  print(paste0("Trajectory rasters built for ", temp_date,"."))
}
```

### 2. Run HYSPLIT and obtain trajectory files

-   Building a matrix of receptor points

Pulling trajectories from multiple points around a target location can
better approximate the movement of an air mass. We construct a matrix of
initiation points around each target location that matches the
resolution of the reanalysis data product used to force HYSPLIT
(\~32km). Data used to force the HYSPLIT model varies extensively. Here we use the most accessible and course dataset to reduce the computational demands associated with downloading it and conducting raster calcuations. The
interactive example below demonstrates a matrix of nine points around
the Medford, OR airport (MFR). Simply respond to the prompt for a directory path for downloading trajectory files. 

```{r visualize_receptors}
# study area locations. Medford Airport is point 5.
locs <- data.frame(
  station = c('1', '2', '3', '4', '5', '6', '7', '8', '9'),
  latitude = c(42.235032, 42.235032, 42.23503, 42.375032, 42.375032, 42.375032, 42.515032, 42.515032, 42.515032),
  longitude = c(-123.017016, -122.877016, -122.7370160, -123.017016, -122.877016, -122.737016, -123.017016, -122.877016, -122.737016)
)

# view the target locations for HYSPLIT models
leaflet(data = locs) %>%
  addTiles() %>%
  addMarkers(
    label = ~station,
    clusterOptions = markerClusterOptions()
  )
```

-   Running HYSPLIT

The example code below generates air mass back trajectories for nine
station locations. Using the matrix method, there are nine air parcel
back trajectories run for each station, four times daily. The processing
time for the 324 back trajectories will vary from 10-20 min. Trajectory
files are stored in "processed_data/example_traj" and can be loaded
instead of downloaded.

```{r run_hysplit, echo = FALSE}

# Create the data frame
df <- tibble::tibble(
  station = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9),
  latitude = c(42.235032, 42.235032, 42.23503, 42.375032, 42.375032, 42.375032, 42.515032, 42.515032, 42.515032,
               42.235032, 42.235032, 42.23503, 42.375032, 42.375032, 42.375032, 42.515032, 42.515032, 42.515032),
  longitude = c(-123.017016, -122.877016, -122.737016, -123.017016, -122.877016, -122.737016,
                -123.017016, -122.877016, -122.737016,
                -123.017016, -122.877016, -122.737016,
                -123.017016, -122.877016, -122.737016,
                -123.017016, -122.877016, -122.737016),
  date = as.Date(c('2021-01-25', '2021-01-25', '2021-01-25', '2021-01-25', '2021-01-25', '2021-01-25',
                   '2021-01-25', '2021-01-25', '2021-01-25',
                   '2021-01-26', '2021-01-26', '2021-01-26', '2021-01-26', '2021-01-26', '2021-01-26',
                   '2021-01-26', '2021-01-26', '2021-01-26'))
)

# Print the data frame
head(df)

# set directory path for where you want the traj files saved
out_dir <- readline(prompt = "Please enter the directory path where you want to store the data, no quotes: ")

# Print the user input to confirm
cat("You entered:", out_dir, "\n")

# simple for-loop that takes in each row of the provided in df
out = for(i in 1:nrow(df)){
  H = df[i,]
  start_time = Sys.time() # log when the model process initiates
  start_date = H$date # define the day being processed
  end_date = H$date + 1 # define the end of the daily trajectory
  out_dir_traj = file.path(out_dir, start_date) # directory to store the trajectory file
  print(out_dir_traj)
  
  if(!dir.exists(out_dir_traj)){dir.create(out_dir_traj, recursive = T)} # create the output location if needed
  
  # run the trajectory model on the date/location combo
  traj_temp <- run_trajectory_model(df =  H,
                                    start_date = start_date,
                                    end_date = end_date)
  end_time = Sys.time()
  
  print(paste("total time per air parcel =", end_time - start_time, "seconds!"))
  
  write_rds(traj_temp, file = file.path(out_dir_traj,paste0("traj_",H$station))) # write the trajectory file with station name into daily folder
}
```

The function `run_trajectory_model` is a wrapper for `trajectory_model`
that provides an application to a study domain. If using parallel
processing or executing on HPC, use the `exec_dir` argument in
`trajectory_model` to create temporary locations for HYSPLIT trajectory
data for each initiation point. It is important to provide the entire
path for this argument so temporary traj files do not overwrite one
another.

\-`start_date`: a string date to begin HYPSLIT models from
('2021-01-25')

\-`end_date`: a string date to end HYPSLIT models on ('2021-01-26')

\-`df`: a dataframe with columns 'location,' 'latitude,' and
'longitude.'

\-`met_dir_in`: optional directory where .gbl files that HYSPLIT runs on
can be downloaded to reduce computation time. Running the function with
this argument defined as 'NULL' results in the .gbl files being
downloaded to the working directory from URL
'<ftp://arlftp.arlhq.noaa.gov/archives/reanalysis>'

### 3. Process trajectory files and create rasters from HYSPLIT indices

-   Daily HYSPLIT indices

We can summarize the meteorological data from HYSPLIT into indices that
estimate air mass conditions. For detailed calculations for each
variable, see the "functions" section above. We run trajectories
for the entire study domain one day at a time in this example. HYSPLIT
trajectories run for a group of nine stations around Medford, OR, over
one day, 2021-01-25, which takes only 10-15 seconds to process.

One function processes the HYSPLIT trajectory files into rasters,
`summarize_to_raster`. It incorporates multiple helper functions that
aggregate the trajectory meteorological data over the entire trajectory
and the portion of the trajectory overland into indices. The result is
one daily value for each HYSPLIT variable. Referring to Section 2, it is
crucial to understand that each daily value is a mean, or otherwise
aggregated value, from nine trajectories surrounding a single station
site. The nine air parcel trajectories are run four times a day, so the
resulting daily air mass variables incorporate 36 trajectory files.

```{r create_rasters, echo = FALSE}

#1: the directory that holds each day of trajectories, choosing the first day
data_in_path <- list.files(out_dir, full.names = T)[1]

#2: the out directory provides a location to produce daily folders with HYSPLIT rasters for the study domain.
# set directory path for where you want the traj files saved
out_path <- readline(prompt = "Please enter the directory path where you want to store the rasters: ")

# Print the user input to confirm
cat("You entered:", out_path, "\n")

# Create the output directory if it doesn't exist
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

print(paste("the data_in_path is ",data_in_path))

print(paste("the out_path is ",out_path))

### get north_america
north_america = get_northAmerica()

### run the day of rasters.
summarize_to_raster(data_in_path = data_in_path,
                    out_path = out_path)

print(paste("rasters saved to", out_path))

```

-   Visualize HYSPLIT Rasters

The rasters of total rainfall are plotted, showing the rainfall history
for air masses associated with each station. The rasters a plotted over
a Stamen background map that shows terrain.

```{r raster_plot, echo = FALSE}

# Load the raster - total rain along the trajectories
raster_file <- paste0(out_path,'/r_2021-01-25_total_rain.tif')
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
  annotate("text", x = -123.25, y = 42.02, 
           label = '© Stadia Maps © Stamen Design © OpenMapTiles © OpenStreetMap contributors',
           hjust = 0, vjust = 0, size = 3, color = "black")+
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
  geom_point(data = locs %>% filter(station==5), aes(longitude, latitude), fill = 'red', col = "white", shape = 21)+
  theme_void()

# Combine the raster plot with the inset plot
final_plot <- ggdraw() +
  draw_plot(main_plot) +
  draw_plot(inset_plot, x = 0.1, y = 0.725, width = 0.2, height = 0.2)  # Adjust position and size as needed

# Display the raster plot
print(final_plot)
```
### 4. Download an elevation raster
An elevation raster is used in this example to obtain a coordinate reference system and resolution. 

```{r elevation, echo = FALSE}
# set boundary for study area to states
sa_bounds <- aoi_get(state = c("OR"))

# pull elevation data
elevation_sa <- elevatr::get_elev_raster(locations = sa_bounds, z = 8, clip = "locations")
names(elevation_sa) <- "elevation_m"

```

### 5. Create a Random Forest model using an example dataset

Only HYSPLIT indices and the elevation predictors are included in this data frame to simplify the workflow in this example. Incorporating other data sources such as gridMET requires following the steps shown here for elevation data. The extent and reference system for the rasters must match. The layer names in the rasters must also match the model predictors. 

```{r rf_model, echo = FALSE}
# using a simplified set of predictors that refer to rasters already produced in this workflow
# example data from roughly 1 month
data <- data.frame(
  lat_120hr = c(56.38302612, 46.16249847, 46.93535995, 61.05019379, 38.10211182, 43.60605621, 51.19655609, 63.22927856, 43.06427765, 54.13580704, 53.26830673, 35.14152908, 40.42441559, 39.84994507, 27.26408386, 34.76025009, 45.0886116, 43.67272186, 43.42355728, 47.14197159, 39.81894302, 51.13405609, 50.68613815, 48.7120018, 49.53575134, 43.20563889, 55.19233322, 51.97366714, 51.56322098),
  long_120hr = c(-149.4461975, -4.685416698, -8.321332932, -142.681839, -159.2973022, -145.81633, -125.2227783, -151.8461456, -163.4849091, -155.2358398, -142.8508301, -144.8025818, -331.6868896, -354.1795554, -155.3153687, -152.7663269, -358.9870278, -348.1305838, -355.9036665, -87.41124725, -95.59505463, -270.587944, -354.9617224, -204.4577179, -276.2919464, -148.4311371, -116.2144699, -303.3214989, -291.6477509),
  elevation_range = c(949.7388916, 536.916687, 579.5888672, 966.3388672, 561.7444458, 522.986084, 554.9611206, 701.4555664, 561.1416626, 668.888916, 575.3166504, 529.1444702, 549.263916, 572.7694702, 546.5861206, 528.7000122, 538.8222046, 573.3778076, 559.6361084, 543.8444214, 541.805542, 524.2666626, 554.6972046, 528.0388794, 612.847229, 495.4222107, 548.3499756, 684.9555664, 788.2471924),
  height_range = c(2232.491699, 983.2999878, 1204.502808, 1842.719482, 781.5861206, 1003.180542, 911.2166748, 1067.06665, 2639.341553, 943.347229, 1907.983276, 950.8721924, 1595.844482, 2633.241699, 997.9277954, 965.6472168, 1261.630615, 1302.047241, 709.7555542, 1062.930542, 2234.875, 1721.891724, 2207.097168, 1264.244385, 1163.636108, 997.4833374, 1005.280579, 914.597229, 1667.077759),
  total_distance = c(2935601, 5862695, 5915071.5, 4018617.75, 4646686.5, 4431209, 4718712.5, 3593392, 5261700, 3575372.75, 2483001, 2785625.25, 6380266, 6085517, 4965872.5, 4783459, 5247175.5, 5013654, 5082132, 4891615.5, 5470483.5, 6674136.5, 5706233, 7334039.5, 5436367.5, 4327792.5, 5009545, 7372173, 7296082.5),
  temp_mean = c(269.872406, 275.7507935, 265.1822205, 270.6975708, 282.8326416, 262.4943848, 275.927124, 270.5460205, 272.8890381, 271.0884399, 259.5909424, 263.987915, 265.8288269, 275.5342407, 287.1350403, 278.198822, 269.4921875, 274.504364, 277.5024414, 273.2569885, 278.1085815, 268.935791, 266.9195251, 273.5153503, 270.0193787, 282.6733398, 271.7858276, 265.0660095, 266.6869507),
  temp_range = c(21.98888969, 14.63333321, 15.25833321, 27.89722252, 11.13333321, 14.45555592, 10.76944447, 17.73333359, 20.5027771, 13.85000038, 17.19166756, 13.53611088, 23.23888969, 31.69166756, 21.21111107, 18.11388969, 11.12222195, 16.5027771, 13.88333321, 13.10000038, 20.28888893, 17.875, 23.19722176, 17.96666718, 15.6916666, 14.33888912, 14.06944466, 20.42499924, 27.89444351),
  total_rain = c(15.572222233, 12.40833378, 21.57777786, 7.066666603, 4.150000095, 16.90833378, 10.93055534, 7.555555344, 8.433333397, 7.566666603, 25.544444561, 19.16944504, 27.06944466, 15.67222214, 16.75555611, 5.81111145, 6.163888931, 6.280555725, 10.16666698, 7.191666603, 6.538888931, 4.541666508, 4.544444561, 6.641666889, 6.402777672, 8.722222328, 8.47222233, 10.51388931, 9.927778244),
  d1hr_rain = c(0.050000001, 0.025, 0, 0.125, 0.200000003, 0.275000006, 0.125, 0.150000006, 0.100000001, 0.025, 0.174999997, 0.550000012, 0.349999994, 0.275000006, 0.550000012, 0.200000003, 0.025, 0.174999997, 0.275000006, 0.025, 0.400000006, 0.050000001, 0.224999994, 0.25, 0.075000003, 0.150000006, 0.100000001, 0.125, 0.075000003),  sat_time = c(1.555555582, 5.5, 7.194444656, 8.527777672, 2.833333254, 0.25, 16.25, 1.055555582, 8.777777672, 26.41666603, 1.75, 22.05555534, 4, 7.805555344, 17.97222137, 39.5, 15.61111069, 7.638888836, 11.91666698, 13.88888931, 9.138889313, 10.86111069, 3.75, 5.638888836, 10.36111069, 0.555555582, 10.58333302, 0.944444418, 2.305555582),
  sat_cycles = c(0.333333343, 0.361111104, 1.138888836, 1.5, 0.888888896, 0.194444448, 1.333333373, 0.583333313, 1.5, 3.111111164, 0, 1.777777791, 1.027777791, 0.5, 1.722222209, 3.027777672, 1.416666627, 0.805555582, 0.611111104, 1.694444418, 0.472222209, 1.222222209, 0.5, 2, 1.25, 0.305555552, 1.444444418, 0.416666657, 0.472222209),
  freeze_time = c(61.16666794, 30.58333397, 31.55555534, 63.61111069, 0, 1.583333373, 27.58333397, 61.58333206, 48.02777863, 68.69444275, 29.44444466, 0.083333336, 35.77777863, 39.66666794, 0.583333313, 27.05555534, 91.66666412, 56.97222137, 19.41666603, 53.22222137, 16.19444466, 84.3611145, 72.08333588, 40.91666794, 73.5, 1.972222209, 68.22222137, 97.1388855, 80.25),
  freeze_cycles = c(0.666666687, 0.583333313, 0.361111104, 0.916666687, 0, 0.222222224, 0.666666687, 0.694444418, 0.972222209, 1, 0.638888896, 0.027777778, 0.805555582, 0.583333313, 0, 0.444444448, 0.361111104, 0.666666687, 0.5, 0.444444448, 0.472222209, 0.777777791, 0.944444418, 1.055555582, 0.972222209, 0.277777791, 0.666666687, 0.944444418, 0.833333313),
  elevation_m = c(130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130),
  d2h_daily = c(-98.4, -66.6, -109.8, -75, -68, -94.8, -53.3, -68.3, -30.1, -34.8, -103.5, -126.3, -114.1, -55.9, -43.6, -43, -43.2, -85.6, -57.3, -62.7, -32.7, -28.4, -21.5, -24.3, -16.4, -50.1, -35.7, -48, -86.7)
)

# model for predicting d2h using only elevation and hysplit indices
rf_model <- randomForest(d2h_daily ~ ., 
                                         data = data, 
                                         ntree = 5000,
                                         keep.forest=TRUE, 
                                         importance=TRUE
)

```
### 6. Prepare the predictor rasters
The raster manipulation shown here produces one daily raster stack for 2021-01-25 that includes layers with the predictors spatially represented. 

```{r raster_prep, echo = FALSE}
# rasters for predictors
# List all the raster files for a specific date
raster_files <- list.files(out_path, pattern = "2021-01-25", full.names = TRUE)

# Create a SpatRaster stack
daily_stack <- rast(raster_files)

# Name the layers in the stack
names(daily_stack) <- gsub(".*r_2021-01-25_|\\.tif", "", raster_files)

# Extract CRS from elevation_sa and convert to character
crs_string <- as.character(crs(elevation_sa))

# Project daily_stack to the CRS of elevation_sa
hysplit_rast <- project(daily_stack, crs_string)

# Convert elevation_sa to SpatRaster
if (inherits(elevation_sa, "RasterLayer")) {
  elevation_sa <- rast(elevation_sa)
}

# crop elevation to match the example extent
elevation_crop <- crop(elevation_sa, hysplit_rast)

# Resample hysplit_rast to match elevation_sa
hysplit_resample <- resample(hysplit_rast, elevation_crop, method = "bilinear")

```
### 7. Make an isoscape

```{r_isoscape, echo = FALSE}
# Combine all predictors into a single SpatRaster
predictors <- c(hysplit_resample, elevation_crop)

# Make predictions using the terra predict function
prediction_raster <- predict(object = predictors, model = rf_model)

# plot the d2h isoscape for 2021-01-25
rasters_df <- as.data.frame(prediction_raster, xy=TRUE)

ggplot(rasters_df) +
  geom_raster(aes(x = x, y = y, fill = lyr1)) +
  geom_point(data = locs, aes(x = longitude, y = latitude),
             shape = 21, color = "white", fill = "red") +
  annotation_north_arrow(location = "tl",
                         style = north_arrow_orienteering()) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  coord_sf(crs = 4326)+
  theme_bw()+
  labs(
    subtitle = expression(delta^2 * "H ‰ VSMOW"),
    fill = "‰",
    y = "",
    x ="")

```
### 4. Conclusion

The presented methodology incorporates selecting a study area and
developing rasterized estimates of daily air mass conditions relevant to
local conditions and the analysis of isotopic variability in daily
precipitation. Rasters enable the projection of predictive relationships
between air mass conditions and isotopic composition into unmeasured
locales. Comparing isotope maps, or isoscapes, with field measurements from
similar locations or broader regions can also highlight how well
isotopic variability is explained by air mass conditions.
