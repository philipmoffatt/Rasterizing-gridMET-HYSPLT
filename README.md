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

```{r packages}

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
  "tidyverse", "splitr", "climateR", "climate", 
  "openair", "maps", "terra", "plyr", "sf", 
  "rnaturalearth", "geosphere", "data.table", 
  "raster", "leaflet", "tictoc", "furrr"
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

```{r functions}
source('R/functions.R')
```

### 2. Run HYSPLIT and obtain trajectory files

-   Building a matrix of receptor points

Pulling trajectories from multiple points around a target location can
better approximate the movement of an air mass. We construct a matrix of
initiation points around each target location that matches the
resolution of the reanalysis data product used to force HYSPLIT
(\~32km). A finer resolution likely would not produce values appreciably
different, whereas a larger resolution sacrifices spatial detail. The
interactive example below demonstrates a matrix of nine points around
the Medford, OR airport (MFR).

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

```{r run_hysplit, eval=FALSE}


# Formatted for use on hpc

args = commandArgs(trailingOnly = T)
#arg are expected in order
#1: the group run file to be iterated over - each line is one location and one day 
#2: the met dir "Hysplt/met_dir" where the .gbl files used for modeling backtrj are stored
#3: the output directory

# example 
args = c("raw_data/locations.csv", "raw_data/met_dir", "processed_data/example_traj")

# assign args
locations_in = read_csv(args[1])
met_dir_in = args[2] 
Out_dir = args[3] 

# check args were correctly loaded
if (!length(args)==3) {
  stop('arguments are expected in order
      1: the group run file to be iterated over - each line is one location and one day
      2: the met dir defalts to "Hysplt/met_dir"
      3: the base output directory Hysplt/processed_data/runs_$sys.date/'
       )
}
print(paste('Locations argument is:', args[1]))
print(paste('Meteorology argument is:', args[2]))
print(paste('Processed data directory argument is:', args[3]))

# simple for-loop that takes in each row of the provided locations_in df
out = for(i in 1:nrow(locations_in)){
    H = locations_in[i,]
    start_time = Sys.time() # log when the model process initiates
    start_date = H$date # define the day being processed
    end_date = H$date + 1 # define the end of the daily trajectory
    Out_dir_traj = file.path(Out_dir, start_date) # directory to store the trajectory file
    print(Out_dir_traj)
    
    if(!dir.exists(Out_dir_traj)){dir.create(Out_dir_traj, recursive = T)} # create the output location if needed
    
    # run the trajectory model on the date/location combo
    traj_temp <- run_trajectory_model(locations_df =  H,
                                      start_date = start_date,
                                      end_date = end_date, 
                                      met_dir_in = met_dir_in)
    end_time = Sys.time()
    
    print(paste("total time per air parcel =", end_time - start_time, "seconds!"))
    
    write_rds(traj_temp, file = file.path(Out_dir_traj,paste0("traj_",H$station))) # write the trajectory file with station name into daily folder
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

\-`location_df`: a dataframe with columns 'location,' 'latitude,' and
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
variable, see the `R/functions.R`. The script example below is formatted
to run on an HPC. The job below can be run on many processes to break up
large study domains and long time spans of interest. We run trajectories
for the entire study domain one day at a time in this example. HYSPLIT
trajectories run for a group of nine stations around Medford, OR, over
one day, 2021-01-25, which takes only 20-40 seconds to process.

One function processes the HYSPLIT trajectory files into rasters,
`summarize_to_raster`. It incorporates multiple helper functions that
aggregate the trajectory meteorological data over the entire trajectory
and the portion of the trajectory overland into indices. The result is
one daily value for each HYSPLIT variable. Referring to Section 2, it is
crucial to understand that each daily value is a mean, or otherwise
aggregated value, from nine trajectories surrounding a single station
site. The nine air parcel trajectories are run four times a day, so the
resulting daily air mass variables incorporate 36 trajectory files.

```{r create_rasters}
# obtain necessary functions
source("R/functions.R")

args = commandArgs(trailingOnly = T)
#arg are expected in order
#1: the directory that holds each day of trajectories -- this will be run from a csv.
# Each job accesses data in one folder that contains the trajectories for a given day for all points in a study domain. 

#2: the out directory provides a location to produce daily folders with HYSPLIT rasters for the study domain.

args = c("processed_data/example_traj/2021-01-25", "processed_data/rasters")

data_in_path = args[1]
print(paste("the data_in_path is ",data_in_path))

out_path = args[2]
print(paste("the out_path is ",out_path))

if (!length(args)==2) {
  stop("arguments are expected in order
      1: the day of trajectories.
      2: the out path for the rasters.")
      }

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

```{r raster_plot}

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
  geom_point(data = locs %>% filter(station==5), aes(longitude, latitude), col = 'red')+
  theme_void()

# Combine the raster plot with the inset plot
final_plot <- ggdraw() +
  draw_plot(main_plot) +
  draw_plot(inset_plot, x = 0.1, y = 0.725, width = 0.2, height = 0.2)  # Adjust position and size as needed

# Display the raster plot
print(final_plot)
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
