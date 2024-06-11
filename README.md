---
editor_options: 
  markdown: 
    wrap: 72
---

# Rasterizing-HYSPLT

This workflow is used for manipulating climate data and creating
atmospheric indices from HYSPLIT trajectories.

Software requirements:\
R-Studio 2023.12.0.369:
<https://posit.co/products/open-source/rstudio/>\
R version 4.3.2 (2023-10-31 ucrt): <https://cran.rstudio.com/>

### Installation

The devtools package is required to install `splitr` from github.

```{r install}
#install.packages("devtools")  
library(devtools)  
#devtools::install_github("rich-iannone/splitr")  
library(splitr)
```

### Pull HYSPLIT Trajectory Data

#### Basic Use

Pull trajectory files with splitR (0.4.0.9000).

```{r, eval = FALSE}
trajectory_model <- create_trajectory_model() %>%
  add_trajectory_params(
    lat = 42,
    lon = -121,
    height = 1000,
    duration = 120,
    days = seq(as.Date('2023-12-15'), as.Date('2023-12-16'), by = "day"),
    daily_hours = "12",
    direction = "backward",
    met_type = "reanalysis",
    model_height = 10000) %>%
  splitr::run_model()
```

`trajectory_model` is a list containing the trajectory data for the
given dates in December 2023. Refer here for more details on the package
and documentation: <https://github.com/rich-iannone/splitr>\
- `lat` the latitude for the HYSPLIT model initiation location

\- `lon` the longitude

\- `height` the above ground level for air parcel initiation in meters

\- `duration` the number of hours to run the HYSPLIT model backward or
forward

\- `days` the sequence of dates to model over

\- `daily_hours` the hours of the each day to initiate a trajectory
model

\- `direction` tells HYSPLIT to run forward or backward from the
beginning datetime

\- `met_type` user selection for atmospheric data downloaded from the
NOAA FTP server to force HYSPLIT (gdas1, narr, reanalysis)

\- `model_height` sets the maximum height HYSPLIT will track air parcels
in meters

The `trajectory_model` converts directly to a dataframe using
splitr::get_output_table(). Plot the trajectories via ggplot to examine
data spatially.

```{r}
library(sf)
library(maps) 
ggplot(data = basic_traj_df) +
    borders("world", xlim = c(-130, -60), ylim = c(20, 50), colour = "gray85", fill = "gray80") + # background map with the extent over North America
    geom_point(aes(x = lon, y = lat, color = traj_dt))+ # air parcel location data, colored by the datetime
    labs(title = "HYSPLIT Hourly Trajectories", subtitle = '120 hours backwards', x = "Longitude",
        y = "Latitude", color = "Datetime") +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 12))+
    theme_minimal()
```

The plot shows two trajectories initiated at 12 PM on December 15 and 16
2023. The trajectories are colored by time, lighter is more recent and
darker as the hour along trajectory approaches a maximum of 120 hours
prior to the imitation time.

#### Matrix Method

Pulling trajectories from multiple points around a target location can
better approximate the movement of an airmass. We construct a matrix of
initiations points around each target location that matches the
resolution of the reanalysis data product used to force HYSPLIT
(\~32km). Finer resolution likely would not produce values appreciably
different, whereas a larger resolution sacrifices spatial detail.

```{r}
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

```{r}
args = commandArgs(trailingOnly = T)
#arg are expected in order
#1: the group run file to be iterated over - each line is one location and one day 
#2: the met dir "Hysplt/met_dir" where the .gbl files used for modeling backtrj are stored
#3: the output directory
#4: temp folder that provides space for parallel jobs to run hysplit models simultaneously

# example 
args = c("raw_data/locations.csv", "raw_data/met_dir", "processed_data/example_traj", paste0("C:/Users/philip.moffatt/Documents/GitHub/Rasterizing-gridMET-HYSPLT/processed_data/exec_dir/runs_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S")))

locations_in = read_csv(args[1]) # "Hysplt/working_dir/traj_group_$SLURM_ARRAY_TASK_ID" array controls # jobs to run
met_dir_in = args[2] # "Hysplt/met_dir"
Out_dir = args[3] # "Hysplt/processed_data/test_runs_2023-09-12"
exec_dir = args[4] # "parent_folder/$SLURM_ARRAY_TASK_ID"

if (!length(args)==4) {
  stop('argumentes are expected in order
      1: the group run file to be iterated over - each line is one location and one day
      2: the met dir defalts to "Hysplt/met_dir"
      3: the base output directory Hysplt/processed_data/runs_$sys.date/
      4: exec_dir where the in/out model files are loaded')
}
print(paste('Locations argument is:', args[1]))
print(paste('Meteorology argument is:', args[2]))
print(paste('Processed data directory argument is:', args[3]))
print(paste('Execution directory argument is:', args[4]))

# Split to list
test_list =(locations_in %>% 
              split(locations_in$station))

out = lapply(
                X = test_list, FUN =  function(H){
                         start_time = Sys.time()
                         start_date = H$Date
                         end_date = H$Date + 1
                         Out_dir = file.path(Out_dir, start_date)
                         
                         if(!dir.exists(Out_dir)){dir.create(Out_dir, recursive = T)}
                         
                         
                         traj_temp <- run_trajectory_model(locations_df =  H,
                                                           start_date = start_date,
                                                           end_date = end_date, 
                                                           met_dir_in = met_dir_in,
                                                           exec_dir = exec_dir)
                         end_time = Sys.time()
                         
                         print(paste("total time per station =", end_time - start_time, "seconds!"))
                         
                         write_rds(traj_temp, file = file.path(Out_dir,paste0("traj_",H$station)))
                       })
```

The function `run_trajectory_model` is a wrapper for `trajectory_model`
that provides an application to a study domain. If using parallel
processing or executing on HPC, use the `exec_dir` argument in
`trajectory_model` to create temporary locations for HYSPLIT trajectory
data for each initiation point. It is important to provide the entire
path for this argument.

\-`start_date`: a string date to begin HYPSLIT models from
('2021-01-25')

\-`end_date`: a string date to end HYPSLIT models on ('2021-01-26')

\-`location_df`: a dataframe with columns 'location,' 'latitude,' and
'longitude.'

\-`met_dir_in`: optional directory where .gbl files that HYSPLIT runs on
can be downloaded to reduce computation time. Running the function with
this argument defined as 'NULL' results in the .gbl files being
downloaded to the working directory from URL
'<ftp://arlftp.arlhq.noaa.gov/archives/reanalysis>

### Process trajectories

#### Daily rasters

```{r}
source("R/traj_agg_2_raster.R")

args = commandArgs(trailingOnly = T)
#arg are expected in order
#1: the directory that holds each day of trajectories -- this will be run from a csv.
# Each job accesses data in one folder that contains the trajectories for a given day for all points in a study domain. 

#2: the out directory provides a location to produce daily folders with HYSPLIT rasters for the study domain.

args = c("/scratch/user/philip.moffatt/20231201_143053/uncompressed_files/2022-03-04", "/scratch/user/philip.moffatt/20231201_143053/rasters")

data_in_path = args[1]
print(paste("the data_in_path is ",data_in_path))

out_path = args[2]
print(paste("the out_path is ",out_path))

if (!length(args)==2) {
  stop("argumentes are expected in order
      1: the day of trajectorys.
      2: the out path for the rasters.")
      }

### get north_america
north_america = get_northAmerica()

### run the day of rasters.
summarize_to_raster(data_in_path = data_in_path,
                    out_path = out_path)
```

#### Cluster analysis

Pull climate data
