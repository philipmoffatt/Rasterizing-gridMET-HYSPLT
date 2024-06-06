# Rasterizing-gridMET-HYSPLT
This workflow is used for manipulating climate data and creating atmospheric indices from HYSPLIT trajectories.

Software requirements:  
R-Studio 2023.12.0.369:  https://posit.co/products/open-source/rstudio/  
R version 4.3.2 (2023-10-31 ucrt):  https://cran.rstudio.com/

### Installation
The devtools package is required to install `splitr` from github.
```{r install}
install.packages("devtools")
library(devtools)
devtools::install_github("rich-iannone/splitr")
library(splitr)
```

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
`trajectory_model` is a list containing the trajectory data for the given dates in December 2023. Refer here for more details on the package and documentation:  https://github.com/rich-iannone/splitr  
-  `lat` the latitude for the HYSPLIT model initiation location
-  `lon` the longitude
-  `height` the above ground level for air parcel initiation in meters
-  `duration` the number of hours to run the HYSPLIT model backward or forward
-  `days` the sequence of dates to model over
-  `daily_hours` the hours of the each day to initiate a trajectory model
-  `direction` tells HYSPLIT to run forward or backward from the beginning datetime
-  `met_type` user selection for atmospheric data downloaded from the NOAA FTP server to force HYSPLIT (gdas1, narr, reanalysis)
-  `model_height` sets the maximum height HYSPLIT will track air parcels in meters

The `trajectory_model` converts directly to a dataframe using splitr::get_output_table() and plotted via ggplot to examine trajectory data spatially. 

```{r}
library(sf)
library(maps) 
ggplot(data = basic_traj_df) +
    borders("world", xlim = c(-130, -60), ylim = c(20, 50), colour = "gray85", fill = "gray80") +
    geom_point(aes(x = lon, y = lat, color = traj_dt))+
    labs(title = "HYSPLIT Hourly Trajectories", subtitle = '120 hours backwards', x = "Longitude",
        y = "Latitude", color = "Datetime") +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 12))+
    theme_minimal()
```


```{r}
start_date <- '2021-01-25' # beginning date for trajectory data
end_date <- '2021-01-26'   # end date for the trajectory data

# study area locations. Medford Airport is point 4.
locs <- data.frame(
  station = c('1', '2', '3', '4', '5', '6', '7', '8', '9'),
  latitude = c(42.235032, 42.235032, 42.23503, 42.375032, 42.375032, 42.375032, 42.515032, 42.515032, 42.515032),
  longitude = c(-123.017016, -122.877016, -122.7370160, -123.017016, -122.877016, -122.737016, -123.017016, -122.877016, -122.737016)
)

example_traj_df <- run_trajectory_model(start_date = start_date,
                                        end_date = end_date,
                                        locations_df = locs,
                                        met_dir_in = NULL)
```
The function `run_trajectory_model` is a wrapper for `trajectory_model` that provides an application to a study domain. If using parallel processing or executing on HPC, use the `exec_dir` argument in `trajectory_model' to create temporary locations for HYSPLIT trajectory data for each initiation point. It is important to provide the entire path for this argument.     
- `start_date`:  a string date to begin HYPSLIT models from ('2021-01-25')  
- `end_date`:  a string date to end HYPSLIT models on ('2021-01-26')  
- `location_df`:  a dataframe with columns 'location,' 'latitude,' and 'longitude.'  
- `met_dir_in`:  optional directory where .gbl files that HYSPLIT runs on can be downloaded to reduce computation time. Running the function with this argument defined as `NULL` results in the .gbl files being downloaded to the working directory from URL 'ftp://arlftp.arlhq.noaa.gov/archives/reanalysis

 

Process trajectories 
  Daily rasters
  Cluster analysis 

Pull climate data


