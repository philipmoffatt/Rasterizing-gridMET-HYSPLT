---
editor_options: 
  markdown: 
    wrap: 72
---

# Rasterizing-HYSPLIT

This workflow is used for manipulating climate data and creating
atmospheric indices from HYSPLIT trajectories.

Software requirements:\
R-Studio 2023.12.0.369:
<https://posit.co/products/open-source/rstudio/>\
R version 4.3.2 (2023-10-31 ucrt): <https://cran.rstudio.com/>

#### Installing and loading packages

```{r}

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

```{r}
source('R/functions.R')
```

#### Matrix method

Pulling trajectories from multiple points around a target location can
better approximate the movement of an airmass. We construct a matrix of
initiations points around each target location that matches the
resolution of the reanalysis data product used to force HYSPLIT
(\~32km). Finer resolution likely would not produce values appreciably
different, whereas a larger resolution sacrifices spatial detail. The
interactive example below demonstrates a matrix of nine points around
the Medford, OR airport (MFR).

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

#### Modeling trajectories with HYSPLIT

The example code below generates air mass back trajectories for nine
station locations. Using the matrix method, there are nine air parcel
back trajectories run for each station, four times daily. The processing
time for the 324 back trajectories will vary from 10-20 min. Traj files
are stored in "processed_data/example_traj"

```{r, eval=FALSE}


# Formatted for use on hpc

args = commandArgs(trailingOnly = T)
#arg are expected in order
#1: the group run file to be iterated over - each line is one location and one day 
#2: the met dir "Hysplt/met_dir" where the .gbl files used for modeling backtrj are stored
#3: the output directory
#4: temp folder that provides space for parallel jobs to run hysplit models simultaneously

# example 
args = c("raw_data/locations.csv", "raw_data/met_dir", "processed_data/example_traj")

locations_in = read_csv(args[1])
met_dir_in = args[2] 
Out_dir = args[3] 

if (!length(args)==3) {
  stop('argumentes are expected in order
      1: the group run file to be iterated over - each line is one location and one day
      2: the met dir defalts to "Hysplt/met_dir"
      3: the base output directory Hysplt/processed_data/runs_$sys.date/'
       )
}
print(paste('Locations argument is:', args[1]))
print(paste('Meteorology argument is:', args[2]))
print(paste('Processed data directory argument is:', args[3]))

# simple for loop that takes in each row of the provided locations_in df
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

### Process trajectories

#### Daily rasters

```{r}
# obtain necessary functions
source("R/functions.R")

args = commandArgs(trailingOnly = T)
#arg are expected in order
#1: the directory that holds each day of trajectories -- this will be run from a csv.
# Each job accesses data in one folder that contains the trajectories for a given day for all points in a study domain. 

#2: the out directory provides a location to produce daily folders with HYSPLIT rasters for the study domain.

args = c("processed_data/example_traj", "processed_data/rasters")

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
