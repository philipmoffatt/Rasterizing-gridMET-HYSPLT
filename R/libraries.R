# necessary libraries for the entire methdology

library(tidyverse)

# HYSPLT trajectories
library(splitr)

# gridMET, GHCN weather
library(climateR)

# sounding data Univerisity of Wyoming
library(climate)

# trajectory cluster analysis 
library(openair)

# maps and figures
library(maps)
library(ggrepel)
library(terra)
library(ggrepel)
library(ggmap)
library(cowplot)

# converting HYSPLIT trajectory data to daily rasters
library(terra)
library(plyr)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(geosphere)
library(data.table)
library(raster)
library(leaflet)
library(tictoc)
library(furrr)

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
