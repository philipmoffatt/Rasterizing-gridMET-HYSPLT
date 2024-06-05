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
    met_type = "reanalysis", model_height = 10000) %>%
  splitr::run_model()
```
`trajectory_model` is a list containing the trajectory data for the given dates in December 2023. Refer here for more details on the package and documentation:  https://github.com/rich-iannone/splitr

```{r}


start_date <- "2020-10-01"
end_date <- '2020-10-02'
locations_df <- read_csv(example_locs.csv)
met_dir_in <- 'Hysplt/met_dir'
exec_dir_in <- 'Hysplt/execute_dir'
```

Process trajectories 
  Daily rasters
  Cluster analysis 

Pull climate data


