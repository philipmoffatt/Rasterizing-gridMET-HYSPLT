

#' @title Run HYSPLIT
#' @description Runs the HYSPLIT model
#' @param start_date the begining of the time period to model back trajectories over (character)
#' @param end_date end of the period (character)
#' @param locations_df a dataframe consisting of one row and four columns:  latitude (num), longitude (num), station (num), date(Date yyyy-mm-dd)
#' @param met_dir_in file path to meterological data used to run HYSPLIT (character), an unnecessary param but it saves considerable computation time. Download .gbl files from 'ftp://arlftp.arlhq.noaa.gov/archives/reanalysis' 
#' @return trajectory dataframe
#' 
run_trajectory_model <- function(start_date, end_date, 
                                 locations_df, 
                                 met_dir_in) {
  # Debugging statements
  print("Starting run_trajectory_model")
  print(paste("Start date:", start_date))
  print(paste("End date:", end_date))
  print(paste("Met directory:", met_dir_in))
  print("Locations DataFrame:")
  print(locations_df)
  
  # Iterate over each location in the locations dataframe
  lat <- locations_df$latitude
  lon <- locations_df$longitude
  
  # Get the StationNumber associated with the current lat and lon
  station <- locations_df$station
  
  # Create trajectory model for the current location
  trajectory_model <- create_trajectory_model() %>%
    # builds a matrix of points
    add_grid(
      lat = lat,
      lon = lon,
      range = c(0.28, 0.28),
      division = c(0.14, 0.14))  %>%
    # load selected parameters
    add_trajectory_params(
      height = 1000, # m agl
      duration = 120, # hours
      days = seq(as.Date(start_date), as.Date(end_date), by = "day"),
      daily_hours = c(0, 6, 12, 18),
      direction = "backward",
      met_type = "reanalysis", 
      met_dir = met_dir_in,
      model_height = 10000) %>% # trajectory upper limit, stops traj before duration
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
  
  test$location = H %>% 
    str_split("_") %>% 
    sapply("[",3) %>% # this may change depending on file paths
    as.numeric()  #
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
  
  test[,dist_test := pmap(l, 
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
  test_sf_list = map(test_it_all, 
                     walk_traj_function)
  # summarizing each trajectory with sum all function
  test_sum_all = map_df(test_sf_list, 
                        sum_all_fn)
  # summarizing each trajectory for land based calcs
  test_sum_land = map_df(test_sf_list, 
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
  
  print(paste("the first file in the data path ",data_in_file[[1]]))
  
  temp_date = data_in_file[[1]] %>% 
    str_split("/") %>% sapply("[",7)
  
  print(paste("the folder being processed is ",temp_date))
  
  tic()
  traj_back = purrr::map(data_in_file,
                         function(H) backtrj_summarize_k(H = H, north_america = north_america))
  toc()
  print("traj_back has run")
  
  #### 4) make raster ######
  
  traj_back_df_filled = traj_back %>% rbind.fill()
  data_raster1 <- rasterFromXYZ(traj_back_df_filled %>%
                                  mutate(x = lon, y = lat) %>%
                                  dplyr::select(x,y, eval(names(traj_back_df_filled)[5:28])),
                                crs = crs("epsg:4326"))
  
  print("traj has been rasterized")
  
  temp_out_path = file.path(out_path, paste0("r_",temp_date))
  terra::writeRaster(data_raster1, paste0(temp_out_path, '.tif'), overwrite = T)
  print(paste0("Trajectory raster built for ", temp_date,".")) # the "." may not work on kamiak
}







