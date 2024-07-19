#' Convert HYSPLIT trajectory data into spatial rasters
#'
#' The function summarizes trajectory data into daily indices and produces daily folders of rasters for a study domain with layers that show HYSPLIT indices
#' @param data_in_path Folder with trajectory files for a given day. The folder is expected to be named for the date (/2021-01-01).
#'   
#'   `sum_all_fn`
#' @param show_hourly An option to show hourly positions and associated data
#'   along trajectories.
#' @param color_scheme Defines the appearance of multiple trajectories in a
#'   single plot. Current options are `cycle_hues` (the default), and
#'   `increasingly_gray`.

#' Example usage:
# args = c("Hysplt/processed_data/2020-10-03", "Hysplt/processed_data/rasters")
# 
# data_in_path = args[1]
# out_path = args[2]
# north_america = get_northAmerica()
# summarize_to_raster(data_in_path = data_in_path,
#                       out_path = out_path)
# # 
# # read and view the first raster from memory
# # Specify the path to the
# grd_file_path <- file.path('Hysplt/processed_data/rasters/r_NA.tif')
# 
# # Read the raster from the .grd file
# first_raster <- terra::rast(grd_file_path)
# 
# # Plot the raster
# plot(first_raster[[1]], main = "Upwind Temp Range")

#### 1) functions #####



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
} ###speed this up! 
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
walk_traj_function= function(traj_in = traj_in){
  traj_in$LC_group =1
  temp = 1
  freeze_cycles = 0
  sat_cycles = 0
  groupslc <- for(i in 2:nrow(traj_in)){
    if(traj_in$landcover[i]!=traj_in$landcover[i-1]){
      temp = temp +1
    }
    
    if(traj_in$air_temp[i]<273 &traj_in$air_temp[i-1]>273){
      freeze_cycles = freeze_cycles+1 
    }
    
    if(traj_in$rh[i]>95 & traj_in$rh[i-1]<95){
      sat_cycles = sat_cycles+1 
    }
    
    traj_in$LC_group[i] = temp  
  }
  traj_in = traj_in %>% mutate(
    freeze_cycles = freeze_cycles, 
    sat_cycles = sat_cycles)
  
  
  return(traj_in)
}

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

#### 3) summarize back trajectory function ######


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
  
  #### 4) make raster ######
  
  traj_back_df_filled = traj_back %>% rbind.fill()
  data_raster1 <- rasterFromXYZ(traj_back_df_filled %>%
                                  mutate(x = lon, y = lat) %>%
                                  dplyr::select(x,y, eval(names(traj_back_df_filled)[5:28])),
                                crs = crs("epsg:4326"))
  
  print("traj has been rasterized")
  
  temp_out_path = file.path(out_path, paste0("r_",temp_date))
  terra::writeRaster(data_raster1, paste0(temp_out_path, '.tif'), overwrite = T)
  print(paste0("Trajectory raster built for ", temp_date,"."))
}
