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
    add_grid(
      lat = lat,
      lon = lon,
      range = c(0.28, 0.28),
      division = c(0.14, 0.14))  %>%
    add_trajectory_params(
      height = 1000,
      duration = 120,
      days = seq(as.Date(start_date), as.Date(end_date), by = "day"),
      daily_hours = c(0, 6, 12, 18),
      direction = "backward",
      met_type = "reanalysis", 
      met_dir = met_dir_in,
      model_height = 10000) %>%
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

# Split the data frame by both `station` and `date`
# test_list <- locations_in %>%
#   group_by(station, date) %>%
#   group_split()
# print(paste('The first station in the test list is ',test_list[[1]]$station))


# apply the `run_trajectory_model` to each element in `test_list` - that is each day and location combo
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
    end_time = Sys.time() # log the end time - should be around 15-20 seconds
    
    print(paste("total time per station =", end_time - start_time, "seconds!"))
    
    write_rds(traj_temp, file = file.path(Out_dir_traj,paste0("traj_",H$station))) # write the trajectory file with station name into daily folder
  }
warn