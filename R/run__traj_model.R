run_trajectory_model <- function(start_date, end_date, 
                                 locations_df, 
                                 met_dir_in) { # add in execu and make traj files in daily folders with line 78
  trajectory_list <- list()
  
  # Iterate over each location in the locations dataframe
  for (j in 1:nrow(locations_df)) {
    lat <- locations_df$latitude
    lon <- locations_df$longitude
    
    # Get the StationNumber associated with the current lat and lon
    station <- locations_df$station[j]
    
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
    # 
    # # Get trajectory output as a data frame for the current location
    trajectory_df <- trajectory_model %>% get_output_tbl()
    # 
    # # Add a column to the dataframe with the corresponding StationNumber
    trajectory_df$station <- station
    # 
    # Append the named trajectory data frame to the list
    trajectory_list[[j]] <- trajectory_df
  }
  
  # Return the list of named trajectory data frames for all locations
  return(trajectory_list)
}


args = commandArgs(trailingOnly = T)
#arg are expected in order
#1: the group run file to be iterated over - each line is one location and one day 
#2: the met dir "Hysplt/met_dir" where the .gbl files used for modeling backtrj are stored
#3: the output directory
#4: temp folder that provides space for parallel jobs to run hysplit models simultaneously

# example 
args = c("raw_data/locations.csv", "raw_data/met_dir", "processed_data/example_traj", paste0("C:/Users/philip.moffatt/Documents/GitHub/Rasterizing-gridMET-HYSPLT/processed_data/exec_dir/runs_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S")))

locations_in = read_csv(args[1])
met_dir_in = args[2] 
Out_dir = args[3] 
exec_dir = args[4]

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