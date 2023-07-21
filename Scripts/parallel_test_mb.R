# ---------------------------------------------------------
# TUM - Technichal University of Munich
#
# Authors:  Aydin Uzun, Taylor Jones
# Date: 2022
# Purpose:  PPB error calculation - This script should be run on the cluster. 
# Run : How to submit a job on the cluster that will run this R script?
#       A sample shell script to submit a job on SLURM is "job_parallel_test_mb.sh" and can be found on
#       '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Scripts'
#       Commands for interacting with SLURM:
#         submit a job
#             sbatch --clusters cm2_tiny job_parallel_test_mb.sh
#         check status
#             squeue --clusters cm2_tiny
#         cancel job
#             scancel --clusters cm2_tiny <job-id>
# Reference: Taylor Jones' transport error example can be found on '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Taylor_Jones_transport_error_example'
# ---------------------------------------------------------


############################ MODIFY THIS PART
############################ Day under investigation HAM_XXXXXXX 
# The working directory is "/dss/dsstumfs01/pn69ki/pn69ki-dss-0000/STILT_Ham/HAM_large/HAM_XXXXXXXX/by-id
# Download this script to the working directory and submit the job 
# XXXXXXXX being e.g. 20180905
# Enter sigma WSPD WDIR for the chosen day
# Declare standard deviations for the wind speed and wind direction uncertainties. 
# These numbers are calculated using the framework Aydin Uzun created
# The values for each day can be found on '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/Weighted_daily_WSPD_WDIR_differences'
# 'weighted_daily_wspd_wdir_differences.csv'
sigma_wind_speed <- 0.57 # meters/second
sigma_wind_dir   <- 12.909 # degrees
# Please remember that you for a single day there are 4 instruments available
# For faster parallel processing it is decided to create 4 different scripts for one day each dealing with the trajectories related to one specific location
# So you need to create 4 scripts for one day [just by changing the instrument_name] and submit 4 different jobs for one day
instrument_name <- "mb" # select the instrument name, it can be "mb","mc","md","me"
range_for_hist_cells <- 10 # define the maximum break point for the histograms cells
break_by <- 1 # define the distance between break points
# e.g. range_for_hist_cells <- 10 break_by <- 1 means that the bin borders are -10 ppb -9 ppb -8 ppb ......9 ppb 10 ppb
############################ MODIFY THIS PART


############################ MODIFICATION TO TAYLOR's ORIGINAL SCRIPT
# We asked him, how he makes sure that the rasterized particle file is in ppb,
# although the footprints we are used to working with are in [ppm*m2*s/umol].
# He approved that the default units of the particle files are ppm /  micro mol m^-2 s^-1
# He stated that we need to divide by the number of particles [in our case 500 particles] 
# after you sum up all the footprint values during the rasterize step. 
# Then we need to multiply by 1000 ( ppm -> ppb) 
############################ MODIFICATION TO TAYLOR's ORIGINAL SCRIPT


# Packages
require(raster)  #inventories and footprints use the raster library
require(geosphere) # just needed for the meters -> degrees lat and lon conversion
require(MASS)      # needed for the awesome 'mvrnorm' function
require(ncdf4) # package for netcdf manipulation
library(parallel) # for parallel processing
# set directory and set the names of the inventories
# The working directory is "/dss/dsstumfs01/pn69ki/pn69ki-dss-0000/STILT_Ham/HAM_large/HAM_XXXXXXXX/by-id
# XXXXXXXX being e.g. 20180905
setwd("./")
# names of the different inventories
# Inventories can be on found on '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/Inventories'
inventory_file_1 <- "HAM_CH4_inventory_city_and_river_smooth.nc"
inventory_file_2 <- "HAM_CH4_inventory_corrected_elev_bothInstruments_city_and_river_smooth.nc"
inventory_file_3 <- "HAM_CH4_inventory_corrected_raw_bothInstruments_city_and_river_smooth.nc"

#load in the 3 different inventories and then rasterize them
inventory_raster_stack_1 <- stack(inventory_file_1)
inventory_raster_stack_2 <- stack(inventory_file_2)
inventory_raster_stack_3 <- stack(inventory_file_3)
inventory_raster_1 <- inventory_raster_stack_1$X1 + inventory_raster_stack_1$X2  
inventory_raster_2 <- inventory_raster_stack_2$X1 + inventory_raster_stack_2$X2 
inventory_raster_3 <- inventory_raster_stack_3$X1 + inventory_raster_stack_3$X2  
# plot(inventory_raster_1)
# plot(inventory_raster_2)
# plot(inventory_raster_3)

# initialization for the directory to save the histograms
# check if there is a direcotry called e.g. tmp_mb exists if not create one to save the histograms
subDir <- paste("tmp_",instrument_name,"/",sep="")
ifelse(!dir.exists(file.path("./", subDir)), dir.create(file.path("./", subDir)), FALSE)

############################# start parallel processing
CORE_COUNT= detectCores()
print(paste("Detected", CORE_COUNT, "cores on the system"))

# find the trajectory files
pattern_for_list <- paste(instrument_name,".*\\_traj.rds$",sep="") # pattern for trajectory files
FILES = list.files(path = "./", pattern = pattern_for_list, recursive = TRUE)
FILES_COUNT = length(FILES)
print(paste("Found", FILES_COUNT, "trajectory files"))

# seperate the files into chunks for parallel processing 
if (FILES_COUNT < CORE_COUNT) { CORE_COUNT = FILES_COUNT }
CHUNKS = split(FILES, cut(seq_along(FILES), CORE_COUNT, labels = FALSE))
print(paste("Processing", FILES_COUNT, "files on", CORE_COUNT, "cores"))

# start calculating ppb error
process_filepath <- function(filepath) {
  xs <- strsplit(filepath, split = "/")[[1]]
  dirname <- xs[1]  #e.g. mb_202107300600_20m
  filename <- xs[2] #e.g. mb_202107300600_20m_traj.rds
  
  # load the particle file as a dataframe and then rasterize
  par <- data.frame(readRDS(filepath)$particle)
  foot_raster_1 <- rasterize(par[c('long', 'lati')], inventory_raster_1, par$foot, sum) 
  #plot(foot_raster_1)
  #plot(foot_raster_1*inventory_raster_1)
  
  
  index = par$indx
  particle_count <- max(index) # particle_count = 500
  # We can then multiply inventories by the footprint to get the gridded enhancement contributions in ppm
  # then divide by the number of particles and do the conversion from ppm to ppb
  # These are the total enhancements in ppb
  # Subscripts indicate that different inventories are used.
  original_ppb_1 <- (cellStats(foot_raster_1 * inventory_raster_1, sum)/particle_count) * 1000
  original_ppb_2 <- (cellStats(foot_raster_1 * inventory_raster_2, sum)/particle_count) * 1000
  original_ppb_3 <- (cellStats(foot_raster_1 * inventory_raster_3, sum)/particle_count) * 1000
  
  # Compute conversion factors from meters to degrees for latitude and longitude
  # Naturally these numbers will change depending where on the globe you are
  center_lat <- mean(par$lati)
  center_lon <- mean(par$long)
  y_meter_offset <- destPoint(c(center_lon, center_lat), b=0, d=1) # north 
  # b =  Bearing (direction) in degrees
  # d = Distance in meters
  # Given a start point, initial bearing (direction), and distance, this function computes the
  # destination point travelling along a the shortest path on an ellipsoid
  lat_deg_per_meter <- y_meter_offset[2] - center_lat # ~=1.469224e-05
  x_meter_offset <- destPoint(c(center_lon, center_lat), b=90, d=1) # east
  lon_deg_per_meter <- x_meter_offset[1] - center_lon # ~= 8.986748e-06
  
  ########################## Create wind errors that vary in time
  # We'll vary the wind direction and wind speed throughout the back trajectory
  # The errors need to change smoothly across timesteps
  # E.g. wind direction errors of −10deg and +8deg are both reasonable if we asssume sigma_wdir=5deg, 
  # but the wind shouldn’t be off by −10deg at t=−5 minutes and then be off by +8deg at t=−6 minutes.
  # It is therefore useful to define correlation time scales for wind speed and wind direction errors.
  # Taylor usually choses a value of 3 hours, which is based on a cursary look at a few wind spectral signiture studies. 
  # These values should be adjusted if you know more about wind variability in your domain. 
  # We decided to keep this value the same.
  ws_error_tscale <- 3 * 60 # minutes
  wd_error_tscale <- 3 * 60 # minutes
  t <- unique(par$time) #note t is also in minutes -1 to -840 decreasing by 1, i.e. t num[1:840] = -1 -2 -3 -4 .... -840
  nt <- length(t) # nt=840
  # initialize wind error covariance matrices.
  ws_error_cov <- matrix(nrow = nt, ncol = nt, data = 0) #of size num[1:840,1:840] hence a large matrix, in Taylor's case num[1:24,1:24]
  wd_error_cov <- matrix(nrow = nt, ncol = nt, data = 0) #of size num[1:840,1:840] hence a large matrix, in Taylor's case num[1:24,1:24]
  # fill wind error covariance matrices
  for (i in 1:nt) {
    for (j in 1:nt) {
      ws_error_cov[i, j] <-
        (sigma_wind_speed ^ 2) * exp(-1 * abs(t[i] - t[j]) / ws_error_tscale)
      wd_error_cov[i, j] <-
        (sigma_wind_dir ^ 2) * exp(-1 * abs(t[i] - t[j]) / wd_error_tscale)
    }
  }
  # image(ws_error_cov)
  # image(wd_error_cov)
  
  simulation_count <- 400 # number of simulations to run. simulation_count should be sufficiently large.
  
  # initialization for ppb errors
  # Subscripts corresponds to different inventories that is used
  ppb_errors_1 <- c() 
  ppb_errors_2 <- c() 
  ppb_errors_3 <- c() 

  
  # 400 simulations for each particle file
  for (i in 1:simulation_count) {
    
    # Wind speed and direction errors that vary with time can be pulled easily
    # thanks to the amazing function MASS::mvrnorm, which does exactly that.
    ws_error <-   mvrnorm(mu = rep(0, nt), Sigma = ws_error_cov) # of size num[1:840]
    wd_error <- mvrnorm(mu = rep(0, nt), Sigma = wd_error_cov) # of size num[1:840]
    # plot(t,ws_error,type="b") # to see how the wind direction error is correlated with time
    # plot(t,wd_error,type="b") # to see how the wind direction error is correlated with time
    
    # compute errors' x and y components
    ws_x_error <- ws_error * sin(wd_error * pi / 180) # m/s
    ws_y_error <- ws_error * cos(wd_error * pi / 180) # m/s
    
    
    #cumulative errors:
    # The cumulative displacement that each particle experience is solely a function of time, 
    # as all particles are experiencing the same wind error within a time step.
    cumm_x_error <-  cumsum(ws_x_error * 60 * c(diff(t)[1] , diff(t))) # of size num[1:840] # meters
    cumm_y_error <-  cumsum(ws_y_error * 60 * c(diff(t)[1] , diff(t))) # of size num[1:840] # meters
    #plot(t,cumm_x_error,type="b") # to see cummulative displacement in the x-direction in meters
    #plot(t,cumm_y_error,type="b") # to see cummulative displacement in the x-direction in meters
    
    # interpolation
    x_errors <- approx(t, cumm_x_error, xout = par$time)$y # of size num[1:420000]
    y_errors <- approx(t, cumm_y_error, xout = par$time)$y # of size num[1:420000]
    # new dataframe
    shifted_part <- par
    # add displacement to every particle to compute new particle positions
    shifted_part$lati <- par$lati + lat_deg_per_meter * y_errors
    shifted_part$long <- par$long + lon_deg_per_meter * x_errors
    # rasterize the new dataframe
    shifted_foot_raster <-rasterize(shifted_part[c('long', 'lati')], inventory_raster_1, shifted_part$foot, sum)
    # multiply by the inventory, cell values summed, divided by the number of particles and then multiplied by 100 to compute the new ppb value
    # ppb enhancements for shifted particle
    ppb_1 <-   (cellStats(shifted_foot_raster * inventory_raster_1, sum)/particle_count )*1000
    ppb_2 <-   (cellStats(shifted_foot_raster * inventory_raster_2, sum)/particle_count )*1000
    ppb_3 <-   (cellStats(shifted_foot_raster * inventory_raster_3, sum)/particle_count )*1000
    # e.g. ppb_1 - original_ppb_1 is now the ppb error for this simulation run
    # but we'll repeat this process a large number of times to get a distribution of ppb errors
    # for each simulation run [i in 1:simulation_count] calculate the ppb error as ppb_1 - original_ppb_1 and in the end concatanate them
    # Again please note that subscripts correspond to the different inventories used 
    ppb_errors_1 <- c(ppb_errors_1, ppb_1 - original_ppb_1) # in the end will be of size num [1:400]
    ppb_errors_2 <- c(ppb_errors_2, ppb_2 - original_ppb_2) # in the end will be of size num [1:400]
    ppb_errors_3 <- c(ppb_errors_3, ppb_3 - original_ppb_3) # in the end will be of size num [1:400]
    
  }
  # now we have the ppb errors for the number of simulations (simulation_count) we did
  # we wanted to convert these simulation_count ppb errors to a distribution
  # e.g. range_for_hist_cells <- 10, break_by <- 1 means that the bin borders are -10 ppb -9 ppb -8 ppb ......9 ppb 10 ppb
  # currently ignoring if ppb_errors > range_for_hist_cells or ppb_errors< -range_for_hist_cells
  hist_ppb_errors_1 <- hist(ppb_errors_1[ppb_errors_1 >= -range_for_hist_cells & ppb_errors_1 < range_for_hist_cells],breaks=seq(-range_for_hist_cells,range_for_hist_cells,by=break_by))
  hist_ppb_errors_2 <- hist(ppb_errors_2[ppb_errors_2 >= -range_for_hist_cells & ppb_errors_2 < range_for_hist_cells],breaks=seq(-range_for_hist_cells,range_for_hist_cells,by=break_by))
  hist_ppb_errors_3 <- hist(ppb_errors_3[ppb_errors_3 >= -range_for_hist_cells & ppb_errors_3 < range_for_hist_cells],breaks=seq(-range_for_hist_cells,range_for_hist_cells,by=break_by))
  
  # save the counts of the histograms in each cell
  # e.g. range_for_hist_cells <- 10, break_by <- 1, the first element of hist_ppb_errors_1$counts corresponds to the number of simulations out of simulation_count simulations which resulted into ppb errors in the range [-10 ppb, -9 ppb]
  saveRDS(hist_ppb_errors_1$counts, paste(subDir, dirname, "_error_1_hist.rds", sep=""))
  saveRDS(hist_ppb_errors_2$counts, paste(subDir, dirname, "_error_2_hist.rds", sep=""))
  saveRDS(hist_ppb_errors_3$counts, paste(subDir, dirname, "_error_3_hist.rds", sep=""))
  
}

# Function to process on chunk
process_chunk <- function(chunk_no) {
  lapply(CHUNKS[[chunk_no]], process_filepath)
}
# Process one chunk per core
print("Running simulations ...")
system.time(
  mclapply(seq(1, CORE_COUNT), process_chunk, mc.cores = CORE_COUNT)
)


print("Merging all histogram files into one ...")
# Initialize the merged histograms
hist_ppb_errors_1_counts <- c()
hist_ppb_errors_2_counts <- c()
hist_ppb_errors_3_counts <- c()
# save the names of the trajectory files
# The name contains the information about the instrument, exact datetime and particle release height
names_of_df <- c()

for (chunk in CHUNKS) {
  for (filepath in chunk) {
    xs <- strsplit(filepath, split = "/")[[1]]
    dirname <- xs[1]  # mb_202107300600_20m
    names_of_df <- c(names_of_df, dirname)
    
    hist_ppb_errors_1_counts <- cbind(
      hist_ppb_errors_1_counts,
      readRDS(paste(subDir, dirname, "_error_1_hist.rds", sep=""))
    )
    hist_ppb_errors_2_counts <- cbind(
      hist_ppb_errors_2_counts,
      readRDS(paste(subDir, dirname, "_error_2_hist.rds", sep=""))
    )
    hist_ppb_errors_3_counts <- cbind(
      hist_ppb_errors_3_counts,
      readRDS(paste(subDir, dirname, "_error_3_hist.rds", sep=""))
    )
  }
}

paste("parallel_",instrument_name,"_hist_ppb_errors_1_counts.rds", sep="")
# save the results we need for further steps
# e.g. hist_ppb_errors_1_counts contains the counts of the histograms of ppb errors for the chosen day and for the chosen instrument
# Again subscript 1,2,3 corresponds to the different inventory files used. 
saveRDS(hist_ppb_errors_1_counts, paste("parallel_",instrument_name,"_hist_ppb_errors_1_counts.rds", sep="")) 
saveRDS(hist_ppb_errors_2_counts, paste("parallel_",instrument_name,"_hist_ppb_errors_2_counts.rds", sep="")) 
saveRDS(hist_ppb_errors_3_counts,paste("parallel_",instrument_name,"_hist_ppb_errors_3_counts.rds", sep="")) 
saveRDS(names_of_df,paste("parallel_",instrument_name,"_particle_files.rds", sep=""))
# These histogramcounts will be later used to calculate a daily histogram
print("Done!")