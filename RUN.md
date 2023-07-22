# Wind Error Analysis in the Inversion Framework

## Lidar Measurements

The lidar measurements used in this study were collected during the Hamburg Campaign of 2021. The data can be accessed at the following location: `./esm/projects/WRFDA_Hamburg/202108_lidar_wettermast_hamburg/DLR85`.

To process the lidar measurements, we use the script [read_lidar_script.m](./Scripts/read_lidar_script.m). This script reads all the lidar measurements and saves the data as matrices for later use. The saved variables can be found at [Lidar](./Data/Lidar).

The lidar measurements have a fine temporal and vertical resolution, with one measurement recorded every 10 minutes. Notably, these measurements do not occur exactly at sharp hours (e.g., xx:55, xx:05, xx:15, xx:25, xx:35, xx:45, etc., where xx is the hour of the day). On the other hand, the available ERA5 model data provides hourly data. To make the lidar measurements compatible with the ERA5 data, we perform data conversion using the [convert_lidar_to_hourly.m](./Scripts/convert_lidar_to_hourly.m) script. The converted data is then saved at: [Lidar hourly](./Data/Lidar_hourly).

## ERA5 Model Data

The next step involves obtaining the model data from ERA5 and combining it.

**ERA5 Hourly Data on Single Levels:**
The ERA5 hourly data on single levels for the time domain of interest is located at: `./esm/projects/WRFDA_Hamburg/ERA5/ERA5-sl`. The script [read_era5_single_levels_script.m](./Scripts/read_era5_single_levels_script.m) reads this data, and the saved variables can be accessed at: [ERA5 sl](./Data/ERA5-sl).

**ERA5 Hourly Data on Pressure Levels:**
Similarly, the ERA5 hourly data on pressure levels for the time domain of interest can be found at: `./esm/projects/WRFDA_Hamburg/ERA5`. The script [read_era5_pressure_levels_script.m](./Scripts/read_era5_pressure_levels_script.m) reads this data, and the saved variables can be found at: [ERA5 pl](./Data/ERA5-pl).

**Combination Run:**
To integrate the ERA5 data from pressure levels with single levels, execute the script [integrate_era5_pressure_levels_with_single_levels.m](./Scripts/integrate_era5_pressure_levels_with_single_levels.m). The output variables are matrices that represent the combined ERA5 model, and they can be accessed at: [ERA5 integrated](./Data/ERA5-integrated).

## Model-Data mismatch (Windspeed and direction)

The script [interpolation_and_error_processing.m](./Scripts/interpolation_and_error_processing.m) takes as input the previously combined ERA5 model data and hourly converted lidar data. Next, since the vertical resolutions of the model data and the measured data are different and we want to make a comparison, they need to look similar, i.e. have the same vertical and vertical resolution. To fix the temporal resolution the range between 06:00 and 18:00 is chosen, since footprints are present at these times. Temporal coverage is the 40-day interval from August 1, 2021 to September 9, 2021. 

Regarding the vertical coverage of lidar measurement, the measurement height above ground level ranges from ~25m to ~2900m with about one sample per 15-25 m. Overall we have WSPD and WDIR values at 149 altitude values, which provides a quite fine resolution. However, this is not the case for the ERA5 model data. WSPD and WDIR values are only available at 39 altitude values. 

First, 13 particle release heights are defined in this script [interpolation_and_error_processing.m](./Scripts/interpolation_and_error_processing.m). Each of these 13 particle release heights defines a layer and the layer boundaries are the average of the two consecutive layers. The lidar measurements and model data (WSPD and WDIR) are interpolated such that there are 10 linearly spaced points in each atmospheric layer. Then the average [in terms of WSPD and WDIR] of these 10 linearly spaced points is taken and said to be representative of that layer. This operation is performed for each atmospheric layer, for each time series and for both lidar measurements and ERA5 data. This gives us as output the representation of the whole atmosphere during the campaign in terms of WSPD and WDIR for both lidar measurements and the ERA5 model. In addition, the mismatches are defined as the distance between the Lidar representation of the atmosphere and the ERA5 model representation of the atmosphere. Both the lidar representations and ERA5 model representations [WSPD and WDIR] and the mismatches between them are saved for later use and can be found at: [Lidar ERA5 representatives](./Data/Lidar_ERA5_representatives).

In addition, one can also generate WSPD or WDIR plots for a specific datetime that visualize the raw data, interpolated data and the calculated representatives. The WSPD plots that visualize the raw data, interpolated data and the representatives for the entire capaign can be found at: [WSPD figures](./Data/WSPD_figures).

The WDIR plots that visualize the raw data, interpolated data and the representatives for the entire capaign can be found at: [WDIR figures](./Data/WDIR_figures)

The script [wspd_wdir_differences_image.m](./Scripts/wspd_wdir_differences_image.m) takes as input the calculated representative WSPD and WDIR mismatches [both 3D matrices of size 40 (representing the days) x 13 (representing the hours from 06:00 to 18:00) x 13 (atmospheric layer)] and then create images with scaled colors for each day to visualize the mismatches between the measurements and the model. The output images can be found at: [WSPD WDIR differences images](.Data/WSPD_WDIR_differences_images).

## How and why to weight the differences?

We know that the receptors at higher altitudes are not sensitive to surface emissions. So mismatches in the higher atmospheric layers should actually be less important.  We need to somehow weight these mismatches that we calculated earlier. The idea was to use the footprint intensities and the vertical scaling factors as weights. 

### Pressure scaling factors

The script [read_scaling_factors.R](./Scripts/read_scaling_factors.R) reads the scaling factors from receptors' .rds documents, which can be found at: `/./esm/data/Footprint_STILT/Hamburg/ERA5/V2_small`

Since there are 13 atmopsheric layers, there are 13 different scaling factors for each day. The output is the 2D scaling matrix for the entire campaign and the daily scaling factors and can be found at: [Scaling factors](.Data/Scaling_factors).

The simplest script [convert_scaling_factors.m](./Scripts/convert_scaling_factors.m) reads the .csv document from the previous R script and normalizes it and saves as a .mat document for later use. The output can again be found in the same folder.

### Combine with footprint intensities to get weighting factors

The footprints exist for each instrument location [mb,mc,md,me], particle release height[13 different values], time series triple and are located on the cluster: `/dss/dsstumfs01/pn69ki/pn69ki-dss-0000/STILT_Ham/HAM_large`

The script [CLUSTER_calculate_footprint_intensities_with_scaling_factors.m](./Scripts/CLUSTER_calculate_footprint_intensities_with_scaling_factors.m) needs to be uploaded and run at this location and takes as input the footprints and the normalized scaling matrix, so the previously calculated normalized scaling matrix needs to be uploaded here as well.

The script [CLUSTER_calculate_footprint_intensities_with_scaling_factors.m](./Scripts/CLUSTER_calculate_footprint_intensities_with_scaling_factors.m) calculates the footprint intensities at sharp hours available and provides as output a csv document where the first column corresponds to the date and time, the second column corresponds to the label of the instrument and the third to last column correspond to the factors for the 13 layers. These factors are calculated by multiplying the available scaling factor (from *.rds files) for the corresponding day by the footprint intensities (basically the sum of all the elements of the 2D footprint image) on each layer. The third column corresponds to the first layer, the fourth column corresponds to the second layer, and so on. By the way, the factors are normalized, so the sum in each row is 1. The output CSV document can be found at [HAM_large_footprint_intensities_with_scaling_factors](./Data/HAM_large_footprint_intensities_with_scaling_factors).

### Visualize weighting factors

The previous CVS document containing the weighting factors for each atmospheric layer, instrument and time series triple is used as input to the script [weighting_factors_image.m](./Scripts/weighting_factors_image.m) to display the weighting factors as images with scaled colors similar to we did it in the script  [wspd_wdir_differences_image.m](./Scripts/wspd_wdir_differences_image.m). Before the images for ech day are generated, weighting factors are averaged over instruments (mb,mc,md,me) and normalized again.

The output images can be found at: [Weighting_factors_images](./Data/Weighting_factors_images).

## Weighting the representatives and their differences

Using the weighting factors we created and visualised earlier, we can weight the mismatches betweem the LiDAR representations and ERA5 model representations for both WSPD and WDIR. Please remember that we have already calculated the mismatches and saved them to [Lidar_ERA5_representatives](./Data/Lidar_ERA5_representatives). The script [plot_weighted_histograms.m](./Scripts/plot_weighted_histograms.m) generates daily weighted histogram plots using the beneficial function [histwc.m](./Scripts/histwc.m). [histwc.m](./Scripts/histwc.m) is a MATLAB function written by [Mehmet Suzen](https://www.mathworks.com/matlabcentral/fileexchange/42493-generate-weighted-histogram) on MATLAB Central File Exchange to generate weighted histograms. In addition, the mean and standard deviation of the daily mismatches are again calculated with the weighting factors and highlighted in the plots. 

The weighted histograms for the entire campaign can be found at: [Weighted_histograms](./Data/Weighted_histograms).

Using the weighting factors created and visualised earlier, we can weight the mismatches betweem the LiDAR representatives and ERA5 model representatives for both WSPD and WDIR in order to find the mean and standard deviation values for the daily mismatch. This must be done for both WSPD and WDIR. The script that does this is called [calculate_weighted_daily_differences.m](./Scripts/calculate_weighted_daily_differences.m). The output document contains the mean WSPD, WDIR mismatches and standard deviation of WSPD, WDIR mismatches for each day and is named [weighted_daily_wspd_wdir_differences.csv](./Data/Weighted_daily_WSPD_WDIR_differences/weighted_daily_wspd_wdir_differences.csv).

Similarly, the WSPD, WDIR lidar and ERA5 representatives are used as input not the mismatches, again saved earlier in [Lidar_ERA5_representatives](./Data/Lidar_ERA5_representatives) with the weighting factors previously created and visualised for the script [calculate_weighted_daily.m](./Scripts/calculate_weighted_daily.m) to calculate the mean and standard deviation values for the daily WSPD and WDIR values for both the lidar measurements and the ERA5 model. The output document is named [weighted_daily_wspd_wdir.csv](./Data/Weighted_daily_WSPD_WDIR/weighted_daily_wspd_wdir.csv).

## Use mean wind direction mismatch to rotate the footprints

The (aggregated) column foortprints can be found at `./esm/campaigns/2021AugHamburg/Data/Footprints/foot`. Using the daily mean WDIR information available in [weighted_daily_wspd_wdir_differences.csv](./Data/Weighted_daily_WSPD_WDIR_differences/weighted_daily_wspd_wdir_differences.csv) (created in the previous part), the script [rotate_footprint.m](./Scripts/rotate_footprint.m) rotates the aggregated footprints. The rotation must be done around the exact position of the instrument [mb,mc,md,me] on the grid. [rotate_footprint.m](./Scripts/rotate_footprint.m) uses the simple function [find_instrument_index_in_mesh.m](./Scripts/find_instrument_index_in_mesh.m) to find the latitude and longitude index of the instrument location using the available latitude and longitude grid and the exact location of the instrument. The output i.e. the latitude and longitude index of the instrument location defines around which point the rotation must be performed. Then the function [rotate_around.m](./Scripts/rotate_around.m) is used to find the rotated footprint image using the rotation angle information (mean WDIR mismatch), the footprint to be rotated and the point around which the rotation needs to be completed. This function is a modified version of the function `rotateAround.m` by [Jan Motl](https://www.mathworks.com/matlabcentral/fileexchange/40469-rotate-an-image-around-a-point) on MATLAB Central File Exchange.

The rotated footprints are saved to: `./esm/campaigns/2021AugHamburg/Data/Footprints/foot_rotated`

## PPB Error calculation 

We have calculated the standard deviations for the wind speed and wind direction uncertainties and please remember they can be found at: [Weighted_daily_WSPD_WDIR_differences](./Data/Weighted_daily_WSPD_WDIR_differences)

Now we will use these values to perform an error simulation to end up finding a transport error according to the Taylor's framework. For Taylor Jones' transport error example, contact [Taylor Jones](https://www.linkedin.com/in/taylor-jones-91568264/).

### Overview of Taylor's transport error calculation method

First he defines standard deviations for the wind speed and wind direction uncertainties. Then the inventory is loaded as a raster layer, then the particle file is loaded as a dataframe, and then the footprint is rasterized. It then multplies the inventory raster by the footprint to get the gridded enhancement contribution in ppb, and these contributions can be summed up to get the total enhancement in ppb. Let's keep this value. 

The error simulation process must be repeated a large number of times, each time selecting new errors. 

The idea was to generate wind errors that vary in time by building wind error covariance matrices. The back trajectory of the particle file and the standard deviations that define the uncertainties are used to create these wind error covariance matrices. At each repretition to select new errors, wind speed and direction errors that vary with time, but are correlation according to these matrices are drawn. Then new particle positions are calculated by adding the displacement (which is a result of the error) to each particle. The new dataframe is created with the new particle positions. Now, similar to the first part, the new dataframe is rasterized and multiplied by the inventory raster to calculate the new total enhancement in ppb. The distance between the new total enhancement in ppb and the original total enhancement in ppb is defined as ppb error. 

As mentioned, this process is repeated a large number of times and he gets a distribution of ppb errors. He defines the standard deviation of this distribution as the transport error statistic.

### Modification to Taylor's original script

We asked him, how he makes sure that the rasterized particle file is in ppb, although the footprints we are used to working with are in [ppm*m2*s/umol]. He approved that the default units of the particle files are ppm /  micro mol m^-2 s^-1. He stated that we need to divide by the number of particles [in our case 500 particles] after we sum up all the footprint values during the rasterize step. Then we need to multiply by 1000 for the conversion from ppm to ppb. 

Apart from that, for one day we have 4 instrument locations [mb,mc,md,me] x 45 time points [06:00 to 17:00 every 15 min] x 13 atm. layers. This equates to more than 2000 particle files, so Taylor's code needs to run on ~2000 particle files. The idea was to seperate the jobs by site [mb,mc,md,me] to allow faster parallel processing.

E.g. [parallel_test_mb.R](./Scripts/parallel_test_mb.R) is the script related to instrument "mb". Since there are a total of ~2000 particle files for a day, ~500 of them relate to the instrument "mb". This means that we have ~500 different distributions as output after running the code, each relating to different time series and atmospheric layers. For standardization and ease of use, these histograms are stored as the number of counts of their bins. The histogram of a single error simulation run [single particle file] is defined by two parameters: [range_for_hist_cells](./Scripts/parallel_test_mb.R), which defines the maximum break point for the histogram cells, and  [break_by](./Scripts/parallel_test_mb.R), which defines the distance between break points. These parameters are freely selectable, but in our work they are set to [range_for_hist_cells=10](./Scripts/parallel_test_mb.R) and [break_by=1](./Scripts/parallel_test_mb.R). With this we predefine the cells of all histograms.


After running the Taylor's error simulation for a day on all particle files, we get ~2000 different histograms i.e. the number of observations falling in the cells we predefined. The goal is to get a single daily histogram. Imagine we have [4 (sites) x 45 (time series) x 13 (atm. layers)] different histograms. The counts corresponding to a pair of site [4], time series [45] can be calculated by combining the counts corresponding to different atm. layers using the pressure scaling factors. This gives us [4 (sites) x 45 (time series)] different histograms and needs to be averaged over time series [45] and then over sites [4] to get a daily distribution. The script [combine_ppb_error_histograms.R](./Scripts/combine_ppb_error_histograms.R) does this job. The saved daily histograms for the days that interest us can be found at: [Daily_histograms](./Data/Daily_histograms/).

In our work, 3 different inventories are used that describe the same area and can be found on: [Inventories](./Data/Inventories)

The subscripts 1,2,3 in the scripts [parallel_test_mb.R](./Scripts/parallel_test_mb.R) and [combine_ppb_error_histograms.R](./Scripts/combine_ppb_error_histograms.R) represent different inventories used: 1 for [HAM_CH4_inventory_city_and_river_smooth.nc](./Data/Inventories/HAM_CH4_inventory_city_and_river_smooth.nc),2 for[HAM_CH4_inventory_corrected_elev_bothInstruments_city_and_river_smooth.nc](./Data/Inventories/HAM_CH4_inventory_corrected_elev_bothInstruments_city_and_river_smooth.nc), 3 for [HAM_CH4_inventory_corrected_raw_bothInstruments_city_and_river_smooth.nc](./Data/Inventories/HAM_CH4_inventory_corrected_raw_bothInstruments_city_and_river_smooth.nc).




### How to run? - [parallel_test_mb.R](./Scripts/parallel_test_mb.R)

The working directory is `/dss/dsstumfs01/pn69ki/pn69ki-dss-0000/STILT_Ham/HAM_large/HAM_XXXXXXXX/by-id` (XXXXXXXX is the day, e.g. 20180905). The jobs will be seperated by the site [mb,mc,md,me] for faster parallel processing. So there should be four R scripts in the working directory, each belonging to a site. E.g. `parallel_test_mb.R` is the one related to the instrument "mb". You need to copy paste this script and create three more named `parallel_test_mc.R` , `parallel_test_md.R`, `parallel_test_me.R`. In the [MODIFY THIS PART](./Scripts/parallel_test_mb.R) section, you need to change the `instrument_name` variable when you create the additional ones. The three different inventories mentioned earlier must also be in the working directory.  Apart from that, one can change the previously mentioned parameters `range_for_hist_cells` and `break_by` but in our work we keep them constant. For the parameters `sigma_wind_speed` and `sigma_wind_dir` please return to [weighted_daily_wspd_wdir_differences.csv](./Data/Weighted_daily_WSPD_WDIR_differences/weighted_daily_wspd_wdir_differences.csv).  You need to manually set these parameters for the desired day, which is XXXXXXXX from HAM_XXXXXXXX.

#### Submit a job on the cluster

An example shell script to submit a job on SLURM is [job_parallel_test_mb.sh](./Scripts/job_parallel_test_mb.sh):
```shell
#!/bin/bash

# ---------------------------------------------------------
# TUM - Technichal University of Munich
#
# Authors:  Aydin Uzun
# Date: 2022
# Purpose:  Run the R script on SLURM calculating the ppb errors according to Taylor's framwork 
# ---------------------------------------------------------

# change this part
# set the name of batch request depending on the instrument name
#SBATCH -J parallel_mb
# set the directory depending on the date
#SBATCH -D /dss/dsstumfs01/pn69ki/pn69ki-dss-0000/STILT_Ham/HAM_large/HAM_20180905/by-id
# set the output logs depending on the date
#SBATCH -o /dss/dsstumfs01/pn69ki/pn69ki-dss-0000/STILT_Ham/HAM_large/HAM_20180905/by-id/logs/myjob.%j.%N.out
# change this part

#SBATCH --export=NONE
#SBATCH --get-user-env	
# use cm2_tiny and only 1 node
#SBATCH --clusters=cm2_tiny 
#SBATCH --partition=cm2_tiny
#SBATCH --nodes=1
#SBATCH --mail-type=end
#SBATCH --mail-user=aydin.uzun@tum.de
#SBATCH --time=03:00:00

# load correct dependencies to compile R libraries
module unload intel-mpi 
module unload intel
module load gcc/11.2.0
# tools for required packages
# The required packages are raster, geosphere, MASS, ncdf4
module load geos/3.9.1-gcc11
module load proj/8.1.0-gcc11
module load sqlite/3.36.0-gcc11
module load netcdf-hdf5-all/4.7_hdf5-1.10-intel21-serial
module load gdal/3.3.3
module load pkgconf/1.8.0
module load slurm_setup
# load the R module
module load r/4.1.2-gcc11-mkl

# R script to run - change this part depending on the instrument name
Rscript parallel_test_mb.R
```



Again we need to create four shell scripts in the working directory, each belonging to a site. E.g. [job_parallel_test_mb.sh](./Scripts/job_parallel_test_mb.sh) is the one related to the instrument "mb". You need to copy paste this script and create three more named `job_parallel_test_mc.sh` , `job_parallel_test_md.sh`, `job_parallel_test_me.sh`. You can set the name of batch request depending on the instrument name. You need to set the directory and the output logs depending on the desired day. All you need to do is changing `HAM_2018090` . The tricky part here was to determine the tools for required packages. For our R script [job_parallel_test_mb.sh](./Scripts/job_parallel_test_mb.sh) we need the packages `raster`, `geosphere`, `MASS` and `ncdf4`. The necessary modules to work properly with these packages can be seen above. Needless to say that you have to change the part `parallel_test_mb.R` for different sites.

Having the shell scripts ready in the working directory `/dss/dsstumfs01/pn69ki/pn69ki-dss-0000/STILT_Ham/HAM_large/HAM_XXXXXXXX/by-id` we are ready to interact with SLURM:
```bash
# submit a job
sbatch --clusters cm2_tiny job_parallel_test_mb.sh
# check status
squeue --clusters cm2_tiny
# cancel job
scancel --clusters cm2_tiny <job-id>
```

### How to run? - [combine_ppb_error_histograms.R](./Scripts/combine_ppb_error_histograms.R)

After running the previous part for the four different locations for the selected day (thus creating the necessary inputs for this part), we need to run the script [combine_ppb_error_histograms.R](./Scripts/combine_ppb_error_histograms.R) to get the daily ppb error histogram. The working directory is again `/dss/dsstumfs01/pn69ki/pn69ki-dss-0000/STILT_Ham/HAM_large/HAM_XXXXXXXX/by-id` (XXXXXXXX is the day, e.g. 20180905). Depending on the desired day, The parameter `Date` must be selected. Another required input is the scaling factor for the desired day, which can be found at: [Scaling factors](./Data/Scaling_factors/).

The corresponding scaling factor document named `HAM_XXXXXXXX_scaling.csv` must be included in the working directory.

Since this operation is not computationally intensive, no parallelization is required. You can just type on the terminal:

```shell
# 1. Load correct dependencies to compile R libraries
module unload intel-mpi
module unload intel
module load gcc
# 2. Load your desired R module 
module load r
# 3. Run your script
Rscript combine_ppb_error_histograms.R
```
The outputs i.e. daily histograms for the ppb errors when using the 3 different inventories are created and can be found for the selected days at: [Daily_histograms](./Data/Daily_histograms/).

