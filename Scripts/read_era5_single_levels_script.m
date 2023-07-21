% ---------------------------------------------------------
% TUM - Technichal University of Munich
%
% Authors:  Aydin Uzun
% Date: 2022
% Purpose: read ERA5 hourly data on single levels
% ---------------------------------------------------------
% Temporal resolution : hourly
% Gridded data - Regular latitude-longitude grid
% Variables of interest:
% 100m u-component of wind
% 100m v-component of wind
% 10m u-component of wind
% 10m v-component of wind
% Geopotential : The (surface) geopotential height can be calculated by dividing the (surface) geopotential by the Earth's gravitational acceleration.
cd '/Volumes/esm/20-Research-Projects/04-Modeling/WRF/WRFDA_Hamburg/ERA5/ERA5-sl';

clear;
close all;

% names of the files 
Jul_name = 'ERA5-2020729-20210731-sl.nc';
Aug_name = 'ERA5-20210801-20210831-sl.nc';
Aug_name_blh = 'ERA5-20210801-20210831-sl_blh.nc';
Sep_name = 'ERA5-2021901-20210910-sl.nc';
Sep_name_blh = 'ERA5-20210901-20210910-sl_blh.nc';

% ncdisp('Aug_name');
 
%% Read data
% longitude Size: 13x1
% latitude Size: 21x1
% time 'hours since 1900-01-01 00:00:00.0'
% u100 longitude,latitude,time 100 metre U wind component
% v100 longitude,latitude,time 100 metre V wind component
% u10 longitude,latitude,time 10 metre U wind component
% v10 longitude,latitude,time 10 metre V wind component
% z longitude,latitude,time Geopotential units = 'm**2 s**-2'

longitude_07 = ncread(Jul_name,'longitude');
latitude_07 = ncread(Jul_name,'latitude');
time_07 = ncread(Jul_name,'time');
z_07 = ncread(Jul_name,'z');
u100_07 = ncread(Jul_name,'u100');
v100_07 = ncread(Jul_name,'v100');
u10_07 = ncread(Jul_name,'u10');
v10_07 = ncread(Jul_name,'v10');

longitude_08 = ncread(Aug_name,'longitude');
latitude_08 = ncread(Aug_name,'latitude');
time_08 = ncread(Aug_name,'time');
z_08 = ncread(Aug_name,'z');
u100_08 = ncread(Aug_name,'u100');
v100_08 = ncread(Aug_name,'v100');
u10_08 = ncread(Aug_name,'u10');
v10_08 = ncread(Aug_name,'v10');


longitude_09 = ncread(Sep_name,'longitude');
latitude_09 = ncread(Sep_name,'latitude');
time_09 = ncread(Sep_name,'time');
z_09 = ncread(Sep_name,'z');
u100_09 = ncread(Sep_name,'u100');
v100_09 = ncread(Sep_name,'v100');
u10_09 = ncread(Sep_name,'u10');
v10_09 = ncread(Sep_name,'v10');

longitude_08_blh = ncread(Aug_name_blh,'longitude');
latitude_08_blh = ncread(Aug_name_blh,'latitude');
time_08_blh = ncread(Aug_name_blh,'time');
z_08_blh = ncread(Aug_name_blh,'z');
u100_08_blh = ncread(Aug_name_blh,'u100');
v100_08_blh = ncread(Aug_name_blh,'v100');
u10_08_blh = ncread(Aug_name_blh,'u10');
v10_08_blh = ncread(Aug_name_blh,'v10');


longitude_09_blh = ncread(Sep_name_blh,'longitude');
latitude_09_blh = ncread(Sep_name_blh,'latitude');
time_09_blh = ncread(Sep_name_blh,'time');
z_09_blh = ncread(Sep_name_blh,'z');
u100_09_blh = ncread(Sep_name_blh,'u100');
v100_09_blh = ncread(Sep_name_blh,'v100');
u10_09_blh = ncread(Sep_name_blh,'u10');
v10_09_blh = ncread(Sep_name_blh,'v10');

%% Find the grid index corresponding to the lidar location
% lidar location
latitude_windLidar   = 53.5191;
longitude_windLidar   = 10.1029;

% find the latitude and longitude index corresponding to the grid which
% includes the lidar instrument
dist_lat = abs(latitude_07-latitude_windLidar);
min_dist_lat = min(dist_lat);
lidar_lat_idx=find(dist_lat==min_dist_lat);
dist_lon = abs(longitude_07-longitude_windLidar);
min_dist_lon = min(dist_lon);
lidar_long_idx=find(dist_lon==min_dist_lon);
%% Create the matrices
z_07_lidar_loc = squeeze(z_07(lidar_long_idx,lidar_lat_idx,:));
z_08_lidar_loc = squeeze(z_08(lidar_long_idx,lidar_lat_idx,:));
z_09_lidar_loc = squeeze(z_09(lidar_long_idx,lidar_lat_idx,1:240)); % only consider the first 10 days, 10x24 hours = 240
z_sl = [z_07_lidar_loc;z_08_lidar_loc;z_09_lidar_loc]; % concatanated
z_mat_sl = reshape(z_sl,24,44); % reshape into matrix form, 24 hours 44 days
z_mat_sl = z_mat_sl'; 
g =9.80665;
gph_mat_sl = z_mat_sl/g; % calculate height from geopotential

% now do the same matrix creation operations for u and v 
u10_07_lidar_loc = squeeze(u10_07(lidar_long_idx,lidar_lat_idx,:));
u10_08_lidar_loc = squeeze(u10_08(lidar_long_idx,lidar_lat_idx,:));
u10_09_lidar_loc = squeeze(u10_09(lidar_long_idx,lidar_lat_idx,1:240));
u10_sl = [u10_07_lidar_loc;u10_08_lidar_loc;u10_09_lidar_loc];
u10_mat_sl = reshape(u10_sl,24,44);
u10_mat_sl = u10_mat_sl';

u100_07_lidar_loc = squeeze(u100_07(lidar_long_idx,lidar_lat_idx,:));
u100_08_lidar_loc = squeeze(u100_08(lidar_long_idx,lidar_lat_idx,:));
u100_09_lidar_loc = squeeze(u100_09(lidar_long_idx,lidar_lat_idx,1:240));
u100_sl = [u100_07_lidar_loc;u100_08_lidar_loc;u100_09_lidar_loc];
u100_mat_sl = reshape(u100_sl,24,44);
u100_mat_sl = u100_mat_sl';

v10_07_lidar_loc = squeeze(v10_07(lidar_long_idx,lidar_lat_idx,:));
v10_08_lidar_loc = squeeze(v10_08(lidar_long_idx,lidar_lat_idx,:));
v10_09_lidar_loc = squeeze(v10_09(lidar_long_idx,lidar_lat_idx,1:240));
v10_sl = [v10_07_lidar_loc;v10_08_lidar_loc;v10_09_lidar_loc];
v10_mat_sl = reshape(v10_sl,24,44);
v10_mat_sl = v10_mat_sl';

v100_07_lidar_loc = squeeze(v100_07(lidar_long_idx,lidar_lat_idx,:));
v100_08_lidar_loc = squeeze(v100_08(lidar_long_idx,lidar_lat_idx,:));
v100_09_lidar_loc = squeeze(v100_09(lidar_long_idx,lidar_lat_idx,1:240));
v100_sl = [v100_07_lidar_loc;v100_08_lidar_loc;v100_09_lidar_loc];
v100_mat_sl = reshape(v100_sl,24,44);
v100_mat_sl = v100_mat_sl';

% to check the corresponding times 
time_sl = [time_07;time_08;time_09(1:240)];
time_mat_sl = reshape(time_sl,24,44);
time_mat_sl = time_mat_sl';
time_double = double(time_mat_sl)./24 + datenum('1900-01-01 00:00:00');
time_vec = datevec(time_double');
time_datetime = datetime(time_vec);

%% Save
cd '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/ERA5-sl';
save('v100.mat','v100_mat_sl'); % [44x24] , 44 days x 24 hours
save('u100.mat','u100_mat_sl'); % [44x24] , 44 days x 24 hours
save('v10.mat','v10_mat_sl'); % [44x24] , 44 days x 24 hours
save('u10.mat','u10_mat_sl'); % [44x24] , 44 days x 24 hours
save('gph_sl.mat','gph_mat_sl'); % [44x24] , 44 days x 24 hours -> also created a matrix for this but we expect this to be scalar.
cd '/Volumes/esm/11-Thesis/03-Scientific Internship/2021 FP Aydin Uzun/Scripts';

