% ---------------------------------------------------------
% TUM - Technichal University of Munich
%
% Authors:  Aydin Uzun
% Date: 2022
% Purpose: read ERA5 hourly data on pressure levels
% ---------------------------------------------------------
% Temporal resolution : hourly
% Gridded data - Regular latitude-longitude grid
% Vertical coverage : 1000 hPa to 1 hPa
% Vertical resolution : 37 pressure levels
%%
%addpath('/Volumes/esm/projects/WRFDA_Hamburg/ERA5');
cd '/Volumes/esm/20-Research-Projects/04-Modeling/WRF/WRFDA_Hamburg/ERA5';

clear;
close all;

files = dir('*.nc');
length_files = length(files);

%% Initialize matrices
% 24 corresponds to the hours
% 37 is the size of levels
% length_files is 44
time_mat = int32(zeros([length_files,24]));
z_mat = double(zeros([length_files,37,24]));
u_mat = double(zeros([length_files,37,24]));
v_mat = double(zeros([length_files,37,24]));
level_mat = int32(zeros([length_files,37]));


for i = 1:length_files

%% Read data
name=files(i).name;
longitude = ncread(name,'longitude');
% Size:       165x1
% Dimensions: longitude
% Datatype:   single
% Attributes:
% units     = 'degrees_east'
% long_name = 'longitude'
latitude = ncread(name,'latitude');
% Size:       173x1
% Dimensions: latitude
% Datatype:   single
% Attributes:
% units     = 'degrees_north'
% long_name = 'latitude'
level = ncread(name,'level');
% Size:       37x1
% Dimensions: level
% Datatype:   int32
% Attributes:
% units     = 'millibars'
% long_name = 'pressure_level'
time = ncread(name,'time');
% Size:       24x1
% Dimensions: time
% Datatype:   int32
% Attributes:
% units     = 'hours since 1900-01-01 00:00:00.0'
% long_name = 'time'
% calendar  = 'gregorian'
z = ncread(name,'z');
% Size:       165x173x37x24
% Dimensions: longitude,latitude,level,time
% Datatype:   int16
% Attributes:
% scale_factor  = 7.2794
% add_offset    = 238320.1669
% _FillValue    = -32767
% missing_value = -32767
% units         = 'm**2 s**-2'
% long_name     = 'Geopotential'
% standard_name = 'geopotential'
u = ncread(name,'u');
% Size:       165x173x37x24
% Dimensions: longitude,latitude,level,time
% Datatype:   int16
% Attributes:
% scale_factor  = 0.0014301
% add_offset    = 2.5814
% _FillValue    = -32767
% missing_value = -32767
% units         = 'm s**-1'
% long_name     = 'U component of wind'
% standard_name = 'eastward_wind'
v = ncread(name,'v');
% Size:       165x173x37x24
% Dimensions: longitude,latitude,level,time
% Datatype:   int16
% Attributes:
% scale_factor  = 0.0011714
% add_offset    = -0.27396
% _FillValue    = -32767
% missing_value = -32767
% units         = 'm s**-1'
% long_name     = 'V component of wind'
% standard_name = 'northward_wind'

%% Find the grid index corresponding to the lidar location
% lidar location
latitude_windLidar   = 53.5191;
longitude_windLidar   = 10.1029;

% find the latitude and longitude index corresponding to the grid which
% includes the lidar instrument
dist_lat = abs(latitude-latitude_windLidar);
min_dist_lat = min(dist_lat);
lidar_lat_idx=find(dist_lat==min_dist_lat);
dist_lon = abs(longitude-longitude_windLidar);
min_dist_lon = min(dist_lon);
lidar_long_idx=find(dist_lon==min_dist_lon);

%% Handle the problematic data
% It is noticed that if i=12, namely for 'ERA5-20210809-pl.nc'
% time is of size 28x1, 
% z is of size 165x173x37x28
% u is of size 165x173x37x28
% v is of size 165x173x37x28
% if i=13, namely for 'ERA5-20210810-pl.nc'
% time is of size 20x1, 
% z is of size 165x173x37x20
% u is of size 165x173x37x20
% v is of size 165x173x37x20
% which means that 'ERA5-20210809-pl.nc' [i=12] includes also the data for the
% first 4 hours of 'ERA5-20210810-pl.nc' [i=13]
% Here I concatanate the the last 4 hours data to the 20 hours data of 'ERA5-20210810-pl.nc'
if i ==12
suffix = time(25:end);
time = time(1:24);
z_suffix = z(:,:,:,25:28);
u_suffix = z(:,:,:,25:28);
v_suffix = z(:,:,:,25:28);
end
if i ==13
time = [suffix;time];
z = cat(4,z_suffix,z);
u = cat(4,u_suffix,u);
v = cat(4,v_suffix,v);
end

%% Assign
level_mat(i,:) = level';
time_mat(i,:) = time';
length_level = length(level);
length_time = length(time);
z_mat(i,1:length_level,1:length_time) = squeeze(z(lidar_long_idx,lidar_lat_idx,:,1:length_time));
u_mat(i,1:length_level,1:length_time) = squeeze(u(lidar_long_idx,lidar_lat_idx,:,1:length_time));
v_mat(i,1:length_level,1:length_time) = squeeze(v(lidar_long_idx,lidar_lat_idx,:,1:length_time));

end
%% Save ER 
cd '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/ERA5-pl';
% Actually levels are the same for each day!
save(strcat('ERA5','_level.mat'),'level_mat'); %Size [44x37x24] , 44 = files (days), 37 = levels
save(strcat('ERA5','_time.mat'),'time_mat'); %Size [44x24] , 44 = files (days), 24 = hours
save(strcat('ERA5','_z.mat'),'z_mat'); %Size [44x37x24] , 44 = files (days), 37 = levels, 24 = hours
save(strcat('ERA5','_u.mat'),'u_mat'); %Size [44x37x24] , 44 = files (days), 37 = levels, 24 = hours
save(strcat('ERA5','_v.mat'),'v_mat'); %Size [44x37x24] , 44 = files (days), 37 = levels, 24 = hours
cd '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Scripts';

