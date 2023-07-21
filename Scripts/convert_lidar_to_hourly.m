% ---------------------------------------------------------
% TUM - Technichal University of Munich
%
% Authors:  Aydin Uzun
% Date: 2022
% Purpose: convert lidar data to hourly data
% ---------------------------------------------------------
clear;
close all;
%%
addpath('/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/Lidar');
% load lidar measurements
load('WSPD.mat');
load('WDIR.mat');
load('UZ.mat');
load('height.mat');
load('datestrings.mat');

addpath('/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/ERA5-pl');
% load ERA5 times 44x24 int32
load('ERA5_time.mat');

%%
% make a column
time_mat=time_mat';
time_mat = reshape(time_mat,1056,1);

% now convert to 1056x1 datetime
time_era5_double = double(time_mat)./24 + datenum('1900-01-01 00:00:00');
time_era5_vec = datevec(time_era5_double);
datetime_era5 = datetime(time_era5_vec); % hourly

% initialize hourly data
WDIR_hourly = zeros(1056,149);
WSPD_hourly = zeros(1056,149);
height_hourly = zeros(1056,149);

% lidar datetime
datetime_lidar = datetime(datestrings,'InputFormat','yyyyMMddHHmmSS');
% convert to datenumber
DateNumber_era5 = datenum(datetime_era5);
DateNumber_lidar = datenum(datetime_lidar);

%% cycle through 44 [days]x24 [hours]=1056 hours
for i=1:1056

    % find the closest timepoint index to the sharp hour
    diff_to_datenumber_era5 =DateNumber_era5(i) - DateNumber_lidar ;
    min_diff_to_datenumber_era5 =min(abs(round(diff_to_datenumber_era5,5)));
    closest_idx = find(abs(round(diff_to_datenumber_era5,5)) ==min_diff_to_datenumber_era5);
    length_closest_idx =length(closest_idx);


    % if there is only 1 closest index to the sharp hour then no further
    % processing
    if length_closest_idx==1
        WSPD_hourly(i,:) = WSPD_mat(closest_idx,:);
        WDIR_hourly(i,:) = WDIR_mat(closest_idx,:);
        height_hourly(i,:) = height_mat(closest_idx,:);
    end

    % if there are 2 closest index to the sharp hour, (e.g. +5 min and -5
    % min)
    if length_closest_idx>1
        % using unit vectors
        % calculate the mean of sin components - sin_comp
        % calculate the mean of cos components - cos_comp
        sin_comp =  mean( sin(WDIR_mat(closest_idx,:)*pi/180),'omitnan');
        cos_comp =  mean( cos(WDIR_mat(closest_idx,:)*pi/180),'omitnan');
        % then find the angle and assign
        wdir=atan2(sin_comp,cos_comp)*180/pi;
        wrap_wdir = wrapTo360(wdir);
        WDIR_hourly(i,:) = wrap_wdir;
        % mean of 2 heights
        height_hourly(i,:) = mean(height_mat(closest_idx,:),'omitnan');
        % mean of wspds
        WSPD_hourly(i,:) = mean(WSPD_mat(closest_idx,:),'omitnan');
    end
end

%% Inıtialize 44 = days, 24 = hours, 149 = altitudes
WSPD_day_by_day_hourly = zeros(44,24,149);
WDIR_day_by_day_hourly = zeros(44,24,149);
height_day_by_day_hourly = zeros(44,24,149);

%% Convert to 3D Mat
% WSPD_hourly 1056 x 149
% WDIR_hourly 1056 x 149
% height_hourly 1056 x 149
for i = 1:44
    for j=1:24
        idx = (24*(i-1))+j;
        WSPD_day_by_day_hourly(i,j,:) = WSPD_hourly(idx,:);
        WDIR_day_by_day_hourly(i,j,:) = WDIR_hourly(idx,:);
        height_day_by_day_hourly(i,j,:) = height_hourly(idx,:);

    end
end
%% Convert datetime to 2D
% datetime_era5 1056x1
reshaped_datetime_era5 = reshape(datetime_era5,[24,44]);
reshaped_datetime_era5 = reshaped_datetime_era5';
%%
cd '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/Lidar_hourly';
save('Lidar_WSPD_hourly.mat','WSPD_day_by_day_hourly');
save('Lidar_WDIR_hourly.mat','WDIR_day_by_day_hourly');
save('datetime_hourly.mat','reshaped_datetime_era5');
save('Lidar_height_hourly.mat','height_day_by_day_hourly');
cd '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Scripts';

