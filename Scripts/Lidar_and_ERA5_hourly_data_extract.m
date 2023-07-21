cd '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/Lidar_hourly';
clear ;
close all;

load('Lidar_WSPD_hourly.mat');
load('Lidar_WDIR_hourly.mat');
load('Lidar_height_hourly.mat');
load('datetime_hourly.mat');

height_hourly = height_day_by_day_hourly(4:end-1,:,1);
WDIR_hourly = WDIR_day_by_day_hourly(4:end-1,:,1);
WSPD_hourly = WSPD_day_by_day_hourly(4:end-1,:,1);
hours_mat = reshaped_datetime_era5(4:end-1,:);
hours = reshape(hours_mat',[],1);
height_hours = reshape(height_hourly',[],1);
wspd_hours = reshape(WSPD_hourly',[],1);
wdir_hours = reshape(WDIR_hourly',[],1);

T_lidar = table(hours,height_hours,wspd_hours,wdir_hours, 'VariableNames', {'Datetime', 'altitudes','wspd','wdir'});

writetable(T_lidar, 'Lidar_timeseries.csv');

%%
cd '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/ERA5-integrated';

load('ERA5_GPH_ext.mat'); % sorted_GPH_ext 44x24x39
load('ERA5_U_ext.mat'); % sorted_U_ext 44x24x39
load('ERA5_V_ext.mat'); % sorted_V_ext 44x24x39

U = sorted_U_ext(4:end-1,:,1);
V = sorted_V_ext(4:end-1,:,1);
GPH = sorted_GPH_ext(4:end-1,:,1);

%% ERA5 RAW DATA - wind speed wind direction all data points - 39 heights available
% initialize 
WSPD_era5_not_processed = zeros(40,24);
degree_mat_not_processed = zeros(40,24);
% we'll calculate wind speed and wind direction using only the provided raw
% era5 data
% need this for visualization purposes
for i =1:40
    for j=1:24
        WSPD_era5_not_processed(i,j) = sqrt(squeeze(U(i,j)).^2 + squeeze(V(i,j)).^2 );
        degree_mat_not_processed(i,j)=atan2(squeeze(-1*U(i,j)),squeeze(-1*V(i,j)))*(180/pi);
    end
end
degree_mat_not_processed = wrapTo360(degree_mat_not_processed);

%% only consider data from 06:00 [idx=7] until 18:00[idx=19] for each day

geopotential_heights = reshape(GPH',[],1);
WSPD_era5 = reshape(WSPD_era5_not_processed',[],1);
degree_era5 = reshape(degree_mat_not_processed',[],1);

T_era5 = table(hours,geopotential_heights,WSPD_era5,degree_era5, 'VariableNames', {'Datetime', 'geopotential heights','wspd','wdir'});
writetable(T_era5, 'ERA5_timeseries.csv');




