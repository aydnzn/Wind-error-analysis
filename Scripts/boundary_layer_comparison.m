cd '/Volumes/esm/20-Research-Projects/04-Modeling/WRF/WRFDA_Hamburg/ERA5/ERA5-sl';
clear;
close all;

Aug_name_blh = 'ERA5-20210801-20210831-sl_blh.nc';
Sep_name_blh = 'ERA5-20210901-20210910-sl_blh.nc';

longitude_08_blh = ncread(Aug_name_blh,'longitude');
latitude_08_blh = ncread(Aug_name_blh,'latitude');
time_08_blh = ncread(Aug_name_blh,'time');
blh_08_blh = ncread(Aug_name_blh,'blh');
z_08_blh = ncread(Aug_name_blh,'z');


longitude_09_blh = ncread(Sep_name_blh,'longitude');
latitude_09_blh = ncread(Sep_name_blh,'latitude');
time_09_blh = ncread(Sep_name_blh,'time');
blh_09_blh = ncread(Sep_name_blh,'blh');
z_09_blh = ncread(Sep_name_blh,'z');

%% Find the grid index corresponding to the lidar location
% lidar location
latitude_windLidar   = 53.5191;
longitude_windLidar   = 10.1029;

% find the latitude and longitude index corresponding to the grid which
% includes the lidar instrument
dist_lat = abs(latitude_08_blh-latitude_windLidar);
min_dist_lat = min(dist_lat);
lidar_lat_idx=find(dist_lat==min_dist_lat);
dist_lon = abs(longitude_08_blh-longitude_windLidar);
min_dist_lon = min(dist_lon);
lidar_long_idx=find(dist_lon==min_dist_lon);

%%
z_08_lidar_loc = squeeze(z_08_blh(lidar_long_idx,lidar_lat_idx,:));
z_09_lidar_loc = squeeze(z_09_blh(lidar_long_idx,lidar_lat_idx,:));
z_sl = [z_08_lidar_loc;z_09_lidar_loc]; % concatanated
z_mat_sl = reshape(z_sl,24,41); % reshape into matrix form, 24 hours 44 days
z_mat_sl = z_mat_sl'; 
g =9.80665;
gph_mat_sl = z_mat_sl/g; % calculate height from geopotential

% now do the same matrix creation operations for u and v 
blh_08_lidar_loc = squeeze(blh_08_blh(lidar_long_idx,lidar_lat_idx,:));
blh_09_lidar_loc = squeeze(blh_09_blh(lidar_long_idx,lidar_lat_idx,:));
blh_sl = [blh_08_lidar_loc;blh_09_lidar_loc];
blh_mat_sl = reshape(blh_sl,24,41);
blh_mat_sl = blh_mat_sl';



% to check the corresponding times 
time_sl = [time_08_blh;time_09_blh(1:240)];
time_mat_sl = reshape(time_sl,24,41);
time_mat_sl = time_mat_sl';
time_double = double(time_mat_sl)./24 + datenum('1900-01-01 00:00:00');
time_vec = datevec(time_double');
time_datetime = datetime(time_vec);

%%
cd '/Volumes/esm/20-Research-Projects/04-Modeling/WRF/WRFDA_Hamburg/ERA5/ERA5-sl';
file = 'norman_blh_data.txt';

% use the 'importdata' function to read the file
data = readtable(file);
datetime_array = data.Var1;
blh_lidar = data.Var2;

% set the year of each date in the data
datetime_array_new = datetime_array;
for i = 1:length(datetime_array)
    datetime_array_new(i) = setfield(datetime_array_new(i), 'Year', 2021);
end


new_date_time = datestr(time_datetime, 'yyyy-mm-dd HH:MM:SS');
%%
cd '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Boundary Layer Figures';
T = table(time_datetime,blh_sl, 'VariableNames', {'Date', 'BLH'});
T_lidar = table(datetime_array_new,blh_lidar, 'VariableNames', {'Date', 'BLH'});

% Export the table to a CSV file
writetable(T, 'era5_blh.csv');
writetable(T_lidar, 'lidar_blh.csv');

%%
close all;
f= figure;
plot(time_datetime,blh_sl,'Linewidth',1.5); hold on;
plot(datetime_array_new,blh_lidar,'Linewidth',1.5);
legend('ERA5 BLH','LiDAR BLH');
% set the x-axis label
xlabel('Date','Fontsize',14);

% set the y-axis label
ylabel('Height in m','Fontsize',14);
set(gcf, 'Position', [100, 100, 1920, 1080]);

%%%%%%%
xlim([time_datetime(529), time_datetime(864)]);
xticks([datetime('2021-08-23'), datetime('2021-08-24'), datetime('2021-08-31'), datetime('2021-09-01'), datetime('2021-09-03'), datetime('2021-09-05')]);
%%%%%%%%

yticks(250:250:2000)
ylim([0,2000]);
% format the x-axis as dates
datetick('x', 'yyyy-mm-dd');
set(gca, 'Fontsize',14);
%%
cd '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Boundary Layer Figures';
print(gcf, '-dpng', 'era_vs_lidar_2');


