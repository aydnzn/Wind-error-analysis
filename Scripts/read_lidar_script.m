% ---------------------------------------------------------
% TUM - Technichal University of Munich
%
% Authors:  Aydin Uzun
% Date: 2022
% Purpose:  Read all lidar measurements and save them
% ---------------------------------------------------------
clear;
close all;
%%
% old location
% cd '/Volumes/esm/projects/WRFDA_Hamburg/202108_lidar_wettermast_hamburg/DLR85';
% new location
% cd '/Volumes/esm/20-Research-Projects/02-Measurement-Campaigns/2021AugHamburg/Data/Wind Measurements/LIDAR Measurements';
cd '/Volumes/esm/20-Research-Projects/02-Measurement-Campaigns/2021AugHamburg/Data/Wind Measurements/LIDAR Measurements/DLR85';

files = dir('*.nc');
%%
% Initialize the outputs
height_mat = zeros(length(files),149);
WSPD_mat = zeros(length(files),149);
WDIR_mat = zeros(length(files),149);
UZ_mat = zeros(length(files),149);
CNR_mat = zeros(length(files),149);
datestrings = cell(length(files),1);

for i =1:length(files) % 7170 measurements
   name=files(i).name;
   
   height = ncread(name,'Z'); % in m
   WSPD = ncread(name,'WSPD'); % hor. wind speed m/s
   WDIR = ncread(name,'WDIR'); % hor. wind direction deg
   UZ = ncread(name,'UZ'); % vert. wind speed m/s
   CNR = ncread(name,'CNR'); % carrier to noise ratio dB
   
   height_mat(i,1:length(height)) = height;
   WSPD_mat(i,1:length(WSPD)) = WSPD;
   WDIR_mat(i,1:length(WDIR)) = WDIR;
   UZ_mat(i,1:length(UZ)) = UZ;
   CNR_mat(i,1:length(CNR)) = CNR;

date = strsplit(name,'_');
date2 = strsplit(date{1,6},'.');
date3 = date2{1,1};
datestrings{i,1} = date3;
end
%%
cd '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/Lidar';
save('height.mat','height_mat'); % 7170x149 double
save('datestrings.mat','datestrings'); % 7170x1 cell
save('WSPD.mat','WSPD_mat'); % 7170x149 double
save('UZ.mat','UZ_mat'); % 7170x149 double
save('CNR.mat','CNR_mat'); % 7170x149 double
save('WDIR.mat','WDIR_mat'); % 7170x149 double
cd '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Scripts';

