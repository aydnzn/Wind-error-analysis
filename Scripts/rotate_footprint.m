% ---------------------------------------------------------
% TUM - Technichal University of Munich
%
% Authors:  Aydin Uzun
% Date: 2022
% Purpose:  Rotate the aggregated footprints using the daily weighted WDIR difference
% ---------------------------------------------------------
clear;
close all;
clc;
%% allocation
% for the functions
addpath("/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Scripts");
% footprint directory
cd '/Volumes/esm/campaigns/2021AugHamburg/Data/Footprints/foot';

% input footprints
inputFiles = dir( fullfile('*bw14.nc') );
fileNames = { inputFiles.name };

% output folder
outputFolder= "/Volumes/esm/campaigns/2021AugHamburg/Data/Footprints/foot_rotated";

% copy the .nc document and recall it with *_rotated.nc and paste it to the
% output folder
for k = 1 : length(inputFiles )
    thisFileName = fileNames{k};
    % Prepare the input filename.
    inputFullFileName = fullfile(pwd, thisFileName);
    % Prepare the output filename.
    outputBaseFileName = sprintf('%s_rotated.nc', thisFileName(1:end-3));
    outputFullFileName = fullfile(outputFolder, outputBaseFileName);
    % Do the copying and renaming all at once.
    copyfile(inputFullFileName, outputFullFileName);
end

% addpath for the daily weighted WDIR differences
addpath('/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/Weighted_daily_WSPD_WDIR_differences');
T=readtable('weighted_daily_wspd_wdir_differences.csv');
datetable = T.date_set;
% read the wdir differences as rotation angles
rot_angles = T.mean_WDIR_diff;
% there are no footprints for 38:40
rot_angles = rot_angles(1:37);

% output footprints i.e. rotated footprints
outputFiles =  dir( fullfile(outputFolder,'*bw14_rotated.nc') );
output_fileNames = { outputFiles.name };
%% do the rotation
% change directory to output location
cd '/Volumes/esm/campaigns/2021AugHamburg/Data/Footprints/foot_rotated';
for j=1: length(inputFiles )

    % read the rotation angle for the corresponding day
    rot_angle_for_day = rot_angles(j);

    % read the footprint names
    original_fp=fileNames{j};
    rotated_fp=output_fileNames{j};

    % load the latitudes and longitudes
    lon = ncread(original_fp,'lon'); % 140x1 double 9.2050 to 10.5950 -> increasing
    lat = ncread(original_fp,'lat'); % 70x1 double 53.1650 to 53.8550 -> increasing
    recep_time = ncread(original_fp,'recep_time'); % 46x1 int32

    % load the footprints that are going to be rotated
    mb_foot = ncread(original_fp,'mb foot'); % 140x70%46
    mc_foot = ncread(original_fp,'mc foot'); % 140x70%46
    md_foot = ncread(original_fp,'md foot'); % 140x70%46
    me_foot = ncread(original_fp,'me foot'); % 140x70%46

    % initialize rotated footprints
    mb_foot_rotated = zeros(size(mb_foot));
    mc_foot_rotated = zeros(size(mc_foot));
    md_foot_rotated = zeros(size(md_foot));
    me_foot_rotated = zeros(size(me_foot));

    % latitudes and longitudes for the instruments "mb" "mc" "md" "me"
    % These locations are the points around which the ration will occur.
    latitude_instrument =[  53.495,53.536,53.421,53.568 ];
    longitude_instrument =[10.200,9.677,9.892,9.974 ]  ;
    % LON = 70x140 double, LAT = 70x140 double -> create the grid
    % LON(1,1) = 9.2050, LON(end,end) = 10.5950
    % LAT(1,1) = 53.8550 LAT(end,end) = 53.1650
    % now we have the right grid for the NORTH EAST hemisphere
    [LON, LAT]=meshgrid(lon,flip(lat));
    % find the latitude and longitude index on the mesh
    % e.g. original latitude mb = 53.495,  original longitude mb = 10.20
    % LAT(latitude_mb,longitude_mb) = 53.4950
    % LON(latitude_mb,longitude_mb) = 10.1950
    [latitude_mb,longitude_mb] =find_instrument_index_in_mesh(LON,LAT,longitude_instrument(1),latitude_instrument(1));
    [latitude_mc,longitude_mc] =find_instrument_index_in_mesh(LON,LAT,longitude_instrument(2),latitude_instrument(2));
    [latitude_md,longitude_md] =find_instrument_index_in_mesh(LON,LAT,longitude_instrument(3),latitude_instrument(3));
    [latitude_me,longitude_me] =find_instrument_index_in_mesh(LON,LAT,longitude_instrument(4),latitude_instrument(4));

    % receptor time length
    time_length = length(recep_time);

    for i = 1:time_length
        % as squeeze(mb_foot(:,:,i)) 140x70 double
        % so now lat is on the x axis increasing to the right, lon is on
        % the y axis increasing to the down
        % if I rotate this 90 deg ccw
        % then lon is on the x axis increasing to the right, lat is on the
        % y axis increasing to the up
        % now we can imagine the footprint as an image on the north east
        % hemisphere
        mb_foot_image_i = rot90(squeeze(mb_foot(:,:,i)));
        mc_foot_image_i = rot90(squeeze(mc_foot(:,:,i)));
        md_foot_image_i = rot90(squeeze(md_foot(:,:,i)));
        me_foot_image_i = rot90(squeeze(me_foot(:,:,i)));

        % rotation
        mb_foot_image_rotated_i =rotate_around(rot_angle_for_day,mb_foot_image_i,latitude_mb,longitude_mb);
        mc_foot_image_rotated_i =rotate_around(rot_angle_for_day,mc_foot_image_i,latitude_mc,longitude_mc);
        md_foot_image_rotated_i =rotate_around(rot_angle_for_day,md_foot_image_i,latitude_md,longitude_md);
        me_foot_image_rotated_i =rotate_around(rot_angle_for_day,me_foot_image_i,latitude_me,longitude_me);

        % rotate 90 deg cw, to get back to the original format
        mb_foot_rotated_i = rot90(mb_foot_image_rotated_i,-1);
        mc_foot_rotated_i = rot90(mc_foot_image_rotated_i,-1);
        md_foot_rotated_i = rot90(md_foot_image_rotated_i,-1);
        me_foot_rotated_i = rot90(me_foot_image_rotated_i,-1);

        % save the footprint at time i
        mb_foot_rotated(:,:,i) = mb_foot_rotated_i;
        mc_foot_rotated(:,:,i) = mc_foot_rotated_i;
        md_foot_rotated(:,:,i) = md_foot_rotated_i;
        me_foot_rotated(:,:,i) = me_foot_rotated_i;
    end

    % now change the footprint variable to the new rotated one on the
    % rotated footprint nc document
    ncwrite(rotated_fp,'mb foot',single(mb_foot_rotated));
    ncwrite(rotated_fp,'mc foot',single(mc_foot_rotated));
    ncwrite(rotated_fp,'md foot',single(md_foot_rotated));
    ncwrite(rotated_fp,'me foot',single(me_foot_rotated));

end
