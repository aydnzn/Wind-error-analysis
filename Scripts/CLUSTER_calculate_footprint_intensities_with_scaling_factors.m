% ---------------------------------------------------------
% TUM - Technichal University of Munich
%
% Authors:  Aydin Uzun
% Date: 2022
% Purpose: calculate footprint intensities 
% Note: Run this script on /dss/dsstumfs01/pn69ki/pn69ki-dss-0000/STILT_Ham/HAM_large
% pwd should be /dss/dsstumfs01/pn69ki/pn69ki-dss-0000/STILT_Ham/HAM_large
% then the output document 'HAM_large_foot_intensity_with_scaling.csv'
% is downloaded manually to NAS
% ---------------------------------------------------------
clear;
close all;
clc;
%% Load normalized scaling factors and addpath the pwd
addpath(genpath(pwd));
load('normalized_scaling.mat');

% list footprints in the directories
FileList = dir(fullfile(pwd, '**', '*foot.nc'));
%% Initialize
instrument_info = cell(length(FileList),1);
date_info = cell(length(FileList),1);
altitude_info = zeros(length(FileList),1);
name_info = cell(length(FileList),1);
lon = zeros(140,1);
lat = zeros(70,1);
time_matrix = zeros(length(FileList),1);
foot_matrix = zeros(length(FileList),140,70);

%% ASSIGN the available data in the directories to matrices
for i = 1: length(FileList)
    % name e.g. mb_202108010600_186m_foot.nc
    name_i = FileList(i).name;
    % split by _ e.g. 4x1 cell "mb" "202108010600" "186m" "foot.nc"
    split_name_i = split(name_i,'_');
    % split the altitude e.g. split 186m to 186
    split_name_i_altitude = split(split_name_i(3),'m');
    split_name_i_altitude = split_name_i_altitude(1);
    % then convert the altitude to double from string and save
    altitude_info(i) = str2double(split_name_i_altitude);
    % save instrument info e.g. "mb"
    instrument_info(i) = split_name_i(1);
    % save date info e.g. "202108010600"
    date_info(i) =  split_name_i(2);
    % save full name info "mb_202108010600_186m_foot.nc"
    name_info(i) = cellstr(name_i);
    
    % no need to save lat and lon again and again
    % as they are the same for all *.nc documents
    % lat is 70x1 double from 53.1650 to 53.8550 increase by 0.01
    % lon is 140x1 double from 9.2050 to 10.5950 increase by 0.01
    if i==1
        lon = ncread(name_info{i},'lon');
        lat = ncread(name_info{i},'lat');
    end
    % read time information directly from the nc document
    % this can be alternatively done using the name of the nc document
    time_matrix(i,:) = ncread(name_info{i},'time');
    % save the footprint
    foot_matrix(i,:,:) = squeeze(ncread(name_info{i},'foot'));
end

% do the time conversion
time_double = double(time_matrix)./86400 + datenum('1970-01-01 00:00:00');
% time_vec of size length(FileList)x6 e.g. 2021 8 1 6 0 0 [year month day hour minute sec]
time_vec = datevec(time_double);
datetime_vec = datetime(time_vec);
% get minutes
minutes = time_vec(:,5);
% find the idx for the FileList for shap hours
idx_for_sharp_hours=find(minutes==0);
%% create a struct for the data at the sharp hour
sharp=struct;
% instruments at the sharp hour idx
sharp.instrument = instrument_info(idx_for_sharp_hours);
% date at the sharp hour idx
sharp.date = date_info(idx_for_sharp_hours);
% altitude at the sharp hour idx
sharp.altitude = altitude_info(idx_for_sharp_hours);
% name at the sharp hour idx
sharp.name = name_info(idx_for_sharp_hours);
% footprint at the sharp hour idx
sharp.foot = foot_matrix(idx_for_sharp_hours,:,:);
% time at the sharp hour idx
sharp.time = datetime_vec(idx_for_sharp_hours,:);

% determine the unique sharp times
% e.g. 
% 01-Aug-2021 06:00:00
% 01-Aug-2021 07:00:00 ...
unique_times = unique(sharp.time);
% determine the day number e.g. for 01-aug-2021 it's 1, for 02-aug-2021
% it's 2
day_number=ceil(days(sharp.time - '01-aug-2021'));
%% Initialize the output document
% remember unique_times only contain sharp hours
% length(unique_times) will be multiplied by the available instrument
% number which is 4, as at each sharp hour and altitude level we have 4
% different footprints corresponding to mb mc md me
sz = length(unique_times)*4;
% sharp hour as datetime, instrument name as string
varTypes = ["datetime","string","double"];
Datetime = NaT(sz,1);
Instrument = strings([sz,1]);
Intensity_and_Scaling=single(zeros(sz,13));
% Intensity_and_Scaling will be the footprint intensities integrated with
% pressure scaling factors, for each sharp hour and instrument pair we have
% 13 numbers defining the atmosphere whose sum will be 1
My_table = table(Datetime,Instrument,Intensity_and_Scaling);
%% Calculate footprint intensities

for j=1: length(unique_times) % for each sharp hour available 
    
    % find the idx corresponding to the sharp hours
    index_of_the_hour = find(sharp.time == unique_times(j));
    % altitudes and instruments 
    altitude_for_the_hour = sharp.altitude(index_of_the_hour);
    instruments_for_the_hour = sharp.instrument(index_of_the_hour);
    % footprints
    foot_for_the_hour = sharp.foot(index_of_the_hour,:,:);
    % day number e.g. 1 for "01-aug-2021", 2 for "02-aug-2021"
    day_for_unique_time = day_number(index_of_the_hour);

    
    % unique instruments for the chosen sharp hour
    % e.g. 4x1 cell array if all the 4 {mb} {mc} {md} {me} are available
    unique_instruments_for_the_hour=unique(instruments_for_the_hour);
    number_of_inst_for_the_hour = length(unique_instruments_for_the_hour);
    
    % initialize fp intensity vector
    fp_intensity_vector = zeros(13,1);
    for k=1:number_of_inst_for_the_hour % number_of_inst_for_the_hour: how many instruments are available at the sharp hour chosen?
        % k=1 "mb", k=2 "mc", k=3 "md", k=4 "md"
        idx_for_inst = find(strcmp(unique_instruments_for_the_hour{k}, instruments_for_the_hour)); % find idx on instruments_for_the_hour for the chosen instrument name
        
        % altitudes and footprints for the chosen sharp hour and for the
        % chosen instrument 
        foot_for_sharp_hour_at_inst = foot_for_the_hour(idx_for_inst,:,:);
        altitude_for_sharp_hour_at_inst = altitude_for_the_hour(idx_for_inst);
        % sort the altitudes 
        % e.g. sorted_alt = [20,186, 355, 527, 701, 879] if there are no
        % available footprints for the higher altitudes
        [sorted_alt,sorted_alt_idx] = sort(altitude_for_sharp_hour_at_inst);
        % sum of the footprints for the chosen sharp hour and for the
        % chosen instrument 
        % e.g. sum_of_fp_for_the_hour_at_inst =
        % [1.4311;0.6974;0.1640;0.0429;0.0089,0.0016] with first one
        % corresponding to 20m altitude the second one corresponding to
        % 186m altitude etc.
        sum_of_fp_for_the_hour_at_inst = sum(foot_for_sharp_hour_at_inst(sorted_alt_idx,:,:),[2 3]);
        % normalize it
        norm_sum_of_fp_for_the_hour_at_inst = sum_of_fp_for_the_hour_at_inst/sum(sum_of_fp_for_the_hour_at_inst);
        % sum(norm_sum_of_fp_for_the_hour_at_inst)=1
        % then assign fp_intensity_vector was of size 13x1
        % it will be all 0s for the layers where footprints are not
        % available
        fp_intensity_vector(1:length(norm_sum_of_fp_for_the_hour_at_inst))=norm_sum_of_fp_for_the_hour_at_inst;
        
        % normalized_scaling 40x13 double [day_number x layers]
        % day_for_unique_time(1) will be day number
        scaling_vec = normalized_scaling(day_for_unique_time(1),:)'; % 13x1 double
        % sum(scaling_vec)=1
        % do the operation integrate the footprint intensity vector with
        % the pressure scaling vector
        fp_intensity_and_scaling =fp_intensity_vector.*scaling_vec;
        % normalize
        fp_intensity_and_scaling= fp_intensity_and_scaling/sum(fp_intensity_and_scaling);

        % index to fill the output vector
        index =(j-1)*4+k;
        My_table.Datetime(index)=unique_times(j);
        My_table.Instrument(index) = unique_instruments_for_the_hour{k};
        My_table.Intensity_and_Scaling(index,:) = single(fp_intensity_and_scaling);

    end


end
% Datetime,Instrument,Intensity_and_Scaling_1,Intensity_and_Scaling_2...
% Intensity_and_Scaling_3,Intensity_and_Scaling_4,Intensity_and_Scaling_5...
% Intensity_and_Scaling_6,Intensity_and_Scaling_7,Intensity_and_Scaling_8...
% Intensity_and_Scaling_9,Intensity_and_Scaling_10,Intensity_and_Scaling_11...
% Intensity_and_Scaling_12,Intensity_and_Scaling_13
% e.g. a line on the table would be "01-Aug-2021
% 06:00:00,mb,0.6090592,0.2978604,0.07023609,0.01831584,0.003827032,0.000701371,0,0,0,0,0,0,0"
writetable(My_table,'HAM_large_foot_intensity_with_scaling.csv')

