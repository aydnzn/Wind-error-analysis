% ---------------------------------------------------------
% TUM - Technichal University of Munich
%
% Authors:  Aydin Uzun
% Date: 2022
% Purpose: use FP intensities with scaling factors to weigh the wspd and wdir differences
% ---------------------------------------------------------
clear;
close all;
%% load wind speed and wind direction difference matrices
addpath('/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/Lidar_ERA5_representatives');
load('WSPD_differences.mat'); % [lidar - era5 model] % 40x13x13 days x hours x layers
load('WDIR_differences.mat'); % [lidar - era5 model] % 40x13x13 days x hours x layers
load('Datetime_hourly.mat');
addpath('/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/HAM_large_footprint_intensities_with_scaling_factors');
% read intensities
T = readtable('HAM_large_foot_intensity_with_scaling.csv');


%% factors going to be used to weigh
factors=[T.Intensity_and_Scaling_1,T.Intensity_and_Scaling_2,T.Intensity_and_Scaling_3,T.Intensity_and_Scaling_4,T.Intensity_and_Scaling_5,T.Intensity_and_Scaling_6, ...
    T.Intensity_and_Scaling_7,T.Intensity_and_Scaling_8,T.Intensity_and_Scaling_9,T.Intensity_and_Scaling_10,T.Intensity_and_Scaling_11,T.Intensity_and_Scaling_12,...
    T.Intensity_and_Scaling_13];
%% Average of intensities at "mb" "mc" "md" "me" at each sharp hour
datetime_factors = T.Datetime;
% unique_datetime_facs 438x1 datetime
% 01-Aug-2021 06:00:00
% 01-Aug-2021 07:00:00
% 01-Aug-2021 08:00:00 ...
unique_datetime_facs = unique(datetime_factors);
% last 3 elements of unique_datetime_facs is NaT
% because on 13-Aug-2021 06:00:00 only footprints at locations "mc" "md" are
% available we get 2 NaTs for "mb" and "me"
% on 07-Sep-2021 07:00:00 only footprints at locations "mb" "mc" "me" are
% available we get 1 NaT for "md"
% at all other sharp hours footprints are available at 4 instrument
% locations
% ignore NaTs
unique_datetime_facs(isnat(unique_datetime_facs))=[];

% initialize
% facs will be the averaged (mean of fp intensities at locations mb mc md
% me) intensity for each sharp hour
facs = zeros(length(unique_datetime_facs),13);

for i=1:length(unique_datetime_facs)

    % get the index for the sharp hour
    % e.g. for 01-Aug-2021 06:00:00 idx_for_date=[1;2;3;4]
    idx_for_date = find(datetime_factors==unique_datetime_facs(i));
    % average of all mb mc md me
    facs_for_date = mean(factors(idx_for_date,:))/sum(mean(factors(idx_for_date,:)));
    % assign the calculated average to a vector
    facs(i,:) =facs_for_date;
end
%% for easier further processing convert the 3D wspd_diff_vec and wdir_diff_vec to 2D matrices
% vectorized datetime vector , 520x1 datetime
datetime_vec = reshape(datetime_mat',[],1) ;
% initialize
size_diff = size(WSPD_diff);
day_sz= size_diff(1); hour_sz= size_diff(2); layer_sz= size_diff(3);
wspd_diff_vec = zeros(day_sz*hour_sz,layer_sz);
wdir_diff_vec = zeros(day_sz*hour_sz,layer_sz);

for i =1:day_sz % days
    for j=1:hour_sz % hours
        my_idx = (i-1)*layer_sz+j;
        % wspd_diff_vec size: 520x13
        wspd_diff_vec(my_idx,:) = WSPD_diff(i,j,:);
        % wdir_diff_vec size: 520x13
        wdir_diff_vec(my_idx,:) = wrapped_wdir_diff(i,j,:);
    end
end
%% Calculate daily values using the intensities
% day_number 435x1 double 1 for 01-aug-2021 etc.
day_number=ceil(days(unique_datetime_facs - '01-aug-2021'));
wspd_mean_vec = single(zeros(day_sz,1));
wdir_mean_vec = single(zeros(day_sz,1));
wspd_std_vec = single(zeros(day_sz,1));
wdir_std_vec = single(zeros(day_sz,1));

for day_no=1:day_sz % days

    % sharp hour idx for the day
    % e.g. idx_day is of size 12x1
    % and datetime_vec_for_day is 01-aug-2021 06:00:00; 01-aug-2021
    % 07:00:00 ... until 01-aug-2021 17:00:00
    idx_day = find(day_number==day_no);
    datetime_vec_for_day = unique_datetime_facs(idx_day);
    % e.g.  datetime_vec also includes 01-aug-2021 18:00:00 but
    % datetime_vec_for_day does not i.e. on 01-aug-2021 at 18:00:00 the
    % footprints do not exist but the measurement and model data are
    % available
    % in this given example idx_dataset contain 12 ones
    idx_dataset = ismember(datetime_vec,datetime_vec_for_day);
    % e.g. for 01-aug-2021 06:00:00
    % wspd_diff_vec = -1.44917573574999	-2.39185820417389	-1.66070837114878	-1.76763578159003	-1.21656380348376	-1.84908797888784	-1.64727520432009	NaN	NaN	NaN	NaN	NaN	NaN
    % facs = 0.698386475141913	0.251896334068268	0.0405900514612419	0.00736567376516425	0.00106018601316670	0.000701279550246064	0	0	0	0	0	0	0
    % then multiply and sum them to get -1.6976
    % do this for all sharp hours for one day and take the average and call
    % it wspd_mu_day
    wspd_mu_day = mean(sum(wspd_diff_vec(idx_dataset,:).*facs(idx_day,:),2,'omitnan'));
    wspd_mean_vec(day_no) = wspd_mu_day;
    % do the same for wdir
    wdir_mu_day = mean(sum(wdir_diff_vec(idx_dataset,:).*facs(idx_day,:),2,'omitnan'));
    wdir_mean_vec(day_no) = wdir_mu_day;
    
    % now convert wspd and wdir differences from 2D to 1D 
    % e.g. there are 12 sharp hours footprint intensities available for
    % 01-aug-2021 then wspd_vec and wdir_vec will be of size 156x1, (12 sharp hours x 13 layers = 156)
    wspd_vec = reshape(wspd_diff_vec(idx_dataset,:)',[],1);
    wdir_vec = reshape(wdir_diff_vec(idx_dataset,:)',[],1);
    % do the same for weight vector
    weight_vec = reshape(facs(idx_day,:)',[],1) ;
    % calculate the standard deviations using the weight_vectors for each
    % day
    wspd_std_vec(day_no) = std(wspd_vec-wspd_mu_day,weight_vec,'omitnan');
    wdir_std_vec(day_no) = std(wdir_vec-wdir_mu_day,weight_vec,'omitnan');
end


%% SAVE
start_date = datetime("01-08-2021",'InputFormat','dd-MM-yyyy');
date_set = start_date + [0:1:39];
date_set = date_set';
% date_set is 40x1 datetime includes the days as datetime
% 01-aug-2021; 02-aug-2021; 03-aug.....
mean_WSPD_diff = round(wspd_mean_vec,3);
std_WSPD_diff = round(wspd_std_vec,3);
mean_WDIR_diff = round(wdir_mean_vec,3) ;
std_WDIR_diff = round(wdir_std_vec,3);

cd '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/Weighted_daily_WSPD_WDIR_differences';

My_table = table(date_set,mean_WSPD_diff,std_WSPD_diff,mean_WDIR_diff,std_WDIR_diff);
writetable(My_table,'weighted_daily_wspd_wdir_differences.csv');

