% ---------------------------------------------------------
% TUM - Technichal University of Munich
%
% Authors:  Aydin Uzun
% Date: 2022
% Purpose: Plot weighted histograms
% Method: use hiswc.m to find weighted histogram count given number of bins
% Reference: Mehmet Suzen (2022). Generate Weighted Histogram (https://www.mathworks.com/matlabcentral/fileexchange/42493-generate-weighted-histogram), 
% MATLAB Central File Exchange. Retrieved July 24, 2022.
% ---------------------------------------------------------
clear;
close all;
%% load wind speed and wind direction difference matrices
addpath('/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/Lidar_ERA5_representatives');
addpath('/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Scripts');
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

cd '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/Weighted histograms';
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
    % do the same for wdir
    wdir_mu_day = mean(sum(wdir_diff_vec(idx_dataset,:).*facs(idx_day,:),2,'omitnan'));

    % now convert wspd and wdir differences from 2D to 1D
    % e.g. there are 12 sharp hours footprint intensities available for
    % 01-aug-2021 then wspd_vec and wdir_vec will be of size 156x1, (12 sharp hours x 13 layers = 156)
    wspd_vec = reshape(wspd_diff_vec(idx_dataset,:)',[],1);
    wdir_vec = reshape(wdir_diff_vec(idx_dataset,:)',[],1);
    % do the same for weight vector
    weight_vec = reshape(facs(idx_day,:)',[],1) ;
    % calculate the standard deviations using the weight_vectors for each
    % day
    wspd_std_day = std(wspd_vec-wspd_mu_day,weight_vec,'omitnan');
    wdir_std_day = std(wdir_vec-wdir_mu_day,weight_vec,'omitnan');
    

    datetime_vec_for_day.Format = 'yyyy-MM-dd';
    f=figure;
    subplot(1,2,1);
    [hist_w_wspd,v_interval_wspd]=histwc(wspd_vec,weight_vec,20);
    bar(v_interval_wspd, hist_w_wspd);
    xline(wspd_mu_day, 'Color', 'g', 'LineWidth', 2); hold on;
    xline( wspd_mu_day-wspd_std_day  , 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--'); hold on;
    xline( wspd_mu_day+ wspd_std_day, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
    set(gca,'Fontsize',16);
    legend('data','\mu','\mu \pm std');
    xlabel('WSPD LiDAR - model [m/s]');
    ylabel('weighted counts');

    subplot(1,2,2);
    [hist_w_wdir,v_interval_wdir]=histwc(wdir_vec,weight_vec,20);
    bar(v_interval_wdir, hist_w_wdir);
    xline(wdir_mu_day, 'Color', 'g', 'LineWidth', 2); hold on;
    xline( wdir_mu_day-wdir_std_day  , 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--'); hold on;
    xline( wdir_mu_day+ wdir_std_day, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
    set(gca,'Fontsize',16);
    legend('data','\mu','\mu \pm std');
    xlabel('WDIR LiDAR - model [deg]');
    ylabel('weighted counts');

    sgtitle(cellstr(datetime_vec_for_day(1,1)),'Fontsize',20);
    set(f, 'PaperPositionMode', 'auto', 'Units', 'Centimeters', 'Position', [0 0 30 15]);
    pause(0.5);
    print(f,string(strcat('weighted_histograms_',cellstr(datetime_vec_for_day(1,1)))),'-dpng');

end

