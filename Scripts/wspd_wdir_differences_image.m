% ---------------------------------------------------------
% TUM - Technichal University of Munich
%
% Authors:  Aydin Uzun
% Date: 2022
% Purpose: Display WSPD and WDIR differences as an image with scaled colors
% ---------------------------------------------------------
clear;
close all;
%% load wind speed and wind direction difference matrices
addpath('/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/Lidar_ERA5_representatives');
load('WSPD_differences.mat'); % [lidar - era5 model] % 40x13x13 days x hours x layers
load('WDIR_differences.mat'); % [lidar - era5 model] % 40x13x13 days x hours x layers
load('Datetime_hourly.mat');
%% 
cd '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/WSPD_WDIR_differences_images';
for j=1:40 % j is the day number

    close all;
    % datetimes for day j
    datetimes_for_day_j=datetime_mat(j,:);
    % hours for day j will be on the x axis
    hours_for_day_j=hour(datetimes_for_day_j);

    f=figure;

    subplot(1,2,1);
    imAlpha=ones(size(squeeze(WSPD_diff(j,:,:))'));
    imAlpha(isnan(squeeze(WSPD_diff(j,:,:))'))=0;
    imagesc(hours_for_day_j,1:13,squeeze(WSPD_diff(j,:,:))','AlphaData',imAlpha)
    ylabel('Layers of the atmosphere');
    c=colorbar;
    caxis([-5,5]);
    ylabel(c,'WSPD differences [lidar - model] m s^{-1}');
    set(gca,'YDir','normal');
    yticks(1:13);
    xticks(hours_for_day_j);
    xlabel('Hours of the day');
    set(gca,'Fontsize',14);
    datetimes_for_day_j.Format = 'yyyy-MM-dd';
    set(gca,'color',0*[1 1 1]);


    subplot(1,2,2);
    imAlpha=ones(size(squeeze(wrapped_wdir_diff(j,:,:))'));
    imAlpha(isnan(squeeze(wrapped_wdir_diff(j,:,:))'))=0;
    imagesc(hours_for_day_j,1:13,squeeze(wrapped_wdir_diff(j,:,:))','AlphaData',imAlpha)
    ylabel('Layers of the atmosphere');
    c=colorbar;
    caxis([-50,50]);
    ylabel(c,'WDIR differences [lidar - model] deg');
    set(gca,'YDir','normal');
    yticks(1:13);
    xticks(hours_for_day_j);
    xlabel('Hours of the day');
    set(gca,'Fontsize',14);
    set(gca,'color',0*[1 1 1]);


    sgtitle(cellstr(datetimes_for_day_j(1,1)),'Fontsize',20);

    set(f, 'PaperPositionMode', 'auto', 'Units', 'Centimeters', 'Position', [0 0 40 20]);
    pause(0.5);
    print(f,string(strcat('wspd_wdir_differences_',cellstr(datetimes_for_day_j(1,1)))),'-dpng');


end