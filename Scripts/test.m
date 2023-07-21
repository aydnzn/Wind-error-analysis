% ---------------------------------------------------------
% TUM - Technichal University of Munich
%
% Authors:  Aydin Uzun
% Date: 2022
% Purpose: test
% ---------------------------------------------------------
clear;
close all;
%% load wind speed and wind direction difference matrices
addpath('/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/Lidar_ERA5_representatives');
addpath('/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/Weighted_daily_WSPD_WDIR');
addpath('/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/Weighted_daily_WSPD_WDIR_differences');
T = readtable('weighted_daily_wspd_wdir_v3.csv');
T_diff = readtable('weighted_daily_wspd_wdir_differences_v3.csv');

dateset =T.date_set;

mu_lidar_wspd =T.mean_lidar_WSPD;
mu_era5_wspd =T.mean_era5_WSPD;

std_lidar_wspd =T.std_lidar_WSPD;
std_era5_wspd =T.std_era5_WSPD;

mu_lidar_wdir =T.mean_lidar_WDIR;
mu_era5_wdir =T.mean_era5_WDIR;

std_lidar_wdir =T.std_lidar_WDIR;
std_era5_wdir =T.std_era5_WDIR;



mu_wspd_diff =T_diff.mean_WSPD_diff;
mu_wdir_diff =T_diff.mean_WDIR_diff;

std_wspd_diff =T_diff.std_WSPD_diff;
std_wdir_diff =T_diff.std_WDIR_diff;


mu_lidar_wspd-mu_era5_wspd-mu_wspd_diff;
mu_lidar_wdir-mu_era5_wdir-mu_wdir_diff;

cd '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/Daily Gaussians';

for i=1:40

    f=figure;
    x=[-15:.1:15];
    lidar_wspd_dist=normpdf( x ,mu_lidar_wspd(i),std_lidar_wspd(i));
    plot(x,lidar_wspd_dist); hold on;
    era5_wspd_dist=normpdf( x ,mu_era5_wspd(i),std_era5_wspd(i));
    plot(x,era5_wspd_dist); hold on;
    diff_wspd_dist=normpdf( x ,mu_wspd_diff(i),std_wspd_diff(i));
    plot(x,diff_wspd_dist);
    set(gca,'Fontsize',14);
    xlabel('Wind speed [m/s]');
    ylabel('Dsitribution');
    legend('Lidar','ERA5','Lidar-ERA5');
    grid on;
    xticks([-15:1:15]);
    set(f, 'PaperPositionMode', 'auto', 'Units', 'Centimeters', 'Position', [0 0 20 20]);
    print(f,strcat('wspd_',datestr(dateset(i))),'-dpng');

    pause(0.5);

    g=figure;
    x2=[-60:0.1:360];
    lidar_wdir_dist=normpdf( x2 ,mu_lidar_wdir(i),std_lidar_wdir(i));
    plot(x2,lidar_wdir_dist); hold on;
    era5_wdir_dist=normpdf( x2 ,mu_era5_wdir(i),std_era5_wdir(i));
    plot(x2,era5_wdir_dist); hold on;
    diff_wdir_dist=normpdf( x2 ,mu_wdir_diff(i),std_wdir_diff(i));
    plot(x2,diff_wdir_dist);
    set(gca,'Fontsize',14);
    xlabel('Wind direction [deg]');
    ylabel('Dsitribution');
    legend('Lidar','ERA5','Lidar-ERA5');
    grid on;
    xticks([-60:30:360]);
    xlim([-60,360]);
    set(g, 'PaperPositionMode', 'auto', 'Units', 'Centimeters', 'Position', [0 0 20 20]);
    print(g,strcat('wdir_',datestr(dateset(i))),'-dpng');


end