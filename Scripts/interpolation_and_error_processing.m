% ---------------------------------------------------------
% TUM - Technichal University of Munich
%
% Authors:  Aydin Uzun
% Date: 2022
% Purpose: Interpolation and create representatives of the layers
% One can also create WSPD and WDIR plots for the entire set of days or for a chosen day-hour pair
% [See the end of the script]
% ---------------------------------------------------------

clear;
close all;
clc;

%% Load the integrated ERA5 data
addpath('/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/ERA5-integrated');
load('ERA5_GPH_ext.mat'); % sorted_GPH_ext 44x24x39
load('ERA5_U_ext.mat'); % sorted_U_ext 44x24x39
load('ERA5_V_ext.mat'); % sorted_V_ext 44x24x39
%% Load hourly Lidar data
addpath('/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/Lidar_hourly');
load('Lidar_WSPD_hourly.mat'); % WSPD_day_by_day_hourly 44x24x149
load('Lidar_WDIR_hourly.mat'); % WDIR_day_by_day_hourly 44x24x149
load('datetime_hourly.mat'); % reshaped_datetime_era5 44x24
load('Lidar_height_hourly.mat'); % height_day_by_day_hourly 44x24x149
%% ignore 29 Jul, 30 Jul, 31 Jul, 10 Sep
% and use easier variable names like U,V and GPH
U = sorted_U_ext(4:end-1,:,:);
V = sorted_V_ext(4:end-1,:,:);
GPH = sorted_GPH_ext(4:end-1,:,:);

%% ignore some irrelevant data
% for this day for the first 6 hours the U and V components from ERA5 model
% were unrealistic
% After discussion it was decided to set them to zero.
U(10,1:6,:) = zeros(6,39);
V(10,1:6,:) = zeros(6,39);

%% again ignore 29 Jul, 30 Jul, 31 Jul, 10 Sep
% and use easier variable names like WDIR,WSPD and height and datetime_mat
WDIR = WDIR_day_by_day_hourly(4:end-1,:,:);
WDIR = wrapTo360(WDIR);
WSPD =WSPD_day_by_day_hourly(4:end-1,:,:);
height =height_day_by_day_hourly(4:end-1,:,:);
datetime_mat = reshaped_datetime_era5(4:end-1,:);
%% only consider data from 06:00 [idx=7] until 18:00[idx=19] for each day
U_new=U(:,7:19,:);
V_new=V(:,7:19,:);
GPH_new = GPH(:,7:19,:);
height_new = height(:,7:19,:);
WDIR_new = WDIR(:,7:19,:);
WSPD_new = WSPD(:,7:19,:);
datetime_mat_new = datetime_mat(:,7:19);
%% back to the same variable names
U = U_new; % 40x13x39 - ERA5
V = V_new; % 40x13x39 - ERA5
GPH = GPH_new; % 40x13x39 - ERA5
height = height_new; % 40x13x149 - LIDAR
WDIR = WDIR_new; % 40x13x149 - LIDAR
WSPD = WSPD_new; % 40x13x149 - LIDAR
datetime_mat = datetime_mat_new; % 40x13
%% INTERPOLATION
% initialize the matrices for interpolation
% ERA5 U wind components layer representatives
U_interp = zeros(40,13,13); % layerwise, 3th dim = 13 as there will be 13 atmospheric layers
% ERA5 V wind components layer representatives
V_interp = zeros(40,13,13); % layerwise, 3th dim = 13 as there will be 13 atmospheric layers
% Lidar WDIR layer representatives
WDIR_interp = zeros(40,13,13); % layerwise, 3th dim = 13 as there will be 13 atmospheric layers
% Lidar WSPD layer representatives
WSPD_interp = zeros(40,13,13); % layerwise, 3th dim = 13 as there will be 13 atmospheric layers

WSPD_era5 = zeros(40,13,13); % ERA5 WSPDS layer representatives
degree_mat = zeros(40,13,13); % ERA5 WDIR layer representatives

% release heights taken from rds document
relative_release = 20 + [0,166,335,507,681,859,1040,1224,1412,1603,1798,1997,2200];
% calculate atmospheric layer boundaries by calculating the mean of each
% relative release heights
relative_release_shifted = zeros(1,13);
relative_release_shifted(1:end-1) = relative_release(2:end);
layers = (relative_release + relative_release_shifted) /2; % layers is of size 1x13
layers(end) = 2933;
ext_layers = [0,layers];

% calculate the heights of interest
% in each atmospheric layer we will have 10 heights
% layer points is of size 1x130, as there are 13 atm. layers
% e.g. 1st relative release height is 20 m
% 2nd relative release height is 186 m
% hence the boundary is at 103 m
% and we are going to introduce linarly spaced 10 points from 0 to 103 m
% 9.3636   18.7273   28.0909   37.4545   46.8182   56.1818   65.5455   74.9091   84.2727   93.6364
% the next task will be to calculate wind speed components at these heights
layer_points = zeros(1,130);
for i=1:13
    b = linspace(ext_layers(i),ext_layers(i+1),12);
    layer_points((i-1)*10+1:i*10)=b(2:end-1);
end
%% create LIDAR WDIR WSPD interpolated - 130 points
% initialize % 40 days, 13 hours, 13o points
sin_comp_WDIR_interp_points = zeros(40,13,130);
cos_comp_WDIR_interp_points = zeros(40,13,130);

% was an other option to consider
% unit_sin_comp_WDIR_interp_points = zeros(40,13,130);
% unit_cos_comp_WDIR_interp_points = zeros(40,13,130);

% calculate sin and cos components using the WSPD and WDIR
sin_comp_WDIR = WSPD.* sin(WDIR*pi/180);
cos_comp_WDIR = WSPD.* cos(WDIR*pi/180);
% unit_sin_comp_WDIR =  sin(WDIR*pi/180);
% unit_cos_comp_WDIR =  cos(WDIR*pi/180);

for i =1:40 % days
    for j=1:13 % hours
        % do the interpolation
        % layer_points are query points
        sin_comp_WDIR_interp_points(i,j,:) = interp1(squeeze(height(i,j,:)),squeeze(sin_comp_WDIR(i,j,:)),layer_points','linear');
        cos_comp_WDIR_interp_points(i,j,:) = interp1(squeeze(height(i,j,:)),squeeze(cos_comp_WDIR(i,j,:)),layer_points','linear');
        % unit_sin_comp_WDIR_interp_points(i,j,:) = interp1(squeeze(height(i,j,:)),squeeze(unit_sin_comp_WDIR(i,j,:)),layer_points','linear');
        % unit_cos_comp_WDIR_interp_points(i,j,:) = interp1(squeeze(height(i,j,:)),squeeze(unit_cos_comp_WDIR(i,j,:)),layer_points','linear');
    end
end

% using interpolated sin and cos components derive the WSPD and WDIR
WSPD_interp_points = sqrt(sin_comp_WDIR_interp_points.^2 + cos_comp_WDIR_interp_points.^2);
WDIR_interp_points = atan2(sin_comp_WDIR_interp_points,cos_comp_WDIR_interp_points)*180/pi;
%% create LIDAR WSPD representative for the layer - 13 layers
for i=1:40 % days
    for j=1:13 % hours
        for k=1:13 % layers
            index = ((k-1)*10+1) : (k*10); % idx for each layer e.g. for the 1st layer 1:10
            is_finite = squeeze(isfinite(WSPD_interp_points(i,j,index))); % 10x1 logical vector; 1 at the idx at which WSPD_interp_points are finite; 0 if NAN
            sum_is_finite = sum(is_finite); % number of finite idx in each layer
            % if there are finite WSPD_interp_points in a layer set their
            % mean as the representative
            if sum_is_finite>0
                finite_idx = index(is_finite);
                WSPD_interp(i,j,k) = mean(squeeze(WSPD_interp_points(i,j,finite_idx)),'omitnan');
            else % if there is no finite WSPD_interp_points in one layer set the representative as NaN
                WSPD_interp(i,j,k) = NaN;
            end

        end
    end
end
%% create LIDAR WDIR representative for the layer - 13 layers
for i=1:40 % days
    for j=1:13 % hours
        for k=1:13 % layers
            index = ((k-1)*10+1) : (k*10); % idx for each layer e.g. for the 1st layer 1:10
            is_finite = squeeze(isfinite(WDIR_interp_points(i,j,index))); % 10x1 logical vector; 1 at the idx at which WDIR_interp_points are finite; 0 if NAN
            sum_is_finite = sum(is_finite);  % number of finite idx in each layer
            % if there are finite WDIR_interp_points in a layer
            if sum_is_finite>0
                finite_idx = index(is_finite);
                % find unit sin and cos components
                unit_sin = sin(squeeze(WDIR_interp_points(i,j,finite_idx))*pi/180);
                unit_cos = cos(squeeze(WDIR_interp_points(i,j,finite_idx))*pi/180);
                % find mean of unit sin component and cos component
                % then calculate the angle
                WDIR_interp(i,j,k) = atan2(mean(unit_sin,'omitnan'),mean(unit_cos,'omitnan'))*180/pi;
                % method 2 - ignored
                % sin_comp_WDIR_for_layer_k =  mean(squeeze(sin_comp_WDIR_interp_points(i,j,finite_idx)));
                % cos_comp_WDIR_for_layer_k =  mean(squeeze(cos_comp_WDIR_interp_points(i,j,finite_idx)));
                % WDIR_interp(i,j,k)=atan2(sin_comp_WDIR_for_layer_k,cos_comp_WDIR_for_layer_k)*180/pi;
            else % if there are no finite WDIR_interp_points in one layer set the representative as NaN
                WDIR_interp(i,j,k) = NaN;
            end

        end
    end
end
%% create U and V interpolated - 130 points
% initialize
U_interp_points = zeros(40,13,130); % 40 days, 13 hours, 130 points
V_interp_points = zeros(40,13,130); % 40 days, 13 hours, 130 points

for i =1:40 % days
    for j=1:13 % hours
        % do the interpolation
        % layer_points are query points
        U_interp_points(i,j,:)=interp1(squeeze(GPH(i,j,:)),squeeze(U(i,j,:)),layer_points','linear');
        V_interp_points(i,j,:)=interp1(squeeze(GPH(i,j,:)),squeeze(V(i,j,:)),layer_points','linear');
    end
end
%% create U representative for the layer - 13 layers

for i=1:40 % days
    for j=1:13 % hours
        for k=1:13 % layers
            index = ((k-1)*10+1) : (k*10); % idx for each layer e.g. for the 1st layer 1:10
            is_finite = squeeze(isfinite(U_interp_points(i,j,index))); % 10x1 logical vector; 1 at the idx at which U_interp_points are finite; 0 if NAN
            sum_is_finite = sum(is_finite); % number of finite idx in each layer
            if sum_is_finite>0
                % if there are finite WSPD_interp_points in a layer set their
                % mean as the representative
                finite_idx = index(is_finite);
                U_interp(i,j,k) = mean(squeeze(U_interp_points(i,j,finite_idx)),'omitnan');
            else %  if there are no finite U_interp_points in one layer set the representative as NaN
                U_interp(i,j,k) = NaN;
            end

        end
    end
end

%% create V representative for the layer - 13 layers
% do the same operations this time for V_interp_points
for i=1:40
    for j=1:13
        for k=1:13
            index = ((k-1)*10+1) : (k*10);
            is_finite = squeeze(isfinite(V_interp_points(i,j,index)));
            sum_is_finite = sum(is_finite);
            if sum_is_finite>0
                finite_idx = index(is_finite);
                V_interp(i,j,k) = mean(squeeze(V_interp_points(i,j,finite_idx)),'omitnan');
            else
                V_interp(i,j,k) = NaN;
            end

        end
    end
end
%% ERA5 - wind speed wind direction layerwise - 13 layers
% calculate representative WSPD [WSPD_era5] and WDIR [degree_mat] using the
% calculated U_interp and V_interp representatives

for i =1:40 % days
    for j=1:13 % hours
        WSPD_era5(i,j,:) = sqrt(squeeze(U_interp(i,j,:)).^2 + squeeze(V_interp(i,j,:)).^2 );
        degree_mat(i,j,:)=atan2(squeeze(-1*U_interp(i,j,:)),squeeze(-1*V_interp(i,j,:)))*(180/pi);
    end
end

%% ERA5 RAW DATA - wind speed wind direction all data points - 39 heights available
% initialize 
WSPD_era5_not_processed = zeros(40,13,39);
degree_mat_not_processed = zeros(40,13,39);
% we'll calculate wind speed and wind direction using only the provided raw
% era5 data
% need this for visualization purposes
for i =1:40
    for j=1:13
        WSPD_era5_not_processed(i,j,:) = sqrt(squeeze(U(i,j,:)).^2 + squeeze(V(i,j,:)).^2 );
        degree_mat_not_processed(i,j,:)=atan2(squeeze(-1*U(i,j,:)),squeeze(-1*V(i,j,:)))*(180/pi);
    end
end
degree_mat_not_processed = wrapTo360(degree_mat_not_processed);

%% ERA5 interpolated DATA - wind speed wind direction - 130 heights available 
% initialize
WSPD_era5_interp_points = zeros(40,13,130);
degree_mat_interp_points = zeros(40,13,130);
% we'll calculate wind speed and wind direction using only the interpolated
% era5 data [before calculating the layer representatives]
% need this for visualization purposes
for i =1:40 % days
    for j=1:13 % hours
        WSPD_era5_interp_points(i,j,:) = sqrt(squeeze(U_interp_points(i,j,:)).^2 + squeeze(V_interp_points(i,j,:)).^2 );
        degree_mat_interp_points(i,j,:)=atan2(squeeze(-1*U_interp_points(i,j,:)),squeeze(-1*V_interp_points(i,j,:)))*(180/pi);
    end
end
degree_mat_interp_points = wrapTo360(degree_mat_interp_points);

%% SAVE
cd '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/Lidar_ERA5_representatives';


WSPD_diff = WSPD_interp - WSPD_era5; % [lidar - era5 model] % 40x13x13 days x hours x layers
WDIR_diff =WDIR_interp - degree_mat; % [lidar - era5 model] % 40x13x13 days x hours x layers
% wrapped_wdir_diff = wrapTo180(WDIR_diff);
wrapped_wdir_diff = WDIR_diff;


% Wrap angle in degrees to [0 360]
% WDIR_interp = wrapTo360(WDIR_interp);
% degree_mat = wrapTo360(degree_mat);

save('WSPD_lidar_repr.mat','WSPD_interp'); % 40x13x13 days x hours x layers
save('WDIR_lidar_repr.mat','WDIR_interp'); % 40x13x13 days x hours x layers
save('WSPD_era5_repr.mat','WSPD_era5'); % 40x13x13 days x hours x layers
save('WDIR_era5_repr.mat','degree_mat'); % 40x13x13 days x hours x layers


save('WSPD_differences.mat','WSPD_diff');
save('WDIR_differences.mat','wrapped_wdir_diff');
save('Datetime_hourly.mat','datetime_mat');
cd '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Scripts';

%% WSPD
%% RUN THIS PART IF YOU WANT TO CREATE A TEST FIGURE FOR A SPECIFIC DATETIME;
close all;
day = 36; % pick the day \in {1,2,....40}
hour = 1; % pick the hour \in {1,2,....13}, 1 corresponds to 06:00
f = figure; 
scatter(squeeze(WSPD(day,hour,:)),squeeze(height(day,hour,:)),'blue','+','LineWidth',2);hold on; % raw measurement
scatter(squeeze(WSPD_interp_points(day,hour,:)),layer_points,'black','+'); hold on; % interpolated measurement
scatter(squeeze(WSPD_interp(day,hour,:)),relative_release','red','+','LineWidth',3); hold on; % measurement layer representatives
scatter(squeeze(WSPD_era5_not_processed(day,hour,:)),squeeze(GPH(day,hour,:)),'blue','Linewidth',2); hold on; % raw model
scatter(squeeze(WSPD_era5_interp_points(day,hour,:)),layer_points,'black'); hold on; % interpolated model
scatter(squeeze(WSPD_era5(day,hour,:)),relative_release','red','LineWidth',3); hold on; % model layer representatives
set(gca,'Fontsize',16);
ylim([0 2933]);
yline(layers);
yticks([relative_release,layers(end)]);
ylabel('Altitude [m]')
xlabel('Horizontal wind speed [ms^{-1}]');
title(datestr(datetime_mat(day,hour)));
legend('Lidar Measurement','interpolated Lidar Measurement','representative for layer - Lidar','ERA5 model','interpolated ERA5 model','representative for layer - ERA5','Location','northwest');
set(f, 'PaperPositionMode', 'auto', 'Units', 'Centimeters', 'Position', [0 0 40 30]);

%% WDIR
%% RUN THIS PART IF YOU WANT TO CREATE A TEST FIGURE FOR A SPECIFIC DATETIME;
close all;
day = 36; % pick the day \in {1,2,....40}
hour = 5; % pick the hour \in {1,2,....13}, 1 corresponds to 06:00
g=figure;
subplot(2,1,1);
polarscatter(squeeze(WDIR(day,hour,:))*(pi/180) ,squeeze(height(day,hour,:)),'magenta','x','LineWidth',2.5); hold on; % raw measurement
polarscatter(squeeze(WDIR_interp_points(day,hour,:))*(pi/180),layer_points,'black','x'); hold on; % interpolated measurement
polarscatter(squeeze(WDIR_interp(day,hour,:))*(pi/180),relative_release',60,'red','x','Linewidth',4); hold on; % measurement layer representatives
polarscatter(squeeze(degree_mat_not_processed(day,hour,:))*(pi/180),squeeze(GPH(day,hour,:)),'filled','magenta','LineWidth',2.5); hold on; % raw model
polarscatter(squeeze(degree_mat_interp_points(day,hour,:))*(pi/180),layer_points,'black'); hold on; % interpolated model
polarscatter(squeeze(degree_mat(day,hour,:))*(pi/180),relative_release',60,'red','LineWidth',4); hold on; % model layer representatives
rlim([0 layers(3)]); % only the first 3 layers
set(gca,'Fontsize',14);
rticks(layers(1:end-1)); % ticks on the radius axis here are layer boundaries
title(strcat('Wind directions - First 3 Layers',{' '},datestr(datetime_mat(day,hour))));
legend('Lidar Measurement','interpolated Lidar Measurement','representative for layer - Lidar','ERA5 model','interpolated ERA5 model','representative for layer - ERA5');
subplot(2,1,2);
polarscatter(squeeze(WDIR(day,hour,:))*(pi/180) ,squeeze(height(day,hour,:)),'blue','x'); hold on;
% commented interpolated part for visual enhancement 
% polarscatter(squeeze(WDIR_interp_points(day,hour,:))*(pi/180),layer_points,'black','+'); hold on;
polarscatter(squeeze(WDIR_interp(day,hour,:))*(pi/180),relative_release',60,'red','x','LineWidth',3); hold on;
polarscatter(squeeze(degree_mat_not_processed(day,hour,:))*(pi/180),squeeze(GPH(day,hour,:)),'blue'); hold on;
% commented interpolated part for visual enhancement 
% polarscatter(squeeze(degree_mat_interp_points(day,hour,:))*(pi/180),layer_points,'black'); hold on;
polarscatter(squeeze(degree_mat(day,hour,:))*(pi/180),relative_release',60,'red','LineWidth',3); hold on;
rlim([0 layers(end-1)]);
title(strcat('Wind directions',{' '},datestr(datetime_mat(day,hour))));
set(gca,'Fontsize',14);
rticks(layers(1:end-1));
legend('Lidar Measurement','representative for layer - Lidar','ERA5 model','representative for layer - ERA5');
set(g, 'PaperPositionMode', 'auto', 'Units', 'Centimeters', 'Position', [0 0 40 50]);

%% PRINT ALL WSPD FIGURES;
cd '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/WSPD_figures';
close all;
cd 
for day=1:40
    for hour = 1:13
        f = figure; 
        scatter(squeeze(WSPD(day,hour,:)),squeeze(height(day,hour,:)),'blue','+','LineWidth',2);hold on;
        scatter(squeeze(WSPD_interp_points(day,hour,:)),layer_points,'black','+'); hold on;
        scatter(squeeze(WSPD_interp(day,hour,:)),relative_release','red','+','LineWidth',3); hold on;
        scatter(squeeze(WSPD_era5_not_processed(day,hour,:)),squeeze(GPH(day,hour,:)),'blue','Linewidth',2); hold on;
        scatter(squeeze(WSPD_era5_interp_points(day,hour,:)),layer_points,'black'); hold on;
        scatter(squeeze(WSPD_era5(day,hour,:)),relative_release','red','LineWidth',3); hold on;
        set(gca,'Fontsize',16);
        ylim([0 2933]);
        yline(layers);
        yticks([relative_release,layers(end)]);
        ylabel('Altitude [m]')
        xlabel('Horizontal wind speed [ms^{-1}]');
        title(datestr(datetime_mat(day,hour)));
        legend('Lidar Measurement','interpolated Lidar Measurement','representative for layer - Lidar','ERA5 model','interpolated ERA5 model','representative for layer - ERA5','Location','northwest');
        set(f, 'PaperPositionMode', 'auto', 'Units', 'Centimeters', 'Position', [0 0 40 30]);
        print(f,strcat('wspd_',datestr(datestr(datetime_mat(day,hour)))),'-dpng');
    end
    pause(0.5);
    close all;
end

%% PRINT ALL WDIR;
cd '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/WDIR_figures';
close all;
for day=1:40
    for hour = 1:13
        g=figure;
        subplot(2,1,1);
        polarscatter(squeeze(WDIR(day,hour,:))*(pi/180) ,squeeze(height(day,hour,:)),'magenta','x','LineWidth',2.5); hold on;
        polarscatter(squeeze(WDIR_interp_points(day,hour,:))*(pi/180),layer_points,'black','x'); hold on;
        polarscatter(squeeze(WDIR_interp(day,hour,:))*(pi/180),relative_release',60,'red','x','Linewidth',4); hold on;
        polarscatter(squeeze(degree_mat_not_processed(day,hour,:))*(pi/180),squeeze(GPH(day,hour,:)),'filled','magenta','LineWidth',2.5); hold on;
        polarscatter(squeeze(degree_mat_interp_points(day,hour,:))*(pi/180),layer_points,'black'); hold on;
        polarscatter(squeeze(degree_mat(day,hour,:))*(pi/180),relative_release',60,'red','LineWidth',4); hold on;
        rlim([0 layers(3)]);
        set(gca,'Fontsize',14);
        rticks(layers(1:end-1));
        title(strcat('Wind directions - First 3 Layers',{' '},datestr(datetime_mat(day,hour))));
        legend('Lidar Measurement','interpolated Lidar Measurement','representative for layer - Lidar','ERA5 model','interpolated ERA5 model','representative for layer - ERA5');
        subplot(2,1,2);
        polarscatter(squeeze(WDIR(day,hour,:))*(pi/180) ,squeeze(height(day,hour,:)),'blue','x'); hold on;
        % polarscatter(squeeze(WDIR_interp_points(day,hour,:))*(pi/180),layer_points,'black','+'); hold on;
        polarscatter(squeeze(WDIR_interp(day,hour,:))*(pi/180),relative_release',60,'red','x','LineWidth',3); hold on;
        polarscatter(squeeze(degree_mat_not_processed(day,hour,:))*(pi/180),squeeze(GPH(day,hour,:)),'blue'); hold on;
        % polarscatter(squeeze(degree_mat_interp_points(day,hour,:))*(pi/180),layer_points,'black'); hold on;
        polarscatter(squeeze(degree_mat(day,hour,:))*(pi/180),relative_release',60,'red','LineWidth',3); hold on;
        rlim([0 layers(end-1)]);
        title(strcat('Wind directions',{' '},datestr(datetime_mat(day,hour))));
        set(gca,'Fontsize',14);
        rticks(layers(1:end-1));
        legend('Lidar Measurement','representative for layer - Lidar','ERA5 model','representative for layer - ERA5');
        set(g, 'PaperPositionMode', 'auto', 'Units', 'Centimeters', 'Position', [0 0 40 50]);
        print(g,strcat('wdir_',datestr(datestr(datetime_mat(day,hour)))),'-dpng');
        pause(1);
        close all;
    end
    close all;
end