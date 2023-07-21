% ---------------------------------------------------------
% TUM - Technichal University of Munich
%
% Authors:  Aydin Uzun
% Date: 2022
% Purpose: Integrate "ERA5 hourly data on pressure levels" with "ERA5 hourly data on single levels"
% ---------------------------------------------------------
clear;
close all;
clc;
%%
addpath('/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/ERA5-sl');
addpath('/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/ERA5-pl');
%% read "ERA5 hourly data on single levels"
load('u100.mat'); % u100_mat_sl 44x24
load('v100.mat');  % v100_mat_sl 44x24
load('u10.mat'); % u10_mat_sl 44x24
load('v10.mat');  % v10_mat_sl 44x24
load('gph_sl.mat'); % gph_mat_sl 44x24

%% read "ERA5 hourly data on pressure levels"
load('ERA5_z.mat'); % z_mat 44x37x24
load('ERA5_v.mat'); % v_mat 44x37x24
load('ERA5_u.mat'); % u_mat 44x37x24
load('ERA5_level.mat'); % level_mat 44x37
g =9.80665;
gph_mat = z_mat/g; % gph_mat 44x37x24

%% Initialize 44 = days, 24 = hours, 37 = levels
U = zeros(44,24,37);
GPH = zeros(44,24,37);
V = zeros(44,24,37);
%% Convert u_mat, gph_mat, v_mat [44x24x37] to U,V,GPH [44x24x37]
for i = 1:44
    for j=1:24
U(i,j,:) = u_mat(i,:,j);
GPH(i,j,:) = gph_mat(i,:,j);
V(i,j,:) = v_mat(i,:,j);

    end
end

%% Inıtialize extended matrices 
GPH_ext = zeros(44,24,39);
U_ext = zeros(44,24,39);
V_ext = zeros(44,24,39);
sorted_GPH_ext = zeros(44,24,39);
sorted_U_ext = zeros(44,24,39);
sorted_V_ext = zeros(44,24,39);

%% Assign extended matrices
U_ext(:,:,1:37) = U;
U_ext(:,:,38) = u100_mat_sl;
U_ext(:,:,39) = u10_mat_sl;

V_ext(:,:,1:37) = V;
V_ext(:,:,38) = v100_mat_sl;
V_ext(:,:,39) = v10_mat_sl;

% add 10m and 100m to single level height
gph_10m = gph_mat_sl+10;
gph_100m = gph_mat_sl+100;

GPH_ext(:,:,1:37) = GPH;
GPH_ext(:,:,38) = gph_100m;
GPH_ext(:,:,39) = gph_10m;

%% Sort according to heights 
for i = 1:44
    for j=1:24
[sorted_GPH_ext_ij, idx ] = sort(squeeze(GPH_ext(i,j,:)));
sorted_GPH_ext(i,j,:) = sorted_GPH_ext_ij;
sorted_U_ext(i,j,:) = U_ext(i,j,idx);
sorted_V_ext(i,j,:) = V_ext(i,j,idx);
    end
end

%% Save
cd '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/ERA5-integrated';
save('ERA5_GPH_ext.mat','sorted_GPH_ext'); % 44x24x39 [1 to 39 increasing heights] 
save('ERA5_U_ext.mat','sorted_U_ext'); % 44x24x39 [1 to 39 increasing heights] 
save('ERA5_V_ext.mat','sorted_V_ext'); % 44x24x39 [1 to 39 increasing heights] 
cd '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Scripts';



