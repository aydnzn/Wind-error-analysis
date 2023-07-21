% ---------------------------------------------------------
% TUM - Technichal University of Munich
%
% Authors:  Aydin Uzun
% Date: 2022
% Purpose: convert "Scaling_factors.csv" to "scaling.mat"
% ---------------------------------------------------------
%%
% scaling.mat [40x13]
% normalized_scaling.mat [40x13]
% 40 days: from 20210801 to 20210909 increasing 1 to 40
% 13 layers: increasing 1[zagl=20],2[zagl=186]..... to 13[zagl=2220]
close all;
clear;

cd '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/Scaling_factors';
name_of_csv = "Scaling_factors.csv";
T = readtable(name_of_csv);
values = table2array(T(:,2:end));
scaling = values;
save('scaling.mat','scaling');

% normalized_scaling
normalized_scaling= zeros(40,13);
for i=1:40
    normalized_scaling(i,:) = values(i,:)/sum(values(i,:));
end
save('normalized_scaling.mat','normalized_scaling');
