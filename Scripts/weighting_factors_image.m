% ---------------------------------------------------------
% TUM - Technichal University of Munich
%
% Authors:  Aydin Uzun
% Date: 2022
% Purpose: Display weighting factors as an image with scaled colors
% ---------------------------------------------------------
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
%% find the days
% day_number 435x1 double 1 for 01-aug-2021 etc.
day_number=ceil(days(unique_datetime_facs - '01-aug-2021'));
%% Create images and save them
cd '/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/Weighting factors images';
for j=1:max(day_number) % j is the day number

    close all;
    % find the idx for day j
    idx_for_day_j = find(day_number==j);
    % these are the intensities for day j
    intensity_image_for_day_j=facs(idx_for_day_j,:);
    % datetimes for day j
    datetimes_for_day_j=unique_datetime_facs(idx_for_day_j,1);
    % hours for day j will be on the x axis
    hours_for_day_j=hour(datetimes_for_day_j);

    f=figure; 
    imagesc(hours_for_day_j,1:13,intensity_image_for_day_j')
    ylabel('Layers of the atmosphere');
    c=colorbar;
    caxis([0,1]);
    ylabel(c,'Weighting factors');
    set(gca,'YDir','normal');
    yticks(1:13);
    xticks(hours_for_day_j);
    xlabel('Hours of the day');
    set(gca,'Fontsize',14);
    datetimes_for_day_j.Format = 'yyyy-MM-dd';
    title(cellstr(datetimes_for_day_j(1,1)));
    set(f, 'PaperPositionMode', 'auto', 'Units', 'Centimeters', 'Position', [0 0 20 15]);
    pause(0.5);
    print(f,string(strcat('weihgting_factors_',cellstr(datetimes_for_day_j(1,1)))),'-dpng');
    

end

