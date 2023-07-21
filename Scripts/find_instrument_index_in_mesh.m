function [latitude_idx_instrument,longitude_idx_instrument] = find_instrument_index_in_mesh(LON,LAT,lon_instrument,lat_instrument)
%   find_instrument_index_in_mesh: This function finds the latitude and
%   longitude idx of the instrument location using the available latitude and
%   longitude grid and its exact location
% ---------------------------------------------------------
% TUM - Technichal University of Munich
%
% Authors:  Aydin Uzun
% Date: 2022
% Purpose: function
% ---------------------------------------------------------
% LAT : meshgrid, 2D matrix
% LON : meshgrid, 2D matrix
% lon_instrument : longitude of instrument (double)
% lat_instrument : latitude of instrument (double)
% distances of the exact instrument location to each grid
dist_mx = sqrt((LON - lon_instrument).^2 + (LAT - lat_instrument).^2);
% find the minimum distance
min_mx = min(min(dist_mx));
% find the index on the grid
[latitude_idx_instrument,longitude_idx_instrument]=find(dist_mx==min_mx);

end