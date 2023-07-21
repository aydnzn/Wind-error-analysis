function im_to_rotate_rotated = rotate_around(angle,im_to_rotate,latitude_mx,longitude_mx)
% ---------------------------------------------------------
% TUM - Technichal University of Munich
%
% Authors:  Aydin Uzun
% Reference to Jan Motl (jan@motl.us)
% Date: 2022
% Purpose:  Rotate the footprint image around a point
% ---------------------------------------------------------
% rotate_around rotates an image.
%   im_to_rotate_rotated=rotate_around(angle,im_to_rotate,latitude_mx,longitude_m) rotates im_to_rotate around
%   the point [latitude_mx,longitude_mx] on the meshgrid by angle degrees. To rotate the image
%   clockwise, specify a negative value for angle.
%
%   Contributed by Jan Motl (jan@motl.us)
%   $Revision: 1.2 $  $Date: 2014/05/01 12:08:01 $

% Initialization.
[imageHeight, imageWidth] = size(im_to_rotate);
centerX = floor(imageWidth/2+1);
centerY = floor(imageHeight/2+1);

% distances to center
dy = centerY-latitude_mx;
dx = centerX-longitude_mx;

% How much would the "rotate around" point shift if the 
% image was rotated about the image center. 
[theta, rho] = cart2pol(-dx,dy);
[newX, newY] = pol2cart(theta+angle*(pi/180), rho);
shiftX = floor(longitude_mx-(centerX+newX));
shiftY = floor(latitude_mx-(centerY-newY));

% Pad the image to preserve the whole image during the rotation.
padX = abs(shiftX);
padY = abs(shiftY);
padded = padarray(im_to_rotate, [padY padX]);

% Rotate the image around the center.
rot = imrotate(padded, angle, 'bilinear', 'crop'); %Bilinear interpolation

% Crop the image.
im_to_rotate_rotated = rot(padY+1-shiftY:end-padY-shiftY, padX+1-shiftX:end-padX-shiftX, :);

end