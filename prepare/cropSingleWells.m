function [wellCenters, radius, img2]=cropSingleWells(img,strName, outputFolder, radiusCrop, tfRotate,tfFlip)
%% cropSingleWells.m — Interactive plate registration and per-well cropping
% Summary:
%   User clicks A1 and H12 corners on a plate image; function infers the grid and crops all 96 wells.
%
% Requirements:
%   Image Processing Toolbox (for strel/imresize).
%
% Dependencies:
%   None.
%
% Inputs:
%   img            double       image
%   strName        char/string  File name of 96 well image
%   outputFolder  char/string   Destination folder for cropped wells.
%   radiusCrop    double        Crop radius by scale factor around each well center.
%   tfRotate      logical       Rotate final crops (default: true).
%   tfFlip        logical       Flip final crops (default: true).
%
% Outputs:
%   wellCenters   double [96×4] [x y row col] center coordinates and indices.
%   radius        double        radius used to crop well
%   img2          double        edited image
%
% Side effects:
%   Writes 96 cropped images with well naming convention.
%
% Usage:
%   centers = cropSingleWells(imgage,'P001.png','output_folder',60,false,false);
%
% Notes:
%   - Click order matters: select A1 (top-left) then H12 (bottom-right).

% house keeping
idRows = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};
idColumns = {'1','2','3','4','5','6','7','8','9','10','11','12'};
wellNames = cell(96,1);
wellCenters = nan(96,4);

if(nargin<6), tfFlip = true; end
if(nargin<5), tfRotate = true; end

% adjust image (rotate and flip)
if(tfFlip), img = fliplr(img); end
if(tfRotate), img = imrotate(img,90); end

img2 = imadjust(img); 
% adjust the image so the short axis is 1000 pixels
shortSide = min(size(img2,1), size(img2,2)); % to get short side because for RGB the min would be # of channels
resizeFactor = 1000 / shortSide;
img2 = imresize(img2,resizeFactor);

% Display the image
h = figure('color','w');
movegui(h,'north')
imshow(img2); hold on;
title('Mark cell of A1 and H12 ...');


% Ask user to select the upper left and bottom right wells
[x, y] = ginput(2); % Get 2 points from the user

% decide on flip procedure

% Define the grid dimensions for a 96-well plate (8 rows, 12 columns)
rows = 8;
cols = 12;

% Calculate the distances between wells
rowDist = (y(2) - y(1)) / (rows - 1);
colDist = (x(2) - x(1)) / (cols - 1);

% Draw the grid
radius = round(radiusCrop*mean([rowDist,colDist]));
from_i_to_rc = nan(96,2);
wellNames = cell(96,1);

for row = 0:rows-1
    for col = 0:cols-1
        i = row*12+col+1;
        wellNames{i} = [idRows{row+1} idColumns{col+1}];
        wellCenters(i,1) = x(1) + col * colDist;
        wellCenters(i,2) = y(1) + row * rowDist;
        wellCenters(i,3) = row+1;
        wellCenters(i,4) = col+1;

        from_i_to_rc(i,[1 2]) = [row+1,col+1]; % record conversion from i index to row and col index
        rectangle('Position', [wellCenters(i,1)-radius, wellCenters(i,2)-radius, 2*radius, 2*radius],'Curvature', [1, 1], 'EdgeColor', "r", 'LineWidth', 2);
    end 
end

pause(2)
close(h)

% saving images 
rowSign = ["A", "B", "C", "D", "E","F","G","H"];
colSign = ["01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12"];
l = radius*2;
mask =strel('disk',radius,0);
for iWell = 1:96
    rect = [wellCenters(iWell,1)-radius, wellCenters(iWell,2)-radius, l, l];
    curWellImage= imcrop(img2, rect);
    if size(curWellImage,3)==3 % Making it gray scale if RGB so crop correctly 
        curWellImage = rgb2gray(curWellImage);
    end
    curWellImage(~mask.Neighborhood) = 0;
    imgName = sprintf("%s/%s_%s%s.png", outputFolder, strName, rowSign(wellCenters(iWell,3)), colSign(wellCenters(iWell,4)));
    imwrite(curWellImage,imgName)
end

end
