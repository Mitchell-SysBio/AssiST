%% prepareTrainingSet.m — Build training set from full plates 
% Summary:
%   Loads full plate images, crops wells,
%
% Requirements:
%   Image Processing Toolbox.
%
% Dependencies:
%   cropSingleWells.m
%
% Inputs:
%   See User Inputs section 
%
% Outputs:
%   Cropped wells for each plate and places into output folder
%
% Usage:
%   Use in parent folder of folder that contains full plate images

%% User Inputs
% Name of folder that contains full plate images
inputFolder = './full_plate_images';

% File type of full plate images 
dataFormat = '.png';

% Name of output folder for cropped well images
outputFolder = './single_wells';

%% housekeeping 

files = dir([inputFolder '/*' dataFormat]); 
nFiles = size(files,1)

% Creating output Folder
if ~isfolder(outputFolder)
    [ok,msg,msgID] = mkdir(outputFolder);
    if ~ok
        error("Could not create folder '%s': %s (%s)",outputFolder, msg, msgID);
    end
end

%% iterate over all files in inputFolder and crop single wells
radiusCrop = 0.325; % factor by which to decrease the radius for cropping the single wells
iFile = 1;
while iFile <=nFiles 
    curFile = files(iFile).name;
    tokens = split(curFile,'.');
    curString = tokens{1};
    img = imread([inputFolder '/' curFile]);
    try
    [wellCenters, radius, img2]=cropSingleWells(img,curString,outputFolder, radiusCrop); % the magic happens here
    iFile = iFile + 1; 
    catch ME
        disp('ERROR: Must click A1 then H12!')
    end
end 


