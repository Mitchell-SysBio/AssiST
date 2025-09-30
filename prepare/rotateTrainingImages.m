%% rotateTrainingImages.m — Generate 90/180/270° augmented training images
% Summary:
%   Walks class folders and writes rotated variants for data augmentation.
%
% Requirements:
%   None (uses built-ins).
%
% Dependencies:
%   None.
%
% Inputs:
%   inputFolder     char/string   Root of original class folders.
%   outputFolder    char/string   Root for augmented images.
%
% Outputs:
%   Creates class subfolders under outputFolder; writes rotated images.
%
% Usage:
%   Use in "prepareTrainingSet" folder
%
% Notes:
%   - Skips non-folder entries


%% definitions
inputFolder = 'training_sets';
outputFolder = 'training_setsRotate';
cp = pwd;
% Creating Output Folder
if ~isfolder(outputFolder)
    [ok,msg,msgID] = mkdir(outputFolder);
    if ~ok
        error("Could not create folder '%s': %s (%s)", outputFolder, msg, msgID);
    end
end
originFolder = pwd();
%% iterate over subfolder, rotate and save images
disp(['Copying all folder from ', inputFolder, ' into ', outputFolder]);

cd(inputFolder); % switch to input folder

subFolders = dir();
nFolders = length(subFolders);

for iFolder = 1:nFolders
    curFolder = subFolders(iFolder).name;
    if strcmp(curFolder, '.') || strcmp(curFolder, '..') || strcmp(curFolder, '.DS_Store')
        continue
    end 
    cd(curFolder);
    imageFiles = dir("*.png");
    nImageFiles = length(imageFiles);
    classFolder = ['../../',outputFolder,'/',curFolder]; 
    if ~exist(classFolder, 'dir')
            mkdir(classFolder)
    end 
    fprintf('Copying %s\n', curFolder)
    for iFile = 1:nImageFiles
        curFile = imageFiles(iFile).name;
        tokens = split(curFile,'.');
        curString = tokens{1};
        img = imread(curFile);
        splits = split(curString,'_');
        for j = [0, 90, 180, 270]
            rotateImg = imrotate(img, j);
            imgName = sprintf("%s/%s/%s/%s_%i.png", cp,outputFolder,curFolder,curString, j);
            imwrite(rotateImg,imgName);
        end
    end
    cd('..');
end

cd(originFolder);

