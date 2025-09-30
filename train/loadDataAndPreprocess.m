function [imdsTrain, imdsValidation] = loadDataAndPreprocess(datasetPath)
%% loadDataAndPreprocess.m — Build training/validation datastores with standard reader
% Summary:
%   Loads images from class folders, splits into train/validation, and assigns a robust ReadFcn.
%
% Requirements:
%   Image Processing Toolbox; Deep Learning Toolbox.
%
% Dependencies:
%   customImageReader.m
%
% Inputs:
%   datasetPath      char/string   Root of class subfolders.
%
% Outputs:
%   imdsTrain, imdsValidation   Split datastores with ReadFcn applied.
%
% Side effects:
%   None.
%
% Usage:
%   [tr,va] = loadDataAndPreprocess('data/wells');
%
% Notes:
%   - Ensures images are [50 50 1] single in [0,1] via customImageReader.

%%
    fullDatasetPath = fullfile(datasetPath);
    imds = imageDatastore(fullDatasetPath, 'IncludeSubfolders', true, 'LabelSource', 'foldernames', 'ReadFcn', @customImageReader);
    [imdsTrain, imdsValidation] = splitEachLabel(imds, 0.8, 'randomized');
end