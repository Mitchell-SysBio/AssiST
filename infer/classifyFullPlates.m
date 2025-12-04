%% classifyFullPlates.m — Batch classify full-plate images and export artifacts
% Summary:
%   Loads a trained network and classifies each well on plate images, saving phenotype,
%   QA, and metric images, plus per-plate outputs.
%
% Requirements:
%   Deep Learning Toolbox; Image Processing Toolbox.
%
% Dependencies:
%   inferPlatePhenotypes.m, inputDialog2.m
%
% Inputs:
%   See User definitions section below 
%
% Outputs:
%   Saved in output mat file: 
%       MICtable, MGCtable: MIC & MGC call for each drug and file in table
%       format
%       drugMIC, drugMGCMaxConc: MIC & MGC call for each drug and file
%       phenotypeImg: classified wells from neural network and user
%       reclassification
%       drug = drugs from plate map
% Images:
%   1) PhenotypeImages: phenotype color overlaid on plate after neural network classification
%   and user reclassification
%   2) QAImages: MIC and MGC overlaid on PhenotypeImages 
%   3) MetricImages: MIC and MGC overlaid on to plate with drug concentrations 
%
% Notes:
%   - Suppress noisy console prints; pre-create all output folders.
%   - Assumes net expects [50 50 1] inputs; direct classify paths must reshape accordingly.

%% User Definitions 
% -- Name of Excel file that contains platemap --
% If you do not want to determine MIC and MGC then put 'none'
    % Drug names must be in B3:M10 (no drug wells must be'PositiveControl')
    % Concentration per well must be in B14:M21
    % Drug type (bacteriostatic or bactericidal) must be in B25:M32 
platemapFileName = 'platemap.xlsx'; 

% Folder location of full plate images
dataFolder = './BMDimages/';

% Name of mat file that contains the best network 
netFileName = 'singleNetwork.mat'; 

% Name of output mat file to hold outputs 
resultsMatName = "results.mat"; % name of output variable 

% True if image needs to be flipped and rotated 90 degrees so A1 is top
% left and H12 is bottom right 
tfFlip = true; 
tfRotate = true;

% Growth classifications for class 
% classes= [1_bubbles, 2_empty, 3_pellet_tiny, 4_cloud_small, 5_pellet_small, 6_cloud, 7_pellet]
numNoGrowth = [1,2,3]; % number of classes corresponding to no growth 
numRestricted = [4,5]; % number of classes corresponding to restricted growth
numFull = [6,7]; % number of classes corresponding to full growth

% OPTIONAL: Flag wells with dubious classification
% This is the confidence threshold for probability of class inference.
% Anything below # will be flagged.
pClassThreshold = nan; % Put # between 0-1.Put "nan" if you want no flag.
 
%% Load data 
% load the neural network classifier
load(netFileName, 'net'); 

% Find images files 
dataFormat = '.png';
filesPng = dir([dataFolder '*' dataFormat]); 
dataFormat = '.jpg';
filesJpg = dir([dataFolder '*' dataFormat]); 
files = [filesPng; filesJpg];
nFiles = size(files,1);
fprintf('# of files: %i\n',nFiles)
% Create output folders 
folderPaths = ["./QAimages", "./PhenotypeImages","./MetricImages"];
for fp = 1:3
if ~isfolder(folderPaths(fp))
    [ok,msg,msgID] = mkdir(folderPaths(fp));
    if ~ok
        error("Could not create folder '%s': %s (%s)", folderPaths(fp), msg, msgID);
    end
end
end 

%% Classify images using trained network   
phenotypeImg = cell(nFiles,6); % preallocate  
drugMGCs = [];
drugMICs = [];
for iFile = 1:(nFiles)
    curFile = [dataFolder files(iFile).name];
    [~,curString,~] = fileparts(files(iFile).name);
    img = imread(curFile);

    % classify using neural network 
    [wellCenters, wellNames, h,h1,phenotypes,phenotypesReclass, reClassWell, reClass, drugMGC, drugMIC, drugs, pClass]=inferPlatePhenotypes(img, curString, platemapFileName, net, tfFlip, tfRotate, numNoGrowth,numRestricted,numFull, pClassThreshold);
    close(h);close(h1);
    phenotypeImg{iFile, 1} = curString; 
    phenotypeImg{iFile, 2} = wellCenters;  
    phenotypeImg{iFile, 3} = phenotypes; % unedited phenotypes
    phenotypeImg{iFile, 4} = phenotypesReclass; % edited phenotypes
    phenotypeImg{iFile, 5} = [reClassWell, reClass]; % storing if user needed to edit the phenotype
    phenotype96 = reshape(phenotypesReclass, 12,8)';
    phenotypeImg{iFile, 6} = phenotype96; % storing phenotypes in 96well format
    phenotypeImg{iFile, 7} = pClass; % posterier probabilites per image
    if ~strcmp(platemapFileName, 'none')
        % MGC: storing concentration with maximum growth
        drugMGCs(:,iFile) = drugMGC;

        % MIC: storing 1st concentration with no growth
        drugMICs(:,iFile)= drugMIC;
    end 
end 

%% Saving variables 
if strcmp(platemapFileName, 'none')
    save(resultsMatName,  "files", "phenotypeImg")
else
    MGCtable = array2table(drugMGCs);
    x1 = strrep(phenotypeImg(:,1),".png","");
    MGCtable.Properties.VariableNames= strrep(x1,".jpg","");
    MGCtable.Properties.RowNames=string(drugs);
    
    fprintf('# of files saved: %i\n',size(drugMGCs,2))

    MICtable = array2table(drugMICs);
    x1 = strrep(phenotypeImg(:,1),".png","");
    MICtable.Properties.VariableNames= strrep(x1,".jpg","");
    MICtable.Properties.RowNames=string(drugs);

    fprintf('# of files saved: %i\n',size(drugMGCs,2))

    save(resultsMatName, "MGCtable","MICtable","drugMGCs","drugMICs", "files","drugs", "phenotypeImg")
end
