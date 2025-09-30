%% constructSingleNetwork.m — Construct, train, and validate a neural network using user given configuration 
% Summary:
%   Constructs the network, trains on training set, and evaluates using
%   validation set
%
% Requirements:
%   Deep Learning Toolbox.
%
% Dependencies:
%   trainEvaluateNetwork.m, loadDataAndPreprocess.m, customImageReader.m, 
%
% Inputs:
%  % See User Inputs section
%
% Outputs:
%   net           DAG/Series    Final trained network.
%   info          struct        Training info for the final run.
%
% Side effects:
%   Writes network (.mat) and training plots to output folder.
%
% Notes:
%   - Uses imageInputLayer([50 50 1]); relies on customImageReader in datastores.
%% User Inputs

% --------------------------------------------------------------------------
% Input desired network configuration 

% numFilters1 <= numFilters2 <= numFilters3
numFilters1 = 8;  % number of filters in layer 1
numFilters2 = 8;  % number of filters in layer 2
numFilters3 = 32; % number of filters in layer 3

filterSize = 5;   % filter size: size of local regions to which neurons connect in the input 
learnRate = 0.01; % learning rate: recommended is 0.01 - if too low training can take a long time, if too high, then training might reach suboptimal result or diverge 

% --------------------------------------------------------------------------

% Folder location of training set 
trainingLocation = "../prepare/training_setsRotate/";

% Name of output mat file that will contain trained network and information
outputFileName = 'singleNetwork_2025-09-30.mat';

%% Load and preprocess the dataset: Creates training and validation sets 
[imdsTrain, imdsValidation] = loadDataAndPreprocess(trainingLocation);

%% Construct Network: define, train, and evaluate the network
if (numFilters1<=numFilters2 && numFilters2<=numFilters3)
    [net, info] = trainEvaluateNetwork(imdsTrain, imdsValidation,filterSize, numFilters1, numFilters2, numFilters3, learnRate);
else 
    disp('FATAL ERROR: numFilters is incorrect, must be numFilters1 <= numFilters2 <= numFilters 3')
end 

% save the trained network
save(outputFileName, 'net','info', 'imdsTrain', "imdsValidation","filterSize", "numFilters1", "numFilters2", "numFilters3", "learnRate");

%% Evaluate Network 
[YPred, ~] = classify(net, imdsValidation);
YTrue = imdsValidation.Labels;

% Compute the confusion matrix
[confMat, order] = confusionmat(YTrue, YPred);
accuracy = mean(YPred == YTrue);

%% Plot Confusion Matrix 
f=figure;
cm = confusionchart(confMat,order,'ColumnSummary','column-normalized','RowSummary','row-normalized');
title(sprintf("Accuracy = %.2f%%",100*accuracy))
saveas(f,"ConfusionChartSingleNetwork.jpg")
