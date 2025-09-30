function [net, info] = trainEvaluateNetwork(imdsTrain, imdsValidation, filterSize, numFilters1, numFilters2, numFilters3, learnRate, varargin)
%% trainEvaluateNetwork.m — Build, train, and validate a CNN on 50×50×1 inputs
% Summary:
%   Constructs a small CNN, trains on imdsTrain, validates on imdsValidation, returns info.
%
% Requirements:
%   Deep Learning Toolbox.
%
% Dependencies:
%   None (expects datastores configured with customImageReader).
%
% Inputs:
%   imdsTrain, imdsValidation   imageDatastore  Training and validation sets.
%   filterSize                  double/vec      e.g., 3.
%   numFilters1,2,3             double          Filters per conv layer (n1 ≤ n2 ≤ n3).
%   learnRate                   double          Initial learning rate.
%   varargin                    struct          Optional: Epochs, MiniBatchSize, ValPatience, etc.
%
% Outputs:
%   net                         DAG/Series      Trained network.
%   info                        struct          Training record incl. accuracy/loss curves.
%
% Side effects:
%   Optional: saves checkpoints if opts.CheckpointPath provided.
%
% Usage:
%   [net,info] = trainEvaluateNetwork(tr,va,[3 3],32,64,128,1e-3,struct('Epochs',15));
%
% Notes:
%   - Input layer is imageInputLayer([50 50 1]); ensure readers produce matching size/type.
%   - ValidationFrequency ~ floor(numObservations/MiniBatchSize) for stable early stopping.

%% Hyperparameter defaults
p = inputParser;
p.addParameter('MiniBatchSize', 128, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('L2', 1e-4, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('MaxEpochs', 30, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('ValidationPatience', 25, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('CheckpointPath', '', @(x)ischar(x)||isstring(x));
p.parse(varargin{:});
mb          = p.Results.MiniBatchSize;
l2          = p.Results.L2;
maxEpochs   = p.Results.MaxEpochs;
valPatience = p.Results.ValidationPatience;
ckptPath    = string(p.Results.CheckpointPath);

% For reproducibility. This controls weight init and shuffling.
rng(0);

%% Define the model
% Input size is 50x50x1 (grayscale), matching customImageReader output.
inputSize = [50 50 1];
numClasses = numel(unique(imdsTrain.Labels));
layers = [
    imageInputLayer(inputSize, 'Name','input')

    % Block 1: Conv -> BN -> ReLU -> MaxPool
    convolution2dLayer(filterSize, numFilters1, 'Padding','same', 'Name','conv1')
    batchNormalizationLayer('Name','bn1')
    reluLayer('Name','relu1')
    maxPooling2dLayer(2, 'Stride',2, 'Name','pool1')

    % Block 2: Conv -> BN -> ReLU -> MaxPool
    convolution2dLayer(filterSize, numFilters2, 'Padding','same', 'Name','conv2')
    batchNormalizationLayer('Name','bn2')
    reluLayer('Name','relu2')
    maxPooling2dLayer(2, 'Stride',2, 'Name','pool2')

    % Block 3: Conv -> BN -> ReLU (no pool to preserve resolution)
    convolution2dLayer(filterSize, numFilters3, 'Padding','same', 'Name','conv3')
    batchNormalizationLayer('Name','bn3')
    reluLayer('Name','relu3')

    % Regularization head: dropout reduces co-adaptation/overfitting
    dropoutLayer(0.5, 'Name','dropout')

    % Classifier head: a hidden FC then final 8-class FC
    fullyConnectedLayer(2*numFilters3, 'Name','fc_hidden')
    reluLayer('Name','relu_hidden')
    fullyConnectedLayer(numClasses, 'Name','fc_out')
    softmaxLayer('Name','softmax')
    classificationLayer('Name','cls')
];

%% Training options
% SGDM is stable for small convnets. We add:
%  - MiniBatchSize: throughput/regularization tradeoff
%  - L2Regularization: weight decay for generalization
%  - ValidationData + ValidationFrequency: enables early stopping
%  - ValidationPatience: early-stop after N checks without improvement
%  - CheckpointPath: optional on-disk snapshots of network states

% Heuristic: validate every ~1 epoch (in steps). If imdsTrain is small, at least every 20 steps.
numTrainObs = numel(imdsTrain.Files);
valFreq = 20; % fallback default
if ~isempty(numTrainObs) && numTrainObs > 0
    stepsPerEpoch = max(1, floor(numTrainObs / mb));
    % Aim to validate ~2x per epoch for good curve granularity
    valFreq = max(10, round(stepsPerEpoch/2));
end

args = {
    'InitialLearnRate',      learnRate,...
    'MiniBatchSize',         mb,...
    'L2Regularization',      l2,...
    'MaxEpochs',             maxEpochs,...
    'Shuffle',               'every-epoch',...
    'ExecutionEnvironment',  'auto',...
    'Verbose',               false,...
    'Plots',                 'none',...
    'ValidationData',        imdsValidation,...
    'ValidationFrequency',   valFreq,...
    'ValidationPatience',    valPatience};

% Add checkpointing if path is given
if strlength(ckptPath) > 0
    if ~exist(ckptPath, 'dir') %#ok<EXIST>
        mkdir(ckptPath);
    end
    args = [args, {'CheckpointPath', char(ckptPath)}];
end

options = trainingOptions('sgdm', args{:});


%% Train
[net, info] = trainNetwork(imdsTrain, layers, options);

end