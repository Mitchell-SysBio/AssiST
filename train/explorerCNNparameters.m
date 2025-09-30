%% Neural Network Architecture Testing Script
% explorerNNparameters.m — tests various CNN configurates to identify the
% best performing architecture on a given dataset
% Summary:
%   Explores a predefined hyperparameter grid (filters, filter sizes, learning rates, etc.),
%   trains small CNNs, evaluates on a validation set, and saves
%   models/metrics to .mat.
%
% Requirements:
%   MATLAB R2021a+; Deep Learning Toolbox; Parallel Computing Toolbox (optional).
%
% Dependencies:
%   trainAndEvaluateNetwork.m, parsave.m, loadDataAndPreprocess.m
%
% Outputs:
%   bestInfo      struct        Summary of best configuration and validation accuracy.
%   Creates ./Networks and subfolders; writes .mat models, training info, and CSV/figures.
%
% Notes:
%   - Expects images to be read via customImageReader → [50 50 1] single in [0,1].
%   - Enforces monotonically non-decreasing filter counts (n1 ≤ n2 ≤ n3).

%% Define parameters for the configurations
numFilters1 = [16, 8, 4,2];
numFilters2 = [32, 16, 8,4];
numFilters3 = [64, 32, 16,8];
filterSizes = [3, 5];
learningRates = [0.01, 0.001];

% Load and preprocess the dataset: Creates training and validation sets 
[imdsTrain, imdsValidation] = loadDataAndPreprocess('../prepare/training_setsRotate/');

% Initialize arrays to store results
numConfigs = length(numFilters1) * length(numFilters2) * length(numFilters3) * length(filterSizes) * length(learningRates);
results = struct('Config', cell(numConfigs, 1), 'TrainingAccuracy', zeros(numConfigs, 1), 'ValidationAccuracy', zeros(numConfigs, 1), 'ExternalAccuracy', zeros(numConfigs, 1), 'ExecutionTime', zeros(numConfigs, 1), 'TrainingLoss', zeros(numConfigs, 1));

% Creating Networks Folder
folderPath = './Networks';
if ~isfolder(folderPath)
    [ok,msg,msgID] = mkdir(folderPath);
    if ~ok
        error("Could not create folder '%s': %s (%s)", folderPath, msg, msgID);
    end
end

% Parallel loop for training and evaluating networks
parfor idx = 1:numConfigs
    fprintf('Training configuration %d...\n', idx);
    tic; % Start timer

    % Extract parameters for the current configuration
    [n1, n2, n3, fs, lr] = extractParameters(idx, numFilters1, numFilters2, numFilters3, filterSizes, learningRates);
    
    strConfig = sprintf('N1:%d, N2:%d, N3:%d, FS:%d LR:%g', n1, n2, n3, fs, lr);
    strFileName = sprintf('./Networks/net_%02d.mat', idx);
    if(n1<=n2 && n2<=n3) % don't run CNNs that more filters in initial layers than in deeper layers
        try
            % Define, train, and evaluate the network
            [net, info] = trainEvaluateNetwork(imdsTrain, imdsValidation,fs, n1, n2, n3, lr);

            % Save the trained network and training information
            parsave(strFileName, net, info); % Assuming parsave is implemented to handle parallel save operations

            % Store results
            results(idx).Config = strConfig;
            results(idx).Configuration = [n1, n2, n3, fs, lr]; 
            results(idx).fileName = strFileName;
            results(idx).Complexity = n1*n2*n3*fs;
            results(idx).TrainingAccuracy = info.TrainingAccuracy(end);
            results(idx).ValidationAccuracy = info.FinalValidationAccuracy;
            results(idx).TrainingLoss = info.TrainingLoss(end);
            
            [YPredVal, ~] = classify(net, imdsValidation);
            accuracyValidation = mean(YPredVal == imdsValidation.Labels);
            results(idx).ValidationAccuracyFull = accuracyValidation;

            [YPredTrain, ~] = classify(net, imdsTrain);
            accuracyTrain = mean(YPredTrain == imdsTrain.Labels);
            results(idx).TrainingAccuracyFull = accuracyTrain;
            
            results(idx).ExecutionTime = toc;
        catch ME
            warning('Failed to train or evaluate configuration %d: %s', idx, ME.message);
            results(idx).Config = sprintf('N1:%d, N2:%d, N3:%d, FS:%d LR:%g - Failed', n1, n2, n3, fs, lr);
            results(idx).ExecutionTime = toc;
        end
    else % network not trained
        % Store results
        results(idx).Config = strConfig;
        results(idx).fileName = 'empty';
        results(idx).Complexity = n1*n2*n3*fs;
        results(idx).Config = sprintf('N1:%d, N2:%d, N3:%d, FS:%d LR:%g', n1, n2, n3, fs, lr);
        results(idx).TrainingAccuracy = nan;
        results(idx).ValidationAccuracy = nan;
        results(idx).TrainingLoss = nan;
        results(idx).ExternalAccuracy = nan;
        results(idx).ExecutionTime = toc;

    end
end
save("results.mat","results")
%% Determine optimal network 
% Optimal networks are networks with top 3 validation accuracy 
[sortVal, valIDX] = sort([results.ValidationAccuracy],'descend', 'MissingPlacement','last');
bestNetIDX = valIDX(1:3);
bestColors = [220 47 2;244 140 6; 255 186 8]./255;
optimalNetInfo = results(valIDX(1));
load(sprintf("./Networks/net_%02d.mat", valIDX(1)))
save("optimalNetwork.mat", "net", "info","optimalNetInfo")

%% Plot metrics
nullCol = [169,169,169]./255;
f = figure;
subplot(1,3,1); hold on;
x = log([results.Complexity]); y = [results.ExecutionTime];
scatter(x,y,20, nullCol, 'filled');
s=[];
for i=1:3
    s=[s, scatter(x(bestNetIDX(i)), y(bestNetIDX(i)), 30,bestColors(i,:),'filled')];
end 
legend(s, {'first','second','third'}, 'location', 'northwest')
xlabel('ln(Complexity)'); ylabel('Execution Time');
box on; grid on;axis square;

subplot(1,3,2); hold on;
scatter(log([results.Complexity]),[results.ValidationAccuracy],20,nullCol,'filled');
scatter(log([results.Complexity]),[results.TrainingAccuracy],20,nullCol);
for i=[3,2,1]
    scatter(log(results(bestNetIDX(i)).Complexity),results(bestNetIDX(i)).ValidationAccuracyFull,30,bestColors(i,:),'filled');
    scatter(log(results(bestNetIDX(i)).Complexity),results(bestNetIDX(i)).TrainingAccuracyFull,30,bestColors(i,:));
end 
xlabel('ln(Complexity)'); ylabel('Accuracy');
legend({'Validation Accuracy','Training Accuracy'}, 'Location','southeast');
box on; grid on; axis square;

subplot(1,3,3); hold on;
scatter([results.TrainingAccuracyFull],[results.ValidationAccuracyFull],20,nullCol,'filled');
plot([.80 100],[.80 100],'--k');
for i=[3,2,1]
    scatter(results(bestNetIDX(i)).TrainingAccuracyFull,results(bestNetIDX(i)).ValidationAccuracyFull,30,bestColors(i,:),'filled');
end 
xlabel('Training Accuracy'); ylabel('Validation Accuracy');
box on; grid on;axis square;
xlim([.79 1]);ylim([.79,1]);


sgtitle(sprintf('First Configuration: %s - Validation Accuracy: %f\nSecond Configuration: %s - Validation Accuracy: %f\nThird Configuration: %s - Validation Accuracy: %f\n', results(bestNetIDX(1)).Config,results(bestNetIDX(1)).ValidationAccuracy, results(bestNetIDX(2)).Config,results(bestNetIDX(2)).ValidationAccuracy, results(bestNetIDX(3)).Config,results(bestNetIDX(3)).ValidationAccuracy))

%%
figure; hold on;
scatter([results.TrainingAccuracyFull],[results.ValidationAccuracyFull],15,log([results(~isnan([results.ValidationAccuracy])).Complexity]),'filled');
colormap(flipud(hot))
plot([.80 100],[.80 100],'--k');
xlabel('Training Accuracy'); ylabel('Validation Accuracy');
box on; grid on;axis square;
xlim([.9 1]);ylim([.9,1]);

%% Display the results
disp(struct2table(results));


