function img = customImageReader(fileName)
%% customImageReader.m — Robust grayscale reader returning [50 50 1] single in [0,1]
% Summary:
%   Standardizes images for the CNN: convert RGB→gray if needed, resize to 50×50, cast to single.
%
% Requirements:
%   Image Processing Toolbox (for rgb2gray, imresize).
%
% Dependencies:
%   None (used as imds.ReadFcn).
%
% Inputs:
%   fileName      char/string   Image path (provided by datastore).
%
% Outputs:
%   img           single [50×50×1]  Normalized grayscale image in [0,1].
%
% Side effects:
%   None.
%
% Usage:
%   imds.ReadFcn = @(f) customImageReader(f);
%
% Notes:
%   - Always reshape to [50 50 1] to match network input layer.
%   - Avoids failures from unexpected 3-channel inputs.

n = 50; % width/length after resize; 

raw = imread(fileName);
if size(raw,3)==3
    raw = rgb2gray(raw); 
end
img = im2single(imresize(raw,[n n]));
img = reshape(img,[n n 1]);


end