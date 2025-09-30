function [a, b, c, d, e] = extractParameters(index, A, B, C, D, E)
%% extractParameters.m — Map linear config index to structured hyperparameters
% Summary:
%   Converts a single integer index into specific hyperparameter values based on predefined lists.
%
% Requirements:
%   None.
%
% Dependencies:
%   None.
%
% Inputs:
%   index         double            Linear index (1-based).
%   A,B,C,D,E     vector(double)    Fields are vectors/cells for each hyperparameter.
%
% Outputs:
%   a,b,c,d,e     double            Selected hyperparameters for training.
%
%
% Usage:
%   [a,b,c,d,e] = extractParameters(i,A,B,C,D,E);
%
% Notes:
%   - Validates index ranges; throws informative error on mismatch.    
%% Validate index
    maxIndex = length(A) * length(B) * length(C) * length(D) * length(E);
    if index < 1 || index > maxIndex
        error('Index out of valid range');
    end
    
    % Adjust index to be zero-based for calculations
    index = index - 1;
    
    % Calculate indices for each parameter array
    lenA = length(A);
    lenB = length(B);
    lenC = length(C);
    lenD = length(D);
    lenE = length(E);
    
    i = floor(index / (lenB * lenC * lenD * lenE)) + 1;
    j = floor(mod(index, lenB * lenC * lenD * lenE) / (lenC * lenD * lenE)) + 1;
    k = floor(mod(index, lenC * lenD * lenE) / (lenD * lenE)) + 1;
    l = floor(mod(index, lenD * lenE) / lenE) + 1;
    m = mod(index, lenE) + 1;
    
    % Extract parameters based on calculated indices
    a = A(i);
    b = B(j);
    c = C(k);
    d = D(l);
    e = E(m);
end