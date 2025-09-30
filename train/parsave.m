function parsave(fname, net, info)
%% parsave.m — Save helper for use inside parfor loops
% Summary:
%   Wraps save() to allow writing variables to .mat from within a parfor loop.
%
% Requirements:
%   None.
%
% Dependencies:
%   None.
%
% Inputs:
%   fname      char/string   Destination .mat path.
%   net, info  variables     Variable names to store.
%
% Outputs:
%   None.
%
% Side effects:
%   Writes a .mat file at filePath with variable varName.
%
% Usage:
%   parsave('model.mat',net,info);
%
% Notes:
%   - Keeps parfor bodies clean and avoids broadcast workspace issues.

    save(fname, 'net', 'info');
end