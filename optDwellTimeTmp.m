%% optDwellTime test

clear all; close all; clc
% add lsramp path
if isempty(strfind(path, 'lsramp'))
    folderExist = exist('lsramp', 'dir');
    if folderExist == 7 % folder
        oldpath = addpath(fullfile(pwd, 'lsramp'));
    end
end

load('data/DWELL_TIME.mat', 'strokeTime', 'dwellTime', 'strkSet');
load('data/OPT_RESULT.mat', 'C', 'Params', 'Rmses', 'nStrks')

leafWidth = 60; % 60 mm
% a_max = 60000; % benchmark test result: 66333 pulse/s^2
timeStep = 0.001; % 1ms
stepAngleDeg = 1.8 / 8;

optDwellTime = nan( size(dwellTime) );

nTiers = numel(strkSet);
for j = 1:nTiers
    for k = strkSet(j).indices
        nw = C(k).num - 1;
        
        dwopt = @(params) dwellopt(params(1:nw), params(nw+(1:nw)), params(2*nw+(1:nw)), params(end-1), params(end), ...
            C(k), strokeTime, stepAngleDeg, leafWidth, timeStep);
        
        [r, optDwellTime_k] = dwopt(Params{k});
        
        if strcmp(strkSet(j).mode, 'DOWN')
            optDwellTime(:, k) = optDwellTime_k(end:-1:1);
        else
            optDwellTime(:, k) = optDwellTime_k;
        end
        
        disp(k)
    end
end

% restore the old path
if exist('oldpath', 'var')
    path(oldpath);
end

% save('data/OPT_RESULT.mat', 'C', 'Params', 'Rmses', 'nStrks', 'optDwellTime')
% disp('saved as data/OPT_RESULT.mat')
