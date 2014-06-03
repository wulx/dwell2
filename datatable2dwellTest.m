%% datatable2dwell
close all; clear all; clc

% load caculated data
load('data/DWELL_TIME.mat', 'strokeTime', 'dwellTime', 'strkSet');
load('data/OPT_RESULT.mat', 'C', 'optDwellTime');
load('data/DATA_TABLE.mat', 'DataTable');

%%
leafWidth = 60; % 60 mm
stepAngleDeg = 1.8 / 8;
timeStep = 0.001; % 1 ms
nScan = ceil(strokeTime / timeStep);

steps = nan(1, nScan);
timeline = linspace(0, strokeTime, nScan);

nsum = @(n, list) sum(list(1:n));

runDwellTime = nan( size(dwellTime) );


nTiers = numel(strkSet);
for j = 1:nTiers
    for k = strkSet(j).indices
        pulseNums = [DataTable(k).pulseNums{:}];
        pulseFreqs = [DataTable(k).pulseFreqs{:}];
        
        timeSeqs = abs(pulseNums) ./ pulseFreqs;
        sumSteps = DataTable(k).crtSteps(1) + arrayfun(@(n) nsum(n, pulseNums), 1:numel(pulseNums));
        
        nStair = numel( pulseNums );
        for i = 1:(nStair-1)
            inds = timeline>=nsum(i-1, timeSeqs) & timeline<=nsum(i, timeSeqs);
            steps(inds) = sumSteps(i);
        end
        steps(timeline>=nsum(nStair-1, timeSeqs)) = sumSteps(nStair);
        
%         timeDiff = strokeTime - sum(timeSeqs)
        
        projWidths = step2width(steps, stepAngleDeg, leafWidth);
        
        runDwellTime1 = timecount(projWidths, strokeTime, C(k).scaleDivs)';
        runDwellTime1 = runDwellTime1(2:end-1); % remove boundaries

        mask = isnan(C(k).ogee);
        
        if strcmp(strkSet(j).mode, 'DOWN')
            runDwellTime1 = runDwellTime1(end:-1:1);
            mask = mask(end:-1:1);
        end

        runDwellTime1(mask) = nan;
        
        figure, plot(runDwellTime1);
        
        runDwellTime(:, k) = runDwellTime1;
    end
end

% save('data/ERROR_MODEL.mat', 'dwellTime', 'optDwellTime', 'runDwellTime')
