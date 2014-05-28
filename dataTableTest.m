%function dataTableTest
% last modified by wulx, 2014/2/25, 2014/5/25

clear all; close all; clc

% add lsramp path
if isempty(strfind(path, 'lsramp'))
    folderExist = exist('lsramp', 'dir');
    if folderExist == 7 % folder
        oldpath = addpath(fullfile(pwd, 'lsramp'));
    end
end

% load stroke time
D = load('data/DWELL_TIME.mat', 'strokeTime');
strokeTime = D.strokeTime;

% load optimization results
O = load('data/OPT_RESULT.mat', 'nStrks', 'C', 'Params');
nStrks = O.nStrks;
crts = O.C;
params = O.Params;

%% datatable1

% use struct array to save objective data
emptyCell = cell(1, nStrks);
DataTable = struct('pulseFreqs', emptyCell, ...
    'pulseNums', emptyCell, ...
    'nRamps', emptyCell, ...
    'crtSteps', emptyCell);

for  n = 1:nStrks
    
    % calculate frequencies and time sequencies
    [F, D, crtSteps] = datatable(params{n}, crts(n), strokeTime);
    
    nRamps = crts(n).num - 1;
    
    pulseNums = cell(1, nRamps);
    pulseFreqs = cell(1, nRamps);
    
    for j = 1:nRamps
        freqs_j = F{j};
        tseqs_j = D{j};
        
        %! using image processing techniques to reduce datatable
        freqs_rshift = [0 freqs_j(1:end-1)];
        
        neighbors = (freqs_j == freqs_rshift);
        
        [L, nLabels] = bwlabel(neighbors);
        
        L_lshift = [L(2:end) 0];
        
        labels = max(L, L_lshift);
        
        reducedFreqs = freqs_j;
        reducedTimeSeqs = tseqs_j;
        for i = 1:nLabels
            inds = find(labels == i);
            
            reducedFreqs(inds(2:end)) = nan;
            reducedTimeSeqs(inds(1)) = sum(tseqs_j(inds));
        end
        
        % clear all NaNs
        nans = isnan(reducedFreqs);
        reducedFreqs(nans) = [];
        reducedTimeSeqs(nans) = [];
        
        pulseNums{j} = round(reducedFreqs .* reducedTimeSeqs);
        pulseFreqs{j} = abs(reducedFreqs);
    end
    
    DataTable(n).pulseFreqs = pulseFreqs;
    DataTable(n).pulseNums = pulseNums;
    DataTable(n).nRamps = nRamps;
    DataTable(n).crtSteps = crtSteps;
    
end

%% datatable2 (reverse datatable1)

DataTable2 = DataTable(end:-1:1);

for m = 1:nStrks
    strk = DataTable2(m);
    
    % reverse critical steps
    DataTable2(m).crtSteps = strk.crtSteps(end:-1:1);

    for i = 1:strk.nRamps
        % reverse pulse numbers
        pulseNum = strk.pulseNums{end+1-i};
        DataTable2(m).pulseNums{i} = -pulseNum(end:-1:1);
        
        % reverse pulse frequencies
        pulseFreq = strk.pulseFreqs{end+1-i};
        DataTable2(m).pulseFreqs{i} = pulseFreq(end:-1:1);
    end

end

%%
% restore the old path
path(oldpath);

% save('data/DATA_TABLE.mat', 'DataTable', 'DataTable2');

