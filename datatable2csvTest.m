%% datatable2csv

format long g

% load caculated data
load('data/DWELL_TIME.mat', 'strokeTime', 'strkSet');
load('data/DATA_TABLE.mat', 'DataTable', 'DataTable2');

%%
nStrkTotal = numel(DataTable);
nStrkIdle = 2;
nStrkRunning = nStrkTotal - nStrkIdle;

RHome2 = DataTable(1).crtSteps(1); % DDT118
YHome2 = 1120; % DDT218

% 23000 pulses <=> 460 mm
targVal1 = 29440; % DDT1004
targVal2 = -targVal1; % DDT2004

freq = floor(targVal1 / strokeTime); % DDT1002 and DDT2002

strokeTimeInMs = 1000 * strokeTime;
delayTime = nStrkIdle * strokeTimeInMs;
roundTime = 2 * nStrkTotal * strokeTimeInMs;
gapTime = 2*nStrkIdle*strokeTimeInMs;
nRoundTrips = 8;

%#TIPS nUpStrks == nDownStrks or nDownStrks+1
nUpStrks = 2 * round(nStrkRunning / 2);
nDownStrks = 2*nStrkRunning - nUpStrks;

keys = {'strokeTime', 'delayTime', 'gapTime', 'roundTime', 'nRoundTrips', 'RHome2', 'YHome2', 'targVal1', 'targVal2', 'nUpStrks', 'nDownStrks', 'upFreq', 'downFreq'};
% keys = [keys, ...
%     arrayfun(@(n) ['nRampsUp' num2str(n)], 1:numel(nRampsUp), 'UniformOutput', false), ...
%     arrayfun(@(n) ['nRampsDown' num2str(n)], 1:numel(nRampsDown), 'UniformOutput', false)];
% seconds ==> milliseconds
vals = num2cell([strokeTimeInMs, delayTime, gapTime, roundTime, nRoundTrips, RHome2, YHome2, targVal1, targVal2, nUpStrks, nDownStrks, freq, freq]);

% params = [keys; vals];
csvwrite_with_headers('data/params.csv', vals, keys, 0, 0, '%8.0f');


%%
TF = @(m, n) [true(1, m) false(1, n)];

newDataTable = [DataTable(~TF(nStrkIdle, nStrkRunning)) DataTable2(TF(nStrkRunning, nStrkIdle))];

nRampsList = [newDataTable.nRamps]; % SV1011 and SV1021

rampLens = cellfun(@(pf) 4*(numel(pf)+1), [newDataTable.pulseFreqs]);
rampLocs = arrayfun(@(n) sum(rampLens(1:n)), 1:numel(rampLens));

csvwrite('data/nramps.csv', nRampsList);
csvwrite('data/lramps.csv', rampLocs);

%% ramps data

nRows = 2*max(cellfun(@(pf) numel(pf), [newDataTable.pulseFreqs])) + 2;
nCols = sum(nRampsList);
rampData = nan(nRows, nCols);

for i = 1:numel(nRampsList)
    for j = 1:nRampsList(i)
        pulFreqs_i = newDataTable(i).pulseFreqs{j};
        pulNums_i =  newDataTable(i).pulseNums{j};
        
        colIdx = sum(nRampsList(1:(i-1))) + j;
        
        if pulNums_i(1)<0
            rampData(1,colIdx) = hex2dec('3');
        else
            rampData(1,colIdx) = hex2dec('2');
        end
        
        nFreqs_i = numel(pulFreqs_i);
        for k = 1:nFreqs_i
           rampData(2*k, colIdx) =  pulFreqs_i(k);
           rampData(2*k+1, colIdx) = pulNums_i(k);
        end
        % tail
        rampData(2*nFreqs_i+2, colIdx) = 0;
    end
end

csvData = rampData(~isnan(rampData));
csvwrite('data/data.csv', csvData');

