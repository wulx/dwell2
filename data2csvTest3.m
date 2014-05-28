%% save as .csv
% last modified by wulx, 2014/3/6, 2014/3/10
% all tests PASS!
format long g

% load caculated data
load('data/DWELL_TIME.mat', 'strokeTime', 'strkSet');
load('data/DATA_TABLE.mat', 'DataTable', 'DataTable2');

% #TMP wulx, 2013/3/7
% select the last 4 up strokes and last 4 down strokes
% DataTable = DataTable(~~[0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1]);

% #TMP wulx, 2013/3/10
% select the first 4 pair (or up and down strokes)
first4 = ~[0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1];
DataTable = [DataTable2(~first4), DataTable(first4)];

% first1 = ~[0 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1];
% last1 = ~[1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 0];
% DataTable = [DataTable2(first1), DataTable(last1)];

first8 = [true(1, 8) false(1, 8)];
newDataTable = DataTable

RHome2 = DataTable(1).crtSteps(1); % DDT118
YHome2 = 1120; % DDT218

% 23000 pulses <=> 460 mm
targVal1 = 29440; % DDT1004
targVal2 = -targVal1; % DDT2004

freq = floor(targVal1 / strokeTime); % DDT1002 and DDT2002

isUpStrk = arrayfun(@(d) isequal(DataTable(d).strkDirs, 'up'), 1:numel(DataTable));
isDownStrk = ~isUpStrk;

upStrokes = DataTable(isUpStrk);
downStrokes = DataTable(isDownStrk);

% cross over
Strokes = DataTable;
Strokes(1:2:end) = upStrokes;
Strokes(2:2:end) = downStrokes;

nStrokes = numel(Strokes);
nUpStrks = sum(isUpStrk); % SV1010
nDownStrks = nStrokes - nUpStrks; % SV1020

% add field: strokeTime, 2014/3/6
% add fields: delayTime, roundTime and nRoundTrips 2014/3/26
strokeTimeInMs = 1000 * strokeTime;
delayTime = 8*strokeTimeInMs;
roundTime = 32*strokeTimeInMs;
gapTime = roundTime - nStrokes*strokeTimeInMs;
nRoundTrips = 10;
keys = {'strokeTime', 'delayTime', 'gapTime', 'roundTime', 'nRoundTrips', 'RHome2', 'YHome2', 'targVal1', 'targVal2', 'nUpStrks', 'nDownStrks', 'upFreq', 'downFreq'};
% keys = [keys, ...
%     arrayfun(@(n) ['nRampsUp' num2str(n)], 1:numel(nRampsUp), 'UniformOutput', false), ...
%     arrayfun(@(n) ['nRampsDown' num2str(n)], 1:numel(nRampsDown), 'UniformOutput', false)];
% seconds ==> milliseconds
vals = num2cell([strokeTimeInMs, delayTime, gapTime, roundTime, nRoundTrips, RHome2, YHome2, targVal1, targVal2, nUpStrks, nDownStrks, freq, freq]);

% params = [keys; vals];

csvwrite_with_headers('params.csv', vals, keys, 0, 0, '%8.0f');

%%
nRampsList = [Strokes.nRamps]; % SV1011 and SV1021

rampLens = cellfun(@(pf) 4*(numel(pf)+1), [Strokes.pulseFreqs]);
rampLocs = arrayfun(@(n) sum(rampLens(1:n)), 1:numel(rampLens));

csvwrite('nramps.csv', nRampsList);
csvwrite('lramps.csv', rampLocs);

%% ramps data

nRows = 2*max(cellfun(@(pf) numel(pf), [Strokes.pulseFreqs])) + 2;
nCols = sum(nRampsList);
rampData = nan(nRows, nCols);

for i = 1:nStrokes
    for j = 1:nRampsList(i)
        pulFreqs_i = Strokes(i).pulseFreqs{j};
        pulNums_i = Strokes(i).pulseNums{j};
        
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
csvwrite('data.csv', csvData');

