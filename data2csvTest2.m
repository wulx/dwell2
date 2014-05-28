%% save as .csv
% last modified by wulx, 2014/1/18
% all tests PASS!

nScan = 2;

% load data table
load('data_table.mat', 'DataTable');
load('data_table2.mat', 'DataTable2');
load('dwell_time.mat', 'strokeTime');

S = [DataTable, DataTable2];
for n = 1:nScan-1
    S = [S, DataTable, DataTable2];
end

RHome2 = S(1).crtSteps(1); % DDT118
YHome2 = 4000; % DDT218, 80 mm

nStrokes = numel(S);

isUpStrk = arrayfun(@(d) isequal(S(d).strkDirs, 'up'), 1:numel(S));
isDownStrk = ~isUpStrk;

nUpStrks = sum(isUpStrk); % SV1010
nDownStrks = nStrokes - nUpStrks; % SV1020

nRampsList = [S.nRamps]; % SV1011 and SV1021

% 23000 pulses <=> 460 mm
targVal1 = 29460; % DDT1004
targVal2 = -targVal1; % DDT2004

freq = round(targVal1 / strokeTime); % DDT1002 and DDT2002

upStrokes = S(isUpStrk);
nRampsUp = nRampsList(isUpStrk);

downStrokes = S(isDownStrk);
nRampsDown = nRampsList(isDownStrk);


keys = {'RHome2', 'YHome2', 'targVal1', 'targVal2', 'nUpStrks', 'nDownStrks', 'upFreq', 'downFreq'};
keys = [keys, ...
    arrayfun(@(n) ['nRampsUp' num2str(n)], 1:numel(nRampsUp), 'UniformOutput', false), ...
    arrayfun(@(n) ['nRampsDown' num2str(n)], 1:numel(nRampsDown), 'UniformOutput', false)];
vals = num2cell([RHome2, YHome2, targVal1, targVal2, nUpStrks, nDownStrks, freq, freq, nRampsList]);

% params = [keys; vals];

csvwrite_with_headers('params.csv', vals, keys);

%% up stroke

nUpRows = 2*max(cellfun(@(pf) numel(pf), [upStrokes.pulseFreqs])) + 2;
nUpCols = sum(nRampsUp);
upData = zeros(nUpRows, nUpCols);

%upHeaders = cell(1, nUpCols);

for upi = 1:nUpStrks
    for upj = 1:nRampsUp(upi)
        pulFreqs_i = upStrokes(upi).pulseFreqs{upj};
        pulNums_i = upStrokes(upi).pulseNums{upj};
        
        colIdx = sum(nRampsUp(1:(upi-1))) + upj;
        
        if pulNums_i(1)<0
            upData(1,colIdx) = hex2dec('3');
        else
            upData(1,colIdx) = hex2dec('2');
        end
        
        nFreqs_i = numel(pulFreqs_i);
        for upk = 1:nFreqs_i
           upData(2*upk, colIdx) =  pulFreqs_i(upk);
           upData(2*upk+1, colIdx) = pulNums_i(upk);
        end
        % tail
        %upData(2*nFreqs_i+2, colIdx) = 0;
        
        %upHeaders{colIdx} = sprintf('us%dr%d', upi, upj);
    end
end

csvwrite('updata.csv', upData');


%% down stroke

nDownRows = 2*max(cellfun(@(pf) numel(pf), [downStrokes.pulseFreqs])) + 2;
nDownCols = sum(nRampsDown);
downData = zeros(nDownRows, nDownCols);

%downHeaders = cell(1, nDownCols);

for di = 1:nDownStrks
    for dj = 1:nRampsDown(di)
        pulFreqs_i = downStrokes(di).pulseFreqs{dj};
        pulNums_i = downStrokes(di).pulseNums{dj};
        
        colIdx = sum(nRampsDown(1:(di-1))) + dj;
        
        if pulNums_i(1)<0
            downData(1,colIdx) = hex2dec('3');
        else
            downData(1,colIdx) = hex2dec('2');
        end
        
        nFreqs_i = numel(pulFreqs_i);
        for dk = 1:nFreqs_i
           downData(2*dk, colIdx) =  pulFreqs_i(dk);
           downData(2*dk+1, colIdx) = pulNums_i(dk);
        end
        % tail
        %downData(2*nFreqs_i+2, colIdx) = 0;
        
        %downHeaders{colIdx} = sprintf('ds%dr%d', di, dj);
    end
end

csvwrite('downdata.csv', downData');
