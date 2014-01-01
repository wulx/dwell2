%% test on rastering and stitching method
close all; clear all; clc


%% load etch-depth map data
load data/ETCH_DEPTH.mat

z = depth; % in nm
x = 0.5:399.5; % mm, 400 points in total
y = 0.5:399.5; % mm
figure, meshc(x, y, z);
axis equal;


%% shear, pad and divide projected area

n = 2; % microsteps per full step, number of periods per ionBeamWidth
nTiers = 2*n;
ionBeamWidth = 60; % mm @ionBeamWidth ------------------------------------@
ionBeamWidth = round(ionBeamWidth/nTiers) * nTiers;
disp(['Ion beam width is revised as ' num2str(ionBeamWidth)])

leafWidth = 60; % mm
subHeight = 400; % mm
subWidth = 400; % mm
vStroke = leafWidth + subHeight; % vertical stroke

halfPeriod = ionBeamWidth / nTiers;

% provided that mean etch rate is 1 nm/s
meanEtchRate = 1;

% stroke sets
%   data: depth map
%   padding: [top right bottom left]'
%   indices: stroke indices
%   nStrks: number of strokes
%   ribbons: stroke profile
%   shear: shear along X-axis
%     [ 1   0   0
%      shx  1   0
%       0   0   1 ];
%   mode: bottom-UP or top-DOWN
%  
strkSet = struct('data', {}, ...
    'padding', {}, ...
    'indices', {}, ...
    'nStrks', {}, ...
    'ribbons', {}, ...
    'shear', {}, ...
    'mode', {});

avgDepth = depth / nTiers;
shx = halfPeriod / vStroke;

% strokes: 1..nTiers, nTiers+1..nTotalStrks
%   start from nTiers, the width of max-tier-area increases halfPeriod with
%   every new stroke adding.
nTotalStrks = nTiers + ceil(subWidth / halfPeriod);

% raster configuring -----------------------------------------------------%
for i = 1:n
    strkSet(2*i-1).mode = 'UP';
    strkSet(2*i-1).shear = -shx;

    strkSet(2*i-1).indices = (2*i-1):nTiers:nTotalStrks;
    strkSet(2*i-1).nStrks = numel(strkSet(2*i-1).indices);
    
    leftPadding = (nTiers - (2*i-1)) * halfPeriod;
    rightPadding = strkSet(2*i-1).nStrks*ionBeamWidth - halfPeriod - subWidth  - leftPadding;

    strkSet(2*i-1).padding = [
        0.5*leafWidth;
        rightPadding;
        0.5*leafWidth;
        leftPadding];
    
    strkSet(2*i).mode = 'DOWN';
    strkSet(2*i).shear = shx;

    strkSet(2*i).indices = (2*i):nTiers:nTotalStrks;
    strkSet(2*i).nStrks =  numel(strkSet(2*i).indices);
    
    leftPadding2 = (nTiers - 2*i) * halfPeriod;
    rightPadding2 = strkSet(2*i).nStrks*ionBeamWidth - halfPeriod - subWidth - leftPadding2;
    
    strkSet(2*i).padding = [
        0.5*leafWidth;
        rightPadding2;
        0.5*leafWidth;
        leftPadding2];
end

depthFitted = zeros(subHeight, subWidth);

for j = 1:nTiers
    % padding ------------------------------------------------------------%
    preLeft = round([strkSet(j).padding(1), strkSet(j).padding(4)]);
    prePad = padarray(avgDepth, preLeft, nan, 'pre');
    postRight = round([strkSet(j).padding(3), strkSet(j).padding(2)]);
    postPad = padarray(prePad, postRight, nan, 'post');
    
    % affine transform ---------------------------------------------------%
    xform = [1, 0, 0; strkSet(j).shear, 1, 0; 0, 0, 1];
    tform = maketform('affine', xform);
    shPad = imtransform(postPad, tform, 'nearest', 'FillValues', nan);
    
    % divide sheared padding area into stroke strips ---------------------%
    [shHeight, shWidth] = size(shPad);
    startIdx = 1:ionBeamWidth:shWidth;
    endIdx = [startIdx(2:end)-1, shWidth];
    
    ribbons = nan(shHeight, strkSet(j).nStrks);
    rasterMap = nan(shHeight, shWidth);

    for n = 1:strkSet(j).nStrks
        colIdx = startIdx(n):endIdx(n);
        for r = 1:shHeight
            row_r = shPad(r, colIdx);
            rowIdx = ~isnan( row_r );
            
            if any(rowIdx) > 0 % any() is better than sum() logically
                ribbons(r, n) = mean( row_r(rowIdx) );
                rasterMap(r, colIdx(rowIdx)) = ribbons(r, n);
            end
        end
    end

    strkSet(j).ribbons = ribbons / meanEtchRate;
    
    xform2 = [1, 0, 0; -strkSet(j).shear, 1, 0; 0, 0, 1];
    tform2 = maketform('affine', xform2);
    
    xdata = strkSet(j).padding(4) + [1 subWidth];
    if strcmp(strkSet(j).mode, 'UP')
        xdata = xdata + halfPeriod;
    end
    
    ydata = strkSet(j).padding(1) + [1 subHeight];
    
    strkSet(j).data = imtransform(rasterMap, tform2, 'nearest', 'FillValues', nan, ...
        'XData', xdata, 'YData', ydata);
    
    figure, mesh(strkSet(j).data);
    axis equal
    
    depthFitted = depthFitted + strkSet(j).data;
end


r = rmse(depth, depthFitted);

figure, mesh( depthFitted );
axis equal
title(['RMSE: ' num2str(r)])


%% convert depth map to time profiles

strips = [strkSet.ribbons];
stripNums = strips(~isnan(strips));

meanEtchTimeG = mean(stripNums); % global mean etch time
% medEtchTimeG = median(stripNums);
maxEtchTimeG = max(stripNums);
minEtchTimeG = min(stripNums);

% max tunable ratio of dwell time
maxTunableRatio = leafWidth/vStroke;

% storke time range:
%   maxEtchTime <= strokeTime <= minEtchTime/(1-maxTunableRatio)
strokeTimeRange = [maxEtchTimeG, minEtchTimeG/(1-maxTunableRatio)];

if strokeTimeRange(1) <= strokeTimeRange(2)
    etchTimeBase = meanEtchTimeG;
    strokeTimeOpt = etchTimeBase / (1 - 0.5*maxTunableRatio);
    
    if strokeTimeOpt>=strokeTimeRange(1) && strokeTimeOpt<=strokeTimeRange(2)
        strokeTime = strokeTimeOpt;
    else
        etchTimeBase = (minEtchTimeG + maxEtchTimeG) / 2;
        strokeTime = etchTimeBase / (1 - 0.5*maxTunableRatio);
    end
end

elapsedTime = nTiers*strokeTime;

figure, hold on;
plot(strips)
title('optimal stroke time')

plot(etchTimeBase*ones(vStroke, 1), 'k--')
xlim([1 vStroke])

plot(maxEtchTimeG*ones(vStroke, 1), 'k:');
plot(minEtchTimeG/(1-maxTunableRatio)*ones(vStroke, 1), 'k:');

plot(strokeTime*ones(vStroke, 1), 'r-')
plot((2*etchTimeBase-strokeTime)*ones(vStroke, 1), 'b:')


%% convert etch time map to dwell time map
% dwellTime = strokeTime - etchTime

dwellTime = nan(vStroke, nTotalStrks);

for k = 1:nTiers
    dwellTime(:, strkSet(k).indices) = strokeTime - strkSet(k).ribbons;
end

figure, hold on;
plot(dwellTime)
title('dwell time of all strokes')
xlim([1 vStroke])

maxDwellTime = maxTunableRatio * strokeTime;
plot(maxDwellTime*ones(vStroke, 1), 'r:')

% save('DWELL_TIME.mat', 'dwellTime', 'strokeTime', 'maxDwellTime', 'strkSet')

