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

n = 3; % microsteps per full step, number of periods per ionBeamWidth
nTiers = 2*n;
ionBeamWidth = 30; % mm @ionBeamWidth ------------------------------------@
ionBeamWidth = round(ionBeamWidth/nTiers) * nTiers;
disp(['Ion beam width is revised as ' num2str(ionBeamWidth)])

leafWidth = 60; % mm
subHeight = 400; % mm
subWidth = 400; % mm
vStroke = leafWidth + subHeight; % vertical stroke

xPeriod = ionBeamWidth / n;

% stroke sets
%   data: depth map
%   padding: [top right bottom left]'
%   shear: shear along X-axis
%     [ 1   0   0
%      shx  1   0
%       0   0   1 ];
%   mode: bottom-UP or top-DOWN
strkSet = struct('data', {}, ...
    'padding', {}, ...
    'indices', {}, ...
    'nStrks', {}, ...
    'shear', {}, ...
    'mode', {});

avgDepth = depth / nTiers;
shx = 0.5*xPeriod / vStroke;

% strokes: 1..nTiers, nTiers+1..nTotalStrks
%   start from nTiers, the width of max-tier-area increases 0.5*xPeriod with
%   every new stroke adding.
nTotalStrks = nTiers + ceil(subWidth / (0.5*xPeriod));

% raster configuring -----------------------------------------------------%
for i = 1:n
    strkSet(2*i-1).mode = 'UP';
    strkSet(2*i-1).shear = -shx;
    %strkSet(2*i-1).data = avgDepth;
    strkSet(2*i-1).indices = (2i-1):nTiers:nTotalStrks;
    strkSet(2*i-1).nStrks = numel(strkSet(2*i-1).indices);
    strkSet(2*i-1).padding = [
        0.5*leafWidth;
        (strkSet(2*i-1).nStrks-1)*ionBeamWidth - subWidth;
        0.5*leafWidth;
        (nTiers-2*i+1)/nTiers * ionBeamWidth];
    
    strkSet(2*i).mode = 'DOWN';
    strkSet(2*i).shear = shx;
    %strkSet(2*i).data = avgDepth;
    strkSet(2*i).indices = 2i:nTiers:nTotalStrks;
    strkSet(2*i).nStrks =  numel(strkSet(2*i).indices);
    strkSet(2*i).padding = [
        0.5*leafWidth;
        (strkSet(2*i).nStrks-1)*ionBeamWidth - subWidth;
        0.5*leafWidth;
        (nTiers-2*i)/nTiers * ionBeamWidth];
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
            
            if sum(rowIdx) > 0
                ribbons(r, n) = mean( row_r(rowIdx) );
                rasterMap(r, colIdx(rowIdx)) = ribbons(r, n);
            end
        end
    end
    
    xform2 = [1, 0, 0; -strkSet(j).shear, 1, 0; 0, 0, 1];
    tform2 = maketform('affine', xform2);
    
    xdata = strkSet(j).padding(4) + [1 subWidth];
    if strcmp(strkSet(j).mode, 'UP')
        xdata = xdata + 0.5*xPeriod;
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



