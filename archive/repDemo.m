%% test on rastering and stitching method
close all; clear all; clc


%% load etch-depth map data
load data/ETCH_DEPTH.mat

z = depth; % in nm
x = 0.5:399.5; % mm, 400 points in total
y = 0.5:399.5; % mm
% figure, meshc(x, y, z);
% axis equal;


%% shear, pad and divide projected area

n = 1; % microsteps per full step, number of periods per ionBeamWidth
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

N = 6;
j = 1;
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
    rasterMap1 = nan(shHeight, shWidth);

    for n = 1:strkSet(j).nStrks
        colIdx = startIdx(n):endIdx(n);
        for r = 1:shHeight
            row_r = shPad(r, colIdx);
            rowIdx = ~isnan( row_r );
            
            if any(rowIdx) > 0 % any() is better than sum() logically
                ribbons(r, n) = mean( row_r(rowIdx) );
                rasterMap(r, colIdx(rowIdx)) = ribbons(r, n);
                
                if n == N
                    rasterMap1(r, colIdx(rowIdx)) = ribbons(r, n);
                end
            end
        end
    end

    
    %#TODO use etch rate vector to repace meanEtchRate
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
    data1 = imtransform(rasterMap1, tform2, 'nearest', 'FillValues', nan, ...
        'XData', xdata, 'YData', ydata);
    
%     figure, mesh(strkSet(j).data);
%     axis equal
    
    depthFitted = depthFitted + strkSet(j).data;


figure, hold on;
mesh(data1-13);
contourf(data1-13, 'EdgeColor','none');
box on
axis tight
lighting phong
% contourf(data1-14)
% axis equal

figure, plot(ribbons(:,N))
axis tight

% r = rmse(depth, depthFitted);
% 
% figure, mesh( depthFitted );
% axis equal
% title(['RMSE: ' num2str(r)])

ogee = ribbons(:,N);

[xx, ys] = spis(ogee);
opts = {'MINPEAKHEIGHT', -1, 'MINPEAKDISTANCE', 15, 'NPEAKS', 5};
[locs, ptype] = findpoints(ys, opts);

figure, hold on;
plot(xx, ys, 'k-')
axis tight

%mark all critical points at the spline
for i = 1:numel(locs)
    li = locs(i);
    ti = ptype(i);
    if ti > 0 % peaks
        plot(xx(li), ys(li), 'ro')
    elseif ti < 0 % valleys
        plot(xx(li), ys(li), 'bo')
    else % inflection points
        plot(xx(li), ys(li), 'ko')
    end
end

sec3 = locs(2):locs(3);
lineX = xx(sec3);
lineY = ys(sec3);
plot(lineX, lineY, 'g-')


