%% transform dwell time into etch depth
%
close all; clear all; clc

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

%%
S1 = load('data/DWELL_TIME.mat', 'dwellTime', 'strokeTime', 'strkSet');
strkSet = S1.strkSet;
strokeTime = S1.strokeTime;

S2 = load('data/OPT_RESULT.mat', 'optDwellTime');
optDwellTime = S2.optDwellTime;

S3 = load('data/ETCH_DEPTH.mat', 'depth');
etchDepth = S3.depth;

nTiers = numel(strkSet);
depthSup = zeros(subHeight, subWidth);
depthSup2 = depthSup;

nanPad = nan(size(strkSet(1).data));

for i = 1:nTiers
    strkSet(i).ribbons = strokeTime - optDwellTime(:, strkSet(i).indices);
    ribbons = meanEtchRate * strkSet(i).ribbons;
    
    depthSup = depthSup + strkSet(i).data;
    
    % padding ------------------------------------------------------------%
    preLeft = round([strkSet(i).padding(1), strkSet(i).padding(4)]);
    prePad = padarray(nanPad, preLeft, nan, 'pre');
    postRight = round([strkSet(i).padding(3), strkSet(i).padding(2)]);
    postPad = padarray(prePad, postRight, nan, 'post');
    
    % affine transform ---------------------------------------------------%
    xform = [1, 0, 0; strkSet(i).shear, 1, 0; 0, 0, 1];
    tform = maketform('affine', xform);
    shPad = imtransform(postPad, tform, 'nearest', 'FillValues', nan);
    
    % divide sheared padding area into stroke strips ---------------------%
    [shHeight, shWidth] = size(shPad);
    startIdx = 1:ionBeamWidth:shWidth;
    endIdx = [startIdx(2:end)-1, shWidth];
    
    rasterMap = nan(shHeight, shWidth);
    
    for n = 1:strkSet(i).nStrks
        for r = 1:shHeight
            rasterMap(r, startIdx(n):endIdx(n)) = ribbons(r, n);
        end
    end
    
    xform2 = [1, 0, 0; -strkSet(i).shear, 1, 0; 0, 0, 1];
    tform2 = maketform('affine', xform2);
    
    xdata = strkSet(i).padding(4) + [1 subWidth];
    if strcmp(strkSet(i).mode, 'UP')
        xdata = xdata + halfPeriod;
    end
    
    ydata = strkSet(i).padding(1) + [1 subHeight];
    
    strkSet(i).data = imtransform(rasterMap, tform2, 'nearest', 'FillValues', nan, ...
        'XData', xdata, 'YData', ydata);
    
    depthSup2 = depthSup2 + strkSet(i).data;
    
end

meanDepth = mean(etchDepth(:));
rmse1 = rmse(etchDepth, depthSup)/meanDepth;
rmse2 = rmse(etchDepth, depthSup2)/meanDepth;

figure

subplot(1, 3, 1)
imshow( mat2gray(etchDepth) )
colormap jet
% set(gca, 'XDir', 'reverse', 'YDir', 'normal', 'Visible', 'off', 'YAxisLocation', 'right')

subplot(1, 3, 2)
imshow( mat2gray(depthSup) )
colormap jet
set(gca, 'XDir', 'reverse', 'YDir', 'normal', 'Visible', 'off', 'YAxisLocation', 'right')
title(['CV(RMSE): ' num2str(rmse1)])

subplot(1, 3, 3)
imshow( mat2gray(depthSup2) )
colormap jet
set(gca, 'XDir', 'reverse', 'YDir', 'normal',  'YAxisLocation', 'right')
title(['CV(RMSE): ' num2str(rmse2)])

