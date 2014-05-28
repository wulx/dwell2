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
nTiers = 2; % microsteps per full step
ionBeamWidth = 10; % mm @ionBeamWidth ------------------------------------@
leafWidth = 60; % mm
subHeight = 400; % mm
subLength = 400; % mm
vStroke = leafWidth + subHeight; % vertical stroke

% padding ----------------------------------------------------------------%
nUpStroke = ceil((ionBeamWidth + subLength) / ionBeamWidth);
upRange = nUpStroke * ionBeamWidth;
nDownStroke = ceil((ionBeamWidth/nTiers + subLength) / ionBeamWidth);
downRange = nDownStroke * ionBeamWidth;

preSize = [0, ionBeamWidth-(ionBeamWidth/nTiers)];
zPadPre = 0.5 * padarray(z, preSize, nan, 'pre');
postSize = [0, upRange-(ionBeamWidth+subLength)];
zPadPost = padarray(zPadPre, postSize, nan, 'post');

zUp = padarray(zPadPost, [leafWidth/2, 0], nan, 'both');

postSize2 = [0, downRange-(ionBeamWidth/nTiers + subLength)];
zPadPost2 = 0.5 * padarray(z, postSize2, nan, 'post');

zDown = padarray(zPadPost2, [leafWidth/2, 0], nan, 'both');

% figure, mesh(zUp);
% axis equal;
% figure, mesh(zDown);
% axis equal;

% shearing ---------------------------------------------------------------%
shx = (ionBeamWidth/nTiers) / vStroke;
% shx = 0.025;

xform = [ 1    0   0
         -shx  1   0
          0    0   1 ];
tform = maketform('affine', xform);
zr = imtransform(zUp, tform, 'nearest', 'FillValues', nan);  % rising strokes

xform2 = [ 1   0   0
          shx  1   0
           0   0   1 ];
tform2 = maketform('affine', xform2);
zf = imtransform(zDown, tform2, 'nearest', 'FillValues', nan);  % falling strokes

% figure, mesh(zr);
% axis equal;
% figure, mesh(zf);
% axis equal;

% divide sheared padding area into stroke strips -------------------------%
% divide rising strokes
rStartIx = 1:ionBeamWidth:upRange-ionBeamWidth+1;
rEndIx = [rStartIx(2:end)-1, upRange];

upRibbons = nan(size(zr,1), nUpStroke);

upDepthMap = nan(size(zr));
upDepthMap2 = upDepthMap;

nSigma = ionBeamWidth / 60 * 0.6;
% beamDist = beamdist(vStroke, ionBeamWidth, nSigma);
for i = 1:nUpStroke
    coli = rStartIx(i):rEndIx(i);
    zri = zr(:, coli);
    for ii = 1:size(zri,1)
        zrii = zri(ii, :);
        idx = ~isnan(zrii);
        c = sum(idx);
        if c > 0
            upRibbons(ii, i) = sum(zrii(idx)) / c;
            
            upDepthMap(ii, coli(idx)) = upRibbons(ii, i);
        end
    end
    
    upDepthMap2(:, coli) = upDepthMap(:, coli) .* beamdist(vStroke, ionBeamWidth, nSigma);
end

% figure;
% hUp = ribbon(upRibbons);
% xlabel('X')
% ylabel('Y')
% zlabel('Z')

zr2 = zr;
zr2(:, rStartIx) = nan;
zr2 = imtransform(zr2, tform2, 'FillValues', nan);

% zr3 = imtransform(zr, tform2, 'nearest', 'FillValues', nan, ...
%     'XData', ionBeamWidth+[1 size(depth,2)], 'YData', 0.5*leafWidth+[1 size(depth,1)]);

upDepthMap = imtransform(upDepthMap, tform2,  'nearest', 'FillValues', nan, ...
    'XData', ionBeamWidth+[1 size(depth,2)], 'YData', 0.5*leafWidth+[1 size(depth,1)]);
upDepthMap2 = imtransform(upDepthMap2, tform2,  'nearest', 'FillValues', nan, ...
    'XData', ionBeamWidth+[1 size(depth,2)], 'YData', 0.5*leafWidth+[1 size(depth,1)]);

figure, mesh(zr2);
title('up strokes')
axis equal

% figure, waterfall(upRibbons');

% divide falling strokes
fStartIx = 1:ionBeamWidth:downRange-ionBeamWidth+1;
fEndIx = [fStartIx(2:end)-1, downRange];

downRibbons = nan(size(zf,1), nDownStroke);

downDepthMap = nan(size(zf));
downDepthMap2 = downDepthMap;
for i = 1:nDownStroke
    zfi = zf(:, fStartIx(i):fEndIx(i));
    for ii = 1:size(zri,1)
        zfii = zfi(ii, :);
        idx = ~isnan(zfii);
        c = sum(idx);
        if c > 0
            downRibbons(ii, i) = sum(zfii(idx)) / c;
            
            coli = fStartIx(i):fEndIx(i);
            downDepthMap(ii, coli(idx)) = downRibbons(ii, i);
        end
    end
    
    downDepthMap2(:, coli) = downDepthMap(:, coli) .* beamdist(vStroke, ionBeamWidth, nSigma);
end

% figure;
% hDown = ribbon(downRibbons);
% xlabel('X')
% ylabel('Y')
% zlabel('Z')

zf2 = zf;
zf2(:, fStartIx) = nan;
zf2 = imtransform(zf2, tform, 'FillValues', nan);

% zf3 = imtransform(zf, tform, 'nearest', 'FillValues', nan, ...
%     'XData', [1 size(depth,2)], 'YData', 0.5*leafWidth+[1 size(depth,1)]);

downDepthMap = imtransform(downDepthMap, tform, 'nearest', 'FillValues', nan, ...
    'XData', [1 size(depth,2)], 'YData', 0.5*leafWidth+[1 size(depth,1)]);
downDepthMap2 = imtransform(downDepthMap2, tform, 'nearest', 'FillValues', nan, ...
    'XData', [1 size(depth,2)], 'YData', 0.5*leafWidth+[1 size(depth,1)]);

figure, mesh(zf2);
title('down strokes')
axis equal

% rastered depth map
% orgDepthMap = zr3 + zf3;
fitDepthMap = upDepthMap + downDepthMap;
fitDepthMap2 = upDepthMap2 + downDepthMap2;
% RMSE
r = rmse(depth, fitDepthMap);
r2 = rmse(depth, fitDepthMap2);

figure, mesh(upDepthMap);
axis equal
figure, mesh(downDepthMap);
axis equal

figure, mesh(fitDepthMap);
axis equal
title(['RMSE: ' num2str(r)])

figure, mesh(fitDepthMap2);
axis equal
title(['RMSE: ' num2str(r2)])



%% convert depth map to time map
% provided that etch rate is 1 nm/s
meanEtchRate = 1; % ------------------------------------------------------@

upStrips = upRibbons / meanEtchRate;
downStrips = downRibbons / meanEtchRate;
strips = [upStrips downStrips];
stripNums = strips(~isnan(strips));

meanEtchTimeG = mean(stripNums); % global mean etch time
maxEtchTimeG = max(stripNums);
minEtchTimeG = min(stripNums);

% max tunable ratio of dwell time
maxTunableRatio = leafWidth/vStroke;

% storke time range:
%   maxEtchTime < strokeTime < minEtchTime/(1-maxTunableRatio)
% elapsedTime = nTiers*strokeTime;
if maxEtchTimeG < minEtchTimeG/(1-maxTunableRatio)
    strokeTimeRange = [maxEtchTimeG, minEtchTimeG/(1-maxTunableRatio)];
    strokeTimeOpt = meanEtchTimeG / (1 - 0.5*maxTunableRatio);

    if strokeTimeOpt<strokeTimeRange(1)
        strokeTime = strokeTimeRange(1);
    elseif strokeTimeOpt>strokeTimeRange(2)
        strokeTime = strokeTimeRange(2);
    else
        strokeTime = strokeTimeOpt;
    end
end

figure, hold on;
plot(strips)
title('opt a stroke time')

plot(meanEtchTimeG*ones(vStroke, 1), 'k--')
xlim([1 vStroke])

plot(maxEtchTimeG*ones(vStroke, 1), 'k:');
plot(minEtchTimeG/(1-maxTunableRatio)*ones(vStroke, 1), 'k:');

plot(strokeTime*ones(vStroke, 1), 'r-')
plot((2*meanEtchTimeG-strokeTime)*ones(vStroke, 1), 'b:')


% convert etch time map to dwell time map
% dwellTime = strokeTime - etchTime
upDwellTime = strokeTime - upStrips;
downDwellTime = strokeTime - downStrips;

maxDwellTime = maxTunableRatio * strokeTime;

figure, hold on;
plot(upDwellTime)
title('dwell time of up strokes')
xlim([1 vStroke])

plot(maxDwellTime*ones(size(upDwellTime, 1), 1), 'r:')

figure, hold on;
plot(downDwellTime)
title('dwell time of down strokes')
xlim([1 vStroke])

plot(maxDwellTime*ones(size(downDwellTime, 1), 1), 'r:')

% save('dwell_time.mat', 'upDwellTime', 'downDwellTime', 'maxDwellTime', 'strokeTime')

