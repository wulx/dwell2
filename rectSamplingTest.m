%% Rectangular Sampling
close all; clear all; clc

load('data/ETCH_DEPTH.mat', 'depth')
load('data/ERROR_MODEL.mat', 'supDepth', 'optSupDepth', 'runSupDepth')

[m, n] = size(depth);
rectMask = false(m, n);

nRows = 50;
nCols = 50;

rowMod = mod(m, nRows-1);
colMod = mod(n, nCols-1);

if ~rowMod
    rowDist = [fix(m/(nRows-1))*ones(1,nRows-1) rowMod];
else
    rowDist = fix(m/nRows)*ones(1,nRows);
end

if ~colMod
    colDist = [fix(m/(nCols-1))*ones(1,nCols-1) rowMod];
else
    colDist = fix(m/nCols)*ones(1,nCols);
end

rectMaskCell = mat2cell(rectMask, rowDist, colDist);

for i = 5:10:50
    for j = 5:10:50
        rectMaskCell{i, j} = true(size(rectMaskCell{i, j}));
    end
end

for i2 = [1 25 50]
    for j2 = [1 25 50]
        rectMaskCell{i2, j2} = true(size(rectMaskCell{i2, j2}));
    end
end

revisedRectMask = cell2mat(rectMaskCell);

sampDepth = depth;
sampDepth(~revisedRectMask) = nan;

figure, imshow(mat2gray(depth), 'colormap', jet);

figure, imshow( mat2gray(revisedRectMask) );

figure, imshow(mat2gray(sampDepth), 'colormap', [jet; 0 0 0]);









