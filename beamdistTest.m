%% test on ion beam density distribution
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
ionBeamWidth = 60; % mm @ionBeamWidth ------------------------------------@
leafWidth = 60; % mm
subHeight = 400; % mm
subLength = 400; % mm
vStroke = leafWidth + subHeight; % vertical stroke
nSig = ionBeamWidth / 60 * 0.6;

nRands = 40;

[exErrMap, errMap, mus] = beamdist(vStroke, ionBeamWidth, nSig);
xx = linspace(1, nRands, 400);
yy = spline(1:nRands, mus, xx);

figure, hold on;
xlim([0 41])
ylim([0 1.1])
plot(xx, yy, 'k-')
stem(mus, 'k:')

plot(xlim, [0.95 0.95], 'r:')
plot(xlim, [1.05, 1.05], 'r:')

figure, hold on;
surf(errMap)

figure, hold on;
meshz(exErrMap)
shading flat

%% minor axis distribution
nSamps = 8;

% standard normal distribution: mu = 0, sigma = 1
stdNormDist = makedist('Normal');

coff = 0.01 * (rand - 0.5);

nSigma = nSig + 0.01 * (rand - 0.5);
    
% sampling points
sampPts = linspace((coff-1)*nSigma, (coff+1)*nSigma, nSamps);

samps = pdf(stdNormDist, sampPts);
samps = samps / mean(samps);

xx = linspace(1, nSamps, ionBeamWidth);
yy = spline(1:nSamps, samps, xx);

figure, plot(1:ionBeamWidth, yy, 'k-')
box on
title(['center offset: ' num2str(coff) ' n-fold sigma: ' num2str(nSigma) ' mean value: 1'])

