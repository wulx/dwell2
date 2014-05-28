function [exErrMap, errMap, mus] = beamdist(nMajors, nMinors, nsig)
%BEAMDIST ion beam density distribution

% ion beam density distribution along major axis

if nargin<3, nsig = 0.6; end

% uncerntainty: +/-5% (lambda)
% ref.: Zhou et al._2007_Amending the uniformity of ion beam current density profile

lambda = 0.1;
nRands = 40;

mus = 1.0 + lambda * (rand(nRands, 1) - 0.5);

coffs = 0.01 * (rand(nRands, 1) - 0.5);

nSigmas = nsig + 0.01 * (rand(nRands, 1) - 0.5);


%% ion beam density distribution along minor axis
nSamps = 8;

errMap = nan(nRands, nSamps);

% standard normal distribution: mu = 0, sigma = 1
stdNormDist = makedist('Normal');

for i = 1:nRands
    mu = mus(i);
    nSigma = nSigmas(i);
    coff = coffs(i); % center offset, range (-1, 1)
    
    % sampling points
    sampPts = linspace((coff-1)*nSigma, (coff+1)*nSigma, nSamps);
    
    samps = pdf(stdNormDist, sampPts);
    
    errMap(i, :) = mu/mean(samps) * samps;
end

[X, Y] = meshgrid(1:nSamps, 1:nRands);
[XI, YI] = meshgrid(linspace(1, nSamps, nMinors), linspace(1, nRands, nMajors));

exErrMap = interp2(X, Y, errMap, XI, YI, 'spline');

% figure, hold on;
% surf(exErrMap)
% shading flat
