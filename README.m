%% dwell2 readme doc
close all; clear all; clc

%% generate free-form optical surface randomly
% based on the following code:
% [f,x,y] = rsgeng2D (N,rL,h,clx,cly) 
% DESCRIPTION: Random rough surface generator of two-dimensional (isotropic and non-isotropic) square surfaces with Gaussian hdf and Gaussian acf.
% INPUT: N-number of surface points, rL-length of surface side, h-rms height, clx-correlation length in x, cly-correlation length in y
% OUTPUT: f-surface heights, x-surface points, y-surface points
% credit: http://www.mysimlabs.com/surface_generation.html

% [height, x, y] = rsgeng2D(400,40,0.4,20,12);
% save('ROUGH_SURF.mat', 'height')


%% conversions
% ROUGH_SURF.mat ==> ETCH_DEPTH.mat
% height ==> etch depth


%% rastering
% test on rastering and stitching method
rasterTest


%% ion beam current density distribution
% Description of second code block

beamdistTest

%%

