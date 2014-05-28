%% dwell time optimization demo
% last modified by wulx, 2014/5/24

clear all; close all; clc
% add lsramp path
if isempty(strfind(path, 'lsramp'))
    folderExist = exist('lsramp', 'dir');
    if folderExist == 7 % folder
        oldpath = addpath(fullfile(pwd, 'lsramp'));
    end
end

% Open pool of MATLAB sessions for parallelcomputation
isOpen = matlabpool('size') > 0;
if ~isOpen
    % fix bug, convert str to double firstly
    nProcessers = str2double(getenv('NUMBER_OF_PROCESSORS'));
    if nProcessers > 3
        matlabpool('local', nProcessers-1);
    else
        matlabpool('local', nProcessers);
    end
end

% settings
stepAngleDeg = 1.8 / 8;
opts = {'MINPEAKHEIGHT', -1, 'MINPEAKDISTANCE', 15, 'NPEAKS', 5};

% struct to save data of critical points
C = struct('ogee', {}, ...
    'scaleDivs', {}, ...
    'locs', {}, ...
    'types', {}, ...
    'num', {}, ...
    'steps', {}, ...
    'strkDir', {});

% load dwell time data
D = load('data/DWELL_TIME.mat');
dwellTime = D.dwellTime;
strokeTime = D.strokeTime;
maxDwellTime = D.maxDwellTime;
strkSet = D.strkSet;

nTiers = numel(strkSet);
nStrks = sum([strkSet.nStrks]);

headsAndTails = nan(2, nStrks);

for i = 1:nTiers
    for ind = strkSet(i).indices
        ogee = dwellTime(:, ind);
        
        % reverse processing directions for down stroks
        if strcmp(strkSet(i).mode, 'DOWN')
        	ogee = ogee(end:-1:1);
        end
        
        C(ind).ogee = ogee;
        
        % spline interpolation (or extrapolation) and smoothing
        [xx, ys] = spis(ogee);
        C(ind).scaleDivs = xx;
        
        % low pass filtering
        lowPass = ys>maxDwellTime;
        if any(lowPass)
            ys(lowPass) = maxDwellTime;
        end
        
        % find critical points
        [locs, ptype] = findpoints(ys, opts);
        
        C(ind).locs = [1 locs numel(xx)]; % all critical locations
        C(ind).types = [0 ptype 0]; % treate boundary points as infletion points
        C(ind).num = numel( C(ind).locs );
        
        crtDwellTimes = ys(C(ind).locs);
        crtDegs = rad2deg( asin(crtDwellTimes / maxDwellTime) );
        crtSteps = round(crtDegs / stepAngleDeg);
        
        headsAndTails(:, ind) = crtSteps([1 end]);
        
        C(ind).steps = crtSteps;
        
        figure, hold on;
        plot(xx, ys, 'k-')
        plot(xx(2:end-1), ogee, 'k:')
        
        % mark all critical points at the spline
        for ni = 1:C(ind).num
            li = C(ind).locs(ni);
            ti = C(ind).types(ni);
            if ti > 0 % peaks
                plot(xx(li), ys(li), 'ro')
            elseif ti < 0 % valleys
                plot(xx(li), ys(li), 'bo')
            else % inflection points
                plot(xx(li), ys(li), 'ko')
            end
        end
        
    end
end

% average and splice heads and tails
headsAndTails(1, 2:end) = round((headsAndTails(1, 2:end) + headsAndTails(2, 1:end-1)) / 2);
headsAndTails(2, 1:end-1) = headsAndTails(1, 2:end);

%%
Params = cell(1, nStrks);
Rmses = cell(1, nStrks);
optDwellTime = nan( size(dwellTime) );

leafWidth = 60; % 60 mm
% a_max = 60000; % benchmark test result: 66333 pulse/s^2
timeStep = 0.001; % 1ms

psopts = psoptimset('Cache', 'on', 'Vectorized','off', 'MaxIter', 20, ...
    'UseParallel', 'always', 'CompletePoll', 'on', 'TolFun', 0.03, ...
    'TolX', 0.05, 'PollMethod', 'GPSPositiveBasis2N', ...
    'MaxMeshSize', 0.2, 'InitialMeshSize', 0.1, 'MeshAccelerator', 'on', ...
    'Display', 'iter');

for j = 1:nTiers
    for k = strkSet(j).indices
        % revise the heads and tails
        C(k).steps([1 end]) = headsAndTails(:, k);
        
        nw = C(k).num - 1;
        
        params = [1/3*ones(1,2*nw), 0.4*ones(1,nw), 0.3, 0.3];
        
        % [a, b, c] = dwopt(params);
        lb = [zeros(1,3*nw), 0.1, 0.1];   % Lower bounds
        ub = [ones(1,3*nw), 0.9, 0.9];  % Upper bounds
        
        % params = [w_a, w_d, w_f, a1, a2]
        dwopt = @(params) dwellopt(params(1:nw), params(nw+(1:nw)), params(2*nw+(1:nw)), params(end-1), params(end), ...
            C(k), strokeTime, stepAngleDeg, leafWidth, timeStep);
        %     'PlotFcns', {@psplotbestf,@psplotmeshsize,@psplotfuncount,@psplotbestx}, ...
        %     'PlotInterval', 1);
        
        % pattern search parallely
        [params1, rmse1] = patternsearch(dwopt, params, [], [], [], [], lb, ub, psopts);
        
        ogee = C(k).ogee;
        xx = C(k).scaleDivs;
        % plot results
        [r, dwellTime] = dwopt(params1);
        isNum = ~isnan(ogee);
        ogee = ogee(isNum)';
        xx = xx([false isNum' false]);
        
        figure, hold on;
        plot(xx, ogee, 'k-');
        plot(xx, dwellTime, 'b-');
        title(['RMSE: ' num2str(rmse1)]);
        
        Params{k} = params1;
        Rmses{k} = rmse1;
        
        if strcmp(strkSet(j).mode, 'DOWN')
            optDwellTime(:, k) = dwellTime(end:-1:1);
        else
            optDwellTime(:, k) = dwellTime;
        end
    end
end

if isOpen
    matlabpool close
end

% restore the old path
path(oldpath);

%save('data/OPT_RESULT.mat', 'C', 'Params', 'Rmses', 'nStrks', 'optDwellTime')
