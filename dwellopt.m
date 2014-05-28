function [r, dwellTime] = dwellopt(w_a, w_d, w_f, a1, a2, Crt, strokeTime, stepAngleDeg, leafWidth, timeStep)
%DWELLOPT motion profile optimization to obtain the minimum RMSE of dwell time
% varargin:
%   @params to be optimized
%   w_a  --  weight factors of number of acceleration steps, (0, 1)
%   w_d  --  weight factors of number of deceleration steps, (0, 1)
%   w_f  --  weight factors of initial frequency, (0, 1)
%   a1   --  amplitude of enlarging peaks or valleys
%   a2   --  amplitude of braodening narrow forks
%   @params to set
%   ogee --  initial dwell time data
%   stroke  --  distance in mm
%   stepAngleDeg  --  step angle in degree
%   a_max  --  maximum acceleration
%   timeStep  --  time step
% varargout:
%   r           --  RMSE of dwell time
%   dwellTime   --  dwell time
%


% copyright (c) wulx, gurdy.woo@gmail.com
% last modified by wulx, 2013/11/21, 2014/5/27


% #1 pre-processing ------------------------------------------------------%

% variable shortcuts
ogee = Crt.ogee;
scaleDivs = Crt.scaleDivs;

crtLocs = Crt.locs;
crtTypes = Crt.types;
crtSteps = Crt.steps;

nCrts = Crt.num;

% steps difference should not be less than 4
narrows = abs(crtSteps(2:end) - crtSteps(1:end-1)) < 4;
broadenMask = [narrows(1) narrows(1:end)] | [narrows(1:end) narrows(end)];

% broaden steps difference
amp1 = 10 * a1 * crtTypes; % @param -------------------------------------------@
amp2 = 10 * a2 * broadenMask .* crtTypes; % @param ----------------------------@

crtSteps = round(crtSteps + amp1 + amp2);

% #2 protocol of dwell time algorithm ------------------------------------%

% pre-allocation for steps and time sequencies
S = cell(1, nCrts-1);
T = cell(1, nCrts-1);

% handles to solve maximum frequency
sf = makeSolvefm;

t_diff = 0;
for j = 2:nCrts
    % total steps and total time
    sn_tot = crtSteps(j) - crtSteps(j-1);
    
    % go backward or forward
    goBackward = false;
    if sn_tot < 0
        sn_tot = -sn_tot;
        goBackward = true;
    end
    
    s_tot = scaleDivs(crtLocs(j)) - scaleDivs(crtLocs(j-1));
    stroke = scaleDivs(end) - scaleDivs(1);
    t_tot = t_diff + strokeTime * s_tot / stroke;
    
    % numbers of steps, [sn_a sn_c sn_d]
    sn_a = round(w_a(j-1) * sn_tot);
    if sn_a < 1
        sn_a = 1;
    elseif sn_a > sn_tot-2
        sn_a = sn_tot - 2;
    end
    
    % make sure not over range
    sn_d = round((1-w_a(j-1))*w_d(j-1) * sn_tot);
    if sn_d < 1
        sn_d = 1;
    elseif sn_d > sn_tot-2
        sn_d = sn_tot - 2;
    end
    
    sn_c = sn_tot - sn_a - sn_d;

    % initial frequency f_i and maximum frequency f_m
    % (0, (sn_tot/t_tot - 2)/(0.25*a_max*t_tot)]
    % f_i = max(floor(sn_tot/t_tot - 0.25*a_max*t_tot*w_f(j-1)), 2)
    f_i = round(w_f(j-1)*(sn_tot/t_tot - 2) + 2);

    %# TODO: maximum frequency may be no solutions.
    % this is not exact solution, but is almost equals to it and greater than it.
    f_m = sf(sn_a, sn_c, sn_d, f_i, t_tot) + 1; % +1

    sn = [sn_a, sn_c, sn_d];
    pf = [f_i, f_m];
    method = 'round';
    [f_list, dt_list] = time_per_step(sn, pf, 1, method);

    while sum(dt_list(:)) > t_tot
        f_m = f_m + 1;
        [f_list, dt_list] = time_per_step(sn, [f_i, f_m], 1, method);
    end
    
    % time difference, positive!
    t_diff = t_tot - sum(dt_list(:));
    
    % time sequencies of every step
    timeSeqs = steptime(f_list, dt_list);
    
    %# tricks: stepper from 0
    steps = 0:numel(timeSeqs);
    
    if goBackward
        steps = crtSteps(j-1) - steps;
    else
        steps = crtSteps(j-1) + steps;
    end
    
    
    if j < nCrts
        S{j-1} = steps(1:end-1);
        T{j-1} = timeSeqs;
    else
        S{j-1} = [steps, steps(end)];
        T{j-1} = [timeSeqs, t_diff]; % time compensation
    end
    
end

% concatenant
T = cell2mat(T);
S = cell2mat(S);

% timeStep = 0.001; % 1 ms
nScan = ceil(strokeTime / timeStep);

stepsamp = zeros(1, nScan);
timeline = linspace(0, strokeTime, nScan);

nStep = numel(T);
for i = 1:nStep
    stepsamp(timeline>=sum(T(1:i-1)) & timeline<=sum(T(1:i))) = S(i);
end

projWidths = step2width(stepsamp, stepAngleDeg, leafWidth);

dwellTime = timecount(projWidths, strokeTime, scaleDivs);


% #3 root-mean-square deviation of dwell time

% filter NaNs
% isNum = ~isnan(ogee);
% ogee = ogee(isNum)';
dwellTime = dwellTime(2:end-1)';
dwellTime(isnan(ogee)) = nan;

r = rmse(ogee, dwellTime);



% solve maximum frequency ------------------------------------------------%
%# approximate solution
    function sf = makeSolvefm
        fm = sym('fm', 'positive');
        
        function f = solvefm(sn_a, sn_c, sn_d, fi, t)
            f = round(double(solve(2*(sn_a+sn_d)/(fi+fm) + sn_c/fm == t, fm)));
        end
        
        sf = @solvefm;
    end
        
end

