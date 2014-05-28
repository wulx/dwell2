function [F, D, crtSteps] = datatable(params, Crt, strokeTime)
%DATATABLE a reduced copy of DWELLOPT

% copyright (c) wulx, gurdy.woo@gmail.com
% last modified by wulx, 2014/5/25


% #1 pre-processing ------------------------------------------------------%

% variable shortcuts
scaleDivs = Crt.scaleDivs;

crtLocs = Crt.locs;
crtTypes = Crt.types;
crtSteps = Crt.steps;

nRamps = Crt.num - 1;

w_a = params(1:nRamps);
w_d = params(nRamps+(1:nRamps));
w_f = params(2*nRamps+(1:nRamps));
a1 = params(end-1);
a2 = params(end);

% steps difference should not be less than 4
narrows = abs(crtSteps(2:end) - crtSteps(1:end-1)) < 4;
broadenMask = [narrows(1) narrows(1:end)] | [narrows(1:end) narrows(end)];

% broaden steps difference
amp1 = 10 * a1 * crtTypes; % @param --------------------------------------@
amp2 = 10 * a2 * broadenMask .* crtTypes; % @param -----------------------@

crtSteps = round(crtSteps + amp1 + amp2);

% #2 protocol of dwell time algorithm ------------------------------------%

% pre-allocation for steps and time sequencies
F = cell(1, nRamps);
D = cell(1, nRamps);

% handles to solve maximum frequency
sf = makeSolvefm;

t_diff = 0;
for j = 1:nRamps
    % total steps and total time
    sn_tot = crtSteps(j+1) - crtSteps(j);
    
    % go backward or forward
    goBackward = false;
    if sn_tot < 0
        sn_tot = -sn_tot;
        goBackward = true;
    end
    
    s_tot = scaleDivs(crtLocs(j+1)) - scaleDivs(crtLocs(j));
    stroke = scaleDivs(end) - scaleDivs(1);
    t_tot = t_diff + strokeTime * s_tot / stroke;
    
    % numbers of steps, [sn_a sn_c sn_d]
    sn_a = round(w_a(j) * sn_tot);
    if sn_a < 1
        sn_a = 1;
    elseif sn_a > sn_tot-2
        sn_a = sn_tot - 2;
    end
    
    % make sure not over range
    sn_d = round((1-w_a(j))*w_d(j) * sn_tot);
    if sn_d < 1
        sn_d = 1;
    elseif sn_d > sn_tot-2
        sn_d = sn_tot - 2;
    end
    
    sn_c = sn_tot - sn_a - sn_d;

    % initial frequency f_i and maximum frequency f_m
    % (0, (sn_tot/t_tot - 2)/(0.25*a_max*t_tot)]
    % f_i = max(floor(sn_tot/t_tot - 0.25*a_max*t_tot*w_f(j-1)), 2)
    f_i = round(w_f(j)*(sn_tot/t_tot - 2) + 2);

    %#TODO: maximum frequency may be no solutions.
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
    
    if goBackward
        F{j} = -f_list;
    else
        F{j} = f_list;
    end
    D{j} = dt_list;
    
end


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

