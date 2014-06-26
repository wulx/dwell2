%% dwell time optimization demo
% last modified by wulx, 2014/5/24

clear all; close all; clc


% load optimization results
S2 = load('data/OPT_RESULT.mat', 'C', 'Rmses', 'optDwellTime');
C = S2.C;
Rmses = S2.Rmses;
optDwellTime = S2.optDwellTime;


ind = 1;

ogee = C(ind).ogee;
% spline interpolation (or extrapolation) and smoothing
[xx, ys] = spis(ogee);

optDwellTime1 = optDwellTime(:, ind);
rmse1 = Rmses{ind};

figure, hold on;
plot(xx, ys, 'LineWidth', 1.2, 'LineStyle', '--', 'Color', 'black')
plot(xx(2:end-1), ogee, 'LineWidth', 1.6, 'LineStyle', '-', 'Color', 'black')
plot(xx(2:end-1), optDwellTime1, 'LineWidth', 1.4, 'LineStyle', ':', 'Color', 'black');
title(['RMSE: ' num2str(rmse1)]);

nums = ogee(~isnan(ogee));
xs = xx([false; ~isnan(ogee); false]);
plot([xs(1) xs(1)], [0, nums(1)], 'LineWidth', 1, 'LineStyle', ':', 'Color', 'black')
plot([xs(end) xs(end)], [0, nums(end)], 'LineWidth', 1, 'LineStyle', ':', 'Color', 'black')

% mark all critical points at the spline
for ni = 1:C(ind).num
    li = C(ind).locs(ni);
    ti = C(ind).types(ni);
    if ti > 0 % peaks
        plot(xx(li), ys(li), 'Color', 'black', 'LineWidth', 2, 'Marker', '^', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'black')
    elseif ti < 0 % valleys
        plot(xx(li), ys(li), 'Color', 'black', 'LineWidth', 2, 'Marker', 'v', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'black')
    else % inflection points
        plot(xx(li), ys(li), 'Color', 'black', 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'black')
    end
end

axis tight





