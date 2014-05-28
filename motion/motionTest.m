%% trapzoidal motion profile
% by wulx, last modified: 2014/5/23

%% clean and add path
close all; clear all; clc

% add path: ./GenTraj/
if isempty(strfind(path, 'GenTraj'))
    folderExist = exist('GenTraj', 'dir');
    if folderExist == 7 % folder
        oldpath = addpath(fullfile(pwd, 'GenTraj'));
    end
end

% add path: ./adjust_quiver_arrowhead_size/
if isempty(strfind(path, 'adjust_quiver_arrowhead_size'))
    folderExist = exist('adjust_quiver_arrowhead_size', 'dir');
    if folderExist == 7 % folder
        addpath(fullfile(pwd, 'adjust_quiver_arrowhead_size'));
    end
end

%%
% position, s
% velocity, v
% acceleration, a and deceleration, d
% jerk, j

% acceleration limited
% A = 1; % max acceleration
% V = 1.2; % max velocity
% P = 3; % total stroke
% [Y,T]=GenTraj(A,V,P);

% figure;
% sp(1)=subplot(3,1,1);plot(T,Y(3,:))
% sp(2)=subplot(3,1,2);plot(T,Y(2,:))
% sp(3)=subplot(3,1,3);plot(T,Y(1,:))
% linkaxes(sp,'x');
% xlim([0 T(end)]) % added by wulx
% ylabel(sp(1),'Position [m]');ylabel(sp(2),'Velocity [m/s]');ylabel(sp(3),'Acceleration [m/s^2]');xlabel(sp(3),'Time [s]')

%% rotation

[Y,T]=GenTraj(1, 0.8, 1.5);

figure;
sp(1)=subplot(3,1,1);plot(T,Y(3,:))
sp(2)=subplot(3,1,2);plot(T,Y(2,:))
sp(3)=subplot(3,1,3);plot(T,Y(1,:))
linkaxes(sp,'x');
xlim([0 T(end)]) % added by wulx
ylabel(sp(1),'Rotation angle [rad]');ylabel(sp(2),'Velocity [rad/s]');ylabel(sp(3),'Acceleration [rad/s^2]');xlabel(sp(3),'Time [s]')

leafWidth = 30;
projWidths = leafWidth * sin(Y(3,:));
strokeTime = T(end);

figure, plot(T, projWidths)
xlim([0 strokeTime])
xlabel('Time [s]')
ylabel('Projected width [mm]')


vScan = 60; % 60 mm/s
scaleDivs = linspace(0, vScan*strokeTime, 1000);
% margins = [leafWidth/2 leafWidth/2];

dwellTime = timecount(projWidths, strokeTime, scaleDivs);

figure, plot(scaleDivs, dwellTime);
xlim([0 scaleDivs(end)])
xlabel('stroke [mm]')
ylabel('dwell time [s]')


filter = [1:1000:(numel(T)-1) numel(T)];
thetas = Y(3, :);
thetas = thetas(filter);
rho = leafWidth/2;

[v, u] = pol2cart(thetas, rho);

figure, hold on;
% plot(thetas)

t = T(filter);
hq1 = quiver(t, zeros(size(filter)), u, v);
hq2 = quiver(t, zeros(size(filter)), -u, -v);
plot(t, zeros(size(filter)), 'ro');
plot(t, zeros(size(filter)), 'k-');
axis equal
axis off

adjust_quiver_arrowhead_size(hq1, 0.3);
adjust_quiver_arrowhead_size(hq2, 0.3);

% restore the old path
path(oldpath);
