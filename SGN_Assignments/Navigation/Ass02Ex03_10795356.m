% Spacecraft Guidance and Navigation (2024/2025)
% Assignment #2: Exercise #3
% Author: Luca Frassinella


% --------------------- EX. 3 - SEQUENTIAL FILTERS ---------------------- %

clearvars; close all; clc;
cspice_kclear()

plotSettings;
addpath('.\kernels\')
addpath('.\sgp4\')
addpath('.\tle\')
addpath('.\mice\src\mice')
addpath('.\mice\lib')
cspice_furnsh('assignment02.tm')

% Load dataset:
[data, constants, settings] = loadSet;

% extract data:
lander = data.lander;
x0 = data.x0;
et0 = data.et0;
etf = data.etf;

%% 1. Check Visibility Window:

% Propagate trajectory in MCI frame (Keplerian motion):
data.tspan = et0 : data.timeStep : etf;
[~, xx] = ode113(@TBP, data.tspan, x0, settings.odeOpt, constants);

% Compute orbiter local coordinates with respect to lunar lander
% (using latitudinal coordinates for lander's state):
[orbiter_lc] = orbLocalCoords(xx, lander, data, 'latrec');

% Check that the lander and orbiter are in relative visibility during the
% entire time interval:
[orbiter_visibleCoords, errFlag] = checkVisibility(orbiter_lc, lander);

% Plot of relative visibility:
plotVisEl(data, orbiter_lc)

% Plot of Lunar Orbiter Trajectory:
figure()
plot3(xx(:, 1), xx(:, 2), xx(:, 3), '-k')
hold on
plot3(xx(1,1), xx(1,2), xx(1,3), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'b')
plot3(xx(end,1), xx(end,2), xx(end,3), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
rMoon = cspice_bodvrd('MOON','RADII',3);
rMoon = rMoon(1);
[X, Y, Z] = sphere(50);
hSurface = surf(X*rMoon, Y*rMoon, Z*rMoon, 'EdgeColor', 'k', 'FaceColor', [0.5 0.5 0.5]);
legend('Trajectory', 'Initial Position', 'Final Position', 'Moon', 'FontSize', 30)
axis equal;
xlabel('x [km]', 'FontSize', 35);
ylabel('y [km]', 'FontSize', 35);
zlabel('z [km]', 'FontSize', 35);
title('Lunar Orbiter Trajectory (@Moon MCI)')


%% 2. Simulate Measurements:

% Simulate the measurements using kernels:
[orbiter_lc_kernels] = orbLocalCoords(xx, lander, data, 'kernels');

% Add the measurements noise:
idealMeasurements.position = xx(:, 1:3);
idealMeasurements.range = orbiter_lc_kernels.range;

% Compute real simulated measurements:
realMeasurements.position = mvnrnd(idealMeasurements.position, data.Rpos);
realMeasurements.range = mvnrnd(idealMeasurements.range', data.sigmaMeas^2);


%%%% PLOTS:
% Plot Ideal position and relative range over time:
figure()
hold on
plot(data.tspan/cspice_spd, idealMeasurements.range, 'DisplayName','Range [km]');
plot(data.tspan/cspice_spd, xx(:,1), 'DisplayName','x Coordinate [km]');
plot(data.tspan/cspice_spd, xx(:,2), 'DisplayName','y Coordinate [km]');
plot(data.tspan/cspice_spd, xx(:,3), 'DisplayName','z Coordinate [km]');
ylabel('Measuremnts [km]', 'FontSize', 35)
xlabel('2024-NOV-18', 'FontSize', 35)
legend('FontSize', 30);
xlim([data.tspan(1) data.tspan(end)]/cspice_spd)
num_ticks = 5;
ax = gca;
ax.FontSize = 25;
tick_indices = round(linspace(1, length(data.tspan), num_ticks)); 
tick_values = data.tspan(tick_indices) / cspice_spd();
tick_labels = cell(num_ticks, 1); 
for i = 1:num_ticks
    utc_full = cspice_et2utc(data.tspan(tick_indices(i)), 'C', 0); 
    tick_labels{i} = utc_full(12:end); 
end
xticks(tick_values);
xticklabels(tick_labels);
xtickangle(0);

% Plot of the error with 3sigma boundaries:
figure()
% Subplot 1: Range error
subplot(4, 1, 1);
plot(data.tspan/cspice_spd, realMeasurements.range - idealMeasurements.range', '-k', 'LineWidth', 1.5);
hold on
yline(3*data.sigmaMeas, '--r', 'LineWidth', 1.5);
yline(-3*data.sigmaMeas, '--r', 'LineWidth', 1.5);
title('Range', 'FontSize', 25);
xlim([data.tspan(1) data.tspan(end)]/cspice_spd);
% Subplot 2: X error
subplot(4, 1, 2);
plot(data.tspan/cspice_spd, realMeasurements.position(:, 1) - xx(:, 1), '-k', 'LineWidth', 1.5);
hold on
yline(3*data.sigmaMeas, '--r', 'LineWidth', 1.5);
yline(-3*data.sigmaMeas, '--r', 'LineWidth', 1.5);
title('x component', 'FontSize', 25);
xlim([data.tspan(1) data.tspan(end)]/cspice_spd);
% Subplot 3: Y error
subplot(4, 1, 3);
plot(data.tspan/cspice_spd, realMeasurements.position(:, 2) - xx(:, 2), '-k', 'LineWidth', 1.5);
hold on
yline(3*data.sigmaMeas, '--r', 'LineWidth', 1.5);
yline(-3*data.sigmaMeas, '--r', 'LineWidth', 1.5);
title('y component', 'FontSize', 25);
xlim([data.tspan(1) data.tspan(end)]/cspice_spd);
% Subplot 4: Z error
subplot(4, 1, 4);
plot(data.tspan/cspice_spd, realMeasurements.position(:, 3) - xx(:, 3), '-k', 'LineWidth', 1.5);
hold on
yline(3*data.sigmaMeas, '--r', 'LineWidth', 1.5);
yline(-3*data.sigmaMeas, '--r', 'LineWidth', 1.5);
title('z', 'Interpreter', 'latex', 'FontSize', 25);
xlim([data.tspan(1) data.tspan(end)]/cspice_spd);
% Adjust x-axis:
subplot(4, 1, 4);
xlim([data.tspan(1) data.tspan(end)]/cspice_spd);
num_ticks = 5;
tick_indices = round(linspace(1, length(data.tspan), num_ticks)); 
tick_values = data.tspan(tick_indices) / cspice_spd();
tick_labels = cell(num_ticks, 1); 
for i = 1:num_ticks
    utc_full = cspice_et2utc(data.tspan(tick_indices(i)), 'C', 0); 
    tick_labels{i} = utc_full(12:end); 
end
xticks(tick_values);
xticklabels(tick_labels);
xtickangle(0);
set(gca, 'FontSize', 25)
% Remove x-ticks from all but the last subplot
for i = 1:3
    subplot(4, 1, i);
    set(gca, 'XTickLabel', []);
end
legend('Error', '$\pm 3\sigma$', 'FontSize', 25, 'Location', 'northeast');
han = axes(gcf, 'Visible', 'off'); 
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
xlabel(han, '2024-NOV-18', 'FontSize', 30);
ylabel(han, 'Errors [km]', 'FontSize', 30);

%% 3. Estimate the lunar orbiter absolute state:

% Mean and Covariance:
data.cov0 = data.P0(1:6, 1:6);
data.mean0 = mvnrnd(data.x0, data.cov0)';
data.n = length(data.x0);

% Initialize Matrices for estimated mean and covariance:
mean_mat = zeros(data.n, length(data.tspan));
cov_mat = zeros(data.n, data.n, length(data.tspan));

mean_mat(:, 1) = data.mean0;
cov_mat(:, :, 1) = data.cov0;

for ii = 2 : length(data.tspan)
    [mean_mat(:, ii), cov_mat(:, :, ii)] = UKF(mean_mat(:, ii-1), cov_mat(:, :, ii-1), ...
                   data.tspan(ii-1), data.tspan(ii), realMeasurements.position(ii, :), data, settings, constants);
end

%% Store results at final time:
results.ex3.mean = mean_mat(:, end);
results.ex3.cov = cov_mat(:, :, end);

% Errors:
errPos = zeros(1,length(data.tspan));
errVel = zeros(1,length(data.tspan));
stdPos = zeros(1,length(data.tspan));
stdVel = zeros(1,length(data.tspan));
for i = 1:length(data.tspan)
    errPos(i) = sqrt((xx(i,1) - mean_mat(1,i))^2 + (xx(i,2) - mean_mat(2,i))^2 +(xx(i,3) - mean_mat(3,i))^2);
    errVel(i) = sqrt((xx(i,4) - mean_mat(4,i))^2 + (xx(i,5) - mean_mat(5,i))^2 +(xx(i,6) - mean_mat(6,i))^2 );
    stdPos(i) = 3 * sqrt(trace(cov_mat(1:3,1:3, i)));
    stdVel(i) = 3 * sqrt(trace(cov_mat(4:6,4:6, i)));
end

% Plot errors:

subplot(1, 2, 1); 
semilogy(data.tspan/cspice_spd, errPos, 'k-', 'LineWidth', 1.8, 'DisplayName', 'Error'); 
hold on;
semilogy(data.tspan/cspice_spd, stdPos, 'r--', 'LineWidth', 1.5, 'DisplayName', '3$\sigma$');
xlabel('2024-NOV-18', 'FontSize', 40);
ylabel('Position Error [km]', 'FontSize', 40);
xlim([data.tspan(1)/cspice_spd, data.tspan(end)/cspice_spd])
% title('Position Error and $3\sigma$');
legend('FontSize', 25, 'Location', 'northeast');
grid on;
set(gca, 'FontSize', 20)
num_ticks = 5;
tick_indices = round(linspace(1, length(data.tspan), num_ticks)); 
tick_values = data.tspan(tick_indices) / cspice_spd();
tick_labels = cell(num_ticks, 1); 
for i = 1:num_ticks
    utc_full = cspice_et2utc(data.tspan(tick_indices(i)), 'C', 0); 
    tick_labels{i} = utc_full(12:end); 
end
xticks(tick_values);
xticklabels(tick_labels);
xtickangle(0);

% Subplot 2: Velocity Error
subplot(1, 2, 2);
semilogy(data.tspan/cspice_spd, errVel, 'k-', 'LineWidth', 1.8, 'DisplayName', 'Error'); 
hold on;
semilogy(data.tspan/cspice_spd, stdVel, 'r--', 'LineWidth', 1.5, 'DisplayName', '3$\sigma$');
xlabel('2024-NOV-18', 'FontSize', 35);
ylabel('Velocity Error [km]', 'FontSize', 35);
xlim([data.tspan(1)/cspice_spd, data.tspan(end)/cspice_spd])
% title('Velocity Error and $3\sigma$');
legend('FontSize', 25, 'Location', 'northeast');
grid on;
set(gca, 'FontSize', 20)
num_ticks = 5;
tick_indices = round(linspace(1, length(data.tspan), num_ticks)); 
tick_values = data.tspan(tick_indices) / cspice_spd();
tick_labels = cell(num_ticks, 1); 
for i = 1:num_ticks
    utc_full = cspice_et2utc(data.tspan(tick_indices(i)), 'C', 0); 
    tick_labels{i} = utc_full(12:end); 
end
xticks(tick_values);
xticklabels(tick_labels);
xtickangle(0);


%% 4. Estimate the lunar lander coordinates:

% Define augmented initial state, mean and covariance:
data.x0_aug = [data.x0; deg2rad(data.lander.lat); deg2rad(data.lander.lon)];
data.cov0_aug = data.P0;
data.mean0_aug = mvnrnd(data.x0_aug, data.cov0_aug);
data.n = length(data.x0_aug);

% Initialize Matrices for estimated mean and covariance:
augMean_mat = zeros(data.n, length(data.tspan));
augCov_mat = zeros(data.n, data.n, length(data.tspan));

augMean_mat(:, 1) = data.mean0_aug;
augCov_mat(:, :, 1) = data.cov0_aug;

for ii = 2 : length(data.tspan)
    [augMean_mat(:, ii), augCov_mat(:, :, ii)] = UKF(augMean_mat(:, ii-1), augCov_mat(:, :, ii-1), ...
                   data.tspan(ii-1), data.tspan(ii), [realMeasurements.position(ii, :), realMeasurements.range(ii)], data, settings, constants);
end

% Compute lander latitude and longitude from kernels:
stateLander = cspice_spkezr(lander.name, data.tspan, 'IAU_MOON', 'NONE', 'MOON');
[alt, real_lon, real_lat] = cspice_reclat(stateLander(1:3, :));

%% Store results at final time:
results.ex4.mean = augMean_mat(:, end);
results.ex4.cov = augCov_mat(:, :, end);
 
% Errors:
errPos = zeros(1,length(data.tspan));
errVel = zeros(1,length(data.tspan));
stdPos = zeros(1,length(data.tspan));
stdVel = zeros(1,length(data.tspan));
errLat = zeros(1, length(data.tspan));
errLon = zeros(1, length(data.tspan));
stdLat = zeros(1, length(data.tspan));
stdLon = zeros(1, length(data.tspan));
for i = 1:length(data.tspan)
    errPos(i) = sqrt((xx(i,1) - mean_mat(1,i))^2 + (xx(i,2) - mean_mat(2,i))^2 +(xx(i,3) - mean_mat(3,i))^2);
    errVel(i) = sqrt((xx(i,4) - mean_mat(4,i))^2 + (xx(i,5) - mean_mat(5,i))^2 +(xx(i,6) - mean_mat(6,i))^2 );
    stdPos(i) = 3 * sqrt(trace(cov_mat(1:3,1:3, i)));
    stdVel(i) = 3 * sqrt(trace(cov_mat(4:6,4:6, i)));
    stdLat(i) = 3 * sqrt(augCov_mat(7,7,i)) * cspice_dpr;
    stdLon(i) = 3 * sqrt(augCov_mat(8,8,i)) * cspice_dpr;
    errLat(i) = abs(real_lat(i)*cspice_dpr - augMean_mat(7,i)*cspice_dpr);
    errLon(i) = abs(real_lon(i)*cspice_dpr - augMean_mat(8,i)*cspice_dpr);
end

% Plot errors:

% Position & Velocity:
figure()
subplot(1, 2, 1); 
semilogy(data.tspan/cspice_spd, errPos, 'k-', 'LineWidth', 1.8, 'DisplayName', 'Error'); 
hold on;
semilogy(data.tspan/cspice_spd, stdPos, 'r--', 'LineWidth', 1.5, 'DisplayName', '3$\sigma$');
xlabel('2024-NOV-18', 'FontSize', 40);
ylabel('Position Error [km]', 'FontSize', 40);
xlim([data.tspan(1)/cspice_spd, data.tspan(end)/cspice_spd])
% title('Position Error and $3\sigma$');
legend('FontSize', 25, 'Location', 'northeast');
grid on;
set(gca, 'FontSize', 20)
num_ticks = 5;
tick_indices = round(linspace(1, length(data.tspan), num_ticks)); 
tick_values = data.tspan(tick_indices) / cspice_spd();
tick_labels = cell(num_ticks, 1); 
for i = 1:num_ticks
    utc_full = cspice_et2utc(data.tspan(tick_indices(i)), 'C', 0); 
    tick_labels{i} = utc_full(12:end); 
end
xticks(tick_values);
xticklabels(tick_labels);
xtickangle(0);
% Subplot 2: Velocity Error
subplot(1, 2, 2);
semilogy(data.tspan/cspice_spd, errVel, 'k-', 'LineWidth', 1.8, 'DisplayName', 'Error'); 
hold on;
semilogy(data.tspan/cspice_spd, stdVel, 'r--', 'LineWidth', 1.5, 'DisplayName', '3$\sigma$');
xlabel('2024-NOV-18', 'FontSize', 35);
ylabel('Velocity Error [km]', 'FontSize', 35);
xlim([data.tspan(1)/cspice_spd, data.tspan(end)/cspice_spd])
% title('Velocity Error and $3\sigma$');
legend('FontSize', 25, 'Location', 'northeast');
grid on;
set(gca, 'FontSize', 20)
num_ticks = 5;
tick_indices = round(linspace(1, length(data.tspan), num_ticks)); 
tick_values = data.tspan(tick_indices) / cspice_spd();
tick_labels = cell(num_ticks, 1); 
for i = 1:num_ticks
    utc_full = cspice_et2utc(data.tspan(tick_indices(i)), 'C', 0); 
    tick_labels{i} = utc_full(12:end); 
end
xticks(tick_values);
xticklabels(tick_labels);
xtickangle(0);

%%

% Latitude and Longitude:
% Plot errors:
figure()
subplot(1, 2, 1); 
plot(data.tspan/cspice_spd, errLat, 'k-', 'LineWidth', 1.8, 'DisplayName', 'Error'); 
hold on;
plot(data.tspan/cspice_spd, stdLat, 'r--', 'LineWidth', 1.5, 'DisplayName', '3$\sigma$');
xlabel('2024-NOV-18', 'FontSize', 40);
ylabel('Latitude Error [deg]', 'FontSize', 40);
xlim([data.tspan(1)/cspice_spd, data.tspan(end)/cspice_spd])
% title('Latitude Error and $3\sigma$');
legend('FontSize', 25, 'Location', 'northeast');
grid on;
set(gca, 'FontSize', 20)
num_ticks = 5;
tick_indices = round(linspace(1, length(data.tspan), num_ticks)); 
tick_values = data.tspan(tick_indices) / cspice_spd();
tick_labels = cell(num_ticks, 1); 
for i = 1:num_ticks
    utc_full = cspice_et2utc(data.tspan(tick_indices(i)), 'C', 0); 
    tick_labels{i} = utc_full(12:end); 
end
xticks(tick_values);
xticklabels(tick_labels);
xtickangle(0);

% Subplot 2: Longitude Error
subplot(1, 2, 2);
plot(data.tspan/cspice_spd, errLon, 'k-', 'LineWidth', 1.8, 'DisplayName', 'Error'); 
hold on;
plot(data.tspan/cspice_spd, stdLon, 'r--', 'LineWidth', 1.5, 'DisplayName', '3$\sigma$');
xlabel('2024-NOV-18', 'FontSize', 35);
ylabel('Longitude Error [deg]', 'FontSize', 35);
xlim([data.tspan(1)/cspice_spd, data.tspan(end)/cspice_spd])
% title('Longitude Error and $3\sigma$');
legend('FontSize', 25, 'Location', 'northeast');
grid on;
set(gca, 'FontSize', 20)
num_ticks = 5;
tick_indices = round(linspace(1, length(data.tspan), num_ticks)); 
tick_values = data.tspan(tick_indices) / cspice_spd();
tick_labels = cell(num_ticks, 1); 
for i = 1:num_ticks
    utc_full = cspice_et2utc(data.tspan(tick_indices(i)), 'C', 0); 
    tick_labels{i} = utc_full(12:end); 
end
xticks(tick_values);
xticklabels(tick_labels);
xtickangle(0);

%% Functions:

function [data, constants, settings] = loadSet()

settings.arcsec2rad = pi/(180*3600);
settings.typerun    = 'u';  % user-provided inputs to SGP4 Matlab function
settings.opsmode    = 'a';  % afspc approach ('air force space command')
settings.whichconst =  72;  % WGS72 constants (radius, gravitational parameter)
settings.odeOpt = odeset('reltol', 1e-12, 'abstol', 1e-12);

% Initial Mean State
data.r0 = [4307.844185282820,-1317.980749248651, 2109.210101634011];
data.v0 = [-0.110997301537882,-0.509392750828585, 0.815198807994189];
data.x0 = [data.r0, data.v0]';

% Initialize data in UTC
data.t0 = '2024-11-18T16:30:00.000';
data.tf = '2024-11-18T20:30:00.000';

% Epoch
data.et0=cspice_str2et(data.t0);
data.etf=cspice_str2et(data.tf);
data.timeStep = 30; % [s]

% Define the measurement noise
data.sigmaMeas = 0.1; %[km]
% Position Noise Covariance matrix:
data.Rpos = diag([data.sigmaMeas^2, data.sigmaMeas^2, data.sigmaMeas^2]);
data.R_aug = diag([data.sigmaMeas^2, data.sigmaMeas^2, data.sigmaMeas^2, data.sigmaMeas^2]);

% Define the initial covatiance
data.P0 = diag([10,1,1,0.001,0.001,0.001,0.00001,0.00001]);

% Define the lander features
data.lander.name = 'MOONLANDER';
data.lander.lat = 78;
data.lander.lon = 15;
data.lander.alt = 0; % [m]
data.lander.minEl = 0; % [deg]

% Constants (Moon data):
constants.mu_moon = cspice_bodvrd('MOON', 'GM',1);
rMoon = cspice_bodvrd('MOON','RADII',3);
constants.rMoon = rMoon(1);

% Lander radius:
data.lander.r = constants.rMoon;

% Unscented Transform Parameters:
data.ut.alpha = 0.01;
data.ut.beta = 2;

end

function [dxx] = TBP(~, xx, constants)
% -------------------------------------------------------------------------
% TBP - RHS of equations of motion for the Two-Body Problem (TBP)
%
%   Inputs:
%       - ~: Time variable (not used).
%       - xx: State vector [6x1]
%       - constants.mu: Gravitational parameter.
%
%   Output:
%       - dxx: Derivative of the state vector [6x1]
% -------------------------------------------------------------------------

x = xx(1);
y = xx(2);
z = xx(3);
r_norm = norm([x;y;z]);

mu = constants.mu_moon;

% Equations of motion
dxx(1:3) = xx(4:6);
dxx(4:6) = -mu/r_norm^3 * [x;y;z];

% Transpose to have a column vector:
dxx = dxx';
    
end

function [localCoords] = orbLocalCoords(xx, lander, data, landerStateEst)
% -------------------------------------------------------------------------
% orbLocalCoords - Function to compute local coordinates of lunar orbiter 
% with respect to moon lander
%
% Inputs:
%   xx      - Propagated state vectors of the lunar orbiter (ECI frame) [nx6]
%   lander - Structure containing lander data: 
%               - (...)
%   constants.mu  - gravitational parameter
%   data    - structure containing epoch data:
%               - data.et0      (initial epoch)
%               - data.etf      (final epoch)
%   landerStateEst  -  String to set the type of estimator for the lander state:
%                       - 'latrec' to estimate the state from latitudinal
%                          coordinates
%                       - 'kernels' to estimate the state from kernels
%
% Outputs:
%   localCoords - Structure containing evolution of range, azimuth, and 
%                 elevation of the lunar orbiter:
%                   - localCoords.range             [1xn]
%                   - localCoords.azimuth           [1xn]
%                   - localCoords.elevation         [1xn]
% -------------------------------------------------------------------------

% Extract Data and Constants:
tspan = data.tspan;
n = length(tspan);

% Compute state of lunar lander using latitudinal coordinates or kernels:
if strcmpi(landerStateEst, 'latrec')
    lander_pos = cspice_latrec(lander.r, lander.lon * cspice_rpd, lander.lat * cspice_rpd);
    lander_vel = [0; 0; 0];
    lander_stateMCMF = [lander_pos; lander_vel];
    % Transformation from MCMF frame to MCI:
    R_MCMF2MCI = zeros(6, 6, n);
    lander_stateMCI = zeros(6, n);
    for i = 1 : n
        R_MCMF2MCI(:, :, i) = cspice_sxform('IAU_MOON', 'J2000', tspan(i));
        lander_stateMCI(:, i) = R_MCMF2MCI(:, :, i) * lander_stateMCMF;
    end
elseif strcmpi(landerStateEst, 'kernels')
    lander_stateMCI = cspice_spkezr(lander.name, tspan, 'J2000', 'NONE', 'MOON');
end
    
% Compute lander-orbiter vector in MCI:
rv_lander_orb_mci = xx' - lander_stateMCI;

% Convert states from MCI into topocentric frame
% Initialize topocentric states array:
topoFrame = [lander.name, '_TOPO'];
orb_topo = zeros(6, n);
R_mci2topo = zeros(6, 6, n);
for i = 1 : n
    R_mci2topo(:, :, i) = cspice_sxform('J2000', topoFrame, tspan(i));
    orb_topo(:, i) = R_mci2topo(:, :, i) * rv_lander_orb_mci(:, i);
end

% Compute range, azimuth, and elevation using cspice_xfmsta
rll_station_sat = cspice_xfmsta(orb_topo, 'RECTANGULAR', 'LATITUDINAL', 'EARTH');
localCoords.range  = rll_station_sat(1, :);
localCoords.azimuth = rll_station_sat(2, :) * cspice_dpr;
localCoords.elevation = rll_station_sat(3, :) * cspice_dpr;

end

function [visibleCoords, errFlag] = checkVisibility(localCoords, lander)
% -------------------------------------------------------------------------
% Function to check that orbiter and lander are in relative visibility
%
% Inputs: localCoords.elevation: elevation evolution over time interval
%         lander.minEl         : minimum elevation constraint for lander
%
% Outputs:
%       visibleCoords          : struct containing visible range, azimuth,
%                                elevation during time interval
%       errFlag                : error flag to check relative visibility
% -------------------------------------------------------------------------

visibilityCond = localCoords.elevation >= lander.minEl;
errFlag = any(~visibilityCond);

if errFlag
    error('lander and orbiter are not in relative visibility during the entire time interval')
end

visibleCoords.range = localCoords.range(visibilityCond);
visibleCoords.azimuth = localCoords.azimuth(visibilityCond);
visibleCoords.elevation = localCoords.elevation(visibilityCond);

end

function [mean, cov] = UKF(mean, cov, ti, tf, measurements, data, settings, constants)

alpha = data.ut.alpha;
beta = data.ut.beta;
lander = data.lander;
% size of state vector to be estimated
stateSize = data.n;

if stateSize == 6

    % realPos = measurements;
    
    %%%%% PREDICTION:
    
    n = stateSize;
    
    % Compute sigma points at time t_{k-1}:
    chi = zeros(n, 2*n+1);
    ll = alpha^2 * n - n;
    mat = sqrtm((n + ll)*cov);
    
    chi(:, 1) = mean;
    for i = 1:n
        chi(:, i + 1)     = mean + mat(:, i);
        chi(:, i + n + 1) = mean - mat(:, i);
    end
    
    % Compute weights:
    weight_mean = zeros(2*n+1, 1);
    weight_cov  = zeros(2*n+1, 1);
    weight_mean(1) = ll/(n+ll);
    weight_cov(1)  = ll/(n+ll) + (1 - alpha^2 + beta);
    weight_mean(2:end) = 1/(2*(n+ll))*ones(2*n, 1);
    weight_cov(2:end)  = 1/(2*(n+ll))*ones(2*n, 1);
    
    % Propagate sigma points to time t_{k} and compute predicted mean and covariance:
    ss_points = zeros(n, 2*n+1);
    predMean = zeros(n, 1);
    predCov = zeros(n, n);
    for i = 1 : 2*n+1
        [~, x] = ode113(@TBP, [ti, tf], chi(:, i), settings.odeOpt, constants);
        ss_points(:, i) = x(end, :)';
        predMean = predMean + weight_mean(i)*ss_points(:, i);
    end
    for i = 1 : 2*n+1
        predCov = predCov + weight_cov(i)*(ss_points(:, i) - predMean)*(ss_points(:, i) - predMean)';
    end
    
    %%%%%% CORRECTION:
    
    % Add measurements from propagated sigma points:
    gamma = ss_points(1:3, :);
    y_mean = zeros(3, 1);
    % y_mean = gamma(:, 1) * weight_mean(1) + gamma(:, 2:end) * weight_mean(2:end)';
    for i = 1 : 2*n+1
        y_mean = y_mean + weight_mean(i) * gamma(:, i);
    end
    
    Pyy = weight_cov(1)*(gamma(:, 1) - y_mean)*(gamma(:, 1) - y_mean)';
    Pxy = weight_cov(1)*(ss_points(:, 1) - predMean)*(gamma(:, 1) - y_mean)';
    for jj = 2 : 2*n+1
        Pyy = Pyy + weight_cov(jj)*(gamma(:, jj) - y_mean)*(gamma(:, jj) - y_mean)';
        Pxy = Pxy + weight_cov(jj)*(ss_points(:, jj) - predMean)*(gamma(:, jj) - y_mean)';
    end
    % Add measurement noise:
    Pyy = Pyy + data.Rpos;
    
    % Kalman Gain:
    Kk = Pxy/Pyy;
    
    % Correction:
    corr = measurements' - y_mean;
    
    % Update mean and covariance:
    mean = predMean + Kk * corr;
    cov = predCov - Kk * Pyy * Kk';

elseif stateSize == 8
    n = stateSize;

    % Compute sigma points at time t_{k-1}:
    chi = zeros(n, 2*n+1);
    ll = alpha^2 * n - n;
    mat = sqrtm((n + ll)*cov);
    
    chi(:, 1) = mean;
    for i = 1:n
        chi(:, i + 1)     = mean + mat(:, i);
        chi(:, i + n + 1) = mean - mat(:, i);
    end
    
    % Compute weights:
    weight_mean = zeros(2*n+1, 1);
    weight_cov  = zeros(2*n+1, 1);
    weight_mean(1) = ll/(n+ll);
    weight_cov(1)  = ll/(n+ll) + (1 - alpha^2 + beta);
    weight_mean(2:end) = 1/(2*(n+ll))*ones(2*n, 1);
    weight_cov(2:end)  = 1/(2*(n+ll))*ones(2*n, 1);
    
    % Propagate sigma points to time t_{k} and compute predicted mean and covariance:
    ss_points = zeros(n, 2*n+1);
    predMean = zeros(n, 1);
    predCov = zeros(n, n);
    range = zeros(1, 2*n+1);
    for i = 1 : 2*n+1
        [~, x] = ode113(@TBP, [ti, tf], chi(1:6, i), settings.odeOpt, constants);
        ss_points(1:6, i) = x(end, :)';
        ss_points(7:8, i) = chi(7:8, i);
        predMean = predMean + weight_mean(i)*ss_points(:, i);
        % Compute range measurements from sigma points:
        lander.lat = ss_points(7, i) * cspice_dpr;
        lander.lon = ss_points(8, i) * cspice_dpr;
        range(i) = computeRange(ss_points(1:6, i)', tf, lander);
    end
    for i = 1 : 2*n+1
        predCov = predCov + weight_cov(i)*(ss_points(:, i) - predMean)*(ss_points(:, i) - predMean)';
    end
    
    %%%%%% CORRECTION:
    
    % Add measurements from propagated sigma points:

    position = ss_points(1:3, :);
    
    gamma = [position; range];
    y_mean = gamma * weight_mean;
        
    Pyy = weight_cov(1)*(gamma(:, 1) - y_mean)*(gamma(:, 1) - y_mean)';
    Pxy = weight_cov(1)*(ss_points(:, 1) - predMean)*(gamma(:, 1) - y_mean)';
    for jj = 2 : 2*n+1
        Pyy = Pyy + weight_cov(jj)*(gamma(:, jj) - y_mean)*(gamma(:, jj) - y_mean)';
        Pxy = Pxy + weight_cov(jj)*(ss_points(:, jj) - predMean)*(gamma(:, jj) - y_mean)';
    end
    % Add measurement noise:
    Pyy = Pyy + data.R_aug;
    
    % Kalman Gain:
    Kk = Pxy/Pyy;
    
    % Correction:
    corr = measurements' - y_mean;
    
    % Update mean and covariance:
    mean = predMean + Kk * corr;
    cov = predCov - Kk * Pyy * Kk';
end
end

function range = computeRange(xx, et, lander)

lander_pos = cspice_latrec(lander.r, lander.lon * cspice_rpd, lander.lat * cspice_rpd);
lander_vel = [0; 0; 0];
lander_stateMCMF = [lander_pos; lander_vel];
% Transformation from MCMF frame to MCI:
R_MCMF2MCI = cspice_sxform('IAU_MOON', 'J2000', et);
lander_stateMCI = R_MCMF2MCI*lander_stateMCMF;

% Compute lander-orbiter vector in MCI:
rv_lander_orb_mci = xx' - lander_stateMCI;

% Convert states from MCI into topocentric frame:
topoFrame = [lander.name, '_TOPO'];
R_mci2topo = cspice_sxform('J2000', topoFrame, et);
orb_topo = R_mci2topo * rv_lander_orb_mci;

% Compute range using cspice_xfmsta:
rll_station_sat = cspice_xfmsta(orb_topo, 'RECTANGULAR', 'LATITUDINAL', 'EARTH');
range  = rll_station_sat(1, :);

end

function plotSettings
%-------------------------------------------------------------------------%
% Settings for figure plots
%-------------------------------------------------------------------------%

% Setting Lines:
set(0, 'defaultLineLineWidth', 1.6);
set(0,'defaultLineMarkerSize', 4) ;
% set(0,'defaultLineMarkerEdgeColor', 'k')
set(0,'defaultLineMarkerFaceColor', 'auto')
% Setting Interpreter:
set(0, 'defaultTextInterpreter', 'latex')
set(0, 'defaultLegendInterpreter', 'latex')
set(0, 'defaultAxesTickLabelInterpreter', 'latex')
% Setting Legends:
set(0, 'defaultLegendLocation','southwest');
set(0, 'defaultLegendOrientation', 'vertical');
set(0, 'defaultLegendFontSize', 12);
% Setting Axes:
set(0, 'defaultAxesXMinorGrid', 'on');
set(0,'defaultAxesYMinorGrid','on');
set(0, 'defaultAxesFontSize', 20);

end

function plotVisEl(data, orbiter_lc)
% Convert epoch times (data.tspan) into time strings (hh:mm:ss)
time_labels = cell(length(data.tspan), 1); 
for i = 1:length(data.tspan)
    utc_full = cspice_et2utc(data.tspan(i), 'C', 0); 
    time_labels{i} = utc_full(12:end);
end

% Select 5 equidistant ticks on x axis
num_ticks = 5;
tick_indices = round(linspace(1, length(data.tspan), num_ticks)); 
tick_values = data.tspan(tick_indices); 
tick_labels = time_labels(tick_indices); 

figure()
plot(data.tspan, orbiter_lc.elevation, 'DisplayName', 'Orbiter Elevation')
hold on
plot(data.tspan, data.lander.minEl * ones(length(data.tspan)), '--', ...
     'HandleVisibility', 'off', 'Color', 'r')
plot(NaN, NaN, '--r', 'DisplayName', 'Minimum Visible Elevation') 
xlabel('NOV-18-2024', 'FontSize', 50)
ylabel('Elevation [deg]', 'FontSize', 50)
legend('FontSize', 35) 
grid on
xticks(tick_values)
xticklabels(tick_labels)
% xtickangle(45)
xlim([data.tspan(1), data.tspan(end)])
ylim([- 5, max(orbiter_lc.elevation) + 5])
set(gca, 'FontSize', 30)


end