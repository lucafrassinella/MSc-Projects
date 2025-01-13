% Spacecraft Guidance and Navigation (2024/2025)
% Assignment #2: Exercise #3
% Author: Luca Frassinella


% --------------------- EX. 3 - SEQUENTIAL FILTERS ---------------------- %

clearvars; close all; clc;
cspice_kclear()
rng default
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

plotVisEl(data, orbiter_lc)
%% 2. Simulate Measurements:

% Simulate the measurements using kernels:
[orbiter_lc_kernels] = orbLocalCoords(xx, lander, data, 'kernels');

% Add the measurements noise:
idealMeasurements.position = xx(:, 1:3);
idealMeasurements.range = orbiter_lc_kernels.range;
realMeasurements.position = mvnrnd(idealMeasurements.position, data.Rpos);
realMeasurements.range = mvnrnd(idealMeasurements.range', data.sigmaMeas^2);

% Plot of the error with 3sigma boundaries:
figure()
plot(data.tspan, realMeasurements.range - idealMeasurements.range');
hold on
yline(3*data.sigmaMeas, 'LineWidth', 1)
yline(-3*data.sigmaMeas, 'LineWidth', 1)

figure()
plot(data.tspan, realMeasurements.position(:, 1) - xx(:, 1))
hold on
plot(data.tspan, realMeasurements.position(:, 2) - xx(:, 2))
plot(data.tspan, realMeasurements.position(:, 3) - xx(:, 3))
yline(3*data.sigmaMeas, 'LineWidth', 1)
yline(-3*data.sigmaMeas, 'LineWidth', 1)

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


stateLander = cspice_spkezr(lander.name, data.tspan, 'IAU_MOON', 'NONE', 'MOON');
[alt, real_lon, real_lat] = cspice_reclat(stateLander(1:3, :));
%% 
for i = 1:length(data.tspan)
    errPos2(i) = sqrt((xx(i,1) - augMean_mat(1,i))^2 + (xx(i,2) - augMean_mat(2,i))^2 +(xx(i,3) - augMean_mat(3,i))^2 );
    errVel2(i) = sqrt((xx(i,4) - augMean_mat(4,i))^2 + (xx(i,5) - augMean_mat(5,i))^2 +(xx(i,6) - augMean_mat(6,i))^2 );
    errLat(i) = sqrt((real_lat(i) - augMean_mat(7, i))^2) * cspice_dpr;
    errLon(i) = sqrt((real_lon(i) - augMean_mat(8, i))^2) * cspice_dpr;
    stdPos2(i) = 3 * sqrt(trace(augCov_mat(1:3,1:3, i)));
    stdVel2(i) = 3 * sqrt(trace(augCov_mat(4:6,4:6, i)));
    stdLat(i) = 3 * sqrt(augCov_mat(7,7,i)) * cspice_dpr;
    stdLon(i) = 3 * sqrt(augCov_mat(8,8,i)) * cspice_dpr;
end

figure
plot(data.tspan, errLat)
hold on
plot(data.tspan, stdLat)

figure
plot(data.tspan, errLon)
hold on
plot(data.tspan, stdLon)
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
xlabel('NOV-18-2024')
ylabel('Elevation [deg]')
legend('Location', 'northeast') 
grid on
xticks(tick_values)
xticklabels(tick_labels)
% xtickangle(45)
xlim([data.tspan(1), data.tspan(end)])
ylim([- 5, max(orbiter_lc.elevation) + 5])

end