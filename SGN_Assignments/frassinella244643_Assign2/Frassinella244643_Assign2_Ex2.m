% Spacecraft Guidance and Navigation (2024/2025)
% Assignment #2: Exercise #2
% Author: Luca Frassinella


% ----------------------- EX. 2 - BATCH FILTERS ------------------------- %

clearvars; close all; clc;
cspice_kclear()
plotSettings;
addpath('.\sgp4\')
addpath('.\tle\')
cspice_furnsh('assignment02.tm')

% Initialize data:
[kourou, troll, svalbard] = loadStationData();
[constants, data, settings] = loadSet();


%% 1. Compute Visibility Windows:

% Extract data from TLE file:
[r_eci, v_eci, data.refEpoch, data.sat] = TLE2car(data, settings);

% Save SMOS cartesian state (in ECI):
results.ex1.refPos = r_eci;
results.ex1.refVel = v_eci;

% Extract reference, initial and final epochs:
refEpoch = data.refEpoch;
et0 = data.et0;
etf = data.etf;

% Propagate from reference TLE epoch to initial epoch:
[~, xx_eci] = ode113(@TBP, [refEpoch et0], [r_eci; v_eci], settings.odeOpt, constants);
x0 = xx_eci(end,:)';

% Perform propagation from local ground stations:
[lc.kourou.TBP, lc.kourou.et_vec, ~] = scLocalCoords(x0, kourou, constants, data, settings, 'TBP', 0);
[lc.troll.TBP, lc.troll.et_vec, ~] = scLocalCoords(x0, troll, constants, data, settings, 'TBP', 0);
[lc.svalbard.TBP, lc.svalbard.et_vec, ~] = scLocalCoords(x0, svalbard, constants, data, settings, 'TBP', 0);

% Compute Visibility Windows:
[visibilityTimes.kourou_TBP, visibleCoords.kourou_TBP, timeWindow.kourou_TBP, visibilityCondition.kourou_TBP] = visibilityWindow(lc.kourou.TBP, kourou, lc.kourou.et_vec);
[visibilityTimes.troll_TBP, visibleCoords.troll_TBP, timeWindow.troll_TBP, visibilityCondition.troll_TBP] = visibilityWindow(lc.troll.TBP, troll, lc.troll.et_vec);
[visibilityTimes.svalbard_TBP, visibleCoords.svalbard_TBP, timeWindow.svalbard_TBP, visibilityCondition.svalbard_TBP] = visibilityWindow(lc.svalbard.TBP, svalbard, lc.svalbard.et_vec);

%%%% Plots:
plotLocalCoords(lc.kourou.TBP, kourou, lc.kourou.et_vec, visibilityTimes.kourou_TBP, visibleCoords.kourou_TBP.azimuth, visibleCoords.kourou_TBP.elevation)
plotLocalCoords(lc.troll.TBP, troll, lc.troll.et_vec, visibilityTimes.troll_TBP, visibleCoords.troll_TBP.azimuth, visibleCoords.troll_TBP.elevation)
plotLocalCoords(lc.svalbard.TBP, svalbard, lc.svalbard.et_vec, visibilityTimes.svalbard_TBP, visibleCoords.svalbard_TBP.azimuth, visibleCoords.svalbard_TBP.elevation)

results.ex1.visibility_windowKOUROU = timeWindow.kourou_TBP;
results.ex1.visibility_windowTROLL = timeWindow.troll_TBP;
results.ex1.visibility_windowSVALBARD = timeWindow.svalbard_TBP;

%% 2. Simulate Measurements:

% 2.a) Compute expected measurements during visibility window with SGP4:
[lc.kourou.SGP, ~, ~] = scLocalCoords(x0, kourou, constants, data, settings, 'SGP4', 0);
[lc.troll.SGP, ~, ~] = scLocalCoords(x0, troll, constants, data, settings, 'SGP4', 0);
[lc.svalbard.SGP, ~, ~] = scLocalCoords(x0, svalbard, constants, data, settings, 'SGP4', 0);

% Compute coordinates over visibility windows identified in Point 1:
visibleCoords.kourou_SGP.azimuth = lc.kourou.SGP.azimuth .* visibilityCondition.kourou_TBP;
visibleCoords.kourou_SGP.elevation = lc.kourou.SGP.elevation .* visibilityCondition.kourou_TBP;
visibleCoords.troll_SGP.azimuth = lc.troll.SGP.azimuth .* visibilityCondition.troll_TBP;
visibleCoords.troll_SGP.elevation = lc.troll.SGP.elevation .* visibilityCondition.troll_TBP;
visibleCoords.svalbard_SGP.azimuth = lc.svalbard.SGP.azimuth .* visibilityCondition.svalbard_TBP;
visibleCoords.svalbard_SGP.elevation = lc.svalbard.SGP.elevation .* visibilityCondition.svalbard_TBP;

% 2.b) Simulate real measurements adding random errors:
% Ideal Measurements:
idealMeasurements.kourou_mat = [lc.kourou.SGP.range', lc.kourou.SGP.azimuth', lc.kourou.SGP.elevation'];
idealMeasurements.kourou.range = idealMeasurements.kourou_mat(:, 1)';
idealMeasurements.kourou.azimuth = idealMeasurements.kourou_mat(:, 2)';
idealMeasurements.kourou.elevation = idealMeasurements.kourou_mat(:, 3)';
idealMeasurements.troll_mat = [lc.troll.SGP.range', lc.troll.SGP.azimuth', lc.troll.SGP.elevation'];
idealMeasurements.troll.range = idealMeasurements.troll_mat(:, 1)';
idealMeasurements.troll.azimuth = idealMeasurements.troll_mat(:, 2)';
idealMeasurements.troll.elevation = idealMeasurements.troll_mat(:, 3)';
idealMeasurements.svalbard_mat = [lc.svalbard.SGP.range', lc.svalbard.SGP.azimuth', lc.svalbard.SGP.elevation'];
idealMeasurements.svalbard.range = idealMeasurements.svalbard_mat(:, 1)';
idealMeasurements.svalbard.azimuth = idealMeasurements.svalbard_mat(:, 2)';
idealMeasurements.svalbard.elevation = idealMeasurements.svalbard_mat(:, 3)';

% Real Measurements using measurement covariance matrix for each station:
realMeasurements.kourou_mat = mvnrnd(idealMeasurements.kourou_mat, kourou.R);
realMeasurements.kourou.range = realMeasurements.kourou_mat(:, 1)';
realMeasurements.kourou.azimuth = realMeasurements.kourou_mat(:, 2)';
realMeasurements.kourou.elevation = realMeasurements.kourou_mat(:, 3)';
realMeasurements.troll_mat = mvnrnd(idealMeasurements.troll_mat, troll.R);
realMeasurements.troll.range = realMeasurements.troll_mat(:, 1)';
realMeasurements.troll.azimuth = realMeasurements.troll_mat(:, 2)';
realMeasurements.troll.elevation = realMeasurements.troll_mat(:, 3)';
realMeasurements.svalbard_mat = mvnrnd(idealMeasurements.svalbard_mat, svalbard.R);
realMeasurements.svalbard.range = realMeasurements.svalbard_mat(:, 1)';
realMeasurements.svalbard.azimuth = realMeasurements.svalbard_mat(:, 2)';
realMeasurements.svalbard.elevation = realMeasurements.svalbard_mat(:, 3)';

% Real visibility windows:
[visibilityTimes.kourou_real, visibleCoords.kourou_real, timeWindow.kourou_real, visibilityCondition.kourou_real] = visibilityWindow(realMeasurements.kourou, kourou, lc.kourou.et_vec);
[visibilityTimes.troll_real, visibleCoords.troll_real, timeWindow.troll_real, visibilityCondition.troll_real] = visibilityWindow(realMeasurements.troll, troll, lc.troll.et_vec);
[visibilityTimes.svalbard_real, visibleCoords.svalbard_real, timeWindow.svalbard_real, visibilityCondition.svalbard_real] = visibilityWindow(realMeasurements.svalbard, svalbard, lc.svalbard.et_vec);

%%%% Plots:
figure()
hold on
plot(visibleCoords.kourou_TBP.azimuth, visibleCoords.kourou_TBP.elevation, 'x', ...
    'MarkerSize', 15, 'MarkerEdgeColor', [0, 0.4470, 0.7410], 'LineWidth', 1.5)
plot(visibleCoords.kourou_real.azimuth, visibleCoords.kourou_real.elevation, '*', ...
    'MarkerSize', 15, 'MarkerEdgeColor', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.5) 
plot(visibleCoords.troll_TBP.azimuth, visibleCoords.troll_TBP.elevation, 'x', ...
    'MarkerSize', 15, 'MarkerEdgeColor', [0.9290, 0.6940, 0.1250], 'LineWidth', 1.5) 
plot(visibleCoords.troll_real.azimuth, visibleCoords.troll_real.elevation, '*', ...
    'MarkerSize', 15, 'MarkerEdgeColor', [0.4940, 0.1840, 0.5560], 'LineWidth', 1.5) 
plot(visibleCoords.svalbard_TBP.azimuth, visibleCoords.svalbard_TBP.elevation, 'x', ...
    'MarkerSize', 15, 'MarkerEdgeColor', [0.4660, 0.6740, 0.1880], 'LineWidth', 1.5) 
plot(visibleCoords.svalbard_real.azimuth, visibleCoords.svalbard_real.elevation, '*', ...
    'MarkerSize', 15, 'MarkerEdgeColor', [0.3010, 0.7450, 0.9330], 'LineWidth', 1.5)
legend('KOUROU (Kepler)', 'KOUROU (SGP4 with noise)', 'TROLL (Kepler)', ...
    'TROLL (SGP4 with noise)', 'SVALBARD (Kepler)', 'SVALBARD (SGP4 with noise)')
axis([-180, 180, 0, 90])
ax = gca;
ax.FontSize = 25;
xlabel('Azimuth [deg]', 'FontSize', 30)
ylabel('Elevation [deg]', 'FontSize', 30)


%% 3. Solve the Navigation Problem:

%%%%%%%% 3.a) Least Squares Solution using KOUROU measurements only:
stations_a{1, 1} = kourou;
stations_a{1, 2} = [visibilityTimes.kourou_real', visibleCoords.kourou_real.range', ...
                  visibleCoords.kourou_real.azimuth', visibleCoords.kourou_real.elevation'];
odeFun = @(t, x) TBP(t, x, constants);
% Propagate to initial guess:
[~, xx_guess] = ode113(@TBP, [refEpoch et0], [r_eci; v_eci], settings.odeOpt, constants);
x0_guess = xx_guess(end,:)';
% Cost function
objective = @(x) costFunction(x, stations_a, odeFun, data, settings);
% Solve the navigation problem using lsqnonlin:
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter', 'StepTolerance', 1e-10);
[navigationSol_a.state, resnorm1, residual1, exitflag1, ~, ~, jacobian1] = lsqnonlin(objective, x0_guess, [], [], options);

% Results:
navigationSol_a.J = full(jacobian1);
navigationSol_a.P = resnorm1/(length(residual1)-length(x0_guess)).*inv(navigationSol_a.J.'*navigationSol_a.J);
navigationSol_a.pos = sqrt(trace(navigationSol_a.P(1:3,1:3)));
navigationSol_a.vel = sqrt(trace(navigationSol_a.P(4:6,4:6)));
% Linear Mapping:
[navigationSol_a.std_a, navigationSol_a.std_i] = linMap(navigationSol_a.state, navigationSol_a.P, constants);

%%%%%%% 3.b) Least Squares Solution using all stations measurements:
stations_b = cell(3, 2);
stations_b{1, 1} = kourou;
stations_b{1, 2} = [visibilityTimes.kourou_real', visibleCoords.kourou_real.range', ...
                  visibleCoords.kourou_real.azimuth', visibleCoords.kourou_real.elevation'];
stations_b{2, 1} = troll;
stations_b{2, 2} = [visibilityTimes.troll_real', visibleCoords.troll_real.range', ...
                  visibleCoords.troll_real.azimuth', visibleCoords.troll_real.elevation'];
stations_b{3, 1} = svalbard;
stations_b{3, 2} = [visibilityTimes.svalbard_real', visibleCoords.svalbard_real.range', ...
                  visibleCoords.svalbard_real.azimuth', visibleCoords.svalbard_real.elevation'];

odeFun = @(t, x) TBP(t, x, constants);
% Propagate to initial guess:
[~, xx_guess] = ode113(@TBP, [refEpoch et0], [r_eci; v_eci], settings.odeOpt, constants);
x0_guess = xx_guess(end,:)';
% Cost function:
objective = @(x) costFunction(x, stations_b, odeFun, data, settings);
% Solve the navigation problem using lsqnonlin:
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter', 'StepTolerance', 1e-10);
[navigationSol_b.state, resnorm2, residual2, exitflag2, ~, ~, jacobian2] = lsqnonlin(objective, x0_guess, [], [], options);

% Results:
navigationSol_b.J = full(jacobian2);
navigationSol_b.P = resnorm2/(length(residual2)-length(x0_guess)).*inv(navigationSol_b.J.'*navigationSol_b.J);
navigationSol_b.pos = sqrt(trace(navigationSol_b.P(1:3,1:3)));
navigationSol_b.vel = sqrt(trace(navigationSol_b.P(4:6,4:6)));
% Linear Mapping:
[navigationSol_b.std_a, navigationSol_b.std_i] = linMap(navigationSol_b.state, navigationSol_b.P, constants);

%%%%%%% 3.c) Least Squares Solution with J2 perturbation:
stations_c = stations_b;
odeFun = @(t, x) PTBP(t, x, constants);
% Propagate to initial guess:
[~, xx_guess] = ode113(@PTBP, [refEpoch et0], [r_eci; v_eci], settings.odeOpt, constants);
x0_guess = xx_guess(end,:)';
% Cost function:
objective = @(x) costFunction(x, stations_c, odeFun, data, settings);
% Solve the navigation problem using lsqnonlin:
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter', 'StepTolerance', 1e-10);
[navigationSol_c.state, resnorm3, residual3, exitflag3, ~, ~, jacobian3] = lsqnonlin(objective, x0_guess, [], [], options);

% Results:
navigationSol_c.J = full(jacobian3);
navigationSol_c.P = resnorm3/(length(residual3)-length(x0_guess)).*inv(navigationSol_c.J.'*navigationSol_c.J);
navigationSol_c.pos = sqrt(trace(navigationSol_c.P(1:3,1:3)));
navigationSol_c.vel = sqrt(trace(navigationSol_c.P(4:6,4:6)));
% Linear Mapping:
[navigationSol_c.std_a, navigationSol_c.std_i] = linMap(navigationSol_c.state, navigationSol_c.P, constants);

results.ex3.a = navigationSol_a;
results.ex3.b = navigationSol_b;
results.ex3.c = navigationSol_c;
%% 4. Trade-Off Analysis:

% set J2-perturbed model for the dynamics:
odeFun = @(t, x) PTBP(t, x, constants);
% set lsqnonlin options:
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter', 'StepTolerance', 1e-10);

% CASE 1: KOUROU - TROLL ground stations
stations_case1 = cell(2, 2);
stations_case1{1, 1} = kourou;
stations_case1{1, 2} = [visibilityTimes.kourou_real', visibleCoords.kourou_real.range', ...
                  visibleCoords.kourou_real.azimuth', visibleCoords.kourou_real.elevation'];
stations_case1{2, 1} = troll;
stations_case1{2, 2} = [visibilityTimes.troll_real', visibleCoords.troll_real.range', ...
                  visibleCoords.troll_real.azimuth', visibleCoords.troll_real.elevation'];
% Cost function:
objective = @(x) costFunction(x, stations_case1, odeFun, data, settings);
% Solve the navigation problem using lsqnonlin:
[navigationSol_case1.state, resnorm_c1, residual_c1, exitflag_c1, ~, ~, jacobian_c1] = lsqnonlin(objective, x0_guess, [], [], options);
navigationSol_case1.J = full(jacobian_c1);
navigationSol_case1.P = resnorm_c1/(length(residual_c1)-length(x0_guess)).*inv(navigationSol_case1.J.'*navigationSol_case1.J);
navigationSol_case1.pos = sqrt(trace(navigationSol_case1.P(1:3,1:3)));
navigationSol_case1.vel = sqrt(trace(navigationSol_case1.P(4:6,4:6)));
% Linear Mapping:
[navigationSol_case1.std_a, navigationSol_case1.std_i] = linMap(navigationSol_case1.state, navigationSol_case1.P, constants);

% CASE 2: KOUROU - SVALBARD ground stations
stations_case2 = cell(2, 2);
stations_case2{1, 1} = kourou;
stations_case2{1, 2} = [visibilityTimes.kourou_real', visibleCoords.kourou_real.range', ...
                  visibleCoords.kourou_real.azimuth', visibleCoords.kourou_real.elevation'];
stations_case2{2, 1} = svalbard;
stations_case2{2, 2} = [visibilityTimes.svalbard_real', visibleCoords.svalbard_real.range', ...
                  visibleCoords.svalbard_real.azimuth', visibleCoords.svalbard_real.elevation'];
% Cost function:
objective = @(x) costFunction(x, stations_case2, odeFun, data, settings);
% Solve the navigation problem using lsqnonlin:
[navigationSol_case2.state, resnorm_c2, residual_c2, exitflag_c2, ~, ~, jacobian_c2] = lsqnonlin(objective, x0_guess, [], [], options);
navigationSol_case2.J = full(jacobian_c2);
navigationSol_case2.P = resnorm_c2/(length(residual_c2)-length(x0_guess)).*inv(navigationSol_case2.J.'*navigationSol_case2.J);
navigationSol_case2.pos = sqrt(trace(navigationSol_case2.P(1:3,1:3)));
navigationSol_case2.vel = sqrt(trace(navigationSol_case2.P(4:6,4:6)));
% Linear Mapping:
[navigationSol_case2.std_a, navigationSol_case2.std_i] = linMap(navigationSol_case2.state, navigationSol_case2.P, constants);

% CASE 3: TROLL - SVALBARD ground stations
stations_case3 = cell(2, 2);
stations_case3{1, 1} = troll;
stations_case3{1, 2} = [visibilityTimes.troll_real', visibleCoords.troll_real.range', ...
                  visibleCoords.troll_real.azimuth', visibleCoords.troll_real.elevation'];
stations_case3{2, 1} = svalbard;
stations_case3{2, 2} = [visibilityTimes.svalbard_real', visibleCoords.svalbard_real.range', ...
                  visibleCoords.svalbard_real.azimuth', visibleCoords.svalbard_real.elevation'];
% Cost function:
objective = @(x) costFunction(x, stations_case3, odeFun, data, settings);
% Solve the navigation problem using lsqnonlin:
[navigationSol_case3.state, resnorm_c3, residual_c3, exitflag_c3, ~, ~, jacobian_c3] = lsqnonlin(objective, x0_guess, [], [], options);
navigationSol_case3.J = full(jacobian_c3);
navigationSol_case3.P = resnorm_c3/(length(residual_c3)-length(x0_guess)).*inv(navigationSol_case3.J.'*navigationSol_case3.J);
navigationSol_case3.pos = sqrt(trace(navigationSol_case3.P(1:3,1:3)));
navigationSol_case3.vel = sqrt(trace(navigationSol_case3.P(4:6,4:6)));
% Linear Mapping:
[navigationSol_case3.std_a, navigationSol_case3.std_i] = linMap(navigationSol_case3.state, navigationSol_case3.P, constants);

results.ex4.Kourou_Troll = navigationSol_case1;
results.ex4.Kourou_Svalbard = navigationSol_case2;
results.ex4.Troll_Svalbard = navigationSol_case3;

%% 5. Long-Term Analysis:

% Set time of observation campaign to reference epoch +/- 5 hours:
data_LTA.et0 = data.refEpoch - 5*60*60;
data_LTA.etf = data.refEpoch + 5*60*60;
data_LTA.refEpoch = data.refEpoch;
data_LTA.sat = data.sat;

[kourou_LTA, tv_kourou, ~] = scLocalCoords(x0, kourou, constants, data_LTA, settings, 'SGP4', 0);
[troll_LTA, tv_troll, ~] = scLocalCoords(x0, troll, constants, data_LTA, settings, 'SGP4', 0);
[svalbard_LTA, tv_svalbard, ~] = scLocalCoords(x0, svalbard, constants, data_LTA, settings, 'SGP4', 0);

% Ideal Measurements:
kourou_LTA_id = [kourou_LTA.range', kourou_LTA.azimuth', kourou_LTA.elevation'];
troll_LTA_id = [troll_LTA.range', troll_LTA.azimuth', troll_LTA.elevation'];
svalbard_LTA_id = [svalbard_LTA.range', svalbard_LTA.azimuth', svalbard_LTA.elevation'];

% Real Measurements using measurement covariance matrix for each station:
kourou_LTA_real = mvnrnd(kourou_LTA_id, kourou.R);
troll_LTA_real = mvnrnd(troll_LTA_id, troll.R);
svalbard_LTA_real = mvnrnd(svalbard_LTA_id, svalbard.R);

% Extract real simulated elevation:
kourou_LTA_el = kourou_LTA_real(:, 3)';
troll_LTA_el = troll_LTA_real(:, 3)';
svalbard_LTA_el = svalbard_LTA_real(:, 3)';

% Plot over augmented time window:
% Kourou:
figure()
plot(tv_kourou / cspice_spd(), kourou_LTA_el, '-k')
hold on
plot(data.refEpoch*ones(100, 1)/cspice_spd, linspace(min(kourou_LTA_el)-10, max(kourou_LTA_el)+10, 100), '--r', 'DisplayName', 'Reference Time')
plot(tv_kourou / cspice_spd(), ones(length(tv_kourou))*kourou.minEl, '--b', 'DisplayName', 'Minimum Elevation')
xlim([tv_kourou(1) tv_kourou(end)]/cspice_spd())
ylim([min(kourou_LTA_el)-10 max(kourou_LTA_el)+10])
num_ticks = 5;
ax = gca;
ax.FontSize = 25;
tick_indices = round(linspace(1, length(tv_kourou), num_ticks)); 
tick_values = tv_kourou(tick_indices) / cspice_spd();
tick_labels = cell(num_ticks, 1); 
for i = 1:num_ticks
    utc_full = cspice_et2utc(tv_kourou(tick_indices(i)), 'C', 0); 
    tick_labels{i} = utc_full(12:end); 
end
xticks(tick_values);
xticklabels(tick_labels);
xtickangle(0);
xlabel('18-NOV-2024', 'FontSize', 40)
ylabel('Elevation [deg]', 'FontSize', 40)
legend('Elevation', 'Reference Time', 'Minimum Elevation', 'FontSize', 35)
title('SMOS elevation over KOUROU (Long-Term Analysis')

% Troll:
figure()
plot(tv_troll / cspice_spd(), troll_LTA_el, '-k')
hold on
plot(data.refEpoch*ones(100, 1)/cspice_spd, linspace(min(troll_LTA_el)-10, max(troll_LTA_el)+10, 100), '--r', 'DisplayName', 'Reference Time')
plot(tv_troll / cspice_spd(), ones(length(tv_troll))*troll.minEl, '--b', 'DisplayName', 'Minimum Elevation')
xlim([tv_troll(1) tv_troll(end)]/cspice_spd())
ylim([min(troll_LTA_el)-10 max(troll_LTA_el)+10])
num_ticks = 5;
ax = gca;
ax.FontSize = 25;
tick_indices = round(linspace(1, length(tv_troll), num_ticks)); 
tick_values = tv_troll(tick_indices) / cspice_spd();
tick_labels = cell(num_ticks, 1); 
for i = 1:num_ticks
    utc_full = cspice_et2utc(tv_troll(tick_indices(i)), 'C', 0); 
    tick_labels{i} = utc_full(12:end); 
end
xticks(tick_values);
xticklabels(tick_labels);
xtickangle(0);
xlabel('18-NOV-2024', 'FontSize', 40)
ylabel('Elevation [deg]', 'FontSize', 40)
legend('Elevation', 'Reference Time', 'Minimum Elevation', 'FontSize', 35)
title('SMOS elevation over TROLL (Long-Term Analysis')

% Svalbard:
figure()
plot(tv_svalbard / cspice_spd(), svalbard_LTA_el, '-k')
hold on
plot(data.refEpoch*ones(100, 1)/cspice_spd, linspace(min(svalbard_LTA_el)-10, max(svalbard_LTA_el)+10, 100), '--r', 'DisplayName', 'Reference Time')
plot(tv_svalbard / cspice_spd(), ones(length(tv_svalbard))*svalbard.minEl, '--b', 'DisplayName', 'Minimum Elevation')
xlim([tv_svalbard(1) tv_svalbard(end)]/cspice_spd())
ylim([min(svalbard_LTA_el)-10 max(svalbard_LTA_el)+10])
num_ticks = 5;
ax = gca;
ax.FontSize = 25;
tick_indices = round(linspace(1, length(tv_svalbard), num_ticks)); 
tick_values = tv_svalbard(tick_indices) / cspice_spd();
tick_labels = cell(num_ticks, 1); 
for i = 1:num_ticks
    utc_full = cspice_et2utc(tv_svalbard(tick_indices(i)), 'C', 0); 
    tick_labels{i} = utc_full(12:end); 
end
xticks(tick_values);
xticklabels(tick_labels);
xtickangle(0);
xlabel('18-NOV-2024', 'FontSize', 40)
ylabel('Elevation [deg]', 'FontSize', 40)
legend('Elevation', 'Reference Time', 'Minimum Elevation', 'FontSize', 35)
title('SMOS elevation over SVALBARD (Long-Term Analysis')

%% Clear Workspace:
clearvars -except data data_LTA idealMeasurements realMeasurements kourou troll svalbard settings results

%% Functions

function [station1, station2, station3] = loadStationData()
% -------------------------------------------------------------------------
% loadStationData - function to load data for the Kourou, Troll and
%                   Svalbard Ground Stations
%
%   Output:
%       - station1: struct containing data of Kourou Ground Station
%       - station2: struct containing data of Troll Ground Station
%       - station3: struct containing data of Svalbard Ground Station
% -------------------------------------------------------------------------


station1.name = 'KOUROU';
station1.lat = 5.25144;       % [deg]
station1.lon = -52.80466;     % [deg]
station1.alt = -14.67;        % [m]
station1.minEl = 6;           % [deg]
station1.measFreq = 60;       % [s]
station1.R = diag([0.01^2; 0.125^2; 0.125^2]); % [km; deg; deg]^2
station1.Wm = inv(sqrtm(station1.R));

station2.name = 'TROLL';
station2.lat = -72.011977;     % [deg]
station2.lon = 2.536103;       % [deg]
station2.alt = 1298;           % [m]
station2.minEl = 0;            % [deg]
station2.measFreq = 30;        % [s]
station2.R = diag([0.01^2; 0.125^2; 0.125^2]); % [km; deg; deg]^2
station2.Wm = inv(sqrtm(station2.R));

station3.name = 'SVALBARD';
station3.lat = 78.229772;   % [deg]
station3.lon = 15.407786;   % [deg]
station3.alt = 458;         % [m]
station3.minEl = 8;         % [deg]
station3.measFreq = 60;     % [s]
station3.R = diag([0.01^2; 0.125^2; 0.125^2]); % [km; deg; deg]^2
station3.Wm = inv(sqrtm(station3.R));

end

function [constants, data, settings] = loadSet()
% -------------------------------------------------------------------------
% loadSet - function to load constants, data and settings
%
%   Output:
%       - constants: struct containing constants
%       - data: struct containing data
%       - settings: struct containing settings
% -------------------------------------------------------------------------

constants.mu = cspice_bodvrd('Earth', 'GM', 1);

% ODE Solver Options
settings.odeOpt = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
settings.arcsec2rad = pi/(180*3600);
% SGP4 Settings:
settings.typerun    = 'u';  
settings.opsmode    = 'a';  
settings.whichconst =  72; 

data.ID = 36036; % Spacecraft ID
% Initialize data in UTC:
data.t0 = '2024-11-18T20:30:00.000';
data.tf = '2024-11-18T22:15:00.000';
% Initial and Final Epochs:
data.et0 = cspice_str2et(data.t0);
data.etf = cspice_str2et(data.tf);

end

function [r_eci, v_eci, refEpoch, sat] = TLE2car(data, settings)
% -------------------------------------------------------------------------
% TLE2car - Function to convert Two-Line Element (TLE) to Cartesian state
%
% Inputs:
%   data         - structure containing data 
%   settings     - structure containing settings
%
% Outputs:
%   r_eci        - extracted position vector (ECI) [3x1]
%   v_eci        - extracted velocity vector (ECI) [3x1]
%   refEpoch     - reference epoch of the TLE data in ET
%   sat          - structure containing extracted TLE data from .3le file
% -------------------------------------------------------------------------

% Extract settings and data:
arcsec2rad = settings.arcsec2rad;
whichconst = settings.whichconst; 
ID = data.ID;

% Compute reference data from TLE file:
sat = read_3LE(ID, 'tle\36036.3le', whichconst);
[year, month, day, hour, min, sec] = invjday(sat.jdsatepoch, sat.jdsatepochf);
sat_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f',[year, month, day, hour, min, sec]);
refEpoch = cspice_str2et(sat_epoch_str);

% Compute position and velocity in TEME reference frame:
[~, rrTEME, vvTEME] = sgp4(sat, 0.0);

% Transformation from TEME to ECI reference frame:
ddpsi = -0.114752*arcsec2rad; %  [rad]
ddeps = -0.007529*arcsec2rad; %  [rad]
ttt = cspice_unitim(refEpoch, 'ET', 'TDT')/cspice_jyear()/100;

[r_eci, v_eci, ~] = teme2eci(rrTEME, vvTEME, [0;0;0], ttt, ddpsi, ddeps);

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

% Equations of motion
dxx(1:3) = xx(4:6);
dxx(4:6) = -constants.mu/r_norm^3 * [x;y;z];

% Transpose to have a column vector:
dxx = dxx';
    
end

function [dxx] = PTBP(t, xx, constants)
% ----------------------------------------------------------------------- %
% Computes the time derivative of the state vector for Keplerian taking
% into account the J2 perturbation motion. 
%
% INPUT:
%   et          : [1, 1] Ephemeris time (TDB) in seconds.
%   xx          : [6, 1] State vector [x, y, z, vx, vy, vz] in ECI where:
%                 x, y, z   - Position components in Cartesian coordinates (km)
%                 vx, vy, vz - Velocity components in Cartesian coordinates (km/s)
%   constants.mu: [1, 1] Standard gravitational parameter of the central body (km^3/s^2)

%
% OUTPUT:
%   dxx          : [6, 1] Time derivative of the state vector, representing:
%                 dxdt = [vx, vy, vz, ax, ay, az], where:
%                 ax, ay, az - Acceleration components due to gravity (km/s^2)
% ----------------------------------------------------------------------- %

mu = constants.mu;
R_EARTH = cspice_bodvrd('EARTH','RADII',3);
RE = R_EARTH(1);
% State :
rECInorm = norm(xx(1:3)) ;
rECI = xx(1:3) ;
M_ROT = cspice_pxform('J2000', 'ITRF93', t);
rECEF = M_ROT * rECI ;
rECEFnorm   = norm(rECEF) ;


J2 = 0.0010826269;
Coeff_J2 = (3/2 * J2 * mu * RE^2/rECEFnorm^5);

aJ2ECEF =   [Coeff_J2*rECEF(1)*(5*rECEF(3)^2/rECEFnorm^2 - 1); 
             Coeff_J2*rECEF(2)*(5*rECEF(3)^2/rECEFnorm^2 - 1);
             Coeff_J2*rECEF(3)*(5*rECEF(3)^2/rECEFnorm^2 - 3)];

aJ2ECI = M_ROT' * aJ2ECEF ;

dxx = zeros(6,1);
dxx(1:3) = xx(4:6);
dxx(4) = -mu/rECInorm^3*xx(1) +  aJ2ECI(1);
dxx(5) = -mu/rECInorm^3*xx(2) +  aJ2ECI(2);
dxx(6) = -mu/rECInorm^3*xx(3) +  aJ2ECI(3);

end

function [station, GS_eci, ROT_eci2topo] = coordsGS(station, et_vec)
% ----------------------------------------------------------------------- %
% coordsGS - Function to compute Ground Station coordinates in ECI frame 
% and rotation matrix to topocentric frame over a give time span
%
% Inputs:
%   station.name    - Ground Station name (e.g., 'KOUROU').
%   et_vec          - Vector of times over which coordinates are computed
%
% Outputs:
%   station.stateECI    - Ground Station coordinates in ECI frame
%   station.eci2topo    - Rotation matrix from ECI to topocentric frame
% ----------------------------------------------------------------------- %

GS_name = station.name;

% Define station topocentric name
GS_topoFrame = [GS_name, '_TOPO'];

% Compute station state in ECI frame
GS_eci = cspice_spkezr(GS_name, et_vec, 'J2000', 'NONE', 'EARTH');

% Transformation from ECI to topocentric frame:
ROT_eci2topo = cspice_sxform('J2000', GS_topoFrame, et_vec);

% Outputs:
station.stateECI = GS_eci;
station.eci2topo = ROT_eci2topo;

end

function [localCoords, et_vec, xx] = scLocalCoords(x0, station, constants, data, settings, propagatorType, J2)
% -------------------------------------------------------------------------
% scLocalCoords - Function to compute local coordinates of a spacecraft 
% with respect to a Ground Station
%
% Inputs:
%   x0      - Initial satellite state vector (ECI frame) [6x1]
%   station - Structure containing station fields: 
%               - stateECI (ECI state of the station)
%               - eci2topo (ECI to topocentric frame rotation matrix)
%               - R (Measurement noise matrix)
%               - minEl (Minimum elevation for visibility)
%               - visibility_time (Time of visibility)
%   constants.mu  - gravitational parameter
%   data    - structure containing epoch data:
%               - data.refEpoch (reference epoch of TLE data)
%               - data.et0      (initial epoch)
%               - data.etf      (final epoch)
%   propagatorType  -  String to set the type of propagator between 'TBP' or 'SGP4'
%   J2      - variable to take into account J2 perturbation, 
%             if = 1 takes into account J2, if = 0 it does not
%
%
% Outputs:
%   localCoords - Structure containing evolution of range, azimuth, and 
%                 elevation of the satellite:
%                   - localCoords.range             [1xN]
%                   - localCoords.azimuth           [1xN]
%                   - localCoords.elevation         [1xN]
%   et_vec      - Vector of discretized epochs      [1xN]
%   xx          - propagated state vectors          [6xN]
% -------------------------------------------------------------------------

% Extract Data and Constants:
refEpoch = data.refEpoch;
et0 = data.et0;
etf = data.etf;

% Define ephemerides time vector starting from GS measurement frequency:
npoints = round((etf-et0)/station.measFreq)+1;
et_vec = linspace(et0, etf, npoints);

% Initialize state vectors:
arcsec2rad = settings.arcsec2rad;
ddpsi = -0.114752*arcsec2rad; %  [rad]
ddeps = -0.007529*arcsec2rad; %  [rad]
rECI = zeros(3, npoints);
vECI = zeros(3, npoints);
xx = zeros(npoints, 6);

% Compute satellite and ground station state vectors in ECI frame:
if strcmpi(propagatorType, 'TBP')
    % Include J2 perturbation or not
    if J2 == 0
       [~, xx] = ode113(@TBP, et_vec, x0, settings.odeOpt, constants);
    elseif J2 == 1
       [~, xx] = ode113(@PTBP, et_vec, x0, settings.odeOpt, constants);
    end

elseif strcmpi(propagatorType, 'SGP4')
  for i = 1:npoints
    sat = data.sat;
    % minutes from TLE epoch:
    tsince = (et_vec(i) - refEpoch)/60.0;
    [~, TEME.r , TEME.v] = sgp4(sat,  tsince);
    
    % Compute centuries from 2000-01-01T00:00:00.00 TDT
    ttt = cspice_unitim(et_vec(i), 'ET', 'TDT')/cspice_jyear()/100;
    
    % TEME to ECI conversion
    [rECI(:,i), vECI(:,i)] = teme2eci(TEME.r, TEME.v, [0.0;0.0;0.0], ttt, ddpsi, ddeps);
    xx(i,:) = [rECI(:,i); vECI(:,i)]';
  end
end

% Station state vectors in ECI frame:
[station] = coordsGS(station, et_vec);

% Compute station-satellite vector in ECI:
rv_station_sat_eci = xx' - station.stateECI;

% Initialize topocentric states array:
n = size(rv_station_sat_eci, 2);
sat_topo = zeros(6, n);

% Convert states from ECI into topocentric frame
for i = 1:n
    sat_topo(:, i) = station.eci2topo(:, :, i) * rv_station_sat_eci(:, i);
end

% Compute range, azimuth, and elevation using cspice_xfmsta
rll_station_sat = cspice_xfmsta(sat_topo, 'RECTANGULAR', 'LATITUDINAL', 'EARTH');
sat_Az = rll_station_sat(2, :) * cspice_dpr;
sat_El = rll_station_sat(3, :) * cspice_dpr;
sat_range = rll_station_sat(1, :);

% Outputs:
localCoords.range = sat_range;
localCoords.azimuth = sat_Az;
localCoords.elevation = sat_El;

end

function [visibilityTimes, visibleCoords, timeWindow, visibilityCondition] = visibilityWindow(localCoords, station, tspan)
% ----------------------------------------------------------------------- %
% visibilityWindow - Function to extract visibility windows and corresponding 
%              coordinates based on ground station's minimum elevation
%              constraint
%
% Inputs:
%   localCoords   - struct containing range, azimuth, elevation from local 
%                   ground station
%   station.minEl - Minimum elevation for visibility from ground station
%   et_vec        - time span of observation campaign
%
% Outputs:
%   visibilityTimes - Extracted visibility times from ground station
%   visibleCoords   - Visible coordinates during time window
%   timeWindow      - Array of initial and final time of visibility window
%   visibilityCondition - logical array containing visible times
% ----------------------------------------------------------------------- %

% Extract visibility condition based on minimum elevation
elevation = localCoords.elevation;
stationMinEl = station.minEl;

visibilityCondition = elevation >= stationMinEl;

% Extract visibility time and coordinates
visibilityTimes = tspan(visibilityCondition);
visibleCoords.range = localCoords.range(:, visibilityCondition);
visibleCoords.azimuth = localCoords.azimuth(:, visibilityCondition);
visibleCoords.elevation = localCoords.elevation(:, visibilityCondition);
timeWindow = [cspice_et2utc(visibilityTimes(1), 'C', 3); cspice_et2utc(visibilityTimes(end), 'C', 3)];

end

function residuals = costFunction(x0, stations, odeFun, data, settings)
% -------------------------------------------------------------------------
% costFUnction - Function to compute the residuals for lsqnonlin
%                optimization
% Inputs:
%   x0         - Initial state vector
%   stations   - Cell array containing information about stations
%   odeFun     - Function handle for the ODE solver (TPB/PTBP)
%   data       - structure containing problem data
%   settings   - structure containing problem settings
%
% Output:
%   residual   - Residuals of the cost function to be minimized
% -------------------------------------------------------------------------

% Collect data from ground stations:
% Initialize index:
ii_old = 0;
% Extract number of stations:
nStations = size(stations, 1);

for i = 1 : nStations
    % Extract number of measurements for the i-th station:
    nMeasures = size(stations{i, 2}, 1);
    ii_new = ii_old + nMeasures;
    % Vector of visibility epochs:
    t_vec(ii_old+1 : ii_new) = stations{i, 2}(:, 1)';
    % Real measurements during visible epochs:
    realMeasurements(ii_old+1 : ii_new, :) = stations{i, 2}(:, 2:4);
    % Measurement weights at each discretized visible epoch:
    measWeights(:, :, ii_old+1 : ii_new) = stations{i, 1}.Wm .* ones(3, 3, ii_new - ii_old);
    % Ground station states (ECI) and transformation to topocentric frame during visible epochs:
    [~, GS_eci(:, ii_old+1:ii_new), ROT_eci2topo(:, :, ii_old+1:ii_new)] = coordsGS(stations{i, 1}, stations{i, 2}(:, 1)');

    % Update index:
    ii_old = ii_new;
end

% Sort measurements from different stations in increasing time: 
t_vec = t_vec(1:ii_new);
[t_vec, index] = sort(t_vec);
t_vec = t_vec';
realMeasurements = realMeasurements(index, :);
measWeights = measWeights(:, :, index);
station.stateECI = GS_eci(:, index);
station.eci2topo = ROT_eci2topo(:, :, index);

% Propagate states at discretized visible epochs:
t_propagation = [data.et0; t_vec];
[~, x] = ode113(odeFun, t_propagation, x0, settings.odeOpt);
propStates = x(2:end, :);

% Convert Propagated States to Range, Azimuth and Elevation:
rv_station_sat_eci = propStates' - station.stateECI;
n = size(rv_station_sat_eci, 2);
sat_topo = zeros(6, n);
% Convert states from ECI into topocentric frame
for i = 1:n
    sat_topo(:, i) = station.eci2topo(:, :, i) * rv_station_sat_eci(:, i);
end
% Compute range, azimuth, and elevation:
rll_station_sat = cspice_xfmsta(sat_topo, 'RECTANGULAR', 'LATITUDINAL', 'EARTH');
range = rll_station_sat(1, :);
Az = rll_station_sat(2, :) * cspice_dpr;
El = rll_station_sat(3, :) * cspice_dpr;

predictedMeas = [range', Az', El'];

% Compute Residuals:
predictedMeas(:, 2) = deg2rad(predictedMeas(:, 2));
predictedMeas(:, 3) = deg2rad(predictedMeas(:, 3));
realMeasurements(:, 2) = deg2rad(realMeasurements(:, 2));
realMeasurements(:, 3) = deg2rad(realMeasurements(:, 3));
residuals = zeros(length(t_vec), 3);
for i = 1 : length(t_vec)
    weightDiff = measWeights(:, :, i) * [predictedMeas(i, 1) - realMeasurements(i, 1), rad2deg(angdiff(predictedMeas(i, 2), realMeasurements(i, 2))), rad2deg(angdiff(predictedMeas(i, 3), realMeasurements(i, 3)))]';
    residuals(i, :) = weightDiff';
end
end

function [std_a, std_i] = linMap(state, P, constants)
% -------------------------------------------------------------------- %
% Compute the linear mapping in terms of standard deviation of semimajor
% axis and inclination.
%
% INPUT:
%   state      : [6x1] Cartesian state vector [x, y, z, vx, vy, vz]
%   Pcar       : [6x6] Covariance matrix of the cartesian state vector
%   constants.mu    : [1x1] Gravitational parameter of the central body [km^3/s^2].
%
% OUTPUT:
%   Pkep       : [6x6] Covariance matrix mapping Cartesian to Keplerian elements.
%   std_a      : [1x1] standard deviation of semimajor axis [km]
%   std_i      : [1x1] standard deviation of inclination [deg]
% -------------------------------------------------------------------- %

mu = constants.mu;
% Initialize Jacobian
J = zeros(6, 6);
r = state(1:3);
v = state(4:6);
 
% Perturb position
for k = 1:3
    epsR = zeros(3,1);
    epsR(k) = sqrt(eps) * max(1, abs(r(k)));

    % Forward perturbation
    [afin, efin, ifin, OMfin, omfin, thfin] = car2kep(r + epsR, v, mu);
    [ain, ein, iin, OMin, omin, thin] = car2kep(r, v, mu);
% Computes the Jacobian from Cartesian to Keplerian elements using finite differences.
    % Numerical derivative
    J(:, k) = ([afin, efin, ifin, OMfin, omfin, thfin]' - [ain, ein, iin, OMin, omin, thin]') ./ epsR(k);
end

% Perturb velocity
for k = 1:3
    epsV = zeros(3,1);
    epsV(k) = sqrt(eps) * max(1, abs(v(k)));

    % Forward and backward perturbation
    [afin, efin, ifin, OMf, omfin, thfin] = car2kep(r, v + epsV, mu);
    [ain, ein, iin, OMin, omin, thin] = car2kep(r, v, mu);
% Computes the Jacobian from Cartesian to Keplerian elements using finite differences.
    % Numerical derivative
    J(:, k+3) = ([afin, efin, ifin, OMf, omfin, thfin]'- [ain, ein, iin, OMin, omin, thin]') ./ epsV(k);
end

% Compute the covariance
Pkep = J * P * J';

std_a = sqrt(Pkep(1,1)); % Standard deviation of semimajor axis [km]
std_i = rad2deg(sqrt(Pkep(3,3))); % Standard deviation of inclination [deg]
end

function [a, e, i, OM, om, th, N] = car2kep(r, v, mu)
% -------------------------------------------------------------------- %
% car2kep -  Converts Cartesian coordinates to Keplerian parameters (Output in radians)
%
% INPUT:
%   r [3x1] : Position vector [km]
%   v [3x1] : Velocity vector [km/s]
%   mu [1x1] : Gravitational parameter [km^3/s^2]
%
% OUTPUT:
%   a  [1x1] : Semi-major axis [km]
%   e  [1x1] : Eccentricity [-]
%   i  [1x1] : Inclination [rad]
%   OM [1x1] : Right Ascension of the Ascending Node (RAAN) [rad]
%   om [1x1] : Argument of periapsis [rad]
%   th [1x1] : True anomaly [rad]
% -------------------------------------------------------------------- %

I = [1 0 0]';
J = [0 1 0]';
K = [0 0 1]';

r_norm = norm(r);
v_norm = norm(v);

h = cross(r,v);
h_norm = norm(h);


i = acos(dot(h,K)/h_norm);


e = 1/mu * ((v_norm^2 - mu/r_norm)*r - (dot(r,v))*v);
e_norm = norm(e);


Eps = 1/2 * v_norm^2 - mu/r_norm;
a = -mu/(2*Eps);

N = cross(K,h);
N_norm = norm(N);


if dot(N,J) >= 0
 OM = acos(dot(N,I)/N_norm);
else
    OM = 2*pi - acos(dot(N,I)/N_norm);
end
if i == 0
    OM = 0;
    N = [1;0;0];
    N_norm = norm(N);
end

if dot(e,K) >= 0
    om = acos(dot(N,e)/(N_norm * e_norm));
else
    om = 2*pi -acos(dot(N,e)/(N_norm * e_norm));
end

if dot(v,r)/r_norm >=0
    th = acos(dot(r,e)/(r_norm * e_norm));
else
    th = 2*pi - acos(dot(r,e)/(r_norm * e_norm));
end

e = e_norm;
end

function plotLocalCoords(lc_station, station, tspan, visibilityTimesVec, visibleAzimuth, visibleElevation)
% -------------------------------------------------------------------- %
% Function to plot azimuth/elevation and visibility windows
% -------------------------------------------------------------------- %

azimuth = lc_station.azimuth;
elevation = lc_station.elevation;
NAME = station.name;

figure()
subplot(1, 2, 1);
hold on
yyaxis left % Left axis for azimuth
plot(tspan / cspice_spd(), azimuth, 'DisplayName', 'Azimuth')
hold on
plot(visibilityTimesVec / cspice_spd(), visibleAzimuth, 'Color',[0, 0, 139] / 255, 'LineStyle','-', 'LineWidth', 3,  'DisplayName', 'Visible Azimuth')
ylabel('Azimuth [deg]','Color','k', 'FontSize', 30) 
ax = gca;
ax.YColor = 'k';
yyaxis right % Right axis for elevation
plot(tspan / cspice_spd(), elevation, 'DisplayName', 'Elevation')
plot(visibilityTimesVec / cspice_spd(), visibleElevation, 'Color',[139, 0, 0] / 255, 'LineStyle','-', 'LineWidth', 3,  'DisplayName', 'Visible Elevation')
grid on
ylabel('Elevation [deg]','Color','k', 'FontSize', 30) 
ax.YColor = 'k';
title(NAME + " station")
% Convert epoch times (tspan) into time strings (hh:mm:ss)
time_labels = cell(length(tspan), 1); 
for i = 1:length(tspan)
    utc_full = cspice_et2utc(tspan(i), 'C', 0); 
    time_labels{i} = utc_full(12:end); % Extract only the time part (hh:mm:ss)
end
% Select 5 equidistant ticks on x axis
num_ticks = 5;
tick_indices = round(linspace(1, length(tspan), num_ticks)); 
tick_values = tspan(tick_indices) / cspice_spd(); % Convert to days
tick_labels = time_labels(tick_indices);
% Set x-axis ticks and labels
xticks(tick_values);
xticklabels(tick_labels);
xtickangle(45);
xlabel('18-NOV-2024', 'FontSize', 30)
legend('Location', 'best', 'FontSize', 20);
hold on
xlim([tspan(1) tspan(end)] / cspice_spd())

subplot(1, 2, 2);
plot(visibleAzimuth, visibleElevation, 'x' , 'MarkerSize', 15)
legend('SMOS Visible Coordinates', 'Location', 'north')
axis([-180, 180, 0, 90])
ax = gca;
ax.FontSize = 25;
xlabel('Azimuth [deg]', 'FontSize', 30)
ylabel('Elevation [deg]', 'FontSize', 30)

end

function plotSettings
%-------------------------------------------------------------------------%
% Settings for figure plots
%-------------------------------------------------------------------------%

% Setting Lines:
set(0, 'defaultLineLineWidth', 1.6);
set(0,'defaultLineMarkerSize', 4) ;
set(0,'defaultLineMarkerEdgeColor', 'k')
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