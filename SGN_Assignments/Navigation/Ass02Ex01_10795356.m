% Spacecraft Guidance and Navigation (2024/2025)
% Assignment #2: Exercise #1
% Author: Luca Frassinella


% ------------------ EX. 1 - UNCERTAINTY PROPAGATION -------------------- %

clearvars; close all; clc;
cspice_kclear()
rng default
format long
plotSettings;
addpath('.\kernels\')
cspice_furnsh('assignment02.tm')

% Load parameters (mu, DU, TU, VU):
parameters = loadConstants();

% Set Initial Conditions:
n = 4;
r0 = [-0.011965533749906; -0.017025663128129];
v0 = [10.718855256727338; 0.116502348513671];
P0 = [1.041e-15, 6.026e-17, 5.647e-16, 4.577e-15;
      6.026e-17, 4.287e-18, 4.312e-17, 1.855e-16;
      5.647e-16, 4.312e-17, 4.432e-16, 1.455e-15;
      4.577e-15, 1.855e-16, 1.455e-15, 2.822e-14];
phi0 = reshape(eye(n), n*n, 1);

parameters.ti = 1.282800225339865;
parameters.tf = 9.595124551366348;
% Options for ode solver:
parameters.odeOptions = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
% Unscented Transform parameters:
parameters.ut.alpha = 1;
parameters.ut.beta = 2;
% Number of time-grid iterations:
N = 5;

%% Pt. 1: LinCov and UT
% Initialize means and covariances for LinCov:
lin.states = [[r0; v0; phi0], zeros(20, N)];
lin.x0 = [r0; v0];
lin.P0 = P0;
lin.P = zeros(n, n, N+1);
lin.P(:, :, 1) = P0;

% Initialize means and covariances for UT:
ut.states = [[r0; v0], zeros(n, N)];
ut.x0 = [r0; v0];
ut.P0 = P0;
ut.P = lin.P;

tvec = linspace(parameters.ti, parameters.tf, N+1);

% Update means and covariances using LinCov and UT:
for i = 2 : N+1
    parameters.tfinal = tvec(i);
    % LinCov
    [lin.states(:, i), lin.P(:, :, i)] = ...
        LinCov(lin.x0, lin.P0, parameters);

    % UT:
    [ut.states(:, i), ut.P(:, :, i), ut.sigma_points] = ...
        UT(ut.x0, ut.P0, parameters);
end

%% 2 - MonteCarlo
% Set the number of simulations
parameters.populationSize = 1000;

% Generate random initial conditions:
mc.x0 = mvnrnd([r0; v0], P0, parameters.populationSize)';

% Initialize mean and covariance matrices:
mc.x_mean = [[r0; v0], zeros(4, N)];
mc.P = zeros(4, 4, N + 1);
mc.P(:, :, 1) = P0;

% Monte Carlo simulations loop
for i = 2 : N + 1
     parameters.tfinal = tvec(i);
    [mc.x_mean(:, i), mc.P(:, :, i), finalStates] = monteCarlo(mc.x0, parameters);
end

%% Outputs:

% Plot of mean and covariance ellipses:
mean_ut = ut.states(:, end);
cov_ut = ut.P(:, :, end);
mean_lin = lin.states(:, end);
cov_lin = lin.P(:, :, end);
mean_mc = mc.x_mean(:, end);
cov_mc = mc.P(:, :, end);

figure()
hold on
covarianceEllipse(mean_ut(1:2), cov_ut(1:2, 1:2), 3)
hold on
covarianceEllipse(mean_lin(1:2), cov_lin(1:2, 1:2), 3)
hold on
covarianceEllipse(mean_mc(1:2), cov_mc(1:2, 1:2), 3)
hold on
plot(finalStates(1, :), finalStates(2, :), 'ko', 'MarkerSize', 0.5, 'MarkerEdgeColor', 'k')
legend('Covariance Ellipse $3\sigma$ (UT)', 'Mean Position (UT)', 'Covariance Ellipse $3\sigma$ (LinCov)', 'Mean Position (LinCov)', ...
    'Covariance Ellipse $3\sigma$ (MC)', 'Mean Position (MC)','Propagated MC Points' ,'FontSize', 23)
xlabel('x [-]', 'FontSize', 30)
ylabel('y [-]', 'FontSize', 30)
title('Comparison Between LinCov, UT and MC')

% Plot of time evolution of  max eigenvalues:
eigMaxPr = zeros(N, 3);
eigMaxPv = zeros(N, 3);
for i = 2 : N+1
    eigMaxPr(i,1) = 3*sqrt(max(eig(lin.P(1:2, 1:2, i))));
    eigMaxPr(i,2) = 3*sqrt(max(eig(ut.P(1:2, 1:2, i))));
    eigMaxPr(i,3) = 3*sqrt(max(eig(mc.P(1:2, 1:2, i))));

    eigMaxPv(i,1) = 3*sqrt(max(eig(lin.P(3:4, 3:4, i))));
    eigMaxPv(i,2) = 3*sqrt(max(eig(ut.P(3:4, 3:4, i))));
    eigMaxPv(i,3) = 3*sqrt(max(eig(mc.P(3:4, 3:4, i))));
end

figure()
plot(tvec, eigMaxPr(:, 1), '--x', 'Color', 'r', 'MarkerEdgeColor', 'r', 'MarkerSize', 8)
hold on
plot(tvec, eigMaxPr(:, 2), '--o', 'Color', 'b', 'MarkerEdgeColor', 'b', 'MarkerSize', 8)
plot(tvec, eigMaxPr(:, 3), '--^', 'Color', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 8)
grid on
ylim([-0.1e-04 20e-04])
hold off
xlabel('TU [-]', 'FontSize', 35)
ylabel('$3\sqrt{\mathrm{max}\lambda_i}P_{rr}$', 'FontSize', 35)
legend('LinCov', 'UT', 'MC', 'FontSize', 40)

figure()
plot(tvec, eigMaxPv(:, 1), '--x', 'Color', 'r', 'MarkerEdgeColor', 'r', 'MarkerSize', 8)
hold on
plot(tvec, eigMaxPv(:, 2), '--o', 'Color', 'b', 'MarkerEdgeColor', 'b', 'MarkerSize', 8)
plot(tvec, eigMaxPv(:, 3), '--^', 'Color', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 8)
grid on
ylim([-0.05 0.45])
hold off
xlabel('TU [-]', 'FontSize', 35)
ylabel('$3\sqrt{\mathrm{max}\lambda_i}P_{vv}$', 'FontSize', 35)
legend('LinCov', 'UT', 'MC', 'FontSize', 40)

% qqplots:
figure()
for i = 1 : 4
    subplot(2, 2, i)
    h = qqplot(finalStates(i, :));
    xlabel('Standard Normal Quantiles', 'Interpreter', 'latex')
    ylabel('Input Quantiles', 'Interpreter', 'latex')
    title('')
    switch i
        case 1
            legend([h(2), h(1)], {'Normal Quantile', 'MC $x$ component'}, 'FontSize', 25)
        case 2
            legend([h(2), h(1)], {'Normal Quantile', 'MC $y$ component'}, 'FontSize', 25)
        case 3
            legend([h(2), h(1)], {'Normal Quantile', 'MC $v_x$ component'}, 'FontSize', 25)
        case 4
            legend([h(2), h(1)], {'Normal Quantile', 'MC $v_y$ component'}, 'FontSize', 25)
    end
end

%% Functions

function [dxdt] = PBRFBP(t, xx, parameters)
x = xx(1);
y = xx(2);
vx = xx(3);
vy = xx(4);

% Constants:
mu = parameters.constants.mu;
ms = parameters.constants.ms;
rho = parameters.constants.rho;
om_s = parameters.constants.om_s;

dOMdx = x - (ms*cos(om_s*t))/rho^2 - (mu*(mu + x - 1))/((mu + x - 1)^2 + y^2)^(3/2) + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(3/2)) - (ms*(2*x - 2*rho*cos(om_s*t)))/(2*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2));
dOMdy = y - (ms*sin(om_s*t))/rho^2 - (mu*y)/((mu + x - 1)^2 + y^2)^(3/2) - (ms*(2*y - 2*rho*sin(om_s*t)))/(2*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2)) + (y*(mu - 1))/((mu + x)^2 + y^2)^(3/2);


if length(xx) == 20
    dxdt = zeros(20, 1);
    phi = reshape(xx(5:end), 4, 4);
    A = zeros(4, 4);
    A(1, 3) = 1;
    A(2, 4) = 1;
    A(3, 1) = (mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - ms/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2) - (3*(mu + x)^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*ms*(x - rho*cos(om_s*t))^2)/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2) + (3*mu*(mu + x - 1)^2)/((mu + x - 1)^2 + y^2)^(5/2) + 1;
    A(3, 2) = (3*ms*(2*x - 2*rho*cos(om_s*t))*(2*y - 2*rho*sin(om_s*t)))/(4*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2)) + (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2));
    A(3, 4) = 2;
    A(4, 1) = (3*ms*(2*x - 2*rho*cos(om_s*t))*(2*y - 2*rho*sin(om_s*t)))/(4*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2)) + (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2));
    A(4, 2) = (mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - ms/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*ms*(y - rho*sin(om_s*t))^2)/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2)^(5/2) + 1;
    A(4, 3) = -2;
    dphi = A*phi;
    dxdt(5:end) = reshape(dphi, 16, 1);
else
    dxdt = zeros(4, 1);
    dxdt(1:2) = xx(3:4);
    dxdt(3) = 2*vy + dOMdx;
    dxdt(4) = -2*vx + dOMdy;
end
end

function [x, P] = LinCov(x0, P0, parameters)
ti = parameters.ti;
tfinal = parameters.tfinal;

% % Solve the Planar Bi-circular Restricted Four Body Problem:
% % Integrate:
% optset = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
% [~, x] = ode78(@(t,xx) PBRFBP_STM(t, xx, parameters), [ti, tfinal], x0, optset);

% % Final state vector:
% x = x(end, :);
mu = parameters.constants.mu;
[xxf, ~, ~, ~, phi] = propagatorSTM(x0, ti, tfinal, mu);

% % Final Covariance Matrix:
% phi = reshape(x(5:end), 4, 4);
x = [xxf; phi(:)];
P = phi * P0 * phi';

end

function [y_mean, P, ss_points] = UT(x0, P0, parameters)
ti = parameters.ti;
tfinal = parameters.tfinal;
odeOptions = parameters.odeOptions;
alpha = parameters.ut.alpha;
beta = parameters.ut.beta;
n = length(x0);

% Initialize y_mean and P:
y_mean = zeros(n, 1);
P = zeros(n);

% Compute sigma points:
chi = zeros(n, 2*n+1);
ll = alpha^2 * n - n;
mat = sqrtm((n + ll)*P0);

chi(:, 1) = x0;
for i = 1:n
    chi(:, i + 1)     = x0 + mat(:, i);
    chi(:, i + n + 1) = x0 - mat(:, i);
end

% Compute weights:
weight_mean = zeros(2*n+1, 1);
weight_cov  = zeros(2*n+1, 1);
weight_mean(1) = ll/(n+ll);
weight_cov(1)  = ll/(n+ll) + (1 - alpha^2 + beta);
weight_mean(2:end) = 1/(2*(n+ll))*ones(2*n, 1);
weight_cov(2:end)  = 1/(2*(n+ll))*ones(2*n, 1);

% Propagate sigma points and compute y_mean and P:
ss_points = zeros(n, 2*n+1);
for i = 1 : 2*n+1
    % Solve the Planar Bi-circular Restricted Four Body Problem:
    [~, x] = ode78(@PBRFBP, [ti, tfinal], chi(:, i), odeOptions, parameters);
    ss_points(:, i) = x(end, :)';
    y_mean = y_mean + weight_mean(i)*ss_points(:, i);
end
for i = 1 : 2*n+1
    P = P + weight_cov(i)*(ss_points(:, i) - y_mean)*(ss_points(:, i) - y_mean)';
end

end

function [x_mean, P, finalStates] = monteCarlo(x0, parameters)

odeOptions = parameters.odeOptions;
ti = parameters.ti;
tfinal = parameters.tfinal;
populationSize = parameters.populationSize;

% Initialize matrix to store all final states:
finalStates = zeros(4, populationSize);

% Run simulations:
parfor i = 1 : populationSize
    [~, x] = ode78(@PBRFBP, [ti, tfinal], x0(:, i), odeOptions, parameters);
    finalStates(:, i) = x(end, :)';
end

% Mean and Covariance of the simulations:
x_mean = mean(finalStates, 2);
P = cov(finalStates');

end

function covarianceEllipse(mean, cov, sigma)

% Build refernce circle:
nPoints = 100;
alpha = 2*pi/nPoints * (0:nPoints);
circle = [cos(alpha); sin(alpha)];

% SVD to extract eigenvalues and first eigenvector of cov matrix:
[R, D] = svd(cov);
eig1 = D(1,1);  
eig2 = D(2,2);  
eigv1 = R(:,1); 

% Orientation angle (angle of the first eigenvector)
theta = atan2(eigv1(2), eigv1(1));

% Scale the circle along the principal axes
circle = [sqrt(eig1) * circle(1, :); sqrt(eig2) * circle(2, :)];

% Rotate the scaled circle using the angle theta:
x = cos(theta) * circle(1, :) - sin(theta) * circle(2, :);
y = sin(theta) * circle(1, :) + cos(theta) * circle(2, :);

% Center ellipse in the mean
x = mean(1) + sigma * x;
y = mean(2) + sigma * y;

% Plot:
hLine = plot(x, y, '--');
hold on;
plot(mean(1), mean(2), 'x', 'MarkerEdgeColor', hLine.Color, 'MarkerSize', 15);
grid on;


end

function parameters = loadConstants()

% Load Constants:
GM_Earth = cspice_bodvrd('EARTH', 'GM', 1);
GM_Moon = cspice_bodvrd('MOON', 'GM', 1);


parameters.constants.mu = GM_Moon/(GM_Earth + GM_Moon); % Earth-Moon system gravity constant
parameters.constants.ms = 3.28900541e05;       % [??]
parameters.constants.rho = 3.88811143e02;      % [??]
parameters.constants.om_s = -9.25195985e-01;   % [??]
parameters.constants.TU = 4.34811305;      % [days]    Time unit
parameters.constants.DU = 3.84405000e05;   % [km]      Distance unit
parameters.constants.VU = 1.02454018;      % [km/s]    Velocity unit
end

function [dxdt] = PBRFBP_STM(t, mu, xx)

% Inputs: xx = [20x1]
% Outputs: dxdt = [20x1]

x = xx(1);
y = xx(2);
vx = xx(3);
vy = xx(4);

PHI = reshape(xx(5: end), 4, 4);

% DEFINE CONSTANTS:
ms = 3.28900541e05;
rho = 3.88811143e02;
om_s = -9.25195985e-01;

dOMdx = x - (ms*cos(om_s*t))/rho^2 - (mu*(mu + x - 1))/((mu + x - 1)^2 + y^2)^(3/2) + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(3/2)) - (ms*(2*x - 2*rho*cos(om_s*t)))/(2*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2));
dOMdy = y - (ms*sin(om_s*t))/rho^2 - (mu*y)/((mu + x - 1)^2 + y^2)^(3/2) - (ms*(2*y - 2*rho*sin(om_s*t)))/(2*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2)) + (y*(mu - 1))/((mu + x)^2 + y^2)^(3/2);

A = zeros(4);
A(1, 3) = 1;
A(2, 4) = 1;
A(3, 1) = (mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - ms/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2) - (3*(mu + x)^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*ms*(x - rho*cos(om_s*t))^2)/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2) + (3*mu*(mu + x - 1)^2)/((mu + x - 1)^2 + y^2)^(5/2) + 1;
A(3, 2) = (3*ms*(2*x - 2*rho*cos(om_s*t))*(2*y - 2*rho*sin(om_s*t)))/(4*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2)) + (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2));
A(3, 4) = 2;
A(4, 1) = (3*ms*(2*x - 2*rho*cos(om_s*t))*(2*y - 2*rho*sin(om_s*t)))/(4*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2)) + (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2));
A(4, 2) = (mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - ms/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*ms*(y - rho*sin(om_s*t))^2)/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2)^(5/2) + 1;
A(4, 3) = -2;

PhiDot = A * PHI;

dxdt = zeros(20, 1);
dxdt(1:2) = xx(3:4);
dxdt(3) = 2*vy + dOMdx;
dxdt(4) = -2*vx + dOMdy;
dxdt(5:end) = PhiDot(:);
end

function [xxf, tf, xx, tt, PHI] = propagatorSTM(xx0, t0, tf, mu)

% Input: xx0 = 4x1
% Output: xx = [:, 4x1] (solo gli stati)
%         xxf = 4x1     (stato finale)
%         tf = 1x1      (tempo finale)
%         PHI = 4x4     (STM(t0, tf))
 
PHI0 = eye(4); % STM at t0
xx0 = [xx0; PHI0(:)]; % new vector of initial conditions (appended with STM)


optset = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[tt, xxv] = ode78(@(t,xx) PBRFBP_STM(t, mu, xx), [t0 tf], xx0, optset);

xx = xxv(:, 1:4)';
PHI = xxv(end, 5:end);
PHI = reshape(PHI, 4, 4);
xxf = xx(:, end);
tf = tt(end);

end
 
function plotSettings
%-------------------------------------------------------------------------%
% Settings for figure plots
%-------------------------------------------------------------------------%
% Setting Lines:
set(0, 'defaultLineLineWidth', 1.6);
set(0,'defaultLineMarkerSize', 4) ;
set(0,'defaultLineMarkerFaceColor', 'auto')
set(0, 'defaultLineMarkerEdgeColor', 'none')
% Setting Interpreter:
set(0, 'defaultTextInterpreter', 'latex')
set(0, 'defaultLegendInterpreter', 'latex')
set(0, 'defaultAxesTickLabelInterpreter', 'latex')
% Setting Legends:
set(0, 'defaultLegendLocation','northwest');
set(0, 'defaultLegendOrientation', 'vertical');
set(0, 'defaultLegendFontSize', 12);
% Setting Axes:
set(groot,'defaultAxesXMinorGrid','on','defaultAxesXMinorGridMode','manual');
set(groot,'defaultAxesYMinorGrid','on','defaultAxesYMinorGridMode','manual');
set(0, 'defaultAxesFontSize', 20);
set(groot, 'defaultAxesFontName', 'Palatino Linotype', 'defaultTextFontName', 'Palatino Linotype');
set(groot, 'defaultAxesFontWeight', 'bold');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
    
end