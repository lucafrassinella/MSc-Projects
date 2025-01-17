% Spacecraft Guidance and Navigation - A.Y. 2024/25
% Assignment #1
% Exercise #2
% Author: Frassinella Luca - 10795356

%% -------------------- EX. 2 - IMPULSIVE GUIDANCE --------------------- %%

clearvars; close all; clc;
cspice_kclear()
format long
plotSettings;
addpath('.\kernels\')
cspice_furnsh('ex02.tm');
tic
% Load constants and data:
[data, constants] = loadDataset();
settings.odeOpt = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% Extract data and constants:
mu = constants.mu;
ti = data.ti;
tf = data.tf;
r0 = data.r0; 
rf = data.rf;
v0 = data.v0;
vf = data.vf;
alpha = data.alpha;

%% Pt. 1: Initial Guess Solution

% Construct initial state vector from data:
xx0 = [r0*cos(alpha) - mu; r0*sin(alpha); -(v0 - r0)*sin(alpha); (v0 - r0)*cos(alpha)];
results.ex1.xx0 = xx0;

% Propagate and Plot initial guess solution:
[tt, xx] = propagate(xx0, ti, tf, constants, settings);

% Plot:
figure()
plot(xx(:, 1), xx(:, 2), 'k', 'DisplayName', 'Trajectory')
grid on
hold on
axis equal
[xPark_i, yPark_i] = circularOrbit(r0, [-mu; 0]); % initial parking orbit
[xPark_f, yPark_f] = circularOrbit(rf, [1-mu; 0]);  % final parking orbit
plot(xPark_i, yPark_i, '--k', 'LineWidth', 0.8, 'DisplayName', 'Initial parking orbit') 
plot(xPark_f, yPark_f, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6, 'DisplayName', 'Final parking orbit')
xlabel('x [-]', 'FontSize', 20)
ylabel('y [-]', 'FontSize', 20)
title('Trajectory from initial guess', 'FontWeight', 'bold', 'FontSize', 20)
subtitle('[@EMB Earth-Moon Rotating Frame]', 'FontSize', 16)
xlim([-1.5 3])
legend('FontSize', 18);

% Rotate to ECI and Plot initial guess solution:
XX = rot2ECI(tt, xx, constants);

% Plot:
figure()
plot(XX(:, 1), XX(:, 2), 'k', 'DisplayName', 'Trajectory')
hold on
[xParkECI, yParkECI] = circularOrbit(r0, [0; 0]);
[xMoon, yMoon] = circularOrbit(1, [0; 0]);
plot(xParkECI, yParkECI, '--k', 'LineWidth', 0.8, 'DisplayName', 'Initial parking orbit')
plot(xMoon, yMoon, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6, 'DisplayName', 'Moon orbit') 
grid on
axis equal
xlabel('x [-]')
ylabel('y [-]')
title('Trajectory from initial guess', 'FontWeight', 'bold', 'FontSize', 20)
subtitle('[@Earth ECI]', 'FontSize', 16)
xlim([-1.25 2.5])
legend('FontSize', 18);

%% Pt. 2: Simple Shooting Optimization:

%%%% a) Simple Shooting without analytical derivatives of the
%%%% objective/constraints:

% Setting options for fmincon:
options = optimoptions('fmincon', 'Display', 'iter-detailed', 'Algorithm', 'active-set','ConstraintTolerance', 1e-10, 'OptimalityTolerance', 1e-10, 'FiniteDifferenceStepSize',1e-12);
% Define Initial Guess:
vars0 = [xx0; ti; tf];

% Cost Function:
J = @(vars) costFunction(vars, data, constants);

% Nonlinear Constraints:
nonlCon = @(vars) constraints(vars, data, constants);

% Solve optimization problem and store results:
[vars_opt, dv_opt] = fmincon(J, vars0, [], [], [], [], [], [], nonlCon, options);

% Compute errors at arrival orbit:
[~, ceq] = nonlCon(vars_opt);
% Dimensionalize errors [m], [m/s]:
errDim = [ceq(3) * constants.DU*1e3; ceq(4) * constants.VU*1e3];

results.ex2a.vars_opt = vars_opt;
results.ex2a.dv_opt = dv_opt;
results.ex2a.errPos = errDim(1);
results.ex2a.errVel = errDim(2);

% Plots:
% Propagate trajectory and rotate it to ECI:
[tt, xx] = propagate(vars_opt(1:4), vars_opt(5), vars_opt(6), constants, settings);
XX = rot2ECI(tt, xx, constants);

% Plot:
figure()
plot(xx(:, 1), xx(:, 2), 'k', 'DisplayName', 'Trajectory')
grid on
hold on
axis equal
[xPark_i, yPark_i] = circularOrbit(r0, [-mu; 0]); % initial parking orbit
[xPark_f, yPark_f] = circularOrbit(rf, [1-mu; 0]);  % final parking orbit
plot(xPark_i, yPark_i, '--k', 'LineWidth', 0.8, 'DisplayName', 'Initial parking orbit') 
plot(xPark_f, yPark_f, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6, 'DisplayName', 'Final parking orbit')
title('Simple Shooting Optimized Trajectory (w/o gradients)', 'FontWeight', 'bold')
subtitle('[@EMB Earth-Moon Rotating Frame]', 'FontSize', 18)
xlim([-1.3 1.45])
ylim([-0.1 1.7])
xlabel('x [-]', 'FontSize', 23)
ylabel('y [-]', 'FontSize', 23)
set(gca, 'FontSize', 16)
legend('FontSize', 18);

figure()
plot(XX(:, 1), XX(:, 2), 'k', 'DisplayName', 'Trajectory')
hold on
[xParkECI, yParkECI] = circularOrbit(r0, [0; 0]);
[xMoon, yMoon] = circularOrbit(1, [0; 0]);
plot(xParkECI, yParkECI, '--k', 'LineWidth', 0.8, 'DisplayName', 'Initial parking orbit')
plot(xMoon, yMoon, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6, 'DisplayName', 'Moon orbit') 
grid on
axis equal
title('Simple Shooting Optimized Trajectory (w/o gradients)', 'FontWeight', 'bold')
subtitle('[@Earth ECI]', 'FontSize', 18)
xlim([-1.7 1.72])
ylim([-1.2 1.22])
xlabel('x [-]', 'FontSize', 23)
ylabel('y [-]', 'FontSize', 23)
set(gca, 'FontSize', 16)
legend('FontSize', 18);

%%
%%%% b) Simple Shooting with analytical derivatives of the
%%%% objective/constraints:

% Update options:
options = optimoptions('fmincon', 'Display', 'iter-detailed', 'Algorithm', 'active-set','ConstraintTolerance', 1e-10, 'OptimalityTolerance', 1e-10, 'FiniteDifferenceStepSize',1e-12, ...
    'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true);

% Cost Function:
J = @(vars) costFunction2(vars, data, constants);

% Nonlinear Constraints:
nonlCon2 = @(vars) constraints2(vars, data, constants);

% Solve optimization problem and store results:
[vars_opt2, dv_opt2] = fmincon(J, vars0, [], [], [], [], [], [], nonlCon2, options);

% Compute errors at arrival orbit:
[~, ceq] = nonlCon2(vars_opt2);
% Dimensionalize errors [m], [m/s]:
errDim = [ceq(3) * constants.DU*1e3; ceq(4) * constants.VU*1e3];

results.ex2b.vars_opt = vars_opt2;
results.ex2b.dv_opt = dv_opt2;
results.ex2b.errPos = errDim(1);
results.ex2b.errVel = errDim(2);

% Plots:
% Propagate trajectory and rotate it to ECI:
[tt, xx] = propagate(vars_opt2(1:4), vars_opt2(5), vars_opt2(6), constants, settings);
XX_ss = rot2ECI(tt, xx, constants);

figure()
plot(xx(:, 1), xx(:, 2), 'k', 'DisplayName', 'Trajectory')
grid on
hold on
axis equal
[xPark_i, yPark_i] = circularOrbit(r0, [-mu; 0]); % initial parking orbit
[xPark_f, yPark_f] = circularOrbit(rf, [1-mu; 0]);  % final parking orbit
plot(xPark_i, yPark_i, '--k', 'LineWidth', 0.8, 'DisplayName', 'Initial parking orbit') 
plot(xPark_f, yPark_f, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6, 'DisplayName', 'Final parking orbit')
xlabel('x [-]')
ylabel('y [-]')
title('Simple Shooting Optimized Trajectory (with gradients)', 'FontWeight', 'bold')
subtitle('[@EMB Earth-Moon Rotating Frame]', 'FontSize', 18)
xlim([-1.3 1.45])
ylim([-0.1 1.7])
xlabel('x [-]', 'FontSize', 23)
ylabel('y [-]', 'FontSize', 23)
set(gca, 'FontSize', 16)
legend('FontSize', 18);

figure()
plot(XX_ss(:, 1), XX_ss(:, 2), 'k', 'DisplayName', 'Trajectory')
hold on
[xParkECI, yParkECI] = circularOrbit(r0, [0; 0]);
[xMoon, yMoon] = circularOrbit(1, [0; 0]);
plot(xParkECI, yParkECI, '--k', 'LineWidth', 0.8, 'DisplayName', 'Initial parking orbit')
plot(xMoon, yMoon, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6, 'DisplayName', 'Moon orbit') 
grid on
axis equal
title('Simple Shooting Optimized Trajectory (with gradients)', 'FontWeight', 'bold')
subtitle('[@Earth ECI]', 'FontSize', 18)
xlim([-1.7 1.72])
ylim([-1.2 1.22])
xlabel('x [-]', 'FontSize', 23)
ylabel('y [-]', 'FontSize', 23)
set(gca, 'FontSize', 16)
legend('FontSize', 18);

%% Pt.3: Multiple Shooting Optimization:

% Number of nodes:
N = 4; 
% Discretized time vector:
t = zeros(N, 1);
for j = 1:N
    t(j) = ti + (j-1)/(N-1) * (tf - ti);
end

% Define Initial guess for Multiple Shooting starting from xx0:
x_ms = zeros(4,N);
x_ms(:,1) = xx0;
settings.odeOpt = odeset('RelTol', 1e-11, 'AbsTol', 1e-11);

for j = 2 : N
 [~,xx] = propagate(x_ms(:,j-1), t(j-1), t(j), constants, settings);
 x_ms(:,j) = xx(end,:)';
end

vars0_ms = zeros(4*N+2,1);
vars0_ms(1:4*N) = reshape(x_ms, 4*N, 1);
vars0_ms(4*N+1) = ti;
vars0_ms(4*N+2) = tf;

% Set options
options = optimoptions('fmincon', 'Algorithm', 'active-set', 'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true, 'Display','iter-detailed', 'MaxFunctionEvaluations',1e6, 'ConstraintTolerance',1e-7, 'OptimalityTolerance',1e-10, 'MaxIterations',10000); 
% Optimization:
[vars_opt_ms, dv_opt_ms] = fmincon(@(vars) costFunction_ms(vars, data, constants), vars0_ms, [],[],[],[],[],[], @(vars) constraints_ms(vars, data, constants), options);

% Compute errors at arrival orbit:
[c, ceq] = constraints_ms(vars_opt_ms, data, constants);
% Dimensionalize errors [m], [m/s]:
errDim = [ceq(end-1) * constants.DU*1e3; ceq(end) * constants.VU*1e3];

results.ex3.vars_opt = vars_opt_ms;
results.ex3.dv_opt = dv_opt_ms;
results.ex3.errPos = errDim(1);
results.ex3.errVel = errDim(2);
% Store dimensionalized errors at internal nodes (1, 2, 3) for position and
% velocity constraints:
results.ex3.errPos_nodes = [norm([ceq(1), ceq(2)]*constants.DU*1e3); ...
                            norm([ceq(5), ceq(6)]*constants.DU*1e3); ...
                            norm([ceq(9), ceq(10)]*constants.DU*1e3);];
results.ex3.errVel_nodes = [norm([ceq(3), ceq(4)]*constants.VU*1e3); ...
                            norm([ceq(7), ceq(8)]*constants.VU*1e3); ...
                            norm([ceq(11), ceq(12)]*constants.VU*1e3);];

% Plots:
plotMultipleShooting(vars_opt_ms, N, constants, data)

%% Pt.4: N-Body Propagation

% Upper/Lower Bounds for et:
tFormat = 'YYYY-MON-DD-HR:MN:SC.####::TDB';
tLstr = 'September 28 00:00:00.000 TDB 2024';
etL = cspice_str2et(tLstr);

% Extract Constants:
ms = constants.ms;
rho = constants.rho;
om_s = constants.om_s;
om_m = constants.om_m;

% Moon Revolution Period:
moonRev = 2*pi/om_m;
etU = etL + moonRev;

% Initial state and time:
xxi = results.ex2b.vars_opt(1:4)';
ti = results.ex2b.vars_opt(5);
tf = results.ex2b.vars_opt(6);

% Time of flight (in seconds):
tof = (tf-ti)*constants.TU*86400;

% Rotate initial state to ECI, append zi = 0, vzi = 0, dimensionalize:
xxi_ECI  = rot2ECI(ti, xxi, constants)';
xxi_ECI = [xxi_ECI(1:2); 0; xxi_ECI(3:4); 0];
xxi_ECI = [xxi_ECI(1:3) .* constants.DU; xxi_ECI(4:6) .* constants.VU]; % [km], [km/s]

results.ex4.xx0_eci = xxi_ECI;

% Define target angle:
thetaTarget = mod(om_s * ti, 2*pi); % target angle [0, 2*pi]

% Initialize Epoch Vector:
etVec = linspace(etL, etU, 10000);

thetaFinderFun = @(et) thetaFinder(et, thetaTarget, 'J2000', 'EMB');
% Solve the zero-finding problem:
etGuess = 781169000;
initialEpoch = fzero(thetaFinderFun, [etGuess etU]); 
results.ex4.initialEpoch = initialEpoch;
results.ex4.initialDateStr = cspice_et2utc(initialEpoch, 'C', 3);

% N-body propagation:
finalEpoch = initialEpoch + tof;
results.ex4.finalDateStr = cspice_et2utc(finalEpoch, 'C', 3);
xxi = xxi_ECI;
labels = {'Sun'; 'Mercury'; 'Venus'; 'Earth'; 'Moon'; 'Mars Barycenter'; 'Jupiter Barycenter';
          'Saturn Barycenter'; 'Uranus Barycenter'; 'Neptune Barycenter'; 'Pluto Barycenter'};

% Save integration frame string for SPICE:
center = 'Earth';
frame = 'J2000';
[bodies] = nBody_init(labels);
[tt, xx] = nBodyPropagator(xxi, initialEpoch, finalEpoch, bodies, frame, center);

figure
plot(xx(:, 1), xx(:, 2), 'Color', [0 0.5 1])
hold on
plot(XX_ss(:, 1)*constants.DU, XX_ss(:, 2)*constants.DU, 'k')
[xParkECI, yParkECI] = circularOrbit(r0, [0; 0]);
[xMoon, yMoon] = circularOrbit(1, [0; 0]);
plot(xParkECI*constants.DU, yParkECI*constants.DU, '--k', 'LineWidth', 0.8, 'DisplayName', 'Initial parking orbit')
plot(xMoon*constants.DU, yMoon*constants.DU, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6, 'DisplayName', 'Moon orbit') 
legend('N-Body', 'PBRFBP', 'Initial Orbit', 'Moon Orbit')
title('Comparison Between N-Body and PBRFBP Propagations')
subtitle('(@Earth ECI)')
grid on
axis equal
xlim([-1.7*constants.DU 1.72*constants.DU])
ylim([-1.2*constants.DU 1.22*constants.DU])
xlabel('x [km]', 'FontSize', 23)
ylabel('y [km]', 'FontSize', 23)
set(gca, 'FontSize', 16)
legend('FontSize', 18);

%% Clear Workspace:
clearvars -except data constants results

%% Functions

function [data, constants] = loadDataset()
% -------------------------------------------------------------------------
% loadSet - function to load constants and data
%
%   Output:
%       - constants: struct containing constants
%       - data: struct containing data
% -------------------------------------------------------------------------

% Load Constants:
GM_Earth = cspice_bodvrd('EARTH', 'GM', 1);
GM_Moon = cspice_bodvrd('MOON', 'GM', 1);
R_Earth = cspice_bodvrd('EARTH', 'RADII', 3);
R_Moon = cspice_bodvrd('MOON', 'RADII', 3);

constants.mu = GM_Moon/(GM_Earth + GM_Moon); % Earth-Moon system gravity constant
constants.TU = 4.34256461;      % [days]    Time unit
constants.DU = 3.84405000e05;   % [km]      Distance unit
constants.VU = 1.02454018;      % [km/s]    Velocity unit

% Earth, Sun & Moon constants:
constants.R_Earth = R_Earth(1);     % [km]      Earth radius
constants.R_Moon = R_Moon(1);       % [km]      Moon radius
constants.ms = 3.28900541e05;       % [-]
constants.rho = 3.88811143e02;      % [-]
constants.om_s = -9.25195985e-01;   % [-]
constants.om_m = 2.661861350e-06;   % [1/s]
constants.l_e = 3.84405e05;         % [km]

% Load Data:
data.alpha = 0.2*pi;            % [rad]     Alpha angle
data.beta = 1.41;               % [rad]     Beta angle
data.ti = 2;                    % [-]       Departure time
data.delta = 4;                 % [-]       Time of flight
data.tf = data.ti + data.delta; % [-]       Arrival time
data.h0 = 167;                  % [km]      Initial parking orbit height
data.hf = 100;                  % [km]      Final parking orbit height

% Compute radius & velocity of parking orbits (@ Earth and @ Moon):
data.r0 = (constants.R_Earth + data.h0)/constants.DU;
data.rf = (constants.R_Moon + data.hf)/constants.DU;
data.v0 = data.beta * sqrt((1-constants.mu)/data.r0);
data.vf = sqrt(constants.mu/data.rf);

end

function [dxdt] = PBRFBP(t, xx, constants, stm)
% ----------------------------------------------------------------------- %
% PBRFBP - Function to compute the RHS of the equations of motion for the 
% Planar Bicircular Restricted Four-Body Problem (PBRFBP) with optional 
% State Transition Matrix (STM) computation.
%
% Inputs:
%   - t: Current time [1, 1]
%   - xx: Current state vector:
%               - [4,1] if STM is not computed
%               - [20, 1] if STM is computed
%   - constants - structure containing constants of the problem:
%               - constants.mu: Mass ratio between the two primary bodies,
%               - constants.ms: Mass of the smaller third body (e.g., a spacecraft),
%               - constants.rho: Distance parameter,
%               - constants.om_s: Angular velocity of the smaller third body
%   - stm: Boolean flag (0 or 1):
%               - 0 if the STM is not computed,
%               - 1 if the STM is computed and included in `xx`.
%
% Outputs:
%   dxdt      - Derivative of the state vector:
% ----------------------------------------------------------------------- %        


% Extract parameters:
mu = constants.mu;
ms = constants.ms;
rho = constants.rho;
om_s = constants.om_s;

% Extract state variables:
x = xx(1);
y = xx(2);
vx = xx(3);
vy = xx(4);

% Compute partial derivatives of the effective potential (OM):
dOMdx = x - (ms*cos(om_s*t))/rho^2 - (mu*(mu + x - 1))/((mu + x - 1)^2 + y^2)^(3/2) + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(3/2)) - (ms*(2*x - 2*rho*cos(om_s*t)))/(2*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2));
dOMdy = y - (ms*sin(om_s*t))/rho^2 - (mu*y)/((mu + x - 1)^2 + y^2)^(3/2) - (ms*(2*y - 2*rho*sin(om_s*t)))/(2*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2)) + (y*(mu - 1))/((mu + x)^2 + y^2)^(3/2);

if stm == 0
    % If STM is not to be computed, output position and velocity EoM:
    dxdt = zeros(4, 1);
    dxdt(1:2) = xx(3:4);
    dxdt(3) = 2*vy + dOMdx;
    dxdt(4) = -2*vx + dOMdy;
elseif stm == 1
    % If STM is to be computed, include dynamics of the STM:

    % Extract STM:
    PHI = reshape(xx(5: end), 4, 4);

    % Assemble the Jacobian Matrix A:
    A = zeros(4);
    A(1, 3) = 1;
    A(2, 4) = 1;
    A(3, 1) = (mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - ms/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2) - (3*(mu + x)^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*ms*(x - rho*cos(om_s*t))^2)/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2) + (3*mu*(mu + x - 1)^2)/((mu + x - 1)^2 + y^2)^(5/2) + 1;
    A(3, 2) = (3*ms*(2*x - 2*rho*cos(om_s*t))*(2*y - 2*rho*sin(om_s*t)))/(4*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2)) + (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2));
    A(3, 4) = 2;
    A(4, 1) = (3*ms*(2*x - 2*rho*cos(om_s*t))*(2*y - 2*rho*sin(om_s*t)))/(4*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2)) + (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2));
    A(4, 2) = (mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - ms/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*ms*(y - rho*sin(om_s*t))^2)/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2)^(5/2) + 1;
    A(4, 3) = -2;

    % Compute derivative of the STM using the variational approach:
    PhiDot = A * PHI;

    % Assemble RHS:
    dxdt = zeros(20, 1);
    dxdt(1:2) = xx(3:4);
    dxdt(3) = 2*vy + dOMdx;
    dxdt(4) = -2*vx + dOMdy;
    dxdt(5:end) = PhiDot(:);
end

end

function [tt, xx] = propagate(xx0, t0, tf, constants, settings)
% ----------------------------------------------------------------------- %
% propagate - Function to propagate the s/c dynamics using the Planar 
% Bicircular Restricted Four-Body Problem (PBRFBP) model (without 
% State Transition Matrix (STM) computation)
%
% Inputs:
%   - xx0: Initial state vector [4,1] 
%   - t0: Initial time [1, 1]
%   - tf: Initial time [1, 1]
%   - constants: structure containing constants of the problem
%   - settings: structure containing settings for ode solver (settings.odeOpt)
%
% Outputs:
%   - xx: integrated discretized states [n, 4] --> i-th row is the state at
%   i-th time
%   - tt: vector of discretized times [n, 1]
% ----------------------------------------------------------------------- % 

% Settings for the ODE solver:
optset = settings.odeOpt;
% Numerical integration
[tt, xx] = ode78(@(t, x) PBRFBP(t, x, constants, 0), [t0 tf], xx0, optset);

end

function [tt, xx, PHIf] = propagate_stm(xx0, t0, tf, constants, settings)
% ----------------------------------------------------------------------- %
% propagate_stm - Function to propagate the s/c dynamics using the Planar 
% Bicircular Restricted Four-Body Problem (PBRFBP) model (including 
% State Transition Matrix (STM) computation through variational approach)
%
% Inputs:
%   - xx0: Initial state vector [4,1] 
%   - t0: Initial time [1, 1]
%   - tf: Initial time [1, 1]
%   - constants: structure containing constants of the problem
%   - settings: structure containing settings for ode solver (settings.odeOpt)
%
% Outputs:
%   - xx: integrated discretized states [n, 20] --> i-th row is the state (and STM)
%   at i-th time
%   - PHIf: STM at final time [4, 4]
%   - tt: vector of discretized times [n, 1]
% ----------------------------------------------------------------------- % 


% Initialize the State Transition Matrix (STM) at the initial time t0
Phi0 = eye(4); % Identity matrix, since STM is identity at t0

% Append the flattened STM (Phi0) to the initial conditions
xx0 = [xx0; Phi0(:)];

% Setings for the ODE solver:
optset = settings.odeOpt;
% Numerical integration:
[tt, xx] = ode78(@(t, x) PBRFBP(t, x, constants, 1), [t0 tf], xx0, optset);

% Extract the final State Transition Matrix (STM) at time tf
% Reshape the Phi vector into a 4x4 matrix:
PHIf = reshape(xx(end, 5:end), 4, 4);

end

function XX = rot2ECI(tt, xx, constants)
% ----------------------------------------------------------------------- %
% rot2eci - Function to rotate the state vector at time tt from the
% Earth-Moon rotating frame to Earth Centered Inertial frame
%
% Inputs:
%   - xx: Matrix of state vectors in the rotating frame [n, 4] 
%   - tt: vector of discretized times [n, 1]
%   - constants.mu: Gravitational parameter, representing the mass ratio for the Earth-Moon system
%
% Outputs:
%   - XX: Matrix of state vectors in the ECI frame [n, 4] 
% ----------------------------------------------------------------------- % 


% Extract mu:
mu = constants.mu;

% Transform position components from rotating to ECI frame
XX(:, 1) = (xx(:, 1) + mu) .* cos(tt) - xx(:, 2) .* sin(tt); 
XX(:, 2) = (xx(:, 1) + mu) .* sin(tt) + xx(:, 2) .* cos(tt); 

% Transform velocity components from rotating to ECI frame
XX(:, 3) = (xx(:, 3) - xx(:, 2)) .* cos(tt) - (xx(:, 4) + xx(:, 1) + mu) .* sin(tt); 
XX(:, 4) = (xx(:, 3) - xx(:, 2)) .* sin(tt) + (xx(:, 4) + xx(:, 1) + mu) .* cos(tt); 

end

function [c, ceq] = constraints(vars, data, constants) 
% ----------------------------------------------------------------------- %
% constraints - function to compute the equality and inequality constraints
% for simple shooting optimization (without gradients computation)
%
% Inputs:
%   - vars: vectors containing the variables of the problems:
%        - vars(1:4): Initial state vector [x0, y0, vx0, vy0]
%        - vars(5): Initial time t0
%        - vars(6): Final time tf
%   - data: struct containing data of the problem:
%        - data.r0: Initial orbit radius
%        - data.rf: Final orbit radius
%   - constants.mu: Mass ratio between the two primary bodies.
%
% Outputs:
%   - c: vector of inequality constraints vector [1,1]
%   - ceq: vector of equality constraints [4, 1]
% ----------------------------------------------------------------------- %


% Extract constants and data:
mu = constants.mu;  
r0 = data.r0;  
rf = data.rf;
    
% Extract Variables:
xx0 = vars(1:4);
ti = vars(5);
tf = vars(6);    
x0 = xx0(1); y0 = xx0(2);
vx0 = xx0(3); vy0 = xx0(4);  

% Propagate from ti a tf to obtain final state:
settings.odeOpt = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[~, xx] = propagate(xx0, ti, tf, constants, settings);
xFinal = xx(end, :);
xf = xFinal(1);
yf = xFinal(2);
vxf = xFinal(3);
vyf = xFinal(4);

% Compute equality constraints:
ceq1 = (x0 + mu)^2 + y0^2 - r0^2;
ceq2 = (x0 + mu)*(vx0 - y0) + y0*(vy0 + x0 + mu);
ceq3 = (xf + mu - 1)^2 + yf^2 - rf^2;
ceq4 = (xf + mu - 1)*(vxf - yf) + yf*(vyf + xf + mu - 1);

% Compute inequality constraint (time):
c = ti - tf;

% Combina tutti i vincoli di uguaglianza
ceq = [ceq1; ceq2; ceq3; ceq4];

end

function [c, ceq, dC, dCeq] = constraints2(vars, data, constants) 
% ----------------------------------------------------------------------- %
% constraints - function to compute the equality and inequality constraints
% for simple shooting optimization (with gradients computation)
%
% Inputs:
%   - vars: vectors containing the variables of the problems:
%        - vars(1:4): Initial state vector [x0, y0, vx0, vy0]
%        - vars(5): Initial time t0
%        - vars(6): Final time tf
%   - data: struct containing data of the problem:
%        - data.r0: Initial orbit radius
%        - data.rf: Final orbit radius
%   - constants.mu: Mass ratio between the two primary bodies.
%
% Outputs:
%   - c: vector of inequality constraints vector [1,1]
%   - ceq: vector of equality constraints [4, 1]
%   - dC: vector of gradients of inequality constraints
%   - dCeq: vector of gradients of equality constraints
% ----------------------------------------------------------------------- %

% Extract constants and data:
mu = constants.mu;  
r0 = data.r0;  
rf = data.rf;
    
% Extract Variables:
xx0 = vars(1:4);
ti = vars(5);
tf = vars(6);    
x0 = xx0(1); y0 = xx0(2);
vx0 = xx0(3); vy0 = xx0(4);  

% Propagate from ti a tf to obtain final state:
settings.odeOpt = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[~, xx, PHI] = propagate_stm(xx0, ti, tf, constants, settings);
xFinal = xx(end, :);
xf = xFinal(1);
yf = xFinal(2);
vxf = xFinal(3);
vyf = xFinal(4);

% Compute equality constraints:
ceq1 = (x0 + mu)^2 + y0^2 - r0^2;
ceq2 = (x0 + mu)*(vx0 - y0) + y0*(vy0 + x0 + mu);
ceq3 = (xf + mu - 1)^2 + yf^2 - rf^2;
ceq4 = (xf + mu - 1)*(vxf - yf) + yf*(vyf + xf + mu - 1);

% Compute inequality constraint (time):
c = ti - tf;

% Combina tutti i vincoli di uguaglianza
ceq = [ceq1; ceq2; ceq3; ceq4];

[dC, dCeq] = gradientConstraints_ss(vars, xFinal, PHI, constants);
end

function [dC, dCeq] = gradientConstraints_ss(vars, xFinal, PHI, constants)
% ----------------------------------------------------------------------- %
% gradientConstraints_ss - function to compute the gradients of the constraints
% for simple shooting optimization
%
% Inputs:
%   - vars: vectors containing the variables of the problems:
%        - vars(1:4): Initial state vector [x0, y0, vx0, vy0]
%        - vars(5): Initial time t0
%        - vars(6): Final time tf
%   - xFinal: final state vector
%   - PHI: final STM matrix [4, 4]
%   - constants.mu: Mass ratio between the two primary bodies.
%
% Outputs:
%   - dC: vector of gradients of inequality constraints
%   - dCeq: vector of gradients of equality constraints
% ----------------------------------------------------------------------- %


xi = vars(1); yi = vars(2); vxi = vars(3); vyi = vars(4);
ti = vars(5); tf = vars(6);
xf = xFinal(1); yf = xFinal(2); vxf = xFinal(3); vyf = xFinal(4);
mu = constants.mu;

% Derivatives of dCeq/dy, with y = {xi, yi, vxi, vyi, ti, tf}
C12_xxi = [2*(xi + mu), 2*yi, 0, 0;
          vxi, vyi, xi + mu, yi]; 
C12_ti = [0; 0];
C12_tf = [0; 0];

C34_xxf = [2*(xf + mu -1), 2*yf, 0, 0;
           vxf, vyf, xf + mu - 1, yf];
C34_xxi = C34_xxf * PHI;

xxi = vars(1:4);
f_ti = PBRFBP(ti, xxi, constants, 0);
C34_ti = -C34_xxf * PHI * f_ti;

f_tf = PBRFBP(tf, xFinal, constants, 0);
C34_tf = C34_xxf * f_tf;

dCeq = [C12_xxi, C12_ti, C12_tf;
        C34_xxi, C34_ti, C34_tf];

dC = [0 0 0 0 1 -1]';
dCeq = dCeq';

end

function dv = costFunction(vars, data, constants)
% ----------------------------------------------------------------------- %
% costFunction - function to compute the cost function for simple shooting
% optimization (without gradients computation)
%
% Inputs:
%   - vars: vectors containing the variables of the problems:
%        - vars(1:4): Initial state vector [x0, y0, vx0, vy0]
%        - vars(5): Initial time t0
%        - vars(6): Final time tf
%   - data: struct containing data of the problem:
%        - data.r0: Initial orbit radius
%        - data.rf: Final orbit radius
%   - constants.mu: Mass ratio between the two primary bodies.
%
% Outputs:
%   - dv: value of the cost function [1,1]
% ----------------------------------------------------------------------- %

% Extract constants and data:
mu = constants.mu; % Mass ratio
r0 = data.r0; % Initial orbital radius
rf = data.rf; % Final orbital radius

% Extract Variables:
xx0 = vars(1:4);
ti = vars(5);
tf = vars(6);    
x0 = xx0(1); y0 = xx0(2);
vx0 = xx0(3); vy0 = xx0(4);

% Calculate the initial deltaV:
dv1 = sqrt((vx0 - y0)^2 + (vy0 + x0 + mu)^2) - sqrt((1 - mu) / r0);

% Propagate from ti a tf to obtain final state:
settings.odeOpt = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[~, xx] = propagate(xx0, ti, tf, constants, settings);
xFinal = xx(end, :);
xf = xFinal(1);
yf = xFinal(2);
vxf = xFinal(3);
vyf = xFinal(4);

% Calculate the final deltaV:
dv2 = sqrt((vxf - yf)^2 + (vyf + xf + mu - 1)^2) - sqrt(mu / rf);

% Calculate the total deltaV:
dv = dv1 + dv2;

end

function [dv, gradDv] = costFunction2(vars, data, constants)
% ----------------------------------------------------------------------- %
% costFunction - function to compute the cost function for simple shooting
% optimization (with gradients computation)
%
% Inputs:
%   - vars: vectors containing the variables of the problems:
%        - vars(1:4): Initial state vector [x0, y0, vx0, vy0]
%        - vars(5): Initial time t0
%        - vars(6): Final time tf
%   - data: struct containing data of the problem:
%        - data.r0: Initial orbit radius
%        - data.rf: Final orbit radius
%   - constants.mu: Mass ratio between the two primary bodies.
%
% Outputs:
%   - dv: value of the cost function [1,1]
%   - gradDv: gradient of the cost function
% ----------------------------------------------------------------------- %

% Extract constants and data:
mu = constants.mu; % Mass ratio
r0 = data.r0; % Initial orbital radius
rf = data.rf; % Final orbital radius

% Extract Variables:
xx0 = vars(1:4);
ti = vars(5);
tf = vars(6);    
x0 = xx0(1); y0 = xx0(2);
vx0 = xx0(3); vy0 = xx0(4);

% Calculate the initial deltaV:
dv1 = sqrt((vx0 - y0)^2 + (vy0 + x0 + mu)^2) - sqrt((1 - mu) / r0);

% Propagate from ti a tf to obtain final state:
settings.odeOpt = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[~, xx, PHI] = propagate_stm(xx0, ti, tf, constants, settings);
xFinal = xx(end, :);
xf = xFinal(1);
yf = xFinal(2);
vxf = xFinal(3);
vyf = xFinal(4);

% Calculate the final deltaV:
dv2 = sqrt((vxf - yf)^2 + (vyf + xf + mu - 1)^2) - sqrt(mu / rf);

% Calculate the total deltaV:
dv = dv1 + dv2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GRADIENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P1 = 1/sqrt((vx0 - y0)^2 + (vy0 + x0 + mu)^2) .* [vy0 + x0 + mu; y0 - vx0; vx0 - y0; vy0 + x0 + mu];
PN = 1/sqrt((vxf - yf)^2 + (vyf + xf + mu - 1)^2) .* [vyf + xf + mu - 1; yf - vxf; vxf - yf; vyf + xf + mu - 1];
 
f_ti = PBRFBP(ti, xx0, constants, 0);
f_tf = PBRFBP(tf, xFinal, constants, 0);


DV_ti = - (PHI' * PN)' * f_ti;
DV_tf = PN' * f_tf;
 
gradDv = [P1 + PHI' * PN;
          DV_ti;
          DV_tf];

end

function [dv, gradDv] = costFunction_ms(vars, data, constants)
% ----------------------------------------------------------------------- %
% costFunction - function to compute the cost function for multiple shooting
% optimization (with gradients computation)
%
% Inputs:
%   - vars: vectors containing the variables of the problems:
%        - vars(1:16): vectors of state vectors at each node
%        - vars(17): Initial time t0
%        - vars(18): Final time tf
%   - data: struct containing data of the problem:
%        - data.r0: Initial orbit radius
%        - data.rf: Final orbit radius
%   - constants.mu: Mass ratio between the two primary bodies.
%
% Outputs:
%   - dv: value of the cost function [1,1]
%   - gradDv: gradient of the cost function
% ----------------------------------------------------------------------- %

% Extract constants and data:
mu = constants.mu;
r0 = data.r0;
rf = data.rf;
v0 = sqrt((1 - mu)/r0);
vf = sqrt(mu/rf);

% Compute number of time instants for multiple shooting
N = (length(vars) - 2) / 4;

% Extract initial state (position and velocity) and final state
xi = vars(1);      % Initial x position
yi = vars(2);      % Initial y position
vxi = vars(3);     % Initial x velocity
vyi = vars(4);     % Initial y velocity

xf = vars(end-5);  % Final x position
yf = vars(end-4);  % Final y position
vxf = vars(end-3); % Final x velocity
vyf = vars(end-2); % Final y velocity

% Compute deltaV:
DV1 = sqrt((vxi - yi)^2 + (vyi + xi + mu)^2) - v0;
DV2 = sqrt((vxf - yf)^2 + (vyf + xf + mu - 1)^2) - vf;
dv = DV1 + DV2;

P1 = 1/sqrt((vxi - yi)^2 + (vyi + xi + mu)^2) * [vyi + xi + mu; yi - vxi; vxi - yi; vyi + xi + mu];
PN = 1/sqrt((vxf - yf)^2 + (vyf + xf + mu - 1)^2) * [vyf + xf + mu - 1; yf - vxf; vxf - yf; vyf + xf + mu - 1];
% Concatenate gradients for initial and final states:
gradDv = [P1; zeros(N * 2, 1); PN; 0; 0];

end

function [c, ceq, dC, dCeq] = constraints_ms(vars, data, constants)
% ----------------------------------------------------------------------- %
% constraints - function to compute the equality and inequality constraints
% for multiple shooting optimization (with gradients computation)
%
% Inputs:
%   - vars: vectors containing the variables of the problems:
%        - vars(1:16): State vectors at each of the 4 nodes
%        - vars(17): Initial time t0
%        - vars(18): Final time tf
%   - data: struct containing data of the problem:
%        - data.r0: Initial orbit radius
%        - data.rf: Final orbit radius
%   - constants.mu: Mass ratio between the two primary bodies.
%
% Outputs:
%   - c: vector of inequality constraints vector 
%   - ceq: vector of equality constraints 
%   - dC: vector of gradients of inequality constraints
%   - dCeq: vector of gradients of equality constraints
% ----------------------------------------------------------------------- %

mu = constants.mu;  % Gravitational parameter
Re = constants.R_Earth;  % Earth radius
Rm = constants.R_Moon;  % Moon radius
DU = constants.DU;
r0 = data.r0;  % Initial orbit radius
rf = data.rf;  % Final orbit radius

% Compute number of time instants for multiple shooting
N = (length(vars) - 2) / 4;

% Extract Variables:
xi = vars(1);      % Initial x position
yi = vars(2);      % Initial y position
vxi = vars(3);     % Initial x velocity
vyi = vars(4);     % Initial y velocity

xf = vars(end-5);  % Final x position
yf = vars(end-4);  % Final y position
vxf = vars(end-3); % Final x velocity
vyf = vars(end-2); % Final y velocity

ti = vars(end-1);
tf = vars(end);

% Linearly interpolate time for each shooting interval
t = linspace(ti, tf, N);  % Interpolated time vector

% Initialize equality constraint vector for positions and velocities
ceq = zeros(4*N, 1);
% Initial orbit constraints:
ceq(end-3) = (xi + mu)^2 + yi^2 - r0^2;  
ceq(end-2) = (xi + mu) * (vxi - yi) + yi * (vyi + xi + mu);  
% Final orbit constraints:
ceq(end-1) = (xf + mu - 1)^2 + yf^2 - rf^2;  
ceq(end) = (xf + mu - 1) * (vxf - yf) + yf * (vyf + xf + mu - 1);  

% Inequality constraints:
c = zeros(2*N + 1, 1);
c(1:2, 1) = [(Re/DU)^2 - (xi + mu)^2 - yi^2; 
            (Rm/DU)^2 - (xi + mu - 1)^2 - yi^2]; 

% Final time inequality constraint:
c(end) = ti - tf;

% Initialize gradient vectors:
Q1 = zeros(4*(N-1), 1);
QN = zeros(4*(N-1), 1);

% Initialize Jacobian matrices:
dCeq = zeros(3*N + 4, 4*N + 2);
dC = zeros(2*N + 1, 4*N + 2);

% Propagate states and apply constraints for each shooting interval
for j = 1:N-1
    % Extract current and next state:
    xxjj = vars(4*j + 1 : 4*j + 4);  % Next state
    xxj = vars(4*j-3 : 4*j);  % Current state

    % Propagate state over the current interval:
    settings.odeOpt = odeset('RelTol', 1e-11, 'AbsTol', 1e-11);
    [~, xx, PHI] = propagate_stm(xxj, t(j), t(j+1), constants, settings);

    % Compute continuity constraint (matching positions and velocities):
    ceq(4*j-3 : 4*j) = xx(end, 1:4)' - xxjj;

    % Update inequality constraints:
    c(2*j+1 : 2*j+2) = [(Re/DU)^2 - (xxjj(1) + mu)^2 - xxjj(2)^2;  
                        (Rm/DU)^2 - (xxjj(1) + mu - 1)^2 - xxjj(2)^2]; 

    % Gradients:
    % Extract current interval position:
    xj = xxj(1);  
    yj = xxj(2);
    % Get RHS of the PBRFBP:
    fj = PBRFBP(t(j), xxj, constants, 0);  
    fjj = PBRFBP(t(j+1), xx(end, :), constants, 0);
    % Compute equality constraint gradients (matching positions and
    % velocities):
    Q1(4*j-3 : 4*j) = -(N-j)/(N-1) * PHI * fj + (N-j-1)/(N-1) * fjj;  
    QN(4*j-3 : 4*j) = -(j-1)/(N-1) * PHI * fj + j/(N-1) * fjj;  

    % Update Jacobians for equality constraints:
    dCeq(4*j-3: 4*j, 4*j+1:  4*j + 4) = -eye(4);  
    dCeq(4*j-3: 4*j, 4*j-3 : 4*j) = PHI;  

    % Update Jacobians for inequality constraints:
    dC(2*j-1: 2*j, 4*j-3:  4*j) = [-2*(xj + mu), -2*yj, 0, 0;  
                                   -2*(xj + mu - 1), -2*yj, 0, 0];
    
end 

% Final inequality constraint:
dC(2*N-1:2*N, 4*N-3:4*N) = [-2*(xf + mu), -2*yf, 0, 0;  % Earth boundary final Jacobian
                            -2*(xf + mu - 1), -2*yf, 0, 0];  % Moon boundary final Jacobian
dC(end, end-1:end) = [1 -1];  % Time constraint Jacobian

% Combine the Jacobian matrices
dC = dC';  % Transpose Jacobian for inequality constraints
dCeq(1:4*(N-1), end-1:end) = [Q1, QN];  % Combine continuity gradients

% Final Jacobian for equality constraints
dPHI1dxi = [2*(xi + mu), 2*yi, 0, 0; vxi, vyi, (xi + mu), yi]; 
dPHI2dxf = [2*(xf + mu - 1), 2*yf, 0, 0; vxf, vyf, (xf + mu - 1), yf];  
dPHI = [dPHI1dxi, zeros(2, 3*N+2); zeros(2, 3*N), dPHI2dxf, zeros(2,2)];  
% Add Jacobian to the equality constraint gradients:
dCeq(end-3:end, :) = dPHI;  

% Transpose the equality constraint Jacobian matrix
dCeq = dCeq';

end

function [bodies] = nBody_init(labels)
%-------------------------------------------------------------------------%
% Initialize planetary data for n-body propagation
%
% INPUTS:
%   labels : [1,n] cell-array with object labels
%
% OUTPUT:
%   bodies : [1,n] cell-array with struct elements containing:
%                  |--bodies{i}.name -> body label
%                  |--bodies{i}.GM   -> gravitational constant [km**3/s**2]
%-------------------------------------------------------------------------%

bodies = cell(size(labels));

for i = 1:length(labels)
    bodies{i}.name = labels{i};
    bodies{i}.GM = cspice_bodvrd(labels{i}, 'GM',1);
end

end

function [dxdt] = nBody_RHS(t, x, bodies, frame, center)
%-------------------------------------------------------------------------%
% nBody_RHS
%
% Function to evaluate the right-hand side (RHS) of an N-body propagator.
% If the system's center is shifted from the Solar System Barycenter (SSB),
% both Keplerian and perturbing terms are considered in the calculation.
%
% INPUTS:
%   t      : [1x1]  Time variable
%   x      : [6x1]  State vector [x; y; z; vx; vy; vz]
%   bodies : [1,n]  Cell-array with struct elements containing:
%                      |-- bodies{i}.name -> Name of the celestial body
%                      |-- bodies{i}.GM   -> Gravitational parameter [km^3/s^2]
%   frame  : [str]   Reference frame (ECLIPJ2000 or J2000)
%   center : [str]   (Optional) System's center (e.g., 'SSB'). If not provided,
%                    the default value is 'SSB'.
%
% OUTPUT:
%   dxdt   : [6x1]  RHS of the N-body propagator (derivatives of state vector)
%
%-------------------------------------------------------------------------%

% Check input reference frame:
if not( strcmpi(frame, 'ECLIPJ2000') || strcmpi(frame, 'J2000') )
    error('Invalid reference frame, select either J2000 or ECLIPJ2000');
end

% Set center to SSB if not provided as input:
if nargin < 5
    center = 'SSB';
end

dxdt = zeros(6,1);   % Initialize RHS
dxdt(1:3) = x(4:6);  % Derivative of the position is object's velocity
rr_obj = x(1:3);     % Extract object's position from state vector

% If center is not 'SSB', compute contribution to acceleration of central body
if ~strcmpi(center, 'SSB')
    GM0 = cspice_bodvrd(center, 'GM', 1);
    dist_center2 = dot(rr_obj, rr_obj);
    dist_center = sqrt(dist_center2);
    aa_grav_center = -GM0 * rr_obj / (dist_center * dist_center2);
    dxdt(4:6) = dxdt(4:6) + aa_grav_center;
end

% Loop over all celestial bodies to compute gravitational accelerations
for i = 1:length(bodies)
    % If center is not 'SSB', skip the central body (no need to compute
    % self-interaction):
        if strcmpi(bodies{i}.name, center)
            continue;
        end
    rv_body = cspice_spkezr(bodies{i}.name, t, frame, 'NONE', center); % i-th body position relative to the center (SSB or other)
    rr_body_obj = rr_obj - rv_body(1:3); % object position relative to i-th body
    d2 = dot(rr_body_obj, rr_body_obj);
    d = sqrt(d2);
    % Compute gravitational acceleration: 
    if strcmpi(center, 'SSB') % Add acceleration to RHS (if center is 'SSB'):
        aa_grav = -bodies{i}.GM * rr_body_obj / (d * d2);
        dxdt(4:6) = dxdt(4:6) + aa_grav;
    else % Add perturbing acceleration (if center is not 'SSB'):
        rho = rv_body(1:3);
        q = dot(rr_obj, (rr_obj - 2*rho))/dot(rho, rho);
        f = (q*(3+3*q+q^2))/(1+(1+q)^1.5);
        aa_pert = -bodies{i}.GM * (rr_obj + rho.*f)/(norm(d)^3);
        dxdt(4:6) = dxdt(4:6) + aa_pert;
    end
end
end

function [tt, xx] = nBodyPropagator(xx0, et0, etf, bodies, frame, center)
% ----------------------------------------------------------------------- %
% nBodyPropagator - Function to propagate the s/c dynamics using the N-Body
% model
%
% Inputs:
%   - xx0: Initial state vector [6,1] 
%   - et0: Initial time [1, 1]
%   - etf: Initial time [1, 1]
%   - bodies : [1,n]  Cell-array with struct elements containing:
%                      |-- bodies{i}.name -> Name of the celestial body
%                      |-- bodies{i}.GM   -> Gravitational parameter [km^3/s^2]
%   - frame  : [str]   Reference frame (ECLIPJ2000 or J2000)
%   - center : [str]   (Optional) System's center (e.g., 'SSB'). If not provided,
%                    the default value is 'SSB'.
%
% Outputs:
%   - xx: integrated discretized states [n, 6] --> i-th row is the state at
%   i-th time
%   - tt: vector of discretized times [n, 1]
% ----------------------------------------------------------------------- % 

% Perform Integration:
options = odeset('reltol', 1e-12, 'abstol', [1e-6*ones(1,3),1e-9*ones(1,3)]);
[tt, xx] = ode113(@(t,x) nBody_RHS(t, x, bodies, frame, center), [et0 etf], xx0, options);
end

function thetaFinderFun = thetaFinder(etVec, thetaTarget, frame, center)
% ----------------------------------------------------------------------- %
% thetaFinder - function to calculate the angular difference (theta) between the Sun 
% and the Moon relative to a central body (Earth or another planet) in a 
% specified reference frame. It then computes the difference between the 
% computed angular difference and a target angular difference (theta_target).
%
% Inputs:
%   - thetaTarget: Target angular difference between the Sun and Moon [rad]
%   - etVec: array of ephemeris times [s]
%   - bodies: Cell array containing the structs for each body (e.g., Sun, Earth, Moon)
%   - frame: Reference frame for the positions (e.g., 'J2000', 'EME2000')
%   - center: Central body for the positions (e.g., 'EARTH')
%
% Outputs:
%   thetaFinderFun - Difference between the computed theta angle and the target angle
% ----------------------------------------------------------------------- %

M_pos = zeros(3, length(etVec));
S_pos = zeros(3, length(etVec));

for i = 1 : length(etVec)
    M_pos(:, i) = cspice_spkpos('MOON', etVec(i), frame, 'none', center);
    S_pos(:, i) = cspice_spkpos('SUN', etVec(i), frame, 'none', center);
end

thetaSun = atan2(S_pos(2, :), S_pos(1, :));
thetaMoon = atan2(M_pos(2, :), M_pos(1, :));

thetaFinderFun = mod(thetaSun - thetaMoon, 2*pi) - thetaTarget;
end

function [x, y] = circularOrbit(radius, center)
% ----------------------------------------------------------------------- %
% Function to compute x,y coordinates for a circular orbit given center and
% radius
% 
% Inputs:
%   - center: center of the orbit [2, 1]
%   - radius: radius of the orbit [1, 1]
%
% Outputs:
%   - x: array of x coordinates [1000, 1]
%   - y: array of y coordinates [1000, 1]
% ----------------------------------------------------------------------- %

% Define theta vector:
theta = linspace(0, 2*pi, 1000); 
% Extract coordinates of center:
xC = center(1);
yC = center(2);
% Compute coordinate vectors of the circular orbit:
x = xC + radius * cos(theta);  
y= yC + radius * sin(theta); 

end

function plotSettings
%-------------------------------------------------------------------------%
% Settings for figure plots
%-------------------------------------------------------------------------%
% Setting Lines:
set(0, 'defaultLineLineWidth', 1.6);
set(0,'defaultLineMarkerSize', 4) ;
set(0,'defaultLineMarkerFaceColor', 'auto')
set(0, 'defaultLineMarkerEdgeColor', 'auto')
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

function plotMultipleShooting(vars_opt_ms, N, constants, data)
% ----------------------------------------------------------------------- %
% Function to plot Multiple Shooting optimal trajectory in the rotating
% frame and ECI frame
% 
% Inputs:
%   - vars_opt_ms: optimal multiple shooting variables [18, 1]
%   - N: number of nodes [1, 1]
%   - data: struct containing data of the problem:
%        - data.r0: Initial orbit radius
%        - data.rf: Final orbit radius
%   - constants.mu: Mass ratio between the two primary bodies.
%
% ----------------------------------------------------------------------- %

mu = constants.mu;
r0 = data.r0;
rf = data.rf;
% Initialize cell arrays for time and state vectors:
time_intervals = cell(N-1, 1); 
state_intervals = cell(N-1, 1); 

% Extract initial and final times from optimization results
ti = vars_opt_ms(end-1);
tf = vars_opt_ms(end);

% Compute time vector for each interval
t = linspace(ti, tf, N);

% Loop over intervals to propagate states and store results
for j = 1:N-1
    % Propagate the state for the current interval
    settings.odeOpt = odeset(RelTol=1e-11, AbsTol=1e-11);
    [tt, xx] = propagate(vars_opt_ms(4*j-3:4*j), t(j), t(j+1), constants, settings); 
    % Store propagated times and states:
    time_intervals{j} = tt;   
    state_intervals{j} = xx(:, 1:4); 
end

% Concatenate results from all intervals:
tvec = vertcat(time_intervals{:}); 
xx_ms = vertcat(state_intervals{:});   

% Figure 1: Optimal trajectory in the rotating frame
figure();
hold on;
grid on;
plot(xx_ms(:, 1), xx_ms(:, 2), 'k');
[xPark_i, yPark_i] = circularOrbit(r0, [-mu; 0]); % initial parking orbit
[xPark_f, yPark_f] = circularOrbit(rf, [1-mu; 0]);  % final parking orbit
plot(xPark_i, yPark_i, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.7) 
plot(xPark_f, yPark_f, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.7)
for j = 1:floor((length(vars_opt_ms)) / 4) - 1
    plot(vars_opt_ms(4 * j + 1), vars_opt_ms(4 * j + 2), 'Marker', 'o', 'MarkerSize', 7, 'LineWidth', 2);
end
plot(vars_opt_ms(1), vars_opt_ms(2), 'Marker', 'o', 'MarkerSize', 7, 'LineWidth', 2);
title('Optimal Multiple Shooting Trajectory (@EMB Earth-Moon Rotating Frame)');
legend('Transfer Orbit', 'Parking Orbit', 'Target Orbit', 'Nodes', 'FontSize', 15);
axis equal;
xlim([-1.3 1.45])
ylim([-0.1 1.7])
xlabel('x [-]', 'FontSize', 23)
ylabel('y [-]', 'FontSize', 23)
set(gca, 'FontSize', 16)
legend('FontSize', 18);

% Rotate to ECI:
xx_ms = xx_ms(:,1:4);
xxECI = rot2ECI(tvec, xx_ms, constants);

% Figure 2: Optimal trajectory in the ECI frame
figure()
hold on;
grid on;
plot(xxECI(:, 1), xxECI(:, 2), 'k', 'LineWidth', 1);
[xParkECI, yParkECI] = circularOrbit(r0, [0; 0]);
[xMoon, yMoon] = circularOrbit(1, [0; 0]);
plot(xParkECI, yParkECI, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.7)
plot(xMoon, yMoon, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.7) 
for i = 1:N - 1
    optStateEciMS = rot2ECI(time_intervals{i}, state_intervals{i}, constants);
    plot(optStateEciMS(:, 1), optStateEciMS(:, 2), 'k');
    plot(optStateEciMS(1, 1), optStateEciMS(1, 2), 'o');
end
title('Optimal Multiple Shooting Trajectory (@Earth ECI)');
legend('Transfer Orbit', 'Parking Orbit', 'Moon Orbit', '', '', '', 'Nodes', 'FontSize', 15);
axis equal;
xlim([-1.7 1.72])
ylim([-1.2 1.22])
xlabel('x [-]', 'FontSize', 23)
ylabel('y [-]', 'FontSize', 23)
set(gca, 'FontSize', 16)
legend('FontSize', 18);

end

