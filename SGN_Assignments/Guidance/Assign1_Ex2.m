% Spacecraft Guidance and Navigation - A.Y. 2024/25
% Assignment #1
% Exercise #2
% Author: Frassinella Luca - 10795356

%% -------------------- EX. 2 - IMPULSIVE GUIDANCE --------------------- %%

clearvars; close all; clc;
cspice_kclear()
rng default
format long
plotSettings;
addpath('.\kernels\')
cspice_furnsh('assignment01.tm')

%% -------------------- Pt.1 - Initial Guess Solution -------------------- 

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

% Construct initial state vector from data:
xx0 = [r0*cos(alpha) - mu; r0*sin(alpha); -(v0 - r0)*sin(alpha); (v0 - r0)*cos(alpha)];

% Propagate and Plot initial guess solution:
[tt, xx] = ode78(@PBRFBP, [ti, tf], xx0, settings.odeOpt, constants);


figure()
plot(xx(:, 1), xx(:, 2), 'k', 'DisplayName', 'Trajectory')
grid on
hold on
axis equal
[xPark_i, yPark_i] = circularOrbit(r0, [-mu; 0]); % initial parking orbit
[xPark_f, yPark_f] = circularOrbit(rf, [1-mu; 0]);  % final parking orbit
plot(xPark_i, yPark_i, '--k', 'LineWidth', 0.8, 'DisplayName', 'Initial parking orbit') 
plot(xPark_f, yPark_f, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6, 'DisplayName', 'Final parking orbit')
% plot(-mu, 0, 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0 0.2 1], 'DisplayName', 'Earth')
% plot(1-mu, 0, 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.5 0.5 0.5], 'DisplayName', 'Moon')
xlabel('x [-]')
ylabel('y [-]')
title('Trajectory from initial guess', 'FontWeight', 'bold')
subtitle('[@EMB Earth-Moon Rotating Frame]', 'FontSize', 18)
xlim([-1.5 3])
legend;

% Rotate to ECI and Plot initial guess solution:
XX = rot2ECI(tt, xx, mu);

figure()
plot(XX(:, 1), XX(:, 2), 'k', 'DisplayName', 'Trajectory')
hold on
[xParkECI, yParkECI] = circularOrbit(r0, [0; 0]);
[xMoon, yMoon] = circularOrbit(1, [0; 0]);
plot(xParkECI, yParkECI, '--k', 'LineWidth', 0.8, 'DisplayName', 'Initial parking orbit')
plot(xMoon, yMoon, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6, 'DisplayName', 'Moon orbit') 
% plot(0, 0, 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0 0.2 1], 'DisplayName', 'Earth')
grid on
axis equal
xlabel('x [-]')
ylabel('y [-]')
title('Trajectory from initial guess', 'FontWeight', 'bold')
subtitle('[@Earth ECI]', 'FontSize', 18)
xlim([-1.25 2.5])
legend;

%% ---------------- Pt.2a - Simple Shooting Optimization ----------------- 

% No gradients of constraints and objective function:
settings.switchGradients = false;

% Initial guess {xi, ti, tf}:
vars0 = [xx0; ti; tf];
% Non-linear constraints:
nonlincon = @(vars) nonlconSS(vars, constants, data, settings);
% Cost function:
J = @(vars) objFunSS(vars, constants, data, settings);

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'active-set', ...
    'ConstraintTolerance', 1e-10, 'FiniteDifferenceType', 'central', .... 
    'SpecifyConstraintGradient', settings.switchGradients, 'SpecifyObjectiveGradient', settings.switchGradients, ...
    'UseParallel', true, 'MaxFunctionEvaluations', 2000, 'MaxIterations', 500);

% Solve optimization problem:
[varsOptSS_A, DV_SS_A] = fmincon(J, vars0, [], [], [], [], [], [], nonlincon, options);
% Check constraints compliance of the solution:
[c_A, ceq_A] = nonlconSS(varsOptSS, costants);


%% ---------------- Pt.2b - Simple Shooting Optimization ----------------- 
clc
% Initial guess {xi, ti, tf}:
vars0 = [xx0; ti; tf];
% Non-linear constraints:
nonlincon = @(vars) nonlconSS(vars, constants, data);
% Cost function:
% J = @(vars) objFunSS(vars, constants, data, settings);
% J = @(vars) objFunSS(vars, constants, data, settings);

% No gradients of constraints and objective function:
settings.switchGradients = true;
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'active-set', ...
    'ConstraintTolerance', 1e-10, 'FiniteDifferenceType', 'central', .... 
    'SpecifyConstraintGradient', settings.switchGradients, 'SpecifyObjectiveGradient', settings.switchGradients, ...
    'UseParallel', true);

% Solve optimization problem:
[varsOptSS_B, DV_SS_B] = fmincon(@(vars) objFunSS(vars, constants, data, settings), vars0, [], [], [], [], [], [], nonlincon, options);
% Check constraints compliance of the solution:
[c_B, ceq_B] = nonlconSS(varsOptSS_B, constants, data);

%% FUNCTIONS:

function [data, constants] = loadDataset()

% Load Constants:
GM_Earth = cspice_bodvrd('EARTH', 'GM', 1);
GM_Moon = cspice_bodvrd('MOON', 'GM', 1);
R_Earth = cspice_bodvrd('EARTH', 'RADII', 3);
R_Moon = cspice_bodvrd('MOON', 'RADII', 3);

constants.mu = GM_Moon/(GM_Earth + GM_Moon); % Earth-Moon system gravity constant
constants.TU = 4.34811305;      % [days]    Time unit
constants.DU = 3.84405000e05;   % [km]      Distance unit
constants.VU = 1.02454018;      % [km/s]    Velocity unit

% Earth, Sun & Moon constants:
constants.R_Earth = R_Earth(1);     % [km]      Earth radius
constants.R_Moon = R_Moon(1);       % [km]      Moon radius
constants.ms = 3.28900541e05;       % [??]
constants.rho = 3.88811143e02;      % [??]
constants.om_s = -9.25195985e-01;   % [??]
constants.om_m = 2.661861350e-06;   % [??]

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

function [x, y] = circularOrbit(radius, center)

% Define theta vector:
theta = linspace(0, 2*pi, 1000); 
% Extract coordinates of center:
xC = center(1);
yC = center(2);
% Compute coordinate vectors of the circular orbit:
x = xC + radius * cos(theta);  
y= yC + radius * sin(theta); 

end

function [tt, xx, STM] = propagator(xx0, ti, tf, constants)

% Integrate:
optset = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[tt, xxv] = ode78(@(t,xx) PBRFBP(t, xx, constants), [ti tf], xx0, optset);

% Extract computations:
if size(xxv, 2) == 20
    STM = reshape(xxv(end, 5:end), 4, 4);
    xx = xxv(:, 1:4);
elseif size(xxv, 2) == 4
    STM = [];
    xx = xxv(:, 1:4);
end

end

function [dxdt] = PBRFBP(t, xx, constants)

% Define Constants:
mu = constants.mu;
ms = constants.ms;
rho = constants.rho;
om_s = constants.om_s;

x = xx(1);
y = xx(2);
vx = xx(3);
vy = xx(4);

% Derivatives of potential function:
dOMdx = x - (ms*cos(om_s*t))/rho^2 - (mu*(mu + x - 1))/((mu + x - 1)^2 + y^2)^(3/2) + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(3/2)) - (ms*(2*x - 2*rho*cos(om_s*t)))/(2*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2));
dOMdy = y - (ms*sin(om_s*t))/rho^2 - (mu*y)/((mu + x - 1)^2 + y^2)^(3/2) - (ms*(2*y - 2*rho*sin(om_s*t)))/(2*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2)) + (y*(mu - 1))/((mu + x)^2 + y^2)^(3/2);

% Derivative of STM
if length(xx) == 20
    dxdt = zeros(20, 1);
    Phi = reshape(xx(5: end), 4, 4);

    % A matrix:
    A = zeros(4);
    A(1, 3) = 1;
    A(2, 4) = 1;
    A(3, 1) = (mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - ms/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2) - (3*(mu + x)^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*ms*(x - rho*cos(om_s*t))^2)/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2) + (3*mu*(mu + x - 1)^2)/((mu + x - 1)^2 + y^2)^(5/2) + 1;
    A(3, 2) = (3*ms*(2*x - 2*rho*cos(om_s*t))*(2*y - 2*rho*sin(om_s*t)))/(4*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2)) + (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2));
    A(3, 4) = 2;
    A(4, 1) = (3*ms*(2*x - 2*rho*cos(om_s*t))*(2*y - 2*rho*sin(om_s*t)))/(4*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2)) + (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2));
    A(4, 2) = (mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - ms/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*ms*(y - rho*sin(om_s*t))^2)/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2)^(5/2) + 1;
    A(4, 3) = -2;

    dPhi = A * Phi;
    dxdt(5:end) = reshape(dPhi, 16, 1);
    % Equations of motion:
    dxdt(1:2) = xx(3:4);
    dxdt(3) = 2*vy + dOMdx;
    dxdt(4) = -2*vx + dOMdy;
elseif length(xx) == 4
    dxdt = zeros(4, 1);
    % Equations of motion:
    dxdt(1:2) = xx(3:4);
    dxdt(3) = 2*vy + dOMdx;
    dxdt(4) = -2*vx + dOMdy;
end

% Transpose to ensure a column vector
% dxdt = dxdt';

end

function XX = rot2ECI(tt, xx, mu)

XX(:, 1) = (xx(:, 1) + mu) .* cos(tt) - xx(:, 2) .* sin(tt);
XX(:, 2) = (xx(:, 1) + mu) .* sin(tt) + xx(:, 2) .* cos(tt);
XX(:, 3) = (xx(:, 3) -xx(:, 2)) .* cos(tt) - (xx(:, 4) + xx(:,1) + mu).* sin(tt);
XX(:, 4) = (xx(:, 3) -xx(:, 2)) .* sin(tt) + (xx(:, 4) + xx(:,1) + mu).* cos(tt);

end

function [J, gradJ] = objFunSS(vars, constants, data, settings)

% Extract Settings:
switchGrads = settings.switchGradients;
% odeOpt = settings.odeOpt;

% Extract initial guess:
xxi = vars(1:4);
ti = vars(5);
tf = vars(6);

% Extract constants:
mu = constants.mu;
r0 = data.r0;
rf = data.rf;
% Initial & Final parking orbits velocities:
vi = sqrt((1 - mu)/r0);
vf = sqrt(mu/rf);

% Compute first Delta V:
xi = xxi(1);
yi = xxi(2);
vxi = xxi(3);
vyi = xxi(4);
DV1 = sqrt((vxi - yi)^2 + (vyi + xi + mu)^2) - vi;

% Propagate orbit to tf:
% [~, xx] = ode78(@PBRFBP, [ti, tf], xxi, odeOpt, constants);
if switchGrads == true
    phi0 = eye(4);
    xx0 = [xxi; phi0(:)];
elseif switchGrads == false
    xx0 = xxi;
end

[~, xx, STM] = propagator(xx0, ti, tf, constants);

finalState = xx(end, :)';
xf = finalState(1);
yf = finalState(2);
vxf = finalState(3);
vyf = finalState(4);
% Compute second Delta V:
DV2 = sqrt((vxf - yf)^2 + (vyf + xf + mu - 1)^2) - vf;

% Total Delta V:
J = DV1 + DV2;

% Gradients:
if switchGrads == true
    Phi = STM;
    P1 = 1/sqrt((vxi - yi)^2 + (vyi + xi + mu)^2) .* [vyi + xi + mu; yi - vxi; vxi - yi; vyi + xi + mu];
    PN = 1/sqrt((vxf - yf)^2 + (vyf + xf + mu - 1)^2) .* [vyf + xf + mu - 1; yf - vxf; vxf - yf; vyf + xf + mu - 1];
    f_ti = PBRFBP(ti, xxi, constants);
    f_tf = PBRFBP(tf, finalState, constants);
    DV_xxi = P1 + Phi' * PN;
    DV_ti = - (Phi' * PN)' * f_ti;
    DV_tf = PN' * f_tf;
     
    gradJ = [DV_xxi;
             DV_ti;
             DV_tf];
elseif switchGrads == false

end

end

function [c, ceq, dC, dCeq] = nonlconSS(vars, constants, data, settings)

% Extract Settings:
switchGrads = settings.switchGradients;
% odeOpt = settings.odeOpt;

% Extract initial guess:
xxi = vars(1:4);
ti = vars(5);
tf = vars(6);

% Extract constants:
mu = constants.mu;
r0 = data.r0;
rf = data.rf;

% Extract initial state:
xi = xxi(1);
yi = xxi(2);
vxi = xxi(3);
vyi = xxi(4);

if switchGrads == true
    phi0 = eye(4);
    xx0 = [xxi; phi0(:)];
elseif switchGrads == false
    xx0 = xxi;
end

% Propagate orbit to tf:
% [~, xx] = ode78(@PBRFBP, [ti, tf], xxi, odeOpt, constants);
[~, xx, STM] = propagator(xx0, ti, tf, constants);

% Extract final state:
finalState = xx(end, :)';
xf = finalState(1);
yf = finalState(2);
vxf = finalState(3);
vyf = finalState(4);

% Equality constraints:
ceq1 = (xi + mu)^2 + yi^2 - r0^2;
ceq2 = (xi + mu)*(vxi - yi) + yi*(vyi + xi + mu);
ceq3 = (xf + mu - 1)^2 + yf^2 - rf^2;
ceq4 = (xf + mu - 1)*(vxf - yf) + yf*(vyf + xf + mu - 1);
ceq = [ceq1; ceq2; ceq3; ceq4];

% Inequality constraint:
c = ti - tf;

% Gradients:
if switchGrads == true
    C12_xxi = [2*(xi + mu), 2*yi, 0, 0;
              vxi, vyi, xi + mu, yi]; 
    C12_ti = [0; 0];
    C12_tf = [0; 0];
    
    C34_xxf = [2*(xf + mu -1), 2*yf, 0, 0;
               vxf, vyf, xf + mu - 1, yf];
    C34_xxi = C34_xxf * STM;
    
    xxi = vars(1:4);
    f_ti = PBRFBP(ti, xxi, constants, false);
    C34_ti = -C34_xxf * STM * f_ti;
    f_tf = PBRFBP(tf, finalState, constants, false);
    C34_tf = C34_xxf * f_tf;
    
    dCeq = [C12_xxi, C12_ti, C12_tf;
            C34_xxi, C34_ti, C34_tf];
    
    dC = [0 0 0 0 1 -1]';
    dCeq = dCeq';
elseif switchGrads == false
end

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
