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
cspice_furnsh('ex02.tm');

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

% Propagate and Plot initial guess solution:
[tt, xx] = propagate(xx0, ti, tf, constants);

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
xlabel('x [-]')
ylabel('y [-]')
title('Trajectory from initial guess', 'FontWeight', 'bold')
subtitle('[@EMB Earth-Moon Rotating Frame]', 'FontSize', 18)
xlim([-1.5 3])
legend;

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
title('Trajectory from initial guess', 'FontWeight', 'bold')
subtitle('[@Earth ECI]', 'FontSize', 18)
xlim([-1.25 2.5])
legend;

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
[tt, xx] = propagate(vars_opt(1:4), vars_opt(5), vars_opt(6), constants);
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
xlabel('x [-]')
ylabel('y [-]')
title('Simple Shooting Optimized Trajectory (w/o gradients)', 'FontWeight', 'bold')
subtitle('[@EMB Earth-Moon Rotating Frame]', 'FontSize', 18)
xlim([-1.5 3])
legend;

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
title('Simple Shooting Optimized Trajectory (w/o gradients)', 'FontWeight', 'bold')
subtitle('[@Earth ECI]', 'FontSize', 18)
xlim([-1.25 2.5])
legend;


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
[c, ceq] = nonlCon2(vars_opt2);
% Dimensionalize errors [m], [m/s]:
errDim = [ceq(3) * constants.DU*1e3; ceq(4) * constants.VU*1e3];

results.ex2b.vars_opt = vars_opt2;
results.ex2b.dv_opt = dv_opt2;
results.ex2b.errPos = errDim(1);
results.ex2b.errVel = errDim(2);

% Plots:
% Propagate trajectory and rotate it to ECI:
[tt, xx] = propagate(vars_opt2(1:4), vars_opt2(5), vars_opt2(6), constants);
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
xlabel('x [-]')
ylabel('y [-]')
title('Simple Shooting Optimized Trajectory (with gradients)', 'FontWeight', 'bold')
subtitle('[@EMB Earth-Moon Rotating Frame]', 'FontSize', 18)
xlim([-1.5 3])
legend;

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
title('Simple Shooting Optimized Trajectory (with gradients)', 'FontWeight', 'bold')
subtitle('[@Earth ECI]', 'FontSize', 18)
xlim([-1.25 2.5])
legend;

%% Functions

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

function [dxdt] = PBRFBP(t,xx, constants, stm)

% Inputs: xx = [20x1]
% Outputs: dxdt = [20x1]

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

function [tt, xx] = propagate(xx0, t0, tf, constants)

% Setings for the ODE solver:
optset = odeset('reltol', 1e-12, 'abstol', 1e-12);
% Numerical integration
[tt, xx] = ode78(@(t, x) PBRFBP(t, x, constants, 0), [t0 tf], xx0, optset);

end

function [tt, xx, PHIf] = propagate_stm(xx0, t0, tf, constants)

% Initialize the State Transition Matrix (STM) at the initial time t0
Phi0 = eye(4); % Identity matrix, since STM is identity at t0

% Append the flattened STM (Phi0) to the initial conditions
xx0 = [xx0; Phi0(:)];

% Setings for the ODE solver:
optset = odeset('reltol', 1e-12, 'abstol', 1e-12);
% Numerical integration:
[tt, xx] = ode78(@(t, x) PBRFBP(t, x, constants, 1), [t0 tf], xx0, optset);

% Extract the final State Transition Matrix (STM) at time tf
% The last row of xx contains [x, y, vx, vy, Phi(:)]
% Reshape the flattened Phi vector into a 4x4 matrix
PHIf = reshape(xx(end, 5:end), 4, 4);

end

function XX = rot2ECI(tt, xx, constants)

% Extract mu:
mu = constants.mu;

% Transform position components from rotating to ECI frame
XX(:, 1) = (xx(:, 1) + mu) .* cos(tt) - xx(:, 2) .* sin(tt); 
XX(:, 2) = (xx(:, 1) + mu) .* sin(tt) + xx(:, 2) .* cos(tt); 

% Transform velocity components from rotating to ECI frame
XX(:, 3) = (xx(:, 3) - xx(:, 2)) .* cos(tt) - (xx(:, 4) + xx(:, 1) + mu) .* sin(tt); 
XX(:, 4) = (xx(:, 3) - xx(:, 2)) .* sin(tt) + (xx(:, 4) + xx(:, 1) + mu) .* cos(tt); 

end
% IN TEORIA BASTA AGGIUNGERE stm COME VARIABILE A TUTTE LE FUNZIONI E UNIRE
% IL CASO CON E IL CASO SENZA GRADIENTI
function [c, ceq] = constraints(vars, data, constants) 
    % constraints Defines the equality and inequality constraints for an
    % optimization problem in the Circular Restricted 4-Body Problem (CR4BP).

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
[~, xx] = propagate(xx0, ti, tf, constants);
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
    % constraints Defines the equality and inequality constraints for an
    % optimization problem in the Circular Restricted 4-Body Problem (CR4BP).

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
[~, xx, PHI] = propagate_stm(xx0, ti, tf, constants);
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
[~, xx] = propagate(xx0, ti, tf, constants);
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
[~, xx, PHI] = propagate_stm(xx0, ti, tf, constants);
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
