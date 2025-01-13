% EXCERCISE 2 - IMPULSIVE GUIDANCE
clearvars; close all; clc;
%% -------------------- EX. 2 - IMPULSIVE GUIDANCE ---------------------- %%

clearvars; close all; clc;
format long

plotSettings;
cspice_kclear()

% -------------------- Pt.1 - Initial Guess Solution -------------------- %
tic
% LOAD KERNELS:
cspice_furnsh('ex02.tm')
% Extract constants from cspice
constants = cspiceConstants('EARTH', 'MOON');
mu = constants.mu;
R_Earth = constants.R_body1;
R_Moon = constants.R_body2;

% DATA:
alpha = 0.2*pi;
beta = 1.41;
ti = 2;
delta = 4; 
tf = ti + delta;
h0 = 167;
hf = 100;
TU = 4.34811305; % [days]
DU = 3.84405000e05; % [km]
VU = 1.02454018;    % [km/s]
l_em = 3.84405e05; % [km]


% Compute radius & velocity of parking orbits (@ Earth and @ Moon):
r0 = (R_Earth + h0)/DU;
rf = (R_Moon + hf)/DU;
v0 = beta * sqrt((1-mu)/r0);
vf = sqrt(mu/rf);


x0 = r0*cos(alpha) - mu;
y0 = r0*sin(alpha);
vx0 = -(v0 - r0)*sin(alpha);
vy0 = (v0 - r0)*cos(alpha);

xx0 = [x0; y0; vx0; vy0];
[tt, xx] = propagator(xx0, ti, tf, mu);

figure()
plot(xx(:, 1), xx(:, 2), 'DisplayName', 'Trajectory')
grid on
hold on
axis equal
% Plot circumferences of r0 and rf centered in P1 (Earth) and P2 (moon)
theta = linspace(0, 2*pi, 801);  
xPark_i = -mu + r0 * cos(theta);  yPark_i = r0 * sin(theta); 
xPark_f = (1 - mu) + rf * cos(theta);  yPark_f = rf * sin(theta); 
plot(xPark_i, yPark_i, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8, 'DisplayName', 'Initial parking orbit') 
plot(xPark_f, yPark_f, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.6, 'DisplayName', 'Final parking orbit')
% plot(-mu, 0, 'o', 'MarkerFaceColor', [0 0.2 1], 'DisplayName', 'Earth')
% plot(1-mu, 0, 'o', 'MarkerFaceColor', [0.5 0.5 0.5], 'DisplayName', 'Moon')
xlabel('x [-]')
ylabel('y [-]')
title('Trajectory from initial guess', 'FontWeight', 'bold')
subtitle('[@EMB Earth-Moon Rotating Frame]', 'FontSize', 15)
legend;

% Rotate to ECI and Plot initial guess solution:
XX = rot2ECI(tt, xx, mu);

figure()
plot(XX(:, 1), XX(:, 2), 'DisplayName', 'Trajectory')
hold on
% Plot Moon Orbit:
theta = linspace(0, 2*pi, 801); rMoon = 1;
xParkECI = r0 * cos(theta);  yParkECI = r0 * sin(theta); 
xMoon = -mu + rMoon*cos(theta); yMoon = rMoon*sin(theta);
plot(xParkECI, yParkECI, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8, 'DisplayName', 'Initial parking orbit')
plot(xMoon, yMoon, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.6, 'DisplayName', 'Moon orbit') 
% plot(0, 0, 'o', 'MarkerFaceColor', [0 0.2 1], 'DisplayName', 'Earth')
grid on
axis equal
xlabel('x [-]')
ylabel('y [-]')
title('Trajectory from initial guess', 'FontWeight', 'bold')
subtitle('[@Earth ECI]', 'FontSize', 15)
legend;
%% 2a) OPTIMIZATION - SIMPLE SHOOTING
clc; close;
mu = constants.mu; constants.r0 = r0; constants.rf = rf;
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'active-set', 'ConstraintTolerance', 1e-10, 'OptimalityTolerance', 1e-10, 'SpecifyConstraintGradient', false, 'SpecifyObjectiveGradient', false, 'FiniteDifferenceStepSize', 1e-12, UseParallel=true);
vars0 = [xx0; ti; tf];

A = [0 0 0 0 1 -1];
b = 0;
% TOGLIERE c DAI NONLINCONS
cost_function = @(vars) costFunction(vars, constants);
nonlincon = @(vars) nonlcon(vars, constants);
[vars_opt, DVtot_singleShoot_2a] = fmincon(cost_function, vars0, A, b, [], [], [], [], nonlincon, options);
[cSS_2a, ceqSS_2a] = nonlincon(vars_opt);

% Plot Solution:
[tt, xx] = propagator(vars_opt(1:4), vars_opt(5), vars_opt(6), mu);

figure()
plot(xx(:, 1), xx(:, 2))
grid on
hold on
axis equal
% Plot circumferences of r0 and rf centered in P1 (Earth) and P2 (moon)
theta = linspace(0, 2*pi, 801);  
xEarth = -mu + r0 * cos(theta);  yEarth = r0 * sin(theta); 
xMoon = (1 - mu) + rf * cos(theta);  yMoon = rf * sin(theta); 
plot(xEarth, yEarth, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8, 'HandleVisibility', 'off')  % Grigio scuro e tratteggiata
plot(xMoon, yMoon, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.6, 'HandleVisibility', 'off')  % Grigio scuro e tratteggiata

% ROTATED PLOT:
XX = rot2ECI(tt, xx, mu);
figure()
hold on
plot(XX(:, 1), XX(:, 2))
% Plot Moon Orbit:
theta = linspace(0, 2*pi, 801); rMoon = l_em/DU;
xMOON = -mu + rMoon*cos(theta); yMOON = rMoon*sin(theta);
plot(xMOON, yMOON, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.6, 'HandleVisibility', 'off')  % Grigio scuro e tratteggiata

grid on
axis equal

%% 2b) OPTIMIZATION - SIMPLE SHOOTING with function gradients:
tic
clc; close;
mu = constants.mu; constants.r0 = r0; constants.rf = rf;
% options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'active-set', 'ConstraintTolerance', 1e-10, 'OptimalityTolerance', 1e-12);
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'active-set', 'ConstraintTolerance', 1e-10, 'OptimalityTolerance', 1e-10, 'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, 'UseParallel', true);

vars0 = [xx0; ti; tf];
vPark = sqrt((1 - mu)/r0);

%
cost_function = @(vars) costFunction(vars, constants);
nonlincon = @(vars) nonlcon(vars, constants);
[vars_opt, DVtot_singleShoot_2b] = fmincon(cost_function, vars0, [], [], [], [], [], [], nonlincon, options);
[cSS_2b, ceqSS_2b] = nonlincon(vars_opt);
toc

%%
% Plot Solution:
[tt, xx] = propagator(vars_opt(1:4), vars_opt(5), vars_opt(6), mu);

figure()
plot(xx(:, 1), xx(:, 2))
grid on
hold on
axis equal
% Plot circumferences of r0 and rf centered in P1 (Earth) and P2 (moon)
theta = linspace(0, 2*pi, 801);  
xEarth = -mu + r0 * cos(theta);  yEarth = r0 * sin(theta); 
xMoon = (1 - mu) + rf * cos(theta);  yMoon = rf * sin(theta); 
plot(xEarth, yEarth, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8, 'HandleVisibility', 'off')  % Grigio scuro e tratteggiata
plot(xMoon, yMoon, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.6, 'HandleVisibility', 'off')  % Grigio scuro e tratteggiata

% ROTATED PLOT:
XX = rot2ECI(tt, xx, mu);
figure()
plot(XX(:, 1), XX(:, 2))
hold on
% Plot Moon Orbit:
theta = linspace(0, 2*pi, 801); rMoon = l_em/DU;
xMOON = -mu + rMoon*cos(theta); yMOON = rMoon*sin(theta);
plot(xMOON, yMOON, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.6, 'HandleVisibility', 'off')  % Grigio scuro e tratteggiata

grid on
axis equal

%% 3) OPTIMIZATION - MULTIPLE SHOOTING
clc;
tic
% Setting Multiple Shooting:
N = 4;

% Create initial guesses at each node:
[vars0] = initialGuessMS(xx0, ti, tf, N, constants);

mu = constants.mu; constants.r0 = r0; constants.rf = rf;
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'active-set', 'ConstraintTolerance', 1e-07, 'MaxFunctionEvaluations', 20000, 'MaxIterations', 6000, 'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, 'UseParallel', true);
cost_functionMS = @(vars) costFunctionMS(vars, constants, N);
[value, error] = checkGradients(cost_functionMS, vars0, 'Display', 'on');
nonlinconMS = @(vars) nonlconMS(vars, constants, N);
% nonlinconMS = @(vars) nonLinearConstraints_MS(vars, constants, N);
[varsOpt, DVtot_multiShoot] = fmincon(cost_functionMS, vars0, [], [], [], [], [], [], nonlinconMS, options);
[cMS, ceqMS] = nonlinconMS(varsOpt);

%% ERRORS AT EACH NODE:
% Position errors in [m]:
errorsPosition = [ceqMS(1:2); ceqMS(5:6); ceqMS(9:10); ceqMS(13); ceqMS(15)] .* DU * 1000; 
% Velocity errors in [m/s]: 
errorsVelocity = [ceqMS(3:4); ceqMS(7:8); ceqMS(11:12); ceqMS(14),; ceqMS(16)] .* VU * 1000; 
%%
toc
% Plot Solution:
plotMultipleShooting(varsOpt, DVtot_multiShoot, constants, N)

%% ROTATED PLOT:
XX = rot2ECI(tt, xx, mu);
figure()
plot(XX(:, 1), XX(:, 2))
hold on
% Plot Moon Orbit:
theta = linspace(0, 2*pi, 801); rMoon = l_em/DU;
xMOON = -mu + rMoon*cos(theta); yMOON = rMoon*sin(theta);
plot(xMOON, yMOON, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.6, 'HandleVisibility', 'off')  % Grigio scuro e tratteggiata

grid on
axis equal
    

%% 4) N-Body PROPAGATION:


% Upper/Lower Bounds for et:
tFormat = 'YYYY-MON-DD-HR:MN:SC.####::TDB';
tLstr = 'September 28 00:00:00.000 TDB 2024';
etL = cspice_str2et(tLstr);

% DEFINE CONSTANTS:
ms = 3.28900541e05;
rho = 3.88811143e02;
om_s = -9.25195985e-01;

om_m = 2.661861350e-06;
moonRev = 2*pi/om_m;
etU = etL + moonRev;

% Initial state and time:
xxi = vars_opt(1:4)';
ti = vars_opt(5);
tf = vars_opt(6);
TOF = (tf-ti)*TU*86400;
% Rotate initial state to ECI, append zi = 0, vzi = 0, dimensionalize:

xxi_ECI  = rot2ECI(ti, xxi, mu)';
xxi_ECI = [xxi_ECI(1:2); 0; xxi_ECI(3:4); 0];
xxi_ECI = [xxi_ECI(1:3) .* DU; xxi_ECI(4:6) .* VU];


% Define target angle
thetaTarget = mod(om_s * ti, 2*pi); % target angle [0, 2*pi]


etVec = linspace(etL, etU, 10000);

figure()
plot(etVec, thetaFinder(etVec, thetaTarget, 'J2000', 'EMB'))
hold on
plot(etVec, zeros(1, length(etVec)), 'k')

thetaFinderFun = @(et) thetaFinder(et, thetaTarget, 'J2000', 'EMB');

% Solve the zero-finding problem:
etGuess = 781169000;
initialEpoch = fzero(thetaFinderFun, [etGuess etU]); 

figure()
plot(etVec, thetaFinder(etVec, thetaTarget, 'J2000', 'EMB'))
hold on
plot(etVec, zeros(1, length(etVec)), 'k')
plot(initialEpoch, thetaFinder(initialEpoch, thetaTarget, 'J2000', 'EMB'), 'o')

% N-body propagation:
finalEpoch = initialEpoch + TOF;
xxi = xxi_ECI;
labels = {'Sun'; 'Mercury'; 'Venus'; 'Earth'; 'Moon'; 'Mars Barycenter'; 'Jupiter Barycenter';
          'Saturn Barycenter'; 'Uranus Barycenter'; 'Neptune Barycenter'; 'Pluto Barycenter'};
% Save integration frame string for SPICE:
center = 'Earth';
frame = 'J2000';
[bodies] = nBody_init(labels);
[tt, xx] = nBodyPropagator(xxi, initialEpoch, finalEpoch, bodies, frame, center);

figure
plot(xx(:, 1), xx(:, 2))
hold on
% Plot Moon Orbit:
theta = linspace(0, 2*pi, 801); rMoon = l_em/DU;
xMOON = -mu + rMoon*cos(theta); yMOON = rMoon*sin(theta);
xMOON = xMOON * DU;
yMOON = yMOON * DU;
plot(xMOON, yMOON, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.6, 'HandleVisibility', 'off')  % Grigio scuro e tratteggiata

grid on
axis equal

toc

%% Functions

function [tt, xx] = propagator(xx0, ti, tf, mu)

optset = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[tt, xx] = ode78(@(t,xx) PBRFBP(t, mu, xx), [ti tf], xx0, optset);

end

function [dxdt] = PBRFBP(t, mu, xx)
x = xx(1);
y = xx(2);
vx = xx(3);
vy = xx(4);

% DEFINE CONSTANTS:
ms = 3.28900541e05;
rho = 3.88811143e02;
om_s = -9.25195985e-01;

dOMdx = x - (ms*cos(om_s*t))/rho^2 - (mu*(mu + x - 1))/((mu + x - 1)^2 + y^2)^(3/2) + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(3/2)) - (ms*(2*x - 2*rho*cos(om_s*t)))/(2*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2));
dOMdy = y - (ms*sin(om_s*t))/rho^2 - (mu*y)/((mu + x - 1)^2 + y^2)^(3/2) - (ms*(2*y - 2*rho*sin(om_s*t)))/(2*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2)) + (y*(mu - 1))/((mu + x)^2 + y^2)^(3/2);

dxdt = zeros(4, 1);
dxdt(1:2) = xx(3:4);
dxdt(3) = 2*vy + dOMdx;
dxdt(4) = -2*vx + dOMdy;
end


function XX = rot2ECI(tt, xx, mu)

% CONSTANTS
% om_em = 2.66186135e-06;
% l_em = 3.84405e08; % [m]
% VU = 1.02454018e03; % [m/s]
% 
% om_scaled = om_em*l_em/VU;
% 
% XX(:, 1) = (xx(:, 1) + mu) .* cos(om_scaled .* tt) - xx(:, 2) .* sin(om_scaled .* tt);
% XX(:, 2) = (xx(:, 1) + mu) .* sin(om_scaled .* tt) + xx(:, 2) .* cos(om_scaled .* tt);
% XX(:, 3) = 0;
% XX(:, 4) = 0;

XX(:, 1) = (xx(:, 1) + mu) .* cos(tt) - xx(:, 2) .* sin(tt);
XX(:, 2) = (xx(:, 1) + mu) .* sin(tt) + xx(:, 2) .* cos(tt);
XX(:, 3) = (xx(:, 3) -xx(:, 2)) .* cos(tt) - (xx(:, 4) + xx(:,1) + mu).* sin(tt);
XX(:, 4) = (xx(:, 3) -xx(:, 2)) .* sin(tt) + (xx(:, 4) + xx(:,1) + mu).* cos(tt);


end

function [DV, gradDV] = costFunction(vars, constants)

xx0 = vars(1:4);
ti = vars(5);
tf = vars(6);
% EXTRACT CONSTANTS:
r0 = constants.r0;
rf = constants.rf;
mu = constants.mu;
v0 = sqrt((1 - mu)/r0);
vf = sqrt(mu/rf);

% Extract input variables:
xi = xx0(1);
yi = xx0(2);
vxi = xx0(3);
vyi = xx0(4);

DV1 = sqrt((vxi - yi)^2 + (vyi + xi + mu)^2) - v0;


% Propagate orbit to tf:
% [~, xx] = propagator(xx0, ti, tf, mu);
[xFinal, ~,  ~, ~, PHI] = propagatorSTM(xx0, ti, tf, mu);
% xFinal = xx(end, :);
xf = xFinal(1);
yf = xFinal(2);
vxf = xFinal(3);
vyf = xFinal(4);

DV2 = sqrt((vxf - yf)^2 + (vyf + xf + mu - 1)^2) - vf;

% DV = abs(DV1) + abs(DV2);
DV = DV1 + DV2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GRADIENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P1 = 1/sqrt((vxi - yi)^2 + (vyi + xi + mu)^2) .* [vyi + xi + mu; yi - vxi; vxi - yi; vyi + xi + mu];
PN = 1/sqrt((vxf - yf)^2 + (vyf + xf + mu - 1)^2) .* [vyf + xf + mu - 1; yf - vxf; vxf - yf; vyf + xf + mu - 1];
 
f_ti = PBRFBP(ti, mu, xx0);
f_tf = PBRFBP(tf, mu, xFinal);


DV_ti = - (PHI' * PN)' * f_ti;
DV_tf = PN' * f_tf;
 
gradDV = [P1 + PHI' * PN;
          DV_ti;
          DV_tf];
end

function [c, ceq, dC, dCeq] = nonlcon(vars, constants)

xx0 = vars(1:4);
ti = vars(5);
tf = vars(6);    
% Estrai le variabili
    x0 = xx0(1); y0 = xx0(2);
    vx0 = xx0(3); vy0 = xx0(4);

    % Parametri orbitali (es. per un problema CR3BP)
    mu = constants.mu;  % Definisci il parametro mu (massa ridotta)
    r0 = constants.r0;  % Definisci il raggio orbitale della parking orbit iniziale
    rf = constants.rf;  % Definisci il raggio orbitale della parking orbit finale

    % Integrazione dal tempo ti a tf per ottenere lo stato finale
    % [~, xx] = propagator(xx0, ti, tf, mu);
    [xFinal, ~,  ~, ~, PHI] = propagatorSTM(xx0, ti, tf, mu);
    % xFinal = xx(end, :);
    xf = xFinal(1);
    yf = xFinal(2);
    vxf = xFinal(3);
    vyf = xFinal(4);

    % 1. Constraint di uguaglianza: (x0 + mu)^2 + y0^2 - r0^2 = 0
    ceq1 = (x0 + mu)^2 + y0^2 - r0^2;

    % 2. Constraint di uguaglianza: (x0 + mu)*(vx0 - y0) + y0*(vy0 + x0 + mu) = 0
    ceq2 = (x0 + mu)*(vx0 - y0) + y0*(vy0 + x0 + mu);

    % 3. Constraint di uguaglianza: (xf + mu - 1)^2 + yf^2 - rf^2 = 0
    ceq3 = (xf + mu - 1)^2 + yf^2 - rf^2;

    % 4. Constraint di uguaglianza: (xf + mu - 1)*(vxf - yf) + yf*(vyf + xf + mu - 1) = 0
    ceq4 = (xf + mu - 1)*(vxf - yf) + yf*(vyf + xf + mu - 1);

    % Non hai vincoli di disuguaglianza
    c = [];

    % Combina tutti i vincoli di uguaglianza
    ceq = [ceq1; ceq2; ceq3; ceq4];

    %%%%%%%%%%%%%%%%%%%% GRADIENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [dC, dCeq] = GradientConstraints(vars, xFinal, PHI, constants);

end

function [dC, dCeq] = GradientConstraints(vars, xFinal, PHI, constants)
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
f_ti = PBRFBP(ti, mu, xxi);
C34_ti = -C34_xxf * PHI * f_ti;

f_tf = PBRFBP(tf, mu, xFinal);
C34_tf = C34_xxf * f_tf;

dCeq = [C12_xxi, C12_ti, C12_tf;
        C34_xxi, C34_ti, C34_tf];

% dC = [0 0 0 0 1 -1]';
dC = [];
dCeq = dCeq';

end


function [vars0] = initialGuessMS(xx0, ti, tf, N, constants)

mu = constants.mu;
% Compute initial guess with internal points:

tvec = zeros(N,1); % time vector
for i = 1 : N
    tvec(i) = ti + (i-1)/(N-1) * (tf - ti);
end

xvecGuesses = zeros(4, N); % Matrix of initial guesses
xvecGuesses(:, 1) = xx0;
for i = 2 : N
    [~, xx] = propagator(xvecGuesses(:, i-1), tvec(i-1), tvec(i), mu);
    xvecGuesses(:, i) = xx(end, :)';
end
vars0 = zeros(4*N + 2, 1); % vector of initial guesses
vars0(1 : 4*N) = reshape(xvecGuesses, 4*N, 1);
vars0(4*N + 1) = ti;
vars0(4*N + 2) = tf;
end

function [DV, gradDV] = costFunctionMS(vars, constants, N)

mu = constants.mu;
r0 = constants.r0;
rf = constants.rf;
v0 = sqrt((1 - mu)/r0);
vf = sqrt(mu/rf);

varsRSP = reshape(vars(1:4*N), 4, N).'; % da vettore delle variabili a matrice 4*N in cui ogni colonna è lo stato al nodo j-esimo

% STEP 2: propagazione dalle variabili e ASSEMBLAMENTO constraint:
xx1 = varsRSP(1, :);
xxN = varsRSP(N, :);

% Extract components from initial and final state vectors:
xi = xx1(1); yi = xx1(2); vxi = xx1(3); vyi = xx1(4);
xf = xxN(1); yf = xxN(2); vxf = xxN(3); vyf = xxN(4);

% Compute deltaV:
DV1 = sqrt((vxi - yi)^2 + (vyi + xi + mu)^2) - v0;
DV2 = sqrt((vxf - yf)^2 + (vyf + xf + mu - 1)^2) - vf;
DV = DV1 + DV2;

P1 = 1/sqrt((vxi - yi)^2 + (vyi + xi + mu)^2) .* [vyi + xi + mu; yi - vxi; vxi - yi; vyi + xi + mu];
PN = 1/sqrt((vxf - yf)^2 + (vyf + xf + mu - 1)^2) .* [vyf + xf + mu - 1; yf - vxf; vxf - yf; vyf + xf + mu - 1];
gradDV = zeros(18, 1);
gradDV(1 : 4, 1) = P1;
gradDV(3*N + 1 : 4*N, 1) = PN;
end

function [c, ceq, dC, dCeq] = nonlconMS(vars, constants, N)

DU = 3.84405000e05;
mu = constants.mu;
Re = constants.R_body1 / DU;
Rm = constants.R_body2 / DU;
r0 = constants.r0;
rf = constants.rf;

varsRSP = reshape(vars(1:4*N), 4, N).';
t1 = vars(4*N+1);
tN = vars(4*N+2);

tvec = zeros(N, 1);
tvec(1) = t1;
for i = 1 : N
    tvec(i) = t1 + (tN - t1) * (i-1)/(N-1);
end
% tvec(N) = tN;

% propagatedStates = zeros(N-1, 4);

% Equality Constraints:
ceq = zeros(4*N, 1);
% def = zeros(N-1, 4);
% Initial position & Velocity:
ceq(end-3) = (varsRSP(1,1) + mu)^2 + varsRSP(1,2)^2 - r0^2;
ceq(end-2) = (varsRSP(1,1) + mu) * (varsRSP(1,3) - varsRSP(1,2)) + varsRSP(1,2) * (varsRSP(1,4) + varsRSP(1,1) + mu);
% Final position & velocity
ceq(end-1) = (varsRSP(end,1) + mu - 1)^2 + varsRSP(end,2)^2 - rf^2;
ceq(end) = (varsRSP(end,1) + mu - 1) * (varsRSP(end,3) - varsRSP(end,2)) + varsRSP(end,2) * (varsRSP(end,4) + varsRSP(end,1) + mu - 1);

 
for j = 1 : N-1
    % Extract current state and next state from variables vector:
    xxj = vars(4*j-3:4*j);
    xxjj = vars(4*j+1 : 4*j+4);
    % [~, xx] = propagator(varsRSP(j, :), tvec(j), tvec(j+1), mu);
    [~, ~, xx, ~, ~] = propagatorSTM(xxj, tvec(j), tvec(j+1), mu);
    % propagatedStates(j, :) = xx(end, :);
    ceq(4*j-3 : 4*j) = xx(1:4, end) - xxjj;
end

% Inequality Constraints:
c = zeros(2*N + 1, 1);
% c(1) = (Re)^2 - (varsRSP(1, 1) + mu)^2 - varsRSP(1, 2)^2;
% c(2) = (Rm)^2 -(varsRSP(1, 1) + mu - 1)^2 - varsRSP(1, 2)^2;
for j = 1 : N
    % xxj = vars(4*j-3:4*j);
    % [~, ~, xx, ~, ~] = propagatorSTM(xxj, tvec(j), tvec(j+1), mu);
    % c(2*j + 1) = (Re)^2 - (xx(1, end) + mu)^2 - xx(2, end)^2;
    % c(2*j + 2) = (Rm)^2 -(xx(1, end) + mu - 1)^2 - xx(1, end)^2;
    c(2*j - 1) = (Re)^2 - (varsRSP(j, 1) + mu)^2 - varsRSP(j, 2)^2;
    c(2*j) = (Rm)^2 -(varsRSP(j, 1) + mu - 1)^2 - varsRSP(j, 2)^2;
end
c(2*N + 1) = t1 - tN;

% ceq(11) = (ceq(11) / 3e-9) - 1;  % Normalizza rispetto a 3*10^-9

% COMPUTE GRADIENTS OF CONSTRAINTS:
% [dC, dCeq] = JacobianCons(vars, constants, N, tvec);
dCeq = zeros(3*N+4, 4*N+2);
dC = zeros(2*N + 1, 4*N+2);

for j = 1 : N - 1
    [xxf, ~, ~, ~, PHIj] = propagatorSTM(varsRSP(j, :)', tvec(j), tvec(j+1), mu);
    dCeq(4*j-3:4*j, 4*j-3:4*j) = PHIj;
    dCeq(4*j-3:4*j, 4*j+1:4*j+4) = - eye(4);
    fj = PBRFBP(tvec(j), mu, varsRSP(j, :)');
    fjProp = PBRFBP(tvec(j+1), mu, xxf);
    Q1j = - (N-j)/(N-1).* PHIj * fj + (N-j-1)/(N-1) .* fjProp;
    QNj = -(j-1)/(N-1) .* PHIj * fj + j/(N-1) .* fjProp;
    dCeq(4*j-3:4*j, 4*N+1) = Q1j;
    dCeq(4*j-3:4*j, 4*N+2) = QNj;
end

xx1 = varsRSP(1, :);
xxN = varsRSP(N, :);

% Extract components from initial and final state vectors:
xi = xx1(1); yi = xx1(2); vxi = xx1(3); vyi = xx1(4);
xf = xxN(1); yf = xxN(2); vxf = xxN(3); vyf = xxN(4);


R1 = [2*(xi + mu), 2*yi, 0, 0;
      vxi, vyi, xi + mu, yi];
RN = [2*(xf + mu - 1), 2*yf, 0, 0;
      vxf, vyf, xf + mu - 1, yf];

dCeq(3*N+1 : 3*N+2, 1:4) = R1;
dCeq(3*N+3 : 3*N+4, 3*N+1:4*N) = RN;

for j = 1 : N
    xj = varsRSP(j, 1);
    yj = varsRSP(j, 2);
    Sj = [-2*(xj + mu), -2*yj, 0, 0;
          -2*(xj + mu - 1), -2*yj, 0, 0];
    dC(2*j-1:2*j, 4*j-3:4*j) = Sj;
end

St = [1 -1];
dC(end, end-1:end) = St;

dCeq = dCeq';
dC = dC';

end

function [constants] = cspiceConstants(body1, body2)
    GM1 = cspice_bodvrd(body1, 'GM', 1);
    GM2 = cspice_bodvrd(body2, 'GM', 1);
    constants.mu = GM2/(GM1 + GM2);

    R_body1 = cspice_bodvrd(body1, 'RADII', 3);
    constants.R_body1 = R_body1(1);
    R_body2 = cspice_bodvrd(body2, 'RADII', 3);
    constants.R_body2 = R_body2(1);
end


function [] = plotMultipleShooting(varsOpt, dV, constants, N)
% Funzione per il plot della traiettoria con metodo multiple shooting
% dV --> da mettere nel titolo del plot
% title('Multiple shooting trajectory, N = 'N', total cost = 'dV')

% CONSTANTS:
mu = constants.mu;
r0 = constants.r0;
rf = constants.rf;

% STEP 1: reshape variables and compute vector of times:
varsRSP = reshape(varsOpt(1:4*N), 4, N); % da vettore delle variabili a matrice 4*N
t1 = varsOpt(4*N + 1);
tN = varsOpt(4*N + 2);

tvec = linspace(t1, tN, N)'; % vettore dei tempi

figure()
grid on
hold on
axis equal

% Plot circumferences of r0 and rf centered in Earth and Moon
theta = linspace(0, 2*pi, 801);  
xEarth = -mu + r0 * cos(theta);  
yEarth = r0 * sin(theta); 
xMoon = (1 - mu) + rf * cos(theta);  
yMoon = rf * sin(theta); 

% Disegna orbite iniziali e finali con leggenda
plot(xEarth, yEarth, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8, 'DisplayName', 'Initial orbit')  
plot(xMoon, yMoon, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.6, 'DisplayName', 'Final orbit')  

% Plot dei nodi e traiettoria tra i nodi
for j = 1 : N - 1
    [~, xx] = propagator(varsRSP(:, j), tvec(j), tvec(j + 1), mu);
    plot(varsRSP(1, j), varsRSP(2, j), 'ro', 'MarkerFaceColor', 'r', 'HandleVisibility', 'off') % Nodo interno
    plot(xx(:, 1), xx(:, 2), 'b', 'HandleVisibility','off') % Traiettoria tra nodo j e j+1 
end

% Ultimo nodo
plot(varsRSP(1, end), varsRSP(2, end), 'ro', 'MarkerFaceColor', 'r', 'HandleVisibility', 'off')

% Disegna Terra e Luna con colori e legenda
plot(-mu, 0, 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Earth') 
plot(1-mu, 0, 'o', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], 'DisplayName', 'Moon') 

% Aggiungi titolo e legenda
title(['Multiple shooting trajectory, N = ', num2str(N), ', $\Delta v_{TOT}$ = ', num2str(dV), ' [km/s]'], 'Interpreter', 'latex')
legend('show', 'Location', 'best')

hold off
end

%%%%%%%% CALCOLI CON STM:

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
 

% function [dC, dCeq] = JacobianCons(vars, constants, N, tvec)
% 
% mu = constants.mu;
% % GRADIENTI:
% dCeq = zeros(3*N+4, 4*N+2);
% dC = zeros(2*N + 1, 4*N+2);
% 
% varsRSP = reshape(vars(1:4*N), 4, N).';
% 
% for j = 1 : N - 1
%     [xxf, ~, ~, ~, PHIj] = propagatorSTM(varsRSP(j, :)', tvec(j), tvec(j+1), mu);
%     dCeq(4*j-3:4*j, 4*j-3:4*j) = PHIj;
%     dCeq(4*j-3:4*j, 4*j+1:4*j+4) = - eye(4);
%     fj = PBRFBP(tvec(j), mu, varsRSP(j, :)');
%     fjProp = PBRFBP(tvec(j+1), mu, xxf);
%     Q1j = - (N-j)/(N-1).* PHIj * fj + (N-j-1)/(N-1) .* fjProp;
%     QNj = -(j-1)/(N-1) .* PHIj * fj + j/(N-1) .* fjProp;
%     dCeq(4*j-3:4*j, 4*N+1) = Q1j;
%     dCeq(4*j-3:4*j, 4*N+2) = QNj;
% end
% 
% % Extract components from initial and final state vectors:
% xi = varsRSP(1, 1); yi = varsRSP(1, 2); vxi = varsRSP(1, 3); vyi = varsRSP(1, 4);
% xf = varsRSP(N, 1); yf = varsRSP(N, 2); vxf = varsRSP(N, 3); vyf = varsRSP(N, 4);
% 
% R1 = [2*(xi + mu), 2*yi, 0, 0;
%       vxi, vyi, xi + mu, yi];
% RN = [2*(xf + mu - 1), 2*yf, 0, 0;
%       vxf, vyf, xf + mu - 1, yf];
% 
% dCeq(3*N+1 : 3*N+2, 1:4) = R1;
% dCeq(3*N+3 : 3*N+4, 3*N+1:4*N) = RN;
% 
% for j = 1 : N
%     xj = varsRSP(j, 1);
%     yj = varsRSP(j, 2);
%     Sj = [-2*(xj + mu), -2*yj, 0, 0;
%           -2*(xj + mu - 1), -2*yj, 0, 0];
%     dC(2*j-1:2*j, 4*j-3:4*j) = Sj;
% end
% 
% St = [1 -1];
% dC(end, end-1:end) = St;
% 
% dCeq = dCeq';
% dC = dC';
% end


% %%%%%%%%%%%%%%%%%%%%%%%%% PT. 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
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
% Perform Integration:
options = odeset('reltol', 1e-12, 'abstol', [1e-6*ones(1,3),1e-9*ones(1,3)]);
[tt, xx] = ode113(@(t,x) nBody_RHS(t, x, bodies, frame, center), [et0 etf], xx0, options);
end

function thetaFinderFun = thetaFinder(etVec, thetaTarget, frame, center)

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