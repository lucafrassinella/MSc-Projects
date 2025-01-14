% Spacecraft Guidance and Navigation - A.Y. 2024/25
% Assignment #1
% Exercise #3
% Author: Frassinella Luca - 10795356

%% -------------------- EX. 3 - CONTINOUS GUIDANCE --------------------- %%

clearvars; close all; clc;
format long
plotSettings;
rng default

% Problem Data:
hi = 800; % [km]
hf = 1000; % [km]
i = deg2rad(0.75); % [rad]
Re = 6378.1366; % [km]
mu = 398600.435; % [km^3/s^2]
rho0 = 750 + Re; % [km]
k1 = 1e-05;
k2 = 1e-04;
m0 = 1000; % [kg]
Tmax = 3e-03; % [kN]
Isp = 3120; % [s]
g0 = 9.81e-3; % [km/s^2]

% Adimensionalization Units:
DU = 7178.1366; % [km]
MU = m0;
TU = sqrt(DU^3/mu);



%% Ex 01: Plot Debris Spatial Density and Compute Initial and Target State:

% Initial and Target States:
initialState = [Re + hi; 0; 0; 0; sqrt(mu/(Re+hi)); 0];
finalState = [Re + hf; 0; 0; 0; sqrt(mu/(Re+hf))*cos(i); sqrt(mu/(Re+hf))*sin(i)];

% Store results:
results.ex1.xx0 = initialState;
results.ex2.xxf = finalState;

% Debris Spatial Density:
q = @(rho) k1./(k2 + ((rho-rho0)/DU).^2);

% Define orbital height and radius vectors:
h_vec = linspace(hi - 100, hf + 100, 1000);
rho_vec = Re + h_vec;

% Plot the debris density:
figure()
plot(rho_vec-Re, q(rho_vec));
hold on
plot(ones(100, 1)*(hi), linspace(0, 0.15, 100), '--k')
plot(ones(100, 1)*(hf), linspace(0, 0.15, 100), '--k')
legend('Debris Spatial Density', 'Initial Orbit Altitude', 'Final Orbit Altitude')
xlabel('Altitude [km]', 'FontSize', 30)
ylabel('Debris Spatial Density (q)', 'FontSize', 30)
grid on
ylim([0, 0.12])
title('Debris Spatial Density over Altitude range of interest')

%% PART 02: Adimensionalization of the problem

% Compute adimensionalized variables and store them in data struct:
data.TU = TU;
data.MU = MU;
data.DU = DU;
data.k1 = k1;
data.k2 = k2;
data.Tmax = Tmax * ((TU^2)/(MU*DU));
data.Isp = Isp/TU;
data.m0 = m0/MU;
data.g0 = g0*(TU^2/DU);
data.mu = mu*(TU^2)/(DU^3);
data.t0 = 0;
data.rho0 = rho0/DU;
% Set ODE options
settings.odeOpt = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% Adimensionalized initial and target state:
x0 = [initialState; m0];
x0 = [x0(1:3)/DU; x0(4:6)*(TU/DU); x0(7)/MU];
xf = [finalState(1:3)/DU; finalState(4:6)*(TU/DU)];

% Store results for point 2 (adimensionalized quantities, adimensionalized 
% initial and target states)
results.ex2.adim_qties = data;
results.ex2.xx0 = x0;
results.ex2.xxf;

%% Ex. 4: Continuous Guidance Solution (T = 3.000 N)
clc; clearvars -except data x0 xf results settings;

% Set Maximum Number Of Iterations to 5:
N_iter = 5;
[sol, err, gradErr] = computeOrbit(N_iter, data, x0, xf, settings);

% Save solutions:
results.ex4.ll0 = sol(1:7);
results.ex4.tf = sol(end);
results.ex4.errPos = norm(err(1:3).*data.DU);
results.ex4.errVel = norm(err(4:6).*data.DU/data.TU);

% Propagation from initial guess:
ll0 = results.ex4.ll0;
tf = results.ex4.tf;
state0 = [x0; ll0];
[tt, xx] = propagate(state0, data.t0, tf, data, settings);

% Compute adimensionalized final mass:
results.ex4.mf = xx(end, 7)*data.MU;

% Check Properties of the Hamiltonian in time-independent problems:
% Initialize hamiltonian vector:
H = zeros(length(tt), 1);
% Compute hamiltonian at each time:
for i = 1 : length(tt)
    % Extract quantities from state and costate vector:
    rr_i = xx(i, 1:3)';
    vv_i = xx(i, 4:6)';
    m_i = xx(i, 7);
    ll_i = xx(i, 8:14)';
    r_i = norm(rr_i);
    lv_i = norm(ll_i(4:6));
    
    % Debris spatial density:
    q_i = data.k1/(data.k2+(r_i - data.rho0)^2);
    
    % Dynamics:
    fi = [vv_i;
         -data.mu/r_i^3*rr_i - data.Tmax/m_i * ll_i(4:6)/lv_i;
         -data.Tmax/(data.Isp*data.g0)];
    
    % Compute Hamiltonian:
    H(i) = q_i + dot(ll_i, fi);
end

% Plot time evolution of the Hamiltonian:
figure()
plot(tt, H, 'k')
xlabel('Time [-]')
ylabel('Hamiltonian [-]')
xlim([tt(1) tt(end)])
grid on
title('Hamiltonian Time Evolution')
subtitle('$T_{max}$=3.000 N')

% Primer vector components evolution in the NTW ref. frame:
alpha = - xx(:, 11:13)' ./ vecnorm(xx(:, 11:13)');
% Rotate to NTW frame:
alphaNTW = zeros(3, length(alpha));
for i = 1 : length(alpha)
    xxECI = xx(i, 1:6);
    alphaECI = alpha(:, i);
    alphaNTW(:, i) = rot2NTW(xxECI, alphaECI);
end
normAlpha = vecnorm(alphaNTW);

figure()
plot(tt, alphaNTW(1, :), 'DisplayName', '$\alpha_N$')
hold on
plot(tt, alphaNTW(2, :), 'DisplayName', '$\alpha_T$')
plot(tt, alphaNTW(3, :), 'DisplayName', '$\alpha_W$')
plot(tt, normAlpha, 'DisplayName', 'norm($\alpha$)')
xlabel('Time [-]', 'Interpreter', 'latex')
ylabel('$\alpha$', 'Interpreter', 'latex')
legend('Interpreter', 'latex')
xlim([tt(1) tt(end)])
grid on
title('Time Evolution of $\alpha$ components in NTW frame')

% Plot Trajectory:
figure()
plot3(xx(:, 1), xx(:, 2), xx(:, 3), 'k', 'LineWidth', 1.5) 
hold on
plot3(x0(1), x0(2), x0(3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'k', 'DisplayName', 'Initial Point') 
plot3(xf(1), xf(2), xf(3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'k', 'DisplayName', 'Target Point') 
plot3(0, 0, 0, 'o', 'MarkerSize', 25, 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'k', 'DisplayName', 'Earth') 
xlabel('x [-]')
ylabel('y [-]')
zlabel('z [-]')
legend('Transfer Trajectory', 'Initial Point', 'Target Point', 'Earth')
grid on
view([-45 15])
title('Orbit Raising Trajectory (@Earth ECI)')
subtitle('$T_{max}$=3.000 N')


%% Ex.5 - Lower Thrust Level solution through numerical continuation:

% initialize Initial Guess:
TmaxVec = linspace(3e-03, 2.86e-03, 10)'.*((data.TU^2)/(data.MU*data.DU));
soll = zeros(8, length(TmaxVec));
soll(:, 1) = sol;
guess = sol;

% Numerical Continuation Loop:
for i = 2 : length(TmaxVec)
    data.Tmax = TmaxVec(i); 
    options = optimoptions('fsolve', 'MaxFunctionEvaluations', 10000, 'MaxIterations', 10000, 'Display','iter-detailed',...
   'SpecifyObjectiveGradient', true, 'FunctionTolerance', 1e-06, 'Algorithm', 'trust-region-dogleg');
    optimizationFun = @(guess) optimization(x0, guess, data, xf, settings);
    [~, gradErr] = checkGradients(optimizationFun, guess);
    [sol, err, exitflag] = fsolve(optimizationFun, guess, options);
    guess = sol;
    soll(:, i) = sol;
end

% Save final results:
results.ex5.ll0 = sol(1:7);
results.ex5.tf = sol(end);
results.ex5.errPos = norm(err(1:3).*data.DU);
results.ex5.errVel = norm(err(4:6).*data.DU/data.TU);

% Propagation from initial guess:
ll0 = sol(1:7);
tf = sol(end);
state0 = [x0; ll0];
[tt, xx] = propagate(state0, data.t0, tf, data, settings);
results.ex5.mf = xx(end, 7)*data.MU;

% Check Properties of the Hamiltonian in time-independent problems:
H = zeros(length(tt), 1);
for i = 1 : length(tt)
    % Hi:
    rr_i = xx(i, 1:3)';
    vv_i = xx(i, 4:6)';
    m_i = xx(i, 7);
    ll_i = xx(i, 8:14)';
    r_i = norm(rr_i);
    lv_i = norm(ll_i(4:6));
    
    q_i = data.k1/(data.k2+(r_i - data.rho0)^2);
    
    % Dynamics:
    fi = [vv_i;
         -data.mu/r_i^3*rr_i - data.Tmax/m_i * ll_i(4:6)/lv_i;
         -data.Tmax/(data.Isp*data.g0)];
    
    H(i) = q_i + dot(ll_i, fi);
end

%%
% Plot time evolution of the Hamiltonian:
figure()
plot(tt, H, 'k')
xlabel('Time [-]')
ylabel('Hamiltonian [-]')
xlim([tt(1) tt(end)])
grid on
title('Hamiltonian Time Evolution')
subtitle('$T_{max}$=2.860 N')

% Plot Trajectory:
figure()
plot3(xx(:, 1), xx(:, 2), xx(:, 3), 'k', 'LineWidth', 1.5) 
hold on
plot3(x0(1), x0(2), x0(3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'k', 'DisplayName', 'Initial Point') 
plot3(xf(1), xf(2), xf(3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'k', 'DisplayName', 'Target Point') 
plot3(0, 0, 0, 'o', 'MarkerSize', 25, 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'k', 'DisplayName', 'Earth') 
xlabel('x [-]')
ylabel('y [-]')
zlabel('z [-]')
legend('Transfer Trajectory', 'Initial Point', 'Target Point', 'Earth')
grid on
view([-45 15])
title('Orbit Raising Trajectory (@Earth ECI)')
subtitle('$T_{max}$=2.860 N')

%% FUNCTIONS

function [dx] = TPBVP(~, x, data)

mu = data.mu;
g0 = data.g0;
Isp = data.Isp;
Tmax = data.Tmax;
k1 = data.k1;
k2 = data.k2;
rho0 = data.rho0;

rr = x(1:3);
vv = x(4:6);
m = x(7);
ll_r = x(8:10);
ll_v = x(11:13);

r = norm(rr);
l_v = norm(ll_v);

dq = [(2*k1*rr(1)*(rho0 - (rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2)))/((k2 + (rho0 - (rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2))^2)^2*(rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2));
      (2*k1*rr(2)*(rho0 - (rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2)))/((k2 + (rho0 - (rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2))^2)^2*(rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2));
      (2*k1*rr(3)*(rho0 - (rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2)))/((k2 + (rho0 - (rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2))^2)^2*(rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2))];

dx = zeros(210, 1);
dx(1:3) = vv;
dx(4:6) = -(mu/r^3)*rr - Tmax/m * ll_v/l_v;
dx(7) = -Tmax/(Isp*g0);
dx(8:10) = - dq - 3*mu/(r^5) * dot(rr, ll_v)*rr + mu/r^3*ll_v;
dx(11:13) = - ll_r;
dx(14) = -l_v*Tmax/m^2;

STM = x(15:end);
% Compute derivative of the STM:
STM = reshape(STM, 14, 14);
F = jacobianDyn(x, data);
STMdot = F*STM;
dx(15:end) = STMdot(:);

end

function [dx] = dynamics(~, x, data)

mu = data.mu;
g0 = data.g0;
Isp = data.Isp;
Tmax = data.Tmax;
k1 = data.k1;
k2 = data.k2;
rho0 = data.rho0;

rr = x(1:3);
vv = x(4:6);
m = x(7);
ll_r = x(8:10);
ll_v = x(11:13);

r = norm(rr);
l_v = norm(ll_v);

dq = [(2*k1*rr(1)*(rho0 - (rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2)))/((k2 + (rho0 - (rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2))^2)^2*(rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2));
      (2*k1*rr(2)*(rho0 - (rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2)))/((k2 + (rho0 - (rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2))^2)^2*(rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2));
      (2*k1*rr(3)*(rho0 - (rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2)))/((k2 + (rho0 - (rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2))^2)^2*(rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2))];

dx = zeros(210, 1);
dx(1:3) = vv;
dx(4:6) = -(mu/r^3)*rr - Tmax/m * ll_v/l_v;
dx(7) = -Tmax/(Isp*g0);
dx(8:10) = - dq - 3*mu/(r^5) * dot(rr, ll_v)*rr + mu/r^3*ll_v;
dx(11:13) = - ll_r;
dx(14) = -l_v*Tmax/m^2;

end

function [tt, xx, STM] = propagate(state0, t0, tf, data, settings)

% Append initial STM to initial state:
STM0 = eye(14);
xv0 = [state0; STM0(:)];

fun = @(t, x) TPBVP(t, x, data);
[tt, xx] = ode113(fun, [t0 tf], xv0, settings.odeOpt);

STMf = xx(end, 15:end);
STM = reshape(STMf, 14, 14);

end

function [residuals, grad] = optimization(x0, guess, data, xf, settings)

mu = data.mu;
g0 = data.g0;
Isp = data.Isp;
Tmax = data.Tmax;
k1 = data.k1;
k2 = data.k2;
rho0 = data.rho0;
t0 = data.t0;

% Propagate s/c dynamics:
tf = guess(end);
ll0 = guess(1:7);
state0 = [x0; ll0];
[tt, xx, STM] = propagate(state0, t0, tf, data, settings);


rr_f = xx(end, 1:3)';
vv_f = xx(end, 4:6)';
m_f = xx(end, 7);
ll_f = xx(end, 8:14)';
r_f = norm(rr_f);
lv_f = norm(ll_f(4:6));
tfinal = tt(end);


q_f = k1/(k2+(r_f - rho0)^2);

% Dynamics:
f = [vv_f;
     -mu/r_f^3*rr_f - Tmax/m_f * ll_f(4:6)/lv_f;
     -Tmax/(Isp*g0)];

H = q_f + dot(ll_f, f);

% Final function to evaluate:
residuals = [rr_f - xf(1:3);
             vv_f - xf(4:6);
             ll_f(7)
             H];


% Compute gradients of objective function:
statefinal = xx(end, 1:14)';
grad = objectiveJacobian(tfinal, statefinal, STM, data);
end

function [sol, err, gradErr] = computeOrbit(N_iter, data, x0, xf, settings)

exitflag = 0;
iter = 0;
tf0 = 20*pi; % initial time guess

while exitflag <= 0 && iter < N_iter
    lambdaRange = [-250 250];
    lambda0 = zeros(7,1);
    lambda0(1:6) = (lambdaRange(2)-lambdaRange(1)).*rand(6,1) +lambdaRange(1);
    lambda0(7) = lambdaRange(2).*rand(1,1);
    guess = [lambda0; tf0];
    
    options = optimoptions('fsolve', 'MaxFunctionEvaluations', 10000, 'MaxIterations', 10000, 'Display','iter-detailed',...
   'SpecifyObjectiveGradient', true, 'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10,  'Algorithm','levenberg-marquardt');
    optimizationFun = @(guess) optimization(x0, guess, data, xf, settings);
    [~, gradErr] = checkGradients(optimizationFun, guess);
    [sol, err, exitflag] = fsolve(optimizationFun, guess, options);
    % Discard solution if the final time is not between 20*pi and 22*pi 
    % (i.e., the transfer requires less than 10 or more than 11 revolutions)
    if sol(8) < 20*pi || sol(8) > 22*pi
        exitflag = 0;
    end
    iter = iter + 1;
end

if iter == N_iter && ~exitflag
    error('Alghoritm did not converge')
end

end


function alphaNTW = rot2NTW(xxECI, alphaECI)

rrECI = xxECI(1:3);
vvECI = xxECI(4:6);

T = vvECI/norm(vvECI);
N = cross(rrECI, vvECI)/norm(cross(rrECI, vvECI));
W = cross(N, T);

RR = [N; T; W];

alphaNTW = RR * alphaECI;

end

function Jacobian = jacobianDyn(state, data)

% function to evaluate the jacobian of the Dynamics to compute the STM
% through variational approach

Tmax = data.Tmax;
mu = data.mu;
k1 = data.k1;
k2 = data.k2;
rho0 = data.rho0;

% Extract state variables:
x = state(1);
y = state(2);
z = state(3);
m = state(7);
lvx = state(11);
lvy = state(12);
lvz = state(13);



% Jacobian of the dynamics:
F = zeros(14, 14);
F(1, 4) = 1;
F(2, 5) = 1;
F(3, 6) = 1;
F(4, 1) = -(mu*(- 2*x^2 + y^2 + z^2))/(x^2 + y^2 + z^2)^(5/2);
F(4, 2) = (3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2);
F(4, 3) = (3*mu*x*z)/(x^2 + y^2 + z^2)^(5/2);
F(4, 7) = (Tmax*lvx)/(m^2*(lvx^2 + lvy^2 + lvz^2)^(1/2));
F(4, 11) = -(Tmax*(lvy^2 + lvz^2))/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2));
F(4, 12) = (Tmax*lvx*lvy)/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2));
F(4, 13) = (Tmax*lvx*lvz)/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2));
F(5, 1) = (3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2);
F(5, 2) = -(mu*(x^2 - 2*y^2 + z^2))/(x^2 + y^2 + z^2)^(5/2);
F(5, 3) = (3*mu*y*z)/(x^2 + y^2 + z^2)^(5/2);
F(5, 7) = (Tmax*lvy)/(m^2*(lvx^2 + lvy^2 + lvz^2)^(1/2));
F(5, 11) = (Tmax*lvx*lvy)/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2));
F(5, 12) = -(Tmax*(lvx^2 + lvz^2))/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2));
F(5, 13) = (Tmax*lvy*lvz)/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2));
F(6, 1) = (3*mu*x*z)/(x^2 + y^2 + z^2)^(5/2);
F(6, 2) = (3*mu*y*z)/(x^2 + y^2 + z^2)^(5/2);
F(6, 3) = -(mu*(x^2 + y^2 - 2*z^2))/(x^2 + y^2 + z^2)^(5/2);
F(6, 7) = (Tmax*lvz)/(m^2*(lvx^2 + lvy^2 + lvz^2)^(1/2));
F(6, 11) = (Tmax*lvx*lvz)/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2));
F(6, 12) = (Tmax*lvy*lvz)/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2));
F(6, 13) = -(Tmax*(lvx^2 + lvy^2))/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2));
F(8, 1) = (15*mu*x^2*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(7/2) - (3*mu*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(5/2) - (2*k1*(rho0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(1/2)) + (2*k1*x^2)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)) - (6*lvx*mu*x)/(x^2 + y^2 + z^2)^(5/2) + (2*k1*x^2*(rho0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(3/2)) - (8*k1*x^2*(rho0 - (x^2 + y^2 + z^2)^(1/2))^2)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2));
F(8, 2) = (2*k1*x*y)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)) - (3*lvx*mu*y)/(x^2 + y^2 + z^2)^(5/2) - (3*lvy*mu*x)/(x^2 + y^2 + z^2)^(5/2) + (15*mu*x*y*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(7/2) + (2*k1*x*y*(rho0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(3/2)) - (8*k1*x*y*(rho0 - (x^2 + y^2 + z^2)^(1/2))^2)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2));
F(8, 3) = (2*k1*x*z)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)) - (3*lvx*mu*z)/(x^2 + y^2 + z^2)^(5/2) - (3*lvz*mu*x)/(x^2 + y^2 + z^2)^(5/2) + (15*mu*x*z*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(7/2) + (2*k1*x*z*(rho0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(3/2)) - (8*k1*x*z*(rho0 - (x^2 + y^2 + z^2)^(1/2))^2)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2));
F(8, 11) = (mu*(- 2*x^2 + y^2 + z^2))/(x^2 + y^2 + z^2)^(5/2);
F(8, 12) = -(3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2);
F(8, 13) = -(3*mu*x*z)/(x^2 + y^2 + z^2)^(5/2);
F(9, 1) = (2*k1*x*y)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)) - (3*lvx*mu*y)/(x^2 + y^2 + z^2)^(5/2) - (3*lvy*mu*x)/(x^2 + y^2 + z^2)^(5/2) + (15*mu*x*y*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(7/2) + (2*k1*x*y*(rho0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(3/2)) - (8*k1*x*y*(rho0 - (x^2 + y^2 + z^2)^(1/2))^2)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2));
F(9, 2) = (15*mu*y^2*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(7/2) - (3*mu*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(5/2) - (2*k1*(rho0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(1/2)) + (2*k1*y^2)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)) - (6*lvy*mu*y)/(x^2 + y^2 + z^2)^(5/2) + (2*k1*y^2*(rho0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(3/2)) - (8*k1*y^2*(rho0 - (x^2 + y^2 + z^2)^(1/2))^2)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2));
F(9, 3) = (2*k1*y*z)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)) - (3*lvy*mu*z)/(x^2 + y^2 + z^2)^(5/2) - (3*lvz*mu*y)/(x^2 + y^2 + z^2)^(5/2) + (15*mu*y*z*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(7/2) + (2*k1*y*z*(rho0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(3/2)) - (8*k1*y*z*(rho0 - (x^2 + y^2 + z^2)^(1/2))^2)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2));
F(9, 11) = -(3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2);
F(9, 12) = (mu*(x^2 - 2*y^2 + z^2))/(x^2 + y^2 + z^2)^(5/2);
F(9, 13) = -(3*mu*y*z)/(x^2 + y^2 + z^2)^(5/2);
F(10, 1) = (2*k1*x*z)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)) - (3*lvx*mu*z)/(x^2 + y^2 + z^2)^(5/2) - (3*lvz*mu*x)/(x^2 + y^2 + z^2)^(5/2) + (15*mu*x*z*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(7/2) + (2*k1*x*z*(rho0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(3/2)) - (8*k1*x*z*(rho0 - (x^2 + y^2 + z^2)^(1/2))^2)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2));
F(10, 2) = (2*k1*y*z)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)) - (3*lvy*mu*z)/(x^2 + y^2 + z^2)^(5/2) - (3*lvz*mu*y)/(x^2 + y^2 + z^2)^(5/2) + (15*mu*y*z*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(7/2) + (2*k1*y*z*(rho0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(3/2)) - (8*k1*y*z*(rho0 - (x^2 + y^2 + z^2)^(1/2))^2)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2));
F(10, 3) = (15*mu*z^2*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(7/2) - (3*mu*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(5/2) - (2*k1*(rho0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(1/2)) + (2*k1*z^2)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)) - (6*lvz*mu*z)/(x^2 + y^2 + z^2)^(5/2) + (2*k1*z^2*(rho0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(3/2)) - (8*k1*z^2*(rho0 - (x^2 + y^2 + z^2)^(1/2))^2)/((k2 + (rho0 - (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2));
F(10, 11) = -(3*mu*x*z)/(x^2 + y^2 + z^2)^(5/2);
F(10, 12) = -(3*mu*y*z)/(x^2 + y^2 + z^2)^(5/2);
F(10, 13) = (mu*(x^2 + y^2 - 2*z^2))/(x^2 + y^2 + z^2)^(5/2);
F(11, 8) = -1;
F(12, 9) = -1;
F(13, 10) = -1;
F(14, 7) = (2*Tmax*(lvx^2 + lvy^2 + lvz^2)^(1/2))/m^3;
F(14, 11) = -(Tmax*lvx)/(m^2*(lvx^2 + lvy^2 + lvz^2)^(1/2));
F(14, 12) = -(Tmax*lvy)/(m^2*(lvx^2 + lvy^2 + lvz^2)^(1/2));
F(14, 13) = -(Tmax*lvz)/(m^2*(lvx^2 + lvy^2 + lvz^2)^(1/2));

Jacobian = F;

end


function grad = objectiveJacobian(tfinal, statefinal, STM, data)


% STM:
STM_rrlr = STM(1:3, 8:10);
STM_rrlv = STM(1:3, 11:13);
STM_rrlm = STM(1:3, 14);


STM_vvlr = STM(4:6, 8:10);
STM_vvlv = STM(4:6, 11:13);
STM_vvlm = STM(4:6, 14);



STM_lmlr = STM(14, 8:10);
STM_lmlv = STM(14, 11:13);
STM_lmlm = STM(14, 14);

% Final Dynamics:
dx = dynamics(tfinal, statefinal, data);
dr = dx(1:3);
dv = dx(4:6);
dlm = dx(14);

% Hamiltonian:
dxx = dx(1:7);
dll = dx(8:14);
STM_xxll = STM(1:7, 8:14);
STM_llll = STM(8:14, 8:14);
dHf = -dll' * STM_xxll + dxx' * STM_llll;

% Assemble Gradients:

grad = [STM_rrlr STM_rrlv STM_rrlm dr;
        STM_vvlr STM_vvlv STM_vvlm dv;
        STM_lmlr STM_lmlv STM_lmlm, dlm
        dHf, 0];

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
set(0, 'defaultLegendLocation','best');
set(0, 'defaultLegendOrientation', 'vertical');
set(0, 'defaultLegendFontSize', 20);
% Setting Axes:
set(0, 'defaultAxesXMinorGrid','on');
set(0,'defaultAxesYMinorGrid','on');
set(0, 'defaultAxesFontSize', 30);


end
