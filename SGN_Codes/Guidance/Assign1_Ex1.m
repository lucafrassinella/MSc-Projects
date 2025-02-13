% Spacecraft Guidance and Navigation - A.Y. 2024/25
% Assignment #1
% Exercise # 1
% Author: Frassinella Luca - 10795356

%% --------------------- EX. 1 - PERIODIC ORBITS ----------------------- %%

clearvars; close all; clc;
format long
plotSettings;

%% -------------------- Pt.1 - Find Lagrangian Points -------------------- 
mu = 0.012150;
constants.mu = mu;
% Compute Lagrangian Points (coordinates and Jacobi Constant):
[lagrangePoints.L1.coords, lagrangePoints.L1.C] = findLagrangePoint(mu, initialGuess('L1'), 'Collinear');
[lagrangePoints.L2.coords, lagrangePoints.L2.C] = findLagrangePoint(mu, initialGuess('L2'), 'Collinear');
[lagrangePoints.L3.coords, lagrangePoints.L3.C] = findLagrangePoint(mu, initialGuess('L3'), 'Collinear');
[lagrangePoints.L4.coords, lagrangePoints.L4.C] = findLagrangePoint(mu, initialGuess('L4'), 'Triangular');
[lagrangePoints.L5.coords, lagrangePoints.L5.C] = findLagrangePoint(mu, initialGuess('L5'), 'Triangular');

% Plots:
plotLagrangePoints(mu, lagrangePoints)

% Results:
results.ex1 = lagrangePoints;

%% ----------------------- Pt.2 - Find Halo Orbit ------------------------ 

% Initial Conditions:
x0 = 1.068792441776;
y0 = 0;
z0 = 0.071093328515;
vx0 = 0; 
vy0 = 0.319422926485;
vz0 = 0; 
data.initialState = [x0; y0; z0; vx0; vy0; vz0];

% Propagate forward initial orbit:
[xx_fwd, ~] = propagate(0, data.initialState, 5, mu, false);
% Propagate backwards initial orbit
[xx_bwd, ~] = propagate(0, data.initialState, -5, mu, false);
% Plot Initial Halo Orbit:
figure();
plot3(xx_fwd(:, 1), xx_fwd(:, 2), xx_fwd(:, 3), 'HandleVisibility', 'off')
hold on
plot3(xx_bwd(:, 1), xx_bwd(:, 2), xx_bwd(:, 3), 'HandleVisibility', 'off')
plot3(lagrangePoints.L2.coords(1), lagrangePoints.L2.coords(2), lagrangePoints.L2.coords(3), ...
    'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r')
text(lagrangePoints.L2.coords(1), lagrangePoints.L2.coords(2), lagrangePoints.L2.coords(3) - 0.012, ...
    'L2', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', 20, 'FontWeight', 'bold')
grid on
xlabel('x [-]', 'FontSize', 30)
ylabel('y [-]', 'FontSize', 30)
zlabel('z [-]', 'FontSize', 30)
view([45 15])
title('Initial Propagated Trajectory (@EMB Earth-Moon Rotating Frame)')

% Define parameters and initialize vectors for Pseudo-Newton Method:
settings.Nmax = 1000;
settings.tol = 1e-12;
data.tfGuess = 5;
data.initialErr = ones(4, 1);

% Target values {y = 0; vx = 0; vz = 0; C = 3.09}:
refC = 3.09;
data.desiredValues = [0; 0; 0; refC];

% Find correct Initial State that satisfies periodicity and energy constraints:
[correctedState, newC, finalErr, tf, iter]  = findHaloOrbit(data, constants, settings);

% Update results structure:
results.ex2.correctedState = correctedState;
results.ex2.tf = tf;

% Propagate corrected orbit:
[xx, tt] = propagate(0, correctedState, tf, mu, true);

% Plot Halo Orbit:
figure();
plot3(xx(:, 1), xx(:, 2), xx(:, 3), 'Color', [0.1043 0.3843 0.6062], 'HandleVisibility', 'off')
hold on
plot3(xx(:, 1), -xx(:, 2), xx(:, 3), 'Color', [0.1043 0.3843 0.6062], 'HandleVisibility', 'off')
plot3(lagrangePoints.L2.coords(1), lagrangePoints.L2.coords(2), lagrangePoints.L2.coords(3), ...
    'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r')
text(lagrangePoints.L2.coords(1), lagrangePoints.L2.coords(2), lagrangePoints.L2.coords(3) - 0.012, ...
    'L2', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', 25, 'FontWeight', 'bold')
grid on
xlabel('x [-]', 'FontSize', 30)
ylabel('y [-]', 'FontSize', 30)
zlabel('z [-]', 'FontSize', 30)
view([40 15])
title('Periodic Halo Orbit (@EMB Earth-Moon Rotating Frame)')

%% ---------------------- Pt.3 - Halo Orbit Family ----------------------- 

% Set vector of target Jacobi constants:
cvec = 3.09 : -0.005 : 3.04;
% Set matrix of target values of {y, vx, vz, C}:
desiredValues_vec = [zeros(3, length(cvec)); cvec];
% Set initial state as the corrected state of Pt.2:
data.initialState = correctedState; 
% Set initial errors and parameter for Halo orbits finder:
data.initialErr = ones(4, 1);
data.tfGuess = 5;
settings.Nmax = 1000;
settings.tol = 1e-12;

corrStates = zeros(6, length(cvec));    % initialize corrected states matrix
newCvec = zeros(length(cvec), 1);       % initialize Jacobi constant vector
errvec = zeros(4, length(cvec));        % initialize error matrix
iterVec = zeros(length(cvec), 1);       % initialize iterations vector


% Begin plot of family of Halo orbits:
figure()
title('Family of Halo orbits (@EMB Earth-Moon Rotating Frame)')
hold on
grid on
xlabel('x [-]', 'FontSize', 30)
ylabel('y [-]', 'FontSize', 30)
zlabel('z [-]', 'FontSize', 30)
plot3(1-mu, 0, 0, 'o', 'MarkerSize', 16, 'MarkerFaceColor', [0.5 0.5 0.5], ...
    'DisplayName', 'Moon')
plot3(lagrangePoints.L2.coords(1), lagrangePoints.L2.coords(2), lagrangePoints.L2.coords(3), ...
    'ro', 'MarkerSize', 14, 'MarkerFaceColor', 'r', 'DisplayName', 'L2')
legend('FontSize', 20);

% Compute Halo orbits exploiting numerical continuation:
for i = 1 : length(cvec)
    data.desiredValues = desiredValues_vec(:, i); % extract vector of target values
    [corrState, newC, err, tfNew, iter]  = findHaloOrbit(data, constants, settings); % compute corrected state for periodic Halo orbit
    corrStates(:, i) = corrState;   % update corrected states matrix
    newCvec(i) = newC;              % update Jacobi constant vector
    errvec(:, i) = err;        % update error matrix
    iterVec(i) = iter;              % update iterations vector

    % propagate and plot Halo orbit:
    [xx, tt] = propagate(0, corrState, tfNew, mu, true);
    
    scatter3(xx(:,1), xx(:,2), xx(:,3), 3, newC*ones(size(xx, 1), 1), 'filled', 'HandleVisibility', 'off');
    scatter3(xx(:,1), -xx(:,2), xx(:,3), 3, newC*ones(size(xx, 1), 1), 'filled', 'HandleVisibility', 'off');
    
    data.initialState = corrState;            % update new initial state as output of last step
end

results.ex3.correctedState = corrState;
results.ex3.tf = tfNew;

% add colorbar to plot:
cb = colorbar;
cb;
clim([min(newCvec) max(newCvec)]);
ylabel(cb, 'Jacobi Constant')
colormap('abyss')
view([35 15])

%% Clear Workspace
clearvars -except constants data settings results

%% Functions:

function [coords, jacobiConst] = findLagrangePoint(mu, guess, pointType)
% ----------------------------------------------------------------------- %
% findLagrangePoint - function to compute coordinates of Earth-Moon
% lagrange points and their associated Jacobi Constant
%
% Inputs:
%       - mu: reduced gravitational parameter                  [1, 1]
%       - guess: initial guess for each of the Lagrange Points [3, 1]
%       - pointType: string to identify type of Lagrange Point
%               - 'Collinear' for L1, L2, L3
%               - 'Triangular' for L4, L5
%
% Outputs:
%       - coords: coordinates of the Lagrange Points in the Earth-Moon
%       rotating frame (position)                              [3, 1]
%       - jacobiConst: Jacobi Constant of each Lagrange Point  [1, 1]
% ----------------------------------------------------------------------- %


if strcmpi(pointType, 'triangular')
    % Define the derivatives wrt x and y:
    dUdx = @(Lxy) Lxy(1) - (1-mu)*(Lxy(1)+mu)./((Lxy(1)+mu).^2 + Lxy(2).^2).^(3/2) - mu*(Lxy(1)-(1-mu))./((Lxy(1)-(1-mu)).^2 + Lxy(2).^2).^(3/2);
    dUdy = @(Lxy) Lxy(2) - (1-mu)*Lxy(2)./((Lxy(1)+mu).^2 + Lxy(2).^2).^(3/2) - mu*Lxy(2)./((Lxy(1)-(1-mu)).^2 + Lxy(2).^2).^(3/2);
    triangularFnctn = @(Lxy) [dUdx(Lxy); dUdy(Lxy)];
    optset = optimoptions('fsolve','FunctionTolerance', 1e-12, 'OptimalityTolerance', 1e-12, 'Display', 'none', 'Algorithm', 'levenberg-marquardt');
    coords = fsolve(triangularFnctn, guess, optset);
elseif strcmpi(pointType, 'collinear')
    % Define the derivatives wrt x:
    dUdx = @(Lx) Lx - (1-mu)*(Lx+mu)./abs(Lx+mu).^3 - mu*(Lx+mu-1)./abs(Lx+mu-1).^3;
    optset = optimoptions('fsolve','FunctionTolerance', 1e-12, 'OptimalityTolerance', 1e-12, 'Display', 'none');
    Lx = fsolve(dUdx, guess(1), optset);
    coords = [Lx; 0; 0]; % append y = 0 and z = 0
end

% Compute State:
state = [coords; 0; 0; 0]; % vx, vy, vz = 0 for eq. point in rotating frame

% Compute Jacobi Constant:
jacobiConst = JacobiConstant(mu, state);

end

function guess = initialGuess(lagrangePoint)
% ----------------------------------------------------------------------- %
% initialGuess - Function to import initial guesses for fsolve solver
%
% Inputs:
%       - lagrangePoint: string to identify Lagrange Point 
%
% Outputs:
%       - guess: initial guess for each Lagrange Point [3, 1]
% ----------------------------------------------------------------------- %

if strcmpi(lagrangePoint, 'L1')
    guess = [0.5; 0; 0];
elseif strcmpi(lagrangePoint, 'L2')
    guess = [1; 0; 0];
elseif strcmpi(lagrangePoint, 'L3')
    guess = [-1; 0; 0];
elseif strcmpi(lagrangePoint, 'L4')
    guess = [0; 1; 0];
elseif strcmpi(lagrangePoint, 'L5')
    guess = [0; -1; 0];
end
end

function C = JacobiConstant(mu, xx)
% ----------------------------------------------------------------------- %
% JacobiConstant - function to compute the Jacobi Constant associated to a
% given state
%
% Inputs:
%       - mu: reduced gravitational parameter                  [1, 1]
%       - xx: actual state (position, velocity)                [6, 1]
%
% Outputs:
%       - C: Jacobi Constant associated to the state xx        [1, 1]
% ----------------------------------------------------------------------- %

% Compute r1, r2:
r1 = norm([xx(1) + mu; xx(2); xx(3)]);
r2 = norm([xx(1) + mu - 1; xx(2); xx(3)]);
% Extract velocity norm:
v = norm(xx(4:6));
% Potential:
OM = 0.5*(xx(1)^2 + xx(2)^2) + (1 - mu)/r1 + mu/r2 + 0.5*mu*(1 - mu); 

% Compute Jacobi Constant:
C = 2*OM - v^2;

end

function [dxdt] = CR3BP(~,mu, xx)
% ----------------------------------------------------------------------- %
% CR3BP - function to evaluate the RHS for the equations of motion of the 
% Circular Restricted 3-Body Problem (CR3BP)
%
% Inputs:
%       - t: current time (here, unitilized)                   [1, 1]
%       - mu: reduced gravitational parameter                  [1, 1]
%       - xx: current state (position, velocity)               [6, 1]
%
% Outputs:
%       - dxdt: Time derivative of the state vector            [6, 1]
% ----------------------------------------------------------------------- %

% Extract variables:
x = xx(1);
y = xx(2);
z = xx(3);
vx = xx(4);
vy = xx(5);

% STM in Matrix Form:
PHI = reshape(xx(7: end), 6, 6);

% Distances from bodies 1 and 2:
r1 = sqrt((x + mu)^2 + y^2 + z^2);
r2 = sqrt((x + mu - 1)^2 + y^2 + z^2);
% Derivative of the Potential:
dUdx = x - (1-mu)/r1^3*(mu+x) + mu/r2^3*(1-mu-x);
dUdy = y - (1-mu)/r1^3*y - mu/r2^3*y;
dUdz = - (1-mu)/r1^3 * z - mu/r2^3 * z;

% Assemble Matrix A = dfdx (6x6):
A = zeros(6);
A(1, 4) = 1;
A(2, 5) = 1;
A(3, 6) = 1;
A(4, 1) = (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + (3*mu*(2*mu + 2*x - 2)*(mu + x - 1))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*(2*mu + 2*x)*(mu + x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2)) + 1;
A(4, 2) = (3*mu*y*(mu + x - 1))/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*(mu + x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2);
A(4, 3) = (3*mu*z*(mu + x - 1))/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*z*(mu + x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2);
A(4, 5) = 2;
A(5, 1) = (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2));
A(5, 2) = (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) + 1;
A(5, 3) = (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2);
A(5, 4) = -2;
A(6, 1) = (3*mu*z*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2));
A(6, 2) = (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2);
A(6, 3) = (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*z^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*z^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2);

% Compute derivative of the STM:
PHIdot = A * PHI;

% Assemble dxdt:
dxdt = zeros(42, 1);
dxdt(1:3) = xx(4: 6);
dxdt(4) = dUdx + 2*vy;
dxdt(5) = dUdy - 2*vx;
dxdt(6) = dUdz;
dxdt(7:end) = PHIdot(:);

end

function [gradJ] = gradientJ(mu, xx)
% ----------------------------------------------------------------------- %
% gradJ - function to evaluate the gradient of the Jacobi Constant given a
% current state xx
%
% Inputs:
%       - mu: reduced gravitational parameter                  [1, 1]
%       - xx: current state (position, velocity)               [6, 1]
%
% Outputs:
%       - gradJ: gradient of the Jacobi Constant J             [6, 1]
% ----------------------------------------------------------------------- %

% Extract Variables:
x = xx(1);
y = xx(2);
z = xx(3);
vx = xx(4);
vy = xx(5);
vz = xx(6);

% Compute gradient:
gradJ = zeros(6, 1);
gradJ(1) = 2*x + ((2*mu + 2*x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2) - (mu*(2*mu + 2*x - 2))/((mu + x - 1)^2 + y^2 + z^2)^(3/2);
gradJ(2) = 2*y - (2*mu*y)/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + (2*y*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2);
gradJ(3) = (2*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2) - (2*mu*z)/((mu + x - 1)^2 + y^2 + z^2)^(3/2);
gradJ(4) = -2*vx;
gradJ(5) = -2*vy;
gradJ(6) = -2*vz;

end

function [xx, tt] = propagate(t0, xx0, tf, mu, evtFlag)
% ----------------------------------------------------------------------- %
% propagate - function to propagate the CR3BP dynamics over a given
% time-span
%
% Inputs:
%       - t0: initial propagation time                         [1, 1]
%       - xx0: initial state (position, velocity)              [6, 1]
%       - tf: final propagation time                           [1, 1]
%       - mu: reduced gravitational parameter                  [1, 1]
%       - evtFlag: boolean vector to stop the integration when event y=0 is
%       reached
%               - true: stops the integration when event is reached
%               - false: doesn't stop the integration when event is reached
%
% Outputs:
%       - xx: propagated state vector and STM at each discrete time [n, 42]
%       - tt: discrete time vector                                  [n, 1]
% ----------------------------------------------------------------------- %

TOF = tf - t0;
PHI0 = eye(6); % STM at t0
xx0 = [xx0; PHI0(:)]; % new vector of initial conditions (with STM)

optset = odeset('RelTol', 3e-14, 'AbsTol', 1e-14, 'Events', @(t, xx)planeCrossing(t, xx, evtFlag));
[tt, xx] = ode78(@(t,xx) CR3BP(t, mu, xx), [0 TOF], xx0, optset);

end

function [value, isterminal, direction] = planeCrossing(~,xx, evtFlag)
% ----------------------------------------------------------------------- %
% planeCrossing - Event function to stop the integration when the xz plane
% is reached
% Inputs:
%   - t: actual time (here, unutilized)              [1, 1]
%   - xx: current state vector                       [6, 1]
%   - evtFlag: [integer] Specifies if the integration should stop when the 
%   event is detected.
%
% Outputs:
%   - value: Value of the y-coordinate, used to detect the crossing of the xz-plane.
%   - isterminal: [integer] Indicates if the integration should terminate (1 = stop, 0 = continue).
%   - direction: [integer] Direction of crossing (0 = any direction, 1 = positive, -1 = negative).
% ----------------------------------------------------------------------- %

    value = xx(2);
    isterminal = evtFlag;
    direction = 0;
end 

function [correctedState, newC, finalErr, tf, iter]  = findHaloOrbit(data, constants, settings)
% ----------------------------------------------------------------------- %
% findHaloOrbit - function to compute periodic Halo Orbit using iterative
% differential correction scheme
%
% Inputs:
%       - data: struct containing data of the problems:
%           - data.initialState -> initial state guess         [6, 1]
%           - data.tfGuess      -> initial time guess          [1, 1]
%           - data.initialErr   -> initialized error vector    [4, 1]
%           - data.desiredValues -> target values              [4, 1]
%       - constants: struct containing constants of the problem (mu)
%       - settings: struct containing settings for the problem (Nmax, tol)
%
% Outputs:
%       - correctedState: Final corrected state vector         [6, 1]
%       - newC: final obtained Jacobi Constant                 [1, 1]
%       - finalErr: final errors wrt target values             [4, 1]
%       - tf: corrected final time                             [1, 1]
%       - iter: number of iterations to reach convergence      [1, 1]
% ----------------------------------------------------------------------- %

% Extract constants, data and settings:
mu = constants.mu;
initialState = data.initialState;
tfGuess = data.tfGuess;
initialErr = data.initialErr;
desiredValues = data.desiredValues;
Nmax = settings.Nmax;
tol = settings.tol;

% Extract desired values:
yRef = desiredValues(1);
vxRef = desiredValues(2);
vzRef = desiredValues(3);
refC = desiredValues(4);

% Initialize iter, newGuess and error to begin the Pseudo-Newton Method:
iter = 0;
initialGuess = [initialState(1); initialState(3); initialState(5); tfGuess];
newGuess = initialGuess;
err = initialErr;
% Initialize new states:
xx0new = [newGuess(1); initialState(2); newGuess(2); initialState(4); newGuess(3); initialState(6); reshape(eye(6), 36, 1)]; % State and STM (42x1)
newState = xx0new(1:6);

% Main Loop:
while (abs(err(1))>tol || abs(err(2))>tol || abs(err(3))>tol || abs(err(4))>tol) && iter < Nmax
     if iter 
         % Assemble reduced 4x4 linear system:
         PHI = reshape(actualSTM, 6, 6);
         A_red = [PHI(2,1), PHI(2, 3), PHI(2, 5), phidot(2);
                  PHI(4,1), PHI(4, 3), PHI(4, 5), phidot(4);
                  PHI(6,1), PHI(6, 3), PHI(6, 5), phidot(6);
                  gradC(1), gradC(3), gradC(5), 0];
         b_red = err;
         corr = A_red\b_red;
         newGuess = newGuess - corr;
         xx0new = [newGuess(1); initialState(2); newGuess(2); initialState(4); newGuess(3); initialState(6); reshape(eye(6), 36, 1)]; % State and STM (42x1)
         newState = xx0new(1:6);
         
     end
    %  Propagate to final point:
    [xx, tt] = propagate(0, newState, newGuess(4), mu, true);
    xf = xx(end, :);
    te = tt(end);
    
    actualState = xf(1:6); % Extract state 
    actualSTM = xf(7:end); % Extract STM 

    newC = JacobiConstant(mu, newState); % compute new Jacobi Constant
    % Update matrices of coefficients:
    [xDot] = CR3BP(te, mu, [actualState'; reshape(eye(6), 36, 1)]); % actual final state + STM
    phidot = xDot(1:6);
    gradC = gradientJ(mu, xx0new);
    % Compute errors
    err = [actualState(2) - yRef; actualState(4) - vxRef; actualState(6) - vzRef; newC - refC];

    iter = iter + 1;
end

% Final Results:
xx0_corr = xx0new;
correctedState = xx0_corr(1: 6);
newC = JacobiConstant(mu, correctedState);
finalErr = err;
tf = newGuess(4);

end

function [] = plotLagrangePoints(mu, lagrangePoints)
% ----------------------------------------------------------------------- %
% plotLagrangePoints - function to plot Lagrange Points od the Earth-Moon
% system
%
% Inputs:
%       - lagrangePoints: structure containing the coordinates of each
%         Lagrange Point 
%       - mu: reduced gravitational parameter                  [1, 1]
%
% Outputs:
%       - []
% ----------------------------------------------------------------------- %

% Begin plot:
figure()
title('Lagrange Points - Earth-Moon System')
hold on
grid on
xlabel('x [-]', 'FontSize', 25)
ylabel('y [-]', 'FontSize', 25)

% Plot circumferences of r = 1 and centered in P1 (Earth) and P2 (moon)
theta = linspace(0, 2*pi, 100);  r = 1;
xEarth = -mu + r * cos(theta);  yEarth = r * sin(theta); 
xMoon = (1 - mu) + r * cos(theta);  yMoon = r * sin(theta); 
plot(xEarth, yEarth, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8, 'HandleVisibility', 'off')
plot(xMoon, yMoon, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.6, 'HandleVisibility', 'off')

% Plot lines of dUdx:
xvec = linspace(-1.1, 2.1, 1000);
dUdx = xvec - (1-mu)*(xvec+mu)./abs(xvec+mu).^3 - mu*(xvec+mu-1)./abs(xvec+mu-1).^3;
dUdx(dUdx > 15 | dUdx < -15) = NaN;   % to avoid vertical lines where the function is discontinuos
plot(xvec, dUdx, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6, 'HandleVisibility', 'off')

% Set coordinates and names of points to plot:
points = {
    -mu, 0, 'Earth';                                                 % Earth
    1-mu, 0, 'Moon';                                                 % Moon
    lagrangePoints.L1.coords(1), lagrangePoints.L1.coords(2), 'L1';  % L1
    lagrangePoints.L2.coords(1), lagrangePoints.L2.coords(2), 'L2';  % L2
    lagrangePoints.L3.coords(1), lagrangePoints.L3.coords(2), 'L3';  % L3
    lagrangePoints.L4.coords(1), lagrangePoints.L4.coords(2), 'L4';  % L4
    lagrangePoints.L5.coords(1), lagrangePoints.L5.coords(2), 'L5'   % L5
};

% Set Colors and markerSizes
colors = {[0 0.2 1], [0.5 0.5 0.5], 'r', 'r', 'r', 'r', 'r'};
markerSizes = [16, 13, 10, 10, 10, 10, 10];

% Plot and add text labels:
for i = 1:2
    plot(points{i, 1}, points{i, 2}, 'o', 'MarkerSize', markerSizes(i), ...
        'MarkerFaceColor', colors{i}, 'DisplayName', points{i, 3});
    text(points{i, 1}, points{i, 2} - 0.1, ['\bf{' points{i, 3} '}'], ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', 16);
end
for i = 3:7
    plot(points{i, 1}, points{i, 2}, 'o', 'MarkerSize', markerSizes(i), ...
        'MarkerFaceColor', colors{i}, 'DisplayName', points{i, 3});
    text(points{i, 1}, points{i, 2} + 0.15, ['\bf{' points{i, 3} '}'], ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', 16);
end
axis equal
xlim([-1.1 2.1])
ylim([-1.2 1.2])
hold off

end

function plotSettings
%-------------------------------------------------------------------------%
% Settings for figure plots
%-------------------------------------------------------------------------%

% Setting Lines:
set(0, 'defaultLineLineWidth', 1.6);
set(0,'defaultLineMarkerSize', 4) ;
set(0,'defaultLineMarkerEdgeColor', 'auto')
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
