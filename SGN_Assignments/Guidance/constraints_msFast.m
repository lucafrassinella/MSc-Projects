function [c, ceq, dC, dCeq] = constraints_msFast(vars, data, constants)

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
    [~, xx, PHI] = propagate_stm(xxj, t(j), t(j+1), constants);

    % Compute continuity constraint (matching positions and velocities):
    ceq(4*j-3 : 4*j) = xx(end, 1:4)' - xxjj;

    % Update inequality constraints:
    xj = xxjj(1); 
    yj = xxjj(2); 
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



%%%%%%% GRADIENTS %%%%%%%%%%


% 
% for j = 1:N-1
%     % Extract current interval state;
%     xxj = vars(4*j-3 : 4*j); 
%     xj = xxj(1);  
%     yj = xxj(2);  
% 
%     % Propagate state and get RHS of the PBRFBP:
%     [~, xx, PHI] = propagate_stm(xxj, t(j), t(j+1), constants); 
%     fj = PBRFBP(t(j), xxj, constants, 0);  
%     fjj = PBRFBP(t(j+1), xx(end, :), constants, 0);  
% 
%     % Compute equality constraint gradients (matching positions and
%     % velocities):
%     Q1(4*j-3 : 4*j) = -(N-j)/(N-1) * PHI * fj + (N-j-1)/(N-1) * fjj;  
%     QN(4*j-3 : 4*j) = -(j-1)/(N-1) * PHI * fj + j/(N-1) * fjj;  
% 
%     % Update Jacobians for equality constraints:
%     dCeq(4*j-3: 4*j, 4*j+1:  4*j + 4) = -eye(4);  
%     dCeq(4*j-3: 4*j, 4*j-3 : 4*j) = PHI;  
% 
%     % Update Jacobians for inequality constraints:
%     dC(2*j-1: 2*j, 4*j-3:  4*j) = [-2*(xj + mu), -2*yj, 0, 0;  
%                                    -2*(xj + mu - 1), -2*yj, 0, 0]; 
% end

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