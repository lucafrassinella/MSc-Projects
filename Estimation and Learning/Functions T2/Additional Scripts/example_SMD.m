clear;
close all;
clc;


%% Test Case: System of two masses connected by springs and dampers

% Define parameters:
m1 = 20; m2 = 20; k1 = 1000; k2 = 2000; d1 = 1; d2 = 5;

% Define time parameters:
time_steps = 300;
t_sample = 0.05;
Ts = 0.05; 
% Define number of inputs and outputs (r, m)
r = 1; m = 1;

% Define continuous-time system
Ac = [0, 1, 0, 0; 
      -(k1+k2)/m1, -(d1+d2)/m1, k2/m1, d2/m1;
      0, 0, 0, 1;
      k2/m2, d2/m2, -k2/m2, -d2/m2];
Bc = [0; 0; 0; 1/m2];
Cc = [1, 0, 0, 0];
Dc = 0;

% Discrete-time system (Backward Euler discretization)
I = eye(size(Ac));
A = inv(I - t_sample * Ac);       % backward Euler A
B = A * t_sample * Bc;            % backward Euler B
C = Cc;
D = 0;


% Input & Initial state (for identification):
rng(42)
x0_id = randn(4,1);
u_id = randn(1, time_steps);

% Input & Initial state (for validation):
rng(99)
x0_val = randn(4,1);
u_val = randn(1, time_steps);

% Simulate discrete-time systems to obtain input-output data for
% identification and validation

[Y_id, X_id] = simulate(A,B,C, u_id, x0_id);
[Y_val, X_val] = simulate(A,B,C, u_val, x0_val);

% Plot:
figure()
plot(u_id)
title('Input sequence (identification)')

figure()
plot(Y_id)
title('true output (identification)')

%% PBSIDopt for system identification

p = 24;
f = 20;
n = 4; 

[Ai, Bi, Ci, Di, Ki, S, X] = PBSID(u_id, Y_id, p, f, n);

figure;
semilogy(S, 'o');


%% Validation:

h = 10;  % window size to estimate initial state

% Estimate the initial state for the validation set
x0_est = estimateInitial(Ai, Bi, Ci, u_val, Y_val, h);

% Simulate the open-loop model
[Y_val_prediction, ~] = simulate(Ai, Bi, Ci, u_val, x0_est);

% Compute validation metrics
[rel_err, VAF, Akaike] = modelError(Y_val, Y_val_prediction, r, m, 30);
fprintf('Final model relative error: %.2f%%\nVAF value: %.2f%%\nAkaike: %.2f\n', rel_err, VAF, Akaike);


% Plot: real vs predicted output (first 100 samples)
figure;
plot(1:100, Y_val(1,1:100), 'k', 'DisplayName', 'Real output'); hold on;
plot(1:100, Y_val_prediction(1,1:100), 'r--', 'DisplayName', 'Prediction');
legend;
xlabel('t');
ylabel('y');
title('Model Validation');
grid on;

%% Check eigs:

eig_est = eig(Ai);
eig_true = eig(A);

figure;
plot(real(eig_true), imag(eig_true), 'bo', 'MarkerSize', 10, 'LineWidth', 2); hold on;
plot(real(eig_est), imag(eig_est), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
legend('Real System', 'PBSIDopt');
xlabel('Re'); ylabel('Im');
title('Eigenvalues Comparison');
axis equal; grid on;

% 1) Crea i modelli SS discreti (assicurati di passare Ts)
sys_estimated = ss(Ai, Bi, Ci, Di, Ts);
sys_true      = ss(A,    B,    C,    D,    Ts);


%%%% STRUCTURING METHOD (hinfstruct):


% PBSID in CT:
I = eye(size(Ai));
A_est_c = (I - inv(Ai)) / Ts;
B_est_c = inv(Ai) * Bi / Ts;
C_est_c = Ci;
D_est_c = Di;

sys_hat_c = ss(A_est_c, B_est_c, C_est_c, D_est_c); 
sys_hat_c.InputName = 'u';
sys_hat_c.OutputName = 'y';

% Known parameters:
m1 = 20; m2 = 20; k1 = 1000; k2 = 2000;

% Unknown parameters:
d1 = realp('d1', 0.9);  
d1.Minimum = 0.2;  d1.Maximum = 3;
d2 = realp('d2', 4);  
d2.Minimum = 4.5;  d2.Maximum = 5.5;


A_c = [  0,        1,         0,       0;
       -(k1+k2)/m1, -(d1+d2)/m1, k2/m1, d2/m1;
        0,         0,         0,       1;
        k2/m2,     d2/m2,    -k2/m2, -d2/m2 ];
B_c = [0; 0; 0; 1/m2];
C_c = [1 0 0 0];
D_c = 0;

sys_th = ss(A_c, B_c, C_c, D_c);
sys_th.InputName = 'u';
sys_th.OutputName = 'y';

M = sys_hat_c - sys_th; 

opt = hinfstructOptions('Display', 'final', 'RandomStart', 10);
[sys_opt, gamma7, info] = hinfstruct(M, opt);


d1_hat = sys_opt.Blocks.d1.Value;
d2_hat = sys_opt.Blocks.d2.Value;

A_hat = [  0,        1,         0,       0;
       -(k1+k2)/m1, -(d1_hat+d2_hat)/m1, k2/m1, d2_hat/m1;
        0,         0,         0,       1;
        k2/m2,     d2_hat/m2,    -k2/m2, -d2_hat/m2 ];
B_hat = [0; 0; 0; 1/m2];
C_hat = [1 0 0 0];
D_hat = 0;



%% Final Plot:

sys_real_c  = ss(Ac, Bc, Cc, Dc);
sys_hat = ss(A_hat, B_hat, C_hat, D_hat);

figure()
bode(sys_real_c, 'k'); hold on;
bode(sys_hat, '--r');
legend('Real','Estimated')


%% Functions:

function [Y, X] = simulate(A, B, C, U, x0)
    % A, B, C: sistema in forma di stato
    % U: input (dimensione r x N)
    % x0: stato iniziale

    simTime = size(U, 2);     
    n = size(A, 1);           
    r = size(C, 1);           % number of outputs

    X = zeros(n, simTime + 1);
    Y = zeros(r, simTime);

    for i = 1:simTime
        if i == 1
            X(:, i) = x0;
            Y(:, i) = C * x0;
            X(:, i + 1) = A * x0 + B * U(:, i);
        else
            Y(:, i) = C * X(:, i);
            X(:, i + 1) = A * X(:, i) + B * U(:, i);
        end
    end
end

function [rel_err_pct, vaf_pct, Akaike_error] = modelError(Ytrue, Ypredicted, r, m, n)
% This function computes the prediction performance of an estimated model
%
% Inputs:
%   Ytrue       - true output (r x N)
%   Ypredicted  - predicted output (r x N)
%   r           - number of outputs
%   m           - number of inputs
%   n           - estimated system order
%
% Outputs:
%   rel_err_pct     - relative error in percentage
%   vaf_pct         - variance accounted for in percentage
%   Akaike_error    - Akaike information criterion estimate

    timeSteps = size(Ytrue, 2);
    total_parameters = n * (n + m + 2 * r);

    % Error matrix
    error_matrix = Ytrue - Ypredicted;

    % Flatten to column vectors (column-major)
    Ytrue_vec = reshape(Ytrue, [], 1);
    Ypred_vec = reshape(Ypredicted, [], 1);
    error_vec = Ytrue_vec - Ypred_vec;

    % Relative error %
    rel_err_pct = norm(error_vec, 2) / norm(Ytrue_vec, 2) * 100;

    % VAF %
    vaf_pct = 100 * (1 - (1/timeSteps) * norm(error_vec)^2 / ((1/timeSteps) * norm(Ytrue_vec)^2));
    vaf_pct = max(vaf_pct, 0);  % clip to zero if negative

    % Akaike error estimate
    cov_matrix = (1 / timeSteps) * (error_matrix * error_matrix');
    Akaike_error = log(det(cov_matrix)) + (2 / timeSteps) * total_parameters;
end

function x0_est = estimateInitial(A, B, C, U, Y, h)
% This function estimates an initial state x0 of the model:
%   x_{k+1} = A x_k + B u_k
%   y_k     = C x_k
%
% using the input and output sequences (u_k, y_k) for k = 0 to h-1
%
% Inputs:
%   A, B, C : system matrices
%   U       : m x N input matrix
%   Y       : r x N output matrix
%   h       : estimation window (number of samples used)
%
% Output:
%   x0_est  : estimated initial state (n x 1)

    n = size(A, 1);  % state dimension
    r = size(C, 1);  % output dimension
    m = size(U, 1);  % input dimension

    % Stack output Y from 0 to h-1 (column-wise)
    Y_0_hm1 = reshape(Y(:,1:h), [], 1);  % (h*r x 1)

    % Stack input U from 0 to h-1 (column-wise)
    U_0_hm1 = reshape(U(:,1:h), [], 1);  % (h*m x 1)

    % Build observability matrix O_{h-1}
    O_hm1 = zeros(h*r, n);
    for i = 0:h-1
        O_hm1(i*r+1:(i+1)*r, :) = C * (A^i);
    end

    % Build input-influence matrix I_{h-1}
    I_hm1 = zeros(h*r, h*m);
    for i = 1:h-1
        for j = 0:i-1
            I_hm1(i*r+1:(i+1)*r, j*m+1:(j+1)*m) = C * (A^(i-j-1)) * B;
        end
    end

    % Estimate initial state
    x0_est = pinv(O_hm1) * (Y_0_hm1 - I_hm1 * U_0_hm1);
end
