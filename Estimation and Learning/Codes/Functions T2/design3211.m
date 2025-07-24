function [ExcitationM, simulation_time, T] = design3211(id, ampl, Ts)
% Function to design a 3211 excitation signal having as time scale half the
% period of the dominant mode of the estimated system
% Inputs:
%   - id: estimated model structure
%   - ampl: amplitude of 3211 signal

% Extract estimated A and C matrices:
A_est = id.A_est;  
C_est = id.C_est; 

c_q = C_est(2,:);  

% Eigs:
[V, D_mat] = eig(A_est);      
poles = diag(D_mat);      

% Modal Contribution:
modal_contrib = abs(c_q * V);  

% Dominant Pole:
[~, idx_dom] = max(modal_contrib);
lambda_dom = poles(idx_dom);

% Natural Frequency and Period:
sigma = real(lambda_dom);
omega_d = imag(lambda_dom);
omega_n = sqrt(sigma^2 + omega_d^2);
Tn = 2*pi / omega_n; 

% Time scale of 3211:
T = Tn / 2;        

% 3211 Signal:
fs = 1 / Ts;                     
t_total = 7*T;                   
t = 0:1/fs:t_total;

u3211 = zeros(size(t));
u3211(t <= 3*T) = +ampl;
u3211(t > 3*T& t <= 5*T) = -ampl;
u3211(t > 5*T & t <= 6*T) = ampl;
u3211(t > 6*T & t <= 7*T) = -ampl;

ExcitationM = zeros(length(t), 2);
ExcitationM(:,1) = t';
ExcitationM(:,2) = u3211';

simulation_time = t(end);
