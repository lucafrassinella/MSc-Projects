%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation and Learning in Aerospace Project A.Y. 24-25 
% TASK 2: Black-Box System Identification
%
% Authors:  Luca Frassinella (luca.frassinella@polimi.it)
%           Gianluca Scudier (gianluca.scudier@polimi.it)                                                  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;
close all;
bdclose('all');
clc;
addpath('common','common/simulator-toolbox','common/simulator-toolbox/attitude_library','common/simulator-toolbox/trajectory_library');
addpath('Functions T2');
format short;
plotSettings;


%% Real Model parameters

% Initial model 
%   - state: longitudinal velocity, pitch rate, pitch angle; 
%   - input: normalised pitching moment; 
%   - outputs: state and longitudinal acceleration;

Xu = -0.1068;
Xq = 0.1192;
Mu = -5.9755;
Mq = -2.6478;
Xd = -10.1647;
Md = 450.71;

A = [Xu, Xq, -9.81; Mu, Mq, 0; 0, 1, 0];
B = [Xd; Md; 0];
C = [1, 0, 0; 0, 1, 0; 0, 0, 1; Xu, Xq, 0]; 
D = [0; 0; 0; Xd];

realModel.A = A;
realModel.B = B;
realModel.C = C;
realModel.D = D;
real_parameters = [Xu; Xq; Mu; Mq; Xd; Md];
realModel.parameters = real_parameters;

% Noise
% noise.Enabler = 0;
noise.Enabler = 1;
noise.pos_stand_dev = noise.Enabler * 0.0011;                            	%[m]
noise.vel_stand_dev = noise.Enabler * 0.01;                               %[m/s]
noise.attitude_stand_dev = noise.Enabler * deg2rad(0.33);                 %[rad]
noise.ang_rate_stand_dev = noise.Enabler * deg2rad(1);                   %[rad/s]

seed.x = 1;
seed.vx = 2;
seed.theta = 3;
seed.q = 4;

% Delays
delay.position_filter = 1;
delay.attitude_filter = 1;
delay.mixer = 1;

% Load controller parameters

ctrl = parameters_controller(); 
Ts = ctrl.sample_time;

% Load Excitation Signal:
load ExcitationM_eq.mat
ExcitationM = ExcitationM_eq;

% selected interval:
t = ExcitationM(:,1);
simulation_time = t(end) - t(1);

sys_real = ss(A,B,C(2,:),D(2));
                  

%% Simulation:

% Simulation of the longitudinal dynamics:
simulation = sim('Simulator_Single_Axis.slx', 'SrcWorkspace', 'current');

% Extract data from simulation:
data.q = simulation.q;
data.Mtot = simulation.Mtot;

% Input-Output data:
u = data.Mtot;  % Input [u(t)]: Total pitch moment (normalized)
y = data.q;
r = ExcitationM(:, 2);
 
% Time
time = 0 : ctrl.sample_time : simulation_time;

%%%%%%%%%%%%%% Task 2: PBSID-H_inf Identification
[G_np, om_low, om_high, f_np] = nonParamFRFest(r,u,y,Ts,sys_real); 

%%% PBSID with (p = p_opt, f = f_opt, n = n_opt):
%   - script to choose {p,f,n}_opt: 'tuning_pbsid.m' 

% Optimal p, f, n:
plotColorMap_pfn();

p_opt = 225;
f_opt = 210;
n_opt = 5;

p = p_opt;
f = f_opt;
n = n_opt;

% Remove Delay from Input-Output:
u = u(1:end-1);
y = y(2:end);


[Ai, Bi, Ci, Di, Ki, S, X] = PBSID(u, y, p, f, n);

% Plot Singular Values:
figure();
semilogy(S, '*b', 'LineWidth', 1); hold on;
semilogy(n, S(n), 'or', 'MarkerSize', 8, 'LineWidth', 1); 
xlabel('Order index');
ylabel('Singular value');
title('Singular Values');
grid on;
xlim([0 f_opt+1])
legend('Singular Values', 'Optimal Model Order');

% D2C of Estimated PBSID Matrices with Tustin method:
sys_PBSID_DT = ss(Ai, Bi, Ci, Di, Ts);
sys_hat_CT = d2c(sys_PBSID_DT, 'tustin');


%%% H-infinity Model Matching (systune + hinfstruct):

% Load and Plot analysis of (om1,om2,N opt)
plotOptimalFilter()

% Optimal filter parameters:
om1_opt = 3.3;      % [rad/s]
om2_opt = 10.6;     % [rad/s]
N_opt = 12;         % [-] Filter Order

om1 = om1_opt;
om2 = om2_opt;
N = N_opt;

[theta_hat, theta0] = modelMatching(sys_hat_CT, om1, om2, N);

% Estimated Model:
g = 9.81;
estModel.A_est = [theta_hat(1) theta_hat(2) -g; theta_hat(3) theta_hat(4) 0; 0 1 0];
estModel.B_est = [theta_hat(5) theta_hat(6) 0]';
estModel.C_est = [1 0 0; 0 1 0; 0 0 1; theta_hat(1), theta_hat(2), 0];
estModel.D_est = [0; 0; 0; theta_hat(3)];

estModel.parameters = theta_hat;

% From MC Analysis on M random unstructured PBSID models, load variances:
%   NB: refer to varianceAnalysis_T2.m
[MC_mean, MC_std] = AsVar_PBSID();

%%% --------------- T2.5: Validation with 3211 signal ----------------- %%%

% Design suitable 3211 sequence:
ampl = 0.1;
[ExcitationM, simulation_time, T] = design3211(estModel, ampl, Ts);

N = length(ExcitationM(:, 1));

% Simulate Real System:
seed.x = 1;
seed.vx = 2;
seed.theta = 3;
seed.q = 4;
val_real = sim('Simulator_Single_Axis.slx', 'SrcWorkspace', 'current');

% Extract data from simulation:
val.q_m = val_real.q;
val.ax_m = val_real.ax;
val.theta_m = val_real.theta;
val.theta0 = val_real.theta0;
val.Mtot = val_real.Mtot;
val.t_vec = val_real.t_vec;

% Simulate Estimated System:
val_est = sim('Validate_Simulator_Single_Axis.slx', 'SrcWorkspace', 'current');

% Extract data from simulation:
val.q_est = val_est.q;
val.ax_est = val_est.ax;
val.theta_est = val_est.theta;

% Validation Results:
[FIT, PEC, VAF] = validation(val, ExcitationM(:, 1), N);


%% Final Plots:
sys_est = ss(estModel.A_est, estModel.B_est, estModel.C_est(2,:), estModel.D_est(2));
w = logspace(log10(0.1), log10(200), 500);

% Uncertain parameters:
pvec = estModel.parameters;      
dcov = 3 * MC_std; 
Xu_u = ureal('Xu',pvec(1),'PlusMinus',[-dcov(1),dcov(1)]);
Xq_u = ureal('Xq',pvec(2),'PlusMinus',[-dcov(2),dcov(2)]);
Mu_u = ureal('Mu',pvec(3),'PlusMinus',[-dcov(3),dcov(3)]);
Mq_u = ureal('Mq',pvec(4),'PlusMinus',[-dcov(4),dcov(4)]);
Xd_u = ureal('Xd',pvec(5),'PlusMinus',[-dcov(5),dcov(5)]);
Md_u = ureal('Md',pvec(6),'PlusMinus',[-dcov(6),dcov(6)]);

% Uncertain Model:
Au = [Xu_u, Xq_u, -9.81; Mu_u, Mq_u, 0; 0, 1, 0];
Bu = [Xd_u; Md_u; 0];
Cu = [0, 1, 0];  
Du = 0;
usys = ss(Au, Bu, Cu, Du);

% Sample +/- 3sigma bounds:
N = 100;
sampled_usys = usample(usys, N);
mag_all = zeros(N, length(w));
pha_all = zeros(N, length(w));

for i = 1:N
    [mag, pha] = bode(sampled_usys(:,:,i), w);
    mag_all(i,:) = squeeze(mag); 
    pha_all(i,:) = squeeze(pha); 
end

mag_min = min(mag_all); mag_max = max(mag_all);
pha_min = min(pha_all); pha_max = max(pha_all);

% Nominal Bode:
[mag_est, pha_est] = bode(sys_est, w);
[mag_real, pha_real] = bode(sys_real, w);
mag_est = squeeze(mag_est); pha_est = squeeze(pha_est);
mag_real = squeeze(mag_real); pha_real = squeeze(pha_real);
% Non Parametric FRF:
mag_np = abs(G_np); phase_np = angle(G_np); 

% Plot:
figure()

subplot(2,1,1);
semilogx(w, 20*log10(mag_real), 'k'); hold on;
semilogx(w, 20*log10(mag_est), '-', 'Color', [0.8 0.1 0.2]);
semilogx(2*pi*f_np, 20*log10(mag_np), '-.', 'Color', [0 0.5 0], 'LineWidth', 0.8);
semilogx(w, 20*log10(mag_min), ':', 'Color', [0.8 0.4 0.5]);
semilogx(w, 20*log10(mag_max), ':', 'Color', [0.8 0.4 0.5]);
xlim([0.5 200])
ylabel('Magnitude [dB]');
grid on;

subplot(2,1,2);
semilogx(w, pha_real, 'k'); hold on;
semilogx(w, pha_est, '-', 'Color', [0.8 0.1 0.2]);
semilogx(2*pi*f_np, rad2deg(phase_np), '-.', 'Color', [0 0.5 0], 'LineWidth', 0.8);
semilogx(w, pha_min, ':', 'Color', [0.8 0.4 0.5]);
semilogx(w, pha_max, ':', 'Color', [0.8 0.4 0.5]);
ylabel('Phase [deg]');
xlabel('Frequency [rad/s]');
xlim([0.5 200])
ylim([-100 10])
grid on;
legend('Real Model', 'Estimated Model', 'Non-Parametric FRF (Welch)', '$3\sigma$');

% Pole/Zero Map:
z_est = tzero(sys_est);
p_est = pole(sys_est);
z_real = tzero(sys_real);
p_real = pole(sys_real);

figure();
plot(real(z_est), imag(z_est), 'o', 'MarkerSize', 10, 'LineWidth', 2, 'Color', [0.8 0.1 0.2]); hold on;
plot(real(p_est), imag(p_est), 'x', 'MarkerSize', 10, 'LineWidth', 2, 'Color', [0.8 0.1 0.2]);
plot(real(z_real), imag(z_real), 'ko', 'MarkerSize', 8, 'LineWidth', 1.5);
plot(real(p_real), imag(p_real), 'kx', 'MarkerSize', 8, 'LineWidth', 1.5);
axis equal;
xlabel('Real Axis [s$^{-1}$])', 'Interpreter','latex');
ylabel('Imaginary Axis [s$^{-1}$])', 'Interpreter','latex');
title('Pole/Zero Map', 'Interpreter','latex');
legend('Zeros - Est.', 'Poles - Est.', 'Zeros - Real', 'Poles - Real', ...
    'Location','best', 'Interpreter','latex');
grid on;
xlim([-4.5 4.5])
ylim([-3.5 3.5])

%% End of code
