%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation and Learning in Aerospace Project A.Y. 24-25 
% TASK 3: Monte Carlo Analysis
%
% Authors:  Luca Frassinella (luca.frassinella@polimi.it)
%           Gianluca Scudier (gianluca.scudier@polimi.it)                                                  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;
clc;
addpath('Mats\');
addpath('Functions T3\')
plotSettings;

%%
% Load Results
load('Results T1.mat')
theta_T1 = est_parameters_ref;
estModel_T1 = estModel;
modelEval_T1 = modelEval;
ssPctg_T1 = ss_pctg;

load('Results T2.mat')
theta_T2 = theta_hat;
estModel_T2 = estModel;

% Clear Workspace:
clear est_parameters est_parameters_ref est_parameters_raw estModel modelEval ss_pctg theta_hat

%% T3.1 - Compare Models Identified in T1 and T2
clc;
close all;
w = logspace(log10(0.1), log10(200), 500);

% Real System:
sys_real = ss(realModel.A, realModel.B, realModel.C(2,:), realModel.D(2));
z_real = tzero(sys_real);
p_real = pole(sys_real);
[mag_real, pha_real] = bode(sys_real, w);
mag_real = squeeze(mag_real);
pha_real = squeeze(pha_real);

% T1:
sys_T1 = ss(estModel_T1.A_est, estModel_T1.B_est, estModel_T1.C_est(2,:), estModel_T1.D_est(2));
z_T1 = tzero(sys_T1);
p_T1 = pole(sys_T1);
[mag_T1, pha_T1] = bode(sys_T1, w);
mag_T1 = squeeze(mag_T1);
pha_T1 = squeeze(pha_T1);

% T2:
sys_T2 = ss(estModel_T2.A_est, estModel_T2.B_est, estModel_T2.C_est(2,:), estModel_T2.D_est(2));
z_T2 = tzero(sys_T2);
p_T2 = pole(sys_T2);
[mag_T2, pha_T2] = bode(sys_T2, w);
mag_T2 = squeeze(mag_T2);
pha_T2 = squeeze(pha_T2);

% Plot ZP:
figure();

plot(real(p_T1), imag(p_T1), 'x', 'MarkerSize', 10, 'LineWidth', 2, 'Color', [0.2 0.1 0.8]); hold on;
plot(real(p_T2), imag(p_T2), 'x', 'MarkerSize', 10, 'LineWidth', 2, 'Color', [0.8 0.1 0.2]);
plot(real(p_real), imag(p_real), 'kx', 'MarkerSize', 8, 'LineWidth', 1.5);
plot(real(z_T1), imag(z_T1), 'o', 'MarkerSize', 10, 'LineWidth', 2, 'Color', [0.2 0.1 0.8]);
plot(real(z_T2), imag(z_T2), 'o', 'MarkerSize', 10, 'LineWidth', 2, 'Color', [0.8 0.1 0.2]);
plot(real(z_real), imag(z_real), 'ko', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Zeros - Real');
xlabel('Re'); ylabel('Im');
xlim([-5 5]);
ylim([-4.25 4.25]);
title('Pole/Zero Map');
legend('Poles (T1)', 'Poles (T2)', 'Real Poles', 'Zeros (T1)', 'Zeros (T2)', 'Real Zeros');
axis equal; grid on;

% Plot FRF:
figure();
% Magnitude
subplot(2,1,1);
semilogx(w, 20*log10(mag_T1), 'LineWidth', 1.2, 'Color', [0.2 0.1 0.8]); hold on;
semilogx(w, 20*log10(mag_T2), 'LineWidth', 1.2, 'Color', [0.8 0.1 0.2]);
semilogx(w, 20*log10(mag_real), 'k', 'LineWidth', 1);
ylabel('Magnitude [dB]');
xlim([0.5 200]);
ylim([-30 50]);
legend('Est. (T1)', 'Est. (T2)', 'Real');
grid on;

% Phase
subplot(2,1,2);
semilogx(w, pha_T1, 'LineWidth', 1.2, 'Color', [0.2 0.1 0.8]); hold on;
semilogx(w, pha_T2, 'LineWidth', 1.2, 'Color', [0.8 0.1 0.2]); hold on;
semilogx(w, pha_real, 'k', 'LineWidth', 1);
xlim([0.5 200]);
ylim([-95 5]);
ylabel('Phase [deg]');
xlabel('Frequency [rad/s]');
grid on;


% Hinf error (only in focus band):
om_low = 2.3;
om_high = 22;
idx = w > om_low & w < om_high;

mag_T1_dB   = 20*log10(mag_T1);
mag_T2_dB   = 20*log10(mag_T2);
mag_real_dB = 20*log10(mag_real);

diff_T1_dB = abs(mag_T1_dB - mag_real_dB);
diff_T2_dB = abs(mag_T2_dB - mag_real_dB);


% H-infinity error in omega range:
Hinf_T1 = max(diff_T1_dB(idx));
Hinf_T2 = max(diff_T2_dB(idx));

%% T3.3 --> Montecarlo results study

load('MC_T1.mat');
load('MC_T2.mat');

%%% Histograms of estimated parameters:
[mu_T1, std_T1, ssPct_T1] = mc_parameters(THETA_T1, 'T1');
[mu_T2, std_T2, ssPct_T2] = mc_parameters(THETA_T2, 'T2');

%%% Histograms of poles/zeros:
[mu_re_p_T1, sigma_re_p_T1, mu_im_p_T1, sigma_im_p_T1, mu_z_T1, sigma_z_T1] = mc_pz(THETA_T1, 'T1');
[mu_re_p_T2, sigma_re_p_T2, mu_im_p_T2, sigma_im_p_T2, mu_z_T2, sigma_z_T2] = mc_pz(THETA_T2, 'T2');

%%% Dispersions of poles/zeros and FRFs:
plotDispersion(THETA_T1, THETA_T2, theta_T1, theta_T2, realModel.parameters);

