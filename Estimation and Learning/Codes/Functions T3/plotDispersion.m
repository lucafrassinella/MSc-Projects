function plotDispersion(THETA_T1, THETA_T2, theta_T1, theta_T2, realParameters)

% Setup:
N = size(THETA_T1, 2);
w = logspace(log10(0.1), log10(200), 500);


% Nominal Systems:
[A_T1,B_T1,C_T1,D_T1] = theta2abcd(theta_T1, 9.81);
sys_nom_T1 = ss(A_T1,B_T1,C_T1,D_T1);
[mag_nom_T1, pha_nom_T1] = bode(sys_nom_T1, w);
mag_nom_T1 = squeeze(mag_nom_T1);
pha_nom_T1 = squeeze(pha_nom_T1);
p_nom_T1 = pole(sys_nom_T1);
z_nom_T1 = tzero(sys_nom_T1);

[A_T2,B_T2,C_T2,D_T2] = theta2abcd(theta_T2, 9.81);
sys_nom_T2 = ss(A_T2,B_T2,C_T2,D_T2);
[mag_nom_T2, pha_nom_T2] = bode(sys_nom_T2, w);
mag_nom_T2 = squeeze(mag_nom_T2);
pha_nom_T2 = squeeze(pha_nom_T2);
p_nom_T2 = pole(sys_nom_T2);
z_nom_T2 = tzero(sys_nom_T2);

% Real System:
[A_r,B_r,C_r,D_r] = theta2abcd(realParameters, 9.81);
sys_r = ss(A_r,B_r,C_r,D_r);
[mag_real, pha_real] = bode(sys_r, w);
mag_real = squeeze(mag_real);
pha_real = squeeze(pha_real);
p_r = pole(sys_r);
z_r = tzero(sys_r);

% Preallocate:
mag_T1_all = zeros(N, length(w));
mag_T2_all = zeros(N, length(w));
pha_T1_all = zeros(N, length(w));
pha_T2_all = zeros(N, length(w));


for i = 1:N
    % T1
    [A,B,C,D] = theta2abcd(THETA_T1(:,i), 9.81);
    sys = ss(A,B,C,D);
    [mag, pha] = bode(sys, w);
    mag_T1_all(i,:) = squeeze(mag);
    pha_T1_all(i,:) = squeeze(pha);

    % T2
    [A,B,C,D] = theta2abcd(THETA_T2(:,i), 9.81);
    sys = ss(A,B,C,D);
    [mag, pha] = bode(sys, w);
    mag_T2_all(i,:) = squeeze(mag);
    pha_T2_all(i,:) = squeeze(pha);
end


% Pole-Zero Dispersion Plots:
poles_MC_T1 = zeros(3*N, 1);
zeros_MC_T1 = zeros(2*N, 1);
poles_MC_T2 = zeros(3*N, 1);
zeros_MC_T2 = zeros(2*N, 1);

% T1, T2:
k = 1;
j = 1;
for i = 1:size(THETA_T1,2)
    [A,B,C,D] = theta2abcd(THETA_T1(:,i), 9.81);
    sys_T1 = ss(A,B,C,D);
    poles_MC_T1(k:(k+2), 1) = pole(sys_T1);
    zeros_MC_T1(j:(j+1), 1) = tzero(sys_T1);

    [A,B,C,D] = theta2abcd(THETA_T2(:,i), 9.81);
    sys_T2 = ss(A,B,C,D);
    poles_MC_T2(k:(k+2), 1) = pole(sys_T2);
    zeros_MC_T2(j:(j+1), 1) = tzero(sys_T2);
    k = k + 3;
    j = j + 2;
end

% Plot Poles/Zeros T1:
figure();
hold on;
plot(real(zeros_MC_T1), imag(zeros_MC_T1), 'o', 'Color', [0.6 0.6 0.6], 'MarkerSize', 8, 'LineWidth', 1);
plot(real(poles_MC_T1), imag(poles_MC_T1), 'x', 'Color', [0.6 0.6 0.6], 'MarkerSize', 8, 'LineWidth', 1);
plot(real(z_nom_T1), imag(z_nom_T1), 'o', 'MarkerSize', 10, 'LineWidth', 1.5, 'Color', [0.2 0.1 0.8]);
plot(real(p_nom_T1), imag(p_nom_T1), 'x', 'MarkerSize', 10, 'LineWidth', 1.5, 'Color', [0.2 0.1 0.8]);
plot(real(z_r), imag(z_r), 'ko', 'MarkerSize', 8, 'LineWidth', 1.5);
plot(real(p_r), imag(p_r), 'kx', 'MarkerSize', 8, 'LineWidth', 1.5);
axis equal;
xlim([-4.5 4.5]);
ylim([-4.25 4.25]);
legend('MC Zeros','MC Poles','Nominal Zeros','Nominal Poles','Real Zeros','Real Poles', ...
    'Location','bestoutside');
xlabel('Real Axis [$s^{-1}$]'); ylabel('Imaginary Axis [$s^{-1}$]');
title('Pole/Zero Disperion (T1)');
grid on; box on;

% Plot Poles/Zeros T2:
figure();
hold on;
plot(real(zeros_MC_T2), imag(zeros_MC_T2), 'o', 'Color', [0.6 0.6 0.6], 'MarkerSize', 8, 'LineWidth', 1);
plot(real(poles_MC_T2), imag(poles_MC_T2), 'x', 'Color', [0.6 0.6 0.6], 'MarkerSize', 8, 'LineWidth', 1);
plot(real(z_nom_T2), imag(z_nom_T2), 'o', 'MarkerSize', 10, 'LineWidth', 1.5, 'Color', [0.8 0.1 0.2]);
plot(real(p_nom_T2), imag(p_nom_T2), 'x', 'MarkerSize', 10, 'LineWidth', 1.5, 'Color', [0.8 0.1 0.2]);
plot(real(z_r), imag(z_r), 'ko', 'MarkerSize', 8, 'LineWidth', 1.5);
plot(real(p_r), imag(p_r), 'kx', 'MarkerSize', 8, 'LineWidth', 1.5);
axis equal;
xlim([-4.5 4.5]);
ylim([-4.25 4.25]);
legend('MC Zeros','MC Poles','Nominal Zeros','Nominal Poles','Real Zeros','Real Poles', ...
    'Location','bestoutside');
xlabel('Real Axis [$s^{-1}$]'); ylabel('Imaginary Axis [$s^{-1}$]');
title('Pole/Zero Disperion (T2)');
grid on; box on;


% Bode Dispersion Plots:

% T1:
figure();
title('Bode Disperion (T1)');
subplot(2,1,1);
semilogx(w, 20*log10(mag_T1_all(1,:)), 'Color', [0.8 0.8 0.8], 'LineWidth', 4); hold on;
for i = 2:size(mag_T1_all,1)
    semilogx(w, 20*log10(mag_T1_all(i,:)), 'Color', [0.8 0.8 0.8],  'LineWidth', 4, 'HandleVisibility', 'off');
end
semilogx(w, 20*log10(mag_nom_T1), 'LineWidth', 1.5, 'Color', [0.2 0.1 0.8]);
semilogx(w, 20*log10(mag_real), '--k', 'LineWidth', 1.2);
ylabel('Magnitude [dB]');
xlim([0.1 200]); ylim([-20 55]);
grid on;

subplot(2,1,2);
semilogx(w, pha_T1_all(1,:), 'Color', [0.8 0.8 0.8], 'LineWidth', 4); hold on;
for i = 2:size(pha_T1_all,1)
    semilogx(w, pha_T1_all(i,:), 'Color', [0.8 0.8 0.8], 'LineWidth', 4, 'HandleVisibility', 'off');
end
semilogx(w, pha_nom_T1, 'LineWidth', 1.5, 'Color', [0.2 0.1 0.8]);
semilogx(w, pha_real, '--k', 'LineWidth', 1.2);
xlabel('Frequency [rad/s]', 'Interpreter', 'latex');
ylabel('Phase [deg]');
xlim([0.1 200]); ylim([-95 5]);
grid on;
legend('Monte Carlo','Nominal (T1)', 'Real System', ...
    'Location','best');

% T2:
figure();
title('Bode Disperion (T2)');
subplot(2,1,1);
semilogx(w, 20*log10(mag_T2_all(1,:)), 'Color', [0.8 0.8 0.8], 'LineWidth', 1.3); hold on;
for i = 2:size(mag_T2_all,1)
    semilogx(w, 20*log10(mag_T2_all(i,:)), 'Color', [0.8 0.8 0.8],  'LineWidth', 1.3, 'HandleVisibility', 'off');
end
semilogx(w, 20*log10(mag_nom_T2), 'LineWidth', 1.5, 'Color', [0.8 0.1 0.2]);
semilogx(w, 20*log10(mag_real), '--k', 'LineWidth', 1.2);
ylabel('Magnitude [dB]');
xlim([0.1 200]); ylim([-20 55]);
grid on;

subplot(2,1,2);
semilogx(w, pha_T2_all(1,:), 'Color', [0.8 0.8 0.8], 'LineWidth', 1.3); hold on;
for i = 2:size(pha_T2_all,1)
    semilogx(w, pha_T2_all(i,:), 'Color', [0.8 0.8 0.8], 'LineWidth', 1.3, 'HandleVisibility', 'off');
end
semilogx(w, pha_nom_T2, 'LineWidth', 1.5, 'Color', [0.8 0.1 0.2]);
semilogx(w, pha_real, '--k', 'LineWidth', 1.2);
xlabel('Frequency [rad/s]', 'Interpreter', 'latex');
ylabel('Phase [deg]');
xlim([0.1 200]); ylim([-95 5]);
grid on;
legend('Monte Carlo','Nominal (T1)', 'Real System', ...
    'Location','best');


end
