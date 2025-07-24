function [G_np, om_low, om_high, f] = nonParamFRFest(r,u,y,Ts,sys_real_red)
% Function to estimate the non-parametric FRF of the real system from time
% domain I/O data

% Define welch estimator parameters:
na = 12;                        % averaging factor
nx = length(r);
window = hanning(floor(nx/na));
noverlap = 1024;
nfft = 2048;                      
fs = 1/Ts; 


% Auto-spectra of excitation signal and output:
[PHI_rr, f] = pwelch(r, window, noverlap, nfft, fs);
[PHI_yy, ~] = pwelch(y, window, noverlap, nfft, fs);

% Cross-spectra:
[PHI_ur, ~] = cpsd(u, r, window, noverlap, nfft, fs);
[PHI_yr, ~] = cpsd(y, r, window, noverlap, nfft, fs);
[PHI_ry, ~] = cpsd(r, y, window, noverlap, nfft, fs);

S = PHI_ur ./ PHI_rr;       
SG = PHI_yr ./ PHI_rr;     
G_np = SG ./ S;            


% Compute model FRF (magnitude and phase)
[H_model, ~] = freqresp(sys_real_red, 2*pi*f);
H_model = squeeze(H_model); 

mag_model = abs(H_model);
phase_model = angle(H_model);  
% Non-parametric estimate
mag_np = abs(G_np);
phase_np = angle(G_np); 


% Coherence Function Estimation:
gamma2_ry = abs(PHI_ry).^2 ./ (PHI_rr .* PHI_yy);

% Find bandwidth with gamma > 0.6:
threshold = 0.6;
idx_band_coh = gamma2_ry > threshold;
f_band_coh = f(idx_band_coh);

% Estimate frequency range
om_low = min(f_band_coh) * 2*pi;
om_high = max(f_band_coh) * 2*pi;


% Plot
figure;
% --- Magnitude subplot ---
subplot(2,1,1);
semilogx(2*pi*f, 20*log10(mag_model), 'b-'); hold on;
semilogx(2*pi*f, 20*log10(mag_np),    '-.', 'Color', [0.9 0.2 0.3]);
grid on;
xlabel('Frequency [rad/s]');
ylabel('Magnitude [dB]');
xlim([2*pi*f(2) 200])
xline(om_low, '--k', 'LineWidth', 1)
xline(om_high,  '--k', 'LineWidth', 1)
legend('$G_{\text{real}(j\omega)$', '$\hat{G}_m(j\omega)$ (Welch Method)', '$\Omega_{\text{range}}$');

% --- Phase subplot ---
subplot(2,1,2);
semilogx(2*pi*f, rad2deg(unwrap(phase_model)), 'b-'); hold on;
semilogx(2*pi*f, rad2deg(unwrap(phase_np)), '-.', 'Color', [0.9 0.2 0.3]);
grid on;
xlabel('Frequency [rad/s]');
ylabel('Phase [deg]');
xline(om_low, '--k', 'LineWidth', 0.8)
xline(om_high,  '--k', 'LineWidth', 0.8)
xlim([2*pi*f(2) 200])
ylim([-200 10])

% Coherence Function:
figure;
semilogx(2*pi*f, gamma2_ry,  '-.s', ...
    'LineWidth', 1, ...
    'MarkerSize', 4, ...
    'MarkerFaceColor', 'w'); hold on;
xline(om_low, '--k', 'LineWidth', 1)
xline(om_high,  '--k', 'LineWidth', 1)
xlabel('Frequency [rad/s]', 'FontSize', 26);
ylabel('$\gamma^2_{\delta,q}$', 'FontSize', 26);
xlim([2*pi*f(2) 5*10^2])
legend('', '$\Omega_{range}$', 'FontSize', 26)
grid on;
