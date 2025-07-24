%%% PBSID with Cross-Validation

clearvars;
close all;
bdclose('all');
clc;
format short;
plotSettings;

%% Real Model parameters
% Initial model 
%   - state: longitudinal velocity, pitch rate, pitch angle; 
%   - input: normalised pitching moment; 
%   - outputs: state and longitudinal acceleration;

% Load controller parameters
ctrl = parameters_controller();
Ts = ctrl.sample_time;

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

% Simulation of the longitudinal dynamics:

load ExcitationM_eq.mat
ExcitationM = ExcitationM_eq;

% Selected interval:
t = ExcitationM(:,1);
simulation_time = t(end) - t(1);

simulation = sim('Simulator_Single_Axis.slx', 'SrcWorkspace', 'current');

% Extract data from simulation:
data.q = simulation.q;
data.Mtot = simulation.Mtot;

% Input-Output data:
u = data.Mtot;  
y = data.q;
r = ExcitationM(:, 2);


% Time
time = 0 : ctrl.sample_time : simulation_time;


% Estimate Non-Parametric FRF and define Frequency Range with Coherence:

sys_real_red = ss(realModel.A, realModel.B, realModel.C(2,:), realModel.D(2));

[G_np, om_low, om_high, fvec] = nonParamFRFest(r,u,y,Ts,sys_real_red);


%% Cross Validation (p, f, n): (First GridSearch --> HP: p = f, n = 1:5);

% Remove delay from outputs:
u = u(1:end-1);
y = y(2:end);

% Ranges
p_range = 10:5:250;
n_range = 1:5;

% Restrict matching of FRF only in valid band (coherence > 0.6):
idx_band = (2*pi*fvec >= 2*om_low) & (2*pi*fvec <= om_high);
om_range = 2*pi*fvec(idx_band);              
G_np_band = G_np(idx_band);            
mag_np_band = abs(G_np_band);

% Configurations (p = f):
configs = [];
for p = p_range
    for n = n_range
        configs = [configs; p, n];
    end 
end


% Preallocate:
numConfigs = size(configs,1);
max_bodemag_errors = zeros(numConfigs,1);
N = numConfigs;

G_all = zeros(length(om_vec), length(p_range));
parfor k = 1 : N
    % Parametri
    p = configs(k,1);
    f = p;
    n = configs(k,2);

    % PBSID
    [Ai, Bi, Ci, Di, ~, ~, ~] = PBSID(u, y, p, f, n);

    
    % Convert Estimated SS system to CT:
    sys_PBSID_DT = ss(Ai, Bi, Ci, Di, Ts);
    sys_PBSID_CT = d2c(sys_PBSID_DT, 'tustin');

    % Store Bodemag of estimated system (for plot):
    [G_pbsid_ct, ~] = freqresp(sys_PBSID_CT, om_vec);
    G_pbsid_ct = squeeze(G_pbsid_ct);
    mag_pbsid_ct = abs(G_pbsid_ct);
    G_all(:,k) = mag_pbsid_ct;  
    
    
    % Compute bodemag of sys_PBSID at the frequencies of om_range:
    [mag_PBSID, ~] = bode(sys_PBSID_CT, om_range);
    mag_PBSID = squeeze(mag_PBSID);  
    

    % Magnitude error:
    err_mag = abs(20*log10(mag_PBSID) - 20*log10(mag_np_band));
    max_bodemag_errors(k) = max(err_mag);
    
end

% Best Configuration:
[minERR, idx] = min(max_bodemag_errors);
best_config = configs(idx,:);
best_bodemag = G_all(:,idx);
fprintf('\nBest config: p=%d, f=%d, n=%d --> minErr=%.2f\n', ...
best_config(1), best_config(1), best_config(2), minERR);

% Fix n to the optimal value computed in the first grid search:
n = best_config(2);

%% Plot Bode Comparisons (Grid Search #1):

% Plot Non-Parametric FRF vs Estimated FRF:
om_vec_range = [0.5 100];
idx_plot = (2*pi*fvec >= om_vec_range(1)) & (2*pi*fvec <= om_vec_range(end));
om_vec = 2*pi*fvec(idx_plot);            % frequenze per il plot in rad/s
G_np_plot = G_np(idx_plot);              % corrispondente FRF stimata in om_vec
mag_np_plot = abs(G_np_plot);

figure();
semilogx(om_vec, 20*log10(G_all(:,1)), '--', 'Color', [0.8 0.8 0.8], 'LineWidth', 1, 'HandleVisibility', 'off');
hold on;
for k = 2:N
    semilogx(om_vec, 20*log10(G_all(:,k)), '--', 'Color', [0.8 0.8 0.8], 'LineWidth', 1, 'HandleVisibility', 'off');
end
% Optimal FRF (PBSID):
semilogx(om_vec, 20*log10(best_bodemag), 'k', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('Optimal PBSID ($p=f=%d$)', best_config(1)));
semilogx(om_vec, 20*log10(mag_np_plot),'-.', 'Color', [0 0.5 0], 'LineWidth', 1.4, 'DisplayName', 'Non-Parametric FRF (Welch)');
semilogx(om_vec, NaN(size(om_vec)), '--', 'Color', [0.8 0.8 0.8], ...
    'LineWidth', 1.2, ...
    'DisplayName', 'PBSID (for $n=5$, $p=f=10{:}250$)');
xline(2*om_low, '--k', 'LineWidth', 1, 'HandleVisibility', 'off')
xline(om_high, '--k', 'LineWidth', 1, 'HandleVisibility', 'off')
grid on;
xlabel('Frequency [rad/s]',  'Interpreter', 'latex');
ylabel('Magnitude [dB]',     'Interpreter', 'latex');
xlim([0.76 100]);
legend('Interpreter', 'latex', 'FontSize', 18, 'Location', 'northeast');


% Plot max error between FRF at different p = f:
figure()
plot(p_range, max_bodemag_errors, 's', ...
    'MarkerSize', 9, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', [0.6 0.5 0.8], ...
    'LineStyle', 'none');
grid on;
xlabel('Past and future horizon $p = f$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\max (|G_{\mathrm{np}}(j\omega) - \hat{G}(j\omega)|)$ [dB]', ...
    'Interpreter', 'latex', 'FontSize', 18);
ylim([0.7 2])
title('Max. Frequency response error (for $n=5$, $p=f=10{:}250$)')


%% Cross Validation (p, f, n): (Second GridSearch --> HP: p =\= f, n = n_opt)
clc;
% Ranges
p_range = 10: 5 : 250;
f_range = 10: 5 : 250;

idx_band = (2*pi*fvec >= 2*om_low) & (2*pi*fvec <= om_high);
om_range = 2*pi*fvec(idx_band);              
G_np_band = G_np(idx_band);              
mag_np_band = abs(G_np_band);

% Combinations: (p <= f)
configs = [];
for p = p_range
    for f = f_range
        if f <= p
            configs = [configs; p, f];
        end
     end
end

% Preallocazione
numConfigs = size(configs,1);
max_bodemag_errors = zeros(numConfigs,1);
N = numConfigs;


parfor k = 1 : N
    % Parameters:
    p = configs(k,1);
    f = configs(k,2);

    % PBSID
    [Ai, Bi, Ci, Di, ~, ~, ~] = PBSID(u, y, p, f, n);
    
    % Estimated SS system:
    sys_PBSID_DT = ss(Ai, Bi, Ci, Di, Ts);
    
    % Compute bodemag of sys_PBSID at the frequencies of om_range:
    [mag_PBSID, ~] = bode(sys_PBSID_DT, om_range);
    mag_PBSID = squeeze(mag_PBSID);  
    
    % Magnitude error:
    err_mag = abs(20*log10(mag_PBSID) - 20*log10(mag_np_band));
    max_bodemag_errors(k) = max(err_mag);
    
end

% Best config
[minERR, idx] = min(max_bodemag_errors);
best_config = configs(idx,:);
fprintf('\nBest config: p=%d, f=%d, n=%d --> minErr=%.2f\n', ...
best_config(1), best_config(2), n, minERR);

p_opt = best_config(1);
f_opt = best_config(2);
n_opt = n;

%% SAVE RESULTS:
% save('Mats/PBSID_FrequencyGridSearch.mat', ...
%     'max_bodemag_errors', 'configs', 'p_range', 'f_range', 'p_opt', 'f_opt', 'n_opt');

%% Colormap (for n = n_opt, p, f in p_range and f_range):

P = length(p_range);
F = length(f_range);
ErrGrid = nan(P, F); 

for i = 1:numConfigs
    p_val = configs(i,1);
    f_val = configs(i,2);

    
    row = find(p_range == p_val);
    col = find(f_range == f_val);

    ErrGrid(row, col) = max_bodemag_errors(i);
end

% NaN where f > p:
for i = 1:P
    for j = 1:F
        if f_range(j) > p_range(i)
            ErrGrid(i,j) = NaN;
        end
    end
end

% Mesh grid
[F_grid, P_grid] = meshgrid(f_range, p_range);

% Plot con pcolor
figure();
h = pcolor(F_grid, P_grid, ErrGrid, 'HandleVisibility', 'off');
set(h, 'EdgeColor', "none");      
set(gca, 'YDir', 'normal');
colormap('parula');
colorbar;
alpha_data = ~isnan(ErrGrid);      
set(h, 'AlphaData', alpha_data);

% colorbar:
c = colorbar;
ylabel(c, '$\max(|G_{\mathrm{err}}|)$ [dB]', 'Interpreter', 'latex', 'FontSize', 20);
xlabel('Future horizon ($f$)', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('Past horizon ($p$)', 'Interpreter', 'latex', 'FontSize', 22);
title('PBSID$_\mathrm{opt}$: Frequency-domain error (focus band)', ...
       'Interpreter', 'latex', 'FontSize', 25);

% Highlight best configuration
hold on;
plot(best_config(2), best_config(1), 'wo', ...
     'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k' , 'MarkerSize', 10, ...
     'DisplayName', sprintf('Optimum ($p=%d$, $f=%d$)', best_config(1), best_config(2)));
legend('show', 'Location', 'northeast', 'Interpreter', 'latex', 'fontsize', 18);