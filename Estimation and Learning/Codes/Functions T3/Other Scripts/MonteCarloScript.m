%% Monte Carlo Task 1:

clearvars;
close all;
bdclose('all');
clc;
addpath('common','common/simulator-toolbox','common/simulator-toolbox/attitude_library','common/simulator-toolbox/trajectory_library');
addpath('Functions T1');
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

% Load Excitation Signal:
load ExcitationM_eq.mat
ExcitationM = ExcitationM_eq;
r = ExcitationM(:,2);

% selected interval:
t = ExcitationM(:,1);
simulation_time = t(end) - t(1);

sys_real = ss(A,B,C(2,:),D(2));

% Load and Set Simulation:
load_system('Simulator_Single_Axis');
set_param('Simulator_Single_Axis',"FastRestart","off");

time = 0 : ctrl.sample_time : simulation_time;

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


%% MONTECARLO SIMULATIONS:

% Number of Monte Carlo Simulations:
N = 1000;

% Pre-allocate matrices for MC dataset:
Nsamples = 26647;

U  = zeros(Nsamples, N);
Q  = zeros(Nsamples, N);
Ax = zeros(Nsamples, N);


simInputs(N,1) = Simulink.SimulationInput('Simulator_Single_Axis');
Ts = ctrl.sample_time;

% Generate N different seeds for different noise realisations:
for k = 1:N

    seed_struct.x     = 1000 + 10 * k;
    seed_struct.vx    = 2000 + 10 * k;
    seed_struct.theta = 3000 + 10 * k;
    seed_struct.q     = 4000 + 10 * k;

    simIn = Simulink.SimulationInput('Simulator_Single_Axis');
    simIn = simIn.setVariable('seed', seed_struct);

    % Noise parameters:
    simIn = simIn.setVariable('noise.Enabler', 1);
    simIn = simIn.setVariable('noise.pos_stand_dev', 0.0011);
    simIn = simIn.setVariable('noise.vel_stand_dev', 0.01);
    simIn = simIn.setVariable('noise.attitude_stand_dev', deg2rad(0.33));
    simIn = simIn.setVariable('noise.ang_rate_stand_dev', deg2rad(1));

    % Other inputs:
    simIn = simIn.setVariable('ctrl', ctrl);
    simIn = simIn.setVariable('simulation_time', simulation_time);
    simIn = simIn.setVariable('ExcitationM', ExcitationM);
    simIn = simIn.setVariable('delay', delay);
    simIn = simIn.setVariable('realModel', realModel);

    simInputs(k) = simIn;
end

% Parallel simulation:
out = parsim(simInputs, 'ShowProgress', 'on', 'TransferBaseWorkspaceVariables','on');

% Extract Data:
for k = 1:N
    U(:,k)  = out(k).Mtot;
    Q(:,k)  = out(k).q;
    Ax(:,k) = out(k).ax;
end


% Save
save('MC_IO_data.mat', 'U', 'Q', 'Ax', 'Ts', '-v7.3');

vars_to_keep = {'U', 'Q', 'Ax', 'Ts', 'ctrl', 'delay', 'realModel', 'N'};
clearvars('-except', vars_to_keep{:});

%% Load Simulated I/O data for MonteCarlo Analysis:
load('MC_IO_data.mat');

%% Monte Carlo Estimates (T1 and T2):
N = 1000;
THETA_T1  = zeros(6,N);
THETA_T2  = zeros(6,N);

% Parameters and Options:
odefun = 'modelClass_SO';
guess = zeros(6,1);
guess(5) = -9;
parameters = {'Xu', guess(1); 'Xq', guess(2); 'Mu', guess(3); ...
              'Mq', guess(4); 'Xd', guess(5); 'Md', guess(6)};
sys_init = idgrey(odefun, parameters, 'c');

% Options:
opt = greyestOptions;
opt.SearchMethod = 'lsqnonlin';
opt.SearchOptions.MaxIterations = 500;
opt.SearchOptions.TolFun = 1e-6;
opt.SearchOptions.TolX = 1e-6;
opt.Regularization.Nominal = 'model';   
opt.Regularization.R = eye(6);          
opt.Regularization.Lambda = 1e-5;   
opt.EnforceStability = false;      
opt.Display = 'off';

p = 225;
f = 210;
n = 5;
om1 = 3.3;      % [rad/s]
om2 = 10.6;     % [rad/s]
N = 12;         % [-] Filter Order

% Parallelized for:
parfor i = 1 : N
    fprintf('Iteration n° %d...\n', i);
    % Extract Input-Outputs:
    u = U(:,i);
    y = Q(:,i);
        
    %%%%%% TASK 1:
    % Estimate Parameters:
    [greyest_data, ~, ~, ~, ~] = prefiltering(y, u, ctrl, r, sys_real);

    % System Identification (Grey-box)
    estModel = struct;
    est_model = greyest(greyest_data, sys_init, opt);
    
    % Outputs
    estModel.parameters = est_model.Report.Parameters.ParVector;

    % Post Processing (using sigma %):
    est_parameters_raw = estModel.parameters;
    threshold = 20; % Set to 0 all parameters with ss% above threshold
    [theta_ref, estModel, ss_pctg] = postProcessing(estModel, threshold);
    
    % Store estimated parameters:
    THETA_T1(:,i)  = theta_ref;

    %%%%%%% TASK 2:
    [G_np, om_low, om_high, f_np] = nonParamFRFest(r,u,y,Ts,sys_real); 
    
    % Remove Delay from Input-Output:
    u = u(1:end-1);
    y = y(2:end);
    
    % PBSID:
    [Ai, Bi, Ci, Di, Ki, S, X] = PBSID(u, y, p, f, n);

    sys_PBSID_DT = ss(Ai, Bi, Ci, Di, Ts);
    sys_hat_CT = d2c(sys_PBSID_DT, 'tustin');

    % Hinf structuring:
    [theta_hat, theta0] = modelMatching(sys_hat_CT, om1, om2, N);

    % Store estimated parameters:
    THETA_T2(:,i)  = theta_hat;
end
 
% save('MC_T1', 'THETA_T1');
% save('MC_T2', 'THETA_T2');

