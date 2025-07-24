function [theta, theta0] = modelMatching(sys_target, om1, om2, N)


%%%% 0) Systune (preprocessor):
g = 9.81;

Xu_sys = realp('Xu', 0); Xu_sys.Minimum = -3; Xu_sys.Maximum = 0;
Xq_sys = realp('Xq', 0);  Xq_sys.Minimum = 0;    Xq_sys.Maximum = 3;
Mu_sys = realp('Mu', 0);    Mu_sys.Minimum = -10; Mu_sys.Maximum = 0;
Mq_sys = realp('Mq', 0);    Mq_sys.Minimum = -6;  Mq_sys.Maximum = 0;
Xd_sys = realp('Xd', -10);  Xd_sys.Minimum = -15;   Xd_sys.Maximum = 0;
Md_sys = realp('Md', 0);  Md_sys.Minimum = 300;     Md_sys.Maximum = 500;

A_th_sys = [Xu_sys Xq_sys -g; Mu_sys Mq_sys 0; 0 1 0];
B_th_sys = [Xd_sys; Md_sys; 0];
C_th_sys = [0 1 0];
D_th_sys = 0;

% Parametric system (structured, CT):
sys_th_sys = ss(A_th_sys, B_th_sys, C_th_sys, D_th_sys);
sys_th_sys.u = 'u'; sys_th_sys.y = 'y';

% Set target system (unstructured, CT):
sys_ns = sys_target;
sys_ns.u = 'u'; sys_ns.y = 'y';

% Error between transfer functions:
Gerr_sys = sys_ns - sys_th_sys;
Gerr_sys.u = 'u'; Gerr_sys.y = 'e';

% systune requirements:
Req = TuningGoal.Gain('u', 'e', 1);  
Req.Focus = [om1 om2]; % Focus matching in a limited band
Req.Stabilize = false; % Remove stability constraint

% Systune:
optTune = systuneOptions('Display', 'final');
optTune.RandomStart = 0;
[sys_tuned, ~, ~, ~] = systune(Gerr_sys, Req, optTune);

% Extract Initial Guess:
init_guess = getBlockValue(sys_tuned);

%%%% 1) Hinfstruct:

% Define initial guess from output of systune:
theta0 = [init_guess.Xu; init_guess.Xq; init_guess.Mu; init_guess.Mq; init_guess.Xd; init_guess.Md;];

% Fix Xu, Xq, Xd (this has been done via trial and error, as fixing these
% parameters to the value estimated by systune yields better results):
Xu_hinf = theta0(1);
Xq_hinf = theta0(2);
Xd_hinf = theta0(5);
% Define variable parameters:
Mu_hinf = realp('Mu', theta0(3));   Mu_hinf.Minimum = -10; Mu_hinf.Maximum = 0;
Mq_hinf = realp('Mq', theta0(4));   Mq_hinf.Minimum = -6;  Mq_hinf.Maximum = 0;
Md_hinf = realp('Md', theta0(6));   Md_hinf.Minimum = 300;     Md_hinf.Maximum = 500;


% Parametric model:
A_th_hinf = [Xu_hinf Xq_hinf -g; Mu_hinf Mq_hinf 0; 0 1 0];
B_th_hinf = [Xd_hinf; Md_hinf; 0];
C_th_hinf = [0 1 0];
D_th_hinf = 0;

sys_theta = ss(A_th_hinf, B_th_hinf, C_th_hinf, D_th_hinf);
sys_theta.u = 'u'; sys_theta.y = 'y';

% Shift systems to ensure stability:
mu = 5;
sys_target_shift = ss(sys_ns.A - mu*eye(size(sys_ns.A, 1)), sys_ns.B, sys_ns.C, sys_ns.D);
sys_theta_shift  = ss(sys_theta.A - mu*eye(3), sys_theta.B, sys_theta.C, sys_theta.D);

% Define Butterworth Filter (pass-band):
[b,a] = butter(N, [om1, om2], 'bandpass', 's');  
GW_c = tf(b,a);
GW = ss(GW_c);

% H-infinity error (with filter):
Gerr_shift = GW*(sys_target_shift - sys_theta_shift);
Gerr_shift.u = 'u'; Gerr_shift.y = 'e';

% hinfstruct:
opts = hinfstructOptions('Display','final','RandomStart',0);
theta_opt = hinfstruct(Gerr_shift, opts);

theta = [theta0(1);
         theta0(2);
         theta_opt.Blocks.Mu.Value;
         theta_opt.Blocks.Mq.Value;
         theta0(5);
         theta_opt.Blocks.Md.Value];
end
