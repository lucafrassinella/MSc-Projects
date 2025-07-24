function [frd_data, G_np, om_low, om_high, fvec] = prefiltering(y, u, ctrl, r, sys_real)
% Prefiltering function that takes as input the initial data and filters around the frequency
% at which the PSD has its peak (-f_low, +f_high):
% Inputs:
%   - u: input data
%   - y: output data
% Output:
%   - frd_data: idfrd input-output data for greyest

f_low = 0.5; 
f_high = 9;

% Remove mean from input and output:
u = u - mean(u);
y = y - mean(y);

% Compute PSD estimate with Welch Estimator:
Ts = ctrl.sample_time;
fs = 1/Ts;

[Pxq, f] = pwelch(y, hanning(1024), 512, 2048, fs);
[~, idx_peak] = max(Pxq);
f_peak = f(idx_peak);

% Filter Design (pass-band filter around f_peak)
f_cut = [max(f_peak - f_low, 0.5), f_peak + f_high];   
[b, a] = butter(4, f_cut/(fs/2));  
q_filt = filtfilt(b, a, y);                 
u_filt = filtfilt(b, a, u);

% FRF estimate (input-output)
[G_np, om_low, om_high, fvec] = nonParamFRFest(r,u_filt,q_filt,Ts,sys_real);

omvec = 2*pi*fvec;

idx_valid = (omvec >= om_low & omvec <= (om_high + 5) & omvec >= 2*pi*f_cut(1) & omvec <= 2*pi*f_cut(2));

Gnp_valid = G_np(idx_valid);
om_valid = omvec(idx_valid);

frd_data = idfrd(Gnp_valid, om_valid, Ts);
