%% Asymptotic Variance of PBSID Estimates:
% Note: to run this script the PBSID Toolbox from TU Delft is needed
% Link: https://www.dcsc.tudelft.nl/~jwvanwingerden/pbsid/pbsidtoolbox_product_page.html 

clearvars; close; clc;
load('Mats\varianceAnalysis.mat');

[S,X,~,U,Zp] = dordvarx(u,y,f,p);
x = dmodx(X,n);
[Ai,Bi,Ci,Di,Ki] = dx2abcdk(x,u,y,f,p);
[P,sigma,dA,dB,dC,dD,dK] = dvar4abcdk(x,u,y,f,p,Ai,Bi,Ci,Di,Ki,U,Zp);


% Save:
% save('AsVar_PBSID', 'P', 'sigma', 'dA', 'dB', 'dC', 'dD', 'dK', 'Ai', 'Bi', 'Ci', 'Di', 'Ki', 'U', 'Zp');


%% MC simulation (20000 Samples)

load('AsVar_PBSID.mat');
g = 9.81;
Ts = 1/250;
om1 = 3.3;
om2 = 10.6;
N = 12;
w = logspace(log10(0.5), log10(200), 500);

% Nominal theta_hat:
sys_nom_CT = d2c(ss(Ai,Bi,Ci,Di,Ts), 'tustin');
[~, ~, ~, ~, thetaNom_hat, ~] = structuring_hinfstruct_systune_pre(sys_nom_CT, om1, om2, N);

nA = numel(Ai);
nB = numel(Bi);
nC = numel(Ci);
nD = numel(Di);
nK = numel(Ki);

% Stack in a single array:
theta_nom = [Ai(:); Bi(:); Ci(:); Di(:); Ki(:)];

% MonteCarlo
M = 20000;
theta_samples = mvnrnd(theta_nom', P, M); 

TH_hat = zeros(6, M);
TH_hat(:,end) = thetaNom_hat;

mag_all = zeros(M, length(w));
phase_all = zeros(M, length(w));

% Nominal Bode:
[magNom,phaseNom] = bode(sys_nom_CT,w);
mag_all(end,:) = 20*log10(squeeze(magNom));
phase_all(end,:) = squeeze(phaseNom);


parfor ii = 1 : M
    fprintf('Iteration: %d \n', ii)
    % Extract random unstructured matrices:
    A_ii = reshape(theta_samples(ii, 1:nA), [5, 5]);
    B_ii = theta_samples(ii, nA+1 : nA+nB);
    B_ii = B_ii';
    C_ii = theta_samples(ii, nA+nB+1: nA+nB+nC);
    D_ii = theta_samples(ii, nA+nB+nC+nD);
    % From discrete to continous time:
    sys_CT_ii = d2c(ss(A_ii,B_ii,C_ii,D_ii,Ts), 'tustin');
    % Estimate physical parameters (structuring method):
    [thetaHat_ii, ~] = modelMatching(sys_CT_ii, om1, om2, N);
    % Store it in TH_hat:
    TH_hat(:,ii) = thetaHat_ii;
    % % Estimated physical system:
    [A,B,C,D] = theta2abcd(thetaHat_ii, g);
    sys_est_ii = ss(A,B,C,D);
    [mag,phase] = bode(sys_est_ii,w);
    mag_all(ii,:) = 20*log10(squeeze(mag));
    phase_all(ii,:) = squeeze(phase);
end


M = 20000;
g = 9.81;
Ts = 1/250;
om1 = 3.3;
om2 = 10.6;
N = 12;
w = logspace(log10(0.5), log10(200), 500);
mag_all = zeros(M, length(w));
phase_all = zeros(M, length(w));

% MC Bodes:
parfor ii = 1 : M
    thetaHat_ii = TH_hat(:,ii);
    % % Estimated physical system:
    [A,B,C,D] = theta2abcd(thetaHat_ii, g);
    sys_est_ii = ss(A,B,C,D);
    [mag,phase] = bode(sys_est_ii,w);
    mag_all(ii,:) = 20*log10(squeeze(mag));
    phase_all(ii,:) = squeeze(phase);
end

% Save results:
% save('MeanCov_PBSID_Hinf_mat', 'TH_hat', 'mag_all', 'phase_all', 'w')
