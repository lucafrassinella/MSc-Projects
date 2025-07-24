function [est_parameters, estModel, ss_pctg] = postProcessing(estModel, threshold)
% Post Processing function to fix to zero the estimated parameters having 
% sigma% over the threshold value (20%).

ss_pctg = zeros(6,1);
for i = 1 : length(estModel.parameters)
    ss_pctg(i) = abs(sqrt(estModel.cov(i,i)) / estModel.parameters(i)) * 100;
    if ss_pctg(i) > threshold
        estModel.parameters(i) = 0;
    end
end

% Post-Processed Estimated Parameters and Matrices:
est_parameters = estModel.parameters;
estModel.A_est = [est_parameters(1), est_parameters(2), -9.81;
                  est_parameters(3), est_parameters(4), 0;
                  0, 1, 0];
estModel.B_est = [est_parameters(5); est_parameters(6); 0];
estModel.C_est = [1 0 0; 0 1 0; 0 0 1; est_parameters(1) est_parameters(2) 0];
estModel.D_est = [0;0;0;est_parameters(5)];