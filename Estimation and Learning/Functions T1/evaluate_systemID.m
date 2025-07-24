function modelEval = evaluate_systemID(estModel, realModel)
% Function to evaluate the estimated model in terms of errors on the
% estimated parameters, asymptotic variances and standard deviations

% Absolute Error:
modelEval.absErr = abs(estModel.parameters - realModel.parameters);

% Relative Error (%):
modelEval.relErr = 100 * modelEval.absErr ./ abs(realModel.parameters);

% Asymptotic Variance:
modelEval.asymptoticVar = diag(estModel.cov);

% Standard Deviations:
modelEval.stdDev = sqrt(modelEval.asymptoticVar);