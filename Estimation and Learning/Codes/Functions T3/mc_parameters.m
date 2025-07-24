function [mu, sigma, ssPct] = mc_parameters(THETA, task)

task = string(task);

xlabels = {
    '$X_u \; [1/s]$', ...
    '$X_q \; [m/(s rad)]$', ...
    '$M_u \; [rad s/m]$', ...
    '$M_q \; [1/s]$', ...
    '$X_\delta \; [m/s^2]$', ...
    '$M_\delta \; [rad/s^2]$' ...
};


if strcmpi(task, ('T1'))
    mu = zeros(1,4);
    sigma = zeros(1,4);
    figure;
    for i = 1:4
        subplot(2,2,i);
        param_raw = THETA(i+2,:);
        % Discard all parameters equal to 0 (happens when ss% is above
        % threshold)
        param = param_raw(param_raw ~= 0);
        
        histogram(param, 'Normalization', 'pdf', 'FaceColor', [0.6 0.6 0.6], ...
            'EdgeColor', [0.3 0.3 0.3], 'NumBins', 50); hold on;
        
        mu(i) = mean(param);
        sigma(i) = std(param);
  
        xline(mu(i), '--b', 'LineWidth', 1.1);
    
        xlabel(xlabels{i+2}, 'Interpreter', 'latex', 'FontSize', 14);
        ylabel('Probability Density', 'FontSize', 14);
        grid on;
        box on;
    end
else
    mu = zeros(1,6);
    sigma = zeros(1,6);
    figure;
    for i = 1:6
        subplot(3,2,i);
        param = THETA(i,:);
        
        histogram(param, 'Normalization', 'pdf', 'FaceColor', [0.6 0.6 0.6], ...
            'EdgeColor', [0.3 0.3 0.3], 'NumBins', 50); hold on;
        
        mu(i) = mean(param);
        sigma(i) = std(param);

        xline(mu(i), '--b', 'LineWidth', 1.1);
    
        xlabel(xlabels{i}, 'Interpreter', 'latex', 'FontSize', 14);
        ylabel('Probability Density', 'FontSize', 14)
        grid on;
        box on;
    end
end

ssPct = 100 * abs(sigma ./ mu);

% Title:
sgtitle(sprintf('Monte Carlo Parameter Distributions (%s)', task));

end