function [MC_mean, MC_std] = AsVar_PBSID()

% Load data: 
load('MeanCov_PBSID_Hinf_mat.mat');

TH_all = TH_hat';  
TH_clean = TH_all;

% Initialize the outlier removal process
iteration = 0;
while true
    is_out = isoutlier(TH_clean, 'quartiles');
    rows_to_remove = any(is_out, 2);
    if ~any(rows_to_remove)
        break;
    end
    TH_clean = TH_clean(~rows_to_remove, :);
    iteration = iteration + 1;
end

% MC statistics
MC_mean = mean(TH_clean)';
MC_std  = std(TH_clean)';
ss_pctg = abs(MC_std ./ MC_mean) * 100;

% Boxplot:
param_labels = {'$X_u$ [1/s]', '$X_q$ [1/s]', '$M_u$ [1/s$^2$]', ...
                '$M_q$ [1/s$^2$]', '$X_{\delta}$ [1/s]', '$M_{\delta}$ [1/s$^2$]'};
xlabel_text = {'$20000$ tests'}; 

figure();
for i = 1:6
    ax = subplot(3,2,i);
    
    boxplot(TH_all(:,i), ...
            'Symbol','r+', ...
            'Colors','k', ...
            'Widths', 0.5, ...
            'Whisker', 1.5);

    set(findobj(ax, 'Tag', 'Box'), 'Color', [0.1 0.4 0.6], 'LineWidth', 1.5);
    set(findobj(ax, 'Tag', 'Whisker'), 'Color', 'k', 'LineWidth', 1.2);
    set(findobj(ax, 'Tag', 'Median'), 'Color', [0.8 0.6 0.5], 'LineWidth', 1.5);
    set(findobj(ax, 'Tag', 'Outliers'), 'MarkerEdgeColor', 'r');

    xticks(1)
    xticklabels(xlabel_text)
    set(gca, ...
        'TickLabelInterpreter', 'latex', ...
        'FontSize', 18, ...
        'XTickLabelRotation', 0);

    ylabel(param_labels{i}, ...
        'Interpreter','latex', ...
        'FontSize', 18, ...
        'FontWeight','bold');
    
    pbaspect([1 1 1])
    if mod(i,2)==1
        ax.Position(1) = ax.Position(1) + 0.09;  
    else
        ax.Position(1) = ax.Position(1) - 0.09;  
    end
    grid on
end

end

