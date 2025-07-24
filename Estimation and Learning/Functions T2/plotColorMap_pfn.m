function plotColorMap_pfn()

addpath('Mats');
load PBSID_FrequencyGridSearch.mat;

% Best config
n= 5;
numConfigs = size(configs, 1);
[minERR, idx] = min(max_bodemag_errors);
best_config = configs(idx,:);
fprintf('\nBest config: p=%d, f=%d, n=%d --> minErr=%.2f\n', ...
best_config(1), best_config(2), n, minERR);

%% ColorMap:
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

for i = 1:P
    for j = 1:F
        if f_range(j) > p_range(i)
            ErrGrid(i,j) = NaN;
        end
    end
end

% Mesh grid
[F_grid, P_grid] = meshgrid(f_range, p_range);


figure();
h = pcolor(F_grid, P_grid, ErrGrid, 'HandleVisibility', 'off');
set(h, 'EdgeColor', "none");       
set(gca, 'YDir', 'normal');
colormap('parula');
colorbar;

alpha_data = ~isnan(ErrGrid);      
set(h, 'AlphaData', alpha_data);

c = colorbar;
ylabel(c, '$\max(|G_{\mathrm{err}}|)$ [dB]', 'Interpreter', 'latex', 'FontSize', 20);
xlabel('Future horizon ($f$)', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('Past horizon ($p$)', 'Interpreter', 'latex', 'FontSize', 22);
title('PBSID$_\mathrm{opt}$: Frequency-domain error (focus band)', ...
       'Interpreter', 'latex', 'FontSize', 25);

hold on;
plot(best_config(2), best_config(1), 'wo', ...
     'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k' , 'MarkerSize', 10, ...
     'DisplayName', sprintf('Optimum ($p=%d$, $f=%d$)', best_config(1), best_config(2)));
legend('show', 'Location', 'northeast', 'Interpreter', 'latex', 'fontsize', 18);
