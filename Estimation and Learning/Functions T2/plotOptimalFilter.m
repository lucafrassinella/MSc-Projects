function plotOptimalFilter()

load FilterOptimizationResults.mat;

numConfigs = size(configs, 1);
X = length(om1_range); Y = length(om2_range);
VAFmap = nan(Y, X); 

for k = 1:numConfigs
    om1 = configs(k,1);
    om2 = configs(k,2);

    idx_x = find(om1_range == om1);
    idx_y = find(om2_range == om2);
    VAFmap(idx_y, idx_x) = VAF_results(k);
end

[OM1_grid, OM2_grid] = meshgrid(om1_range, om2_range);
subMatrix = VAFmap(10:20, 10:20);

[~, linear_idx] = max(subMatrix(:));
[row_sub, col_sub] = ind2sub(size(subMatrix), linear_idx);

row_idx = row_sub + 9;
col_idx = col_sub + 9;  

idx_om1_opt = row_idx;
idx_om2_opt = col_idx;

om1_opt = om1_range(idx_om1_opt);
om2_opt = om2_range(idx_om2_opt);
VAF_opt = VAFmap(idx_om1_opt, idx_om2_opt);


figure();
surf(OM1_grid, OM2_grid, VAFmap, 'EdgeColor', [0.4 0.4 0.4], 'EdgeLighting', 'gouraud', 'HandleVisibility', 'off');
colormap('parula');
xlabel('$\omega_1$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\omega_2$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 18);
zlabel('VAF [\%]', 'Interpreter', 'latex', 'FontSize', 18);
view(45, 30);
c = colorbar;
ylabel(c, 'VAF [\%]', 'Interpreter', 'latex', 'FontSize', 18);
hold on;
plot3(om1_opt, om2_opt, VAF_opt, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
h_fake = plot3(nan, nan, nan, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
legend(h_fake, ...
    {sprintf('Optimum: $\\omega_1$ = %.1f rad/s, $\\omega_2$ = %.1f rad/s', om1_opt, om2_opt)}, ...
    'Interpreter', 'latex', 'FontSize', 18, 'Location', 'best');
title('$N=12$', 'Interpreter', 'latex', 'FontSize', 18)
