function [mu_re_p, sigma_re_p, mu_im_p, sigma_im_p, mu_z, sigma_z] = mc_pz(THETA, task)

N = size(THETA, 2);
g = 9.81;
task = string(task);

% Preallocate:
all_poles_real = NaN(3, N);
all_poles_imag = NaN(3, N);
all_zeros_real = NaN(2, N);

% Extract from montecarlo:
for i = 1:N
    theta = THETA(:,i);
    if any(theta(3:end) == 0)
        continue
    end
    [A, B, C, D] = theta2abcd(theta, g);
    sys = ss(A,B,C,D);
    p = pole(sys);
    z = tzero(sys);
    [~, idx] = sort(real(p), 'descend');  
    p_sorted = p(idx);
    all_poles_real(1:length(p_sorted), i) = real(p_sorted);
    all_poles_imag(1:length(p_sorted), i) = imag(p_sorted);
    all_zeros_real(1:length(z), i) = real(z);
end


mu_re_p = zeros(3,1);
mu_im_p = zeros(3,1);
sigma_re_p = zeros(3,1);
sigma_im_p = zeros(3,1);
figure();
for i = 1:3
    subplot(3,2,2*(i-1)+1);
    data = all_poles_real(i,:);
    histogram(data, 'Normalization', 'pdf', ...
        'FaceColor', [0.6 0.6 0.6], 'EdgeColor', [0.3 0.3 0.3], 'NumBins', 40); hold on;
    mu_re_p(i) = mean(data, 'omitnan'); sigma_re_p(i) = std(data, 'omitnan');
    xline(mu_re_p(i), '--b', 'LineWidth', 1.1);
    xlabel(['Re($p_' num2str(i) '$)'], 'FontSize', 18);
    ylabel('Probability Density', 'FontSize', 14);
    grid on; box on;

    subplot(3,2,2*i);
    data = all_poles_imag(i,:);
    histogram(data, 'Normalization', 'pdf', ...
        'FaceColor', [0.6 0.6 0.6], 'EdgeColor', [0.3 0.3 0.3], 'NumBins', 40); hold on;
    mu_im_p(i) = mean(data, 'omitnan'); sigma_im_p(i) = std(data, 'omitnan');
    xline(mu_im_p(i), '--b', 'LineWidth', 1.1);
     xlabel(['Im($p_' num2str(i) '$)'], 'FontSize', 18);
    ylabel('Probability Density', 'FontSize', 14);
    grid on; box on;
end

sgtitle(sprintf('Monte Carlo Pole Distributions (%s)', task));

% Zeros:
mu_z = zeros(2,1);
sigma_z = zeros(2,1);

figure();
for i = 1:2
    subplot(1,2,i);
    data = all_zeros_real(i,:);
    if i == 1
        histogram(data, 'Normalization', 'pdf', ...
            'FaceColor', [0.6 0.6 0.6], 'EdgeColor', [0.3 0.3 0.3], 'BinWidth', 0.01); hold on;
        xlim([-0.2 0.2])
    else
        histogram(data, 'Normalization', 'pdf', ...
            'FaceColor', [0.6 0.6 0.6], 'EdgeColor', [0.3 0.3 0.3], 'NumBins', 40); hold on;
    end
    mu_z(i) = mean(data, 'omitnan'); sigma_z(i) = std(data, 'omitnan');
    xline(mu_z(i), '--b', 'LineWidth', 1.1);
    xlabel(['Re($z_' num2str(i) '$)'], 'FontSize', 18);
    ylabel('Probability Density', 'FontSize', 14);
    grid on; box on;
end

sgtitle(sprintf('Monte Carlo Zero Distributions (%s)', task));

end