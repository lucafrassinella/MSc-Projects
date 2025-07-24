function [FIT, PEC, VAF] = validation(val, t_vec, N)


% FIT:
FIT.q = max((1 - norm(val.q_m - val.q_est)^2 / norm(val.q_m - mean(val.q_m))^2), 0) * 100;
FIT.ax = max((1 - norm(val.ax_m - val.ax_est)^2 / norm(val.ax_m - mean(val.ax_m))^2), 0) * 100;
FIT.theta = max((1 - norm(val.theta_m - val.theta_est)^2 / norm(val.theta_m - mean(val.theta_m))^2), 0) * 100;

% PEC:
PEC.q = 1/sqrt(N) * norm(val.q_m - val.q_est)^2;
PEC.ax = 1/sqrt(N) * norm(val.ax_m - val.ax_est)^2;
PEC.theta = 1/sqrt(N) * norm(val.theta_m - val.theta_est)^2;

% VAF (Variance Accounted For)
VAF.q = 100 * max(1 - var(val.q_m - val.q_est)/var(val.q_m), 0);
VAF.ax = 100 * max(1 - var(val.ax_m - val.ax_est)/var(val.ax_m), 0);
VAF.theta = 100 * max(1 - var(val.theta_m - val.theta_est)/var(val.theta_m), 0);

figure();

% Subplot 1: pitch rate q
subplot(2,2,1);
plot(t_vec, val.q_m, '--r', 'LineWidth', 1)
hold on
plot(t_vec, val.q_est, 'k', 'LineWidth', 1.2)
legend('Real System', 'Estimated System')
title(['$VAF_{q} = ', num2str(VAF.q, '%.2f'), '\%$'], ...
      'Interpreter','latex', 'FontWeight', 'bold');
xlabel('$t$ [s]')
ylabel('$q$ [rad/s]')
grid on;
xlim([t_vec(1) t_vec(end)])
ylim([-1.2*max(abs(val.q_est)) +1.2*max(abs(val.q_est))])

subplot(2,2,2);
plot(t_vec, val.ax_m, '--r', 'LineWidth', 1)
hold on
plot(t_vec, val.ax_est, 'k', 'LineWidth', 1.2)
title(['$VAF_{a_x} = ', num2str(VAF.ax, '%.2f'), '\%$'], ...
      'Interpreter','latex', 'FontWeight', 'bold');
xlabel('$t$ [s]')
ylabel('$a_x$ [m/s$^2$]')
grid on;
xlim([t_vec(1) t_vec(end)])
ylim([-1.2*max(abs(val.ax_m)) +1.2*max(abs(val.ax_m))])

subplot(2,2,3);
plot(t_vec, val.theta_m, '--r', 'LineWidth', 1)
hold on
plot(t_vec, val.theta_est, 'k', 'LineWidth', 1.2)
title(['$VAF_{\theta} = ', num2str(VAF.theta, '%.2f'), '\%$'], ...
      'Interpreter','latex', 'FontWeight', 'bold');
xlabel('$t$ [s]')
ylabel('$\theta$ [rad]')
grid on;
xlim([t_vec(1) t_vec(end)])
ylim([-1.2*max(abs(val.theta_est)) +1.2*max(abs(val.theta_est))])

subplot(2,2,4);
axis off