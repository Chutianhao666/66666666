% 滞后校正器参数选择（可调）
beta = 10;           % 滞后因子，>1
tau = 1;             % 零点位置参数
Kc = 10;             % 增益调节因子

Gc_lag = Kc * tf([tau 1], [beta*tau 1]);

% 校正后开环系统
G_open_lag = Gc_lag * G_open;

% Bode图
figure;
margin(G_open_lag);
grid on;
title('加入滞后校正器后的 Bode 图');

% 新性能指标
[~, phase_margin, wcg, ~] = margin(G_open_lag);
Kv_new = dcgain(s * G_open_lag);
fprintf('滞后校正后 Kv = %.2f\n', Kv_new);
fprintf('相位裕度 = %.2f°\n', phase_margin);
fprintf('交叉频率 = %.2f rad/s\n', wcg);
G_cl_lag = feedback(G_open_lag, 1);
figure;
step(G_cl_lag);
title('加入滞后校正器后的闭环阶跃响应');
grid on;
