% 开环传递函数
G_open = tf(1, [1, 0.716, 0]);

% 绘制Bode图
bode(G_open);
grid on;
margin(G_open);