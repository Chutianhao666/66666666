% Zadoff-Chu (ZC) 序列互相关测试

% 参数设置
N = 64;     % ZC 序列的长度
m1 = 7;     % ZC 序列1的常数 m1，满足 gcd(m1, N) = 1
m2 = 13;     % ZC 序列2的常数 m2，满足 gcd(m2, N) = 1
k_max = 64; % 测试的最大延迟

% 生成两个 ZC 序列
z1 = zeros(1, N);
z2 = zeros(1, N);
for n = 0:N-1
    z1(n+1) = exp(1j * pi * m1 * n * (n+1) / N);  % ZC 序列1
    z2(n+1) = exp(1j * pi * m2 * n * (n+1) / N);  % ZC 序列2
end

% 计算互相关函数
R_cross = zeros(1, k_max);
for k = 0:k_max-1
    % 计算互相关值 R_cross(k)
    R_cross(k+1) = sum(z1 .* conj(circshift(z2, [0, k])));  % 循环移位并计算互相关
end

% 绘制互相关结果
figure;
plot(abs(R_cross));
