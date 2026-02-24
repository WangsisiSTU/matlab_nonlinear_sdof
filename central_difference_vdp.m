% ---------------------------------------------------------
% 中心差分法：求解非线性 SDOF (Van der Pol 形式)
% 方程: (1 + x^2)*xdd + (c0 + x^2)*xd + x = 0
% 参数: c0=2, x0=-0.02, v0=0
% ---------------------------------------------------------
clear; clc; close all;

%% 1) 参数设置
c0 = 2;
dt = 0.01;          % 中心差分法建议取较小时间步
T_total = 100;
time = 0:dt:T_total;
N = numel(time);

x = zeros(1, N);    % 位移
v = zeros(1, N);    % 速度
a = zeros(1, N);    % 加速度

x(1) = -0.02;
v(1) = 0;

% 动力学方程提供的加速度
acc = @(disp, vel) -((c0 + disp.^2).*vel + disp) ./ (1 + disp.^2);
a(1) = acc(x(1), v(1));

%% 2) 起步：用泰勒展开得到 x(2)
% x(t+dt) = x(t) + dt*v(t) + 0.5*dt^2*a(t)
x(2) = x(1) + dt*v(1) + 0.5*dt^2*a(1);

%% 3) 中心差分主循环
% 在 t_n 处采用:
%   xdd_n ≈ (x_{n+1} - 2x_n + x_{n-1}) / dt^2
%   xd_n  ≈ (x_{n+1} - x_{n-1}) / (2dt)
% 代入方程并整理得显式更新:
%   A_n*x_{n+1} + B_n = 0  ->  x_{n+1} = -B_n / A_n
for n = 2:N-1
    m_n = 1 + x(n)^2;
    c_n = c0 + x(n)^2;

    A_n = m_n/dt^2 + c_n/(2*dt);
    B_n = -2*m_n*x(n)/dt^2 + m_n*x(n-1)/dt^2 ...
        - c_n*x(n-1)/(2*dt) + x(n);

    x(n+1) = -B_n / A_n;
end

%% 4) 后处理：由中心差分恢复 v, a
% 内点用中心差分，端点用单边差分
for n = 2:N-1
    v(n) = (x(n+1) - x(n-1)) / (2*dt);
    a(n) = (x(n+1) - 2*x(n) + x(n-1)) / dt^2;
end

v(N) = (x(N) - x(N-1)) / dt;
a(N) = acc(x(N), v(N));

%% 5) 绘图
figure('Color','w');
subplot(3,1,1);
plot(time, x, 'b-', 'LineWidth', 1.2); grid on;
ylabel('Displacement x'); xlim([0, T_total]);
title(sprintf('Central Difference Method (c0=%.1f, dt=%.4f)', c0, dt));

subplot(3,1,2);
plot(time, v, 'r-', 'LineWidth', 1.2); grid on;
ylabel('Velocity xdot'); xlim([0, T_total]);

subplot(3,1,3);
plot(time, a, 'k-', 'LineWidth', 1.2); grid on;
xlabel('Time (s)'); ylabel('Acceleration xddot'); xlim([0, T_total]);
