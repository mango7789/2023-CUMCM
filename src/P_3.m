clear;clc;
L = 4 * 1852; % 海域的长度
S = 2 * 1852; % 海域的宽度
theta = 2/3 * pi; % 多波束换能器的开角
alpha = 1.5/180 * pi; % 坡度
D_center = 110; % 海域中心处的海水深度
D_0 = D_center + L / 2 * tan(alpha); % 最深处的海水深度
eta = 0.1; % 覆盖率

n = 5000; % 选取的点数
beta = linspace(atan(0.5), pi / 2, n); % 方向角
l = zeros(1, n); % 存储测线的总长度

for i = 1:n
    beta_now = beta(i);
    gamma = acos(1 / sqrt(1 + sin(beta_now) ^ 2 * tan(alpha) ^ 2));
    % l(x),r(x),w(x)的系数
    k_l = 1 + tan(alpha) * sin(theta / 2) * cos(gamma) / (cos(theta / 2 + gamma));
    c_l = -D_0 * sin(theta / 2) * cos(gamma) / (cos(theta / 2 + gamma));

    k_r = 1 - tan(alpha) * sin(theta / 2) * cos(gamma) / (cos(theta / 2 - gamma));
    c_r = D_0 * sin(theta / 2) * cos(gamma) / (cos(theta / 2 - gamma));

    k_w = -tan(alpha) * sin(theta / 2) * cos(gamma) * (1 / cos(theta / 2 + gamma) + 1 / cos(theta / 2 - gamma));
    c_w = D_0 * sin(theta / 2) * cos(gamma) * (1 / cos(theta / 2 + gamma) + 1 / cos(theta / 2 - gamma));

    x_0 = -c_l / k_l; % 初始最优点
    l_std = S / sin(beta_now); % 一条测线的标准长度
    flag = 0; % 用于判断测线是否越过右边界

    while (1)
        right = k_r * (x_0 - l_std * cos(beta_now)) + c_r;
        % 测线探测的范围已覆盖整个海域，退出循环
        if (right >= L)
            break
        end

        % 测线越过右边界
        if (x_0 > L)
            flag = 1;
            l(i) = l(i) + l_std - (x_0 - L) / cos(beta_now);
            % 测线越过左边界
        elseif (x_0 - l_std * cos(beta_now) < 0)
            l(i) = l(i) + x_0 / cos(beta_now);
            % 整条侧线都在边界范围内
        else
            l(i) = l(i) + l_std;
        end

        if (flag == 0)
            x_0 = (k_r * x_0 - k_l * x_0 - k_w * eta * x_0 + c_r - c_l - c_w * eta) / ((k_w * eta + k_l) * sin(beta_now) ^ 2) + x_0; % 测线未越界时的最优点
        else
            x_0 = x_0 - (k_w * L * eta + c_w * eta - k_r * L + k_l * L + c_l - c_r) / (k_r * sin(beta_now) ^ 2); % 测线越界时的最优点
        end

    end

end

semilogy(beta, l, "r-.");
hold on;
plot(beta(end), l(end), "b*");
hold off;
title("测线总长度与方向角beta的关系")
xlabel("beta"); ylabel("测线总长度/m");
legend("测线总长度--beta", "最低点");
