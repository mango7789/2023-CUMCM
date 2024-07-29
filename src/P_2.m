clear;clc;
%% 已知变量
d = 0:0.3:2.1;
d = d * 1852; % 距海域中心点距离，单位转化为m
beta = 0:45:315;
beta = beta / 180 * pi; % 测线方向夹角β，单位转化为弧度制
theta = 120/180 * pi; % 多波束测量器开角
alpha = 1.5/180 * pi; % 坡度角
D_0 = 120; % 海域中心处海水深度

%% 覆盖宽度计算
W = zeros(8, 8);

for i = 1:8
    for j = 1:8
        beta_now = beta(i);
        d_now = d(j);
        gamma = acos(1 / sqrt(1 + sin(beta_now) ^ 2 * tan(alpha) ^ 2)); % 计算与测线垂直的平面的坡度角
        D = D_0 + d_now * cos(beta_now) * tan(alpha); % 沿beta方向运动d距离后的水深
        W(i, j) = D * sin(theta / 2) * cos(gamma) * (1 / cos(theta / 2 + gamma) + 1 / cos(theta / 2 - gamma)); % 覆盖宽度
    end

end

%% 将结果写入result2.xlsx文件
writematrix(W, "result2.xlsx", "Range", "C3:J10")
