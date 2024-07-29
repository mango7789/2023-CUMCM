clear;clc;
theta = pi / 3 * 2; % 多波束测量器开角
alpha = 1.5/180 * pi; % 坡度角
D_0 = 70; % 海域中心处海水深度
b = 200; % 每个测量点的水平间隔

x = -800:b:800;
D = D_0 - x * tan(alpha); % 海水深度
W = D * sin(theta / 2) * cos(alpha) * (1 / cos(theta / 2 + alpha) + 1 / cos(theta / 2 - alpha)); % 覆盖宽度

seq = -4:1:4;
l = seq * b - D * sin(theta / 2) / cos(theta / 2 + alpha) * cos(alpha); % 每个测量点左边界坐标
r = seq * b + D * sin(theta / 2) / cos(theta / 2 - alpha) * cos(alpha); % 每个测量点右边界坐标
r = r(1:end - 1);
l = l(2:end);
eta = (r - l) ./ W(2:end); % 与前一条测线重叠率

%% 将计算结果写入result1.xlsx
writematrix(D, "result1.xlsx", "Range", "B2:J2");
writematrix(W, "result1.xlsx", "Range", "B3:J3");
writematrix(eta, "result1.xlsx", "Range", "C4:J4");
