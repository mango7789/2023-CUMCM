clear;clc;
A=readmatrix("附件.xlsx");
A=-A(2:end,3:end);
[m,n]=size(A);
theta=pi*2/3;
eta=0.1;

%% 分割A1,A2,A3,A4
O_x=28;O_y=28; % 海水最浅处
P_1y=228;P_2x=178;P_3x=121; % 分割点,基于Excel数据的分析
A1=A(P_1y:m,1:P_2x);
A2=A(P_1y:m,P_2x:n);
A3=A(1:P_1y,1:P_3x);
A4=A(1:P_1y,P_3x:n);
% 将A3,A4边界以外的深度设置为-Inf
% for i = 1:57
%     A3(1:4*i,P_3x+i)=-Inf;
%     A4(4*i-2:P_1y,i)=-Inf;
% end

%% 处理A1
A1=mean(A1,2);num_A1=m-P_1y+1;
X=ones(num_A1,2);
Y=zeros(num_A1,1);
for i=1:num_A1
    X(i,2)=(P_1y-2+i)*0.02*1582;
    Y(i)=A1(i);
end
K1=(X'*X)\X'*Y;
num=100;
x=linspace(P_1y-1,m-1,num)'*0.02*1582;
y=[ones(num,1),x]*K1;
hold on;
plot(x,y,"r-");plot((P_1y-1:1:m-1)*0.02*1582,A1,"b*");
legend("最小二乘拟合直线","数据点");
xlabel("沿南-北方向");ylabel("海水深度");
hold off;
alpha=atan(-K1(2));gamma=alpha;
D_0=-A1(end);
L=(m-P_1y+1)*0.02*1582;
k_l=1+tan(alpha)*sin(theta/2)*cos(alpha)/(cos(theta/2+alpha));
c_l=-D_0*sin(theta/2)*cos(alpha)/(cos(theta/2+alpha));
k_r=1-tan(alpha)*sin(theta/2)*cos(alpha)/(cos(theta/2-alpha));
c_r=D_0*sin(theta/2)*cos(alpha)/(cos(theta/2-alpha));
k_w=-tan(alpha)*sin(theta/2)*cos(gamma)*(1/cos(theta/2+gamma)+1/cos(theta/2-gamma));
c_w=D_0*sin(theta/2)*cos(gamma)*(1/cos(theta/2+gamma)+1/cos(theta/2-gamma));
x_0=-c_l/k_l;
total_l=x_0; 
r=k_r*total_l+c_r;
l=k_l*total_l+c_l;
count_l1=1; % 所需测线的条数
while(r<L)
    count_l1=count_l1+1;
    total_l=(k_r*total_l+c_r-c_l-eta*c_w)/(eta*k_w+k_l); % 每次以最小重叠率选取下一条测线
    r=k_r*total_l+c_r;
    l=k_l*total_l+c_l;
end
l_A1=count_l1*P_2x*0.02*1582;
out_A1=r-L; % A1测线超出部分的长度


%% 处理A2
num_A2=(m-P_1y+1)*(n-P_2x+1);
Z=zeros(num_A2,1);
I=ones(num_A2,2);
x_A2=zeros(num_A2,1);
y_A2=zeros(num_A2,1);
for i=1:m-P_1y+1
    for j=1:n-P_2x+1
        count=(i-1)*(m-P_1y+1)+j;
        x_A2(count)=sqrt((i-(n-P_2x+1))^2+(j-(m-P_1y+1))^2)*0.02*1852;
        y_A2(count)=-A(i-1+P_1y,j-1+P_2x);
        Z(count)=y_A2(count);
        I(count,2)=x_A2(count);
    end
end 
K2=(I'*I)\I'*Z;
hold on;
x=linspace(0,1200,1000)';
plot(x,[ones(1000,1),x]*K2,"r-");
plot(x_A2,y_A2,"bx");
legend("最小二乘拟合直线","数据点");
xlabel("沿半径方向");ylabel("海水深度");
hold off;

alpha=atan(K2(2));gamma=alpha;
D_0=-A2(m-P_1y+1,1);
L=(n-P_2x+1)*0.02*1582;
k_l=1+tan(alpha)*sin(theta/2)*cos(alpha)/(cos(theta/2+alpha));
c_l=-D_0*sin(theta/2)*cos(alpha)/(cos(theta/2+alpha));
k_r=1-tan(alpha)*sin(theta/2)*cos(alpha)/(cos(theta/2-alpha));
c_r=D_0*sin(theta/2)*cos(alpha)/(cos(theta/2-alpha));
k_w=-tan(alpha)*sin(theta/2)*cos(gamma)*(1/cos(theta/2+gamma)+1/cos(theta/2-gamma));
c_w=D_0*sin(theta/2)*cos(gamma)*(1/cos(theta/2+gamma)+1/cos(theta/2-gamma));
total_l=0; 
r=k_r*total_l+c_r;
count_l2=0; % 测线的长度
while(r<L)
    count_l2=count_l2+pi/2*(D_0-total_l*tan(alpha)-K2(1))/K2(2);
    total_l=(k_r*total_l+c_r-c_l-eta*c_w)/(eta*k_w+k_l); % 每次以最小重叠率选取下一条测线
    r=k_r*total_l+c_r;
end

%% 处理A3
num_A3=P_1y*P_3x;
x_A3=zeros(num_A3,1);
y_A3=zeros(num_A3,1);
Z=zeros(num_A3,1);
I=ones(num_A3,5);
count=0;
for i=1:P_1y
    for j=1:P_3x
        count=count+1;
        x_A3(count)=sqrt((i-27)^2/4+(j-27)^2)*0.02*1852;
        y_A3(count)=log(-A(i,j));
        Z(count)=y_A3(count);
        I(count,2)=(j-1)*0.2*1852;I(count,3)=I(count,2)^2;
        I(count,4)=(i-1)*0.2*1852;I(count,5)=I(count,4)^2;
    end
end
plot(x_A3,y_A3,"bx");
K3=(I'*I)\I'*Z;

%% 处理A4
num_A4=P_1y*(n-P_3x+1);
x_A4=zeros(num_A4,1);
y_A4=zeros(num_A4,1);
Z=zeros(num_A4,1);
I=ones(num_A4,3);
count=0;
for i=1:P_1y
    for j=1:n-P_3x+1
        count=count+1;
        Z(count)=-A4(i,j);
        I(count,2)=(j-n)*0.2*1852;
        I(count,3)=(i-1)*0.2*1852;

    end
end
K4=(I'*I)\I'*Z;
alpha=atan(K4(2)^2+K4(3)^2);
D_0=A4(1,n-P_3x+1);

