%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   正弦反应性，不改变反应性大小，支持代码中改变
%   配合函数文件 sineCalculateNeutronWithRhoArray.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%常量 常量数组 定义
T = 800;   %总时间
dt = 0.001; %时间步长
N = T/dt;  %总步数
time = 0;  %时间积累
age = 1e-4;%中子寿命
lambda = [0.0124,0.0305,0.111,0.301,1.14,3.01]; %六组先驱核衰变常数
beta = [2.15e-4,1.424e-3,1.274e-3,2.568e-3,7.48e-4,2.73e-4]; %六组先驱核份额
betaAll = sum(beta); %总 beta
neutron_change = zeros(4,N+1);
rho_0=[0.03,8e-4,9e-4,1e-3];

%rho_0 = 1e-3;%设定反应性大小
rho = rho_0(1)*ones(1,N+1);      %阶跃正反应性，此后不随时间变化
Lambda = age*(1-rho(1));   %中子代时间
F1 = dt*rho(1)/Lambda;         %阶跃反应性不随时间变的时候，中子数计算式系数
F2 = -0.5*dt*F1; %阶跃反应性不随时间变的时候，中子数计算式系数
neutron_change(1,:)=sineCalculateNeutronWithRhoArray(T,dt,rho,lambda,beta,age);


figure; % 创建新的图形窗口

% 绘制四条中子数量变化曲线
hold on;
t = 0:dt:T;
plot(t, neutron_change(1, :), 'b');


% 设置坐标轴标签和标题
xlabel('时间（s）');
ylabel('中子密度');
title('正弦变化正反应性下中子密度随时间变化曲线');

% 设置横轴范围（可根据需求调整）
xlim([0, T]); 

% 设置纵轴范围（可根据实际数据情况调整）
ylim([0.8, max(neutron_change(1,:))]); 

% 添加图例
legend('rho=0.03sin(t*pi/5)','Location', 'best'); % 将图例位置设置为最佳位置，避免遮挡曲线

% 显示网格
grid on;