%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   阶跃正反应性，并验证理论计算曲线
%   配合函数文件 calculateNeutronWithRhoArray.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%常量 常量数组 定义
T = 100;   %总时间
dt = 0.001; %时间步长
N = T/dt;  %总步数
time = 0;  %时间积累
age = 1e-4;%中子寿命
rho_0 = 1e-3;%设定反应性大小
rho = rho_0*ones(1,N+1);      %阶跃正反应性，此后不随时间变化
Lambda = age*(1-rho(1));   %中子代时间
F1 = dt*rho(1)/Lambda;         %阶跃反应性不随时间变的时候，中子数计算式系数
F2 = -0.5*dt*F1; %阶跃反应性不随时间变的时候，中子数计算式系数
lambda = [0.0124,0.0305,0.111,0.301,1.14,3.01]; %六组先驱核衰变常数
beta = [2.15e-4,1.424e-3,1.274e-3,2.568e-3,7.48e-4,2.73e-4]; %六组先驱核份额
betaAll = sum(beta); %总 beta

%调用函数，推进时间计算
neutron_1 = calculateNeutronWithRhoArray(T,dt,rho,lambda,beta,age);
%验证程序
fprintf("100 秒后中子密度 neutron(N+1)为：%d\n",neutron(N+1));
% 计算理论曲线数据
t = 0:dt:T;
exact_solution = 1.0 * (1.446 * exp(0.0182 * t) - 0.0359 * exp(-0.0136 * t)...
    - 0.140 * exp(-0.0598 * t) - 0.0637 * exp(-0.183 * t) - 0.0205 * ...
    exp(-1.005 * t) - 0.00767 * exp(-2.875 * t) - 0.179 * exp(-55.6 * t));

% 计算误差
error = neutron_1 - exact_solution;
% 绘制理论、数值计算曲线和误差曲线在同一张图上
figure; % 创建新的图形窗口
yyaxis left; % 使用左侧纵轴绘制数值计算和理论曲线
% 'b' 表示蓝色线条绘制数值计算曲线，'r' 表示红色线条绘制理论曲线
plot(t, neutron_1, 'b', t, exact_solution, 'r'); 
xlabel('时间（s）');
ylabel('中子密度');
title('理论与数值计算中子密度随时间变化曲线');
xlim([0, T]);
yyaxis right; % 使用右侧纵轴绘制误差曲线
plot(t, error, 'g'); % 'g' 表示绿色线条绘制误差曲线
ylabel('误差');
grid on;
legend('数值计算', '理论计算','计算减理论','Location','northwest'); % 添加图例