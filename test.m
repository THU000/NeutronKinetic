%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   缓冲文件，用于测试和DeBug
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%常量 常量数组 定义
T = 100;   %总时间
dt = 0.001; %时间步长
N = T/dt;  %总步数
time = 0;  %时间积累
age = 1e-4;%中子寿命
lambda = [0.0124,0.0305,0.111,0.301,1.14,3.01]; %六组先驱核衰变常数
beta = [2.15e-4,1.424e-3,1.274e-3,2.568e-3,7.48e-4,2.73e-4]; %六组先驱核份额
betaAll = sum(beta); %总 beta
neutron_change = zeros(4,N+1);
rho_0=[1e-3,2e-3,3e-3,4e-3];
for i=1:4
    %rho_0 = 1e-3;%设定反应性大小
    rho = rho_0(i)*ones(1,N+1);      %阶跃正反应性，此后不随时间变化
    Lambda = age*(1-rho(1));   %中子代时间
    F1 = dt*rho(1)/Lambda;         %阶跃反应性不随时间变的时候，中子数计算式系数
    F2 = -0.5*dt*F1; %阶跃反应性不随时间变的时候，中子数计算式系数
    neutron_change(i,:)=calculateNeutronWithRhoArray(T,dt,rho,lambda,beta,age);
end

figure; % 创建新的图形窗口

% 绘制四条中子数量变化曲线，使用对数纵坐标
semilogy(t, neutron_change(1, :), 'LineWidth', 1.5, 'DisplayName',['rho = ', num2str(rho_0(1))]);
hold on;
for i = 2:4
    semilogy(t, neutron_change(i, :), 'LineWidth', 1.5, 'DisplayName',['rho = ', num2str(rho_0(i))]);
end

% 设置坐标轴标签和标题
xlabel('时间（s）');
ylabel('中子密度（对数坐标）');
title('不同阶跃正反应性下中子密度随时间变化曲线');
% 设置横轴范围（可根据需求调整）
xlim([0, T]); 
% 设置纵轴范围（可根据实际数据情况调整，这里先找到所有数据中的最大值和最小值来设置范围）
%ylim([min(min(neutron_change)), max(max(neutron_change))]); 
ylim([min(min(neutron_change)), 1e3]); 

% 添加图例
legend('Location', 'best'); % 将图例位置设置为最佳位置，避免遮挡曲线
% 显示网格
grid on;