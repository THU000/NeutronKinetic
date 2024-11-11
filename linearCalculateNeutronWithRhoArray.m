function neutron = linearCalculateNeutronWithRhoArray(T, dt, rhoArray, lambda, beta, age)
% 计算不同反应性情况下点堆中子动力学方程中的中子数量随时间变化
% 输入：
%   T: 总时间
%   dt: 时间步长
%   rhoArray: 不同反应性值的数组，每个元素代表一个时间步长下的反应性
%   lambda: 六组先驱核衰变常数数组
%   beta: 六组先驱核份额数组
%   age: 中子寿命
% 输出：
%   neutron: 初始 t=0 到 t=T 下的中子数数组

N = T / dt;
betaAll = sum(beta);

neutron = zeros(1, N + 1);
neutron_1 = zeros(1, N + 1);
G1 = zeros(1, 6);
G2 = zeros(1, 6);
C = zeros(6, N + 1);

neutron(1) = 1.0;
C(1, 1) = 173.3871;
C(2, 1) = 466.8852;
C(3, 1) = 114.7748;
C(4, 1) = 85.31561;
C(5, 1) = 6.561404;
C(6, 1) = 0.9069767;
for i = 1:6
    G1(i) = (1 - exp(-lambda(i) * dt)) / lambda(i);
    G2(i) = (exp(-lambda(i) * dt) * (dt * lambda(i) + 1) - 1) / lambda(i) / lambda(i);
end

for n = 2:(N + 1)
    rho = rhoArray(n)*n*dt;
    Lambda = age * (1 - rho);
    F1 = dt * rhoArray(n)*(2*n-1)*dt / 2*Lambda;
    F2 = -0.5 * dt * F1;

    m = zeros(1, 10);
    for i = 1:6
        m(1) = m(1) + C(i, n - 1);
        m(2) = m(2) + exp(-dt * lambda(i)) * C(i, n - 1);
        m(3) = m(3) + beta(i) * G2(i) / Lambda;
        m(4) = m(4) + lambda(i) * exp(-dt * lambda(i)) * C(i, n - 1);
        m(5) = m(5) + lambda(i) * beta(i) * G2(i) / Lambda;
        m(6) = m(6) + beta(i) * G1(i) / Lambda;
        m(7) = m(7) + lambda(i) * beta(i) * G1(i) / Lambda;
    end
    neutron(n) = (neutron(n - 1) + m(1) - m(2) + (F2 - m(3)) * m(4) / (1 - m(5))) /...
        (1 - F1 + m(6) - (F2 - m(3)) * ((rho - betaAll) / Lambda + m(7)) / (1 - m(5)));
    neutron_1(n) = (neutron(n) * (m(7) + (rho - betaAll) / Lambda) + m(4)) / (1 - m(5));
    for i = 1:6
        C(i, n) = C(i, n - 1) * exp(-dt * lambda(i)) + (neutron(n) * G1(i) -...
            neutron_1(n) * G2(i)) * beta(i) / Lambda;
    end
end
end