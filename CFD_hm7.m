%%  参数设置
clear;
gamma = 1.4;       % 绝热指数
x_left = 0;     % 计算域左边界
x_right = 1;     % 计算域右边界
t_end = 0.2;      % 计算终止时间
Mx = 200;           % 网格数
delta_x = (x_right - x_left)/Mx; % 空间步长
Nt = 200;
delta_t = t_end/Nt;
U = zeros(3, Mx + 4); % 3行（ρ, ρu, E），列对应单元
F = zeros(3, Mx + 4);

% 填充实际单元初始值
for i = 1:Mx+4
    [rho,u,p] = InitialConditions3(x_left+(i-3)*delta_x);
    U(1,i) = rho;
    U(2,i) = rho*u;
    U(3,i) = rho*(p/(rho*(gamma-1))+0.5*u^2);
end

%% Roe格式
s = 0;   
for t = 1:Nt
    % 计算Roe通量
    F = zeros(3, Mx+4);
    for i = 3:Mx+3
        % 左右状态
        UL = U(:,i);
        UR = U(:,i+1);
        
        % 计算Roe平均
        rhoL = UL(1);
        uL = UL(2)/UL(1);
        pL = UL(1)*(gamma-1)*(UL(3)/UL(1)-0.5*(UL(2)/UL(1))^2);
        rhoR = UR(1);
        uR = UR(2)/UR(1);
        pR = UR(1)*(gamma-1)*(UR(3)/UR(1)-0.5*(UR(2)/UR(1))^2);
        %HL = gamma*pL/((gamma-1)*rhoL);
        %HR = gamma*pR/((gamma-1)*rhoR);
        HL = (UL(3)+pL)/rhoL;
        HR = (UR(3)+pR)/rhoR;

        sqrt_rhoL = sqrt(rhoL);
        sqrt_rhoR = sqrt(rhoR);
        rho_roe = sqrt_rhoL*sqrt_rhoR;
        u_roe = (sqrt_rhoL*uL + sqrt_rhoR*uR)/(sqrt_rhoL + sqrt_rhoR);
        H_roe = (sqrt_rhoL*HL + sqrt_rhoR*HR)/(sqrt_rhoL + sqrt_rhoR);
        a_roe = sqrt((gamma-1)*(H_roe - 0.5*u_roe^2));
        
        %计算通量
        lamda = [u_roe-a_roe, u_roe, u_roe+a_roe];
        alpha = [0.5/a_roe^2*(pR-pL-rho_roe*a_roe*(uR-uL)), rhoR-rhoL-(pR-pL)/a_roe^2, 0.5/a_roe^2*(pR-pL+rho_roe*a_roe*(uR-uL))];
        r = [1 1 1;
             u_roe-a_roe u_roe u_roe+a_roe;
             H_roe-u_roe*a_roe 0.5*u_roe^2 H_roe+u_roe*a_roe];


        % 熵修正
        epsilon = 0.1;
        for k = 1:3
            if abs(lamda(k)) < epsilon
                lamda(k) = (lamda(k)^2 + epsilon^2)/(2*epsilon);
            end
        end

        FL = [rhoL*uL; rhoL*uL^2+pL; (UL(3)+pL)*uL];
        FR = [rhoR*uR; rhoR*uR^2+pR; (UR(3)+pR)*uR];
        F(:,i) = 0.5*(FL + FR) - 0.5*(abs(lamda(1))*alpha(1)*r(:,1)+abs(lamda(2))*alpha(2)*r(:,2)+abs(lamda(3))*alpha(3)*r(:,3));
    end
    
    % 更新守恒变量
    for i = 4:Mx+2
        U(:,i) = U(:,i) - delta_t/delta_x*(F(:,i) - F(:,i-1));
    end

    s = s+1;
    disp(s);
end

%% 作图
x = x_left:delta_x:x_right;
UU = zeros(3,Mx+1);
for i = 3:Mx+3
    UU(:,i-2) = U(:,i);
end
for i = 1:3
    plot(x,UU(i,:),'LineWidth',1.5);
    hold on;
end
legend('\rho','\rhou','\rhoE');
xlabel('x');
ylabel('U的各分量');
title(['Roe格式(t = ',num2str(t_end),'s)']);