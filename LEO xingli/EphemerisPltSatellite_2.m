

 % =========================================================================
%
%                      卫星位置速度解算
% 
%
% =========================================================================
%
%　(C)2019-2020 广州海格通信有限公司
%   版本：V1.2
%   日期：2019年7月23日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能:  1.复制的1的，考虑地球引力摄动
%        2. 
%        3. 
%        4. 
%        5.
%        6.
%--------------------------------------------------------------------------
clear all;
close all;
%% 六根数
% ----随便获得一组轨道六根数-------------
[oe,epoch,yr,M,E,satname] = TLE2oe('xingli.txt');    % 两行星历到六根数
a = oe(1);      % a 半长轴
e = oe(2);      % e 偏心率
i = oe(3);      % i 轨道平面倾角
% i = 1.0996;
% i =  0.9425;
Om0 = oe(4);     % Om 升交点赤经
om = oe(5);     % om 近升角距，近地点角距，或者叫近地点幅角
nu = oe(6);     % nu 时刻（单位是s）
% 实际上这里有个问题就是，这里并不是时刻（卫星近地点时刻）
% 应该是与下面的 sita 值是一回事

%% % ---GPS六根数--------------------------------
% a = 0.515373216438e4^2;
% e = 0.202130551916e-1;
% i = 0.93408487355;
% Om = -0.586481950333;
% om = -0.211637251436e1;

%%



%%

% ------参数----------------------------
aG = 40/180*pi;         % 格林尼治恒星时角
tw = 0;                 % 书里的参量
GASTw = aG;             % 格林尼治恒星时角
toe = 0;                % 卫星近地点时刻 (单位：s)
                        % toe 星历参考时刻(周积秒)
t  = toe;                % none，初始赋值
tspan = 1;              % 时间步进
coutmax =  floor(4.3080e+05); 
coutmax = 86400;
% -------------信关站（经纬度坐标系下）-------------------------------------
Hx = 6371e3;        % 高度
Bx = 40/180*pi;     % 纬度
Lx = 60/180*pi;     % 经度

% ------------地球-----------------------
eE2    = 0.00669437999013;  % 地球离心率的平方
eE     = 0.0818191910428;   % 离心率  
rE     = 6378.137*1e3;      % 地球半径
owE    = 7.29211567e-5;     % 地球自转的角速度
% owE = 0;
% J2     = -0.4841667798797018e-3;    % CGCS2000球谐函数模型系数(应该不准确)
J2     = 1082.63e-6;        % 论文中的模型

% ------经纬度坐标系-->信关站大地坐标系-------------------
N  = rE / (sqrt(1-eE2*sin(Bx)^2));
Xx = (N+Hx)*cos(Bx)*cos(Lx);
Yx = (N+Hx)*cos(Bx)*sin(Lx);
Zx = (N*(1-eE2)+Hx)*sin(Bx);
rX = [Xx;Yx;Zx];            % 信关站的位置
% -------------------------------------
dotOm = 0;
for ii = 1:coutmax+1
    miu = 3.986004415e14;   % 地球常数
    n0 = sqrt(miu/(a)^3);   % 角速度
    t  = t + tspan;
    tk = t - toe;
    %----------2. 计算归化时间tk----------
    if tk > 302400
        tk = tk - 604800;
    elseif tk < -302400
        tk = tk + 604800;
    end
    % 更新 M
    if ii ==1
        M = n0*tk;
    else
        M = M + dotM*tspan;          % 计算平近点角（应该是不断去更新的） 
    end

    % GAST = GASTw + (dotOm - owE)*tk - owE*toe;      % 更新格林时角
    
    % ---------------------------------------------
    sigema = 0.0001;      % 精度
    E0 = 1.2;           % 随便设置的初值，计算开普勒方程
    for inx = 1:1000
        E1 = M + e*sin(E0);
        Etp(inx) = E1;
        if abs(E1 - E0)<sigema
            E = E1;
            break;
        end
        E0 = E1;
    end
    r = a*(1-e*cos(E));
    sita = 2*atan( (sqrt((1+e)/(1-e)))*tan(E/2) );  
    uk = sita + om;         % 书里面的做法，求取位置时暂时不考虑，相当于将om的变化已经考虑进去了
    
    % 与书中的有点不一样，注意一下
    x_a = r*cos(sita);      % 轨道坐标系下的三轴位置，原点是地球球心
    y_a = r*sin(sita);
    z_a = 0;
    if ii ==1
        savedata1 = [x_a;y_a];
    elseif ii == coutmax
        savedata2 = [x_a;y_a];
    end
    
    % -------速度------------------------
    Ekdot  = n0/(1-e*sin(E));
    Faidot = sqrt((1+e)/(1-e))*(cos(sita/2)^2)/(cos(E/2)^2)*Ekdot;
    ukdot  = Faidot;
    rkdot  = Ekdot*a*e*sin(E);
    vx = rkdot*cos(uk) - ukdot*y_a;
    vy = rkdot*sin(uk) + ukdot*x_a;
    vz = 0;
    
    % ----------位置转换到大地坐标系-----------------------------------
    lamdak = Om0 + (dotOm - owE)*tk - owE*toe;
    
    
    R1 = [cos(om),-sin(om),0 ; sin(om),cos(om),0 ; 0 ,0 , 1];
    R2 = [1 ,0 ,0 ; 0 ,cos(i),-sin(i) ; 0 ,sin(i),cos(i)];
    R3 = [cos(lamdak),-sin(lamdak),0 ; sin(lamdak),cos(lamdak),0 ; 0 ,0 ,1];    % Lk
    r_bar(:,ii) = R1*R2*R3*[ x_a; y_a; z_a ];     % 卫星在大地坐标系，应该就是WGS-84坐标系
    v_bar(:,ii) = R2*R3*[ vx ; vy ; vz ];         % 卫星在大地坐标系

    % -----------------更新----------------------------------
    K = n0*J2/2*( rE/(a*(1-e^2)) )^2;
    dotOm = -3*K*cos(i);
    Om0 = Om0 + dotOm * tspan;
    dotom = 3*K*(2 - 5/2*sin(i)^2);
    om = om + dotom * tspan;
    dotM   = n0 + 3/4*n0*J2*( 3*cos(i)^2-1 )/ (1-e^2)^(3/2) * (rE/a)^2;

    lamdakSave(ii) = lamdak;
end
lastP = r_bar(:,end);
% ------------------------------------

fprintf('开始时刻：%3.0f , 结束时刻：%3.0f , \n',toe,t);
% -----------画图----------------------------
figure(1);plot3(r_bar(1,:),r_bar(2,:),r_bar(3,:),'r','LineWidth',1);
hold on;plot3(0,0,0,'*');grid on;
R_satellite = r_bar - rX;               % 卫星-->信关站 矢量


%%
hold on;
[Xe,Ye,Ze] = sphere(24);
Xe = rE*Xe;
Ye = rE*Ye;
Ze = rE*Ze;
surf(Xe,Ye,Ze)
axis equal;
% -------------北京-模拟信关站的位置---------------------------------
longitude_beijing = (39+48/60)*pi/180; 
latitude_beijing = (116+28/60)*pi/180;
% -------------坐标系转换（大地坐标系？）-------------------
x_beijing=rE*cos(longitude_beijing)*cos(latitude_beijing);
y_beijing= rE*cos(longitude_beijing)*sin(latitude_beijing);
z_beijing= rE*sin(longitude_beijing);

plot3(x_beijing,y_beijing,z_beijing,'pr');
legend('卫星运行轨迹','地球球心','地球','信关站位置');

figure;plot(toe:tspan:toe+tspan*coutmax,lamdakSave)
%% 
norm(savedata1-savedata2);      

% figure(4);
% plot(t0:tspan:t0+tspan*coutmax, v_bar(1,:) , 'r');  hold on;
% plot(t0:tspan:t0+tspan*coutmax, v_bar(2,:) , 'k');  hold on;
% plot(t0:tspan:t0+tspan*coutmax, v_bar(3,:) , 'b');  hold on;
% grid on;
% legend('x轴速度','y轴速度','z轴速度');
% title('大地坐标系下卫星状态变化')
% xlabel('Time s')
% ylabel('v m/s')
