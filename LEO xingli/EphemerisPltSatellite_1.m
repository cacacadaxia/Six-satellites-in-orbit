

 % =========================================================================
%
%                  卫星位置速度解算
% 
%
% =========================================================================
%
%　(C)2019-2020 广州海格通信有限公司
%   版本：V1.1
%   日期：2019年7月18日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能:  1.读取二行星历信息转换成轨道为六根数信息
%        2. 利用轨道六根数计算卫星的速度和位置
%        3. 坐标系的转换，具体看文档
%        4. 后期需要加入"地球引力摄动"的影响，需要将轨道变成不断变化的
%        5.
%        6.
%--------------------------------------------------------------------------
clear all;
close all;

% ----随便获得一组轨道六根数-------------
[oe,epoch,yr,M,E,satname] = TLE2oe('xingli.txt');    % 两行星历到六根数
a = oe(1);      % a 半长轴
e = oe(2);      % e 偏心率
i = oe(3);      % i 轨道平面倾角
Om = oe(4);     % Om 升交点赤经
om = oe(5);     % om 近升角距，近地点角距，或者叫近地点幅角
nu = oe(6);     % nu 时刻（单位是s）
% 实际上这里有个问题就是，这里并不是时刻（卫星近地点时刻）
% 应该是与下面的 sita 值是一回事


% ------参数----------------------------
aG = 40/180*pi;         % 格林尼治恒星时角
tw = 0;                 % 书里的参量
GASTw = aG;             % 格林尼治恒星时角

t0 = 0;                 % 卫星近地点时刻 (单位：s)
t  = t0;                % none，初始赋值
tspan = 1;
coutmax = 6.117e+03;
% -------------信关站-------------------------------------
Hx = 6371e3;        % 高度
Bx = 40/180*pi;     % 纬度
Lx = 60/180*pi;     % 经度

% ------------地球-----------------------
eE2    = 0.00669437999013;  % 地球离心率的平方
rE     = 6378.137*1e3;      % 地球半径
owE    = 7.29211567e-5;     % 地球自转的角速度
% ------经纬度坐标系-->信关站大地坐标系-------------------
N  = rE / (sqrt(1-eE2*sin(Bx)^2));
Xx = (N+Hx)*cos(Bx)*cos(Lx);
Yx = (N+Hx)*cos(Bx)*sin(Lx);
Zx = (N*(1-eE2)+Hx)*sin(Bx);
rX = [Xx;Yx;Zx];            % 信关站的位置
% -------------------------------------
for ii = 1:coutmax+1
    miu = 3.986004415e14;   % 地球常数
    n0 = sqrt(miu/(a)^3);   % 角速度
    t  = t + tspan;
    M = n0*(t-t0);          % 计算平近点角（应该是不断去更新的）
    
    GAST = GASTw + owE*(t-t0);
    % ---------------------------------------------
    sigema = 0.01;      % 精度
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
    % figure;plot(Etp)
    r = a*(1-e*cos(E));
    sita = 2*atan( (sqrt((1+e)/(1-e)))*tan(E/2) );  % 
    uk = sita + om;         % 书里面的做法，求取位置时暂时不考虑，相当于将om的变化已经考虑进去了
    
    x_a(ii) = r*cos(sita);  % 轨道坐标系下的三轴位置，原点是地球球心
    y_a(ii) = r*sin(sita);
    z_a(ii) = 0;
    % -------速度------------------------
    Ekdot  = n0/(1-e*sin(E));
    Faidot = sqrt((1+e)/(1-e))*(cos(sita/2)^2)/(cos(E/2)^2)*Ekdot;
    ukdot  = Faidot;
    rkdot  = Ekdot*a*e*sin(E);
    vx(ii) = rkdot*cos(uk) - ukdot*y_a(ii);
    vy(ii) = rkdot*sin(uk) + ukdot*x_a(ii);
    vz(ii) = 0;
end
% ------------------------------------
R1 = [cos(om),-sin(om),0 ; sin(om),cos(om),0 ; 0 ,0 , 1];
R2 = [1 ,0 ,0 ; 0 ,cos(i),-sin(i) ; 0 ,sin(i),cos(i)];
R3 = [cos(Om-GAST),-sin(Om-GAST),0 ; sin(Om-GAST),cos(Om-GAST),0 ; 0 ,0 ,1];    % Lk
r_bar = R1*R2*R3*[ x_a; y_a; z_a ];     % 卫星在大地坐标系，应该就是WGS-84坐标系
v_bar = R2*R3*[ vx ; vy ; vz ];         % 卫星在大地坐标系



fprintf('开始时刻：%3.0f , 结束时刻：%3.0f , \n',t0,t);
% -----------画图----------------------------
figure(1);plot3(r_bar(1,:),r_bar(2,:),r_bar(3,:),'r','LineWidth',1);
hold on;plot3(0,0,0,'*');grid on;
legend('卫星运行轨迹','地球球心');

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

%% 

% figure(4);
% plot(t0:tspan:t0+tspan*coutmax, v_bar(1,:) , 'r');  hold on;
% plot(t0:tspan:t0+tspan*coutmax, v_bar(2,:) , 'k');  hold on;
% plot(t0:tspan:t0+tspan*coutmax, v_bar(3,:) , 'b');  hold on;
% grid on;
% legend('x轴速度','y轴速度','z轴速度');
% title('大地坐标系下卫星状态变化')
% xlabel('Time s')
% ylabel('v m/s')
