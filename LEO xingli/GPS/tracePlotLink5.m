
% =========================================================================
%
%                      导航卫星位置速度解算
%
%
% =========================================================================
%
%　(C)2019-2020 广州海格通信有限公司
%   版本：V1.5
%   日期：2019年7月25日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能:  1.复制的Link4
%        2.利用其他的一些卫星参数进行比较，尤其是低轨的卫星轨道，发现与之前有很大的出入
%        3.多普勒频偏
%        4.经过验证，模型无误
%        5.
%        6.
%--------------------------------------------------------------------------

clear all;  close all;
%---------------------读取RINEX格式n文件--------------------------------
data = RinexNreader('BR.16n','G01'); % 注意读取路径和卫星编号
%---------------------计算测量日的周积秒---------------------------------
% [JD,FOD,GPSW,SOW,DOY,DOW] = GCtoGPS(data(1,1),data(2,1),data(3,1),data(4,1),data(5,1),data(6,1));
%----------计算n文件中参考历元的周积秒----------
for j = 1:size(data,2) % size(data，2) - 返回data矩阵的列数;
    [JD,FOD,GPSW,SOW(j),DOY,DOW] = GCtoGPS(data(1,j),data(2,j),data(3,j),data(4,j),data(5,j),data(6,j));
    % 将指定卫星的 年、月、日、时、分、秒 → 儒略日(整数部分)、儒略日(小数部分)、GPS周、周积秒、年积日、星期数
end

% 地球参数
GM     = 3986005e8; % μ=GM
we     = 7.2921151467e-5; % 地球自转角速度
eE2    = 0.00669437999013;  % 地球离心率的平方
eE     = 0.0818191910428;   % 离心率
rE     = 6378.137*1e3;      % 地球半长轴
J2     = 1082.63e-6;        % 论文中的模型

% -------------北京-模拟信关站的位置 (东经120，北纬37)---------------------------------
B_bj = (39+48/60)*pi/180;       % 纬度
L_bj = (0 +28/60 - 180)*pi/180;      % 经度
hX = 0;                         % 地面高度
N  = rE / sqrt( 1 - eE2*sin( B_bj )^2 );
% -------------坐标系转换（大地坐标系？）-------------------
x_beijing = (N+hX)*cos(B_bj)*cos(L_bj);
y_beijing = (N+hX)*cos(B_bj)*sin(L_bj);
z_beijing = (N*(1-eE2)+hX)*sin(B_bj);
rXGZ = [x_beijing;y_beijing;z_beijing];     % 从地心到信关站的矢量
%-----------------------不是很重要的参数------------------------------------
inx = 1;
tspan = 1;
dotM2 = 0;
dotOm = 0;
delta_f = [];

%% 
index   = 1;
% delta_n = data(12,index);       % detn 平近点角的长期变化(近地点参数)
% sqrta   = data(17,index);       % sqrta 长半轴的平方根
toe       = data(18,index);       % t0e 星历参考时刻(周积秒)
% M0      = data(13,index);       % M0 参考时刻的平近点角
% e       = data(15,index);       % e 偏心率
% om      = data(24,index);       % omega 近地点角距
% Om0     = data(20,index);       % OMEGA0 非参考时刻toe升交点赤经，是始于格林尼治子午圈到卫星轨道升交点所谓准经度
% % dotOm   = data(25,index);       % OMEGADOT 升交点赤经在赤道平面中的长期变化
% i0      = data(22,index);       % i0 参考时刻轨道倾角
% doti    = data(26,index);       % IDOT 轨道倾角变化率
    % p.s.
% 86400         一天
% 4.3077e+04    一个周期
% -------------------------------------------
sqrta = sqrt(7229175.73885051);
e     = 0.000394000000000000;
om    = 0.982207577030007;
Om0   = 5.41244749270632;
i0    = 1.23689040627885+1;
doti  = 0;
M0    = 0;
delta_n = 0;            % 很小，可以忽略不计

% -----------仿真参数----------------------
TimeEnd = 13;
tEnd = 40000;
% for t = toe:tspan:SOW(1,TimeEnd)    % 一个运行周期  
for t = toe:tspan:toe + tEnd
    %% 按步骤计算卫星轨道
    %----------1. 计算卫星运动的平均角速度n---------
    n = sqrt(GM) / sqrta^3 + delta_n;
    %----------2. 计算归化时间tk----------
    % 604800为一周的时间
    tk = t - toe;
    if tk > 302400
        tk = tk - 604800;
    elseif tk < -302400
        tk = tk + 604800;
    end
    %----------3. 计算观测时刻的平近点角Mk----------
    Mk = M0 + n * tk + dotM2*tk;
    %----------4. 计算观测时刻的偏近点角Ek----------
    Ek = Mk;
    Etemp = 10e6;
    while 1
        if abs(Ek - Etemp) < 1e-10
            break;
        end
        Etemp = Ek;
        Ek = Mk + e * sin(Ek);
    end
    %----------5. 计算观测时刻的真近点角Vk------------
    Vkc = (cos(Ek)- e )/(1-e*cos(Ek));
    Vks = ((1-(e)^2)^0.5*sin(Ek))/(1-e*cos(Ek));
    Vk = atan2(Vks,Vkc);
    %----------6. 计算升交角距----------
    phk = Vk + om;
    %----------7. 计算摄动改正项----------
    %     setau = data(14,index) * cos(2*phk) + data(16,index) * sin(2*phk);
    %     setar = data(23,index) * cos(2*phk) + data(11,index) * sin(2*phk);
    %     setai = data(19,index) * cos(2*phk) + data(21,index) * sin(2*phk);
    setau = 0;  setar = 0;  setai = 0;
    %----------8.计算经改正后的升交角距、卫星矢径、轨道倾角----------
    uk = phk + setau;
    rk = sqrta^2 * (1 - e * cos(Ek)) + setar;
    ik = i0 + setai + doti * tk;              % 轨道平面倾角
    %----------9. 计算卫星在轨道坐标系的位置----------
    xECI = rk * cos(uk);
    yECI = rk * sin(uk);
    
    %----------10. 计算观测时刻t的升交点赤径----------
    lamdak = Om0 + (dotOm - we) * tk - we * toe;
    %----------11. 计算卫星在WGS-84坐标系的位置----------
    xECEF = xECI * cos(lamdak) - yECI * cos(ik) * sin(lamdak);
    yECEF = xECI * sin(lamdak) + yECI * cos(ik) * cos(lamdak);
    zECEF = yECI * sin(ik);
    %----------11. 输出位置-------------------------
    satPosECEF(:,inx) = [xECEF;yECEF;zECEF];    % 在WGS-84坐标系下的位置，因为有地球自转的存在，所以看起来比较怪
    rECI(:,inx)       = [xECI;yECI;0];          % 轨道坐标系下的位置
    
    %% ---------12. 计算速度--------------------------
    Ekdot  = n/(1-e*sin(Ek));
    phkdot = sqrt((1+e)/(1-e))*(cos(Vk/2)^2)/(cos(Ek/2)^2)*Ekdot;
    %         setau_dot = -2*phkdot* [data(14,index) * sin(2*phk) - data(16,index) * cos(2*phk)];
    %         setar_dot = -2*phkdot* [data(23,index) * sin(2*phk) - data(11,index) * cos(2*phk)];
    %         setai_dot = -2*phkdot* [data(19,index) * sin(2*phk) - data(21,index) * cos(2*phk)];
    setau_dot =0;   setar_dot =0;   setai_dot =0;
    % -----------12.2--------
    dotuk = phkdot + setau_dot;
    dotrk = Ekdot*sqrta^2*e*sin(Ek) + setar_dot;
    dotik = setai_dot + doti;
    % ------------12.3 -------------
    dotLk = dotOm - we;
    xk_dot = dotrk*cos(uk) - dotuk*yECI;
    yk_dot = dotrk*sin(uk) + dotuk*xECI;
    vECI(:,inx) = [xk_dot; yk_dot; 0];

    % ------------12.4 转换坐标系 --------------------
    VxECEF = xk_dot*cos(lamdak) - yk_dot*cos(dotik)*sin(lamdak) + dotik*zECEF*sin(lamdak) - dotLk*yECEF;
    VyECEF = xk_dot*sin(lamdak) + yk_dot*cos(dotik)*cos(lamdak) - dotik*zECEF*cos(lamdak) + dotLk*xECEF;
    VzECEF = yk_dot*sin(ik) + dotik*yECI*cos(dotik);
    vECEF(:,inx)  = [VxECEF;VyECEF;VzECEF];          % 在WGS-84坐标系下的卫星的速度

    % 换一种 (错误的，和坐标系相关，但是计算多普勒是正确的，目前来看也不是很正确。)
    %     r1 = [yECI/xECI;-1;0];
    %     r2 = r1/norm(r1);
    %     vECI(:,inx) = 7425*r2;      % 7425约等于平均速度
    %     vECEF(1,inx) = vECI(1,inx) * cos(lamdak) - vECI(2,inx) * cos(ik) * sin(lamdak);
    %     vECEF(2,inx) = vECI(1,inx) * sin(lamdak) + vECI(2,inx) * cos(ik) * cos(lamdak);
    %     vECEF(3,inx) = vECI(2,inx) * sin(ik);
    %% ------------13. 摄动考虑 -------------------------
    K = n*J2/2*( rE/(sqrta^2*(1-e^2)) )^2;
    dotOm2 = -3*K*cos(ik);
    dotOm  = 0 + dotOm2 ;               % 使用摄动计算出来的
    dotom2 = 3*K*(2 - 5/2*sin(ik)^2);
    om = om + dotom2 * tspan;
    dotM2   = 3/4*n*J2*( 3*cos(ik)^2-1 )/ (1-e^2)^(3/2) * (rE/sqrta^2)^2;
    % --------------14.1 波束选择--------------------------------
    alpha = 80/180*pi;      % 天线的角度
    % --------------14 多普勒速度计算 ---------------------------------
    r_xinToStal = rXGZ - satPosECEF(:,inx);
    r_xinToStal = r_xinToStal/norm(r_xinToStal);
    r_StalToXin =  - r_xinToStal;       % 矢量：信关站指向卫星
    RefCosAlpha = r_StalToXin'*rXGZ / ( norm(r_StalToXin)*norm(rXGZ) );
    
    if ( RefCosAlpha >= cos(alpha) ) && RefCosAlpha>0      % 在天线的波束范围之内
        vUE = [100;100;0];
        vRef = vECEF(:,inx)'*r_StalToXin;
        fs = 19.3e9;   %  载波频率
        c = 3e8;       %  光速
        delta_f(inx) = vRef/c*fs;
        sita(inx) = acos(RefCosAlpha);      % 保存仰角大小
    else                                    % 不在天线的波束范围之内
        sita(inx) = acos(RefCosAlpha);      % 保存仰角大小
    end

    % 更新i
    inx = inx + 1;
    % others
    lamdakSave(inx) = lamdak;   
end
rECILast = satPosECEF(:,end);
% -----------画图 1：频偏------------------------
figure;plot(0:tspan:(0+tEnd),sita/pi*180);
if ~isempty(delta_f);   figure;plot(0:tspan:length(delta_f)-1,delta_f);grid on;xlabel('时间 t s');ylabel('频偏 \deltaf Hz');    end
% -----------画图 2：地球----------------------------
rE     = 6378.137*1e3;      % 地球半径
r_bar = satPosECEF;
figure;plot3(r_bar(1,:),r_bar(2,:),r_bar(3,:),'r','LineWidth',1);
hold on;plot3(rXGZ(1),rXGZ(2),rXGZ(3),'*');grid on;hold on;
[Xe,Ye,Ze] = sphere(24);
Xe = rE*Xe;     Ye = rE*Ye;     Ze = rE*Ze;
surf(Xe,Ye,Ze); axis equal;     hold on;
plot3([0,1.8*rXGZ(1)],[0,1.8*rXGZ(2)],[0,1.8*rXGZ(3)],'LineWidth',2.5,'Color','b');

% -----------画图 3：速度 ----------------------------
% figure;
% plot([0:tspan:(0+tEnd)]/3600,vECEF(1,:));hold on;
% plot([0:tspan:(0+tEnd)]/3600,vECEF(2,:));hold on;
% plot([0:tspan:(0+tEnd)]/3600,vECEF(3,:));
%-----------画图 4：画经纬度轨迹 ----------------------------
% for i = 1:length(satPosECEF)
%     [B_ECEF(i),L_ECEF(i)] = XYZtoBLH(satPosECEF(1,i),satPosECEF(2,i),satPosECEF(3,i));
% end
% figure;
% plot(L_ECEF,B_ECEF,'r.'); % 绘制坐标
% grid on;
% title('卫星飞行轨迹');
% xlabel('经度 rad');
% ylabel('纬度 rad');
% hold on;
% --------------画图 5： -----------------------------




