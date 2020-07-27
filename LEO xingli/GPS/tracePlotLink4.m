
% =========================================================================
%
%                      导航卫星位置速度解算
%
%
% =========================================================================
%
%　(C)2019-2020 广州海格通信有限公司
%   版本：V1.4
%   日期：2019年7月25日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能:  1.使用北斗的星历数据，画出北斗卫星的飞行轨迹
%        2.整合到一个文件里来，增加速度计算
%        3.摄动考虑
%        4.切入到低轨卫星，比较准确星历与计算之间的误差
%        5.经过验证，模型无误
%        6.
%--------------------------------------------------------------------------

clear all;  close all;
%---------------------读取RINEX格式n文件--------------------------------
data = RinexNreader('BR.16n','G22'); % 注意读取路径和卫星编号
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
% we = 0;
eE2    = 0.00669437999013;  % 地球离心率的平方
eE     = 0.0818191910428;   % 离心率
rE     = 6378.137*1e3;      % 地球半径
J2     = 1082.63e-6;        % 论文中的模型
%-----------------------计算卫星坐标------------------------------------
inx = 1;
tspan = 1;
dotM2 = 0;
dotOm = 0;


%% 首先进行赋值
index   = 1;
delta_n = data(12,index);       % detn 平近点角的长期变化(近地点参数)
sqrta   = data(17,index);       % sqrta 长半轴的平方根
toe     = data(18,index);       % t0e 星历参考时刻(周积秒)
M0      = data(13,index);       % M0 参考时刻的平近点角
e       = data(15,index);       % e 偏心率
om      = data(24,index);       % omega 近地点角距
Om0     = data(20,index);       % OMEGA0 非参考时刻toe升交点赤经，是始于格林尼治子午圈到卫星轨道升交点所谓准经度
% dotOm   = data(25,index);       % OMEGADOT 升交点赤经在赤道平面中的长期变化
i0      = data(22,index);       % i0 参考时刻轨道倾角
doti    = data(26,index);       % IDOT 轨道倾角变化率

% 86400         一天
% 4.3077e+04    一个周期
% -----------仿真参数----------------------
TimeEnd = 13;
for t = toe:tspan:SOW(1,TimeEnd)    % 一个运行周期  

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
    %         setau = data(14,index) * cos(2*phk) + data(16,index) * sin(2*phk);
    %         setar = data(23,index) * cos(2*phk) + data(11,index) * sin(2*phk);
    %         setai = data(19,index) * cos(2*phk) + data(21,index) * sin(2*phk);
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
%                 setau_dot = -2*phkdot* [data(14,index) * sin(2*phk) - data(16,index) * cos(2*phk)];
%                 setar_dot = -2*phkdot* [data(23,index) * sin(2*phk) - data(11,index) * cos(2*phk)];
%                 setai_dot = -2*phkdot* [data(19,index) * sin(2*phk) - data(21,index) * cos(2*phk)];
        setau_dot =0;setar_dot =0;setai_dot =0;
        % -----------12.2--------
        dotuk = phkdot + setau_dot;
        dotrk = Ekdot*sqrta^2*e*sin(Ek) + setar_dot;
        dotik = setai_dot + doti;
        % ------------12.3 -------------
        dotLk = dotOm - we;
        xk_dot = dotrk*cos(uk) - dotuk*yECI;
        yk_dot = dotrk*sin(uk) + dotuk*xECI;
        % ------------12.4 转换坐标系 --------------------
        VxECEF = xk_dot*cos(lamdak) - yk_dot*cos(dotik)*sin(lamdak) + dotik*zECEF*sin(lamdak) - dotLk*yECEF;
        VyECEF = xk_dot*sin(lamdak) + yk_dot*cos(dotik)*cos(lamdak) - dotik*zECEF*cos(lamdak) + dotLk*xECEF;
        VzECEF = yk_dot*sin(ik) + dotik*yECI*cos(dotik);
        vECEF(:,inx)  = [VxECEF;VyECEF;VzECEF];          % 在WGS-84坐标系下的速度
    
        %% ------------13. 摄动考虑 -------------------------
        K = n*J2/2*( rE/(sqrta^2*(1-e^2)) )^2;
        dotOm2 = -3*K*cos(ik);
        dotOm  = 0 + dotOm2   ;               % 使用摄动计算出来的
        dotom2 = 3*K*(2 - 5/2*sin(ik)^2);
%         om = om + dotom2*tspan ;               % 为什么这个量出现了没有什么用
        dotM2   = 3/4*n*J2*( 3*cos(ik)^2-1 )/ (1-e^2)^(3/2) * (rE/sqrta^2)^2;
        omsave(inx) = om;
%         
    % 更新i
    inx = inx + 1;
    % others
    lamdakSave(inx) = lamdak;   
end

rECILast = satPosECEF(:,end);
% -----------画图----------------------------
rE     = 6378.137*1e3;      % 地球半径
r_bar = satPosECEF;
figure(1);plot3(r_bar(1,:),r_bar(2,:),r_bar(3,:),'r','LineWidth',1);
hold on;plot3(0,0,0,'*');grid on;
hold on;
[Xe,Ye,Ze] = sphere(24);
Xe = rE*Xe;
Ye = rE*Ye;
Ze = rE*Ze;
surf(Xe,Ye,Ze)
axis equal;
% ---------------------------------------------
% figure;
% plot(rECI(1,:),rECI(2,:))

figure;plot([toe:tspan:SOW(1,TimeEnd)]/3600,vECEF(1,:));hold on;
plot([toe:tspan:SOW(1,TimeEnd)]/3600,vECEF(2,:));hold on;
plot([toe:tspan:SOW(1,TimeEnd)]/3600,vECEF(3,:));
% figure;plot(0:1*60:(0+86400)+60,lamdakSave)



%% 计算下一时刻（从SOW里选择）的位置
[satPos,rECI] = orbitDetermine(data,SOW(TimeEnd));
fprintf('经历时间为：%4.4f h , 星历计算出来的和推演的轨道相距位置为：%4.4fkm\n',...
    (SOW(1,TimeEnd)-toe)/3600,norm([rECILast - satPos'])/1e3)

errR = rECILast - satPos';
fprintf('在速度方向上的误差为：%4.4fkm \n\n',...
errR'*vECEF(:,end)/norm(vECEF(:,end))/1000)

%% 测试
figure;plot([toe:tspan:SOW(1,TimeEnd)]/3600 , omsave)



