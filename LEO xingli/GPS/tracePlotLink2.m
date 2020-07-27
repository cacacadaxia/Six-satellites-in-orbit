
% =========================================================================
%
%                      导航卫星位置速度解算
%
%
% =========================================================================
%
%　(C)2019-2020 广州海格通信有限公司
%   版本：V1.2
%   日期：2019年7月25日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能:  1.使用北斗的星历数据，画出北斗卫星的飞行轨迹
%        2.整合到一个文件里来，增加速度计算
%        3.
%        4.
%        5.
%        6.
%--------------------------------------------------------------------------

clear all;
close all;
%---------------------读取RINEX格式n文件--------------------------------
data = RinexNreader('BR.16n','G0'); % 注意读取路径和卫星编号
%---------------------计算测量日的周积秒---------------------------------
[JD,FOD,GPSW,SOW,DOY,DOW] = GCtoGPS(data(1,1),data(2,1),data(3,1),data(4,1),data(5,1),data(6,1));
t0 = SOW; % SOW-周积秒
%-----------------------计算卫星坐标------------------------------------
inx = 1;
for t = t0:1*60:(t0+86400)
    % 从 t0 到 t0+86400秒(1天) 间隔 600秒(10min)，决定轨迹疏密
    % [satPosECEF(:,i) , rECI(:,i) ] = orbitDetermine(data,t);
    GM = 3986005e8; % μ=GM
    we = 7.2921151467e-5; % 地球自转角速度
%     we = 0;
    %----------计算n文件中参考历元的周积秒----------
    for j = 1:size(data,2) % size(data，2) - 返回data矩阵的列数;
        [JD,FOD,GPSW,SOW(j),DOY,DOW] = GCtoGPS(2000+data(1,j),data(2,j),data(3,j),data(4,j),data(5,j),data(6,j));
        % 将指定卫星的 年、月、日、时、分、秒 → 儒略日(整数部分)、儒略日(小数部分)、GPS周、周积秒、年积日、星期数
    end
    %----------找出与t最近的那列数据----------
    index = 1;      % 暂时先让 i==1
    if t >= SOW(size(data,2))
        index = size(data,2);
    elseif t < SOW(1)
    else
        while t >= SOW(index+1)     % 如果时间超过SOW(i)，那么就找SOW(i+1)的时间
            index = index + 1;
        end
    end
    %% 首先进行赋值
    delta_n = data(12,index);
    a       = data(17,index);
    toe     = data(18,index);
    M0      = data(13,index);
    e       = data(15,index);
    om      = data(24,index);
    Om0     = data(20,index);       % OMEGA0 非参考时刻toe升交点赤经，是始于格林尼治子午圈到卫星轨道升交点所谓准经度
    dotOm   = data(25,index);
    i0      = data(22,index);
    doti    = data(26,index);
    %% 按步骤计算卫星轨道
    %----------1. 计算卫星运动的平均角速度n---------
    n = sqrt(GM) / a^3 + delta_n;
    %----------2. 计算归化时间tk----------
    % 604800为一周的时间
    tk = t - toe;
    if tk > 302400
        tk = tk - 604800;
    elseif tk < -302400
        tk = tk + 604800;
    end
    %----------3. 计算观测时刻的平近点角Mk----------
    Mk = M0 + n * tk;
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
    setau = data(14,index) * cos(2*phk) + data(16,index) * sin(2*phk);
    setar = data(23,index) * cos(2*phk) + data(11,index) * sin(2*phk);
    setai = data(19,index) * cos(2*phk) + data(21,index) * sin(2*phk);
    %----------8.计算经改正后的升交角距、卫星矢径、轨道倾角----------
    uk = phk + setau;
    rk = a^2 * (1 - e * cos(Ek)) + setar;
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
    setau_dot = -2*phkdot* [data(14,index) * sin(2*phk) - data(16,index) * cos(2*phk)];
    setar_dot = -2*phkdot* [data(23,index) * sin(2*phk) - data(11,index) * cos(2*phk)];
    setai_dot = -2*phkdot* [data(19,index) * sin(2*phk) - data(21,index) * cos(2*phk)];    
    % -----------12.2--------
    dotuk = phkdot + setau_dot;
    dotrk = Ekdot*a*e*sin(Ek) + setar_dot;
    dotik = setai_dot + doti;
    % ------------12.3 -------------
    dotLk = dotOm - we;
    xk_dot = dotrk*cos(uk) - dotuk*yECI;
    yk_dot = dotrk*sin(uk) + dotuk*xECI;
    % ------------12.4 转换坐标系 --------------------
    VxECEF = xk_dot*cos(lamdak) - yk_dot*cos(dotik)*sin(lamdak) + dotik*zECEF*sin(lamdak) - dotLk*yECEF;
    VyECEF = xk_dot*sin(lamdak) + yk_dot*cos(dotik)*cos(lamdak) - dotik*zECEF*cos(lamdak) + dotLk*xECEF;
    VzECEF = yk_dot*sin(ik) + dotik*yECI*cos(dotik);
    
    vECI(:,inx)  = [VxECEF;VyECEF;VzECEF];          % 在WGS-84坐标系下的速度
    % 更新i
    inx = inx + 1;
    % others
    lamdakSave(inx) = lamdak;
end

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

figure;plot(vECI(1,:))
figure;plot(0:1*60:(0+86400)+60,lamdakSave)