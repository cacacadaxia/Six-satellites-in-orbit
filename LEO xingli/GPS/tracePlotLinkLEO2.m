
% =========================================================================
%
%                      导航卫星位置速度解算
%
%
% =========================================================================
%
%　(C)2019-2020 广州海格通信有限公司
%   版本：V1.2
%   日期：2019年8月1日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能:  1.复制的 LEO1
%        2. 多颗卫星画图(同一轨道面)
%        3. 找到所有的卫星，并且保存数据，后期修改就修改此文件
%        4.
%        5.
%        6.
%--------------------------------------------------------------------------

clear all;  close all;
rng('default');
tic;
%---------------------读取RINEX格式n文件--------------------------------
data = RinexNreader('BR.16n','G01'); % 注意读取路径和卫星编号
%---------------------计算测量日的周积秒---------------------------------
% [JD,FOD,GPSW,SOW,DOY,DOW] = GCtoGPS(data(1,1),data(2,1),data(3,1),data(4,1),data(5,1),data(6,1));
%----------计算n文件中参考历元的周积秒----------
for j = 1:size(data,2) % size(data，2) - 返回data矩阵的列数;
    [JD,FOD,GPSW,SOW(j),DOY,DOW] = GCtoGPS(data(1,j),data(2,j),data(3,j),data(4,j),data(5,j),data(6,j));
    % 将指定卫星的 年、月、日、时、分、秒 → 儒略日(整数部分)、儒略日(小数部分)、GPS周、周积秒、年积日、星期数
end
% -----------------地球参数---------------------------
GM     = 3986005e8;         % μ=GM
we     = 7.2921151467e-5;   % 地球自转角速度
% we     = 0;
eE2    = 0.00669437999013;  % 地球离心率的平方
eE     = 0.0818191910428;   % 离心率
rE     = 6378.137*1e3;      % 地球半长轴
J2     = 1082.63e-6;        % 论文中的模型
% ------------------卫星参数-----------------------------
alphaEar = 2000e3 /rE/2;
r3tp  = sqrt(rE^2 + (rE+1175e3)^2 - 2* cos(alphaEar) *rE*(rE+1175e3) );
betaState = asin( sin(alphaEar) *rE/r3tp);            % 用不到，对应文档

% -------------北京-模拟信关站的位置 (东经120，北纬37)---------------------------------
B_bj = (39+48/60)*pi/180;            % 纬度
L_bj = (120 +28/60 )*pi/180;      % 经度
hX = 0;                         % 地面高度
N  = rE / sqrt( 1 - eE2*sin( B_bj )^2 );
% -------------大地地心坐标系&&WGS坐标系-------------------
x_bj = -(N+hX)*cos(B_bj)*cos(L_bj);
y_bj = -(N+hX)*cos(B_bj)*sin(L_bj);
z_bj = (N*(1-eE2)+hX)*sin(B_bj);
rXGZ = [x_bj;y_bj;z_bj];     % 从地心到信关站的矢量
% ---------------地面坐标系&&WGS坐标系--------------------------
wgs_surface = [cos(L_bj)*sin(B_bj)   sin(L_bj)*sin(B_bj)    -cos(B_bj);
                -sin(L_bj)               cos(L_bj)           0        ;
               cos(L_bj)*cos(B_bj)   sin(L_bj)*cos(B_bj)    sin(B_bj)]';
wgs_surface = wgs_surface*diag([-1,1,1]);           % 南天东坐标系-->北天东坐标系
% ----------------本体坐标系&&地面坐标系-----------------------
theta = 0;          % 俯仰
fai   = 0;          % 偏航
gamma = 0;          % 滚降
ben_surface = [cos(theta)*cos(fai)                                        sin(theta)             -cos(theta)*sin(fai);
                -sin(theta)*cos(fai)*cos(gamma)+sin(fai)*sin(gamma)     cos(theta)*cos(gamma)    sin(theta)*sin(fai)*cos(gamma)+cos(fai)*sin(gamma)
                sin(theta)*cos(fai)*sin(gamma)                          -cos(theta)*sin(gamma)   -sin(theta)*sin(fai)*sin(gamma)+cos(fai)*cos(gamma)];
wgs_ben     = wgs_surface*ben_surface';
antDirec    = [0;0;1];      % 天线在本体坐标系中的方向
antDirecWgs = wgs_ben*antDirec;

%% 卫星轨道参数
index = 1;
toe   = data(18,index);         % t0e 星历参考时刻(周积秒)
Om0   = 120/180*pi;             % 暂时先随便设定一个值
om    = 0.982207577030007;      % 近地点幅角
% -----------仿真参数----------------------
TimeEnd = 13;                   % 这里没有作用
tEnd = 6532;                    % 一圈
NumStall = 18;                  % 一共18个轨道面
% for staInxFlat = 1:NumStall
%     [satPosECEF(:,:,:,staInxFlat), vECEF(:,:,:,staInxFlat) ] = pltState( tEnd, Om0+2*pi/36*(staInxFlat-1) ,toe , om+2*pi/18*(staInxFlat-1) );
%     disp(staInxFlat);
%     % satPosECEF 3*tEnd*48*18
% end
% save('ds/dataTrac','satPosECEF','vECEF');
%%
load('ds/dataTrac');


%% 
% ******************终端还有其他的一些参数***************************
vUE = [0;0;277];        % NED坐标系
[ BdegUe,LdegUe ] = XYZsettoBLH(rXGZ);
Berr = 10;      Lerr = 10;              % 终端与卫星可以允许的经纬度差
% ******************选星****************************
tspan = 1;
time = 0:tspan:tEnd;
tAcess = [];
for inxTim = 1000       % 具体的某一个时间
    t = time(inxTim);
    rStallWgs = satPosECEF(:,inxTim,1,:);     % 18颗卫星的位置，第一个轨道面
    vStallWgs = vECEF(:,inxTim,1,:);          % 18颗卫星的速度，第一个轨道面
    % -------- 选轨道面 ---------------------------------
    for stallinx = 1:NumStall
        [ BdegSta(stallinx) , LdegSta(stallinx) ] = XYZsettoBLH( rStallWgs(:,stallinx) );
    end
    if sum( ( (abs(BdegSta - BdegUe)<Berr ) + (abs(LdegSta - LdegUe)<Lerr) )==2 )
        %        fprintf('可以接入这个星了\n'); 
        tAcess = [tAcess , t];
    end
    % ---------- 选一个轨道面的卫星 ---------------------------
    
    
end
q(tAcess) = 1;
figure;plot(q,'.');

%% 波束选择
% --------------14.1 波束选择--------------------------------
% alpha = 80/180*pi;      % 天线的角度，80度是非常大
%     % --------------14 多普勒速度计算 ---------------------------------
%     r_xinToStal = rXGZ - satPosECEF(:,inx);
%     r_xinToStal = r_xinToStal/norm(r_xinToStal);
%     r_StalToXin =  - r_xinToStal;       % 矢量：信关站指向卫星
%     RefCosAlpha = r_StalToXin'*antDirecWgs / ( norm(r_StalToXin)*norm(antDirecWgs) );     % 计算天线方向角 theta
%
%     if ( RefCosAlpha >= cos(alpha) ) && RefCosAlpha>0      % 在天线的波束范围之内
%         alphaRef = acos( rXGZ'*satPosECEF(:,inx)/( norm(rXGZ)*norm(satPosECEF(:,inx)) ) );      % 卫星和终端相对地心的角
%         if alphaRef <= alphaEar
%             vUE = [0;0;277];                    % 对于车辆来说，只有xy轴有速度，在地球表面坐标系
%             vUEwgs = wgs_surface*vUE;
%             vRef = ( vECEF(:,inx) + vUEwgs )'*r_StalToXin;
%             fs = 19.3e9;                        % 载波频率
%             c = 3e8;                            % 光速
%             delta_f(inx) = vRef/c*fs;           % 保存频偏
%             sita(inx) = acos(RefCosAlpha);      % 保存天线方向角大小
%         end
%     else                                        % 不在天线的波束范围之内
%         % sita(inx) = acos(RefCosAlpha);        % 保存天线方向角大小
%         sita(inx) = 0;
%     end

%%
% figure;plot(Vkgot);

% rECILast = satPosECEF(:,end);       % 最终的位置
% sita_dot = sita(61:end)-sita(1:end-60);
% % -----------画图 0：天线指向的变化率 ------------------------
% figure;plot(0:tspan:length(sita_dot)-1 , sita_dot./pi*180 ,'r.');grid on;xlabel('时间t s');ylabel('\theta变化率 度数/min');
% % -----------画图 1：频偏------------------------
% figure;plot(0:tspan:(0+tEnd),sita/pi*180);
% if ~isempty(delta_f);   figure;plot(0:tspan:length(delta_f)-1,delta_f);grid on;xlabel('时间 t s');ylabel('频偏 \deltaf Hz');    end
% -----------画图 2：地球----------------------------
figure;
rE     = 6378.137*1e3;      % 地球半径
[Xe,Ye,Ze] = sphere(24);    Xe = rE*Xe;     Ye = rE*Ye;     Ze = rE*Ze;     cla reset;   load topo;
s = surface(Xe,Ye,Ze,'FaceColor','texturemap','CData',topo); axis equal;colormap(topomap1);brighten(0.9);hold on;
for i = 1:18
    plot3(satPosECEF(1,:,1,i),satPosECEF(2,:,1,i),satPosECEF(3,:,1,i),'r-','LineWidth',1.3);hold on;
end
plot3(rXGZ(1),rXGZ(2),rXGZ(3),'*');grid on;hold on;
plot3([0,1.4*rXGZ(1)],[0,1.4*rXGZ(2)],[0,1.4*rXGZ(3)],'LineWidth',2.5,'Color','b'); 

% -----------画图 3：所有的卫星 ----------------------------
figure;
rE  = 6378.137*1e3;      % 地球半径
[Xe,Ye,Ze] = sphere(24);    Xe = rE*Xe;     Ye = rE*Ye;     Ze = rE*Ze;     cla reset;   load topo;
s = surface(Xe,Ye,Ze,'FaceColor','texturemap','CData',topo); axis equal;colormap(topomap1);brighten(0.9);hold on;
nowTime = 10;
h = plot3(0,0,0,'k.','LineWidth',10);
for nowTime = 1:10:6300
    set( h,'XData',reshape(satPosECEF(1,nowTime,:,:),1,48*18),'YData',reshape(satPosECEF(2,nowTime,:,:),1,48*18),'ZData',...
        reshape(satPosECEF(3,nowTime,:,:),1,48*18) );
    drawnow
    pause(0.1);
end


toc



%% 计算所有的卫星
function [satPosECEF ,vECEF] = pltState( tEnd, Om0 ,toe ,om)
% -----------------地球参数---------------------------
GM     = 3986005e8;         % μ=GM
we     = 7.2921151467e-5;   % 地球自转角速度
eE2    = 0.00669437999013;  % 地球离心率的平方
eE     = 0.0818191910428;   % 离心率
rE     = 6378.137*1e3;      % 地球半长轴
J2     = 1082.63e-6;        % 论文中的模型
%-----------------------不是很重要的参数------------------------------------
inx = 1;
tspan = 1;
dotM2 = 0;
dotOm = 0;
delta_f = [];
% -------------------------------------------
sqrta = sqrt(1175e3 + rE);      % 目前低轨的参数
e     = 0.000394000000000000;   % 离心率
% om    = 0.982207577030007;      % 近地点幅角
% Om0   = 130/180*pi;
i0    = 86.7/180*pi;            % 轨道倾角，确定值
doti  = 0;
M0    = 0;
delta_n = 0;            % 很小，可以忽略不计
% -----------仿真参数----------------------
% tEnd = 15000;
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
    Vk  = atan2(Vks,Vkc);
    % ----------5.1 等分真近角点-----------------
    num1Surf = 48;
    VkSav = Vk + [0:2*pi/num1Surf:2*pi-2*pi/num1Surf];
    % -----------5.2 一个轨道面分成48颗卫星-----------------------------
    for FacInxS = 1:num1Surf
        %----------6. 计算升交角距----------
        phk = VkSav( FacInxS ) + om;
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
        %----------11. 计算卫星的位置--WGS-84坐标系----------
        xECEF = xECI * cos(lamdak) - yECI * cos(ik) * sin(lamdak);
        yECEF = xECI * sin(lamdak) + yECI * cos(ik) * cos(lamdak);
        zECEF = yECI * sin(ik);
        %----------11. 输出位置 -------------------------
        satPosECEF(:,inx,FacInxS) = [xECEF;yECEF;zECEF];    % 在WGS-84坐标系下的位置，因为有地球自转的存在，所以看起来比较怪
        %         rECI(:,inx)       = [xECI;yECI;0];          % 轨道坐标系下的位置
        %% ---------12. 计算速度--------------------------
        Ekdot  = n/(1-e*sin(Ek));
        phkdot = sqrt((1+e)/(1-e))*(cos(Vk/2)^2)/(cos(Ek/2)^2)*Ekdot;
        %         setau_dot = -2*phkdot* [data(14,index) * sin(2*phk) - data(16,index) * cos(2*phk)];
        %         setar_dot = -2*phkdot* [data(23,index) * sin(2*phk) - data(11,index) * cos(2*phk)];
        %         setai_dot = -2*phkdot* [data(19,index) * sin(2*phk) - data(21,index) * cos(2*phk)];
        setau_dot =0;   setar_dot =0;   setai_dot =0;
        % -----------12.2 计算变化率-------------------------------------
        dotuk = phkdot + setau_dot;
        dotrk = Ekdot*sqrta^2*e*sin(Ek) + setar_dot;
        dotik = setai_dot + doti;
        % ------------12.3 计算速度--轨道坐标系---------------------------
        dotLk = dotOm - we;
        xk_dot = dotrk*cos(uk) - dotuk*yECI;
        yk_dot = dotrk*sin(uk) + dotuk*xECI;
        % vECI(:,inx) = [xk_dot; yk_dot; 0];
        
        % ------------12.4 转换坐标系--WGS坐标系 --------------------------
        VxECEF = xk_dot*cos(lamdak) - yk_dot*cos(dotik)*sin(lamdak) + dotik*zECEF*sin(lamdak) - dotLk*yECEF;
        VyECEF = xk_dot*sin(lamdak) + yk_dot*cos(dotik)*cos(lamdak) - dotik*zECEF*cos(lamdak) + dotLk*xECEF;
        VzECEF = yk_dot*sin(ik) + dotik*yECI*cos(dotik);
        vECEF(:,inx,FacInxS)  = [VxECEF;VyECEF;VzECEF];          % 在WGS-84坐标系下的卫星的速度
    end
    %% ------------13. 摄动考虑 -------------------------
    K = n*J2/2*( rE/(sqrta^2*(1-e^2)) )^2;
    dotOm2 = -3*K*cos(ik);
    dotOm  = 0 + dotOm2 ;               % 使用摄动计算出来的
    dotom2 = 3*K*(2 - 5/2*sin(ik)^2);       % 那这个就没啥用
    %     om     = om + dotom2 * tspan;     % 暂时这一步先不要算
    dotM2  = 3/4*n*J2*( 3*cos(ik)^2-1 )/ (1-e^2)^(3/2) * (rE/sqrta^2)^2;
    % 更新i
    inx = inx + 1;
    % others
    % lamdakSave(inx) = lamdak;
    
end
end





