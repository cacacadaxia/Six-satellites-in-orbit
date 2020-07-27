
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
%  功能:  1. 复制的 LEO2
%        2. 算天线的对准算法
%        3. 求解可以对的卫星颗数
%        4. 波束范围覆盖
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
eE2    = 0.00669437999013;  % 地球离心率的平方
eE     = 0.0818191910428;   % 离心率
rE     = 6378.137*1e3;      % 地球半长轴
J2     = 1082.63e-6;        % 论文中的模型
% ------------------卫星参数-----------------------------
alphaEar = 1250e3 /rE/2;
r3tp  = sqrt(rE^2 + (rE+1175e3)^2 - 2* cos(alphaEar) *rE*(rE+1175e3) );
betaState = asin( sin(alphaEar) *rE/r3tp);            % 用不到，对应文档

% -------------北京-模拟信关站的位置 (东经120，北纬37)---------------------------------
B_bj = (20 +48/60)*pi/180;            % 纬度
L_bj = (118 +28/60-180 )*pi/180;      % 经度
hX = 0;                         % 地面高度
N  = rE / sqrt( 1 - eE2*sin( B_bj )^2 );
% -------------大地地心坐标系&&WGS坐标系-------------------
x_bj = (N+hX)*cos(B_bj)*cos(L_bj);
y_bj = (N+hX)*cos(B_bj)*sin(L_bj);
z_bj = (N*(1-eE2)+hX)*sin(B_bj);
rXGZ = [x_bj;y_bj;z_bj];     % 从地心到信关站的矢量
% ---------------地面坐标系&&WGS坐标系--------------------------
wgs_surface = [cos(L_bj)*sin(B_bj)   sin(L_bj)*sin(B_bj)    -cos(B_bj);
                     -sin(L_bj)               cos(L_bj)           0   ;
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
% -----------------仿真参数-----------------------------
TimeEnd  = 13;                  % 这里没有作用
tEnd     = 6532;                % 一圈，周期时间
NumStall = 18;                  % 一共18个轨道面

%
%% load data of satellite
load('ds/dataTrac');

%% 
% ******************终端还有其他的一些参数***************************
vUE = [0;0;277];                        % NED坐标系
[ BdegUe,LdegUe ] = XYZsettoBLH(rXGZ);  % UE经纬度，度数
vUEwgs = wgs_surface*vUE;

Berr = 5;      Lerr = 5;              % 终端与卫星可以允许的经纬度差
% ******************选星****************************
tspan = 1;
time = 0:tspan:tEnd;
tAcess = [];
LenSet = [];
for inxTim = 1:1:1000       % 具体的某一个时间
    t = time(inxTim);
    % -------- 选星 1 ---------------------------------
    %         cnt = 1;
    %         for faceinx = 1:18
    %             for trainx = 1:48
    %                 rStallWgs = satPosECEF(:,inxTim,trainx,faceinx);     % 18颗卫星的位置，第一个轨道面
    %                 vStallWgs = vECEF(:,inxTim,trainx,faceinx);          % 18颗卫星的速度，第一个轨道面
    %                 [ BdegSta(trainx , faceinx) , LdegSta(trainx , faceinx) ] = XYZsettoBLH( rStallWgs );
    %                 if  (abs(BdegSta(trainx , faceinx) - BdegUe)<Berr ) && (abs(LdegSta(trainx , faceinx) - LdegUe)<Lerr)
    %                     stallSettp(cnt,:) = [ trainx , faceinx ];
    %                     cnt = cnt + 1;
    %                 end
    %             end
    %         end
    % ---------- 选星 2 ---------------------------
    for faceinx = 1:18
        for trainx = 1:48
            rStallWgstp = satPosECEF(:,inxTim,trainx,faceinx);
            alphaRef( trainx , faceinx ) =  rXGZ'*rStallWgstp/( norm(rXGZ)*norm(rStallWgstp) );
        end
    end
    % 选星
    stallSet = intersect( find( alphaRef >=cos( alphaEar ) ) , find( alphaRef>0 ) );        % 同时满足两个条件
    vReftp = [];
    if isempty(stallSet)
        fprintf('没有可以选择的星\n');
        stallSet = [];
    else
        stallSettp = [];
        stallSettp(:,2) = floor((stallSet -0.1)/48)+1;                      % 48里选1个
        stallSettp(:,1) = stallSet - ( stallSettp(:,2)-1 ).*48;             % 18里选1个
        for index = 1:length(stallSet)
            vStall      =      vECEF(:,inxTim,stallSettp(index,1) ,stallSettp(index,2) );   %
            rStallWgstp = satPosECEF(:,inxTim,stallSettp(index,1) ,stallSettp(index,2) );   %
            r_xinToStal = rXGZ - rStallWgstp;
            r_xinToStal = r_xinToStal/norm(r_xinToStal);
            r_StalToXin =  - r_xinToStal;                                   % 矢量: 信关站指向卫星
            vReftp( index ) = (vStall - vUEwgs)'*r_xinToStal;               % 相对速度，用来计算多普勒频偏
        end
        inxReftp = find(vReftp>0);
        if isempty(inxReftp)
            fprintf(' 相对速度，没有大于0的 ');
        end
        statelSect = stallSettp( find(vReftp == max( vReftp ) ), : );       % 找最大的一个
        fprintf('%3.0f\n',length( stallSet ));
    end
    LenSet = [LenSet, length( stallSet )];
    
end
fprintf('波束范围%4.0fkm\n',alphaEar*2*rE/1e3)
[f,xi] = ksdensity(LenSet);
figure;plot(xi,f);xlabel('可选卫星的颗数');ylabel('概率');

% ********************波束范围覆盖************************
% rbeam = alphaEar*rE;
% poinum = 17;
% sitaPoi = 0:pi/8:2*pi;
% TpSet8PoitIQ = exp(1j*sitaPoi);
% TpSet8Poit = [real(exp(1j*sitaPoi));imag(exp(1j*sitaPoi));zeros(1,17)].*rbeam;      % 在地表附近的一个圈，在地球表面坐标系下
% 
% inxTim = 100;
% cout = 1;
% for faceinx = 1:18
%     for trainx = 1:48
%         rStallWgstp = satPosECEF(:,inxTim,trainx,faceinx);
%         alphaRef( trainx , faceinx ) =  rXGZ'*rStallWgstp/( norm(rXGZ)*norm(rStallWgstp) );
%         % 坐标系转换
%         [ BstateTp , LstateTp , ~ ] = XYZsettoBLH(rStallWgstp,1);
%         LstateTp = LstateTp-pi;         % 不知道为什么
%         wgs_surface = [cos(LstateTp)*sin(BstateTp)   sin(LstateTp)*sin(BstateTp)    -cos(BstateTp);
%                                    -sin(LstateTp)               cos(LstateTp)           0   ;
%                      cos(LstateTp)*cos(BstateTp)   sin(LstateTp)*cos(BstateTp)    sin(BstateTp)]';
%         wgs_surface = wgs_surface*diag([-1,1,1]);
%         % 转到WGS坐标系下
%         x_bj = (N + 0) *cos(BstateTp)*cos(LstateTp);
%         y_bj = (N + 0) *cos(BstateTp)*sin(LstateTp);
%         z_bj = (N*(1-eE2)+0)*sin(BstateTp);
%         rStallWgstpSur = [x_bj;y_bj;z_bj];
%         Point8Wgs = wgs_surface*TpSet8Poit + repmat(rStallWgstpSur,1,poinum);
%         poiTpPlt{cout} = Point8Wgs;
%         cout = cout + 1;
%     end
% end

% -----------画图 1.1：卫星波束覆盖范围 ----------------------------
% figure;
% for inx = 1:864
%     b=[];l=[];
%     for i=1:poinum
%         [b(i),l(i),~] = XYZsettoBLH(poiTpPlt{ inx }(:,i));
%     end
%     if min(abs(b))<70&&min(abs(l))<150
%         plot(b,l);hold on;
%     end
% end
% xlabel('纬度');ylabel('经度');title('卫星波束覆盖示意图')

% -----------画图 1.2：地球表面画个圈 (具体没有很大的作用) ----------------------------
% figure;
% rE     = 6378.137*1e3;      % 地球半径
% [Xe,Ye,Ze] = sphere(24);    Xe = rE*Xe;     Ye = rE*Ye;     Ze = rE*Ze;     cla reset;   load topo;
% s = surface(Xe,Ye,Ze,'FaceColor','texturemap','CData',topo); axis equal;colormap(topomap1);brighten(0.9);hold on;
% plot3(poiTpPlt{1}(1,:),poiTpPlt{1}(2,:),poiTpPlt{1}(3,:),'r','LineWidth',1.4);  hold on;
% rStallWgstp = satPosECEF(:,inxTim,1,1);
% plot3(rStallWgstp(1),rStallWgstp(2),rStallWgstp(3),'b*');hold on;
% plot3(rStallWgstpSur(1),rStallWgstpSur(2),rStallWgstpSur(3),'k*')

%% 波束选择 （）
% % --------------14.1 波束选择--------------------------------
% alpha = 80/180*pi;      % 天线的角度，80度是非常大
% % --------------14 多普勒速度计算 ---------------------------------
% r_xinToStal = rXGZ - satPosECEF(:,inx);
% r_xinToStal = r_xinToStal/norm(r_xinToStal);
% r_StalToXin =  - r_xinToStal;       % 矢量：信关站指向卫星
% RefCosAlpha = r_StalToXin'*antDirecWgs / ( norm(r_StalToXin)*norm(antDirecWgs) );     % 计算天线方向角 theta
% 
% if ( RefCosAlpha >= cos(alpha) ) && RefCosAlpha>0      % 在天线的波束范围之内
%     alphaRef = acos( rXGZ'*satPosECEF(:,inx)/( norm(rXGZ)*norm(satPosECEF(:,inx)) ) );      % 卫星和终端相对地心的角
%     if alphaRef <= alphaEar
%         vUE = [0;0;277];                    % 对于车辆来说，只有xy轴有速度，在地球表面坐标系
%         vUEwgs = wgs_surface*vUE;
%         vRef = ( vECEF(:,inx) + vUEwgs )'*r_xinToStal;
%         fs = 19.3e9;                        % 载波频率
%         c = 3e8;                            % 光速
%         delta_f(inx) = vRef/c*fs;           % 保存频偏
%         sita(inx) = acos(RefCosAlpha);      % 保存天线方向角大小
%     end
% else                                        % 不在天线的波束范围之内
%     % sita(inx) = acos(RefCosAlpha);        % 保存天线方向角大小
%     sita(inx) = 0;
% end

%% 画图
% figure;plot(Vkgot);

% rECILast = satPosECEF(:,end);       % 最终的位置
% sita_dot = sita(61:end)-sita(1:end-60);
% % -----------画图 0：天线指向的变化率 ------------------------
% figure;plot(0:tspan:length(sita_dot)-1 , sita_dot./pi*180 ,'r.');grid on;xlabel('时间t s');ylabel('\theta变化率 度数/min');
% % -----------画图 1：频偏------------------------
% figure;plot(0:tspan:(0+tEnd),sita/pi*180);
% if ~isempty(delta_f);   figure;plot(0:tspan:length(delta_f)-1,delta_f);grid on;xlabel('时间 t s');ylabel('频偏 \deltaf Hz');    end
% ----------- 画图 2：地球----------------------------
% figure;
% rE     = 6378.137*1e3;      % 地球半径
% [Xe,Ye,Ze] = sphere(24);    Xe = rE*Xe;     Ye = rE*Ye;     Ze = rE*Ze;     cla reset;   load topo;
% s = surface(Xe,Ye,Ze,'FaceColor','texturemap','CData',topo); axis equal;colormap(topomap1);brighten(0.9);hold on;
% for i = 1:18
%     plot3(satPosECEF(1,:,1,i),satPosECEF(2,:,1,i),satPosECEF(3,:,1,i),'r--','LineWidth',0.1);hold on;
% end
% plot3(rXGZ(1),rXGZ(2),rXGZ(3),'*');grid on;hold on;
% plot3([0,1.2*rXGZ(1)],[0,1.2*rXGZ(2)],[0,1.2*rXGZ(3)],'LineWidth',2.5,'Color','b'); 

% -----------画图 3：所有的卫星动态的运动轨迹 (动图) ----------------------------
% figure;
% rE  = 6378.137*1e3;      % 地球半径
% [Xe,Ye,Ze] = sphere(24);    Xe = rE*Xe;     Ye = rE*Ye;     Ze = rE*Ze;     cla reset;   load topo;
% s = surface(Xe,Ye,Ze,'FaceColor','texturemap','CData',topo); axis equal;colormap(topomap1);brighten(0.9);hold on;
% nowTime = 10;
% h = plot3(0,0,0,'k.','LineWidth',10);
% for nowTime = 1:10:6300
%     set( h,'XData',reshape(satPosECEF(1,nowTime,:,:),1,48*18),'YData',reshape(satPosECEF(2,nowTime,:,:),1,48*18),'ZData',...
%         reshape(satPosECEF(3,nowTime,:,:),1,48*18) );
%     drawnow
%     pause(0.1);
% end

% -----------画图 4：选择卫星 ----------------------------
% figure;
% rE  = 6378.137*1e3;      % 地球半径
% [Xe,Ye,Ze] = sphere(24);    Xe = rE*Xe;     Ye = rE*Ye;     Ze = rE*Ze;     cla reset;   load topo;
% s = surface(Xe,Ye,Ze,'FaceColor','texturemap','CData',topo); axis equal;colormap(topomap1);brighten(0.9);hold on;
% h = plot3(0,0,0,'k.','LineWidth',10);
% set( h,'XData',reshape(satPosECEF(1,inxTim,:,:),1,48*18),'YData',reshape(satPosECEF(2,inxTim,:,:),1,48*18),'ZData',...
%     reshape(satPosECEF(3,inxTim,:,:),1,48*18) );hold on;
% plot3(rXGZ(1),rXGZ(2),rXGZ(3),'*');grid on;hold on;
% plot3([0,1.4*rXGZ(1)],[0,1.4*rXGZ(2)],[0,1.4*rXGZ(3)],'LineWidth',2.5,'Color','b'); hold on;
% for i = 1: size(stallSettp,1)
%     rStallWgstp = satPosECEF(:,inxTim,stallSettp(i,1) ,stallSettp(i,2) );
% %     plot3([rXGZ(1),rStallWgstp(1)] , [rXGZ(2) ,rStallWgstp(2)], [rXGZ(3),rStallWgstp(3)] ,'LineWidth',1,'Color','r');
%     plot3([ rStallWgstp(1)] , [rStallWgstp(2)], [ rStallWgstp(3)] ,'r*');           % 被选择的卫星标记出来
%     hold on;
% end
% rStallWgstp = satPosECEF(:,inxTim,statelSect(1,1) ,statelSect(1,2) );
% plot3( rStallWgstp(1),rStallWgstp(2),rStallWgstp(3) , 'ko' )

% -----------画图 5：波束范围 ----------------------------




