
% =========================================================================
%
%                      ��������λ���ٶȽ���
%
%
% =========================================================================
%
%��(C)2019-2020 ���ݺ���ͨ�����޹�˾
%   �汾��V1.2
%   ���ڣ�2019��8��1��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ����:  1.���Ƶ� LEO1
%        2. ������ǻ�ͼ(ͬһ�����)
%        3. �ҵ����е����ǣ����ұ������ݣ������޸ľ��޸Ĵ��ļ�
%        4.
%        5.
%        6.
%--------------------------------------------------------------------------

clear all;  close all;
rng('default');
tic;
%---------------------��ȡRINEX��ʽn�ļ�--------------------------------
data = RinexNreader('BR.16n','G01'); % ע���ȡ·�������Ǳ��
%---------------------��������յ��ܻ���---------------------------------
% [JD,FOD,GPSW,SOW,DOY,DOW] = GCtoGPS(data(1,1),data(2,1),data(3,1),data(4,1),data(5,1),data(6,1));
%----------����n�ļ��вο���Ԫ���ܻ���----------
for j = 1:size(data,2) % size(data��2) - ����data���������;
    [JD,FOD,GPSW,SOW(j),DOY,DOW] = GCtoGPS(data(1,j),data(2,j),data(3,j),data(4,j),data(5,j),data(6,j));
    % ��ָ�����ǵ� �ꡢ�¡��ա�ʱ���֡��� �� ������(��������)��������(С������)��GPS�ܡ��ܻ��롢����ա�������
end
% -----------------�������---------------------------
GM     = 3986005e8;         % ��=GM
we     = 7.2921151467e-5;   % ������ת���ٶ�
% we     = 0;
eE2    = 0.00669437999013;  % ���������ʵ�ƽ��
eE     = 0.0818191910428;   % ������
rE     = 6378.137*1e3;      % ����볤��
J2     = 1082.63e-6;        % �����е�ģ��
% ------------------���ǲ���-----------------------------
alphaEar = 2000e3 /rE/2;
r3tp  = sqrt(rE^2 + (rE+1175e3)^2 - 2* cos(alphaEar) *rE*(rE+1175e3) );
betaState = asin( sin(alphaEar) *rE/r3tp);            % �ò�������Ӧ�ĵ�

% -------------����-ģ���Ź�վ��λ�� (����120����γ37)---------------------------------
B_bj = (39+48/60)*pi/180;            % γ��
L_bj = (120 +28/60 )*pi/180;      % ����
hX = 0;                         % ����߶�
N  = rE / sqrt( 1 - eE2*sin( B_bj )^2 );
% -------------��ص�������ϵ&&WGS����ϵ-------------------
x_bj = -(N+hX)*cos(B_bj)*cos(L_bj);
y_bj = -(N+hX)*cos(B_bj)*sin(L_bj);
z_bj = (N*(1-eE2)+hX)*sin(B_bj);
rXGZ = [x_bj;y_bj;z_bj];     % �ӵ��ĵ��Ź�վ��ʸ��
% ---------------��������ϵ&&WGS����ϵ--------------------------
wgs_surface = [cos(L_bj)*sin(B_bj)   sin(L_bj)*sin(B_bj)    -cos(B_bj);
                -sin(L_bj)               cos(L_bj)           0        ;
               cos(L_bj)*cos(B_bj)   sin(L_bj)*cos(B_bj)    sin(B_bj)]';
wgs_surface = wgs_surface*diag([-1,1,1]);           % ���춫����ϵ-->���춫����ϵ
% ----------------��������ϵ&&��������ϵ-----------------------
theta = 0;          % ����
fai   = 0;          % ƫ��
gamma = 0;          % ����
ben_surface = [cos(theta)*cos(fai)                                        sin(theta)             -cos(theta)*sin(fai);
                -sin(theta)*cos(fai)*cos(gamma)+sin(fai)*sin(gamma)     cos(theta)*cos(gamma)    sin(theta)*sin(fai)*cos(gamma)+cos(fai)*sin(gamma)
                sin(theta)*cos(fai)*sin(gamma)                          -cos(theta)*sin(gamma)   -sin(theta)*sin(fai)*sin(gamma)+cos(fai)*cos(gamma)];
wgs_ben     = wgs_surface*ben_surface';
antDirec    = [0;0;1];      % �����ڱ�������ϵ�еķ���
antDirecWgs = wgs_ben*antDirec;

%% ���ǹ������
index = 1;
toe   = data(18,index);         % t0e �����ο�ʱ��(�ܻ���)
Om0   = 120/180*pi;             % ��ʱ������趨һ��ֵ
om    = 0.982207577030007;      % ���ص����
% -----------�������----------------------
TimeEnd = 13;                   % ����û������
tEnd = 6532;                    % һȦ
NumStall = 18;                  % һ��18�������
% for staInxFlat = 1:NumStall
%     [satPosECEF(:,:,:,staInxFlat), vECEF(:,:,:,staInxFlat) ] = pltState( tEnd, Om0+2*pi/36*(staInxFlat-1) ,toe , om+2*pi/18*(staInxFlat-1) );
%     disp(staInxFlat);
%     % satPosECEF 3*tEnd*48*18
% end
% save('ds/dataTrac','satPosECEF','vECEF');
%%
load('ds/dataTrac');


%% 
% ******************�ն˻���������һЩ����***************************
vUE = [0;0;277];        % NED����ϵ
[ BdegUe,LdegUe ] = XYZsettoBLH(rXGZ);
Berr = 10;      Lerr = 10;              % �ն������ǿ�������ľ�γ�Ȳ�
% ******************ѡ��****************************
tspan = 1;
time = 0:tspan:tEnd;
tAcess = [];
for inxTim = 1000       % �����ĳһ��ʱ��
    t = time(inxTim);
    rStallWgs = satPosECEF(:,inxTim,1,:);     % 18�����ǵ�λ�ã���һ�������
    vStallWgs = vECEF(:,inxTim,1,:);          % 18�����ǵ��ٶȣ���һ�������
    % -------- ѡ����� ---------------------------------
    for stallinx = 1:NumStall
        [ BdegSta(stallinx) , LdegSta(stallinx) ] = XYZsettoBLH( rStallWgs(:,stallinx) );
    end
    if sum( ( (abs(BdegSta - BdegUe)<Berr ) + (abs(LdegSta - LdegUe)<Lerr) )==2 )
        %        fprintf('���Խ����������\n'); 
        tAcess = [tAcess , t];
    end
    % ---------- ѡһ������������ ---------------------------
    
    
end
q(tAcess) = 1;
figure;plot(q,'.');

%% ����ѡ��
% --------------14.1 ����ѡ��--------------------------------
% alpha = 80/180*pi;      % ���ߵĽǶȣ�80���Ƿǳ���
%     % --------------14 �������ٶȼ��� ---------------------------------
%     r_xinToStal = rXGZ - satPosECEF(:,inx);
%     r_xinToStal = r_xinToStal/norm(r_xinToStal);
%     r_StalToXin =  - r_xinToStal;       % ʸ�����Ź�վָ������
%     RefCosAlpha = r_StalToXin'*antDirecWgs / ( norm(r_StalToXin)*norm(antDirecWgs) );     % �������߷���� theta
%
%     if ( RefCosAlpha >= cos(alpha) ) && RefCosAlpha>0      % �����ߵĲ�����Χ֮��
%         alphaRef = acos( rXGZ'*satPosECEF(:,inx)/( norm(rXGZ)*norm(satPosECEF(:,inx)) ) );      % ���Ǻ��ն���Ե��ĵĽ�
%         if alphaRef <= alphaEar
%             vUE = [0;0;277];                    % ���ڳ�����˵��ֻ��xy�����ٶȣ��ڵ����������ϵ
%             vUEwgs = wgs_surface*vUE;
%             vRef = ( vECEF(:,inx) + vUEwgs )'*r_StalToXin;
%             fs = 19.3e9;                        % �ز�Ƶ��
%             c = 3e8;                            % ����
%             delta_f(inx) = vRef/c*fs;           % ����Ƶƫ
%             sita(inx) = acos(RefCosAlpha);      % �������߷���Ǵ�С
%         end
%     else                                        % �������ߵĲ�����Χ֮��
%         % sita(inx) = acos(RefCosAlpha);        % �������߷���Ǵ�С
%         sita(inx) = 0;
%     end

%%
% figure;plot(Vkgot);

% rECILast = satPosECEF(:,end);       % ���յ�λ��
% sita_dot = sita(61:end)-sita(1:end-60);
% % -----------��ͼ 0������ָ��ı仯�� ------------------------
% figure;plot(0:tspan:length(sita_dot)-1 , sita_dot./pi*180 ,'r.');grid on;xlabel('ʱ��t s');ylabel('\theta�仯�� ����/min');
% % -----------��ͼ 1��Ƶƫ------------------------
% figure;plot(0:tspan:(0+tEnd),sita/pi*180);
% if ~isempty(delta_f);   figure;plot(0:tspan:length(delta_f)-1,delta_f);grid on;xlabel('ʱ�� t s');ylabel('Ƶƫ \deltaf Hz');    end
% -----------��ͼ 2������----------------------------
figure;
rE     = 6378.137*1e3;      % ����뾶
[Xe,Ye,Ze] = sphere(24);    Xe = rE*Xe;     Ye = rE*Ye;     Ze = rE*Ze;     cla reset;   load topo;
s = surface(Xe,Ye,Ze,'FaceColor','texturemap','CData',topo); axis equal;colormap(topomap1);brighten(0.9);hold on;
for i = 1:18
    plot3(satPosECEF(1,:,1,i),satPosECEF(2,:,1,i),satPosECEF(3,:,1,i),'r-','LineWidth',1.3);hold on;
end
plot3(rXGZ(1),rXGZ(2),rXGZ(3),'*');grid on;hold on;
plot3([0,1.4*rXGZ(1)],[0,1.4*rXGZ(2)],[0,1.4*rXGZ(3)],'LineWidth',2.5,'Color','b'); 

% -----------��ͼ 3�����е����� ----------------------------
figure;
rE  = 6378.137*1e3;      % ����뾶
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



%% �������е�����
function [satPosECEF ,vECEF] = pltState( tEnd, Om0 ,toe ,om)
% -----------------�������---------------------------
GM     = 3986005e8;         % ��=GM
we     = 7.2921151467e-5;   % ������ת���ٶ�
eE2    = 0.00669437999013;  % ���������ʵ�ƽ��
eE     = 0.0818191910428;   % ������
rE     = 6378.137*1e3;      % ����볤��
J2     = 1082.63e-6;        % �����е�ģ��
%-----------------------���Ǻ���Ҫ�Ĳ���------------------------------------
inx = 1;
tspan = 1;
dotM2 = 0;
dotOm = 0;
delta_f = [];
% -------------------------------------------
sqrta = sqrt(1175e3 + rE);      % Ŀǰ�͹�Ĳ���
e     = 0.000394000000000000;   % ������
% om    = 0.982207577030007;      % ���ص����
% Om0   = 130/180*pi;
i0    = 86.7/180*pi;            % �����ǣ�ȷ��ֵ
doti  = 0;
M0    = 0;
delta_n = 0;            % ��С�����Ժ��Բ���
% -----------�������----------------------
% tEnd = 15000;
% for t = toe:tspan:SOW(1,TimeEnd)    % һ����������
for t = toe:tspan:toe + tEnd
    %% ������������ǹ��
    %----------1. ���������˶���ƽ�����ٶ�n---------
    n = sqrt(GM) / sqrta^3 + delta_n;
    %----------2. ����黯ʱ��tk----------
    % 604800Ϊһ�ܵ�ʱ��
    tk = t - toe;
    if tk > 302400
        tk = tk - 604800;
    elseif tk < -302400
        tk = tk + 604800;
    end
    %----------3. ����۲�ʱ�̵�ƽ�����Mk----------
    Mk = M0 + n * tk + dotM2*tk;
    %----------4. ����۲�ʱ�̵�ƫ�����Ek----------
    Ek = Mk;
    Etemp = 10e6;
    while 1
        if abs(Ek - Etemp) < 1e-10
            break;
        end
        Etemp = Ek;
        Ek = Mk + e * sin(Ek);
    end
    %----------5. ����۲�ʱ�̵�������Vk------------
    Vkc = (cos(Ek)- e )/(1-e*cos(Ek));
    Vks = ((1-(e)^2)^0.5*sin(Ek))/(1-e*cos(Ek));
    Vk  = atan2(Vks,Vkc);
    % ----------5.1 �ȷ�����ǵ�-----------------
    num1Surf = 48;
    VkSav = Vk + [0:2*pi/num1Surf:2*pi-2*pi/num1Surf];
    % -----------5.2 һ�������ֳ�48������-----------------------------
    for FacInxS = 1:num1Surf
        %----------6. ���������Ǿ�----------
        phk = VkSav( FacInxS ) + om;
        %----------7. �����㶯������----------
        %     setau = data(14,index) * cos(2*phk) + data(16,index) * sin(2*phk);
        %     setar = data(23,index) * cos(2*phk) + data(11,index) * sin(2*phk);
        %     setai = data(19,index) * cos(2*phk) + data(21,index) * sin(2*phk);
        setau = 0;  setar = 0;  setai = 0;
        %----------8.���㾭������������Ǿࡢ����ʸ����������----------
        uk = phk + setau;
        rk = sqrta^2 * (1 - e * cos(Ek)) + setar;
        ik = i0 + setai + doti * tk;              % ���ƽ�����
        %----------9. ���������ڹ������ϵ��λ��----------
        xECI = rk * cos(uk);
        yECI = rk * sin(uk);
        %----------10. ����۲�ʱ��t��������ྶ----------
        lamdak = Om0 + (dotOm - we) * tk - we * toe;
        %----------11. �������ǵ�λ��--WGS-84����ϵ----------
        xECEF = xECI * cos(lamdak) - yECI * cos(ik) * sin(lamdak);
        yECEF = xECI * sin(lamdak) + yECI * cos(ik) * cos(lamdak);
        zECEF = yECI * sin(ik);
        %----------11. ���λ�� -------------------------
        satPosECEF(:,inx,FacInxS) = [xECEF;yECEF;zECEF];    % ��WGS-84����ϵ�µ�λ�ã���Ϊ�е�����ת�Ĵ��ڣ����Կ������ȽϹ�
        %         rECI(:,inx)       = [xECI;yECI;0];          % �������ϵ�µ�λ��
        %% ---------12. �����ٶ�--------------------------
        Ekdot  = n/(1-e*sin(Ek));
        phkdot = sqrt((1+e)/(1-e))*(cos(Vk/2)^2)/(cos(Ek/2)^2)*Ekdot;
        %         setau_dot = -2*phkdot* [data(14,index) * sin(2*phk) - data(16,index) * cos(2*phk)];
        %         setar_dot = -2*phkdot* [data(23,index) * sin(2*phk) - data(11,index) * cos(2*phk)];
        %         setai_dot = -2*phkdot* [data(19,index) * sin(2*phk) - data(21,index) * cos(2*phk)];
        setau_dot =0;   setar_dot =0;   setai_dot =0;
        % -----------12.2 ����仯��-------------------------------------
        dotuk = phkdot + setau_dot;
        dotrk = Ekdot*sqrta^2*e*sin(Ek) + setar_dot;
        dotik = setai_dot + doti;
        % ------------12.3 �����ٶ�--�������ϵ---------------------------
        dotLk = dotOm - we;
        xk_dot = dotrk*cos(uk) - dotuk*yECI;
        yk_dot = dotrk*sin(uk) + dotuk*xECI;
        % vECI(:,inx) = [xk_dot; yk_dot; 0];
        
        % ------------12.4 ת������ϵ--WGS����ϵ --------------------------
        VxECEF = xk_dot*cos(lamdak) - yk_dot*cos(dotik)*sin(lamdak) + dotik*zECEF*sin(lamdak) - dotLk*yECEF;
        VyECEF = xk_dot*sin(lamdak) + yk_dot*cos(dotik)*cos(lamdak) - dotik*zECEF*cos(lamdak) + dotLk*xECEF;
        VzECEF = yk_dot*sin(ik) + dotik*yECI*cos(dotik);
        vECEF(:,inx,FacInxS)  = [VxECEF;VyECEF;VzECEF];          % ��WGS-84����ϵ�µ����ǵ��ٶ�
    end
    %% ------------13. �㶯���� -------------------------
    K = n*J2/2*( rE/(sqrta^2*(1-e^2)) )^2;
    dotOm2 = -3*K*cos(ik);
    dotOm  = 0 + dotOm2 ;               % ʹ���㶯���������
    dotom2 = 3*K*(2 - 5/2*sin(ik)^2);       % �������ûɶ��
    %     om     = om + dotom2 * tspan;     % ��ʱ��һ���Ȳ�Ҫ��
    dotM2  = 3/4*n*J2*( 3*cos(ik)^2-1 )/ (1-e^2)^(3/2) * (rE/sqrta^2)^2;
    % ����i
    inx = inx + 1;
    % others
    % lamdakSave(inx) = lamdak;
    
end
end





