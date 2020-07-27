
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
%  ����:  1. ���Ƶ� LEO2
%        2. �����ߵĶ�׼�㷨
%        3. �����ԶԵ����ǿ���
%        4. ������Χ����
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
eE2    = 0.00669437999013;  % ���������ʵ�ƽ��
eE     = 0.0818191910428;   % ������
rE     = 6378.137*1e3;      % ����볤��
J2     = 1082.63e-6;        % �����е�ģ��
% ------------------���ǲ���-----------------------------
alphaEar = 1250e3 /rE/2;
r3tp  = sqrt(rE^2 + (rE+1175e3)^2 - 2* cos(alphaEar) *rE*(rE+1175e3) );
betaState = asin( sin(alphaEar) *rE/r3tp);            % �ò�������Ӧ�ĵ�

% -------------����-ģ���Ź�վ��λ�� (����120����γ37)---------------------------------
B_bj = (20 +48/60)*pi/180;            % γ��
L_bj = (118 +28/60-180 )*pi/180;      % ����
hX = 0;                         % ����߶�
N  = rE / sqrt( 1 - eE2*sin( B_bj )^2 );
% -------------��ص�������ϵ&&WGS����ϵ-------------------
x_bj = (N+hX)*cos(B_bj)*cos(L_bj);
y_bj = (N+hX)*cos(B_bj)*sin(L_bj);
z_bj = (N*(1-eE2)+hX)*sin(B_bj);
rXGZ = [x_bj;y_bj;z_bj];     % �ӵ��ĵ��Ź�վ��ʸ��
% ---------------��������ϵ&&WGS����ϵ--------------------------
wgs_surface = [cos(L_bj)*sin(B_bj)   sin(L_bj)*sin(B_bj)    -cos(B_bj);
                     -sin(L_bj)               cos(L_bj)           0   ;
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
% -----------------�������-----------------------------
TimeEnd  = 13;                  % ����û������
tEnd     = 6532;                % һȦ������ʱ��
NumStall = 18;                  % һ��18�������

%
%% load data of satellite
load('ds/dataTrac');

%% 
% ******************�ն˻���������һЩ����***************************
vUE = [0;0;277];                        % NED����ϵ
[ BdegUe,LdegUe ] = XYZsettoBLH(rXGZ);  % UE��γ�ȣ�����
vUEwgs = wgs_surface*vUE;

Berr = 5;      Lerr = 5;              % �ն������ǿ�������ľ�γ�Ȳ�
% ******************ѡ��****************************
tspan = 1;
time = 0:tspan:tEnd;
tAcess = [];
LenSet = [];
for inxTim = 1:1:1000       % �����ĳһ��ʱ��
    t = time(inxTim);
    % -------- ѡ�� 1 ---------------------------------
    %         cnt = 1;
    %         for faceinx = 1:18
    %             for trainx = 1:48
    %                 rStallWgs = satPosECEF(:,inxTim,trainx,faceinx);     % 18�����ǵ�λ�ã���һ�������
    %                 vStallWgs = vECEF(:,inxTim,trainx,faceinx);          % 18�����ǵ��ٶȣ���һ�������
    %                 [ BdegSta(trainx , faceinx) , LdegSta(trainx , faceinx) ] = XYZsettoBLH( rStallWgs );
    %                 if  (abs(BdegSta(trainx , faceinx) - BdegUe)<Berr ) && (abs(LdegSta(trainx , faceinx) - LdegUe)<Lerr)
    %                     stallSettp(cnt,:) = [ trainx , faceinx ];
    %                     cnt = cnt + 1;
    %                 end
    %             end
    %         end
    % ---------- ѡ�� 2 ---------------------------
    for faceinx = 1:18
        for trainx = 1:48
            rStallWgstp = satPosECEF(:,inxTim,trainx,faceinx);
            alphaRef( trainx , faceinx ) =  rXGZ'*rStallWgstp/( norm(rXGZ)*norm(rStallWgstp) );
        end
    end
    % ѡ��
    stallSet = intersect( find( alphaRef >=cos( alphaEar ) ) , find( alphaRef>0 ) );        % ͬʱ������������
    vReftp = [];
    if isempty(stallSet)
        fprintf('û�п���ѡ�����\n');
        stallSet = [];
    else
        stallSettp = [];
        stallSettp(:,2) = floor((stallSet -0.1)/48)+1;                      % 48��ѡ1��
        stallSettp(:,1) = stallSet - ( stallSettp(:,2)-1 ).*48;             % 18��ѡ1��
        for index = 1:length(stallSet)
            vStall      =      vECEF(:,inxTim,stallSettp(index,1) ,stallSettp(index,2) );   %
            rStallWgstp = satPosECEF(:,inxTim,stallSettp(index,1) ,stallSettp(index,2) );   %
            r_xinToStal = rXGZ - rStallWgstp;
            r_xinToStal = r_xinToStal/norm(r_xinToStal);
            r_StalToXin =  - r_xinToStal;                                   % ʸ��: �Ź�վָ������
            vReftp( index ) = (vStall - vUEwgs)'*r_xinToStal;               % ����ٶȣ��������������Ƶƫ
        end
        inxReftp = find(vReftp>0);
        if isempty(inxReftp)
            fprintf(' ����ٶȣ�û�д���0�� ');
        end
        statelSect = stallSettp( find(vReftp == max( vReftp ) ), : );       % ������һ��
        fprintf('%3.0f\n',length( stallSet ));
    end
    LenSet = [LenSet, length( stallSet )];
    
end
fprintf('������Χ%4.0fkm\n',alphaEar*2*rE/1e3)
[f,xi] = ksdensity(LenSet);
figure;plot(xi,f);xlabel('��ѡ���ǵĿ���');ylabel('����');

% ********************������Χ����************************
% rbeam = alphaEar*rE;
% poinum = 17;
% sitaPoi = 0:pi/8:2*pi;
% TpSet8PoitIQ = exp(1j*sitaPoi);
% TpSet8Poit = [real(exp(1j*sitaPoi));imag(exp(1j*sitaPoi));zeros(1,17)].*rbeam;      % �ڵر�����һ��Ȧ���ڵ����������ϵ��
% 
% inxTim = 100;
% cout = 1;
% for faceinx = 1:18
%     for trainx = 1:48
%         rStallWgstp = satPosECEF(:,inxTim,trainx,faceinx);
%         alphaRef( trainx , faceinx ) =  rXGZ'*rStallWgstp/( norm(rXGZ)*norm(rStallWgstp) );
%         % ����ϵת��
%         [ BstateTp , LstateTp , ~ ] = XYZsettoBLH(rStallWgstp,1);
%         LstateTp = LstateTp-pi;         % ��֪��Ϊʲô
%         wgs_surface = [cos(LstateTp)*sin(BstateTp)   sin(LstateTp)*sin(BstateTp)    -cos(BstateTp);
%                                    -sin(LstateTp)               cos(LstateTp)           0   ;
%                      cos(LstateTp)*cos(BstateTp)   sin(LstateTp)*cos(BstateTp)    sin(BstateTp)]';
%         wgs_surface = wgs_surface*diag([-1,1,1]);
%         % ת��WGS����ϵ��
%         x_bj = (N + 0) *cos(BstateTp)*cos(LstateTp);
%         y_bj = (N + 0) *cos(BstateTp)*sin(LstateTp);
%         z_bj = (N*(1-eE2)+0)*sin(BstateTp);
%         rStallWgstpSur = [x_bj;y_bj;z_bj];
%         Point8Wgs = wgs_surface*TpSet8Poit + repmat(rStallWgstpSur,1,poinum);
%         poiTpPlt{cout} = Point8Wgs;
%         cout = cout + 1;
%     end
% end

% -----------��ͼ 1.1�����ǲ������Ƿ�Χ ----------------------------
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
% xlabel('γ��');ylabel('����');title('���ǲ�������ʾ��ͼ')

% -----------��ͼ 1.2��������滭��Ȧ (����û�кܴ������) ----------------------------
% figure;
% rE     = 6378.137*1e3;      % ����뾶
% [Xe,Ye,Ze] = sphere(24);    Xe = rE*Xe;     Ye = rE*Ye;     Ze = rE*Ze;     cla reset;   load topo;
% s = surface(Xe,Ye,Ze,'FaceColor','texturemap','CData',topo); axis equal;colormap(topomap1);brighten(0.9);hold on;
% plot3(poiTpPlt{1}(1,:),poiTpPlt{1}(2,:),poiTpPlt{1}(3,:),'r','LineWidth',1.4);  hold on;
% rStallWgstp = satPosECEF(:,inxTim,1,1);
% plot3(rStallWgstp(1),rStallWgstp(2),rStallWgstp(3),'b*');hold on;
% plot3(rStallWgstpSur(1),rStallWgstpSur(2),rStallWgstpSur(3),'k*')

%% ����ѡ�� ����
% % --------------14.1 ����ѡ��--------------------------------
% alpha = 80/180*pi;      % ���ߵĽǶȣ�80���Ƿǳ���
% % --------------14 �������ٶȼ��� ---------------------------------
% r_xinToStal = rXGZ - satPosECEF(:,inx);
% r_xinToStal = r_xinToStal/norm(r_xinToStal);
% r_StalToXin =  - r_xinToStal;       % ʸ�����Ź�վָ������
% RefCosAlpha = r_StalToXin'*antDirecWgs / ( norm(r_StalToXin)*norm(antDirecWgs) );     % �������߷���� theta
% 
% if ( RefCosAlpha >= cos(alpha) ) && RefCosAlpha>0      % �����ߵĲ�����Χ֮��
%     alphaRef = acos( rXGZ'*satPosECEF(:,inx)/( norm(rXGZ)*norm(satPosECEF(:,inx)) ) );      % ���Ǻ��ն���Ե��ĵĽ�
%     if alphaRef <= alphaEar
%         vUE = [0;0;277];                    % ���ڳ�����˵��ֻ��xy�����ٶȣ��ڵ����������ϵ
%         vUEwgs = wgs_surface*vUE;
%         vRef = ( vECEF(:,inx) + vUEwgs )'*r_xinToStal;
%         fs = 19.3e9;                        % �ز�Ƶ��
%         c = 3e8;                            % ����
%         delta_f(inx) = vRef/c*fs;           % ����Ƶƫ
%         sita(inx) = acos(RefCosAlpha);      % �������߷���Ǵ�С
%     end
% else                                        % �������ߵĲ�����Χ֮��
%     % sita(inx) = acos(RefCosAlpha);        % �������߷���Ǵ�С
%     sita(inx) = 0;
% end

%% ��ͼ
% figure;plot(Vkgot);

% rECILast = satPosECEF(:,end);       % ���յ�λ��
% sita_dot = sita(61:end)-sita(1:end-60);
% % -----------��ͼ 0������ָ��ı仯�� ------------------------
% figure;plot(0:tspan:length(sita_dot)-1 , sita_dot./pi*180 ,'r.');grid on;xlabel('ʱ��t s');ylabel('\theta�仯�� ����/min');
% % -----------��ͼ 1��Ƶƫ------------------------
% figure;plot(0:tspan:(0+tEnd),sita/pi*180);
% if ~isempty(delta_f);   figure;plot(0:tspan:length(delta_f)-1,delta_f);grid on;xlabel('ʱ�� t s');ylabel('Ƶƫ \deltaf Hz');    end
% ----------- ��ͼ 2������----------------------------
% figure;
% rE     = 6378.137*1e3;      % ����뾶
% [Xe,Ye,Ze] = sphere(24);    Xe = rE*Xe;     Ye = rE*Ye;     Ze = rE*Ze;     cla reset;   load topo;
% s = surface(Xe,Ye,Ze,'FaceColor','texturemap','CData',topo); axis equal;colormap(topomap1);brighten(0.9);hold on;
% for i = 1:18
%     plot3(satPosECEF(1,:,1,i),satPosECEF(2,:,1,i),satPosECEF(3,:,1,i),'r--','LineWidth',0.1);hold on;
% end
% plot3(rXGZ(1),rXGZ(2),rXGZ(3),'*');grid on;hold on;
% plot3([0,1.2*rXGZ(1)],[0,1.2*rXGZ(2)],[0,1.2*rXGZ(3)],'LineWidth',2.5,'Color','b'); 

% -----------��ͼ 3�����е����Ƕ�̬���˶��켣 (��ͼ) ----------------------------
% figure;
% rE  = 6378.137*1e3;      % ����뾶
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

% -----------��ͼ 4��ѡ������ ----------------------------
% figure;
% rE  = 6378.137*1e3;      % ����뾶
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
%     plot3([ rStallWgstp(1)] , [rStallWgstp(2)], [ rStallWgstp(3)] ,'r*');           % ��ѡ������Ǳ�ǳ���
%     hold on;
% end
% rStallWgstp = satPosECEF(:,inxTim,statelSect(1,1) ,statelSect(1,2) );
% plot3( rStallWgstp(1),rStallWgstp(2),rStallWgstp(3) , 'ko' )

% -----------��ͼ 5��������Χ ----------------------------




