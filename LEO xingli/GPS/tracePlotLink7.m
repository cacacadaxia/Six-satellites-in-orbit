
% =========================================================================
%
%                      ��������λ���ٶȽ���
%
%
% =========================================================================
%
%��(C)2019-2020 ���ݺ���ͨ�����޹�˾
%   �汾��V1.7
%   ���ڣ�2019��7��31��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ����:  1.���Ƶ�Link6
%        2.������һ��������Χ����Ҫ����
%        3.���������
%        4.
%        5.
%        6.
%--------------------------------------------------------------------------

clear all;  close all;
rng('default');
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
GM     = 3986005e8; % ��=GM
we     = 7.2921151467e-5; % ������ת���ٶ�
eE2    = 0.00669437999013;  % ���������ʵ�ƽ��
eE     = 0.0818191910428;   % ������
rE     = 6378.137*1e3;      % ����볤��
J2     = 1082.63e-6;        % �����е�ģ��
% ------------------���ǲ���-----------------------------
alphaEar = 1000e3 /rE/2;
r3tp  = sqrt(rE^2 + (rE+1175e3)^2 - 2* cos(alphaEar) *rE*(rE+1175e3) );
betaState = asin( sin(alphaEar) *rE/r3tp);            % �ò���


% -------------����-ģ���Ź�վ��λ�� (����120����γ37)---------------------------------
B_bj = (25+48/60)*pi/180;            % γ��
L_bj = (-64+28/60 - 180)*pi/180;      % ����
hX = 0;                         % ����߶�
N  = rE / sqrt( 1 - eE2*sin( B_bj )^2 );
% -------------��ص�������ϵ&&WGS����ϵ-------------------
x_bj = (N+hX)*cos(B_bj)*cos(L_bj);
y_bj = (N+hX)*cos(B_bj)*sin(L_bj);
z_bj = (N*(1-eE2)+hX)*sin(B_bj);
rXGZ = [x_bj;y_bj;z_bj];     % �ӵ��ĵ��Ź�վ��ʸ��
% ---------------��������ϵ&&WGS����ϵ--------------------------
wgs_surface = [cos(L_bj)*sin(B_bj)   sin(L_bj)*sin(B_bj)    -cos(B_bj);
                -sin(L_bj)               cos(L_bj)           0        ;
               cos(L_bj)*cos(B_bj)   sin(L_bj)*cos(B_bj)    sin(B_bj)]';
wgs_surface = wgs_surface*diag([-1,1,1]);           % ���춫����ϵ-->���춫����ϵ
% ----------------��������ϵ&&��������ϵ-----------------------
theta = 0 /180*pi ;          % ���� ��ʱ����
fai   = 0 /180*pi ;          % ƫ��
gamma = 45 /180*pi ;          % ����
ben_surface = [cos(theta)*cos(fai)                                        sin(theta)             -cos(theta)*sin(fai);
                -sin(theta)*cos(fai)*cos(gamma)+sin(fai)*sin(gamma)     cos(theta)*cos(gamma)    sin(theta)*sin(fai)*cos(gamma)+cos(fai)*sin(gamma)
                sin(theta)*cos(fai)*sin(gamma)                          -cos(theta)*sin(gamma)   -sin(theta)*sin(fai)*sin(gamma)+cos(fai)*cos(gamma)];
wgs_ben     = wgs_surface*ben_surface';
antDirec    = [0; 0; 1];      % �����ڱ�������ϵ�еķ���
antDirecWgs = wgs_ben*antDirec;
%-----------------------���Ǻ���Ҫ�Ĳ���------------------------------------
inx = 1;
tspan = 1;
dotM2 = 0;
dotOm = 0;
delta_f = [];
inxAnt = 1;

%% ���ǹ������
index   = 1;
% delta_n = data(12,index);       % detn ƽ����ǵĳ��ڱ仯(���ص����)
% sqrta   = data(17,index);       % sqrta �������ƽ����
toe       = data(18,index);       % t0e �����ο�ʱ��(�ܻ���)
% M0      = data(13,index);       % M0 �ο�ʱ�̵�ƽ�����
% e       = data(15,index);       % e ƫ����
% om      = data(24,index);       % omega ���ص�Ǿ�
% Om0     = data(20,index);       % OMEGA0 �ǲο�ʱ��toe������ྭ����ʼ�ڸ�����������Ȧ�����ǹ����������ν׼����
% % dotOm   = data(25,index);       % OMEGADOT ������ྭ�ڳ��ƽ���еĳ��ڱ仯
% i0      = data(22,index);       % i0 �ο�ʱ�̹�����
% doti    = data(26,index);       % IDOT �����Ǳ仯��
    % p.s.
% 86400         һ��
% 4.3077e+04    һ������
% -------------------------------------------
sqrta = sqrt(1175e3 + rE);      % Ŀǰ�͹�Ĳ���
e     = 0.000394000000000000;
om    = 0.982207577030007;
Om0   = 5.41244749270632;
i0    = 86.7/180*pi ;
doti  = 0;
M0    = 0;
delta_n = 0;            % ��С�����Ժ��Բ���
% -----------�������----------------------
TimeEnd = 13;
tEnd = 3600;
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
    Vk = atan2(Vks,Vkc);
    %----------6. ���������Ǿ�----------
    phk = Vk + om;
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
    satPosECEF(:,inx) = [xECEF;yECEF;zECEF];    % ��WGS-84����ϵ�µ�λ�ã���Ϊ�е�����ת�Ĵ��ڣ����Կ������ȽϹ�
    rECI(:,inx)       = [xECI;yECI;0];          % �������ϵ�µ�λ��
    
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
    vECI(:,inx) = [xk_dot; yk_dot; 0];

    % ------------12.4 ת������ϵ--WGS����ϵ --------------------------
    VxECEF = xk_dot*cos(lamdak) - yk_dot*cos(dotik)*sin(lamdak) + dotik*zECEF*sin(lamdak) - dotLk*yECEF;
    VyECEF = xk_dot*sin(lamdak) + yk_dot*cos(dotik)*cos(lamdak) - dotik*zECEF*cos(lamdak) + dotLk*xECEF;
    VzECEF = yk_dot*sin(ik) + dotik*yECI*cos(dotik);
    vECEF(:,inx)  = [VxECEF;VyECEF;VzECEF];          % ��WGS-84����ϵ�µ����ǵ��ٶ�

    %% ------------13. �㶯���� -------------------------
    K = n*J2/2*( rE/(sqrta^2*(1-e^2)) )^2;
    dotOm2 = -3*K*cos(ik);              
    dotOm  = 0 + dotOm2 ;               % ʹ���㶯���������
    dotom2 = 3*K*(2 - 5/2*sin(ik)^2);
    om     = om + dotom2 * tspan;
    dotM2  = 3/4*n*J2*( 3*cos(ik)^2-1 )/ (1-e^2)^(3/2) * (rE/sqrta^2)^2;
    % --------------14.1 ����ѡ��--------------------------------
    alpha = 80/180*pi;      % ���ߵĽǶȣ�80���Ƿǳ���
    % --------------14 �������ٶȼ��� ---------------------------------
    r_xinToStal = rXGZ - satPosECEF(:,inx);
    r_xinToStal = r_xinToStal/norm(r_xinToStal);
    r_StalToXin =  - r_xinToStal;       % ʸ�����Ź�վָ������
    RefCosAlpha = r_StalToXin'*antDirecWgs / ( norm(r_StalToXin)*norm(antDirecWgs) );     % �������߷���� theta
    
    if ( RefCosAlpha >= cos(alpha) ) && RefCosAlpha>0      % �����ߵĲ�����Χ֮��
        alphaRef = acos( rXGZ'*satPosECEF(:,inx)/( norm(rXGZ)*norm(satPosECEF(:,inx)) ) );      % ���Ǻ��ն���Ե��ĵĽ�
        if alphaRef <= alphaEar 
            inxAnt = inxAnt + 1;
            % ���߷�λ�Ǳ仯
            rECIAnt = satPosECEF(:,inx) - rXGZ;
            
            Rb(:,inxAnt) = wgs_ben'*rECIAnt;
            lam(inxAnt) = asin( Rb(3,inxAnt)/norm(rECIAnt) );
            vflow(inxAnt) = atan( Rb(1,inxAnt)/Rb(2,inxAnt) );
            
            % ����Ƶƫ
            vUE = [0;0;277];                    % ���ڳ�����˵��ֻ��xy�����ٶȣ��ڵ����������ϵ
            vUEwgs = wgs_surface*vUE;
            vRef = ( vECEF(:,inx)  -  vUEwgs )'*r_xinToStal;        % ���������ٶ��Ǽ�
            fs = 19.3e9;                        % �ز�Ƶ��
            c = 3e8;                            % ����
            delta_f(inx) = vRef/c*fs;           % ����Ƶƫ
            sita(inx) = acos(RefCosAlpha);      % �������߷���Ǵ�С
        end
        alphaRef_sav(inx) = alphaRef;
    else                                        % �������ߵĲ�����Χ֮��
        % sita(inx) = acos(RefCosAlpha);        % �������߷���Ǵ�С
        sita(inx) = 0;
    end

    % ����i
    inx = inx + 1;
    % others
    lamdakSave(inx) = lamdak;   
end
rECILast = satPosECEF(:,end);       % ���յ�λ��
sita_dot = sita(61:end)-sita(1:end-60);

% figure;plot(alphaRef_sav/pi*180);
% -----------��ͼ 0.1������ָ��仯 ------------------------
figure;plot(1:length(lam),lam/pi*180,'r');hold on
plot(1:length(lam) , vflow/pi*180,'k');
legend('\lambda','\nu');xlabel('t s');ylabel('rad');grid on
% -----------��ͼ 0������ָ��ı仯�� ------------------------
figure;plot(0:tspan:length(sita_dot)-1 , sita_dot./pi*180 ,'r.');grid on;xlabel('ʱ��t s');ylabel('\theta�仯�� ����/min');
% -----------��ͼ 1��Ƶƫ------------------------
figure;plot(0:tspan:(0+tEnd),sita/pi*180);
if ~isempty(delta_f);   figure;plot(0:tspan:length(delta_f)-1,delta_f);grid on;xlabel('ʱ�� t s');ylabel('Ƶƫ \deltaf Hz');    end
% -----------��ͼ 2������----------------------------
figure;
rE     = 6378.137*1e3;      % ����뾶
r_bar = satPosECEF;     
[Xe,Ye,Ze] = sphere(24);
Xe = rE*Xe;     Ye = rE*Ye;     Ze = rE*Ze;
cla reset;   load topo;
s = surface(Xe,Ye,Ze,'FaceColor','texturemap','CData',topo); axis equal;     
view(0,100);
colormap(topomap1);
brighten(0.8);
hold on;
plot3([0,1.4*rXGZ(1)],[0,1.4*rXGZ(2)],[0,1.4*rXGZ(3)],'LineWidth',2.5,'Color','b');hold on
plot3(r_bar(1,:),r_bar(2,:),r_bar(3,:),'r','LineWidth',1);hold on;
plot3(rXGZ(1),rXGZ(2),rXGZ(3),'*');grid on;hold on;
% ------���߽Ƕ�-------
% antDirecWgs
rFin = rXGZ + antDirecWgs*5000e3;
plot3([rXGZ(1) , rFin(1) ] , [rXGZ(2) ,rFin(2) ] , [rXGZ(3) , rFin(3) ] ,'k','LineWidth',2.5)

    % hold on;
    % v1 = vUE;
    % v2 = wgs_surface*v1*1e5;
    % vtp = v2  + rXGZ;
    % plot3([rXGZ(1),vtp(1)],[rXGZ(2),vtp(2)],[rXGZ(3),vtp(3)],'LineWidth',2.5,'Color','b');      % �����ٶȵķ���
% -----------��ͼ 3���ٶ� ----------------------------
% figure;
% plot([0:tspan:(0+tEnd)]/3600,vECEF(1,:));hold on;
% plot([0:tspan:(0+tEnd)]/3600,vECEF(2,:));hold on;
% plot([0:tspan:(0+tEnd)]/3600,vECEF(3,:));
%-----------��ͼ 4������γ�ȹ켣 ----------------------------
% for i = 1:length(satPosECEF)
%     [B_ECEF(i),L_ECEF(i)] = XYZtoBLH(satPosECEF(1,i),satPosECEF(2,i),satPosECEF(3,i));
% end
% figure;
% plot(L_ECEF,B_ECEF,'r.'); % ��������
% grid on;
% title('���Ƿ��й켣');
% xlabel('���� rad');
% ylabel('γ�� rad');
% hold on;
% --------------��ͼ 5�� -----------------------------

%% ����Ա�����
% save('ds/datasave','delta_f');
