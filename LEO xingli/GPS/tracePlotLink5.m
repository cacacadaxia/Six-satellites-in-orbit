
% =========================================================================
%
%                      ��������λ���ٶȽ���
%
%
% =========================================================================
%
%��(C)2019-2020 ���ݺ���ͨ�����޹�˾
%   �汾��V1.5
%   ���ڣ�2019��7��25��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ����:  1.���Ƶ�Link4
%        2.����������һЩ���ǲ������бȽϣ������ǵ͹�����ǹ����������֮ǰ�кܴ�ĳ���
%        3.������Ƶƫ
%        4.������֤��ģ������
%        5.
%        6.
%--------------------------------------------------------------------------

clear all;  close all;
%---------------------��ȡRINEX��ʽn�ļ�--------------------------------
data = RinexNreader('BR.16n','G01'); % ע���ȡ·�������Ǳ��
%---------------------��������յ��ܻ���---------------------------------
% [JD,FOD,GPSW,SOW,DOY,DOW] = GCtoGPS(data(1,1),data(2,1),data(3,1),data(4,1),data(5,1),data(6,1));
%----------����n�ļ��вο���Ԫ���ܻ���----------
for j = 1:size(data,2) % size(data��2) - ����data���������;
    [JD,FOD,GPSW,SOW(j),DOY,DOW] = GCtoGPS(data(1,j),data(2,j),data(3,j),data(4,j),data(5,j),data(6,j));
    % ��ָ�����ǵ� �ꡢ�¡��ա�ʱ���֡��� �� ������(��������)��������(С������)��GPS�ܡ��ܻ��롢����ա�������
end

% �������
GM     = 3986005e8; % ��=GM
we     = 7.2921151467e-5; % ������ת���ٶ�
eE2    = 0.00669437999013;  % ���������ʵ�ƽ��
eE     = 0.0818191910428;   % ������
rE     = 6378.137*1e3;      % ����볤��
J2     = 1082.63e-6;        % �����е�ģ��

% -------------����-ģ���Ź�վ��λ�� (����120����γ37)---------------------------------
B_bj = (39+48/60)*pi/180;       % γ��
L_bj = (0 +28/60 - 180)*pi/180;      % ����
hX = 0;                         % ����߶�
N  = rE / sqrt( 1 - eE2*sin( B_bj )^2 );
% -------------����ϵת�����������ϵ����-------------------
x_beijing = (N+hX)*cos(B_bj)*cos(L_bj);
y_beijing = (N+hX)*cos(B_bj)*sin(L_bj);
z_beijing = (N*(1-eE2)+hX)*sin(B_bj);
rXGZ = [x_beijing;y_beijing;z_beijing];     % �ӵ��ĵ��Ź�վ��ʸ��
%-----------------------���Ǻ���Ҫ�Ĳ���------------------------------------
inx = 1;
tspan = 1;
dotM2 = 0;
dotOm = 0;
delta_f = [];

%% 
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
sqrta = sqrt(7229175.73885051);
e     = 0.000394000000000000;
om    = 0.982207577030007;
Om0   = 5.41244749270632;
i0    = 1.23689040627885+1;
doti  = 0;
M0    = 0;
delta_n = 0;            % ��С�����Ժ��Բ���

% -----------�������----------------------
TimeEnd = 13;
tEnd = 40000;
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
    %----------11. ����������WGS-84����ϵ��λ��----------
    xECEF = xECI * cos(lamdak) - yECI * cos(ik) * sin(lamdak);
    yECEF = xECI * sin(lamdak) + yECI * cos(ik) * cos(lamdak);
    zECEF = yECI * sin(ik);
    %----------11. ���λ��-------------------------
    satPosECEF(:,inx) = [xECEF;yECEF;zECEF];    % ��WGS-84����ϵ�µ�λ�ã���Ϊ�е�����ת�Ĵ��ڣ����Կ������ȽϹ�
    rECI(:,inx)       = [xECI;yECI;0];          % �������ϵ�µ�λ��
    
    %% ---------12. �����ٶ�--------------------------
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

    % ------------12.4 ת������ϵ --------------------
    VxECEF = xk_dot*cos(lamdak) - yk_dot*cos(dotik)*sin(lamdak) + dotik*zECEF*sin(lamdak) - dotLk*yECEF;
    VyECEF = xk_dot*sin(lamdak) + yk_dot*cos(dotik)*cos(lamdak) - dotik*zECEF*cos(lamdak) + dotLk*xECEF;
    VzECEF = yk_dot*sin(ik) + dotik*yECI*cos(dotik);
    vECEF(:,inx)  = [VxECEF;VyECEF;VzECEF];          % ��WGS-84����ϵ�µ����ǵ��ٶ�

    % ��һ�� (����ģ�������ϵ��أ����Ǽ������������ȷ�ģ�Ŀǰ����Ҳ���Ǻ���ȷ��)
    %     r1 = [yECI/xECI;-1;0];
    %     r2 = r1/norm(r1);
    %     vECI(:,inx) = 7425*r2;      % 7425Լ����ƽ���ٶ�
    %     vECEF(1,inx) = vECI(1,inx) * cos(lamdak) - vECI(2,inx) * cos(ik) * sin(lamdak);
    %     vECEF(2,inx) = vECI(1,inx) * sin(lamdak) + vECI(2,inx) * cos(ik) * cos(lamdak);
    %     vECEF(3,inx) = vECI(2,inx) * sin(ik);
    %% ------------13. �㶯���� -------------------------
    K = n*J2/2*( rE/(sqrta^2*(1-e^2)) )^2;
    dotOm2 = -3*K*cos(ik);
    dotOm  = 0 + dotOm2 ;               % ʹ���㶯���������
    dotom2 = 3*K*(2 - 5/2*sin(ik)^2);
    om = om + dotom2 * tspan;
    dotM2   = 3/4*n*J2*( 3*cos(ik)^2-1 )/ (1-e^2)^(3/2) * (rE/sqrta^2)^2;
    % --------------14.1 ����ѡ��--------------------------------
    alpha = 80/180*pi;      % ���ߵĽǶ�
    % --------------14 �������ٶȼ��� ---------------------------------
    r_xinToStal = rXGZ - satPosECEF(:,inx);
    r_xinToStal = r_xinToStal/norm(r_xinToStal);
    r_StalToXin =  - r_xinToStal;       % ʸ�����Ź�վָ������
    RefCosAlpha = r_StalToXin'*rXGZ / ( norm(r_StalToXin)*norm(rXGZ) );
    
    if ( RefCosAlpha >= cos(alpha) ) && RefCosAlpha>0      % �����ߵĲ�����Χ֮��
        vUE = [100;100;0];
        vRef = vECEF(:,inx)'*r_StalToXin;
        fs = 19.3e9;   %  �ز�Ƶ��
        c = 3e8;       %  ����
        delta_f(inx) = vRef/c*fs;
        sita(inx) = acos(RefCosAlpha);      % �������Ǵ�С
    else                                    % �������ߵĲ�����Χ֮��
        sita(inx) = acos(RefCosAlpha);      % �������Ǵ�С
    end

    % ����i
    inx = inx + 1;
    % others
    lamdakSave(inx) = lamdak;   
end
rECILast = satPosECEF(:,end);
% -----------��ͼ 1��Ƶƫ------------------------
figure;plot(0:tspan:(0+tEnd),sita/pi*180);
if ~isempty(delta_f);   figure;plot(0:tspan:length(delta_f)-1,delta_f);grid on;xlabel('ʱ�� t s');ylabel('Ƶƫ \deltaf Hz');    end
% -----------��ͼ 2������----------------------------
rE     = 6378.137*1e3;      % ����뾶
r_bar = satPosECEF;
figure;plot3(r_bar(1,:),r_bar(2,:),r_bar(3,:),'r','LineWidth',1);
hold on;plot3(rXGZ(1),rXGZ(2),rXGZ(3),'*');grid on;hold on;
[Xe,Ye,Ze] = sphere(24);
Xe = rE*Xe;     Ye = rE*Ye;     Ze = rE*Ze;
surf(Xe,Ye,Ze); axis equal;     hold on;
plot3([0,1.8*rXGZ(1)],[0,1.8*rXGZ(2)],[0,1.8*rXGZ(3)],'LineWidth',2.5,'Color','b');

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




