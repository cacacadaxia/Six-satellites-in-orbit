

 % =========================================================================
%
%                      ����λ���ٶȽ���
% 
%
% =========================================================================
%
%��(C)2019-2020 ���ݺ���ͨ�����޹�˾
%   �汾��V1.2
%   ���ڣ�2019��7��23��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ����:  1.���Ƶ�1�ģ����ǵ��������㶯
%        2. 
%        3. 
%        4. 
%        5.
%        6.
%--------------------------------------------------------------------------
clear all;
close all;
%% ������
% ----�����һ����������-------------
[oe,epoch,yr,M,E,satname] = TLE2oe('xingli.txt');    % ����������������
a = oe(1);      % a �볤��
e = oe(2);      % e ƫ����
i = oe(3);      % i ���ƽ�����
% i = 1.0996;
% i =  0.9425;
Om0 = oe(4);     % Om ������ྭ
om = oe(5);     % om �����Ǿ࣬���ص�Ǿ࣬���߽н��ص����
nu = oe(6);     % nu ʱ�̣���λ��s��
% ʵ���������и�������ǣ����ﲢ����ʱ�̣����ǽ��ص�ʱ�̣�
% Ӧ����������� sita ֵ��һ����

%% % ---GPS������--------------------------------
% a = 0.515373216438e4^2;
% e = 0.202130551916e-1;
% i = 0.93408487355;
% Om = -0.586481950333;
% om = -0.211637251436e1;

%%



%%

% ------����----------------------------
aG = 40/180*pi;         % �������κ���ʱ��
tw = 0;                 % ����Ĳ���
GASTw = aG;             % �������κ���ʱ��
toe = 0;                % ���ǽ��ص�ʱ�� (��λ��s)
                        % toe �����ο�ʱ��(�ܻ���)
t  = toe;                % none����ʼ��ֵ
tspan = 1;              % ʱ�䲽��
coutmax =  floor(4.3080e+05); 
coutmax = 86400;
% -------------�Ź�վ����γ������ϵ�£�-------------------------------------
Hx = 6371e3;        % �߶�
Bx = 40/180*pi;     % γ��
Lx = 60/180*pi;     % ����

% ------------����-----------------------
eE2    = 0.00669437999013;  % ���������ʵ�ƽ��
eE     = 0.0818191910428;   % ������  
rE     = 6378.137*1e3;      % ����뾶
owE    = 7.29211567e-5;     % ������ת�Ľ��ٶ�
% owE = 0;
% J2     = -0.4841667798797018e-3;    % CGCS2000��г����ģ��ϵ��(Ӧ�ò�׼ȷ)
J2     = 1082.63e-6;        % �����е�ģ��

% ------��γ������ϵ-->�Ź�վ�������ϵ-------------------
N  = rE / (sqrt(1-eE2*sin(Bx)^2));
Xx = (N+Hx)*cos(Bx)*cos(Lx);
Yx = (N+Hx)*cos(Bx)*sin(Lx);
Zx = (N*(1-eE2)+Hx)*sin(Bx);
rX = [Xx;Yx;Zx];            % �Ź�վ��λ��
% -------------------------------------
dotOm = 0;
for ii = 1:coutmax+1
    miu = 3.986004415e14;   % ������
    n0 = sqrt(miu/(a)^3);   % ���ٶ�
    t  = t + tspan;
    tk = t - toe;
    %----------2. ����黯ʱ��tk----------
    if tk > 302400
        tk = tk - 604800;
    elseif tk < -302400
        tk = tk + 604800;
    end
    % ���� M
    if ii ==1
        M = n0*tk;
    else
        M = M + dotM*tspan;          % ����ƽ����ǣ�Ӧ���ǲ���ȥ���µģ� 
    end

    % GAST = GASTw + (dotOm - owE)*tk - owE*toe;      % ���¸���ʱ��
    
    % ---------------------------------------------
    sigema = 0.0001;      % ����
    E0 = 1.2;           % ������õĳ�ֵ�����㿪���շ���
    for inx = 1:1000
        E1 = M + e*sin(E0);
        Etp(inx) = E1;
        if abs(E1 - E0)<sigema
            E = E1;
            break;
        end
        E0 = E1;
    end
    r = a*(1-e*cos(E));
    sita = 2*atan( (sqrt((1+e)/(1-e)))*tan(E/2) );  
    uk = sita + om;         % ���������������ȡλ��ʱ��ʱ�����ǣ��൱�ڽ�om�ı仯�Ѿ����ǽ�ȥ��
    
    % �����е��е㲻һ����ע��һ��
    x_a = r*cos(sita);      % �������ϵ�µ�����λ�ã�ԭ���ǵ�������
    y_a = r*sin(sita);
    z_a = 0;
    if ii ==1
        savedata1 = [x_a;y_a];
    elseif ii == coutmax
        savedata2 = [x_a;y_a];
    end
    
    % -------�ٶ�------------------------
    Ekdot  = n0/(1-e*sin(E));
    Faidot = sqrt((1+e)/(1-e))*(cos(sita/2)^2)/(cos(E/2)^2)*Ekdot;
    ukdot  = Faidot;
    rkdot  = Ekdot*a*e*sin(E);
    vx = rkdot*cos(uk) - ukdot*y_a;
    vy = rkdot*sin(uk) + ukdot*x_a;
    vz = 0;
    
    % ----------λ��ת�����������ϵ-----------------------------------
    lamdak = Om0 + (dotOm - owE)*tk - owE*toe;
    
    
    R1 = [cos(om),-sin(om),0 ; sin(om),cos(om),0 ; 0 ,0 , 1];
    R2 = [1 ,0 ,0 ; 0 ,cos(i),-sin(i) ; 0 ,sin(i),cos(i)];
    R3 = [cos(lamdak),-sin(lamdak),0 ; sin(lamdak),cos(lamdak),0 ; 0 ,0 ,1];    % Lk
    r_bar(:,ii) = R1*R2*R3*[ x_a; y_a; z_a ];     % �����ڴ������ϵ��Ӧ�þ���WGS-84����ϵ
    v_bar(:,ii) = R2*R3*[ vx ; vy ; vz ];         % �����ڴ������ϵ

    % -----------------����----------------------------------
    K = n0*J2/2*( rE/(a*(1-e^2)) )^2;
    dotOm = -3*K*cos(i);
    Om0 = Om0 + dotOm * tspan;
    dotom = 3*K*(2 - 5/2*sin(i)^2);
    om = om + dotom * tspan;
    dotM   = n0 + 3/4*n0*J2*( 3*cos(i)^2-1 )/ (1-e^2)^(3/2) * (rE/a)^2;

    lamdakSave(ii) = lamdak;
end
lastP = r_bar(:,end);
% ------------------------------------

fprintf('��ʼʱ�̣�%3.0f , ����ʱ�̣�%3.0f , \n',toe,t);
% -----------��ͼ----------------------------
figure(1);plot3(r_bar(1,:),r_bar(2,:),r_bar(3,:),'r','LineWidth',1);
hold on;plot3(0,0,0,'*');grid on;
R_satellite = r_bar - rX;               % ����-->�Ź�վ ʸ��


%%
hold on;
[Xe,Ye,Ze] = sphere(24);
Xe = rE*Xe;
Ye = rE*Ye;
Ze = rE*Ze;
surf(Xe,Ye,Ze)
axis equal;
% -------------����-ģ���Ź�վ��λ��---------------------------------
longitude_beijing = (39+48/60)*pi/180; 
latitude_beijing = (116+28/60)*pi/180;
% -------------����ϵת�����������ϵ����-------------------
x_beijing=rE*cos(longitude_beijing)*cos(latitude_beijing);
y_beijing= rE*cos(longitude_beijing)*sin(latitude_beijing);
z_beijing= rE*sin(longitude_beijing);

plot3(x_beijing,y_beijing,z_beijing,'pr');
legend('�������й켣','��������','����','�Ź�վλ��');

figure;plot(toe:tspan:toe+tspan*coutmax,lamdakSave)
%% 
norm(savedata1-savedata2);      

% figure(4);
% plot(t0:tspan:t0+tspan*coutmax, v_bar(1,:) , 'r');  hold on;
% plot(t0:tspan:t0+tspan*coutmax, v_bar(2,:) , 'k');  hold on;
% plot(t0:tspan:t0+tspan*coutmax, v_bar(3,:) , 'b');  hold on;
% grid on;
% legend('x���ٶ�','y���ٶ�','z���ٶ�');
% title('�������ϵ������״̬�仯')
% xlabel('Time s')
% ylabel('v m/s')
