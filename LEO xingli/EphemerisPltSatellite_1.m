

 % =========================================================================
%
%                  ����λ���ٶȽ���
% 
%
% =========================================================================
%
%��(C)2019-2020 ���ݺ���ͨ�����޹�˾
%   �汾��V1.1
%   ���ڣ�2019��7��18��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ����:  1.��ȡ����������Ϣת���ɹ��Ϊ��������Ϣ
%        2. ���ù���������������ǵ��ٶȺ�λ��
%        3. ����ϵ��ת�������忴�ĵ�
%        4. ������Ҫ����"���������㶯"��Ӱ�죬��Ҫ�������ɲ��ϱ仯��
%        5.
%        6.
%--------------------------------------------------------------------------
clear all;
close all;

% ----�����һ����������-------------
[oe,epoch,yr,M,E,satname] = TLE2oe('xingli.txt');    % ����������������
a = oe(1);      % a �볤��
e = oe(2);      % e ƫ����
i = oe(3);      % i ���ƽ�����
Om = oe(4);     % Om ������ྭ
om = oe(5);     % om �����Ǿ࣬���ص�Ǿ࣬���߽н��ص����
nu = oe(6);     % nu ʱ�̣���λ��s��
% ʵ���������и�������ǣ����ﲢ����ʱ�̣����ǽ��ص�ʱ�̣�
% Ӧ����������� sita ֵ��һ����


% ------����----------------------------
aG = 40/180*pi;         % �������κ���ʱ��
tw = 0;                 % ����Ĳ���
GASTw = aG;             % �������κ���ʱ��

t0 = 0;                 % ���ǽ��ص�ʱ�� (��λ��s)
t  = t0;                % none����ʼ��ֵ
tspan = 1;
coutmax = 6.117e+03;
% -------------�Ź�վ-------------------------------------
Hx = 6371e3;        % �߶�
Bx = 40/180*pi;     % γ��
Lx = 60/180*pi;     % ����

% ------------����-----------------------
eE2    = 0.00669437999013;  % ���������ʵ�ƽ��
rE     = 6378.137*1e3;      % ����뾶
owE    = 7.29211567e-5;     % ������ת�Ľ��ٶ�
% ------��γ������ϵ-->�Ź�վ�������ϵ-------------------
N  = rE / (sqrt(1-eE2*sin(Bx)^2));
Xx = (N+Hx)*cos(Bx)*cos(Lx);
Yx = (N+Hx)*cos(Bx)*sin(Lx);
Zx = (N*(1-eE2)+Hx)*sin(Bx);
rX = [Xx;Yx;Zx];            % �Ź�վ��λ��
% -------------------------------------
for ii = 1:coutmax+1
    miu = 3.986004415e14;   % ������
    n0 = sqrt(miu/(a)^3);   % ���ٶ�
    t  = t + tspan;
    M = n0*(t-t0);          % ����ƽ����ǣ�Ӧ���ǲ���ȥ���µģ�
    
    GAST = GASTw + owE*(t-t0);
    % ---------------------------------------------
    sigema = 0.01;      % ����
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
    % figure;plot(Etp)
    r = a*(1-e*cos(E));
    sita = 2*atan( (sqrt((1+e)/(1-e)))*tan(E/2) );  % 
    uk = sita + om;         % ���������������ȡλ��ʱ��ʱ�����ǣ��൱�ڽ�om�ı仯�Ѿ����ǽ�ȥ��
    
    x_a(ii) = r*cos(sita);  % �������ϵ�µ�����λ�ã�ԭ���ǵ�������
    y_a(ii) = r*sin(sita);
    z_a(ii) = 0;
    % -------�ٶ�------------------------
    Ekdot  = n0/(1-e*sin(E));
    Faidot = sqrt((1+e)/(1-e))*(cos(sita/2)^2)/(cos(E/2)^2)*Ekdot;
    ukdot  = Faidot;
    rkdot  = Ekdot*a*e*sin(E);
    vx(ii) = rkdot*cos(uk) - ukdot*y_a(ii);
    vy(ii) = rkdot*sin(uk) + ukdot*x_a(ii);
    vz(ii) = 0;
end
% ------------------------------------
R1 = [cos(om),-sin(om),0 ; sin(om),cos(om),0 ; 0 ,0 , 1];
R2 = [1 ,0 ,0 ; 0 ,cos(i),-sin(i) ; 0 ,sin(i),cos(i)];
R3 = [cos(Om-GAST),-sin(Om-GAST),0 ; sin(Om-GAST),cos(Om-GAST),0 ; 0 ,0 ,1];    % Lk
r_bar = R1*R2*R3*[ x_a; y_a; z_a ];     % �����ڴ������ϵ��Ӧ�þ���WGS-84����ϵ
v_bar = R2*R3*[ vx ; vy ; vz ];         % �����ڴ������ϵ



fprintf('��ʼʱ�̣�%3.0f , ����ʱ�̣�%3.0f , \n',t0,t);
% -----------��ͼ----------------------------
figure(1);plot3(r_bar(1,:),r_bar(2,:),r_bar(3,:),'r','LineWidth',1);
hold on;plot3(0,0,0,'*');grid on;
legend('�������й켣','��������');

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

%% 

% figure(4);
% plot(t0:tspan:t0+tspan*coutmax, v_bar(1,:) , 'r');  hold on;
% plot(t0:tspan:t0+tspan*coutmax, v_bar(2,:) , 'k');  hold on;
% plot(t0:tspan:t0+tspan*coutmax, v_bar(3,:) , 'b');  hold on;
% grid on;
% legend('x���ٶ�','y���ٶ�','z���ٶ�');
% title('�������ϵ������״̬�仯')
% xlabel('Time s')
% ylabel('v m/s')
