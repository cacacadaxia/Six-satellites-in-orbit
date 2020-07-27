
% =========================================================================
%
%                      ��������λ���ٶȽ���
%
%
% =========================================================================
%
%��(C)2019-2020 ���ݺ���ͨ�����޹�˾
%   �汾��V1.2
%   ���ڣ�2019��7��25��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ����:  1.ʹ�ñ������������ݣ������������ǵķ��й켣
%        2.���ϵ�һ���ļ������������ٶȼ���
%        3.
%        4.
%        5.
%        6.
%--------------------------------------------------------------------------

clear all;
close all;
%---------------------��ȡRINEX��ʽn�ļ�--------------------------------
data = RinexNreader('BR.16n','G0'); % ע���ȡ·�������Ǳ��
%---------------------��������յ��ܻ���---------------------------------
[JD,FOD,GPSW,SOW,DOY,DOW] = GCtoGPS(data(1,1),data(2,1),data(3,1),data(4,1),data(5,1),data(6,1));
t0 = SOW; % SOW-�ܻ���
%-----------------------������������------------------------------------
inx = 1;
for t = t0:1*60:(t0+86400)
    % �� t0 �� t0+86400��(1��) ��� 600��(10min)�������켣����
    % [satPosECEF(:,i) , rECI(:,i) ] = orbitDetermine(data,t);
    GM = 3986005e8; % ��=GM
    we = 7.2921151467e-5; % ������ת���ٶ�
%     we = 0;
    %----------����n�ļ��вο���Ԫ���ܻ���----------
    for j = 1:size(data,2) % size(data��2) - ����data���������;
        [JD,FOD,GPSW,SOW(j),DOY,DOW] = GCtoGPS(2000+data(1,j),data(2,j),data(3,j),data(4,j),data(5,j),data(6,j));
        % ��ָ�����ǵ� �ꡢ�¡��ա�ʱ���֡��� �� ������(��������)��������(С������)��GPS�ܡ��ܻ��롢����ա�������
    end
    %----------�ҳ���t�������������----------
    index = 1;      % ��ʱ���� i==1
    if t >= SOW(size(data,2))
        index = size(data,2);
    elseif t < SOW(1)
    else
        while t >= SOW(index+1)     % ���ʱ�䳬��SOW(i)����ô����SOW(i+1)��ʱ��
            index = index + 1;
        end
    end
    %% ���Ƚ��и�ֵ
    delta_n = data(12,index);
    a       = data(17,index);
    toe     = data(18,index);
    M0      = data(13,index);
    e       = data(15,index);
    om      = data(24,index);
    Om0     = data(20,index);       % OMEGA0 �ǲο�ʱ��toe������ྭ����ʼ�ڸ�����������Ȧ�����ǹ����������ν׼����
    dotOm   = data(25,index);
    i0      = data(22,index);
    doti    = data(26,index);
    %% ������������ǹ��
    %----------1. ���������˶���ƽ�����ٶ�n---------
    n = sqrt(GM) / a^3 + delta_n;
    %----------2. ����黯ʱ��tk----------
    % 604800Ϊһ�ܵ�ʱ��
    tk = t - toe;
    if tk > 302400
        tk = tk - 604800;
    elseif tk < -302400
        tk = tk + 604800;
    end
    %----------3. ����۲�ʱ�̵�ƽ�����Mk----------
    Mk = M0 + n * tk;
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
    setau = data(14,index) * cos(2*phk) + data(16,index) * sin(2*phk);
    setar = data(23,index) * cos(2*phk) + data(11,index) * sin(2*phk);
    setai = data(19,index) * cos(2*phk) + data(21,index) * sin(2*phk);
    %----------8.���㾭������������Ǿࡢ����ʸ����������----------
    uk = phk + setau;
    rk = a^2 * (1 - e * cos(Ek)) + setar;
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
    % ------------12.4 ת������ϵ --------------------
    VxECEF = xk_dot*cos(lamdak) - yk_dot*cos(dotik)*sin(lamdak) + dotik*zECEF*sin(lamdak) - dotLk*yECEF;
    VyECEF = xk_dot*sin(lamdak) + yk_dot*cos(dotik)*cos(lamdak) - dotik*zECEF*cos(lamdak) + dotLk*xECEF;
    VzECEF = yk_dot*sin(ik) + dotik*yECI*cos(dotik);
    
    vECI(:,inx)  = [VxECEF;VyECEF;VzECEF];          % ��WGS-84����ϵ�µ��ٶ�
    % ����i
    inx = inx + 1;
    % others
    lamdakSave(inx) = lamdak;
end

% -----------��ͼ----------------------------
rE     = 6378.137*1e3;      % ����뾶
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