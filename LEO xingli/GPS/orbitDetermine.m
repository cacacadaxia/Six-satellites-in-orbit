function [satPos,rECI] = orbitDetermine(data,t)
% ����:���ǹ������
% ����:���Ǳ�� data ��ʱ�� t
% ���:����WGS-84����ϵx��y��z���� - satPos(1)��satPos(2)��satPos(3)
% ����----------------------------------------
GM = 3986005e8; % ��=GM
we = 7.2921151467e-5; % ������ת���ٶ�

%----------����n�ļ��вο���Ԫ���ܻ���----------
for j = 1:size(data,2) % size(data��2) - ����data���������;
    [JD,FOD,GPSW,SOW(j),DOY,DOW] = GCtoGPS(data(1,j),data(2,j),data(3,j),data(4,j),data(5,j),data(6,j));
    % ��ָ�����ǵ� �ꡢ�¡��ա�ʱ���֡��� �� ������(��������)��������(С������)��GPS�ܡ��ܻ��롢����ա�������
end
%----------�ҳ���t�������������----------
i = 1;      % ��ʱ���� i==1
if t >= SOW(size(data,2))
    i = size(data,2);
elseif t < SOW(1)
    i = 1;
else
    while t >= SOW(i+1)     % ���ʱ�䳬��SOW(i)����ô����SOW(i+1)��ʱ��
        i = i + 1;
    end
end
% disp(i);

%% ������������ǹ��
%----------1. ���������˶���ƽ�����ٶ�n---------
n = sqrt(GM) / data(17,i)^3 + data(12,i);
%----------2. ����黯ʱ��tk----------
% 604800Ϊһ�ܵ�ʱ��
tk = t - data(18,i);
if tk > 302400
    tk = tk - 604800;
elseif tk < -302400
    tk = tk + 604800;
end
%----------3. ����۲�ʱ�̵�ƽ�����Mk----------
Mk = data(13,i) + n * tk;
%----------4. ����۲�ʱ�̵�ƫ�����Ek----------
Ek = Mk;
Etemp = 10e6;
while 1
    if abs(Ek - Etemp) < 1e-10
        break;
    end
    Etemp = Ek;
    Ek = Mk + data(15,i) * sin(Ek);
end
%----------5. ����۲�ʱ�̵�������Vk------------
Vkc = (cos(Ek)-data(15,i))/(1-data(15,i)*cos(Ek));
Vks = ((1-(data(15,i))^2)^0.5*sin(Ek))/(1-data(15,i)*cos(Ek));
Vk = atan2(Vks,Vkc);
%----------6. ���������Ǿ�----------
phk = Vk + data(24,i);
%----------7. �����㶯������----------
setau = data(14,i) * cos(2*phk) + data(16,i) * sin(2*phk);
setar = data(23,i) * cos(2*phk) + data(11,i) * sin(2*phk);
setai = data(19,i) * cos(2*phk) + data(21,i) * sin(2*phk);
%----------8.���㾭������������Ǿࡢ����ʸ����������----------
uk = phk + setau;
rk = data(17,i)^2 * (1 - data(15,i) * cos(Ek)) + setar;
ik = data(22,i) + setai + data(26,i) * tk;              % ���ƽ�����
%----------9. ���������ڹ������ϵ��λ��----------
xECI = rk * cos(uk);
yECI = rk * sin(uk);

rECI = [xECI;yECI;0];       % �������ϵ

%----------10. ����۲�ʱ��t��������ྶ----------
lamdak = data(20,i) + (data(25,i) - we) * tk - we * data(18,i);
%----------11. ����������WGS-84����ϵ��λ��----------
xECEF = xECI * cos(lamdak) - yECI * cos(ik) * sin(lamdak);
yECEF = xECI * sin(lamdak) + yECI * cos(ik) * cos(lamdak);
zECEF = yECI * sin(ik);
%----------������----------
satPos(1) = xECEF;
satPos(2) = yECEF;
satPos(3) = zECEF;

end
