


function out = LEOtraCa( t,in, a)
[oe] = TLE2oe('xingli.txt');
a = oe(1);      % �볤��
e = oe(2);
i = oe(3);
Om = oe(4);
om = oe(5);
nu = oe(6);
% --------------------------
mu = 3.986004415e14; % ������
n0 = sqrt(mu/a^3);  % ���ٶ�
t0 = nu;            % ʱ��
M = n0*(t-t0);      % ����ƽ�����
% ---------
sigema = 0.01;      % ����
E0 = 1.2;       % ������õĳ�ֵ
% e  = 0.1;
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
sita = 2*atan( (sqrt((1+e)/(1-e)))*tan(E/2) );

x_a = r*cos(sita);  % ���ٶ�
y_a = r*sin(sita);  % ���ٶ�
z_a = 0;            % ���ٶ�
accel = [x_a,y_a,z_a];   % ���ٶ�
a = accel;
% ---------------------

dx = in(4);
dy = in(5);
dz = in(6);
dvx = a(1);
dvy = a(2);
dvz = a(3);

out = [dx;dy;dz;dvx;dvy;dvz];


end