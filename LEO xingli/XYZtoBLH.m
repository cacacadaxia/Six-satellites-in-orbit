function [B,L,H] = XYZtoBLH(X,Y,Z) % GPS����WGS-84����Ĳ���
a = 6378137;                % ����İ뾶
f = 1/298.257223563;        % ����ı���
e2 = 2 * f - f^2;
%���ؾ��ȣ�����������
L = atan(Y/X);
if X < 0
   L = L + pi;
elseif Y < 0 && X > 0
L = L + 2*pi;
end
L= L * 180 / pi - 180;
%���γ�ȣ����ֵ������
B = atan(Z / sqrt(X^2 + Y^2));
Btemp = 20;
while abs(B - Btemp) > 1e-11
   N = a / sqrt(1 - e2 * (sin(B))^2);
   H = sqrt(X^2 + Y^2) / cos(B) - N;
   Btemp = B;
   B = atan(Z / sqrt(X^2 + Y^2) / (1 - e2 * N / (N + H)));
end
H = Z / sin(B) - N * (1 - e2);
B = B * 180 / pi;
end