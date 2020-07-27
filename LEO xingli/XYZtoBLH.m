function [B,L,H] = XYZtoBLH(X,Y,Z) % GPS采用WGS-84坐标的参数
a = 6378137;                % 地球的半径
f = 1/298.257223563;        % 地球的扁率
e2 = 2 * f - f^2;
%算大地经度，反正切修正
L = atan(Y/X);
if X < 0
   L = L + pi;
elseif Y < 0 && X > 0
L = L + 2*pi;
end
L= L * 180 / pi - 180;
%大地纬度，设初值，迭代
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