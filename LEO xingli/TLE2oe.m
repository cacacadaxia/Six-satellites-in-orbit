%      function [oe,epoch,yr,M,E,satname] = TLE2oe(fname);
%      fname is a filename string for a file containing
%            a two-line element set (TLE)
%      oe is a 1/6 matrix containing the orbital elements
%            [a e i Om om nu]
%      yr is the two-digit year
%      M  is the mean anomaly at epoch
%      E  is the eccentric anomaly at epoch
%      satname is the satellite name
%
% Calls Newton iteration function file EofMe.m

function [oe,epoch,yr,M,E,satname] = TLE2oe(fname)

% Open the file up and scan in the elements

fid = fopen(fname, 'r');
A = fscanf(fid,'%13c%*s',1);
B = fscanf(fid,'%d%6d%*c%5d%*3c%2d%f%f%5d%*c%*d%5d%*c%*d%d%5d',[1,10]);
C = fscanf(fid,'%d%6d%f%f%f%f%f%f',[1,8]);
fclose(fid);
satname=A;

% The value of mu is for the earth
% mu = 3.986004415e5;     % 单位？
miu = 3.986004415e14;     % 单位？
%  Calculate 2-digit year (Oh no!, look out for Y2K bug!)

yr = B(1,4);

% Calculate epoch in julian days
epoch = B(1,5);
%ndot = B(1,6);
% n2dot = B(1,7);

% Assign variables to the orbital elements
i = C(1,3)*pi/180;          % inclination
Om = C(1,4)*pi/180;         % Right Ascension of the Ascending Node
e = C(1,5)/1e7;             % Eccentricity偏心率
om = C(1,6)*pi/180;         % Argument of periapsis
M = C(1,7)*pi/180;          % Mean anomaly平近点角
n = C(1,8)*2*pi/(24*3600);  % 平均角速度


% Calculate the semi-major axis
a = (miu/n^2)^(1/3);     

% Calculate the eccentric anomaly using mean anomaly
E = EofMe(M,e,1e-10);

% Calculate true anomaly from eccentric anomaly
cosnu = (e-cos(E)) / (e*cos(E)-1);
sinnu = ((a*sqrt(1-e*e)) / (a*(1-e*cos(E))))*sin(E);
nu = atan2(sinnu,cosnu);        % 这个是什么？
if (nu<0)
    nu=nu+2*pi; 
end

% Return the orbital elements in a 1x6 matrix
oe = [a e i Om om nu];

% a 轨道椭圆半长轴
% e 偏心率
% i 轨道平面倾角
% Om 升交点赤经
% om 近升角距
% nu 时刻？不对




 