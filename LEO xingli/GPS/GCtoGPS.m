function [JD,FOD,GPSW,SOW,DOY,DOW] = GCtoGPS(Y,M,D,h,m,s) %将测量日 格里历→GPS时 % JD-儒略日整数，FOD-儒略日小数，GPSW-GPS周数，SOW-GPS周积秒，DOY-年积日，DOW-星期数
JD = GCtoJD(Y,M,D); % 调用函数: GCtoJD
FOD = h / 24 + m / 1440 + s / 86400;
GPSW = floor((JD - GCtoJD(1980,1,6))/7); % floor - 取整
DOW = floor(rem(JD+1.5,7)); % rem ? 取余，r=rem(a,b)表示a除以b取余数 
SOW = DOW * 86400 + h * 3600 + m * 60 + s; % 周积秒
L = [31,28,31,30,31,30,31,31,30,31,30,31];
if Y == 4 * floor(Y/4) % 判断Y是不是闰年
    L(2) = 29; % 如果是闰年，将L中第二个数改为29天
end
DOY = 0;
for i = 1:M-1
    DOY = DOY + L(i);
end
DOY = DOY + D;
end