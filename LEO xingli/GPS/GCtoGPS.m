function [JD,FOD,GPSW,SOW,DOY,DOW] = GCtoGPS(Y,M,D,h,m,s) %�������� ��������GPSʱ % JD-������������FOD-������С����GPSW-GPS������SOW-GPS�ܻ��룬DOY-����գ�DOW-������
JD = GCtoJD(Y,M,D); % ���ú���: GCtoJD
FOD = h / 24 + m / 1440 + s / 86400;
GPSW = floor((JD - GCtoJD(1980,1,6))/7); % floor - ȡ��
DOW = floor(rem(JD+1.5,7)); % rem ? ȡ�࣬r=rem(a,b)��ʾa����bȡ���� 
SOW = DOW * 86400 + h * 3600 + m * 60 + s; % �ܻ���
L = [31,28,31,30,31,30,31,31,30,31,30,31];
if Y == 4 * floor(Y/4) % �ж�Y�ǲ�������
    L(2) = 29; % ��������꣬��L�еڶ�������Ϊ29��
end
DOY = 0;
for i = 1:M-1
    DOY = DOY + L(i);
end
DOY = DOY + D;
end