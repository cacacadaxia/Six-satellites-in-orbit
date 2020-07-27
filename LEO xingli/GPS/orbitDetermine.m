function [satPos,rECI] = orbitDetermine(data,t)
% 功能:卫星轨道计算
% 输入:卫星编号 data 、时间 t
% 输出:卫星WGS-84坐标系x、y、z坐标 - satPos(1)、satPos(2)、satPos(3)
% 常量----------------------------------------
GM = 3986005e8; % μ=GM
we = 7.2921151467e-5; % 地球自转角速度

%----------计算n文件中参考历元的周积秒----------
for j = 1:size(data,2) % size(data，2) - 返回data矩阵的列数;
    [JD,FOD,GPSW,SOW(j),DOY,DOW] = GCtoGPS(data(1,j),data(2,j),data(3,j),data(4,j),data(5,j),data(6,j));
    % 将指定卫星的 年、月、日、时、分、秒 → 儒略日(整数部分)、儒略日(小数部分)、GPS周、周积秒、年积日、星期数
end
%----------找出与t最近的那列数据----------
i = 1;      % 暂时先让 i==1
if t >= SOW(size(data,2))
    i = size(data,2);
elseif t < SOW(1)
    i = 1;
else
    while t >= SOW(i+1)     % 如果时间超过SOW(i)，那么就找SOW(i+1)的时间
        i = i + 1;
    end
end
% disp(i);

%% 按步骤计算卫星轨道
%----------1. 计算卫星运动的平均角速度n---------
n = sqrt(GM) / data(17,i)^3 + data(12,i);
%----------2. 计算归化时间tk----------
% 604800为一周的时间
tk = t - data(18,i);
if tk > 302400
    tk = tk - 604800;
elseif tk < -302400
    tk = tk + 604800;
end
%----------3. 计算观测时刻的平近点角Mk----------
Mk = data(13,i) + n * tk;
%----------4. 计算观测时刻的偏近点角Ek----------
Ek = Mk;
Etemp = 10e6;
while 1
    if abs(Ek - Etemp) < 1e-10
        break;
    end
    Etemp = Ek;
    Ek = Mk + data(15,i) * sin(Ek);
end
%----------5. 计算观测时刻的真近点角Vk------------
Vkc = (cos(Ek)-data(15,i))/(1-data(15,i)*cos(Ek));
Vks = ((1-(data(15,i))^2)^0.5*sin(Ek))/(1-data(15,i)*cos(Ek));
Vk = atan2(Vks,Vkc);
%----------6. 计算升交角距----------
phk = Vk + data(24,i);
%----------7. 计算摄动改正项----------
setau = data(14,i) * cos(2*phk) + data(16,i) * sin(2*phk);
setar = data(23,i) * cos(2*phk) + data(11,i) * sin(2*phk);
setai = data(19,i) * cos(2*phk) + data(21,i) * sin(2*phk);
%----------8.计算经改正后的升交角距、卫星矢径、轨道倾角----------
uk = phk + setau;
rk = data(17,i)^2 * (1 - data(15,i) * cos(Ek)) + setar;
ik = data(22,i) + setai + data(26,i) * tk;              % 轨道平面倾角
%----------9. 计算卫星在轨道坐标系的位置----------
xECI = rk * cos(uk);
yECI = rk * sin(uk);

rECI = [xECI;yECI;0];       % 轨道坐标系

%----------10. 计算观测时刻t的升交点赤径----------
lamdak = data(20,i) + (data(25,i) - we) * tk - we * data(18,i);
%----------11. 计算卫星在WGS-84坐标系的位置----------
xECEF = xECI * cos(lamdak) - yECI * cos(ik) * sin(lamdak);
yECEF = xECI * sin(lamdak) + yECI * cos(ik) * cos(lamdak);
zECEF = yECI * sin(ik);
%----------输出结果----------
satPos(1) = xECEF;
satPos(2) = yECEF;
satPos(3) = zECEF;

end
