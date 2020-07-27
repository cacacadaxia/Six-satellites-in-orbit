
function data = RinexNreader(file,PRN)
fid = fopen(file,'r+');

% 判断文件是否打开
if fid == -1
    %     file!
    error('Cannot open this file!')
end

% 找到 END OF HEADER
while 1
    tline = fgets(fid);
    % 输入file ? 文件地址;RPN ? 卫星编号
    % fopen - 打开文件;r+ - 打开文件，允许读、写
    % 判断fid 是否为有误-1，如果有误则输出 Cannnot open this
    % while - 条件为 true即1 时重复执行的 while 循环 % fgets - 读取文件中的行，并保留换行符
    if tline(61:73) == 'END OF HEADER' % 判断第61-73列是否为 E O H,本句也可用 contains(tline,'END OF HEADER')
        break;
    end
end
% ---------------------------------------
% 文件标识符无效
i = 1;
while 1
    tline = fgets(fid);
    if tline == -1
        break
    end
    if tline(1:3) == PRN % PRN所指示的卫星编号
        % size - 数组大小，其中size(A,1)只求行数;size(A,2)只求列数 % fgets - 读取文件中的行，并保留换行符
        % 判断读取的行是否有数据，没有就跳出
        % str2double ? 将字符串转换为双精度值， 判断第1-2位是不是
        data(1,i) = str2double(tline([3:5+2]+1));
        data(2,i) = str2double(tline(10:11));
        data(3,i) = str2double(tline(13:14));
        data(4,i) = str2double(tline(16:17));
        data(5,i) = str2double(tline(19:20));
%         tline([38,57,76]) = 'E';
        data(6,i) = str2double(tline(22:23));
        data(7,i) = str2double(tline([23:41]+1));
        data(8,i) = str2double(tline([42:60]+1));
        data(9,i) = str2double(tline([61:79]+1));
        % 年 % 月 % 日
        %时
        %分
        % 将数据文件中的第38列、第57列、第76列直接变成E
        %秒
        % a0 卫星钟差常数项
        % a1 卫星钟差漂移项
        % a2 卫星钟差漂移速率
        % 新行读取
        % 将数据文件中的第19列、第38列、第57列、第76列直接变成E
        tline = fgets(fid);
%         tline([19,38,57,76]) = 'E';
        data(10,i) = str2double(tline([4:22]+1));
        data(11,i) = str2double(tline([23:41]+1));
        data(12,i) = str2double(tline([42:60]+1));
        data(13,i) = str2double(tline([61:79]+1));
        tline = fgets(fid);
%         tline([19,38,57,76]) = 'E';
        data(14,i) = str2double(tline([4:22]+1));
        data(15,i) = str2double(tline([23:41]+1));
        data(16,i) = str2double(tline([42:60]+1));
        data(17,i) = str2double(tline([61:79]+1));
        tline = fgets(fid); 
%         tline([19,38,57,76]) = 'E';
        data(18,i) = str2double(tline([4:22]+1));
        data(19,i) = str2double(tline([23:41]+1));
        data(20,i) = str2double(tline([42:60]+1));% OMEGA0 非参考时刻toe升交点赤经，是始于格林尼治子午圈到卫星轨道升交点所谓准经度
        data(21,i) = str2double(tline([61:79]+1));
        tline = fgets(fid);
%         tline([19,38,57,76]) = 'E';
        data(22,i) = str2double(tline([4:22]+1));
        data(23,i) = str2double(tline([23:41]+1));
        data(24,i) = str2double(tline([42:60]+1));
        % IODE 卫星星历的数据龄期
        % Crs 星历参考时刻在轨道径向方向上周期改正余弦项的振幅 % detn 平近点角的长期变化(近地点参数)
        % M0 参考时刻的平近点角的振幅
        % t0e 星历参考时刻(周积秒)
        % Cic 星历参考时刻轨道倾角(近似于法向)周期改正余弦项
        % OMEGA0 非参考时刻toe升交点赤经，是始于格林尼治子午 %Cis 星历参考时刻轨道倾角(近似于法向)周期改正余弦项的振幅
        % 新行读取
        % 将数据文件中的第19列、第38列、第57列、第76列直接变成E
        % Cuc 星历参考时刻在轨道延迹方向上周期改正余弦项的振幅 % e 偏心率
        % Cus 星历参考时刻在轨道延迹方向上周期改正余弦项的振幅 % sqrta 长半轴的平方根
        % 新行读取
        % 将数据文件中的第19列、第38列、第57列、第76列直接变成E
        % 新行读取
        % 将数据文件中的第19列、第38列、第57列、第76列直接变成E
        % i0 参考时刻轨道倾角
        % Crc 星历参考时刻在轨道径向方向上周期改正余弦项的振幅 % omega 近地点角距
        data(25,i) = str2double(tline([61:79]+1)); % OMEGADOT 升交点赤经在赤道平面中的长期变化
        tline = fgets(fid);
%         tline([19,38,57,76]) = 'E';
        data(26,i) = str2double(tline([4:22]+1));
        data(27,i) = str2double(tline([23:41]+1));
        data(28,i) = str2double(tline([42:60]+1));
        data(29,i) = str2double(tline([61:79]+1));
        tline = fgets(fid);
%         tline([19,38,57,76]) = 'E';
        data(30,i) = str2double(tline([4:22]+1));
        data(31,i) = str2double(tline([23:41]+1));
        data(32,i) = str2double(tline([42:60]+1));
        data(33,i) = str2double(tline([61:79]+1));
        tline = fgets(fid); % 新行读取
%         tline([19,38]) = 'E';
        data(34,i) = str2double(tline([4:22]+1));
        % 将数据文件中的第19列、第38列直接变成E
        % TTOM信息传送时间(与接收机对接收到的卫星信号解码有关)
        data(35,i) = str2double(tline([23:41]+1));
        i = i + 1; % 读第二个
        % unknow 星历拟合区间标志
    end
end
    fclose(fid); % 关闭文件
    
end