
function data = RinexNreader(file,PRN)
fid = fopen(file,'r+');

% �ж��ļ��Ƿ��
if fid == -1
    %     file!
    error('Cannot open this file!')
end

% �ҵ� END OF HEADER
while 1
    tline = fgets(fid);
    % ����file ? �ļ���ַ;RPN ? ���Ǳ��
    % fopen - ���ļ�;r+ - ���ļ����������д
    % �ж�fid �Ƿ�Ϊ����-1�������������� Cannnot open this
    % while - ����Ϊ true��1 ʱ�ظ�ִ�е� while ѭ�� % fgets - ��ȡ�ļ��е��У����������з�
    if tline(61:73) == 'END OF HEADER' % �жϵ�61-73���Ƿ�Ϊ E O H,����Ҳ���� contains(tline,'END OF HEADER')
        break;
    end
end
% ---------------------------------------
% �ļ���ʶ����Ч
i = 1;
while 1
    tline = fgets(fid);
    if tline == -1
        break
    end
    if tline(1:3) == PRN % PRN��ָʾ�����Ǳ��
        % size - �����С������size(A,1)ֻ������;size(A,2)ֻ������ % fgets - ��ȡ�ļ��е��У����������з�
        % �ж϶�ȡ�����Ƿ������ݣ�û�о�����
        % str2double ? ���ַ���ת��Ϊ˫����ֵ�� �жϵ�1-2λ�ǲ���
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
        % �� % �� % ��
        %ʱ
        %��
        % �������ļ��еĵ�38�С���57�С���76��ֱ�ӱ��E
        %��
        % a0 �����Ӳ����
        % a1 �����Ӳ�Ư����
        % a2 �����Ӳ�Ư������
        % ���ж�ȡ
        % �������ļ��еĵ�19�С���38�С���57�С���76��ֱ�ӱ��E
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
        data(20,i) = str2double(tline([42:60]+1));% OMEGA0 �ǲο�ʱ��toe������ྭ����ʼ�ڸ�����������Ȧ�����ǹ����������ν׼����
        data(21,i) = str2double(tline([61:79]+1));
        tline = fgets(fid);
%         tline([19,38,57,76]) = 'E';
        data(22,i) = str2double(tline([4:22]+1));
        data(23,i) = str2double(tline([23:41]+1));
        data(24,i) = str2double(tline([42:60]+1));
        % IODE ������������������
        % Crs �����ο�ʱ���ڹ�������������ڸ������������� % detn ƽ����ǵĳ��ڱ仯(���ص����)
        % M0 �ο�ʱ�̵�ƽ����ǵ����
        % t0e �����ο�ʱ��(�ܻ���)
        % Cic �����ο�ʱ�̹�����(�����ڷ���)���ڸ���������
        % OMEGA0 �ǲο�ʱ��toe������ྭ����ʼ�ڸ����������� %Cis �����ο�ʱ�̹�����(�����ڷ���)���ڸ�������������
        % ���ж�ȡ
        % �������ļ��еĵ�19�С���38�С���57�С���76��ֱ�ӱ��E
        % Cuc �����ο�ʱ���ڹ���Ӽ����������ڸ������������� % e ƫ����
        % Cus �����ο�ʱ���ڹ���Ӽ����������ڸ������������� % sqrta �������ƽ����
        % ���ж�ȡ
        % �������ļ��еĵ�19�С���38�С���57�С���76��ֱ�ӱ��E
        % ���ж�ȡ
        % �������ļ��еĵ�19�С���38�С���57�С���76��ֱ�ӱ��E
        % i0 �ο�ʱ�̹�����
        % Crc �����ο�ʱ���ڹ�������������ڸ������������� % omega ���ص�Ǿ�
        data(25,i) = str2double(tline([61:79]+1)); % OMEGADOT ������ྭ�ڳ��ƽ���еĳ��ڱ仯
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
        tline = fgets(fid); % ���ж�ȡ
%         tline([19,38]) = 'E';
        data(34,i) = str2double(tline([4:22]+1));
        % �������ļ��еĵ�19�С���38��ֱ�ӱ��E
        % TTOM��Ϣ����ʱ��(����ջ��Խ��յ��������źŽ����й�)
        data(35,i) = str2double(tline([23:41]+1));
        i = i + 1; % ���ڶ���
        % unknow ������������־
    end
end
    fclose(fid); % �ر��ļ�
    
end