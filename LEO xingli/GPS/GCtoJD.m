function JD = GCtoJD(Year,Month,DAY) % ������ת��Ϊ������
JD = 367 * Year - floor(7 / 4 * (Year + floor((Month + 9) / 12))) + floor(275 * Month / 9) + DAY + 1721013.5; % floor - ȡ�� 
end