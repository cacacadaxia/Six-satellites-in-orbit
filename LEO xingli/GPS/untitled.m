%%
%���ñ��������ķ�������̬�Ļ��㣬���Ҷ�̬�ı�����ϵ
% t,m ��Ϊһ�� �����Ҳ���Ϊ����
t = 0;
m = 0;
n = 0;
p = plot3(0,0,0,'*');
x=-1.5*pi;
axis([x x+2*pi -1.5 1.5]);
grid on;

for i=1:1000
    t=0.1*i;       %������������׷��
    m=sin(0.1*i);
    n = 0;
    set(p,'XData',t,'YData',m,'ZData',n)
    x=x+0.1;
    drawnow
    axis([x x+2*pi -1.5 1.5]);
    pause(0.1);
end
