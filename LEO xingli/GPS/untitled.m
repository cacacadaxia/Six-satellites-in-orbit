%%
%采用背景擦除的方法，动态的划点，并且动态改变坐标系
% t,m 均为一行 ，并且不能为多行
t = 0;
m = 0;
n = 0;
p = plot3(0,0,0,'*');
x=-1.5*pi;
axis([x x+2*pi -1.5 1.5]);
grid on;

for i=1:1000
    t=0.1*i;       %两个变量均不追加
    m=sin(0.1*i);
    n = 0;
    set(p,'XData',t,'YData',m,'ZData',n)
    x=x+0.1;
    drawnow
    axis([x x+2*pi -1.5 1.5]);
    pause(0.1);
end
