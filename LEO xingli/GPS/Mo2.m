

% =========================================================================
%
%                      球谐函数
%
%
% =========================================================================
%
%　(C)2019-2020 广州海格通信有限公司
%   版本：V1.0
%   日期：2019年8月2日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能:  1.天线角度
%        2.
%        3.
%        4.
%        5.
%        6.
%--------------------------------------------------------------------------
clear;
Nx = 16;
Ny = 32;
lamda = 1;
k   = 2*pi/lamda;
dx  = 0.5*lamda;
dy  = dx;

I = 1;
% sita = 30/180*pi;
% fai  = 60/180*pi;
sitainx0 = 0;
for sitainx = -100:100
    sitainx0 = sitainx0 + 1;
    faiinx0 = 0;
    for faiinx = -100:100
        sita = sitainx/180*pi;
        fai  = faiinx/180*pi;
        S = 0;
        for m = 0:Nx - 1
            for n = 0:Ny - 1
                S = S + I*exp( 1j*k*( m*dx*cos(fai)+n*dy*sin(fai) )*sin(sita) );
            end
        end
        faiinx0 = faiinx0 + 1;
        SS(sitainx0,faiinx0) = S;
    end
end
X = -100:100;
Y = -100:100;
% R = sqrt(X.^2 + Y.^2) + eps;
% Z = sin(R)./R;
C = del2(20*log10(abs(SS)));
mesh(X,Y, 20*log10(abs(SS)),C,'FaceLighting','gouraud','LineWidth',0.3);
% colormap
% pcolor(-100:100, -90:90, 20*log10(abs(SS)));
% colorbar
axis([-100,100,-100,100,0,50])


