



clear all;
addpath('ds');

load datasave01;
d1 = delta_f;
load datasave02;
d2 = delta_f;
load datasave03;
d3 = delta_f;
load datasave04;
d4 = delta_f;

%%
tspan = 1;
t = 0:tspan:length(delta_f)-1;

figure;
plot(t ,d1,'LineWidth',1);hold on;
plot(t ,d2,'LineWidth',1);hold on;
plot(t ,d3,'LineWidth',1);hold on;
plot(t ,d4,'LineWidth',1);hold on;
legend('UE���ٶ�','UE���ž����˶�','UE����γ���˶�','UE��������ĳ�����˶�')
xlabel('ʱ��t s')
ylabel('Ƶƫ\deltaf Hz')



