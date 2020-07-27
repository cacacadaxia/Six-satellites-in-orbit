% clear;
% figure
% grs80=almanac('earth','grs80');
% ax=axesm('globe','Geoid',grs80,'Grid','on','GLineStyle','-','Gcolor','yellow');
% set(ax,'Position',[0 0 1 1]);
% view(3);
% axis equal off vis3d;
% set(gcf,'Renderer','opengl');
% load topo
% geoshow(topo,topolegend,'DisplayType','texturemap');
% demcmap(topo);


cla reset;
load topo;
[x,y,z] = sphere(24);
s = surface(x,y,z,'FaceColor','texturemap','CData',topo);
colormap(topomap1);
view(0,50)
% Brighten the colormap for better annotation visibility:
brighten(.8)
% Create and arrange the camera and lighting for better visibility:
% campos([1.3239? -14.4250? 9.4954]);
% campos([0 0 0]);
% % camlight;
% lighting gouraud;
% axis off vis3d;
% % Set the x- and y-coordinates of the textarrow object:
% x = [0.7698 0.5851];
% y = [0.3593 0.5492];
% % Create the textarrow object:?
% txtar =annotation('textarrow',x,y,'String','This is China.','color','red','FontSize',14);
