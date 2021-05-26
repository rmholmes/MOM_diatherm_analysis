
% Load one year of data:
load('MOM025_kb3seg_output120_VertInt_T22p5C.mat');
load('MOM025_kb3seg_output120_BaseVars.mat');

figure;
set(gcf,'Position',get(0,'ScreenSize'));
subplot(2,2,1);
pcolor(lon,lat,FlM);
shading flat;
caxis([-120 120]);
colormap(redblue);
cb = colorbar;
ylabel(cb,'Wm$^{-2}$');
set(gca,'color','k');
title('Total vertical mixing flux across $22.5^\circ$C isotherm [One year MOM025] - Fig 4a');

subplot(2,2,2);
pcolor(lon,lat,FlMkppish);
shading flat;
caxis([-120 120]);
colormap(redblue);
cb = colorbar;
ylabel(cb,'Wm$^{-2}$');
set(gca,'color','k');
title('Shear-driven vertical mixing flux across $22.5^\circ$C isotherm [One year MOM025] - Fig. 4c');

