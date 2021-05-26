fname = '~/data/ocean_wmass.ncra.nc';
yt = ncread(fname,'grid_yt_ocean');
T = ncread(fname,'neutral');
Te = ncread(fname,'neutralrho_edges');
rho0 = 1035;

% load and zonal sum:
tmp = squeeze(nansum(ncread(fname,'ty_trans_nrho')/rho0,1));
[yL,TL] = size(tmp);
tmp(isnan(tmp)) = 0;
psi = squeeze(cat(2,zeros(yL,1),cumsum(tmp,2)));

tmp = squeeze(nansum(ncread(fname,'ty_trans_nrho_gm')/rho0,1));
tmp(isnan(tmp)) = 0;
psi_gm = tmp;
tmp = squeeze(nansum(ncread(fname,'ty_trans_nrho_submeso')/rho0,1));
tmp(isnan(tmp)) = 0;
psi_submeso = tmp;

[Xe,Ye] = ndgrid(yt,Te);
[X,Y] = ndgrid(yt,T);

figure;
subplot(2,2,1);
pcolPlot(Xe,Ye,psi/1e6);
caxis([-20 20]);
xlabel('Latitude ($^\circ$E)');
ylabel('Temperature ($^\circ$C)');
cb = colorbar;
ylabel(cb,'Sv');
title('Resolved overturning');
subplot(2,2,2);
pcolPlot(X,Y,psi_gm/1e6);
caxis([-0.2 0.2]);
xlabel('Latitude ($^\circ$E)');
ylabel('Temperature ($^\circ$C)');
cb = colorbar;
ylabel(cb,'Sv');
title('GM overturning');
subplot(2,2,3);
pcolPlot(X,Y,psi_submeso/1e6);
caxis([-0.05 0.05]);
xlabel('Latitude ($^\circ$E)');
ylabel('Temperature ($^\circ$C)');
cb = colorbar;
ylabel(cb,'Sv');
title('Submeso overturning');
colormap(redblue);


