addpath(genpath('/short/e14/rmh561/software/matlab-utilities/'));
startup;
name = 'ocean_wmass.nc';

Cp = 3992.10322329649; % J kg-1 degC-1
rho0 = 1035; % kgm-3
tsc = 1e9;

Te = ncread(name,'neutralrho_edges');
T = ncread(name,'neutral');
lonv = ncread(name,'xt_ocean');
latv = ncread(name,'yt_ocean');
[tmp ln1] = min(abs(lonv+90));
[tmp ln2] = min(abs(lonv-45));
[tmp lt1] = min(abs(latv+80));
[tmp lt2] = min(abs(latv+20));

lon = ncread(name,'geolon_t');
lat = ncread(name,'geolat_t');
lon = lon(ln1:ln2,lt1:lt2);
lat = lat(ln1:ln2,lt1:lt2);
lonv = lonv(ln1:ln2);
latv = latv(lt1:lt2);
[xL,yL] = size(lon)
TL = length(Te)-1;

txtrans = zeros(xL,yL,TL);
hxtrans = zeros(xL,yL,TL);
for ti=1:12
    txtrans = txtrans+ncread(name,'tx_trans_nrho',[ln1 lt1 1 ti],[xL yL 74 1])*tsc/rho0;
    hxtrans = txtrans+ncread(name,'temp_xflux_adv_on_nrho',[ln1 lt1 1 ti],[xL yL 74 1]);
    ti
end
txtrans = txtrans/12;
txtrans(isnan(txtrans)) = 0;
hxtrans = hxtrans/12;
hxtrans(isnan(hxtrans)) = 0;

% $$$ tytrans = zeros(xL,yL,TL);
% $$$ for ti=1:12
% $$$     ti
% $$$     tytrans = tytrans+ncread(name,'ty_trans_nrho',[ln1 lt1 1 ti],[xL yL 74 1])*tsc/rho0;
% $$$ end
% $$$ tytrans = tytrans/12;
% $$$ tytrans(isnan(tytrans)) = 0;

% Warm Route:
[tmp ind] = min(abs(lonv-20));
PSI = cumsum(cat(2,zeros(yL,1),squeeze(cumsum(txtrans(ind,:,:),3))),1,'reverse');
A = cumsum(cat(2,zeros(yL,1),squeeze(cumsum(hxtrans(ind,:,:),3))),1,'reverse');
AE = rho0*Cp*PSI.*repmat(Te',[yL 1]);
AI = A-AE;
ind2 = find(PSI(:,end)<0,1,'first');
AI(ind2,end)/1e15
A(ind2,end)/1e15
AE(ind2,end)/1e15
xlims = [-60 -30];

% Cold Route:
[tmp ind] = min(abs(lonv+68));
PSI = cumsum(cat(2,zeros(yL,1),squeeze(cumsum(txtrans(ind,:,:),3))),1,'reverse');
A = cumsum(cat(2,zeros(yL,1),squeeze(cumsum(hxtrans(ind,:,:),3))),1,'reverse');
AE = rho0*Cp*PSI.*repmat(Te',[yL 1]);
AI = A-AE;
[tmp ind2] = min(abs(latv+60));
AI(ind2,end)/1e15
A(ind2,end)/1e15
AE(ind2,end)/1e15
xlims = [-70 -50];


[X,Y] = ndgrid(latv,Te);

%figure;
clf;
subplot(3,2,1);
VIX = sum(txtrans(1:4:end,1:4:end,:),3);
VIX(VIX==0) = NaN;
pcolPlot(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),VIX/1e6);
hold on;
plot([lonv(ind) lonv(ind)],[latv(1) latv(end)],'--k');
cb = colorbar;
ylabel(cb,'Sv/cell');
caxis([-5 5]);
title('Vertically-integrated zonal transport');
set(gca,'color','k');

subplot(3,2,2);
pcolPlot(X,Y,PSI/1e6);
hold on;
plot([latv(ind2) latv(ind2)],[Te(1) Te(end)],'--k');
caxis([-200 200]);
cb = colorbar;
ylabel(cb,'Sv');
colormap(redblue);
ylim([0 34]);
title('$\Psi(y,\Theta)$');
xlabel('Latitude ($^\circ$N)');
ylabel('Temperature ($^\circ$C)');
xlim(xlims);

subplot(3,2,3);
pcolPlot(X,Y,A/1e15);
hold on;
plot([latv(ind2) latv(ind2)],[Te(1) Te(end)],'--k');
cb = colorbar;
ylabel(cb,'PW');
%caxis([-0.2 0.2]);
caxis([-0.4 0.4]);
colormap(redblue);
ylim([0 34]);
title('$\mathcal{A}(y,\Theta)$');
xlabel('Latitude ($^\circ$N)');
ylabel('Temperature ($^\circ$C)');
xlim(xlims);

subplot(3,2,4);
pcolPlot(X,Y,AE/1e15);
hold on;
plot([latv(ind2) latv(ind2)],[Te(1) Te(end)],'--k');
cb = colorbar;
ylabel(cb,'PW');
caxis([-20 20]);
colormap(redblue);
ylim([0 34]);
title('$\mathcal{A}_E(y,\Theta)$');
xlabel('Latitude ($^\circ$N)');
ylabel('Temperature ($^\circ$C)');
xlim(xlims);

subplot(3,2,5);
pcolPlot(X,Y,AI/1e15);
hold on;
plot([latv(ind2) latv(ind2)],[Te(1) Te(end)],'--k');
cb = colorbar;
ylabel(cb,'PW');
caxis([-20 20]);
colormap(redblue);
ylim([0 34]);
title('$\mathcal{A}_I(y,\Theta)$');
xlabel('Latitude ($^\circ$N)');
ylabel('Temperature ($^\circ$C)');
xlim(xlims);

subplot(3,2,6);
plot(A(ind2,:)/1e15,Te,'-k');
hold on;
plot(AI(ind2,:)/1e15,Te,'-r');
plot(AE(ind2,:)/1e15,Te,'-b');
plot(A(ind2+1,:)/1e15,Te,'--k');
plot(AI(ind2+1,:)/1e15,Te,'--r');
plot(AE(ind2+1,:)/1e15,Te,'--b');
plot(A(ind2-1,:)/1e15,Te,':k');
plot(AI(ind2-1,:)/1e15,Te,':r');
plot(AE(ind2-1,:)/1e15,Te,':b');
lat2 = sprintf('%5.1f',latv(ind2));
lat2p1 = sprintf('%5.1f',latv(ind2+1));
lat2m1 = sprintf('%5.1f',latv(ind2-1));
legend(['A(' lat2 ')'], ...
       ['AI(' lat2 ')'], ...
       ['AE(' lat2 ')'], ...
       ['A(' lat2p1 ')'], ...
       ['AI(' lat2p1 ')'], ...
       ['AE(' lat2p1 ')'], ...
       ['A(' lat2m1 ')'], ...
       ['AI(' lat2m1 ')'], ...
       ['AE(' lat2m1 ')']);
xlabel('Westward Transport (PW)');
ylabel('Temperature ($^\circ$C)');
ylim([0 34]);

