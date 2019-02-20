% A quick script to calculate eddy volume transports across 5N in
% temperature intervals and their potential contribution to WWV
% changes.

% $$$ %%% Get data and save in .mat files (run on raijin):
% $$$ 
% $$$ base = '/g/data/e14/rmh561/mom/';
% $$$ model = 'MOM_HeatDiag';
% $$$ switchtransfilename = 0;
% $$$ 
% $$$ base = '/g/data/e14/mv7494/mom/archive/';
% $$$ model = 'pnEXP2_restart000_windstress';
% $$$ switchtransfilename = 1;
% $$$ 
% $$$ % $$$ for output=1:4
% $$$ for output=0:3
% $$$     output
% $$$ 
% $$$ 
% $$$     fold = [base model '/' sprintf('output%03d',output) '/'];
% $$$     fname = [fold 'ocean.nc'];
% $$$     wname = [fold 'ocean_wmass.nc'];
% $$$     gname = [fold 'ocean_grid.nc'];
% $$$ 
% $$$     lon = ncread(gname,'geolon_t');lat = ncread(gname,'geolat_t');
% $$$     lonu = ncread(gname,'geolon_c');latu = ncread(gname,'geolat_c');
% $$$     area = ncread(gname,'area_t');[xL,yL] = size(lon);
% $$$     lonv_t = ncread(gname,'xt_ocean');lonv_u = ncread(gname,'xu_ocean');
% $$$     latv_t = ncread(gname,'yt_ocean');latv_u = ncread(gname,'yu_ocean');
% $$$     z = ncread(fname,'st_ocean');zL = length(z);
% $$$ 
% $$$     mask = ncread(fname,'temp',[1 1 1 1],[xL yL zL 1]);
% $$$     mask(~isnan(mask)) = 1; mask(isnan(mask)) = 0;
% $$$     mask = mask == 1;
% $$$ 
% $$$     time = ncread(fname,'time');
% $$$     tL = length(time);
% $$$ 
% $$$ % $$$ dys = [31 28 31 30 31 30 31 31 30 31 30 31];
% $$$ % $$$ Cp = 3992.1; % J kg-1 degC-1
% $$$     rho0 = 1035; % kgm-3
% $$$     T = ncread(wname,'neutral');
% $$$     Te = ncread(wname,'neutralrho_edges');
% $$$     TL = length(T);dT = T(2)-T(1);
% $$$ 
% $$$     % Calculate eddy transports at 5N:
% $$$ 
% $$$     [tmp ltind_u] = min(abs(latv_u-5));
% $$$     [tmp ltind_t] = min(abs(latv_t-5));
% $$$ 
% $$$     % online:
% $$$     if (~switchtransfilename)
% $$$         tname = wname;
% $$$     else
% $$$         tname = fname;
% $$$     end
% $$$     ty_trans_theta_on = squeeze(ncread(tname,'ty_trans_nrho',[1 ltind_u 1 1],[xL 1 TL tL]))*1e9/rho0;
% $$$ 
% $$$     % offline (XXX: This calculation is not exactly copying the model
% $$$     % because of the lack of the exact thickness dzt on u-points):
% $$$     temp = squeeze(ncread(fname,'temp',[1 ltind_t 1 1],[xL 2 zL tL]));
% $$$     dzt = squeeze(ncread(fname,'dzt',[1 ltind_t 1 1],[xL 2 zL tL]));
% $$$     temp = squeeze(nanmean(temp,2));
% $$$     dzt = squeeze(nanmean(dzt,2));
% $$$     v = squeeze(ncread(fname,'v',[1 ltind_u 1 1],[xL 1 zL tL]));
% $$$     v = cat(1,(v(1,:,:)+v(end,:,:))/2,(v(2:end,:,:)+v(1:(end-1),:,:))/2);
% $$$ 
% $$$     Re = 6371000;
% $$$     dx = (area(:,ltind_t)+area(:,ltind_t-1))/2./(Re*sin((latv_u(ltind_t)-latv_u(ltind_t-1))*pi/180));
% $$$ 
% $$$     ty_trans_theta_off = zeros(size(ty_trans_theta_on));
% $$$     for ti = 1:tL
% $$$         ti
% $$$         for xi=1:xL
% $$$             for i = 2:TL
% $$$                 inds = temp(xi,:,ti) > Te(i-1) & temp(xi,:,ti) <= Te(i);
% $$$                 ty_trans_theta_off(xi,i-1,ti) = nansum(v(xi,inds,ti).*dzt(xi,inds,ti))*dx(xi);
% $$$             end
% $$$         end
% $$$     end
% $$$ 
% $$$     [tmp lnind1] = min(abs(lonv_u+260));
% $$$     [tmp lnind2] = min(abs(lonv_u+70));
% $$$ 
% $$$     ty_trans_theta_on = ty_trans_theta_on(lnind1:lnind2,:,:);
% $$$     ty_trans_theta_off = ty_trans_theta_off(lnind1:lnind2,:,:);
% $$$     lonv = lonv_u(lnind1:lnind2);
% $$$     time = [1:12];
% $$$ 
% $$$     save(['/short/e14/rmh561/trans5N_' model '_' sprintf('output%03d',output) '.mat'],'ty_trans_theta_on','ty_trans_theta_off',...
% $$$          'lonv','time','T');
% $$$ end

%% Plotting (on linux desktop):

base = '/srv/ccrc/data03/z3500785/mom/5Ntrans/';

% Control:
model = 'MOM_HeatDiag';
outputs = [15:19];
load([base 'trans5N_' model '_' sprintf('output%03d',outputs(1)) '.mat']);
ty_ona = ty_trans_theta_on;
ty_offa = ty_trans_theta_off;
for i = 2:length(outputs)
    load([base 'trans5N_' model '_' sprintf('output%03d',outputs(i)) '.mat']);
    ty_ona = ty_ona + ty_trans_theta_on;
    ty_offa = ty_offa + ty_trans_theta_off;
end
ty_on_CT = ty_ona/length(outputs);
ty_off_CT = ty_offa/length(outputs);
clear ty_ona ty_offa ty_on ty_off;

% El Nino:
model = 'pnEXP1_restart000_windstress';
outputs = [0:3];
load([base 'trans5N_' model '_' sprintf('output%03d',outputs(1)) '.mat']);
ty_on_EN = ty_trans_theta_on;
ty_off_EN = ty_trans_theta_off;
for i = 2:length(outputs)
    load([base 'trans5N_' model '_' sprintf('output%03d',outputs(i)) '.mat']);
    ty_on_EN = cat(3,ty_on_EN,ty_trans_theta_on);
    ty_off_EN = cat(3,ty_off_EN,ty_trans_theta_off);
end

% La Nina:
model = 'pnEXP2_restart000_windstress';
outputs = [0:3];
load([base 'trans5N_' model '_' sprintf('output%03d',outputs(1)) '.mat']);
ty_on_LN = ty_trans_theta_on;
ty_off_LN = ty_trans_theta_off;
for i = 2:length(outputs)
    load([base 'trans5N_' model '_' sprintf('output%03d',outputs(i)) '.mat']);
    ty_on_LN = cat(3,ty_on_LN,ty_trans_theta_on);
    ty_off_LN = cat(3,ty_off_LN,ty_trans_theta_off);
end

ndays = [31 28 31 30 31 30 31 31 30 31 30 31]';
Re = 6371000;
dx = 2*pi*Re*(0.25)/360*cos(5/180*pi);
dT = diff(T);
dT = dT(1);

% Longitudinal structure of control:
[X,Y] = ndgrid(lonv,T);
figure;
caxs = [-100 -2.5:0.05:2.5 100];
subplot(2,1,1);
contourf(X,Y,monmean(ty_off_CT(:,:,7:12),3,ndays(7:12))/dx/dT,caxs,'linestyle','none');
xlim([-240 -80]);
ylim([10 32]);
caxis([-2.5 2.5]);
cb = colorbar;
ylabel(cb,'m$^2$ s$^{-1}$ $^\circ$ C$^{-1}$');
title(['Monthly-mean transport through $5^\circ$N, July-December ' ...
       'MOM025 Control']);
xlabel('Longitude ($^\circ$E)');
ylabel('Temperature ($^\circ$C)');
subplot(2,1,2);
contourf(X,Y,monmean(ty_on_CT(:,:,7:12)-ty_off_CT(:,:,7:12),3,ndays(7:12))/dx/dT,caxs,'linestyle','none');
xlim([-240 -80]);
ylim([10 32]);
caxis([-2.5 2.5]);
cb = colorbar;
ylabel(cb,'m$^2$ s$^{-1}$ $^\circ$ C$^{-1}$');
title(['Sub-monthly transport through $5^\circ$N, July-December ' ...
       'MOM025 Control']);
xlabel('Longitude ($^\circ$E)');
ylabel('Temperature ($^\circ$C)');
colormap(redblue);

[tmp lnind1] = min(abs(lonv+150));
[tmp lnind2] = min(abs(lonv+120));

% $$$ [X,Y] = ndgrid(time,T);
[X,Y] = ndgrid(1:(4*12),T);
figure;
caxs = [-100 -5:0.1:5 100];
subplot(2,1,1);
contourf(X,Y,squeeze(nansum(ty_off_LN(lnind1:lnind2,:,:),1)/1e6/dT)',caxs,'linestyle','none');
% $$$ xlim([0.5 12.5]);
xlim([0.5 24.5]);
ylim([10 32]);
caxis([-5 5]);
cb = colorbar;
ylabel(cb,'Sv $^\circ$ C$^{-1}$');
title(['Monthly-mean transport through $5^\circ$N, $150^\circ$W-$120^\circ$W, ' ...
       'MOM025 Control']);
xlabel('Month');
ylabel('Temperature ($^\circ$C)');
subplot(2,1,2);
contourf(X,Y,squeeze(nansum(ty_on_LN(lnind1:lnind2,:,:)-ty_off_LN(lnind1:lnind2,:,:),1)/1e6/dT)',caxs,'linestyle','none');
% $$$ xlim([0.5 12.5]);
xlim([0.5 24.5]);
ylim([10 32]);
caxis([-5 5]);
cb = colorbar;
ylabel(cb,'Sv $^\circ$ C$^{-1}$');
title(['Sub-monthly transport through $5^\circ$N, $150^\circ$W-$120^\circ$W, ' ...
       'MOM025 Control']);
xlabel('Month');
ylabel('Temperature ($^\circ$C)');
colormap(redblue)

[tmp lnind1] = min(abs(lonv+150));
[tmp lnind2] = min(abs(lonv+120));
[tmp lnind1] = min(abs(lonv+180));
[tmp lnind2] = min(abs(lonv+100));
% $$$ lnind1 = 1;
% $$$ lnind2 = length(lonv);
TEMP = 20;
[tmp ind] = min(abs(T-TEMP));
filt = 3;
figure;
plot(1:(4*12),filter_field(repmat(squeeze(nansum(nansum(ty_off_CT(lnind1:lnind2,ind:end,:),1),2)/1e6)',[1 4]),filt,'-t'),'--k');
hold on;
plot(1:(4*12),filter_field(repmat(squeeze(nansum(nansum(ty_on_CT(lnind1:lnind2,ind:end,:)-ty_off_CT(lnind1:lnind2,ind:end,:),1),2)/1e6)',[1 4]),filt,'-t'),':k');
% $$$ plot(1:(4*12),filter_field(repmat(squeeze(nansum(nansum(ty_on_CT(lnind1:lnind2,ind:end,:),1),2)/1e6)',[1 4]),filt,'-t'),'-k');
plot(1:(4*12),filter_field(squeeze(nansum(nansum(ty_off_EN(lnind1:lnind2,ind:end,:),1),2)/1e6)',filt,'-t'),'--r');
plot(1:(4*12),filter_field(squeeze(nansum(nansum(ty_on_EN(lnind1:lnind2,ind:end,:)-ty_off_EN(lnind1:lnind2,ind:end,:),1),2)/1e6)',filt,'-t'),':r');
% $$$ plot(1:(4*12),filter_field(squeeze(nansum(nansum(ty_on_EN(lnind1:lnind2,ind:end,:),1),2)/1e6)',filt,'-t'),'-r');
plot(1:(4*12),filter_field(squeeze(nansum(nansum(ty_off_LN(lnind1:lnind2,ind:end,:),1),2)/1e6)',filt,'-t'),'--b');
plot(1:(4*12),filter_field(squeeze(nansum(nansum(ty_on_LN(lnind1:lnind2,ind:end,:)-ty_off_LN(lnind1:lnind2,ind:end,:),1),2)/1e6)',filt,'-t'),':b');
% $$$ plot(1:(4*12),filter_field(squeeze(nansum(nansum(ty_on_LN(lnind1:lnind2,ind:end,:),1),2)/1e6)',filt,'-t'),'-b');

% $$$ legend('Control Mean','Control Eddy','Control Total', ...
% $$$        'El Nino Mean','El Nino Eddy','El Nino Total', ...
% $$$        'La Nina Mean','La Nina Eddy','La Nina Total');
legend('Control Mean','Control Eddy', ...
       'El Nino Mean','El Nino Eddy', ...
       'La Nina Mean','La Nina Eddy');
xlabel('Month');
ylabel('Transport through $5^\circ$N (Sv)');
title(['Mean and Eddy transports $180^\circ$-$100^\circ$W above ' ...
       num2str(TEMP) ' $^\circ$C']);
xlim([1 24]);
 
