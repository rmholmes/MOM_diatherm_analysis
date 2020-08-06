% This script makes plots of the spatial structure of the
% temperature-integrated heat budget in the MOM simulations down to a
% particular temperature (or to -3C for the vertically-integrated
% budget) to check closures and look at the roles of skew-diffusion
% etc.

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/ncdata/';

% $$$ % MOM025 Control:
% $$$ wname = [base 'ocean_wmass.out111-120.mom025cntrl.ncra.nc'];
% $$$ gname = [base 'ocean_grid.mom025.nc'];
% $$$ tsc = 1e9;
% $$$ haveRedi = 0;
% $$$ haveGM = 0;

% ACCCESS-OM2-025-RG:
wname = [base 'ocean_wmass.aom2rg.out7683.ncra.nc'];
% $$$ wname = [base 'ocean_wmass.aom2.out083lap.ncra.nc'];
%wname = [base 'ocean_wmass.aom2.out084bih.ncra.nc'];
gname = [base 'ocean_grid.aom2025.nc'];
tsc = 1;
haveRedi = 1;
haveGM = 1;

area = ncread(gname,'area_t');
lon = ncread(gname,'geolon_t');
lat = ncread(gname,'geolat_t');
% Use the following fro masking over land:
%load(['/srv/ccrc/data03/z3500785/mom/mat_data/MOM025_kb3seg_output090_BaseVars.mat'],'lon','lat');

Te = ncread(wname,'neutralrho_edges');
T = ncread(wname,'neutral');
TL = length(T);
[xL,yL] = size(area);
Cp = 3992.10322329649; % J kg-1 degC-1
rho0 = 1035; % kgm-3

% sum down to this temperature:
Tmin = -3;%15;
[tmp Tmini] = min(abs(Te-Tmin));

ovars = {'temp_tendency_on_nrho','temp_advection_on_nrho',...
          'temp_submeso_on_nrho'};
diffvars = {'temp_vdiffuse_diff_cbt_on_nrho',...
          'temp_nonlocal_KPP_on_nrho','temp_vdiffuse_sbc_on_nrho',...
          'sfc_hflux_pme_on_nrho','frazil_on_nrho','temp_eta_smooth_on_nrho',...
          'sw_heat_on_nrho','temp_rivermix_on_nrho'};
tx_vars = {'tx_trans_nrho'};
ty_vars = {'ty_trans_nrho'};
qx_vars = {'temp_xflux_adv_on_nrho'};
qy_vars = {'temp_yflux_adv_on_nrho'};
qxS_vars = {'temp_xflux_submeso_on_nrho'};
qyS_vars = {'temp_yflux_submeso_on_nrho'};
Js_var = {'mass_pmepr_on_nrho'};
ten_vars = {'dHdt','dVdt'};

% Redi/GM:
if (haveRedi)
diffvars = {diffvars{:},'temp_vdiffuse_k33_on_nrho','neutral_diffusion_on_nrho_temp'};
end
if (haveGM)
qxS_vars = {qxS_vars{:},'temp_xflux_gm_on_nrho'};
qyS_vars = {qyS_vars{:},'temp_yflux_gm_on_nrho'};
ovars = {ovars{:},'neutral_gm_on_nrho_temp'};
end

allvars = {ovars{:},diffvars{:},tx_vars{:},ty_vars{:},qx_vars{:},qy_vars{:}, ...
           Js_var{:},ten_vars{:},qxS_vars{:},qyS_vars{:}};

for ii=1:length(allvars)
    eval([allvars{ii} ' = nansum(ncread(wname,''' allvars{ii} ''',[1 1 Tmini 1],[xL yL TL+1-Tmini 1]),3);' ...
                        '']);
end

% Load pre-processed num-mix:
ndifa = ncread(wname,'temp_numdiff_heat_on_nrho',[1 1 Tmini 1],[xL yL 1 ...
                    1]);

% Load skew-diffusion associated terms:
if (Tmini ~= 1)
JxS = avg(ncread(wname,'tx_trans_nrho_submeso', ...
                          [1 1 Tmini-1 1],[xL yL 2 1])*tsc/rho0,3);
JyS = avg(ncread(wname,'ty_trans_nrho_submeso', ...
                          [1 1 Tmini-1 1],[xL yL 2 1])*tsc/rho0,3);
if (haveGM)
    JxS = JxS + avg(ncread(wname,'tx_trans_nrho_gm', ...
                          [1 1 Tmini-1 1],[xL yL 2 1])*tsc/rho0,3);
    JyS = JyS + avg(ncread(wname,'ty_trans_nrho_gm', ...
                          [1 1 Tmini-1 1],[xL yL 2 1])*tsc/rho0,3);
end
else 
    JxS = zeros(size(tx_trans_nrho));
    JyS = zeros(size(ty_trans_nrho));
end

QxS = zeros(size(tx_trans_nrho));
for ii=1:length(qxS_vars)
    eval(['QxS = QxS + ' qxS_vars{ii} ';']);
end
QyS = zeros(size(ty_trans_nrho));
for ii=1:length(qyS_vars)
    eval(['QyS = QyS + ' qyS_vars{ii} ';']);
end

% calculate convergence:
Jx = zeros(size(tx_trans_nrho));
for ii=1:length(tx_vars)
    eval(['Jx = Jx + ' tx_vars{ii} '*tsc/rho0;']);
end
Jy = zeros(size(ty_trans_nrho));
for ii=1:length(ty_vars)
    eval(['Jy = Jy + ' ty_vars{ii} '*tsc/rho0;']);
end
Qx = zeros(size(tx_trans_nrho));
for ii=1:length(qx_vars)
    eval(['Qx = Qx + ' qx_vars{ii} ';']);
end
Qy = zeros(size(ty_trans_nrho));
for ii=1:length(qy_vars)
    eval(['Qy = Qy + ' qy_vars{ii} ';']);
end

rcT = rho0*Cp*Te(Tmini);

JI = zeros(size(temp_tendency_on_nrho));
QI = JI;
JI(2:end,2:end) = JI(2:end,2:end) + (Jx(1:(end-1),2:end) - Jx(2:end,2:end) ...
                                     +Jy(2:end,1:(end-1)) - Jy(2:end,2:end))./area(2:end,2:end);
JI(1,2:end) = JI(1,2:end) + (Jx(end,2:end) - Jx(1,2:end) ...
                             +Jy(1,1:(end-1)) - Jy(1,2:end))./area(1,2:end);        
QI(2:end,2:end) = QI(2:end,2:end) + (Qx(1:(end-1),2:end) - Qx(2:end,2:end) ...
                                     +Qy(2:end,1:(end-1)) - Qy(2:end,2:end))./area(2:end,2:end);
QI(1,2:end) = QI(1,2:end) + (Qx(end,2:end) - Qx(1,2:end) ...
                             +Qy(1,1:(end-1)) - Qy(1,2:end))./area(1,2:end);        

JIS = zeros(size(temp_tendency_on_nrho));
QIS = JIS;
JIS(2:end,2:end) = JIS(2:end,2:end) + (JxS(1:(end-1),2:end) - JxS(2:end,2:end) ...
                                     +JyS(2:end,1:(end-1)) - JyS(2:end,2:end))./area(2:end,2:end);
JIS(1,2:end) = JIS(1,2:end) + (JxS(end,2:end) - JxS(1,2:end) ...
                             +JyS(1,1:(end-1)) - JyS(1,2:end))./area(1,2:end);        
QIS(2:end,2:end) = QIS(2:end,2:end) + (QxS(1:(end-1),2:end) - QxS(2:end,2:end) ...
                                     +QyS(2:end,1:(end-1)) - QyS(2:end,2:end))./area(2:end,2:end);
QIS(1,2:end) = QIS(1,2:end) + (QxS(end,2:end) - QxS(1,2:end) ...
                             +QyS(1,1:(end-1)) - QyS(1,2:end))./area(1,2:end);        

JS = mass_pmepr_on_nrho./area/rho0;
dVdt = dVdt*1e9/rho0./area;
dHdt = dHdt./area;

dift = zeros(size(temp_vdiffuse_diff_cbt_on_nrho));
for ii = 1:length(diffvars)
    eval(['dift = dift + ' diffvars{ii} ';']);
end

xvec = 760:1:980;
yvec = 560:1:730; % +15 -> +50
xvec = 1:4:xL;
yvec = 1:4:yL;

poss = [0.0886    0.69   0.3375    0.2581; ...
        0.4503    0.69   0.3375    0.2581; ...
        0.0886    0.39   0.3375    0.2581; ...
        0.4503    0.39   0.3375    0.2581; ...
        0.0886    0.09   0.3375    0.2581; ...
        0.4503    0.09   0.3375    0.2581];

% $$$ %%% Checking skew-diffusion:
% $$$ % $$$ subndifa = 1;
% $$$ % $$$ caxs1 = [-1 1];
% $$$ % $$$ caxsn = [-1e-4 1e-4];
% $$$ subndifa = 1;
% $$$ caxs1 = [-125 125];
% $$$ caxsn = [-125 125];
% $$$ figure;
% $$$ set(gcf,'defaulttextfontsize',15);
% $$$ set(gcf,'defaultaxesfontsize',15);
% $$$ ii = 1;
% $$$ % 1) standard (QIS yes but no JIS)
% $$$ Gt = (dVdt - JI - JS)*rcT;
% $$$ ndif = dHdt - Gt - dift - QI - QIS;
% $$$ subplot(2,2,ii);
% $$$ var = ndif;
% $$$ var(var == 0) = NaN;
% $$$ pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
% $$$ title('$\mathcal{I}$ not including $\Psi_{GM}$');
% $$$ ylabel('Latitude ($^\circ$N)');
% $$$ caxis(caxs1);
% $$$ % $$$ set(gca,'xticklabel',[]);
% $$$ xlabel('Longitude ($^\circ$E)');
% $$$ colormap(redblue);
% $$$ colorbar off;
% $$$ %set(gca,'Position',poss(ii,:));
% $$$ set(gca,'color','k');
% $$$ cb = colorbar;
% $$$ ylabel(cb,'Wm$^{-2}$');
% $$$ % Check several different ways:
% $$$ % 2) including vol flux (QIS and JIS)
% $$$ Gt = (dVdt - JI - JS - JIS)*rcT;
% $$$ ndif = dHdt - Gt - dift - QI - QIS;
% $$$ subplot(2,2,ii+1);
% $$$ var = ndif;
% $$$ if (subndifa); var = var-ndifa; end;
% $$$ var(var == 0) = NaN;
% $$$ pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
% $$$ title('$\mathcal{I}$ including $\Psi_{GM}$');
% $$$ ylabel('Latitude ($^\circ$N)');
% $$$ caxis(caxsn);
% $$$ %set(gca,'xticklabel',[]);
% $$$ xlabel('Longitude ($^\circ$E)');
% $$$ colormap(redblue);
% $$$ %set(gca,'Position',poss(ii+1,:));
% $$$ set(gca,'color','k');
% $$$ cb = colorbar;
% $$$ ylabel(cb,'Wm$^{-2}$');
% $$$ %colorbar off;
% $$$ 
% $$$ % 3) Purely by 3D convergence (no QIS or JIS, but add submeso and
% $$$ % gm 3D convergences to dift).
% $$$ Gt = (dVdt - JI - JS)*rcT;
% $$$ ndif = dHdt - Gt - dift - QI - temp_submeso_on_nrho - neutral_gm_on_nrho_temp;
% $$$ subplot(2,2,ii+2);
% $$$ var = ndif;
% $$$ if (subndifa); var = var-ndifa; end;
% $$$ var(var == 0) = NaN;
% $$$ pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
% $$$ title('$\mathcal{I}$ using 3D GM heat flux convergence');
% $$$ ylabel('Latitude ($^\circ$N)');
% $$$ caxis(caxsn);
% $$$ %set(gca,'xticklabel',[]);
% $$$ xlabel('Longitude ($^\circ$E)');
% $$$ colormap(redblue);
% $$$ %set(gca,'Position',poss(ii+2,:));
% $$$ set(gca,'color','k');
% $$$ cb = colorbar;
% $$$ ylabel(cb,'Wm$^{-2}$');
% $$$ %colorbar off;
% $$$ subplot(2,2,ii+3);
% $$$ Gt = (dVdt - JI - JS - JIS)*rcT;
% $$$ ndif = dHdt - Gt - dift - QI - temp_submeso_on_nrho - neutral_gm_on_nrho_temp;
% $$$ var = ndif;%-ndifa;
% $$$ if (subndifa); var = var-ndifa; end;
% $$$ var(var == 0) = NaN;
% $$$ pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
% $$$ title('$\mathcal{I}$ using 3D GM heat flux convergence and $\Psi_{GM}$');
% $$$ % $$$ title('ndif 3D convergence');
% $$$ ylabel('Latitude ($^\circ$N)');
% $$$ caxis(caxsn);
% $$$ %set(gca,'xticklabel',[]);
% $$$ xlabel('Longitude ($^\circ$E)');
% $$$ colormap(redblue);
% $$$ %set(gca,'Position',poss(ii+3,:));
% $$$ set(gca,'color','k');
% $$$ cb = colorbar;
% $$$ ylabel(cb,'Wm$^{-2}$');
% $$$ %colorbar off;
% $$$ % 4) Any others? To maintain the meridional fluxes properly I need
% $$$ % to include one and only one of QIS or the 3D convergences - but
% $$$ % what about the volume budget? The volume budget is made up -
% $$$ % through G? Where does this come from?
% $$$ 
% $$$ % Volume budget:
% $$$ caxs = [-1 1]*1e-7;
% $$$ figure;
% $$$ set(gcf,'defaulttextfontsize',15);
% $$$ set(gcf,'defaultaxesfontsize',15);
% $$$ ii = 1;
% $$$ subplot(2,2,ii);
% $$$ var = dVdt;
% $$$ var(var == 0) = NaN;
% $$$ pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
% $$$ title('dV/dt');%$\mathcal{I}$ not including $\Psi_{GM}$');
% $$$ ylabel('Latitude ($^\circ$N)');
% $$$ caxis(caxs);
% $$$ % $$$ set(gca,'xticklabel',[]);
% $$$ xlabel('Longitude ($^\circ$E)');
% $$$ colormap(redblue);
% $$$ colorbar off;
% $$$ %set(gca,'Position',poss(ii,:));
% $$$ set(gca,'color','k');
% $$$ cb = colorbar;
% $$$ ylabel(cb,'ms$^{-1}$');
% $$$ subplot(2,2,ii+1);
% $$$ var = JI;
% $$$ var(var==0) = NaN;
% $$$ pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
% $$$ title('$J_I$');%$\mathcal{I}$ including $\Psi_{GM}$');
% $$$ ylabel('Latitude ($^\circ$N)');
% $$$ caxis(caxs);
% $$$ %set(gca,'xticklabel',[]);
% $$$ xlabel('Longitude ($^\circ$E)');
% $$$ colormap(redblue);
% $$$ %set(gca,'Position',poss(ii+1,:));
% $$$ set(gca,'color','k');
% $$$ cb = colorbar;
% $$$ ylabel(cb,'ms$^{-1}$');
% $$$ subplot(2,2,ii+2);
% $$$ var = JS;%-ndifa;
% $$$ var(var==0) = NaN;
% $$$ pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
% $$$ title('$J_S$');%$\mathcal{I}$ using 3D GM heat flux convergence');
% $$$ ylabel('Latitude ($^\circ$N)');
% $$$ caxis(caxs);
% $$$ %set(gca,'xticklabel',[]);
% $$$ xlabel('Longitude ($^\circ$E)');
% $$$ colormap(redblue);
% $$$ %set(gca,'Position',poss(ii+2,:));
% $$$ set(gca,'color','k');
% $$$ cb = colorbar;
% $$$ ylabel(cb,'ms$^{-1}$');
% $$$ %colorbar off;
% $$$ subplot(2,2,ii+3);
% $$$ var = dVdt-JI-JS;%ndif;%-ndifa;
% $$$ var(var==0) = NaN;
% $$$ pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
% $$$ title('$dV/dt - J_I - J_S$');%$\mathcal{I}$ using 3D GM heat flux convergence and $\Psi_{GM}$');
% $$$ % $$$ title('ndif 3D convergence');
% $$$ ylabel('Latitude ($^\circ$N)');
% $$$ caxis(caxs);
% $$$ %set(gca,'xticklabel',[]);
% $$$ xlabel('Longitude ($^\circ$E)');
% $$$ colormap(redblue);
% $$$ %set(gca,'Position',poss(ii+3,:));
% $$$ set(gca,'color','k');
% $$$ cb = colorbar;
% $$$ ylabel(cb,'ms$^{-1}$');

    
%%% Vertically-integrated vol/heat and internal heat budget residuals:
xvec = 1:1:xL;
yvec = 1:1:yL;
caxs = [-5 5];
caxsV = [-2 2]*1e-7;
figure;
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);
ii = 2;
subplot(2,2,ii);
var = dHdt-dift-QI-QIS;
var(var == 0) = NaN;
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
title('$\frac{\partial\mathcal{H}}{\partial t}-\mathcal{M}-\mathcal{R}-\mathcal{F}-\mathcal{P}-\mathcal{Q}$');
% $$$ ylabel('Latitude ($^\circ$N)');
caxis(caxs);
set(gca,'xticklabel',[]);
colormap(redblue);
% $$$ colorbar off;
cb = colorbar;
ylabel(cb,'Wm$^{-2}$');
set(gca,'Position',poss(ii,:));
set(gca,'color','k');
subplot(2,2,ii+2);
var = dVdt - JI - JS;
var(var == 0) = NaN;
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
title('$\frac{\partial\mathcal{V}}{\partial t}-\mathcal{J}-\Psi$');
caxis(caxsV);
% $$$ ylabel('Latitude ($^\circ$N)');
xlabel('Longitude ($^\circ$E)');
% $$$ set(gca,'xticklabel',[]);
colormap(redblue);
set(gca,'Position',poss(ii+2,:));
set(gca,'color','k');
cb = colorbar;
ylabel(cb,'ms$^{-1}$');
% $$$ colorbar off;
% $$$ subplot(2,2,ii+4);
% $$$ Gt = dVdt - JI - JS;
% $$$ ndif = dHdt - Gt - dift - QI - QIS;
% $$$ var(var == 0) = NaN;
% $$$ pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
% $$$ title('Residual $\mathcal{I}$');
% $$$ ylabel('Latitude ($^\circ$N)');
% $$$ caxis(caxs);
% $$$ colormap(redblue);
% $$$ set(gca,'Position',poss(ii+4,:));
% $$$ xlabel('Longitude ($^\circ$E)');
% $$$ set(gca,'color','k');
% $$$ colorbar off;

% $$$ subplot(3,3,2);
% $$$ var = Gt;
% $$$ pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
% $$$ title('-(dVdt-JI-JS)*rho0*Cp*(-3C)');
% $$$ subplot(3,3,3);
% $$$ var = dHdt;
% $$$ pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
% $$$ title('dHdt');
% $$$ 
% $$$ 
% $$$ % $$$ for ii=1:length(vars)
% $$$ % $$$     subplot(4,4,ii);
% $$$ % $$$     eval(['pcolPlot(lon(xvec,yvec),lat(xvec,yvec),' vars{ii} '(xvec,yvec));']);
% $$$ % $$$     title(strrep(vars{ii},'_',' '));
% $$$ % $$$ end
% $$$ subplot(3,3,1);
% $$$ var = dift;
% $$$ pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
% $$$ title('Surface Forcing + vertical mixing');
% $$$ subplot(3,3,2);
% $$$ var = Gt;
% $$$ pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
% $$$ title('-(dVdt-JI-JS)*rho0*Cp*(-3C)');
% $$$ subplot(3,3,3);
% $$$ var = dHdt;
% $$$ pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
% $$$ title('dHdt');
% $$$ subplot(3,3,4);
% $$$ var = QI;
% $$$ pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
% $$$ title('QI');
% $$$ subplot(3,3,5);
% $$$ var = ndif;
% $$$ pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
% $$$ title('RES');
% $$$ subplot(3,3,6);
% $$$ var = dHdt-dift-QI;
% $$$ pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
% $$$ title('dHdt-dift-QI');
% $$$ subplot(3,3,7);
% $$$ var = ndifa;
% $$$ pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
% $$$ title('ndifa');%dHdt-dift-QI');
