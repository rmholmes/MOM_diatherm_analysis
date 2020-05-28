% This script makes plots of the vertically-/temperature-integrated
% heat budget in the MOM simulations to check closure.

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/ncdata/';

% MOM025:
wname = [base 'ocean_wmass.out111-120.mom025cntrl.ncra.nc'];
tsc = 1e9;

% ACCCESS-OM2:
wname = [base 'ocean_wmass.aom2025.ncra.nc'];
tsc = 1;


Te = ncread(wname,'neutralrho_edges');%')
Cp = 3992.10322329649; % J kg-1 degC-1
rho0 = 1035; % kgm-3

ovars = {'temp_tendency_on_nrho','temp_advection_on_nrho',...
          'temp_submeso_on_nrho'};
diffvars = {'temp_vdiffuse_diff_cbt_on_nrho',...
          'temp_nonlocal_KPP_on_nrho','temp_vdiffuse_sbc_on_nrho',...
          'sfc_hflux_pme_on_nrho','frazil_on_nrho','temp_eta_smooth_on_nrho',...
          'sw_heat_on_nrho','temp_rivermix_on_nrho'};
tx_vars = {'tx_trans_nrho'};
ty_vars = {,'ty_trans_nrho'};
qx_vars = {'temp_xflux_adv_on_nrho',...
          'temp_xflux_submeso_on_nrho'};
qy_vars = {'temp_yflux_adv_on_nrho',...
          'temp_yflux_submeso_on_nrho'};
Js_var = {'mass_pmepr_on_nrho'};
ten_vars = {'dHdt','dVdt'};

% $$$ % Redi/GM:
% $$$ diffvars = {diffvars{:},'temp_vdiffuse_k33_on_nrho','neutral_diffusion_on_nrho_temp'};
% $$$ qx_vars = {qx_vars{:},'temp_xflux_gm_on_nrho'};
% $$$ qy_vars = {qy_vars{:},'temp_yflux_gm_on_nrho'};


allvars = {ovars{:},diffvars{:},tx_vars{:},ty_vars{:},qx_vars{:},qy_vars{:}, ...
           Js_var{:},ten_vars{:}};
for ii=1:length(allvars)
    eval([allvars{ii} ' = sum(ncread(wname,''' allvars{ii} '''),3);' ...
                        '']);
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

area = ncread([base 'ocean_grid_mom025.nc'],'area_t');
load(['/srv/ccrc/data03/z3500785/mom/mat_data/MOM025_kb3seg_output090_BaseVars.mat'],'lon','lat');
% $$$ lon = ncread([base 'ocean_grid_mom025.nc'],'geolon_t');
% $$$ lat = ncread([base 'ocean_grid_mom025.nc'],'geolat_t');
[xL,yL] = size(area);
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

JI = JI*rho0*Cp*Te(1);
JS = (mass_pmepr_on_nrho./area/rho0)*rho0*Cp*Te(1);
dVdt = (dVdt*1e9/rho0./area)*rho0*Cp*Te(1);
dHdt = dHdt./area;

ndifa = ncread(wname,'temp_numdiff_heat_on_nrho',[1 1 1 1],[xL yL 1 ...
                    1]);
dift = zeros(size(temp_vdiffuse_diff_cbt_on_nrho));
for ii = 1:length(diffvars)
    eval(['dift = dift + ' diffvars{ii} ';']);
end

Gt = dVdt - JI - JS;
ndif = dHdt - Gt - dift - QI;

vars = {diffvars{:},'JI','JS','dVdt','dHdt','QI','dift','ndif','ndifa'};
xvec = 1:1:xL;
yvec = 1:1:yL;

    poss = [0.0886    0.69   0.3375    0.2581; ...
            0.4503    0.69   0.3375    0.2581; ...
            0.0886    0.39   0.3375    0.2581; ...
            0.4503    0.39   0.3375    0.2581; ...
            0.0886    0.09   0.3375    0.2581; ...
            0.4503    0.09   0.3375    0.2581];

figure;
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);
ii = 1;
subplot(3,2,ii);
var = dHdt-dift-QI;
var(var == 0) = NaN;
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
title('$\frac{\partial\mathcal{H}}{\partial t}-\mathcal{M}-\mathcal{R}-\mathcal{F}-\mathcal{P}-\mathcal{Q}_I)$');
ylabel('Latitude ($^\circ$N)');
caxis([-1 1]);
set(gca,'xticklabel',[]);
colormap(redblue);
colorbar off;
set(gca,'Position',poss(ii,:));
set(gca,'color','k');
colorbar off;
subplot(3,2,ii+2);
var = -Gt;
var(var == 0) = NaN;
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
title('$-\rho_0C_p\Theta\left(\frac{\partial\mathcal{V}}{\partial t}-\mathcal{J}_S-\mathcal{J}_A\right)$');
caxis([-1 1]);
set(gca,'xticklabel',[]);
colormap(redblue);
set(gca,'Position',poss(ii+2,:));
set(gca,'color','k');
colorbar off;
subplot(3,2,ii+4);
var = ndif;
var(var == 0) = NaN;
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
title('Residual $\mathcal{I}$');
ylabel('Latitude ($^\circ$N)');
caxis([-1 1]);
colormap(redblue);
set(gca,'Position',poss(ii+4,:));
xlabel('Longitude ($^\circ$E)');
set(gca,'color','k');
colorbar off;

subplot(3,3,2);
var = Gt;
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
title('-(dVdt-JI-JS)*rho0*Cp*(-3C)');
subplot(3,3,3);
var = dHdt;
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
title('dHdt');


% $$$ for ii=1:length(vars)
% $$$     subplot(4,4,ii);
% $$$     eval(['pcolPlot(lon(xvec,yvec),lat(xvec,yvec),' vars{ii} '(xvec,yvec));']);
% $$$     title(strrep(vars{ii},'_',' '));
% $$$ end
subplot(3,3,1);
var = dift;
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
title('Surface Forcing + vertical mixing');
subplot(3,3,2);
var = Gt;
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
title('-(dVdt-JI-JS)*rho0*Cp*(-3C)');
subplot(3,3,3);
var = dHdt;
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
title('dHdt');
subplot(3,3,4);
var = QI;
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
title('QI');
subplot(3,3,5);
var = ndif;
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
title('RES');
subplot(3,3,6);
var = dHdt-dift-QI;
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
title('dHdt-dift-QI');
subplot(3,3,7);
var = ndifa;
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),var(xvec,yvec));
title('ndifa');%dHdt-dift-QI');
