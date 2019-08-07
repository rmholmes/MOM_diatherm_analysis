% This script extracts specified data from MOM025 and MOM01 runs 

baseL = '/short/e14/rmh561/mom/archive/';
% $$$ baseL = '/g/data/e14/rmh561/mom/';
% $$$ baseL = '/short/e14/rmh561/access-om2/archive/';
% $$$ baseL = '/srv/ccrc/data03/z3500785/';
% $$$ types = {'kds50','gfdl50','kds75','kds100','kds135'};

model = 'MOM025_kb3seg';
baseD = [baseL 'MOM_HeatDiag_kb3seg/'];

outD = [baseD];

output = 95;
% $$$ for output=86:90
    if (output==0)
        restart=0;
    else
        restart = output-1;
    end

%% file-names and grid properties:
base = [baseD sprintf('output%03d/',output)];
hname = [base 'ocean_heat.nc'];
if (strfind(baseD,'01'))
    fname = [base 'ocean_month.nc'];
    m3name = [base 'ocean.nc'];
else
    fname = [base 'ocean.nc'];
end
gname = [base 'ocean_grid.nc'];
sname = [base 'ocean_snap.nc'];
wname = [base 'ocean_wmass.nc'];
         
lon = ncread(gname,'geolon_t');lat = ncread(gname,'geolat_t');
area = ncread(gname,'area_t');[xL,yL] = size(lon);
lonu = ncread(gname,'geolon_c');latu = ncread(gname,'geolat_c');

z = ncread(fname,'st_ocean');zL = length(z);
zw = ncread(fname,'sw_ocean');

time = ncread(fname,'time');
tL = length(time);

Cp = 3992.1; % J kg-1 degC-1
rho0 = 1035; % kgm-3

% $$$ T = ncread(wname,'neutral');
% $$$ Te = ncread(wname,'neutralrho_edges');
% $$$ TL = length(T);dT = T(2)-T(1);
% $$$ 
% $$$ %% Get lat-depth slice of variables:
% $$$ for lonsl=[-110 -140]
% $$$ [tmp lt1] = min(abs(lat(1,:)+20));
% $$$ [tmp lt2] = min(abs(lat(1,:)-20));
% $$$ 
% $$$ [tmp lnind] = min(abs(lon(:,round(mean([lt1 lt2])))-lonsl));
% $$$ 
% $$$ temp = squeeze(ncread(fname,'temp',[lnind lt1 1 1],[1 lt2-lt1+1 zL tL]));
% $$$ 
% $$$ if (strfind(baseD,'01'))
% $$$     u = squeeze(ncread(m3name,'u',[lnind lt1 1 1],[1 lt2-lt1+1 zL 1]));
% $$$     v = squeeze(ncread(m3name,'v',[lnind lt1 1 1],[1 lt2-lt1+1 zL 1]));
% $$$ else
% $$$     u = squeeze(ncread(fname,'u',[lnind lt1 1 1],[1 lt2-lt1+1 zL tL]));
% $$$     v = squeeze(ncread(fname,'v',[lnind lt1 1 1],[1 lt2-lt1+1 zL tL]));
% $$$ end
% $$$ 
% $$$ kappa = squeeze(ncread(fname,'diff_cbt_t',[lnind lt1 1 1],[1 lt2-lt1+1 zL tL]));
% $$$ taux = squeeze(ncread(fname,'tau_x',[lnind lt1 1],[1 lt2-lt1+1 tL]));
% $$$ tauy = squeeze(ncread(fname,'tau_y',[lnind lt1 1],[1 lt2-lt1+1 tL]));
% $$$ mld = squeeze(ncread(fname,'mld',[lnind lt1 1],[1 lt2-lt1+1 tL]));
% $$$ vdif = squeeze(ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[lnind lt1 1 1],[1 lt2-lt1+1 TL tL]));
% $$$ vnlc = squeeze(ncread(wname,'temp_nonlocal_KPP_on_nrho',[lnind lt1 1 1],[1 lt2-lt1+1 TL tL]));
% $$$ ndif = squeeze(ncread(wname,'temp_numdiff_heat_on_nrho',[lnind lt1 1 1],[1 lt2-lt1+1 TL tL]));
% $$$ pmer = squeeze(ncread(wname,'sfc_hflux_pme_on_nrho',[lnind lt1 1 1],[1 lt2-lt1+1 TL tL]) + ...
% $$$                ncread(wname,'temp_rivermix_on_nrho',[lnind lt1 1 1],[1 lt2-lt1+1 TL tL]));
% $$$ sufc = squeeze(ncread(wname,'temp_vdiffuse_sbc_on_nrho',[lnind lt1 1 1],[1 lt2-lt1+1 TL tL]) + ...
% $$$                ncread(wname,'frazil_on_nrho',[lnind lt1 1 1],[1 lt2-lt1+1 TL tL]) + ...
% $$$                ncread(wname,'temp_eta_smooth_on_nrho',[lnind lt1 1 1],[1 lt2-lt1+1 TL tL]) + ...
% $$$                ncread(wname,'sw_heat_on_nrho',[lnind lt1 1 1],[1 lt2-lt1+1 TL tL]));
% $$$ swrd = squeeze(ncread(wname,'sw_heat_on_nrho',[lnind lt1 1 1],[1 lt2-lt1+1 TL tL]));
% $$$ 
% $$$ [Yt,Zt] = ndgrid(lat(lnind,lt1:lt2),z);
% $$$ [Yu,Zu] = ndgrid(latu(lnind,lt1:lt2),z);
% $$$ 
% $$$ name = [outD 'mat_data/' model sprintf('_output%03d',output) '_varsat_' num2str(-lonsl) 'W.mat']
% $$$ save(name,'Yt','Zt','Yu','Zu','temp','kappa','taux','tauy','mld', ...
% $$$      'vdif','vnlc','pmer','sufc','swrd','ndif');
% $$$ end

%% Get lon-depth slices of variables:
% $$$ % Equatorial Pacific:
% $$$ rname = 'EqPM2';
% $$$ [tmp lt1] = min(abs(lat(1,:)+2));
% $$$ [tmp lt2] = min(abs(lat(1,:)-2));
% $$$ [tmp ln1] = min(abs(lon(:,lt1)+240));
% $$$ [tmp ln2] = min(abs(lon(:,lt1)+70));

% Gulf Stream:
rname = 'GulfSt_42pm0p5';
[tmp lt1] = min(abs(lat(1,:)-41.5));
[tmp lt2] = min(abs(lat(1,:)-42.5));
[tmp ln1] = min(abs(lon(:,lt1)+78));
[tmp ln2] = min(abs(lon(:,lt1)+8));

temp = squeeze(mean(ncread(fname,'temp',[ln1 lt1 1 1],[ln2-ln1+1 lt2-lt1+1 zL tL]),2));
if (strfind(baseD,'01'))
    u = squeeze(mean(ncread(m3name,'u',[ln1 lt1 1 1],[ln2-ln1+1 lt2-lt1+1 zL 1]),2));
    v = squeeze(mean(ncread(m3name,'v',[ln1 lt1 1 1],[ln2-ln1+1 lt2-lt1+1 zL 1]),2));
else
    u = squeeze(mean(ncread(fname,'u',[ln1 lt1 1 1],[ln2-ln1+1 lt2-lt1+1 zL tL]),2));
    v = squeeze(mean(ncread(fname,'v',[ln1 lt1 1 1],[ln2-ln1+1 lt2-lt1+1 zL tL]),2));
    u_sq = squeeze(mean(ncread(fname,'u_sq',[ln1 lt1 1 1],[ln2-ln1+1 lt2-lt1+1 zL tL]),2));
    v_sq = squeeze(mean(ncread(fname,'v_sq',[ln1 lt1 1 1],[ln2-ln1+1 lt2-lt1+1 zL tL]),2));
end
kappa = squeeze(mean(ncread(fname,'diff_cbt_t',[ln1 lt1 1 1],[ln2-ln1+1 lt2-lt1+1 zL tL]),2));
taux = squeeze(mean(ncread(fname,'tau_x',[ln1 lt1 1],[ln2-ln1+1 lt2-lt1+1 tL]),2));
tauy = squeeze(mean(ncread(fname,'tau_y',[ln1 lt1 1],[ln2-ln1+1 lt2-lt1+1 tL]),2));
mld = squeeze(mean(ncread(fname,'mld',[ln1 lt1 1],[ln2-ln1+1 lt2-lt1+1 tL]),2));
w = squeeze(mean(ncread(fname,'wt',[ln1 lt1 1 1],[ln2-ln1+1 lt2-lt1+1 zL tL]),2));
w_sq = squeeze(mean(ncread(fname,'wt_sq',[ln1 lt1 1 1],[ln2-ln1+1 lt2-lt1+1 zL tL]),2));
Tdxsq = squeeze(mean(ncread(fname,'temp_dxsq',[ln1 lt1 1 1],[ln2-ln1+1 lt2-lt1+1 zL tL]),2));
Tdysq = squeeze(mean(ncread(fname,'temp_dysq',[ln1 lt1 1 1],[ln2-ln1+1 lt2-lt1+1 zL tL]),2));
Tdzsq = squeeze(mean(ncread(fname,'temp_dzsq',[ln1 lt1 1 1],[ln2-ln1+1 lt2-lt1+1 zL tL]),2));

% On isotherm terms:
vdif = squeeze(mean(ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[ln1 lt1 1 1],[ln2-ln1+1 lt2-lt1+1 TL tL]),2));
ndif = squeeze(mean(ncread(wname,'temp_numdiff_heat_on_nrho',[ln1 lt1 1 1],[ln2-ln1+1 lt2-lt1+1 TL tL]),2));
vnlc = squeeze(mean(ncread(wname,'temp_nonlocal_KPP_on_nrho',[ln1 lt1 1 1],[ln2-ln1+1 lt2-lt1+1 TL tL]),2));
pmer = squeeze(mean(ncread(wname,'sfc_hflux_pme_on_nrho',[ln1 lt1 1 1],[ln2-ln1+1 lt2-lt1+1 TL tL]) + ...
               ncread(wname,'temp_rivermix_on_nrho',[ln1 lt1 1 1],[ln2-ln1+1 lt2-lt1+1 TL tL]),2));
sufc = squeeze(mean(ncread(wname,'temp_vdiffuse_sbc_on_nrho',[ln1 lt1 1 1],[ln2-ln1+1 lt2-lt1+1 TL tL]) + ...
               ncread(wname,'frazil_on_nrho',[ln1 lt1 1 1],[ln2-ln1+1 lt2-lt1+1 TL tL]) + ...
               ncread(wname,'temp_eta_smooth_on_nrho',[ln1 lt1 1 1],[ln2-ln1+1 lt2-lt1+1 TL tL]) + ...
               ncread(wname,'sw_heat_on_nrho',[ln1 lt1 1 1],[ln2-ln1+1 lt2-lt1+1 TL tL]),2));
swrd = squeeze(mean(ncread(wname,'sw_heat_on_nrho',[ln1 lt1 1 1],[ln2-ln1+1 lt2-lt1+1 TL tL]),2));

[Xt,Zt] = ndgrid(mean(lon(ln1:ln2,lt1:lt2),2),z);
[Xu,Zu] = ndgrid(mean(lonu(ln1:ln2,lt1:lt2),2),z);
[Xw,Zw] = ndgrid(mean(lon(ln1:ln2,lt1:lt2),2),zw);

name = [outD 'mat_data/' model sprintf('_output%03d',output) ...
        '_varsat_' rname '.mat']
% Non-isotherm only:
save(name,'Xt','Zt','Xu','Zu','Xw','Zw','temp','mld','u','v','u_sq','v_sq','w','w_sq','Tdxsq','Tdysq','Tdzsq');
% Isotherm only:
% $$$ save(name,'Xt','Zt','Xu','Zu','temp','kappa','taux','tauy','mld', ...
% $$$      'vdif','vnlc','pmer','sufc','swrd','ndif');
end

