% This script extracts specified data from MOM025 and MOM01 runs 

% $$$ model = 'MOM01';
model = 'MOM025_kb3seg';
baseD = '/short/e14/rmh561/mom/archive/MOM_HeatDiag_kb3seg/'; %Data Directory.
% $$$ baseD = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/'; %Data Directory.
% $$$ baseD = '/g/data/e14/rmh561/MOM01_HeatDiag/';
% $$$ baseD = '/g/data/e14/rmh561/MOM_HeatDiag/';
% $$$ baseD = '/short/e14/mv7494/mom_perturbations/EXP1_and_EXP2_restart000_windstress/archive/';
outD = '/short/e14/rmh561/mom/archive/MOM_HeatDiag_kb3seg/';

for output=75:79
% $$$ for output=[0 1 2 3 266 267 268 269]
    if (output==0)
        restart=0;
    else
        restart = output-1;
    end

%% file-names and grid properties:
base = [baseD sprintf('output%03d/',output)];
if (output==0)
    baser = '/short/e14/mv7494/mom_control/archive/restart000/';
else
    baser = [baseD sprintf('restart%03d/',restart)];
end
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
rnameT = [baser 'ocean_temp_salt.res.nc'];
rnameZ = [baser 'ocean_thickness.res.nc'];
rnametime = [baser 'coupler.res'];
         
lon = ncread(hname,'geolon_t');lat = ncread(hname,'geolat_t');
area = ncread(gname,'area_t');[xL,yL] = size(lon);
lonu = ncread(hname,'geolon_c');latu = ncread(hname,'geolat_c');

z = ncread(hname,'st_ocean');zL = length(z);

time = ncread(hname,'time');
dys = [31 28 31 30 31 30 31 31 30 31 30 31];
C = textread(rnametime, '%s','delimiter', '\n');
rtime = [str2num(C{3}(1:3)) str2num(C{3}(8:9)) str2num(C{3}(14:15)) str2num(C{3}(20:21)) str2num(C{3}(26:27)) str2num(C{3}(32:33))];
time_snap = [(rtime(1)-1)*365+sum(dys(1:(rtime(2)-1)))+(rtime(3)-1)+rtime(4)/24+rtime(5)/24/60+rtime(6)/24/60/60;
             ncread(sname,'time')];
if (time_snap(end) == 0) time_snap(end) = 365;end
tL = length(time);

Cp = 3992.1; % J kg-1 degC-1
rho0 = 1035; % kgm-3

T = ncread(wname,'neutral');
Te = ncread(wname,'neutralrho_edges');
TL = length(T);dT = T(2)-T(1);

%% Get lat-depth slice of variables:
for lonsl=[-110 -140]
[tmp lt1] = min(abs(lat(1,:)+20));
[tmp lt2] = min(abs(lat(1,:)-20));

[tmp lnind] = min(abs(lon(:,round(mean([lt1 lt2])))-lonsl));

temp = squeeze(ncread(fname,'temp',[lnind lt1 1 1],[1 lt2-lt1+1 zL tL]));
if (strfind(baseD,'01'))
    u = squeeze(ncread(m3name,'u',[lnind lt1 1 1],[1 lt2-lt1+1 zL 1]));
    v = squeeze(ncread(m3name,'v',[lnind lt1 1 1],[1 lt2-lt1+1 zL 1]));
else
    u = squeeze(ncread(fname,'u',[lnind lt1 1 1],[1 lt2-lt1+1 zL tL]));
    v = squeeze(ncread(fname,'v',[lnind lt1 1 1],[1 lt2-lt1+1 zL tL]));
end

kappa = squeeze(ncread(fname,'diff_cbt_t',[lnind lt1 1 1],[1 lt2-lt1+1 zL tL]));
taux = squeeze(ncread(fname,'tau_x',[lnind lt1 1],[1 lt2-lt1+1 tL]));
tauy = squeeze(ncread(fname,'tau_y',[lnind lt1 1],[1 lt2-lt1+1 tL]));
mld = squeeze(ncread(fname,'mld',[lnind lt1 1],[1 lt2-lt1+1 tL]));
vdif = squeeze(ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[lnind lt1 1 1],[1 lt2-lt1+1 TL tL]));
vnlc = squeeze(ncread(wname,'temp_nonlocal_KPP_on_nrho',[lnind lt1 1 1],[1 lt2-lt1+1 TL tL]));
pmer = squeeze(ncread(wname,'sfc_hflux_pme_on_nrho',[lnind lt1 1 1],[1 lt2-lt1+1 TL tL]) + ...
               ncread(wname,'temp_rivermix_on_nrho',[lnind lt1 1 1],[1 lt2-lt1+1 TL tL]));
sufc = squeeze(ncread(wname,'temp_vdiffuse_sbc_on_nrho',[lnind lt1 1 1],[1 lt2-lt1+1 TL tL]) + ...
               ncread(wname,'frazil_on_nrho',[lnind lt1 1 1],[1 lt2-lt1+1 TL tL]) + ...
               ncread(wname,'temp_eta_smooth_on_nrho',[lnind lt1 1 1],[1 lt2-lt1+1 TL tL]) + ...
               ncread(wname,'sw_heat_on_nrho',[lnind lt1 1 1],[1 lt2-lt1+1 TL tL]));
swrd = squeeze(ncread(wname,'sw_heat_on_nrho',[lnind lt1 1 1],[1 lt2-lt1+1 TL tL]));

[Yt,Zt] = ndgrid(lat(lnind,lt1:lt2),z);
[Yu,Zu] = ndgrid(latu(lnind,lt1:lt2),z);

name = [outD 'mat_data/' model sprintf('_output%03d',output) '_varsat_' num2str(-lonsl) 'W.mat']
save(name,'Yt','Zt','Yu','Zu','temp','u','v','kappa','taux','tauy','mld', ...
     'vdif','vnlc','pmer','sufc','swrd');
end

%% Get equatorial slices of variables:
latsl = 0;
[tmp ltind] = min(abs(lat(1,:)-latsl));
[tmp ln1] = min(abs(lon(:,ltind)+240));
[tmp ln2] = min(abs(lon(:,ltind)+70));


temp = squeeze(ncread(fname,'temp',[ln1 ltind 1 1],[ln2-ln1+1 1 zL tL]));
if (strfind(baseD,'01'))
    u = squeeze(ncread(m3name,'u',[ln1 ltind 1 1],[ln2-ln1+1 1 zL 1]));
    v = squeeze(ncread(m3name,'v',[ln1 ltind 1 1],[ln2-ln1+1 1 zL 1]));
else
    u = squeeze(ncread(fname,'u',[ln1 ltind 1 1],[ln2-ln1+1 1 zL tL]));
    v = squeeze(ncread(fname,'v',[ln1 ltind 1 1],[ln2-ln1+1 1 zL tL]));
end
kappa = squeeze(ncread(fname,'diff_cbt_t',[ln1 ltind 1 1],[ln2-ln1+1 1 zL tL]));
taux = squeeze(ncread(fname,'tau_x',[ln1 ltind 1],[ln2-ln1+1 1 tL]));
tauy = squeeze(ncread(fname,'tau_y',[ln1 ltind 1],[ln2-ln1+1 1 tL]));
mld = squeeze(ncread(fname,'mld',[ln1 ltind 1],[ln2-ln1+1 1 tL]));
vdif = squeeze(ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]));
vnlc = squeeze(ncread(wname,'temp_nonlocal_KPP_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]));
pmer = squeeze(ncread(wname,'sfc_hflux_pme_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]) + ...
               ncread(wname,'temp_rivermix_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]));
sufc = squeeze(ncread(wname,'temp_vdiffuse_sbc_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]) + ...
               ncread(wname,'frazil_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]) + ...
               ncread(wname,'temp_eta_smooth_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]) + ...
               ncread(wname,'sw_heat_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]));
swrd = squeeze(ncread(wname,'sw_heat_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]));

[Xt,Zt] = ndgrid(lon(ln1:ln2,ltind),z);
[Xu,Zu] = ndgrid(lonu(ln1:ln2,ltind),z);

name = [outD 'mat_data/' model sprintf('_output%03d',output) '_varsat_Eq.mat']
save(name,'Xt','Zt','Xu','Zu','temp','u','v','kappa','taux','tauy','mld', ...
     'vdif','vnlc','pmer','sufc','swrd');

end

