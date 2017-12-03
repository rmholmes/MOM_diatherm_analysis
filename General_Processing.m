% This script extracts specified data from MOM025 and MOM01 runs 

% $$$ model = 'MOM01';
model = 'MOM025';
% $$$ baseD = '/short/e14/rmh561/mom/archive/MOM_HeatDiag/'; %Data Directory.
baseD = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/'; %Data Directory.

for output=2:6
    restart = output-1;

%% file-names and grid properties:
base = [baseD sprintf('output%03d/',output)];
baser = [baseD sprintf('restart%03d/',restart)];
hname = [base 'ocean_heat.nc'];
if (strfind(baseD,'01'))
    fname = [base 'ocean_month.nc'];
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


% $$$ %% Get slices of variables at 140W:
% $$$ lonsl = -140;
% $$$ [tmp lt1] = min(abs(lat(1,:)+20));
% $$$ [tmp lt2] = min(abs(lat(1,:)-20));
% $$$ 
% $$$ [tmp lnind] = min(abs(lon(:,round(mean([lt1 lt2])))-lonsl));
% $$$ 
% $$$ temp = squeeze(ncread(fname,'temp',[lnind lt1 1 1],[1 lt2-lt1+1 zL tL]));
% $$$ u = squeeze(ncread(fname,'u',[lnind lt1 1 1],[1 lt2-lt1+1 zL tL]));
% $$$ v = squeeze(ncread(fname,'v',[lnind lt1 1 1],[1 lt2-lt1+1 zL tL]));
% $$$ kappa = squeeze(ncread(fname,'diff_cbt_t',[lnind lt1 1 1],[1 lt2-lt1+1 zL tL]));
% $$$ taux = squeeze(ncread(fname,'tau_x',[lnind lt1 1],[1 lt2-lt1+1 tL]));
% $$$ tauy = squeeze(ncread(fname,'tau_y',[lnind lt1 1],[1 lt2-lt1+1 tL]));
% $$$ 
% $$$ [Y,Z] = ndgrid(lat(lnind,:),z);
% $$$ 
% $$$ save([model '_varsat_' num2str(-lonsl) 'W.mat'],'Y','Z','temp','u','v','kappa','taux','tauy');

%% Get equatorial slices of variables:
latsl = 0;
[tmp ln1] = min(abs(lon(:,1)+240));
[tmp ln2] = min(abs(lon(:,1)+70));

[tmp ltind] = min(abs(lat(1,:)-latsl));

output
temp = squeeze(ncread(fname,'temp',[ln1 ltind 1 1],[ln2-ln1+1 1 zL tL]));
u = squeeze(ncread(fname,'u',[ln1 ltind 1 1],[ln2-ln1+1 1 zL tL]));
v = squeeze(ncread(fname,'v',[ln1 ltind 1 1],[ln2-ln1+1 1 zL tL]));
kappa = squeeze(ncread(fname,'diff_cbt_t',[ln1 ltind 1 1],[ln2-ln1+1 1 zL tL]));
taux = squeeze(ncread(fname,'tau_x',[ln1 ltind 1],[ln2-ln1+1 1 tL]));
tauy = squeeze(ncread(fname,'tau_y',[ln1 ltind 1],[ln2-ln1+1 1 tL]));
mld = squeeze(ncread(fname,'mld',[ln1 ltind 1],[ln2-ln1+1 1 tL]));
vdif = squeeze(ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]));
vnlc = squeeze(ncread(wname,'temp_nonlocal_KPP_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]));

output
[Xt,Zt] = ndgrid(lon(ln1:ln2,1),z);
[Xu,Zu] = ndgrid(lonu(ln1:ln2,1),z);

save([baseD 'mat_data/' model sprintf('_output%03d',output) '_varsat_Eq.mat'],'X','Z','temp','u','v','kappa','taux','tauy','vdif','vnlc','mld');

end
