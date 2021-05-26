close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/';

hname = [base 'ncdata/ocean_heat.aom2025.out086.ncra.nc'];
gname = [base 'ncdata/ocean_grid.aom2025.nc'];

Cp = 3992.10322329649; % J kg-1 degC-1
rho0 = 1035; % kgm-3

lon = ncread(gname,'geolon_t');
lat = ncread(gname,'geolat_t');
area = ncread(gname,'area_t');
ht = ncread(gname,'ht');

sq_adv = ncread(hname,'temp_sq_advection')/Cp; % W degK m^-2
adv_sq = ncread(hname,'temp_advection_sq')/Cp/rho0; % W degK m^-1 s^-1

Dt = 365*2*86400;

tot = nansum(nansum(nansum(adv_sq,3).*area./ht,2),1)*Dt/1e15



V = 1.3e18;
tot = nansum(nansum(nansum((adv.*repmat(area,[1 1 50])/rho0/Cp).^2,1),2),3)*rho0*Cp*Dt/V;






