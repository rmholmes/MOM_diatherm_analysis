% This script checks the closure of the Eulerian heat budget

base = '/g/data/p66/cm2704/archive/bj594/history/ocn/';
fname = [base 'ocean_month.nc-19820630'];

area = ncread(fname,'area_t');
[xL,yL] = size(area);
zL = 50;

vars2D = {'sfc_hflux_pme','sfc_hflux_coupler', ...
          'sfc_hflux_from_runoff', 'sfc_hflux_from_water_evap', ...
          'sfc_hflux_from_water_prec', 'frazil_2d'};
vars3D = {'frazil_3d','temp_vdiffuse_sbc','temp_rivermix', ...
          'temp_tendency'};


for vi = 1:length(vars2D)
    tmp = nansum(nansum(ncread(fname,vars2D{vi},[1 1 1],[xL yL 1]).*area,1),2);
    eval([vars2D{vi} ' = tmp;']);
end

for vi = 1:length(vars3D)
    tmp = nansum(nansum(nansum(ncread(fname,vars3D{vi},[1 1 1 1],[xL yL zL 1]),3).*area,1),2);
    eval([vars3D{vi} ' = tmp;']);
end

for vi = 1:length(vars2D)
    eval(['var = ' vars2D{vi} ';']);
    fprintf(['Total Flux = %6.6e, ' vars2D{vi} '\n'],var);
end

for vi = 1:length(vars3D)
    eval(['var = ' vars3D{vi} ';']);
    fprintf(['Total Flux = %6.6e, ' vars3D{vi} '\n'],var);
end

total_sfc_from2D = frazil_3d+sfc_hflux_coupler+sfc_hflux_pme+sfc_hflux_from_runoff;
fprintf(['Total Flux = %6.6e, total_sfc_from2D = frazil_3d+sfc_hflux_coupler+sfc_hflux_pme+sfc_hflux_from_runoff \n'],total_sfc_from2D);

total_sfc_from3D = temp_vdiffuse_sbc+temp_rivermix+frazil_3d+sfc_hflux_pme;
fprintf(['Total Flux = %6.6e, total_sfc_from3D = temp_vdiffuse_sbc+temp_rivermix+frazil_3d+sfc_hflux_pme \n'],total_sfc_from3D);

hfds_prior = sfc_hflux_from_runoff + ...
    sfc_hflux_coupler + ...
    sfc_hflux_from_water_evap + ...
    sfc_hflux_from_water_prec + ...
    frazil_2d;
fprintf(['Total Flux = %6.6e, hfds_prior = sfc_hflux_from_runoff+sfc_hflux_coupler+sfc_hflux_from_water_evap+sfc_hflux_from_water_prec+frazil_2d \n'],hfds_prior);

hfds_prior_with_frazil3d = hfds_prior-frazil_2d+frazil_3d;
fprintf(['Total Flux = %6.6e, hfds_prior_with_frazil3d = hfds_prior ' ...
         '-frazil_2d+frazil_3d \n'],hfds_prior_with_frazil3d);

hfds_new  = sfc_hflux_from_runoff + ...
    sfc_hflux_coupler + ...
    sfc_hflux_pme + ...
    frazil_3d;
fprintf(['Total Flux = %6.6e, hfds_new = sfc_hflux_from_runoff+sfc_hflux_coupler+sfc_hflux_pm+frazil_3d \n'],hfds_new);



