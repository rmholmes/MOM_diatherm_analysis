% Get and plot the vertical grids from the 1-degree simulations.

base = '/scratch/e14/rmh561/access-om2/archive/';

names = {{'GFDL50','1deg_jra55_ryf_gfdl50'}, ...
         {'KDS50','1deg_jra55_ryf_kds50'}, ...
         {'KDS75','1deg_jra55_ryf_kds75'}, ...
         {'KDS100','1deg_jra55_ryf_kds100'}, ...
         {'KDS150','1deg_jra55_ryf_kds135'}};

for i = 1:length(names)
    out{i}{1} = names{i}{1};
    gname = [base names{i}{2} '/output031/ocean/ocean.nc'];
    z = ncread(gname,'st_ocean');
    dz = diff(ncread(gname,'st_edges_ocean'));
    out{i}{2} = z;%names{i}{1};
    out{i}{3} = dz;%names{i}{1};
end
save('vertical_grids.mat','out');



