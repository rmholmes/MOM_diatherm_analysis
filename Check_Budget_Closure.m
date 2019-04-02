% This script checks the closure of the Eulerian heat budget

close all;
clear all;

addpath(genpath('/short/e14/rmh561/software/matlab-utilities'));

% Load Base Variables:
base = '/short/e14/rmh561/access-om2/archive/1deg_jra55_ryf8485_kds50_july/';
model = 'ACCESS-OM2_1deg_jra55_ryf8485_kds50_july';
outputs = [36];

load([base 'mat_data/' model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
if (~exist('ndays'))
    ndays = diff(time_snap);
end
region = 'Global';

outputs = [37];

haveRedi = 1; %Redi on
haveGM = 1; %GM on
haveMD = 1; %mixdownslope on
haveSG = 1; % sigma-diff
havef3D = 1; %frazil_3d

%%%%% Eulerian budget:
% $$$ hname = '/srv/ccrc/data03/z3500785/MOM_wombat/output1978/ocean_heat.nc';
hname = sprintf([base 'output%03d/ocean/ocean_heat.nc'],outputs);
vars3D = {'temp_tendency','temp_advection','temp_submeso', ...
        'temp_vdiffuse_diff_cbt','temp_nonlocal_KPP','sw_heat','temp_vdiffuse_sbc', ...
        'temp_rivermix'};
if (haveRedi)
    vars3D = {vars3D{:},'temp_vdiffuse_k33','neutral_diffusion_temp'};
end
if (haveGM)
    vars3D = {vars3D{:},'neutral_gm_temp'};
end
if (haveMD)
    vars3D = {vars3D{:},'mixdownslope_temp'};
end
if (haveSG)
    vars3D = {vars3D{:},'temp_sigma_diff'};
end
if (havef3D)
    vars3D = {vars3D{:},'frazil_3d'};
end

for zi=1:50
ti = 1;
if (zi == 1)
    vars = {vars3D{:},'sfc_hflux_pme','temp_eta_smooth'};
else
    vars = vars3D;
end

amps = zeros(length(vars)+1,1);
res = zeros(xL,yL);
for i=1:length(vars)
    try
        var = ncread(hname,vars{i},[1 1 zi ti],[xL yL 1 1]);
    catch
        var = ncread(hname,vars{i},[1 1 ti],[xL yL 1]);
        ['caught ' vars{i}]
    end
    amps(i) = nansum(nansum(var.^2,1),2);
    if (i==1)
        res = -var;
    else
        res = res + var;
    end
end
amps(end) = nansum(nansum(res.^2,1),2);

vars{i+1} = 'residual';
fprintf(['Depth' num2str(zi) '\n ---- \n']);
for i=1:length(vars)
    fprintf(['Sq amp = %6.1e,' vars{i} '\n'],amps(i))
end
fprintf('\n');
end

%%%% T-binned budget:
wname = sprintf([base 'output%03d/ocean/ocean_wmass.nc'],outputs);
vars3D = {'temp_tendency_on_nrho','temp_advection_on_nrho','temp_submeso_on_nrho', ...
        'temp_vdiffuse_diff_cbt_on_nrho','temp_nonlocal_KPP_on_nrho','sw_heat_on_nrho','temp_vdiffuse_sbc_on_nrho', ...
        'temp_rivermix_on_nrho','temp_eta_smooth_on_nrho','sfc_hflux_pme_on_nrho'};
if (haveRedi)
    vars3D = {vars3D{:},'temp_vdiffuse_k33_on_nrho','neutral_diffusion_on_nrho_temp'};
end
if (haveGM)
    vars3D = {vars3D{:},'neutral_gm_on_nrho_temp'};
end
if (haveMD)
    vars3D = {vars3D{:},'mixdownslope_temp_on_nrho'};
end
if (haveSG)
    vars3D = {vars3D{:},'temp_sigma_diff_on_nrho'};
end
if (havef3D)
    vars3D = {vars3D{:},'frazil_on_nrho'};
end

for Ti=1:TL
ti = 1;
vars = vars3D;

amps = zeros(length(vars)+1,1);
res = zeros(xL,yL);
for i=1:length(vars)
    var = ncread(wname,vars{i},[1 1 Ti ti],[xL yL 1 1]);
    amps(i) = nansum(nansum(var.^2,1),2);
    if (i==1)
        res = -var;
    else
        res = res + var;
    end
end
amps(end) = nansum(nansum(res.^2,1),2);

vars{i+1} = 'residual';
fprintf(['Level ' num2str(Ti) '\n ---- \n']);
for i=1:length(vars)
    fprintf(['Sq amp = %6.1e,' vars{i} '\n'],amps(i))
end
fprintf('\n');
end
