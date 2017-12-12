% This script checks the closure of the Eulerian heat budget

close all;
clear all;

% Load Base Variables:
% $$$ base = '/srv/ccrc/data03/z3500785/MOM_wombat/mat_data/';
% $$$ model = 'MOM025';
% $$$ outputs = [1978];
% $$$ base = '/short/e14/rmh561/access-om2/control/1deg_jra55_ryf/archive/mat_data/';
% $$$ model = 'ACCESS-OM2_1deg_jra55_ryf';
% $$$ outputs = [4];
base = '/short/e14/rmh561/access-om2/control/025deg_jra55_ryf8485/archive/';
model = 'ACCESS-OM2_025deg_jra55_ryf8485';
outputs = [078];

load([base 'mat_data/' model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
ndays = diff(time_snap);
region = 'Global';

haveRedi = 0; %Redi on
haveGM = 0; %GM on
haveMD = 0; %mixdownslope on
havef3D = 1; %frazil_3d

%Eulerian budget:
% $$$ hname = '/srv/ccrc/data03/z3500785/MOM_wombat/output1978/ocean_heat.nc';
hname = sprintf([base 'output%03d/ocean/ocean_heat.nc'],outputs);
zi = 5;
ti = 1;
residual = ...
    ncread(hname,'temp_tendency',[1 1 zi ti],[xL yL 1 1]) - ...
    ncread(hname,'temp_advection',[1 1 zi ti],[xL yL 1 1]) - ...
    ncread(hname,'temp_submeso',[1 1 zi ti],[xL yL 1 1]) - ...
    ncread(hname,'temp_vdiffuse_diff_cbt',[1 1 zi ti],[xL yL 1 1]) - ...
    ncread(hname,'temp_nonlocal_KPP',[1 1 zi ti],[xL yL 1 1]) - ...
    ncread(hname,'sw_heat',[1 1 zi ti],[xL yL 1 1]) - ...
    ncread(hname,'temp_vdiffuse_sbc',[1 1 zi ti],[xL yL 1 1]) - ...
    ncread(hname,'temp_rivermix',[1 1 zi ti],[xL yL 1 1]);
    if (haveRedi)
        residual = residual - ...
            ncread(hname,'temp_vdiffuse_k33',[1 1 zi ti],[xL yL 1 1])- ...
            ncread(hname,'neutral_diffusion_temp',[1 1 zi ti],[xL yL 1 1]);
    end
    if (haveGM)
        residual = residual - ...
            ncread(hname,'neutral_gm_temp',[1 1 zi ti],[xL yL 1 1]);
    end
    if (haveMD)
        residual = residual - ...
            ncread(hname,'mixdownslope_temp',[1 1 zi ti],[xL yL 1 1]);
    end
    if (zi == 1)
        residual = residual - ...
            ncread(hname,'sfc_hflux_pme',[1 1 ti],[xL yL 1]) - ...
            ncread(hname,'temp_eta_smooth',[1 1 ti],[xL yL 1]);
    end
    if (havef3D)
        residual = residual - ...
            ncread(hname,'frazil_3d',[1 1 zi ti],[xL yL 1 1]);
    else
        if (zi == 1)
            residual = residual -  ...
            ncread(hname,'frazil_2d',[1 1 ti],[xL yL 1]);
        end
    end
            
figure;
pcolPlot(lon,lat,residual);
pcolPlot(lon,lat,ncread(hname,'frazil_3d',[1 1 zi ti],[xL yL 1 1]));

%T-binned budget:
% $$$ wname = '/srv/ccrc/data03/z3500785/MOM_wombat/output1978/ocean_wmass.nc';
wname = sprintf([base 'output%03d/ocean/ocean_wmass.nc'],outputs);

%all: 
[tmp ii] = min(abs(Te - 22));
% $$$ ii = ;
TLL = 1;
ti = 1;
residual = ...
    ncread(wname,'temp_tendency_on_nrho',[1 1 ii ti],[xL yL TLL 1]) - ...
    ncread(wname,'temp_advection_on_nrho',[1 1 ii ti],[xL yL TLL 1]) - ...
    ncread(wname,'temp_submeso_on_nrho',[1 1 ii ti],[xL yL TLL 1]) - ...
    ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL TLL 1]) - ...
    ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL TLL 1]) - ...
    ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL TLL 1]) - ...
    ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL TLL 1]) - ...
    ncread(wname,'temp_rivermix_on_nrho',[1 1 ii ti],[xL yL TLL 1]) - ...
    ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 ii ti],[xL yL TLL 1]) - ...
    ncread(wname,'frazil_on_nrho',[1 1 ii ti],[xL yL TLL 1]) - ...
    ncread(wname,'temp_eta_smooth_on_nrho',[1 1 ii ti],[xL yL TLL 1]);
    if (haveRedi)
        residual = residual - ...
            ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 ii ti],[xL yL TLL 1])- ...
            ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 ii ti],[xL yL TLL 1]);
    end
    if (haveGM)
        residual = residual - ...
            ncread(wname,'neutral_gm_on_nrho_temp',[1 1 ii ti],[xL yL TLL 1]);
    end
    if (haveMD)
        residual = residual - ...
            ncread(wname,'mixdownslope_temp_on_nrho',[1 1 ii ti],[xL yL TLL 1]);
    end
            
figure;
pcolPlot(lon,lat,residual);
pcolPlot(lon,lat,ncread(wname,'mixdownslope_temp_on_nrho',[1 1 ii ti],[xL yL 1 1]));
pcolPlot(lon,lat,ncread(wname,'frazil_on_nrho',[1 1 ii ti],[xL yL 1 1]));

