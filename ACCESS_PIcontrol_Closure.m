% This script checks the closure of the Eulerian heat budget

base = '/g/data/p66/cm2704/archive/bi889/history/ocn/';

fname = [base 'ocean_month.nc-08500630'];

area = ncread(fname,'area_t');
[xL,yL] = size(area);
zL = 50;
rho0 = 1035;
Cp = 3992.10322329649; % J kg-1 degC-1

% $$$ %%% HEAT BUDGET CLOSURE
% $$$ 
% $$$ hfds_prior = [];
% $$$ hfds_new = [];
% $$$ temp_tend = [];
% $$$ OHC = [];
% $$$ time = [];
% $$$ DT_A = [];
% $$$ 
% $$$ vars2D = {'sfc_hflux_pme','sfc_hflux_coupler', ...
% $$$           'sfc_hflux_from_runoff', 'sfc_hflux_from_water_evap', ...
% $$$           'sfc_hflux_from_water_prec', 'frazil_2d'};
% $$$ vars3D = {'frazil_3d','temp_tendency'};
% $$$ 
% $$$ files = dir(base);
% $$$ 
% $$$ for fi = 1:length(files)
% $$$     if (strfind(files(fi).name,'month'))
% $$$ 
% $$$         fname = [base files(fi).name];
% $$$         sprintf('Doing %03d of %03d',fi,length(files))
% $$$         time_t = ncread(fname,'time');
% $$$         DT_A_t = ncread(fname,'average_DT')*86400;
% $$$         
% $$$         time = cat(1,time,time_t);
% $$$         DT_A = cat(1,DT_A,DT_A_t);
% $$$ 
% $$$         tL = length(time_t);
% $$$         
% $$$         for vi = 1:length(vars2D)
% $$$             tmp = nansum(nansum(ncread(fname,vars2D{vi},[1 1 1],[xL yL tL]).*area,1),2);
% $$$             eval([vars2D{vi} ' = squeeze(tmp);']);
% $$$         end
% $$$ 
% $$$         for vi = 1:length(vars3D)
% $$$             tmp = nansum(nansum(nansum(ncread(fname,vars3D{vi},[1 1 1 1],[xL yL zL tL]),3).*area,1),2);
% $$$             eval([vars3D{vi} ' = squeeze(tmp);']);
% $$$         end
% $$$ 
% $$$         OHC_t = squeeze(Cp*nansum(nansum(nansum(ncread(fname,'temp_rhodzt',[1 1 1 1],[xL yL zL ...
% $$$                             tL]),3).*area,1),2));
% $$$ 
% $$$         hfds_prior_t = sfc_hflux_from_runoff + ...
% $$$             sfc_hflux_coupler + ...
% $$$             sfc_hflux_from_water_evap + ...
% $$$             sfc_hflux_from_water_prec + ...
% $$$             frazil_2d;
% $$$         hfds_prior_t = hfds_prior_t.*DT_A_t;
% $$$         
% $$$         hfds_new_t  = sfc_hflux_from_runoff + ...
% $$$             sfc_hflux_coupler + ...
% $$$             sfc_hflux_pme + ...
% $$$             frazil_3d;
% $$$         hfds_new_t = hfds_new_t.*DT_A_t;
% $$$ 
% $$$         temp_tend_t = temp_tendency.*DT_A_t;
% $$$ 
% $$$         OHC = cat(1,OHC,OHC_t);
% $$$         hfds_prior = cat(1,hfds_prior,hfds_prior_t);
% $$$         hfds_new = cat(1,hfds_new,hfds_new_t);
% $$$         temp_tend = cat(1,temp_tend,temp_tend_t);
% $$$ 
% $$$         if (mod(fi,20)==0)
% $$$             save('ACCESS-CM2_PIcontrol_Closure.mat','OHC','time','DT_A','temp_tend','hfds_new','hfds_prior');
% $$$         end
% $$$     end
% $$$ end

%%% MASS BUDGET CLOSURE

OM = [];
OMeta = [];
time = [];
DT_A = [];

vars2D = {'pme_river','pme_mass','pme_net',...
          'wfimelt','wfiform'};
for vi=1:length(vars2D)
    eval([vars2D{vi} ' = [];']);
end

files = dir(base);

for fi = 1:length(files)
    if (strfind(files(fi).name,'month'))

        fname = [base files(fi).name];
        sprintf('Doing %03d of %03d',fi,length(files))
        time_t = ncread(fname,'time');
        DT_A_t = ncread(fname,'average_DT')*86400;
        
        time = cat(1,time,time_t);
        DT_A = cat(1,DT_A,DT_A_t);

        tL = length(time_t);
        
        for vi = 1:length(vars2D)
            tmp = nansum(nansum(ncread(fname,vars2D{vi},[1 1 1],[xL yL tL]).*area,1),2);
            eval([vars2D{vi} ' = cat(1,' vars2D{vi} ',squeeze(tmp).*DT_A_t);']);
        end

        OM_t = squeeze(nansum(nansum(nansum(ncread(fname,'rho_dzt',[1 1 1 1],[xL yL zL ...
                            tL]),3).*area,1),2));
        OMeta_t = rho0*squeeze(nansum(nansum(ncread(fname,'eta_t',[1 ...
                            1 1],[xL yL tL]).*area,1),2));
        OM = cat(1,OM,OM_t);
        OMeta = cat(1,OMeta,OMeta_t);

        if (mod(fi,20)==0)
            save('ACCESS-CM2_PIcontrol_Mass_Closure.mat','OM','time','DT_A','OMeta');
            for vi=1:length(vars2D)
                save('ACCESS-CM2_PIcontrol_Mass_Closure.mat',vars2D{vi},'-append');
            end
        end
    end
end


        


% $$$ for vi = 1:length(vars2D)
% $$$     eval(['var = ' vars2D{vi} ';']);
% $$$     fprintf(['Total Flux = %6.6e, ' vars2D{vi} '\n'],var);
% $$$ end
% $$$ 
% $$$ for vi = 1:length(vars3D)
% $$$     eval(['var = ' vars3D{vi} ';']);
% $$$     fprintf(['Total Flux = %6.6e, ' vars3D{vi} '\n'],var);
% $$$ end
% $$$ 
% $$$ total_sfc_from2D = frazil_3d+sfc_hflux_coupler+sfc_hflux_pme+sfc_hflux_from_runoff;
% $$$ fprintf(['Total Flux = %6.6e, total_sfc_from2D = frazil_3d+sfc_hflux_coupler+sfc_hflux_pme+sfc_hflux_from_runoff \n'],total_sfc_from2D);
% $$$ 
% $$$ % $$$ total_sfc_from3D = temp_vdiffuse_sbc+temp_rivermix+frazil_3d+sfc_hflux_pme;
% $$$ % $$$ fprintf(['Total Flux = %6.6e, total_sfc_from3D = temp_vdiffuse_sbc+temp_rivermix+frazil_3d+sfc_hflux_pme \n'],total_sfc_from3D);
% $$$ 
% $$$ hfds_prior = sfc_hflux_from_runoff + ...
% $$$     sfc_hflux_coupler + ...
% $$$     sfc_hflux_from_water_evap + ...
% $$$     sfc_hflux_from_water_prec + ...
% $$$     frazil_2d;
% $$$ fprintf(['Total Flux = %6.6e, hfds_prior = sfc_hflux_from_runoff+sfc_hflux_coupler+sfc_hflux_from_water_evap+sfc_hflux_from_water_prec+frazil_2d \n'],hfds_prior);
% $$$ 
% $$$ hfds_prior_with_frazil3d = hfds_prior-frazil_2d+frazil_3d;
% $$$ fprintf(['Total Flux = %6.6e, hfds_prior_with_frazil3d = hfds_prior ' ...
% $$$          '-frazil_2d+frazil_3d \n'],hfds_prior_with_frazil3d);
% $$$ 
% $$$ hfds_new  = sfc_hflux_from_runoff + ...
% $$$     sfc_hflux_coupler + ...
% $$$     sfc_hflux_pme + ...
% $$$     frazil_3d;
% $$$ fprintf(['Total Flux = %6.6e, hfds_new = sfc_hflux_from_runoff+sfc_hflux_coupler+sfc_hflux_pm+frazil_3d \n'],hfds_new);
% $$$ 
% $$$ 
% $$$ 
