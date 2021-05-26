% Process surface forcing and parameterized mixing terms from
% ACCESS-CM2 runs binned into temperature coordinates.

% $$$ # Surface heat fluxes (not including surface volume flux terms):
% $$$ SFCH = temp_vdiffuse_sbc + sw_heat + frazil_3d + # 3D vars
% $$$       temp_eta_smooth; # 2D vars
% $$$       
% $$$ # Surface heat fluxes from surface volume fluxes
% $$$ SFCV = temp_rivermix + # 3D vars
% $$$        sfc_hflux_pme # 2D vars
% $$$        
% $$$ # Shortwave redistribution
% $$$ SWR = sw_heat # 3D vars
% $$$ 
% $$$ # Vertical mixing
% $$$ VMIX = temp_vdiffuse_diff_cbt + temp_nonlocal_KPP # 3D vars
% $$$ 
% $$$ # Miscellaneous mixing
% $$$ SMIX = mixdownslope_temp + temp_sigma_diff # 3D vars
% $$$ 
% $$$ # Neutral diffusion
% $$$ RMIX = temp_vdiffuse_k33 + neutral_diffusion_temp # 3D vars
% $$$ 
% $$$ # Total tendency
% $$$ TEN = temp_tendency # 3D vars
% $$$ 
% $$$ # Total external surface forcing
% $$$ SFC = SFCH + SFCV
% $$$ 
% $$$ # Total internal surface forcing
% $$$ SFCI = SFC - rho0*Cp*THETA*SVF
% $$$ 
% $$$ # where SVF is the accumulated total surface volume flux above the
% $$$ # temperature THETA in m3s-1, rho0 = 1035.0 kgm-3 and 
% $$$ # Cp = 3992.10322329649 J kg-1 degC-1
% $$$ 
% $$$ # Total explicit mixing
% $$$ MIX = VMIX+SMIX+RMIX
% $$$ 
% $$$ # Numerical mixing (by residual)
% $$$ NMIX = dHI/dt - SFCI - MIX
% $$$ 
% $$$ # where dHI/dt is the internal heat content tendency (which will
% $$$ # have to be calculated from the monthly average heat content as
% $$$ # you've already done, although you could also try using TEN above).

clear all;

plot_only = 0;
PI_or_his = 1; % 1 = PI-control, 0 = historical simualtion
mname = 'ACCESS_PIcontrol_TVP.mat';
% $$$ PI_or_his = 0; % 1 = PI-control, 0 = historical simualtion
% $$$ mname = 'ACCESS_SpecificHeat_historical.mat';

if (~plot_only)

    if (PI_or_his)
        base = '/g/data/p66/cm2704/archive/bi889/history/ocn/';
        name = 'PIcontrol';
        fname = [base 'ocean_month.nc-08500630'];
    else
        base = '/g/data/p66/cm2704/archive/bj594/history/ocn/';
        name = 'historical';
        fname = [base 'ocean_month.nc-18500630'];
    end

    outD = '/scratch/e14/rmh561/access-cm2/';

    % Constants:
    Cp = 3992.10322329649; % J kg-1 degC-1
    rho0 = 1035; % kgm-3

    % Define a temperature grid:
    dT = 0.5;
    Te = -3:dT:34;
    T = (Te(2:end)+Te(1:end-1))/2;
    TL = length(T);
    
    % Define a percentile grid:
    dP = 0.5;
    Pe = 0:dP:100;
    P = (Pe(2:end)+Pe(1:end-1))/2;
    PL = length(P);

    % Constant grid parameters:
    area = ncread(fname,'area_t');
% $$$     lon = ncread(fname,'geolon_t');
% $$$     lat = ncread(fname,'geolat_t');
    [xL,yL] = size(area);
    zL = 50;

    % Depth grid:
    Z = ncread(fname,'st_ocean');
    Ze = ncread(fname,'st_ocean_edges');
    zL = length(Z);
    
    % 3D mask:
    temp = ncread(fname,'temp',[1 1 1 1],[xL yL zL 1]);
    mask = ~isnan(temp);
    
    % A(z):
    A = zeros(zL,1);
    for zi=1:zL
        A(zi) = nansum(nansum(area(mask(:,:,zi))));
    end
    
    % Save grid info:
    save(mname,'Cp','rho0','dT','Te','T','TL','dP','Pe','P','PL',...
               'xL','yL','zL','Z','Ze','zL','A');
    
    % Initialize variables:
    Hz = []; % Heat content at depth z (J)
    Vz = []; % Volume at depth z (m3)
    HT = []; % Heat content at temperature T (J)
    VT = []; % Volume at temperature T (m3)
    time = []; % time axis
    DT_A = []; % averaging time
    
    % Start file loop:
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
            
            temp = ncread(fname,'temp');
            temp(~mask) = NaN;
            if (max(max(temp))>120); temp=temp-273.15;end;
            V = ncread(fname,'dzt').*repmat(area,[1 1 zL tL]);
            V(isnan(V)) = 0;
            
            % Depth space calculations:
            Hz_t = squeeze(nansum(nansum(rho0*Cp*temp.*V,1),2));
            Vz_t = squeeze(nansum(nansum(V,1),2));
            
            Hz = cat(2,Hz,Hz_t);
            Vz = cat(2,Vz,Vz_t);
            
            % Temp space calculations:
            VT_t = zeros(TL,tL);
            HT_t = zeros(TL,tL);            
            for ti=1:tL
                for Ti=1:TL
                    %Accumulate sums:
                    inds = temp(:,:,:,ti)>=Te(Ti) & temp(:,:,:,ti)<Te(Ti+1);
                    Vt_t(Ti,ti) = Vt_t(Ti,ti) + nansum(nansum(nansum(V(:,:,:,ti).*inds,1),2),3);
                    Ht_t(Ti,ti) = Ht_t(Ti,ti) + nansum(nansum(nansum(rho0*Cp*temp(:,:,:,ti).*V(:,:,:,ti).*inds,1),2),3);
                end
            end
            HT = cat(2,HT,HT_t);
            VT = cat(2,VT,VT_t);
            
            if (mod(fi,5)==0)
                save(mname,'time','DT_A','Hz','Vz','HT','VT','-append');%, ...
            end
        end
    end
end

% $$$ % term options:
% $$$ haveRedi = 1; % 1 = Redi diffusion is on, 0 = off
% $$$ haveMDS = 1; % 1 = MDS is on, 0 = off;
% $$$ haveSIG = 1; % 1 = SIG is on, 0 = off;
% $$$ 
% $$$ % files:
% $$$ filenames = dir(baseL);
% $$$ 
% $$$ for fi=1:length(filenames)
% $$$     fname = [baseL filenames(fi).name];
% $$$     if (length(strfind(fname,'ocean_month.nc'))>0)
% $$$         dstr = fname(end-8:end);
% $$$ 
% $$$         % grids:
% $$$         area = ncread(fname,'area_t');
% $$$         time = ncread(fname,'time');
% $$$         [xL,yL] = size(area);
% $$$         tL = length(time);
% $$$ 
% $$$ GWBmon.SWH    = zeros(TL+1,tL);GWBmon.VDS    = zeros(TL+1,tL);
% $$$ GWBmon.RMX    = zeros(TL+1,tL);GWBmon.PME    = zeros(TL+1,tL);
% $$$ GWBmon.FRZ    = zeros(TL+1,tL);GWBmon.ETS    = zeros(TL+1,tL);
% $$$ GWBmon.VDF    = zeros(TL+1,tL);GWBmon.KNL    = zeros(TL+1,tL);
% $$$ GWBmon.TEN    = zeros(TL+1,tL);
% $$$ if (haveRedi)
% $$$ GWBmon.K33    = zeros(TL+1,tL);
% $$$ GWBmon.RED    = zeros(TL+1,tL);
% $$$ end
% $$$ if (haveMDS)
% $$$ GWBmon.MDS   = zeros(TL+1,tL);
% $$$ end
% $$$ if (haveSIG)
% $$$ GWBmon.SIG   = zeros(TL+1,tL);
% $$$ end
% $$$ 
% $$$ for zi=1:zL
% $$$     for ti=1:tL
% $$$         sprintf('Calculating MON/AN binned time %03d of %03d, depth %02d of %02d',ti,tL,zi,zL)
% $$$ 
% $$$         temp = ncread(fname,'temp',[1 1 zi ti],[xL yL 1 1]);
% $$$         if (max(max(temp))>120);temp = temp-273.15;end;
% $$$         temp(abs(temp)>100) = NaN;
% $$$         
% $$$         TEN = area.*ncread(fname,'temp_tendency',[1 1 zi ti],[xL yL 1 1]);
% $$$         RMX = area.*ncread(fname,'temp_rivermix',[1 1 zi ti],[xL yL 1 1]);
% $$$         VDS = area.*ncread(fname,'temp_vdiffuse_sbc',[1 1 zi ti],[xL yL 1 1]);
% $$$         SWH = area.*ncread(fname,'sw_heat',[1 1 zi ti],[xL yL 1 1]);
% $$$         VDF = area.*ncread(fname,'temp_vdiffuse_diff_cbt',[1 1 zi ti],[xL yL 1 1]);
% $$$         KNL = area.*ncread(fname,'temp_nonlocal_KPP',[1 1 zi ti],[xL yL 1 1]);
% $$$         FRZ = area.*ncread(fname,'frazil_3d',[1 1 zi ti],[xL yL 1 1]);
% $$$         if (haveRedi)
% $$$             K33 = area.*ncread(fname,'temp_vdiffuse_k33',[1 1 zi ti],[xL yL 1 1]);
% $$$             RED = area.*ncread(fname,'neutral_diffusion_temp',[1 1 zi ti],[xL yL 1 1]);
% $$$         end
% $$$         if (haveMDS)
% $$$             MDS = area.*ncread(fname,'mixdownslope_temp',[1 1 zi ti],[xL yL 1 1]);
% $$$         end
% $$$         if (haveSIG)
% $$$             SIG = area.*ncread(fname,'temp_sigma_diff',[1 1 zi ti],[xL yL 1 1]);
% $$$         end
% $$$ 
% $$$         if (zi == 1)
% $$$             ETS = area.*ncread(fname,'temp_eta_smooth',[1 1 ti],[xL yL 1]);
% $$$             PME = area.*ncread(fname,'sfc_hflux_pme',[1 1 ti],[xL yL 1]);
% $$$         end
% $$$                 
% $$$         for Ti=1:TL
% $$$             %Accumulate sums:
% $$$             inds = find(temp>=Te(Ti) & temp<Te(Ti+1));
% $$$             GWBmon.TEN(Ti,ti) = GWBmon.TEN(Ti,ti)+nansum(TEN(inds));
% $$$             GWBmon.RMX(Ti,ti) = GWBmon.RMX(Ti,ti)+nansum(RMX(inds));
% $$$             GWBmon.VDS(Ti,ti) = GWBmon.VDS(Ti,ti)+nansum(VDS(inds));
% $$$             GWBmon.SWH(Ti,ti) = GWBmon.SWH(Ti,ti)+nansum(SWH(inds));
% $$$             GWBmon.VDF(Ti,ti) = GWBmon.VDF(Ti,ti)+nansum(VDF(inds));
% $$$             GWBmon.KNL(Ti,ti) = GWBmon.KNL(Ti,ti)+nansum(KNL(inds));
% $$$             GWBmon.FRZ(Ti,ti) = GWBmon.FRZ(Ti,ti)+nansum(FRZ(inds));
% $$$             if (haveRedi)
% $$$                 GWBmon.K33(Ti,ti) = GWBmon.K33(Ti,ti)+nansum(K33(inds));
% $$$                 GWBmon.RED(Ti,ti) = GWBmon.RED(Ti,ti)+nansum(RED(inds));
% $$$             end
% $$$             if (haveMDS)
% $$$                 GWBmon.MDS(Ti,ti) = GWBmon.MDS(Ti,ti)+nansum(MDS(inds));
% $$$             end
% $$$             if (haveSIG)
% $$$                 GWBmon.SIG(Ti,ti) = GWBmon.SIG(Ti,ti)+nansum(SIG(inds));
% $$$             end
% $$$             
% $$$             if (zi == 1)
% $$$                 GWBmon.ETS(Ti,ti) = GWBmon.ETS(Ti,ti)+nansum(ETS(inds));
% $$$                 GWBmon.PME(Ti,ti) = GWBmon.PME(Ti,ti)+nansum(PME(inds));
% $$$             end
% $$$         end
% $$$         inds = find(temp>=Te(TL+1));
% $$$         GWBmon.TEN(TL+1,ti) = GWBmon.TEN(TL+1,ti)+nansum(TEN(inds));
% $$$         GWBmon.RMX(TL+1,ti) = GWBmon.RMX(TL+1,ti)+nansum(RMX(inds));
% $$$         GWBmon.VDS(TL+1,ti) = GWBmon.VDS(TL+1,ti)+nansum(VDS(inds));
% $$$         GWBmon.SWH(TL+1,ti) = GWBmon.SWH(TL+1,ti)+nansum(SWH(inds));
% $$$         GWBmon.VDF(TL+1,ti) = GWBmon.VDF(TL+1,ti)+nansum(VDF(inds));
% $$$         GWBmon.KNL(TL+1,ti) = GWBmon.KNL(TL+1,ti)+nansum(KNL(inds));
% $$$         GWBmon.FRZ(TL+1,ti) = GWBmon.FRZ(TL+1,ti)+nansum(FRZ(inds));
% $$$         if (haveRedi)
% $$$             GWBmon.K33(TL+1,ti) = GWBmon.K33(TL+1,ti)+nansum(K33(inds));
% $$$             GWBmon.RED(TL+1,ti) = GWBmon.RED(TL+1,ti)+nansum(RED(inds));
% $$$         end
% $$$         if (haveMDS)
% $$$             GWBmon.MDS(TL+1,ti) = GWBmon.MDS(TL+1,ti)+nansum(MDS(inds));
% $$$         end
% $$$         if (haveSIG)
% $$$             GWBmon.SIG(TL+1,ti) = GWBmon.SIG(TL+1,ti)+nansum(SIG(inds));
% $$$         end
% $$$ 
% $$$         if (zi == 1)
% $$$             GWBmon.ETS(TL+1,ti) = GWBmon.ETS(TL+1,ti)+nansum(ETS(inds));
% $$$             GWBmon.PME(TL+1,ti) = GWBmon.PME(TL+1,ti)+nansum(PME(inds));
% $$$         end
% $$$     end
% $$$ end
% $$$ 
% $$$ % Integrate to get to T'>T:
% $$$ names = fieldnames(GWBmon);
% $$$ for i=1:length(names)
% $$$     for ti=1:tL
% $$$         eval(['GWBmon.' names{i} '(:,ti) = flipud(cumsum(flipud(GWBmon.' ...
% $$$               names{i} '(:,ti))));']);
% $$$     end
% $$$ end
% $$$ 
% $$$ save([outD model dstr '_GlobalHBud.mat'],'GWBmon','T','Te','-v7.3');
% $$$     
% $$$     end
% $$$ end
% $$$ 
