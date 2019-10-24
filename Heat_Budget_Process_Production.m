% This script processes the heat budget and associated variables in
% MOM025 or MOM01 simulations and save's into .mat files

% $$$ baseL = '/short/e14/rmh561/mom/archive/';
baseL = '/g/data/e14/rmh561/access-om2/archive/';
% $$$ baseL = '/short/e14/rmh561/access-om2/archive/';
% $$$ baseL = '/srv/ccrc/data03/z3500785/';

% MOM-SIS025:
% $$$ model = 'MOM025_kb3seg';
% $$$ baseD = [baseL 'MOM_HeatDiag_kb3seg/']; %Data Directory.
% $$$ ICdir = [baseL 'MOM_HeatDiag_kb3seg/restart100/'];
% $$$ % ACCESS-OM2:
model = 'ACCESS-OM2_025deg_jra55_iaf';
baseD = [baseL '025deg_jra55_iaf_Maurice/']; %Data Directory.
ICdir = [baseL '025deg_jra55_iaf_Maurice/'];
% $$$ % MOM-SIS01:
% $$$ model = 'MOM01';
% $$$ baseD = [baseL 'MOM01_HeatDiag/']; %Data Directory.

outD = [baseD 'mat_data/'];
rstbaseD = baseD;

post = 'ocean/'; % For ACCESS-OM2 output coulpled;
% $$$ post = ''; % For MOM-SIS.

% term options:
haveRedi = 1; % 1 = Redi diffusion is on, 0 = off
haveGM = 1; % 1 = GM is on, 0 = off;
haveSUB = 1; % 1 = submeso is on, 0 = off;
haveMDS = 0; % 1 = MDS is on, 0 = off;
haveSIG = 0; % 1 = SIG is on, 0 = off;
haveMIX = 0; % 1 = Do mixing components (vdiffuse_diff_cbt_*), 0 = don't. 

% Processing options:
doBASE     = 0; % 1 = save BaseVars.mat file
dodVdtdHdt = 0; % 1 = calculate dVdt/dHdt and save into .nc file
doNUMDIF   = 0; % 1 = calculate tempdiff x,y,T,t and save into .nc file
doSGMviac  = 0; % 1 = calculate SUB/GM influence via binned
                % convergence (otherwise uses lateral flux). The
                % better option is to use the lateral flux -> This
                % calculates the net numerical mixing (i.e. it
                % includes the numerical mixing associated with the GM
                % and SUB schemes).

doGWB      = 0; % 1 = calculate global online budget terms
doXY       = 0; % 1 = calculate spatial fluxes-on-an-isotherm
doWMT      = 1; % 1 = calculate WMT volume fluxes spatial structure
doSURF     = 0; % 1 = calculate surface flux field and SST
doZA       = 0; % 1 = calculate zonal average budget

doHND      = 0; % 1 = calculate global online numdif
doTENMON   = 0; % 1 = do monthly eulerian tendency binning
doMONANN   = 0; % 1 = calculate monthly and annually binned eulerian global budget
doXYall    = 0; % 1 = do all XY calcs (most not used)
doXYtran   = 0; % 1 = calculate vertically-integrated heat
                % transports below given isotherm/s.

% scaling constant on the transports:
if (strcmp(model(1),'A')) %ACCESS-OM2, transport in kg/s
    tsc = 1;
else % MOM-SIS, transport in 1e9 kg/s
    tsc = 1e9;
end

% $$$ for output = 86;
restart = output-1;

% $$$ region = 'Global';

% file-names -----------------------------------------
base = [baseD sprintf('output%03d/',output) post];
basem1 = [baseD sprintf('output%03d/',output-1) post];
if (output==0)
    baser = ICdir;
else
    baser = [rstbaseD sprintf('restart%03d/',restart) post];
end
hname = [base 'ocean_heat.nc'];
if (strfind(baseD,'01'))
    fname = [base 'ocean_month.nc'];
else
    fname = [base 'ocean.nc'];
end
fname2 = [base 'ocean_month.nc'];
gname = [base 'ocean_grid.nc'];
sname = [base 'ocean_snap.nc'];
wname = [base 'ocean_wmass.nc'];
tname = [base 'time_stamp.out'];
         
% Horizontal Grid  -----------------------------------------
lon = ncread(gname,'geolon_t');lat = ncread(gname,'geolat_t');
lonu = ncread(gname,'geolon_c');latu = ncread(gname,'geolat_c');
area = ncread(gname,'area_t');[xL,yL] = size(lon);
xt = ncread(gname,'xt_ocean');xu = ncread(gname,'xu_ocean');
yt = ncread(gname,'yt_ocean');yu = ncread(gname,'yu_ocean');

% Vertical grid  -----------------------------------------
z = ncread(fname,'st_ocean');zL = length(z);
if (output ==0)
    % Initial dzt accounting for partial bottom cells:
    ht = ncread(gname,'ht');ze = ncread(fname,'st_edges_ocean');
    kmt = ncread(gname,'kmt');
    dztI = repmat(permute(diff(ze),[3 2 1]),[size(ht) 1]);
    for ii = 1:xL
        for jj = 1:yL
            if (kmt(ii,jj)>1)
                dztI(ii,jj,kmt(ii,jj)) = ht(ii,jj) - ze(kmt(ii,jj));
            end
        end
    end    
end

% 3D mask ------------------------------------------------
mask = ncread(fname,'temp',[1 1 1 1],[xL yL zL 1]);
mask(~isnan(mask)) = 1; mask(isnan(mask)) = 0;
mask = mask == 1;

if (~strcmp(region,'Global'))
    [tmaskREG,umaskREG,~] = Heat_Budget_Mask(region,gname,wname,outD,model);
else
    tmaskREG = ones(xL,yL);
    umaskREG = ones(xL,yL);
end

% Time  -----------------------------------------
time = ncread(fname,'time');
ndays = ncread(fname,'average_DT');
tL = length(time);

if (exist(baser))
    found_rst = 1;rstti = 1;
    rnameT = [baser 'ocean_temp_salt.res.nc'];
    rnameZ = [baser 'ocean_thickness.res.nc'];
else
    found_rst = 0;rstti = tL;
    rnameT = [basem1 'ocean_snap.nc'];
    rnameZ = [basem1 'ocean_snap.nc'];
end

Cp = 3992.10322329649; % J kg-1 degC-1
rho0 = 1035; % kgm-3

T = ncread(wname,'neutral');
Te = ncread(wname,'neutralrho_edges');
TL = length(T);dT = T(2)-T(1);

if (doBASE)
save([outD model sprintf('_output%03d',output) '_BaseVars.mat'], ...
     'T','Te','TL','dT','Cp','rho0','time','ndays','tL', ...
     'z','zL','lon','lat','area','xL','yL','yt','yu', 'xt','xu','lonu','latu','-v7.3');
end

%% Calculate dVdt, dHdt and save back into wmass file:
%% Also, calculate TENMON
if (dodVdtdHdt)
% Create variables:
ncid = netcdf.open(wname,'NC_WRITE');
try
    id = netcdf.inqVarID(ncid,'dHdt');
    not_there = 0;
catch
    not_there = 1;
end
%If variable not there, add it:
if (not_there)
    xid = netcdf.inqDimID(ncid,'grid_xt_ocean');yid = netcdf.inqDimID(ncid,'grid_yt_ocean');zid = netcdf.inqDimID(ncid,'neutral');tid = netcdf.inqDimID(ncid,'time');
    netcdf.reDef(ncid);
    dVdtID = netcdf.defVar(ncid,'dVdt','NC_FLOAT',[xid yid zid tid]);
    netcdf.putAtt(ncid,dVdtID,'long_name','Change in time of volume within temperature bin');
    netcdf.putAtt(ncid,dVdtID,'units','Sv (10^9 kg/s)');
% $$$     netcdf.putAtt(ncid,dVdtID,'_FillValue',single(-1e20));
    dHdtID = netcdf.defVar(ncid,'dHdt','NC_FLOAT',[xid yid zid tid]);
    netcdf.putAtt(ncid,dHdtID,'long_name','Change in time of heat content within temperature bin');
    netcdf.putAtt(ncid,dHdtID,'units','Watts');
% $$$     netcdf.putAtt(ncid,dHdtID,'_FillValue',single(-1e20));
    netcdf.endDef(ncid);
else
    dVdtID = netcdf.inqVarID(ncid,'dVdt');
    dHdtID = netcdf.inqVarID(ncid,'dHdt');
end

Vsnap = zeros(xL,yL,TL);
Hsnap = zeros(xL,yL,TL);

%Do IC for Vsnap and Hsnap:
for zi = 1:zL
    sprintf('Calculating Vsnap and Hsnap IC, depth %02d of %02d',zi,zL)
    %Temperature snapshot:
    tempsnap = ncread(rnameT,'temp',[1 1 zi rstti],[xL yL 1 1]);
    tempsnap(~mask(:,:,zi)) = NaN;
    if (max(max(tempsnap))>120);tempsnap = tempsnap-273.15;end;
    
    if (found_rst)
        if (output == 0) % Initial dzt:
            Volsnap = dztI(:,:,zi).*area;
        else
            Volsnap = ncread(rnameZ,'rho_dzt',[1 1 zi rstti],[xL yL 1 1]).*area/rho0;
        end
    else
        Volsnap = ncread(rnameT,'dzt',[1 1 zi rstti],[xL yL 1 1]).*area;
    end
    Volsnap(isnan(Volsnap)) = 0;
    
    for Ti=1:TL
        %Accumulate sums:
        inds = tempsnap>=Te(Ti) & tempsnap<Te(Ti+1);
        Vsnap(:,:,Ti) = Vsnap(:,:,Ti) + Volsnap.*inds;
        Hlay = Volsnap.*tempsnap.*inds*rho0*Cp;
        Hlay(isnan(Hlay)) = 0;
        Hsnap(:,:,Ti) = Hsnap(:,:,Ti) + Hlay;
    end
end

VsnapM = Vsnap;
HsnapM = Hsnap;
Vsnap = zeros(xL,yL,TL);
Hsnap = zeros(xL,yL,TL);

if (doTENMON)
    TENMON = zeros(TL+1,tL);
end

%Do other times for Vsnap and Hsnap:
for ti=1:tL
    for zi=1:zL
        sprintf('Calculating Vsnap and Hsnap later months time %03d of %03d, depth %02d of %02d',ti,tL,zi,zL)
        
        if (doTENMON)
            temp = ncread(fname,'temp',[1 1 zi ti],[xL yL 1 1]);
            temp(~mask(:,:,zi)) = NaN;
            if (max(max(temp))>120);temp = temp-273.15;end;
        end

        tempsnap = ncread(sname,'temp',[1 1 zi ti],[xL yL 1 1]);
        tempsnap(~mask(:,:,zi)) = NaN;
        if (max(max(tempsnap))>120);tempsnap = tempsnap-273.15;end;
        Volsnap = ncread(sname,'dzt',[1 1 zi ti],[xL yL 1 1]).*area;
        Volsnap(isnan(Volsnap)) = 0;

        if (doTENMON)
            TENf = area.*ncread(hname,'temp_tendency',[1 1 zi ti],[xL ...
                                yL 1 1]);
        end
        
        for Ti=1:TL
            %Accumulate sums:
            inds = tempsnap>=Te(Ti) & tempsnap<Te(Ti+1);
            Vsnap(:,:,Ti) = Vsnap(:,:,Ti) + Volsnap.*inds;
            Hlay = Volsnap.*tempsnap.*inds*rho0*Cp;
            Hlay(isnan(Hlay)) = 0;
            Hsnap(:,:,Ti) = Hsnap(:,:,Ti) + Hlay;
            if (doTENMON)
                inds = find(temp>=Te(Ti) & temp<Te(Ti+1));
                TENMON(Ti,ti) = TENMON(Ti,ti)+nansum(TENf(inds));
            end
        end
        if (doTENMON)
            inds = find(temp>=Te(TL+1));
            TENMON(TL+1,ti) = TENMON(TL+1,ti)+nansum(TENf(inds));
        end
    end
    
    if (doTENMON)
        % Integrate to get to T'>T:
        TENMON(:,ti) = flipud(cumsum(flipud(TENMON(:,ti))));
    end
    
    for Ti=1:TL
    netcdf.putVar(ncid,dVdtID,[0 0 Ti-1 ti-1],[xL yL 1 1],(Vsnap(:,:,Ti)-VsnapM(:,:,Ti)) ...
                        /ndays(ti)/86400*rho0/1e9);
    netcdf.putVar(ncid,dHdtID,[0 0 Ti-1 ti-1],[xL yL 1 1],(Hsnap(:,:,Ti)-HsnapM(:,:,Ti)) ...
                        /ndays(ti)/86400);
    end
    VsnapM = Vsnap;
    HsnapM = Hsnap;
    Vsnap = zeros(xL,yL,TL);
    Hsnap = zeros(xL,yL,TL);
end
netcdf.close(ncid);

end

% $$$ % Calculate non-snap zonal-average annual-average volumes:
% $$$ V = zeros(yL,TL,tL);
% $$$ H = zeros(yL,TL,tL);
% $$$ for ti=1:tL
% $$$     for zi=1:zL
% $$$         sprintf('Calculating Vsnap and Hsnap later months time %03d of %03d, depth %02d of %02d',ti,tL,zi,zL)
% $$$ 
% $$$         temp = ncread(fname,'temp',[1 1 zi ti],[xL yL 1 1]);
% $$$         temp(~mask(:,:,zi)) = NaN;
% $$$         if (max(max(temp))>120);temp = temp-273.15;end;
% $$$         Vol = ncread(fname,'dzt',[1 1 zi ti],[xL yL 1 1]).*area;
% $$$         Vol(isnan(Vol)) = 0;
% $$$         
% $$$         for Ti=1:TL
% $$$             %Accumulate sums:
% $$$             inds = temp>=Te(Ti) & temp<Te(Ti+1);
% $$$             V(:,Ti,ti) = V(:,Ti,ti) + nansum(Vol.*inds,1)';
% $$$             Hlay = Vol.*temp.*inds*rho0*Cp;
% $$$             Hlay(isnan(Hlay)) = 0;
% $$$             H(:,Ti,ti) = H(:,Ti,ti) + nansum(Hlay,1)';
% $$$         end
% $$$     end
% $$$ end

%% Calculate numerical mixing by residual from heat budget and save
%% back into netcdf file:
if (doNUMDIF)
% Create variables:
ncid = netcdf.open(wname,'NC_WRITE');
varname = 'temp_numdiff_heat_on_nrho';
try
    id = netcdf.inqVarID(ncid,varname);
    not_there = 0;
catch
    not_there = 1;
end
%If variable not there, add it:
if (not_there)
    xid = netcdf.inqDimID(ncid,'grid_xt_ocean');yid = netcdf.inqDimID(ncid,'grid_yt_ocean');zid = netcdf.inqDimID(ncid,'neutralrho_edges');tid = netcdf.inqDimID(ncid,'time');
    netcdf.reDef(ncid);
    ndifID = netcdf.defVar(ncid,varname,'NC_FLOAT',[xid yid zid tid]);
    netcdf.putAtt(ncid,ndifID,'long_name',['Diffusion of heat due ' ...
                        'to numerical mixing estimated from heat ' ...
                        'fluxes binned to neutral density']);
    netcdf.putAtt(ncid,ndifID,'units','Watts/m^2');
    %    netcdf.putAtt(ncid,ndifID,'_FillValue',single(-1e20));
    netcdf.endDef(ncid);
else
    ndifID = netcdf.inqVarID(ncid,varname);
end


for ti=1:tL
    
    JI = zeros(xL,yL);
    QI = zeros(xL,yL);
    JS = zeros(xL,yL);
    dVdt = zeros(xL,yL);
    dHdt = zeros(xL,yL);
    dift = zeros(xL,yL);
    netcdf.putVar(ncid,ndifID,[0 0 TL ti-1],[xL yL 1 1],zeros(xL,yL));

    for Ti=TL:-1:1
        sprintf('Calculating numdif time %03d of %03d, Temp %03d of %03d',ti,tL,Ti,TL)
        JS = JS + ncread(wname,'mass_pmepr_on_nrho',[1 1 Ti ti],[xL yL 1 1])./area/rho0;
        dVdt = dVdt + double(ncread(wname,'dVdt',[1 1 Ti ti],[xL yL 1 1]))*1e9/rho0./area;
        dHdt = dHdt + double(ncread(wname,'dHdt',[1 1 Ti ti],[xL yL 1 1]))./area;
        dift = dift + ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
               ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
               ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
               ncread(wname,'sw_heat_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
               ncread(wname,'frazil_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
               ncread(wname,'temp_eta_smooth_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
               ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
               ncread(wname,'temp_rivermix_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        if (haveRedi)
            dift = dift + ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
                   ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 Ti ti],[xL yL 1 1]);
        end
        if (haveMDS)
            dift = dift + ncread(wname,'mixdownslope_temp_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        end    
        if (haveSIG)
            dift = dift + ncread(wname,'temp_sigma_diff_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        end    
    
        txtrans = ncread(wname,'tx_trans_nrho',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
        tytrans = ncread(wname,'ty_trans_nrho',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
        % NOTE: THE tx_trans_submeso contribution SHOULD NOT be
        % added if it is implemented via skew-diffusion.
        % NOTE: THE tx_trans_gm contribution SHOULD NOT be
        % added if it is implemented via skew-diffusion.

        qxtrans = ncread(wname,'temp_xflux_adv_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        qytrans = ncread(wname,'temp_yflux_adv_on_nrho',[1 1 Ti ti],[xL yL 1 1]);

        % Submesoscale and GM: two options:
        if (haveSUB)
            if (doSGMviac)
                dift = dift + ncread(wname,'temp_submeso_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            else
                qxtrans = qxtrans + ncread(wname,'temp_xflux_submeso_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
                qytrans = qytrans + ncread(wname,'temp_yflux_submeso_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            end
        end
        if (haveGM)
            if (doSGMviac)
                dift = dift + ncread(wname,'neutral_gm_on_nrho_temp',[1 1 Ti ti],[xL yL 1 1]);
            else
                qxtrans = qxtrans + ncread(wname,'temp_xflux_gm_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
                qytrans = qytrans + ncread(wname,'temp_yflux_gm_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            end
        end

        JI(2:end,2:end) = JI(2:end,2:end) + (txtrans(1:(end-1),2:end) - txtrans(2:end,2:end) ...
                                             +tytrans(2:end,1:(end-1)) - tytrans(2:end,2:end))./area(2:end,2:end);
        JI(1,2:end) = JI(1,2:end) + (txtrans(end,2:end) - txtrans(1,2:end) ...
                                     +tytrans(1,1:(end-1)) - tytrans(1,2:end))./area(1,2:end);        
        QI(2:end,2:end) = QI(2:end,2:end) + (qxtrans(1:(end-1),2:end) - qxtrans(2:end,2:end) ...
                                             +qytrans(2:end,1:(end-1)) - qytrans(2:end,2:end))./area(2:end,2:end);
        QI(1,2:end) = QI(1,2:end) + (qxtrans(end,2:end) - qxtrans(1,2:end) ...
                                     +qytrans(1,1:(end-1)) - qytrans(1,2:end))./area(1,2:end);        

        ndif = dHdt - (dVdt - JI - JS)*rho0*Cp*Te(Ti) - dift - QI;
    
        netcdf.putVar(ncid,ndifID,[0 0 Ti-1 ti-1],[xL yL 1 1],ndif);
    end
end
netcdf.close(ncid);

end

%% Calculate Monthly and Annual binned volume integrated budget:

if (doMONANN)
GWBmon.SWH    = zeros(TL+1,tL);GWBmon.VDS    = zeros(TL+1,tL);
GWBmon.RMX    = zeros(TL+1,tL);GWBmon.PME    = zeros(TL+1,tL);
GWBmon.FRZ    = zeros(TL+1,tL);GWBmon.ETS    = zeros(TL+1,tL);
GWBmon.SUB    = zeros(TL+1,tL);GWBmon.VDF    = zeros(TL+1,tL);
GWBmon.KNL    = zeros(TL+1,tL);GWBmon.ADV    = zeros(TL+1,tL);
GWBmon.TEN    = zeros(TL+1,tL);GWBmon.SFW    = zeros(TL+1,tL);

GWBann = GWBmon;

for zi=1:zL
    tempAN = zeros(xL,yL);
    SWHAN = zeros(xL,yL);VDSAN = zeros(xL,yL);RMXAN = zeros(xL,yL);
    SUBAN = zeros(xL,yL);VDFAN = zeros(xL,yL);KNLAN = zeros(xL,yL);
    ADVAN = zeros(xL,yL);TENAN = zeros(xL,yL);SFWAN = zeros(xL,yL);
    PMEAN = zeros(xL,yL);FRZAN = zeros(xL,yL);ETSAN = zeros(xL,yL);

    for ti=1:tL
        sprintf('Calculating MON/AN binned time %03d of %03d, depth %02d of %02d',ti,tL,zi,zL)

        temp = ncread(fname,'temp',[1 1 zi ti],[xL yL 1 1]);
        temp(~mask(:,:,zi)) = NaN;
        if (max(max(temp))>120);temp = temp-273.15;end;

        tempAN = tempAN + temp;
        
        TEN = area.*ncread(hname,'temp_tendency',[1 1 zi ti],[xL yL 1 1]);
        ADV = area.*ncread(hname,'temp_advection',[1 1 zi ti],[xL yL 1 1]);
        SUB = area.*ncread(hname,'temp_submeso',[1 1 zi ti],[xL yL 1 1]);
        RMX = area.*ncread(hname,'temp_rivermix',[1 1 zi ti],[xL yL 1 1]);
        VDS = area.*ncread(hname,'temp_vdiffuse_sbc',[1 1 zi ti],[xL yL 1 1]);
        SWH = area.*ncread(hname,'sw_heat',[1 1 zi ti],[xL yL 1 1]);
        VDF = area.*ncread(hname,'temp_vdiffuse_diff_cbt',[1 1 zi ti],[xL yL 1 1]);
        KNL = area.*ncread(hname,'temp_nonlocal_KPP',[1 1 zi ti],[xL yL 1 1]);

        if (zi == 1)
            FRZ = area.*ncread(hname,'frazil_2d',[1 1 ti],[xL yL 1]);
            ETS = area.*ncread(hname,'temp_eta_smooth',[1 1 ti],[xL yL 1]);
            PME = area.*ncread(hname,'sfc_hflux_pme',[1 1 ti],[xL yL 1]);
        end
        
        TENAN = TENAN + TEN;
        ADVAN = ADVAN + ADV;
        SUBAN = SUBAN + SUB;
        RMXAN = RMXAN + RMX;
        VDSAN = VDSAN + VDS;
        SWHAN = SWHAN + SWH;
        VDFAN = VDFAN + VDF;
        KNLAN = KNLAN + KNL;
        
        if (zi == 1)
            PMEAN = PMEAN + PME;
            FRZAN = FRZAN + FRZ;
            ETSAN = ETSAN + ETS;
        end
        
        for Ti=1:TL
            %Accumulate sums:
            inds = find(temp>=Te(Ti) & temp<Te(Ti+1));
            GWBmon.TEN(Ti,ti) = GWBmon.TEN(Ti,ti)+nansum(TEN(inds));
            GWBmon.ADV(Ti,ti) = GWBmon.ADV(Ti,ti)+nansum(ADV(inds));
            GWBmon.SUB(Ti,ti) = GWBmon.SUB(Ti,ti)+nansum(SUB(inds));
            GWBmon.RMX(Ti,ti) = GWBmon.RMX(Ti,ti)+nansum(RMX(inds));
            GWBmon.VDS(Ti,ti) = GWBmon.VDS(Ti,ti)+nansum(VDS(inds));
            GWBmon.SWH(Ti,ti) = GWBmon.SWH(Ti,ti)+nansum(SWH(inds));
            GWBmon.VDF(Ti,ti) = GWBmon.VDF(Ti,ti)+nansum(VDF(inds));
            GWBmon.KNL(Ti,ti) = GWBmon.KNL(Ti,ti)+nansum(KNL(inds));
            
            if (zi == 1)
                GWBmon.FRZ(Ti,ti) = GWBmon.FRZ(Ti,ti)+nansum(FRZ(inds));
                GWBmon.ETS(Ti,ti) = GWBmon.ETS(Ti,ti)+nansum(ETS(inds));
                GWBmon.PME(Ti,ti) = GWBmon.PME(Ti,ti)+nansum(PME(inds));
            end
        end
        inds = find(temp>=Te(TL+1));
        GWBmon.TEN(TL+1,ti) = GWBmon.TEN(TL+1,ti)+nansum(TEN(inds));
        GWBmon.ADV(TL+1,ti) = GWBmon.ADV(TL+1,ti)+nansum(ADV(inds));
        GWBmon.SUB(TL+1,ti) = GWBmon.SUB(TL+1,ti)+nansum(SUB(inds));
        GWBmon.RMX(TL+1,ti) = GWBmon.RMX(TL+1,ti)+nansum(RMX(inds));
        GWBmon.VDS(TL+1,ti) = GWBmon.VDS(TL+1,ti)+nansum(VDS(inds));
        GWBmon.SWH(TL+1,ti) = GWBmon.SWH(TL+1,ti)+nansum(SWH(inds));
        GWBmon.VDF(TL+1,ti) = GWBmon.VDF(TL+1,ti)+nansum(VDF(inds));
        GWBmon.KNL(TL+1,ti) = GWBmon.KNL(TL+1,ti)+nansum(KNL(inds));

        if (zi == 1)
            GWBmon.FRZ(TL+1,ti) = GWBmon.FRZ(TL+1,ti)+nansum(FRZ(inds));
            GWBmon.ETS(TL+1,ti) = GWBmon.ETS(TL+1,ti)+nansum(ETS(inds));
            GWBmon.PME(TL+1,ti) = GWBmon.PME(TL+1,ti)+nansum(PME(inds));
        end
    end
    
    ti = 1;
    for Ti=1:TL
        %Accumulate sums:
        inds = find(tempAN/tL>=Te(Ti) & tempAN/tL<Te(Ti+1));
        GWBann.TEN(Ti,ti) = GWBann.TEN(Ti,ti)+nansum(TENAN(inds)/tL);
        GWBann.ADV(Ti,ti) = GWBann.ADV(Ti,ti)+nansum(ADVAN(inds)/tL);
        GWBann.SUB(Ti,ti) = GWBann.SUB(Ti,ti)+nansum(SUBAN(inds)/tL);
        GWBann.RMX(Ti,ti) = GWBann.RMX(Ti,ti)+nansum(RMXAN(inds)/tL);
        GWBann.VDS(Ti,ti) = GWBann.VDS(Ti,ti)+nansum(VDSAN(inds)/tL);
        GWBann.SWH(Ti,ti) = GWBann.SWH(Ti,ti)+nansum(SWHAN(inds)/tL);
        GWBann.VDF(Ti,ti) = GWBann.VDF(Ti,ti)+nansum(VDFAN(inds)/tL);
        GWBann.KNL(Ti,ti) = GWBann.KNL(Ti,ti)+nansum(KNLAN(inds)/tL);
    
        if (zi == 1)
            GWBann.FRZ(Ti,ti) = GWBann.FRZ(Ti,ti)+nansum(FRZAN(inds)/tL);
            GWBann.ETS(Ti,ti) = GWBann.ETS(Ti,ti)+nansum(ETSAN(inds)/tL);
            GWBann.PME(Ti,ti) = GWBann.PME(Ti,ti)+nansum(PMEAN(inds)/tL);
        end
    end
    inds = find(temp>=Te(TL+1));
    GWBann.TEN(TL+1,ti) = GWBann.TEN(TL+1,ti)+nansum(TENAN(inds)/tL);
    GWBann.ADV(TL+1,ti) = GWBann.ADV(TL+1,ti)+nansum(ADVAN(inds)/tL);
    GWBann.SUB(TL+1,ti) = GWBann.SUB(TL+1,ti)+nansum(SUBAN(inds)/tL);
    GWBann.RMX(TL+1,ti) = GWBann.RMX(TL+1,ti)+nansum(RMXAN(inds)/tL);
    GWBann.VDS(TL+1,ti) = GWBann.VDS(TL+1,ti)+nansum(VDSAN(inds)/tL);
    GWBann.SWH(TL+1,ti) = GWBann.SWH(TL+1,ti)+nansum(SWHAN(inds)/tL);
    GWBann.VDF(TL+1,ti) = GWBann.VDF(TL+1,ti)+nansum(VDFAN(inds)/tL);
    GWBann.KNL(TL+1,ti) = GWBann.KNL(TL+1,ti)+nansum(KNLAN(inds)/tL);

    if (zi == 1)
        GWBann.FRZ(TL+1,ti) = GWBann.FRZ(TL+1,ti)+nansum(FRZAN(inds)/tL);
        GWBann.ETS(TL+1,ti) = GWBann.ETS(TL+1,ti)+nansum(ETSAN(inds)/tL);
        GWBann.PME(TL+1,ti) = GWBann.PME(TL+1,ti)+nansum(PMEAN(inds)/tL);
    end
end

% Integrate to get to T'>T:
names = fieldnames(GWBmon);
for i=1:length(names)
    for ti=1:tL
        eval(['GWBmon.' names{i} '(:,ti) = flipud(cumsum(flipud(GWBmon.' ...
              names{i} '(:,ti))));']);
    end
    ti = 1;
    eval(['GWBann.' names{i} '(:,ti) = flipud(cumsum(flipud(GWBann.' ...
          names{i} '(:,ti))));']);
end

save([outD model sprintf('_output%03d',output) '_' region '_HBud_MonAnBin.mat'],'GWBmon','GWBann','-v7.3');

end


%% Calculate volume integrated budget from online T-binned values -----------------------------------------------------------------------------------------------------------
if (doGWB)
    
GWB.dVdt   = zeros(TL+1,tL); % Sv 
try
    GWB.TENMON = TENMON; % Tendency from monthly averages (from previous code block)
catch
    GWB.TENMON = GWB.dVdt;
end
GWB.dHdt   = zeros(TL+1,tL); % Total W
GWB.SWH    = zeros(TL+1,tL); % W due to SW redistribution
GWB.VDS    = zeros(TL+1,tL); % W due to vdiffuse_sbc.
GWB.RMX    = zeros(TL+1,tL); % W due to rivermix.
GWB.PME    = zeros(TL+1,tL); % W due to P-E.
GWB.FRZ    = zeros(TL+1,tL); % W due to frazil.
GWB.ETS    = zeros(TL+1,tL); % W due to eta_smoothing.
GWB.VDF    = zeros(TL+1,tL); % W due to vdiffusion
GWB.KNL    = zeros(TL+1,tL); % W due to KPP non-local
if (haveSUB)
    GWB.SUB    = zeros(TL+1,tL); % W due to submesoscale.
end
if (haveMIX)
    GWB.VDFkppiw = zeros(TL+1,tL); % W due to KPP Internal Wave
    GWB.VDFkppish = zeros(TL+1,tL); % W due to KPP Shear
    GWB.VDFkppicon = zeros(TL+1,tL); % W due to KPP Convection
    GWB.VDFkppbl = zeros(TL+1,tL); % W due to KPP Boundary Layer
    GWB.VDFkppdd = zeros(TL+1,tL); % W due to KPP Double-Diffusion
    GWB.VDFwave = zeros(TL+1,tL); % W due to Simmons Internal Wave
end    
if (haveRedi)
    GWB.K33    = zeros(TL+1,tL); % W due to K33
    GWB.RED    = zeros(TL+1,tL); % W due to Redi diffusion
end
if (haveGM)
    GWB.NGM    = zeros(TL+1,tL); % W due to GM
end
if (haveMDS)
    GWB.MDS    = zeros(TL+1,tL); % W due to mixdownslope
end
if (haveSIG)
    GWB.SIG    = zeros(TL+1,tL); % W due to sigma diff
end
GWB.ADV    = zeros(TL+1,tL); % W due to advection
GWB.TEN    = zeros(TL+1,tL); % W due to tendency
GWB.SFW    = zeros(TL+1,tL); % surface volume flux into ocean (m3s-1)

if (doHND)
    GWB.NUM    = zeros(TL+1,tL); % W due to numerical mixing from heat budget
    % Get NUM:
    for ti=1:tL
        for Ti = (TL+1):-1:1
            sprintf('Calculating NUMH heat budget time %03d of %03d, temp %03d of %03d',ti,tL,Ti,TL)
            GWB.NUM(Ti,ti) = nansum(nansum(area.*ncread(wname, ...
                                                        'temp_numdiff_heat_on_nrho',[1 1 Ti ti],[xL yL 1 1]),1),2);
        end
    end
end

for ti=1:tL
    ii = TL;
    sprintf('Calculating global water-mass heat budget time %03d of %03d, temp %03d of %03d',ti,tL,ii,TL)
    GWB.dVdt(ii,ti) = nansum(nansum(tmaskREG.*ncread(wname,'dVdt',[1 1 ii ti],[xL yL 1 1])*1e9/rho0,1),2);
    GWB.dHdt(ii,ti) = nansum(nansum(tmaskREG.*ncread(wname,'dHdt',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.TEN(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_tendency_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.ADV(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_advection_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.PME(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.RMX(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_rivermix_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDS(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.SWH(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDF(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.KNL(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    if (haveSUB)
    GWB.SUB(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveMIX)
    GWB.VDFkppiw(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppiw_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFkppish(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppish_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFkppicon(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppicon_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFkppbl(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppbl_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFkppdd(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppdd_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFwave(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_wave_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    end        
    if (haveRedi)
    GWB.K33(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.RED(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveGM)
    GWB.NGM(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'neutral_gm_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveMDS)
    GWB.MDS(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'mixdownslope_temp_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveSIG)
    GWB.SIG(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_sigma_diff_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    GWB.FRZ(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'frazil_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.ETS(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_eta_smooth_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.SFW(ii,ti) = nansum(nansum(tmaskREG.*ncread(wname,'mass_pmepr_on_nrho',[1 1 ii ti],[xL yL 1 1])/rho0,1),2);
for ii=TL-1:-1:1
    sprintf('Calculating global water-mass heat budget time %03d of %03d, temp %03d of %03d',ti,tL,ii,TL)
    GWB.dVdt(ii,ti) = GWB.dVdt(ii+1,ti) + nansum(nansum(tmaskREG.*ncread(wname,'dVdt',[1 1 ii ti],[xL yL 1 1])*1e9/rho0,1),2);
    GWB.dHdt(ii,ti) = GWB.dHdt(ii+1,ti) + nansum(nansum(tmaskREG.*ncread(wname,'dHdt',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.TEN(ii,ti) = GWB.TEN(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_tendency_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.ADV(ii,ti) = GWB.ADV(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_advection_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.PME(ii,ti) = GWB.PME(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.RMX(ii,ti) = GWB.RMX(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_rivermix_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDS(ii,ti) = GWB.VDS(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.SWH(ii,ti) = GWB.SWH(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDF(ii,ti) = GWB.VDF(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.KNL(ii,ti) = GWB.KNL(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    if (haveSUB)
    GWB.SUB(ii,ti) = GWB.SUB(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveMIX)
    GWB.VDFkppiw(ii,ti)   = GWB.VDFkppiw(ii+1,ti)   + nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppiw_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFkppish(ii,ti)  = GWB.VDFkppish(ii+1,ti)  + nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppish_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFkppicon(ii,ti) = GWB.VDFkppicon(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppicon_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFkppbl(ii,ti)   = GWB.VDFkppbl(ii+1,ti)   + nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppbl_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFkppdd(ii,ti)   = GWB.VDFkppdd(ii+1,ti)   + nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppdd_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFwave(ii,ti)    = GWB.VDFwave(ii+1,ti)    + nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_wave_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);        
    end
    if (haveRedi)
    GWB.K33(ii,ti) = GWB.K33(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.RED(ii,ti) = GWB.RED(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveGM)
    GWB.NGM(ii,ti) = GWB.NGM(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'neutral_gm_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveMDS)
    GWB.MDS(ii,ti) = GWB.MDS(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'mixdownslope_temp_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveSIG)
    GWB.SIG(ii,ti) = GWB.SIG(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_sigma_diff_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    GWB.FRZ(ii,ti) = GWB.FRZ(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'frazil_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.ETS(ii,ti) = GWB.ETS(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_eta_smooth_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.SFW(ii,ti) = GWB.SFW(ii+1,ti) + nansum(nansum(tmaskREG.*ncread(wname,'mass_pmepr_on_nrho',[1 1 ii ti],[xL yL 1 1])/rho0,1),2);
end
end
save([outD model sprintf('_output%03d',output) '_' region '_HBud.mat'],'GWB','-v7.3');

end

%% Vertical Integrate down to level from online T-binned values -----------------------------------------------------------------------------------------------------------
if (doXY)
Tls = [0 5 10 15 20 22.5 25 27.5];
Nremain = length(Tls);
Ti = TL;

FlM = zeros(xL,yL,tL); % vdiffuse and nonlocal_KPP
FlSP = zeros(xL,yL,tL); % solar penetration
FlI = zeros(xL,yL,tL); % numerical mixing
if (haveMIX)
    FlMdif  = zeros(xL,yL,tL); % vdiffuse
    FlMkppiw  = zeros(xL,yL,tL);
    FlMkppish  = zeros(xL,yL,tL);
    FlMkppicon  = zeros(xL,yL,tL);
    FlMkppbl  = zeros(xL,yL,tL);
    FlMkppdd  = zeros(xL,yL,tL);
    FlMwave  = zeros(xL,yL,tL);
end
if (haveSUB)
    FlSUB = zeros(xL,yL,tL);
end
if (haveGM)
    FlGM  = zeros(xL,yL,tL);
end
if (doXYall)
    FlF = zeros(xL,yL,tL); % surface forcing
    FlP = zeros(xL,yL,tL); % P-E+R
    if (haveRedi)
        FlK = zeros(xL,yL,tL); % K33
        FlR = zeros(xL,yL,tL); % Redi
    end
end
 
while (Nremain > 0 & Ti >= 1)
    Tl = Te(Ti);

    for ti=1:tL
        sprintf(['Calculating water-mass heat budget time %03d of ' ...
                 '%03d, temp %2.2f, going down to %2.2f'],ti,tL,Te(Ti),min(Tls))
        FlM(:,:,ti) = FlM(:,:,ti)+ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
            ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        if (haveMDS)
            FlM(:,:,ti) = FlM(:,:,ti)+ncread(wname,'mixdownslope_temp_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        end
        if (haveSIG)
            FlM(:,:,ti) = FlM(:,:,ti)+ncread(wname,'temp_sigma_diff_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        end
        if (haveMIX)
            FlMdif(:,:,ti) = FlMdif(:,:,ti)+ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            FlMkppiw(:,:,ti) = FlMkppiw(:,:,ti)+ncread(wname,'temp_vdiffuse_diff_cbt_kppiw_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            FlMkppish(:,:,ti) = FlMkppish(:,:,ti)+ncread(wname,'temp_vdiffuse_diff_cbt_kppish_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            FlMkppicon(:,:,ti) = FlMkppicon(:,:,ti)+ncread(wname,'temp_vdiffuse_diff_cbt_kppicon_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            FlMkppbl(:,:,ti) = FlMkppbl(:,:,ti)+ncread(wname,'temp_vdiffuse_diff_cbt_kppbl_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            FlMkppdd(:,:,ti) = FlMkppdd(:,:,ti)+ncread(wname,'temp_vdiffuse_diff_cbt_kppdd_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            FlMwave(:,:,ti) = FlMwave(:,:,ti)+ncread(wname,'temp_vdiffuse_diff_cbt_wave_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        end            
        FlSP(:,:,ti) = FlSP(:,:,ti)+ncread(wname,'sw_heat_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        FlI(:,:,ti) = ncread(wname,'temp_numdiff_heat_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        if (haveSUB)
            FlSUB(:,:,ti) = FlSUB(:,:,ti)+ncread(wname,'temp_submeso_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        end
        if (haveGM)
            FlGM(:,:,ti) = FlGM(:,:,ti)+ncread(wname,'neutral_gm_on_nrho_temp',[1 1 Ti ti],[xL yL 1 1]);
        end
        if (doXYall)
            if (haveRedi)
                FlK(:,:,ti) = FlK(:,:,ti)+ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
                FlR(:,:,ti) = FlR(:,:,ti)+ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 Ti ti],[xL yL 1 1]);
            end
            FlP(:,:,ti) = FlP(:,:,ti)+ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
                ncread(wname,'temp_rivermix_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            FlF(:,:,ti) = FlF(:,:,ti)+ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
                ncread(wname,'sw_heat_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
                ncread(wname,'frazil_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
                ncread(wname,'temp_eta_smooth_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        end
    end

    % Save heat flux terms:
    [sp,ind] = min(abs(Tls-Tl));
    if (abs(sp) <= dT/4)
        name = [outD model sprintf('_output%03d',output) '_VertInt_T' strrep(num2str(Tls(ind)),'.','p') 'C.mat']

        save(name,'FlM','FlSP','FlI','Tl','-v7.3');
        if (haveMIX)
            save(name,'FlMdif','FlMkppiw','FlMkppish', ...
                 'FlMkppicon','FlMkppbl','FlMkppdd','FlMwave','-append');
        end
        if (haveSUB)
            save(name,'FlSUB','-append');
        end
        if (haveGM)
            save(name,'FlGM','-append');
        end
        if (doXYall)
            save(name,'FlF','FlP','-append');
            if (haveRedi)
                save(name,'FlK','FlR','-append');
            end
        end
        
        Nremain = Nremain-1;
    end

    Ti = Ti-1;
end

end

%% Calculate WMT due to different (resolved) terms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (doWMT)
Tls = [0:2.5:27.5]-0.25;

for ii = 1:length(Tls)
    Tl = Tls(ii);

    WMTM = NaN*zeros(xL,yL,tL); % vdiffuse and nonlocal_KPP
    WMTF = NaN*zeros(xL,yL,tL); % surface forcing
    WMTI = NaN*zeros(xL,yL,tL); % numerical
    WMTP = NaN*zeros(xL,yL,tL); % P-E+R
    WMTSP = NaN*zeros(xL,yL,tL); % solar penetration
    if (haveRedi)
        WMTK = NaN*zeros(xL,yL,tL); % K33
        WMTR = NaN*zeros(xL,yL,tL); % Redi
    end
    T = ncread(wname,'neutral');
    Te = ncread(wname,'neutralrho_edges');
    [tmp Ti] = min(abs(T-Tl));
    
    for ti=1:tL
        sprintf('Calculating WMT time %03d of %03d, temp %03d of %03d',ti,tL,ii,length(Tls))
        WMTP(:,:,ti) = 1/rho0/Cp/dT*(ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
                                      ncread(wname,'temp_rivermix_on_nrho',[1 1 Ti ti],[xL yL 1 1]));
        WMTM(:,:,ti) = 1/rho0/Cp/dT*(ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
                                      ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 Ti ti],[xL yL 1 1]));
        WMTF(:,:,ti) = 1/rho0/Cp/dT*(ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
                                      ncread(wname,'sw_heat_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
                                      ncread(wname,'frazil_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
                                      ncread(wname,'temp_eta_smooth_on_nrho',[1 1 Ti ti],[xL yL 1 1]));
        WMTSP(:,:,ti) = 1/rho0/Cp/dT*(ncread(wname,'sw_heat_on_nrho',[1 1 Ti ti],[xL yL 1 1]));
        WMTI(:,:,ti) = 1/rho0/Cp/dT*(ncread(wname,'temp_numdiff_heat_on_nrho',[1 1 Ti ti],[xL yL 1 1])-...
                                     ncread(wname,'temp_numdiff_heat_on_nrho',[1 1 Ti+1 ti],[xL yL 1 1]));
        if (haveRedi)
            WMTK(:,:,ti) = 1/rho0/Cp/dT*(ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 Ti ti],[xL yL 1 1]));
            WMTR(:,:,ti) = 1/rho0/Cp/dT*(ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 Ti ti],[xL yL 1 1]));
        end
        if (haveMDS)
            WMTM(:,:,ti) = WMTM(:,:,ti) + 1/rho0/Cp/dT* ...
                ncread(wname,'mixdownslope_temp_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        end
    end
    save([outD model sprintf('_output%03d',output) '_WMT_T' strrep(num2str(Tl),'.','p') 'C.mat'],'WMTM','WMTSP','WMTP','WMTF','WMTI','Tl','-v7.3');
    if (haveRedi)
        save([outD model sprintf('_output%03d',output) '_WMT_T' strrep(num2str(Tl),'.','p') 'C.mat'],'WMTK','WMTR','-append');
    end
end

end
% $$$ 
% $$$ 
% $$$ %% Save isotherm depths -------------------------------------------------------------------------------------
% $$$ Tls = [22 22.5 23];
% $$$ 
% $$$ for ii = 1:length(Tls)
% $$$     Tl = Tls(ii);
% $$$ 
% $$$     ziso = NaN*zeros(xL,yL,tL);
% $$$ 
% $$$     for ti=1:tL
% $$$         sprintf('Calculating isotherm depth time %03d of %03d, temp %03d of %03d',ti,tL,ii,length(Tls))
% $$$         temp = reshape(ncread(fname,'temp',[1 1 1 ti],[xL yL zL 1]),[xL*yL zL 1]);
% $$$         tmax = max(temp,[],2);
% $$$         inds = tmax >= Tl;
% $$$         temp = temp(inds,:);
% $$$         temp(isnan(temp)) = -1000;
% $$$         temp = temp -0.001*repmat(1:zL,[length(temp(:,1)) 1]);
% $$$         zisot = zeros(length(temp(:,1)),1);
% $$$         for jj=1:length(temp(:,1))
% $$$             zisot(jj) = interp1(temp(jj,:),z,Tl,'linear');
% $$$         end
% $$$         tmp = NaN*zeros(xL,yL);
% $$$         tmp(inds) = zisot;
% $$$         ziso(:,:,ti) = tmp;
% $$$     end
% $$$     save([outD model sprintf('_output%03d',output) '_Ziso_T' strrep(num2str(Tl),'.','p') 'C.mat'],'ziso');
% $$$ end
% $$$ 

%% Save surface heat flux, wind stress, SST, meridional heat flux:
if (doSURF)
try
    shflux = ncread(fname,'net_sfc_heating',[1 1 1],[xL yL tL]);
catch
    shflux = ncread(fname2,'net_sfc_heating',[1 1 1],[xL yL tL]);
end    
SST = squeeze(ncread(fname,'temp',[1 1 1 1],[xL yL 1 tL]));
% $$$ taux = ncread(fname,'tau_x',[1 1 1],[xL yL tL]);
% $$$ tauy = ncread(fname,'tau_y',[1 1 1],[xL yL tL]);

save([outD model sprintf('_output%03d',output) '_SurfaceVars.mat'],'shflux','SST');%,'taux','tauy');
end

%% Save vertically-integrated heat transports:
if (doXYtran)
Tls = [10 12.5 15 20 34];
Nremain = length(Tls);
Ti = 1;

xflux = zeros(xL,yL); % vdiffuse and nonlocal_KPP
yflux = zeros(xL,yL); % solar penetration

while (Nremain > 0 & Ti <= TL)
    Tl = Te(Ti+1);

    qxtrans = zeros(xL,yL);
    qytrans = zeros(xL,yL);
    for ti=1:tL
        sprintf(['Calculating vertically-integrated heat fluxes time %03d of ' ...
                 '%03d, temp %2.2f, going up to %2.2f'],ti,tL,Te(Ti),max(Tls))

        qxtrans = qxtrans+ncread(wname,'temp_xflux_adv_on_nrho',[1 1 Ti ti],[xL yL 1 1])*ndays(ti);
        qytrans = qytrans+ncread(wname,'temp_yflux_adv_on_nrho',[1 1 Ti ti],[xL yL 1 1])*ndays(ti);

        % Submesoscale and GM: two options:
        if (haveSUB)
            qxtrans = qxtrans + ncread(wname,'temp_xflux_submeso_on_nrho',[1 1 Ti ti],[xL yL 1 1])*ndays(ti);
            qytrans = qytrans + ncread(wname,'temp_yflux_submeso_on_nrho',[1 1 Ti ti],[xL yL 1 1])*ndays(ti);
        end
        if (haveGM)
            qxtrans = qxtrans + ncread(wname,'temp_xflux_gm_on_nrho',[1 1 Ti ti],[xL yL 1 1])*ndays(ti);
            qytrans = qytrans + ncread(wname,'temp_yflux_gm_on_nrho',[1 1 Ti ti],[xL yL 1 1])*ndays(ti);
        end
    end
    xflux = xflux + qxtrans/sum(ndays);
    yflux = yflux + qytrans/sum(ndays);

    % Save heat flux terms:
    [sp,ind] = min(abs(Tls-Tl));
    if (abs(sp) <= dT/4)
        name = [outD model sprintf('_output%03d',output) '_XYtrans_T' strrep(num2str(Tls(ind)),'.','p') 'C.mat']

        save(name,'xflux','yflux','Tl','-v7.3');
        Nremain = Nremain-1;
    end

    Ti = Ti+1;
end
end

%% Zonally-averaged fluxes -------------------------------------------------------------
if (doZA)
ZA.F = zeros(yL,TL+1,tL); % Wdeg-1 due to F
ZA.P = zeros(yL,TL+1,tL); % Wdeg-1 due to P
ZA.M = zeros(yL,TL+1,tL); % Wdeg-1 due to M
ZA.KPPNL = zeros(yL,TL+1,tL); % Wdeg-1 due to KPP non-local
ZA.dVdt = zeros(yL,TL+1,tL); % Wdeg-1 dVdt
ZA.dHdt = zeros(yL,TL+1,tL); % Wdeg-1 dHdt
ZA.SWH = zeros(yL,TL+1,tL); % Wdeg-1 due to SW redistribution
if (haveMIX)
    ZA.Mkppiw = zeros(yL,TL+1,tL); %Wdeg-1 due to kppiw;
    ZA.Mkppish = zeros(yL,TL+1,tL); %Wdeg-1 due to kppiw;
    ZA.Mwave = zeros(yL,TL+1,tL); %Wdeg-1 due to kppiw;
    ZA.Mkppbl = zeros(yL,TL+1,tL); %Wdeg-1 due to kppiw;
    ZA.Moth   = zeros(yL,TL+1,tL); %Wdeg-1 due to kppicon+kppdd+KPPnloc;
end
if (haveRedi)
    ZA.K33    = zeros(yL,TL+1,tL); % W due to K33
    ZA.RED    = zeros(yL,TL+1,tL); % W due to Redi diffusion
    ZA.AHDR   = zeros(yL,TL+1,tL); % Meridional heat flux due to
                                   % Redi diffusion
end
ZA.JS = zeros(yL,TL+1,tL); % m3s-1deg-1 due to surface volume flux
ZA.PSI = zeros(yL,TL+1,tL); % m3s-1 northward transport
ZA.AHD = zeros(yL,TL+1,tL); % W A direct using heat fluxes
if (haveSUB)
    ZA.SUB    = zeros(yL,TL+1,tL);
    ZA.PSISUB = zeros(yL,TL+1,tL);
    ZA.AHDSUB = zeros(yL,TL+1,tL);
end
if (haveGM)
    ZA.GM    = zeros(yL,TL+1,tL);
    ZA.PSIGM = zeros(yL,TL+1,tL);
    ZA.AHDGM = zeros(yL,TL+1,tL);
end
if (doHND)
    ZA.NUM   = zeros(yL,TL+1,tL);
end
if (haveMDS)
    ZA.MDS = zeros(yL,TL+1,tL); % W
end
if (haveSIG)
    ZA.SIG = zeros(yL,TL+1,tL); % W
end
% ignoring tri-polar for now.
for ti=1:tL
    if (doHND)
        for ii = 1:(TL+1)
            ZA.NUM(:,ii,ti) = nansum(tmaskREG.*area.*ncread(wname,'temp_numdiff_heat_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)';
        end
    end

    for ii = 1:TL
    sprintf('Calculating zonally-averaged water-mass heat budget time %03d of %03d, temp %03d of %03d',ti,tL,ii,TL)
    ZA.F(:,ii+1,ti) = ZA.F(:,ii,ti) + (nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL 1 1]),1) + ...
                    nansum(tmaskREG.*area.*ncread(wname,'frazil_on_nrho',[1 1 ii ti],[xL yL 1 1]),1) + ...
                    nansum(tmaskREG.*area.*ncread(wname,'temp_eta_smooth_on_nrho',[1 1 ii ti],[xL yL 1 1]),1))';
    ZA.P(:,ii+1,ti) = ZA.P(:,ii,ti) + (nansum(tmaskREG.*area.*ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 ii ti],[xL yL 1 1]),1) + ...
                    nansum(tmaskREG.*area.*ncread(wname,'temp_rivermix_on_nrho',[1 1 ii ti],[xL yL 1 1]),1))';
    ZA.M(:,ii+1,ti) = ZA.M(:,ii,ti) + nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)';
    ZA.KPPNL(:,ii+1,ti) = ZA.KPPNL(:,ii,ti) + nansum(tmaskREG.*area.*ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)';
    if (haveMDS)
        ZA.MDS(:,ii+1,ti) = ZA.MDS(:,ii,ti) + nansum(tmaskREG.*area.*ncread(wname,'mixdownslope_temp_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)';
    end
    if (haveSIG)
        ZA.SIG(:,ii+1,ti) = ZA.SIG(:,ii,ti) + nansum(tmaskREG.*area.*ncread(wname,'temp_sigma_diff_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)';
    end
    ZA.dVdt(:,ii+1,ti) = ZA.dVdt(:,ii,ti) + nansum(tmaskREG.*ncread(wname,'dVdt',[1 1 ii ti],[xL yL 1 1])*1e9/rho0,1)';
    ZA.dHdt(:,ii+1,ti) = ZA.dHdt(:,ii,ti) + nansum(tmaskREG.*ncread(wname,'dHdt',[1 1 ii ti],[xL yL 1 1]),1)';
    ZA.PSI(:,ii+1,ti)  = ZA.PSI(:,ii,ti) + nansum(umaskREG.*ncread(wname,'ty_trans_nrho',[1 1 ii ti],[xL yL 1 1]))'*tsc/rho0;
    ZA.AHD(:,ii+1,ti)  = ZA.AHD(:,ii,ti) + nansum(umaskREG.*ncread(wname,'temp_yflux_adv_on_nrho',[1 1 ii ti],[xL yL 1 1]))';
    if (haveSUB)
        ZA.SUB(:,ii+1,ti) = ZA.SUB(:,ii,ti) + nansum(tmaskREG.*ncread(wname,'temp_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)';
        ZA.PSISUB(:,ii+1,ti) = ZA.PSISUB(:,ii,ti) + nansum(umaskREG.*ncread(wname,'ty_trans_nrho_submeso',[1 1 ii ti],[xL yL 1 1]),1)'*tsc/rho0;
        ZA.AHDSUB(:,ii+1,ti) = ZA.AHDSUB(:,ii,ti) + nansum(umaskREG.*ncread(wname,'temp_yflux_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)';
    end
    if (haveGM) 
        ZA.GM(:,ii+1,ti) = ZA.GM(:,ii,ti) + nansum(tmaskREG.*ncread(wname,'neutral_gm_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1)';
        ZA.PSIGM(:,ii+1,ti) = ZA.PSIGM(:,ii,ti) + nansum(umaskREG.*ncread(wname,'ty_trans_nrho_gm',[1 1 ii ti],[xL yL 1 1]),1)'*tsc/rho0;
        ZA.AHDGM(:,ii+1,ti) = ZA.AHDGM(:,ii,ti) + nansum(umaskREG.*ncread(wname,'temp_yflux_gm_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)';
    end

    ZA.SWH(:,ii+1,ti) = ZA.SWH(:,ii,ti) + nansum(tmaskREG.*area.*ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)';
    ZA.F(:,ii+1,ti) = ZA.F(:,ii+1,ti) + nansum(tmaskREG.*area.*ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)';
    if (haveMIX)
        ZA.Mkppiw(:,ii+1,ti) = ZA.Mkppiw(:,ii,ti) + (nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppiw_on_nrho',[1 1 ii ti],[xL yL 1 1]),1))';
        ZA.Mkppish(:,ii+1,ti) = ZA.Mkppish(:,ii,ti) + (nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppish_on_nrho',[1 1 ii ti],[xL yL 1 1]),1))';
        ZA.Mwave(:,ii+1,ti) = ZA.Mwave(:,ii,ti) + (nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_wave_on_nrho',[1 1 ii ti],[xL yL 1 1]),1))';
        ZA.Mkppbl(:,ii+1,ti) = ZA.Mkppbl(:,ii,ti) + (nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppbl_on_nrho',[1 1 ii ti],[xL yL 1 1]),1))';
        ZA.Moth(:,ii+1,ti) = ZA.Moth(:,ii,ti) + (nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppicon_on_nrho',[1 1 ii ti],[xL yL 1 1]),1) + ...
                            nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppdd_on_nrho',[1 1 ii ti],[xL yL 1 1]),1) + ...
                            nansum(tmaskREG.*area.*ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]),1))';
    end
    if (haveRedi)
        ZA.K33(:,ii+1,ti) = ZA.K33(:,ii,ti) + (nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 ii ti],[xL yL 1 1]),1))';
        ZA.RED(:,ii+1,ti) = ZA.RED(:,ii,ti) + (nansum(tmaskREG.*area.*ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1))';
        ZA.AHDR(:,ii+1,ti) = ZA.AHDR(:,ii,ti) + nansum(umaskREG.*ncread(wname,'temp_yflux_ndiffuse_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)';
    end        
    ZA.JS(:,ii+1,ti) = ZA.JS(:,ii,ti) + (nansum(tmaskREG.*ncread(wname,'mass_pmepr_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)/rho0)';
    end
    save([outD model sprintf('_output%03d',output) '_' region '_ZAHBud.mat'],'ZA','yu','yt','-v7.3');
end
end % End doZA


% $$$ %% Zonally-averaged diffusivity (rough) --------------------------------------------------
% $$$ 
% $$$ rho0 = 1035; % kgm-3
% $$$ area = ncread('ocean_grid.nc','area_t');
% $$$ area(isnan(area)) = 0;
% $$$ [xL,yL] = size(area);
% $$$ T = ncread('ocean.nc','neutral');
% $$$ yt = ncread('ocean.nc','yt_ocean');
% $$$ TL = length(T);
% $$$ tL = 12;
% $$$ diff_cbt_T_G = zeros(yL,TL,tL);
% $$$ diff_cbt_T_P = zeros(yL,TL,tL);
% $$$ diff_cbt_T_A = zeros(yL,TL,tL);
% $$$ 
% $$$ [tmaskIP,~,~] = Heat_Budget_Mask('IndoPacific2BAS','ocean_grid.nc','ocean_wmass.nc','../mat_data/','MOM025_kb3seg');
% $$$ [tmaskA,~,~] = Heat_Budget_Mask('Atlantic2BAS','ocean_grid.nc','ocean_wmass.nc','../mat_data/','MOM025_kb3seg');
% $$$ 
% $$$ for i=1:tL
% $$$     for ii=1:TL
% $$$         [i tL ii TL]
% $$$         diff = ncread('ocean.nc','diff_cbt_t_on_nrho',[1 1 ii i],[xL ...
% $$$                             yL 1 1]);
% $$$         mass = ncread('ocean.nc','mass_t_on_nrho',[1 1 ii i],[xL ...
% $$$                             yL 1 1]);
% $$$         diff(isnan(diff)) = 0;
% $$$         mass(isnan(mass)) = 0;
% $$$         
% $$$         diff_cbt_T_G(:,ii,i) = sum(diff.*area,1)./sum(mass/rho0,1);
% $$$         diff_cbt_T_P(:,ii,i) = sum(diff.*area.*tmaskIP,1)./sum(mass.*tmaskIP/rho0,1);
% $$$         diff_cbt_T_A(:,ii,i) = sum(diff.*area.*tmaskA,1)./sum(mass.*tmaskA/rho0,1);
% $$$     end
% $$$ end
% $$$ save('../mat_data/MOM025_kb3seg_output121_diff_cbt_T','diff_cbt_T_G','diff_cbt_T_A','diff_cbt_T_P','T','yt');
% $$$ SST_G = squeeze(ncread('ocean.nc','temp',[1 1 1 1],[xL yL 1 tL]));
% $$$ SST_P = squeeze(ncread('ocean.nc','temp',[1 1 1 1],[xL yL 1 ...
% $$$                     tL])).*repmat(tmaskIP,[1 1 tL]);
% $$$ SST_A = squeeze(ncread('ocean.nc','temp',[1 1 1 1],[xL yL 1 tL])).*repmat(tmaskA,[1 1 tL]);
% $$$ save('../mat_data/MOM025_kb3seg_output121_SurfaceVars.mat','SST_G','SST_A','SST_P');

% $$$ %% Annual average/max, zonal-average/max depth of isotherms:
% $$$ 
% $$$ tempZA = zeros(yL,zL);
% $$$ saltZA = zeros(yL,zL);
% $$$ tempZM = -100*zeros(yL,zL);
% $$$ tempZMA = -100*zeros(yL,zL);
% $$$ 
% $$$ for ti=1:tL
% $$$     sprintf('Calculating zonally-averaged temperature %03d of %03d',ti,tL)
% $$$     temp = ncread(fname,'temp',[1 1 1 ti],[xL yL zL 1]);
% $$$     salt = ncread(fname,'salt',[1 1 1 ti],[xL yL zL 1]);
% $$$     areaNAN = repmat(area,[1 1 zL]).*(~isnan(temp));
% $$$     
% $$$     tempZ = squeeze(nansum(areaNAN.*temp,1))./squeeze(nansum(areaNAN,1));
% $$$     tempZA = tempZA + tempZ;
% $$$ 
% $$$     saltZ = squeeze(nansum(areaNAN.*salt,1))./squeeze(nansum(areaNAN,1));
% $$$     saltZA = saltZA + saltZ;
% $$$ 
% $$$     maxt = squeeze(max(temp,[],1));
% $$$     tempZM = max(tempZM,maxt);
% $$$     tempZMA = max(tempZMA,tempZ);
% $$$ end
% $$$ tempZA = tempZA/tL;
% $$$ saltZA = saltZA/tL;
% $$$ 
% $$$ % Depth of isotherms:
% $$$ ZAtemp = zeros(yL,TL+1);
% $$$ ZMtemp = zeros(yL,TL+1);
% $$$ ZMAtemp = zeros(yL,TL+1);
% $$$ 'Calculating zonally-averaged isotherms'
% $$$ for yi=1:yL
% $$$     tvec = squeeze(tempZA(yi,:));
% $$$     zvec = -z;
% $$$     tvec(isnan(tvec)) = -1000;
% $$$     tvec = tvec - 0.01*(1:zL);
% $$$     ZAtemp(yi,:) = interp1(tvec,zvec,Te,'linear');
% $$$     ind = find(~isnan(ZAtemp(yi,:)),1,'last');
% $$$     ZAtemp(yi,(ind+1):end) = max(zvec);
% $$$ 
% $$$     tvec = squeeze(tempZM(yi,:));
% $$$     tvec(isnan(tvec)) = -1000;
% $$$     tvec = tvec - 0.01*(1:zL);
% $$$     ZMtemp(yi,:) = interp1(tvec,zvec,Te,'linear');
% $$$     ind = find(~isnan(ZMtemp(yi,:)),1,'last');
% $$$     ZMtemp(yi,(ind+1):end) = max(zvec);
% $$$ 
% $$$     tvec = squeeze(tempZMA(yi,:));
% $$$     tvec(isnan(tvec)) = -1000;
% $$$     tvec = tvec - 0.01*(1:zL);
% $$$     ZMAtemp(yi,:) = interp1(tvec,zvec,Te,'linear');
% $$$     ind = find(~isnan(ZMAtemp(yi,:)),1,'last');
% $$$     ZMAtemp(yi,(ind+1):end) = max(zvec);
% $$$ end
% $$$ 
% $$$ save([outD model sprintf('_output%03d',output) '_ZAHBud.mat'],'tempZA','saltZA','ZAtemp','tempZM','ZMtemp','tempZMA','ZMAtemp','z','Te','latv','-append');

% $$$ end



% $$$ end
% $$$ %% Swap in non-NaN'd lon/lat:
% $$$ base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
% $$$ model = 'MOM025';
% $$$ 
% $$$ load([base model sprintf('_output%03d_BaseVars.mat',8)]);
% $$$ region = 'Global';
% $$$ 
% $$$ base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
% $$$ model = 'MOM025_kb3seg';
% $$$ for output = [87:90]
% $$$     save([base model sprintf('_output%03d_BaseVars.mat',output)], ...
% $$$          'lon','lat','lonu','latu','area','-append');
% $$$ end
% $$$ 
