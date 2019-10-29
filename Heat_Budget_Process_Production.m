% This script processes the heat budget and associated variables in
% MOM025 or MOM01 simulations and save's into .mat files

baseL = '/short/e14/rmh561/mom/archive/';

% MOM-SIS025:
model = 'MOM025_kb3seg';
baseD = [baseL 'MOM_HeatDiag_kb3seg/']; %Data Directory.
ICdir = [baseL 'MOM_HeatDiag_kb3seg/restart100/'];

outD = [baseD 'mat_data/'];
rstbaseD = baseD;

post = ''; % For MOM-SIS.

% term options:
haveRedi = 1; % 1 = Redi diffusion is on, 0 = off
haveGM = 1; % 1 = GM is on, 0 = off;
haveSUB = 1; % 1 = submeso is on, 0 = off;
haveMDS = 0; % 1 = MDS is on, 0 = off;
haveSIG = 0; % 1 = SIG is on, 0 = off;
haveMIX = 0; % 1 = Do mixing components (vdiffuse_diff_cbt_*), 0 = don't. 

% Processing options:
doBASE     = 1; % 1 = save BaseVars.mat file
dodVdtdHdt = 1; % 1 = calculate dVdt/dHdt and save into .nc file
doGWB      = 1; % 1 = calculate global online budget terms
doSURF     = 1; % 1 = calculate surface flux field and SST
doZA       = 1; % 1 = calculate zonal average budget

% scaling constant on the transports:
if (strcmp(model(1),'A')) %ACCESS-OM2, transport in kg/s
    tsc = 1;
else % MOM-SIS, transport in 1e9 kg/s
    tsc = 1e9;
end

for output = 101:110;
restart = output-1;

% CHOOSE REGION HERE!!!!
region = 'Global';
%region = 'IndoPacific2BAS';
%region = 'Atlantic2BAS';

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

end

