% This script processes the heat budget and associated variables in
% MOM025 or MOM01 simulations and save's into .mat files

baseL = '/short/e14/rmh561/access-om2/archive/';
% ACCESS-OM2:
model = 'ACCESS-OM2_1deg_jra55_ryf8485_kds50_may';
baseD = [baseL '1deg_jra55_ryf8485_kds50_may/']; %Data Directory.
ICdir = '/g/data1/ua8/MOM/initial_conditions/WOA/10_KDS50/';

outD = [baseD 'mat_data/'];
rstbaseD = baseD;

post = 'ocean/'; % For ACCESS-OM2 output coulpled;

haveRedi = 0; % 1 = Redi diffusion is on, 0 = off
haveGM = 0; % 1 = GM is on, 0 = off;
haveMDS = 0; % 1 = MDS is on, 0 = off;

output=86;
restart = output-1;

region = 'Global';

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
if (exist(baser))
    found_rst = 1;rstti = 1;
    rnameT = [baser 'ocean_temp_salt.res.nc'];
    rnameZ = [baser 'ocean_thickness.res.nc'];
else
    found_rst = 0;rstti = 12;
    rnameT = [basem1 'ocean_snap.nc'];
    rnameZ = [basem1 'ocean_snap.nc'];
end
         
% Horizontal Grid  -----------------------------------------
lon = ncread(gname,'geolon_t');lat = ncread(gname,'geolat_t');
lonu = ncread(gname,'geolon_c');latu = ncread(gname,'geolat_c');
area = ncread(gname,'area_t');[xL,yL] = size(lon);
lonv_t = ncread(gname,'xt_ocean');lonv_u = ncread(gname,'xu_ocean');
latv_t = ncread(gname,'yt_ocean');latv_u = ncread(gname,'yu_ocean');

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
mask = ncread(fname,'temp',[1 1 1 rstti],[xL yL zL 1]);
mask(~isnan(mask)) = 1; mask(isnan(mask)) = 0;
mask = mask == 1;

if (~strcmp(region,'Global'))
    [maskREG,~,~,~,~,~,~] = Heat_Budget_Mask(region,gname,fname,wname,outD,model);
else
    maskREG = ones(xL,yL);
end

% Time  -----------------------------------------
time = ncread(fname,'time');

dys = [31 28 31 30 31 30 31 31 30 31 30 31];
C = textread(tname, '%s','delimiter', '\n');
C = strsplit(C{1});
rtime = [str2num(C{1}) str2num(C{2}) str2num(C{3}) str2num(C{4}) str2num(C{5}) str2num(C{6})];
time_snap = [(rtime(1)-1)*365+sum(dys(1:(rtime(2)-1)))+(rtime(3)-1)+rtime(4)/24+rtime(5)/24/60+rtime(6)/24/60/60;
             ncread(sname,'time')];

tL = length(time);

Cp = 3992.1; % J kg-1 degC-1
rho0 = 1035; % kgm-3

T = ncread(wname,'neutral');
Te = ncread(wname,'neutralrho_edges');
TL = length(T);dT = T(2)-T(1);

%latitude vector for heat function:
latv = max(lat,[],1);
late = [-90 (latv(2:end)+latv(1:(end-1)))/2 90];

save([outD model sprintf('_output%03d',output) '_BaseVars.mat'], ...
     'T','Te','TL','dT','Cp','rho0','time','time_snap','tL', ...
     'z','zL','lon','lat','area','xL','yL','latv','late', ...
     'lonu','latu','-v7.3');

%% Calculate dVdt, dHdt and save back into wmass file:
%% Also, calculate TENMON

% Create variables:
ncid = netcdf.open(wname,'NC_WRITE');
try
    id = netcdf.inqVarID(ncid,'dVdt');
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
    netcdf.putAtt(ncid,dVdtID,'_FillValue',single(-1e20));
    dHdtID = netcdf.defVar(ncid,'dHdt','NC_FLOAT',[xid yid zid tid]);
    netcdf.putAtt(ncid,dHdtID,'long_name','Change in time of heat content within temperature bin');
    netcdf.putAtt(ncid,dHdtID,'units','Watts');
    netcdf.putAtt(ncid,dHdtID,'_FillValue',single(-1e20));
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

TENMON = zeros(TL+1,tL);

%Do other times for Vsnap and Hsnap:
for ti=1:tL
    for zi=1:zL
        sprintf('Calculating Vsnap and Hsnap later months time %03d of %03d, depth %02d of %02d',ti,tL,zi,zL)

        temp = ncread(fname,'temp',[1 1 zi ti],[xL yL 1 1]);
        temp(~mask(:,:,zi)) = NaN;
        if (max(max(temp))>120);temp = temp-273.15;end;

        tempsnap = ncread(sname,'temp',[1 1 zi ti],[xL yL 1 1]);
        tempsnap(~mask(:,:,zi)) = NaN;
        if (max(max(tempsnap))>120);tempsnap = tempsnap-273.15;end;
        Volsnap = ncread(sname,'dzt',[1 1 zi ti],[xL yL 1 1]).*area;
        Volsnap(isnan(Volsnap)) = 0;
        
        TENf = area.*ncread(hname,'temp_tendency',[1 1 zi ti],[xL ...
                            yL 1 1]);

        for Ti=1:TL
            %Accumulate sums:
            inds = tempsnap>=Te(Ti) & tempsnap<Te(Ti+1);
            Vsnap(:,:,Ti) = Vsnap(:,:,Ti) + Volsnap.*inds;
            Hlay = Volsnap.*tempsnap.*inds*rho0*Cp;
            Hlay(isnan(Hlay)) = 0;
            Hsnap(:,:,Ti) = Hsnap(:,:,Ti) + Hlay;
            inds = find(temp>=Te(Ti) & temp<Te(Ti+1));
            TENMON(Ti,ti) = TENMON(Ti,ti)+nansum(TENf(inds));
        end
        inds = find(temp>=Te(TL+1));
        TENMON(TL+1,ti) = TENMON(TL+1,ti)+nansum(TENf(inds));
    end
    
    % Integrate to get to T'>T:
    TENMON(:,ti) = flipud(cumsum(flipud(TENMON(:,ti))));
    
    netcdf.putVar(ncid,dVdtID,[0 0 0 ti-1],[xL yL TL 1],(Vsnap-VsnapM) ...
                  /(time_snap(ti+1)-time_snap(ti))/86400*rho0/1e9);
    netcdf.putVar(ncid,dHdtID,[0 0 0 ti-1],[xL yL TL 1],(Hsnap-HsnapM) ...
                  /(time_snap(ti+1)-time_snap(ti))/86400);
    VsnapM = Vsnap;
    HsnapM = Hsnap;
    Vsnap = zeros(xL,yL,TL);
    Hsnap = zeros(xL,yL,TL);
end
netcdf.close(ncid);

%% Calculate volume integrated budget from online T-binned values -----------------------------------------------------------------------------------------------------------

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
GWB.SUB    = zeros(TL+1,tL); % W due to submesoscale.
GWB.VDF    = zeros(TL+1,tL); % W due to vdiffusion
GWB.KNL    = zeros(TL+1,tL); % W due to KPP non-local
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
GWB.ADV    = zeros(TL+1,tL); % W due to advection
GWB.TEN    = zeros(TL+1,tL); % W due to tendency
GWB.SFW    = zeros(TL+1,tL); % surface volume flux into ocean (m3s-1)

if (haveHND)
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
    GWB.dVdt(ii,ti) = nansum(nansum(maskREG.*ncread(wname,'dVdt',[1 1 ii ti],[xL yL 1 1])*1e9/rho0,1),2);
    GWB.dHdt(ii,ti) = nansum(nansum(maskREG.*ncread(wname,'dHdt',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.TEN(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_tendency_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.ADV(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_advection_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.SUB(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.PME(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.RMX(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_rivermix_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDS(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.SWH(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDF(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.KNL(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    if (haveMIX)
    GWB.VDFkppiw(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppiw_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFkppish(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppish_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFkppicon(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppicon_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFkppbl(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppbl_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFkppdd(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppdd_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFwave(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_wave_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    end        
    if (haveRedi)
    GWB.K33(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.RED(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveGM)
    GWB.NGM(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'neutral_gm_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveMDS)
    GWB.MDS(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'mixdownslope_temp_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    GWB.FRZ(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'frazil_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.ETS(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_eta_smooth_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.SFW(ii,ti) = nansum(nansum(maskREG.*ncread(wname,'mass_pmepr_on_nrho',[1 1 ii ti],[xL yL 1 1])/rho0,1),2);
for ii=TL-1:-1:1
    sprintf('Calculating global water-mass heat budget time %03d of %03d, temp %03d of %03d',ti,tL,ii,TL)
    GWB.dVdt(ii,ti) = GWB.dVdt(ii+1,ti) + nansum(nansum(maskREG.*ncread(wname,'dVdt',[1 1 ii ti],[xL yL 1 1])*1e9/rho0,1),2);
    GWB.dHdt(ii,ti) = GWB.dHdt(ii+1,ti) + nansum(nansum(maskREG.*ncread(wname,'dHdt',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.TEN(ii,ti) = GWB.TEN(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'temp_tendency_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.ADV(ii,ti) = GWB.ADV(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'temp_advection_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.SUB(ii,ti) = GWB.SUB(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'temp_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.PME(ii,ti) = GWB.PME(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.RMX(ii,ti) = GWB.RMX(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'temp_rivermix_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDS(ii,ti) = GWB.VDS(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.SWH(ii,ti) = GWB.SWH(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDF(ii,ti) = GWB.VDF(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.KNL(ii,ti) = GWB.KNL(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    if (haveMIX)
    GWB.VDFkppiw(ii,ti)   = GWB.VDFkppiw(ii+1,ti)   + nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppiw_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFkppish(ii,ti)  = GWB.VDFkppish(ii+1,ti)  + nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppish_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFkppicon(ii,ti) = GWB.VDFkppicon(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppicon_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFkppbl(ii,ti)   = GWB.VDFkppbl(ii+1,ti)   + nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppbl_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFkppdd(ii,ti)   = GWB.VDFkppdd(ii+1,ti)   + nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppdd_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFwave(ii,ti)    = GWB.VDFwave(ii+1,ti)    + nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_wave_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);        
    end
    if (haveRedi)
    GWB.K33(ii,ti) = GWB.K33(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.RED(ii,ti) = GWB.RED(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveGM)
    GWB.NGM(ii,ti) = GWB.NGM(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'neutral_gm_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveMDS)
    GWB.MDS(ii,ti) = GWB.MDS(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'mixdownslope_temp_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    GWB.FRZ(ii,ti) = GWB.FRZ(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'frazil_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.ETS(ii,ti) = GWB.ETS(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'temp_eta_smooth_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.SFW(ii,ti) = GWB.SFW(ii+1,ti) + nansum(nansum(maskREG.*ncread(wname,'mass_pmepr_on_nrho',[1 1 ii ti],[xL yL 1 1])/rho0,1),2);
end
end
save([outD model sprintf('_output%03d',output) '_' region 'HBud.mat'],'GWB','-v7.3');

%% Vertical Integrate down to level from online T-binned values -----------------------------------------------------------------------------------------------------------
Tls = [0:2.5:27.5];
Nremain = length(Tls);
Ti = TL;

FlM = zeros(xL,yL,tL); % vdiffuse and nonlocal_KPP
FlF = zeros(xL,yL,tL); % surface forcing
FlP = zeros(xL,yL,tL); % P-E+R
% $$$ FlA = zeros(xL,yL,tL); % advection + submeso + GM
if (haveMIX)
    FlMdif  = zeros(xL,yL,tL); % vdiffuse
    FlMkppiw  = zeros(xL,yL,tL);
    FlMkppish  = zeros(xL,yL,tL);
    FlMkppicon  = zeros(xL,yL,tL);
    FlMkppbl  = zeros(xL,yL,tL);
    FlMkppdd  = zeros(xL,yL,tL);
    FlMwave  = zeros(xL,yL,tL);
end
if (haveRedi)
    FlK = zeros(xL,yL,tL); % K33
    FlR = zeros(xL,yL,tL); % Redi
end
% $$$ if (haveGM)
% $$$     FlG = zeros(xL,yL,tL); % GM
% $$$ end
% $$$ FlT = zeros(xL,yL,tL); % tendency
FlSP = zeros(xL,yL,tL); % solar penetration
FlI = zeros(xL,yL,tL); % numerical mixing
 
while (Nremain > 0 & Ti >= 1)
    Tl = Te(Ti);

    for ti=1:tL
        sprintf(['Calculating water-mass heat budget time %03d of ' ...
                 '%03d, temp %2.2f, going down to %2.2f'],ti,tL,Te(Ti),min(Tls))
% $$$         FlT(:,:,ti) = FlT(:,:,ti)+ncread(wname,'temp_tendency_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        if (haveRedi)
            FlK(:,:,ti) = FlK(:,:,ti)+ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            FlR(:,:,ti) = FlR(:,:,ti)+ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 Ti ti],[xL yL 1 1]);
        end
% $$$         FlA(:,:,ti) = FlA(:,:,ti)+ncread(wname,'temp_advection_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
% $$$             ncread(wname,'temp_submeso_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
% $$$         if (haveGM)
% $$$             FlG(:,:,ti) = FlG(:,:,ti)+ncread(wname,'neutral_gm_on_nrho_temp',[1 1 Ti ti],[xL yL 1 1]);
% $$$             FlA(:,:,ti) = FlA(:,:,ti)+ncread(wname,'neutral_gm_on_nrho_temp',[1 1 Ti ti],[xL yL 1 1]);
% $$$         end
        FlP(:,:,ti) = FlP(:,:,ti)+ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
            ncread(wname,'temp_rivermix_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        FlM(:,:,ti) = FlM(:,:,ti)+ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
            ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        if (haveMDS)
            FlM(:,:,ti) = FlM(:,:,ti)+ncread(wname,'mixdownslope_temp_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
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
        FlF(:,:,ti) = FlF(:,:,ti)+ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
            ncread(wname,'sw_heat_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
            ncread(wname,'frazil_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
            ncread(wname,'temp_eta_smooth_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        FlSP(:,:,ti) = FlSP(:,:,ti)+ncread(wname,'sw_heat_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        FlI(:,:,ti) = ncread(wname,'temp_numdiff_heat_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
    end

    % Save heat flux terms:
    [sp,ind] = min(abs(Tls-Tl));
    if (abs(sp) <= dT/4)
        name = [outD model sprintf('_output%03d',output) '_VertInt_T' strrep(num2str(Tls(ind)),'.','p') 'C.mat']

        save(name,'FlM','FlSP','FlF','FlP','FlI','Tl','-v7.3');
% $$$         save(name,'FlM','FlSP','FlF','FlT','FlA','FlP','FlI','FlIH','Tl','-v7.3');
        if (haveMIX)
            save(name,'FlMdif','FlMkppiw','FlMkppish', ...
                 'FlMkppicon','FlMkppbl','FlMkppdd','FlMwave','-append');
        end
        if (haveRedi)
            save(name,'FlK','FlR','-append');
        end
% $$$         if (haveGM)
% $$$             save(name,'FlG','-append');
% $$$         end
        
        Nremain = Nremain-1;
    end

    Ti = Ti-1;
end

%% Save surface heat flux, wind stress, SST, meridional heat flux:
try
    shflux = ncread(fname,'net_sfc_heating',[1 1 1],[xL yL tL]);
catch
    shflux = ncread(fname2,'net_sfc_heating',[1 1 1],[xL yL tL]);
end    
SST = squeeze(ncread(fname,'temp',[1 1 1 1],[xL yL 1 tL]));
% $$$ taux = ncread(fname,'tau_x',[1 1 1],[xL yL tL]);
% $$$ tauy = ncread(fname,'tau_y',[1 1 1],[xL yL tL]);

save([outD model sprintf('_output%03d',output) '_SurfaceVars.mat'],'shflux','SST');%,'taux','tauy');

% Do meridional heat flux:
mhflux = zeros(yL,tL);
for ti = 1:tL
    for Ti = 1:TL
        sprintf('Calculating meridional heat flux time %03d of %03d, temp %03d of %03d',ti,tL,Ti,TL)
        tytrans = ncread(wname,'ty_trans_nrho',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0 + ...
                  ncread(wname,'ty_trans_nrho_submeso',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
        if (haveGM)
            tytrans = tytrans + ncread(wname,'ty_trans_nrho_gm',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
        end
        mhflux(:,ti) = mhflux(:,ti) + rho0*Cp*T(Ti)*nansum(tytrans,1)';
    end
end
save([outD model sprintf('_output%03d',output) '_SurfaceVars.mat'],'mhflux','-append');

%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%

% % Plotting Code:

% This script makes plots of the heat budget in the MOM
% simulations.

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
% $$$ base = 'archive/mat_data/';

RUNS = { ...
         {'ACCESS-OM2_1deg_jra55_ryf8485_kds50_may',[37]}, ...
       };

rr = 1;
for rr = 1:length(RUNS);
    outputs = RUNS{rr}{2};
    model = RUNS{rr}{1};

    clearvars -except base RUNS rr outputs model leg legh;
    
    load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
    ndays = diff(time_snap);
    ndays = ndays(1:12);
    if (ndays(end) <= 0); ndays(end) = 365-ndays(end);end;
    region = 'Global';
% $$$     region = 'IndoPacific';
    nyrs = tL/12;szTe = [TL+1 12 nyrs];szT  = [TL 12 nyrs];
    yrs = 1:nyrs;
    months = 1:12;
    
    ycur = 1;

    %% Global Calculations:
    for i=1:length(outputs)
        
        load([base model sprintf('_output%03d_',outputs(i)) region 'HBud.mat']);
        
        % Fluxes:
        P(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.PME+GWB.RMX,szTe); % PME effective heat flux (W)
        F(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SWH+GWB.VDS+GWB.FRZ+GWB.ETS,szTe); % Surface heat flux (W)
        Ffz(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.FRZ,szTe); % Surface heat flux (W)
        Fsw(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SWH,szTe); % Surface heat flux (W)
        Fsh(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDS,szTe); % Surface heat flux (W)
        Fet(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.ETS,szTe); % Surface heat flux (W)
        M(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDF+GWB.KNL,szTe); % Vertical mixing flux (W)
        if (isfield(GWB,'RED')) % Redi Diffusion
            R(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.RED+GWB.K33,szTe); % Redi diffusion (W)
        else
            R(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
        end
        if (isfield(GWB,'NGM')) % GM parameterization
            GM(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.NGM,szTe); % GM (W)
        else
            GM(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
        end    
        if (isfield(GWB,'MDS')) % Mix-downslope
            MD(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.MDS,szTe);; % GM (W)
            M(:,:,ycur:(ycur+nyrs-1)) = M(:,:,ycur:(ycur+nyrs-1)) + reshape(GWB.MDS,szTe); %ADD TO VERTICAL MIXING, but it's small...
        else
            MD(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
        end    
        if (isfield(GWB,'NUM')) % Pre-calculated numerical mixing
            NUM(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.NUM,szTe); % NUM (W)
        else
            NUM(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
        end    
        D(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.TEN-GWB.ADV-GWB.SUB,szTe)-GM(:,:,ycur:(ycur+nyrs-1)); % Material derivative of T (W)
        SW(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SWH,szTe); % Short-wave heat
        JS(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SFW,szTe); % Surface Volume Flux
        SUB(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SUB,szTe);
        % Snapshot fields:
        dVdt(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.dVdt,szTe); % V Change (m3s-1)
        dHdt(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.dHdt,szTe); % H Change (W)

        % Water-mass transformation:
        G(:,:,ycur:(ycur+nyrs-1)) = dVdt(:,:,ycur:(ycur+nyrs-1)) - JS(:,:,ycur:(ycur+nyrs-1)); %Water-mass transformation (m3s-1)

        % Surface Volume flux base flux (not P!)
        JSH(:,:,ycur:(ycur+nyrs-1)) = JS(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 12 nyrs])*rho0*Cp;

        % Interior heat source P:
        PI(:,:,ycur:(ycur+nyrs-1)) = P(:,:,ycur:(ycur+nyrs-1)) - JSH(:,:,ycur:(ycur+nyrs-1));

        % Across-isotherm advective heat flux:
        CIA(:,:,ycur:(ycur+nyrs-1)) = G(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 12 nyrs])*rho0*Cp;

        % External HC Tendency:
        EHC(:,:,ycur:(ycur+nyrs-1)) = dVdt(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 12 nyrs])*rho0*Cp;

        % Internal HC Tendency:
        N(:,:,ycur:(ycur+nyrs-1)) = dHdt(:,:,ycur:(ycur+nyrs-1)) - EHC(:,:,ycur:(ycur+nyrs-1));

        % Implicit mixing:
        I(:,:,ycur:(ycur+nyrs-1)) = N(:,:,ycur:(ycur+nyrs-1)) - F(:,:,ycur:(ycur+nyrs-1)) - P(:,:,ycur:(ycur+nyrs-1)) - M(:,:,ycur:(ycur+nyrs-1)) - R(:,:,ycur:(ycur+nyrs-1)) + JSH(:,:,ycur:(ycur+nyrs-1));

        % Non-advective flux into volume:
        B(:,:,ycur:(ycur+nyrs-1)) = F(:,:,ycur:(ycur+nyrs-1))+M(:,:,ycur:(ycur+nyrs-1))+I(:,:,ycur:(ycur+nyrs-1))+R(:,:,ycur:(ycur+nyrs-1));

        ycur = ycur+nyrs;
    end
    months = [1:length(P(1,:,1))];
    yrs = [1:length(P(1,1,:))];

%%%%Heat Flux: ---------------------------------------------------------------------------------------------
% Production fields:
fields = { ...
          {N(:,months,yrs), 'Internal HC Tendency $\mathcal{N}$','m',2,'-'}, ...
          {F(:,months,yrs)+PI(:,months,yrs), 'Surface Forcing $\mathcal{F}+\mathcal{P}_I$','k',2,'-'}, ...
          {M(:,months,yrs), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
          {R(:,months,yrs), 'Redi Mixing $\mathcal{R}$',[0 0.5 0],2,'-'}, ...
          {I(:,months,yrs), 'Numerical Mixing $\mathcal{I}$','b',2,'-'}, ...
          };

Fscale = 1/1e15;

yrtyps = {'-','--','-.',':'}; % line-types for different years

%Fluxes only:
% $$$ figure;
% $$$ set(gcf,'Position',[207          97        1609         815]);
% $$$ leg = {};
% $$$ legh = [];
for i=1:length(fields)
    hold on;
    if (length(fields{i}{1}(:,1)) == length(Te))
        x = Te;
    else
        x = T;
    end
    
    % Plot years from a single run separately:
% $$$     for j=1:length(yrs) 
% $$$         h = plot(Te,monmean(fields{i}{1}(:,:,yrs(j)),2,ndays(months))*Fscale,yrtyps{j}, 'color',fields{i}{3} ...
% $$$              ,'linewidth',3);
% $$$         if (j == 1)
% $$$             legh(i) = h;
% $$$         end
% $$$     end
% $$$     leg{i} = fields{i}{2};

% $$$     % Average years together for a single run:
% $$$     legh(i) = plot(x,mean(monmean(fields{i}{1},2,ndays(months))*Fscale,3),fields{i}{5}, 'color',fields{i}{3} ...
% $$$          ,'linewidth',fields{i}{4});
% $$$     leg{i} = fields{i}{2};
    
    % Average years together for multiple runs:
    tmp = plot(x,mean(monmean(fields{i}{1},2,ndays(months))*Fscale,3),yrtyps{rr}, 'color',fields{i}{3} ...
         ,'linewidth',fields{i}{4});
    if i==1
        leg{rr} = strrep(RUNS{rr}{1},'_',' ');
        legh(rr) = tmp;
    end
end
ylim([-1.5 1.5]);
xlim([-3 31]);
box on; 
grid on;
ylabel('Heat flux into fluid warmer than $\Theta$ (PW)');
xlabel('Temperature $\Theta$ ($^\circ$C)');
lg = legend(legh,leg);
% $$$ set(lg,'Position',[0.5881    0.5500    0.2041    0.2588]);

