 % This script processes the heat budget and associated variables in
% MOM025 or MOM01 simulations and save's into .mat files

% $$$ baseL = '/short/e14/rmh561/mom/archive/';
% $$$ baseL = '/short/e14/rmh561/access-om2/;
baseL = '/srv/ccrc/data03/z3500785/';

% MOM-SIS025-WOMBAT:
% $$$ model = 'MOM025';
% $$$ baseD = [baseL 'MOM_wombat/']; %Data Directory.
% MOM-SIS025:
model = 'MOM025';
baseD = [baseL 'MOM_HeatDiag/']; %Data Directory.
% ACCESS-OM2:
% $$$ model = 'ACCESS-OM2_025deg_jra55_ryf8485';
% $$$ baseD = [baseL 'control/025deg_jra55_ryf8485/archive/']; %Data Directory.

outD = [baseD 'mat_data/'];
rstbaseD = baseD;

% $$$ post = 'ocean/'; % For ACCESS-OM2 output coulpled;
post = ''; % For MOM-SIS.

haveRedi = 0; % 1 = Redi diffusion is on, 0 = off
haveGM = 0; % 1 = GM is on, 0 = off;
haveMDS = 0; % 1 = MDS is on, 0 = off;

for output = 8:12;
% $$$ output = 8;
restart = output-1;

% file-names -----------------------------------------
base = [baseD sprintf('output%03d/',output) post];
basem1 = [baseD sprintf('output%03d/',output-1) post];
baser = [rstbaseD sprintf('restart%03d/',restart) post];
hname = [base 'ocean_heat.nc'];
if (strfind(baseD,'01'))
    fname = [base 'ocean_month.nc'];
else
    fname = [base 'ocean.nc'];
end
gname = [base 'ocean_grid.nc'];
sname = [base 'ocean_snap.nc'];
wname = [base 'ocean_wmass.nc'];
if (exist(baser))
    found_rst = 1;rstti = 1;
    rnameT = [baser 'ocean_temp_salt.res.nc'];
    rnameZ = [baser 'ocean_thickness.res.nc'];
    rnametime = [baser 'coupler.res'];
    if (~exist(rnametime))
        rnametime = [baser 'ocean_solo.res'];
    end
else
    found_rst = 0;rstti = 12;
    rnameT = [basem1 'ocean_snap.nc'];
    rnameZ = [basem1 'ocean_snap.nc'];
    rnametime = [basem1 'ocean_snap.nc'];
end    
         
% Horizontal Grid  -----------------------------------------
lon = ncread(gname,'geolon_t');lat = ncread(gname,'geolat_t');
lonu = ncread(gname,'geolon_c');latu = ncread(gname,'geolat_c');
area = ncread(gname,'area_t');[xL,yL] = size(lon);
lonv_t = ncread(gname,'xt_ocean');lonv_u = ncread(gname,'xu_ocean');
latv_t = ncread(gname,'yt_ocean');latv_u = ncread(gname,'yu_ocean');

% Vertical grid  -----------------------------------------
z = ncread(hname,'st_ocean');zL = length(z);

% 3D mask ------------------------------------------------
mask = ncread(fname,'temp',[1 1 1 rstti],[xL yL zL 1]);
mask(~isnan(mask)) = 1; mask(isnan(mask)) = 0;
mask = mask == 1;

% Time  -----------------------------------------
time = ncread(hname,'time');

if (found_rst)
    dys = [31 28 31 30 31 30 31 31 30 31 30 31];
    C = textread(rnametime, '%s','delimiter', '\n');
    C = strsplit(C{3});
    rtime = [str2num(C{1}) str2num(C{2}) str2num(C{3}) str2num(C{4}) str2num(C{5}) str2num(C{6})];
    time_snap = [(rtime(1)-1)*365+sum(dys(1:(rtime(2)-1)))+(rtime(3)-1)+rtime(4)/24+rtime(5)/24/60+rtime(6)/24/60/60;
                 ncread(sname,'time')];
else
    time_snapl = ncread(rnametime,'time');
    time_snap = [time_snapl(end); ncread(sname,'time')];
end

time_snap = mod(time_snap,365);
if (time_snap(end) == 0) time_snap(end) = 365;end

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
    dHdtID = netcdf.defVar(ncid,'dHdt','NC_FLOAT',[xid yid zid tid]);
    netcdf.putAtt(ncid,dHdtID,'long_name','Change in time of heat content within temperature bin');
    netcdf.putAtt(ncid,dHdtID,'units','Watts');
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
        Volsnap = ncread(rnameZ,'rho_dzt',[1 1 zi rstti],[xL yL 1 1]).*area/rho0;
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

GWB.TENMON = TENMON; % Tendency from monthly averages (from
                     % previous code block
GWB.dVdt   = zeros(TL+1,tL); % Sv 
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

for ti=1:tL
    ii = TL;
    sprintf('Calculating global water-mass heat budget time %03d of %03d, temp %03d of %03d',ti,tL,ii,TL)
    GWB.dVdt(ii,ti) = nansum(nansum(ncread(wname,'dVdt',[1 1 ii ti],[xL yL 1 1])*1e9/rho0,1),2);
    GWB.dHdt(ii,ti) = nansum(nansum(ncread(wname,'dHdt',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.TEN(ii,ti) = nansum(nansum(area.*ncread(wname,'temp_tendency_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.ADV(ii,ti) = nansum(nansum(area.*ncread(wname,'temp_advection_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.SUB(ii,ti) = nansum(nansum(area.*ncread(wname,'temp_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.PME(ii,ti) = nansum(nansum(area.*ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.RMX(ii,ti) = nansum(nansum(area.*ncread(wname,'temp_rivermix_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDS(ii,ti) = nansum(nansum(area.*ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.SWH(ii,ti) = nansum(nansum(area.*ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDF(ii,ti) = nansum(nansum(area.*ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.KNL(ii,ti) = nansum(nansum(area.*ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    if (haveRedi)
    GWB.K33(ii,ti) = nansum(nansum(area.*ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.RED(ii,ti) = nansum(nansum(area.*ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveGM)
    GWB.NGM(ii,ti) = nansum(nansum(area.*ncread(wname,'neutral_gm_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveMDS)
    GWB.MDS(ii,ti) = nansum(nansum(area.*ncread(wname,'mixdownslope_temp_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    GWB.FRZ(ii,ti) = nansum(nansum(area.*ncread(wname,'frazil_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.ETS(ii,ti) = nansum(nansum(area.*ncread(wname,'temp_eta_smooth_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.SFW(ii,ti) = nansum(nansum(ncread(wname,'mass_pmepr_on_nrho',[1 1 ii ti],[xL yL 1 1])/rho0,1),2);
for ii=TL-1:-1:1
    sprintf('Calculating global water-mass heat budget time %03d of %03d, temp %03d of %03d',ti,tL,ii,TL)
    GWB.dVdt(ii,ti) = GWB.dVdt(ii+1,ti) + nansum(nansum(ncread(wname,'dVdt',[1 1 ii ti],[xL yL 1 1])*1e9/rho0,1),2);
    GWB.dHdt(ii,ti) = GWB.dHdt(ii+1,ti) + nansum(nansum(ncread(wname,'dHdt',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.TEN(ii,ti) = GWB.TEN(ii+1,ti) + nansum(nansum(area.*ncread(wname,'temp_tendency_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.ADV(ii,ti) = GWB.ADV(ii+1,ti) + nansum(nansum(area.*ncread(wname,'temp_advection_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.SUB(ii,ti) = GWB.SUB(ii+1,ti) + nansum(nansum(area.*ncread(wname,'temp_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.PME(ii,ti) = GWB.PME(ii+1,ti) + nansum(nansum(area.*ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.RMX(ii,ti) = GWB.RMX(ii+1,ti) + nansum(nansum(area.*ncread(wname,'temp_rivermix_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDS(ii,ti) = GWB.VDS(ii+1,ti) + nansum(nansum(area.*ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.SWH(ii,ti) = GWB.SWH(ii+1,ti) + nansum(nansum(area.*ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDF(ii,ti) = GWB.VDF(ii+1,ti) + nansum(nansum(area.*ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.KNL(ii,ti) = GWB.KNL(ii+1,ti) + nansum(nansum(area.*ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    if (haveRedi)
    GWB.K33(ii,ti) = GWB.K33(ii+1,ti) + nansum(nansum(area.*ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.RED(ii,ti) = GWB.RED(ii+1,ti) + nansum(nansum(area.*ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveGM)
    GWB.NGM(ii,ti) = GWB.NGM(ii+1,ti) + nansum(nansum(area.*ncread(wname,'neutral_gm_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveMDS)
    GWB.MDS(ii,ti) = GWB.MDS(ii+1,ti) + nansum(nansum(area.*ncread(wname,'mixdownslope_temp_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    GWB.FRZ(ii,ti) = GWB.FRZ(ii+1,ti) + nansum(nansum(area.*ncread(wname,'frazil_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.ETS(ii,ti) = GWB.ETS(ii+1,ti) + nansum(nansum(area.*ncread(wname,'temp_eta_smooth_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.SFW(ii,ti) = GWB.SFW(ii+1,ti) + nansum(nansum(ncread(wname,'mass_pmepr_on_nrho',[1 1 ii ti],[xL yL 1 1])/rho0,1),2);
end
end
save([outD model sprintf('_output%03d',output) '_GlobalHBud.mat'],'GWB','-v7.3');

%% Vertical Integrate down to level from online T-binned values -----------------------------------------------------------------------------------------------------------
Tls = [5 10 15:2.5:27.5];
Nremain = length(Tls);
Ti = TL;

FlM = zeros(xL,yL,tL); % vdiffuse and nonlocal_KPP
FlF = zeros(xL,yL,tL); % surface forcing
FlP = zeros(xL,yL,tL); % P-E+R
FlA = zeros(xL,yL,tL); % advection + submeso + GM
if (haveRedi)
    FlK = zeros(xL,yL,tL); % K33
    FlR = zeros(xL,yL,tL); % Redi
end
if (haveGM)
    FlG = zeros(xL,yL,tL); % GM
end
FlT = zeros(xL,yL,tL); % tendency
FlSP = zeros(xL,yL,tL); % solar penetration

while (Nremain > 0 & Ti >= 1)
    Tl = Te(Ti);

    for ti=1:tL
        sprintf(['Calculating water-mass heat budget time %03d of ' ...
                 '%03d, temp %2.2f, going down to %2.2f'],ti,tL,Te(Ti),min(Tls))
        FlT(:,:,ti) = FlT(:,:,ti)+ncread(wname,'temp_tendency_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        if (haveRedi)
            FlK(:,:,ti) = FlK(:,:,ti)+ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            FlR(:,:,ti) = FlR(:,:,ti)+ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 Ti ti],[xL yL 1 1]);
        end
        FlA(:,:,ti) = FlA(:,:,ti)+ncread(wname,'temp_advection_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
            ncread(wname,'temp_submeso_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        if (haveGM)
            FlG(:,:,ti) = FlG(:,:,ti)+ncread(wname,'neutral_gm_on_nrho_temp',[1 1 Ti ti],[xL yL 1 1]);
            FlA(:,:,ti) = FlA(:,:,ti)+ncread(wname,'neutral_gm_on_nrho_temp',[1 1 Ti ti],[xL yL 1 1]);
        end
        FlP(:,:,ti) = FlP(:,:,ti)+ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
            ncread(wname,'temp_rivermix_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        FlM(:,:,ti) = FlM(:,:,ti)+ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
            ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        FlF(:,:,ti) = FlF(:,:,ti)+ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
            ncread(wname,'sw_heat_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
            ncread(wname,'frazil_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
            ncread(wname,'temp_eta_smooth_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        FlSP(:,:,ti) = FlSP(:,:,ti)+ncread(wname,'sw_heat_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
    end

    % Save heat flux terms:
    [sp,ind] = min(abs(Tls-Tl));
    if (abs(sp) <= dT/4)
        name = [outD model sprintf('_output%03d',output) '_VertInt_T' strrep(num2str(Tls(ind)),'.','p') 'C.mat']

        save(name,'FlM','FlSP','FlF','FlT','FlA','FlP','Tl','-v7.3');
        if (haveRedi)
            save(name,'FlK','FlR','-append');
        end
        if (haveGM)
            save(name,'FlG','-append');
        end
        
        Nremain = Nremain-1;
    end

    Ti = Ti-1;
    
end

%% Calculate WMT due to different (resolved) terms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tls = [5 10 15:2.5:27.5]-0.25;

for ii = 1:length(Tls)
    Tl = Tls(ii);

    WMTM = NaN*zeros(xL,yL,tL); % vdiffuse and nonlocal_KPP
    WMTF = NaN*zeros(xL,yL,tL); % surface forcing
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
        if (haveRedi)
            WMTK(:,:,ti) = 1/rho0/Cp/dT*(ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 Ti ti],[xL yL 1 1]));
            WMTR(:,:,ti) = 1/rho0/Cp/dT*(ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 Ti ti],[xL yL 1 1]));
        end
    end
    save([outD model sprintf('_output%03d',output) '_WMT_T' strrep(num2str(Tl),'.','p') 'C.mat'],'WMTM','WMTSP','WMTP','WMTF','Tl','-v7.3');
    if (haveRedi)
        save([outD model sprintf('_output%03d',output) '_WMT_T' strrep(num2str(Tl),'.','p') 'C.mat'],'WMTK','WMTR','-append');
    end
end

%% Add horizontally-resolved volume fluxes for implicit mixing residual:
WMTTls = [5 10 15:2.5:27.5]-0.25;
FLTls = [5 10 15:2.5:27.5];

Nremain = length(WMTTls)+length(FLTls);
Ti = TL+1;

dVdtM = zeros(xL,yL,tL);
JIM = zeros(xL,yL,tL);
JSM = zeros(xL,yL,tL);

FldVdt = zeros(xL,yL,tL); % rho0*Cp*\int_theta^infty dVdt dtheta
FlJI = zeros(xL,yL,tL); % 
FlJS = zeros(xL,yL,tL); %

JSP = zeros(xL,yL,tL);
JIP = zeros(xL,yL,tL);
dVdtP = zeros(xL,yL,tL);

while (Nremain > 0 & Ti >= 1)
    
    for ti=1:tL
        sprintf(['Calculating JS, JI, dVdt time %03d of %03d, temp ' ...
        '%2.2f, going down to %2.2f'],ti,tL,Te(Ti),min([WMTTls FLTls]))
        JSM(:,:,ti) = JSM(:,:,ti) + ncread(wname,'mass_pmepr_on_nrho',[1 1 Ti-1 ti],[xL yL 1 1])./area/rho0;
        dVdtM(:,:,ti) = dVdtM(:,:,ti) + ncread(wname,'dVdt',[1 1 Ti-1 ti],[xL yL 1 1])*1e9/rho0./area;
        txtrans = ncread(wname,'tx_trans_nrho',[1 1 Ti-1 ti],[xL yL 1 1])*1e9/rho0 + ...
                  ncread(wname,'tx_trans_nrho_submeso',[1 1 Ti-1 ti],[xL yL 1 1])*1e9/rho0;
        tytrans = ncread(wname,'ty_trans_nrho',[1 1 Ti-1 ti],[xL yL 1 1])*1e9/rho0 + ...
                  ncread(wname,'ty_trans_nrho_submeso',[1 1 Ti-1 ti],[xL yL 1 1])*1e9/rho0;
        if (haveGM)
            txtrans = txtrans + ncread(wname,'tx_trans_nrho_gm',[1 1 Ti-1 ti],[xL yL 1 1])*1e9/rho0;
            tytrans = tytrans + ncread(wname,'ty_trans_nrho_gm',[1 1 Ti-1 ti],[xL yL 1 1])*1e9/rho0;
        end
            
        JIM(2:end,2:end,ti) = JIM(2:end,2:end,ti)+(txtrans(1:(end-1),2:end) - txtrans(2:end,2:end) ...
                +tytrans(2:end,1:(end-1)) - tytrans(2:end,2:end))./area(2:end,2:end);
        JIM(1,2:end,ti) = JIM(1,2:end,ti)+(txtrans(end,2:end) - txtrans(1,2:end) ...
            +tytrans(1,1:(end-1)) - tytrans(1,2:end))./area(1,2:end);        
    end    
    
    dVdt = (dVdtM + dVdtP)/2;
    JS = (JSM + JSP)/2;
    JI = (JIM + JIP)/2;
    
    FldVdt = FldVdt + rho0*Cp*dVdt*dT;
    FlJS = FlJS + rho0*Cp*JS*dT;
    FlJI = FlJI + rho0*Cp*JI*dT;
        
    % Save WMT terms:
    WMT_temp = Te(Ti)-dT/2;
    sp = min(abs(WMT_temp-WMTTls));
    if (abs(sp) <= dT/4)
        name = [outD model sprintf('_output%03d',output) '_WMT_T' strrep(num2str(WMT_temp),'.','p') 'C.mat']
        save(name,'dVdt','JS','JI','-append');

        % Calculate implicit mixing by residual:
        load(name,'WMTM','WMTF');
        WMTI = dVdt-WMTM-WMTF-JI-JS;
        if (haveRedi)
            load(name,'WMTK','WMTR');
            WMTI = WMTI-WMTK-WMTR;
        end
        save(name,'WMTI','-append');
        
        Nremain = Nremain-1;
    end

    % Save heat flux terms:
    Fl_temp = Te(Ti-1);
    sp = min(abs(Fl_temp-FLTls));
    if (abs(sp) <= dT/4)
        name = [outD model sprintf('_output%03d',output) '_VertInt_T' strrep(num2str(Fl_temp),'.','p') 'C.mat']
        save(name,'FldVdt','FlJS','FlJI','-append');

        % Calculate implicit mixing by residual:
        load(name,'FlM','FlF');
        FlI = FldVdt-FlM-FlF-FlJI-FlJS;
        if (haveRedi)
            load(name,'FlK','FlR');
            FlI = FlI-FlK-FlR;
        end
        save(name,'FlI','-append');
    
        Nremain = Nremain-1;
    end

    JSP = JSM;
    JIP = JIM;
    dVdtP = dVdtM;
    Ti = Ti-1;
end

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

%% Save surface heat flux, wind stress, SST, meridional heat flux:
shflux = ncread(fname,'net_sfc_heating',[1 1 1],[xL yL tL]);
SST = squeeze(ncread(fname,'temp',[1 1 1 1],[xL yL 1 tL]));
taux = ncread(fname,'tau_x',[1 1 1],[xL yL tL]);
tauy = ncread(fname,'tau_y',[1 1 1],[xL yL tL]);

save([outD model sprintf('_output%03d',output) '_SurfaceVars.mat'],'shflux','SST','taux','tauy');
p
% Do meridional heat flux:
mhflux = zeros(yL,tL);
for ti = 1:tL
    for Ti = 1:TL
        sprintf('Calculating meridional heat flux time %03d of %03d, temp %03d of %03d',ti,tL,Ti,TL)
        tytrans = ncread(wname,'ty_trans_nrho',[1 1 Ti ti],[xL yL 1 1])*1e9/rho0 + ...
                  ncread(wname,'ty_trans_nrho_submeso',[1 1 Ti ti],[xL yL 1 1])*1e9/rho0;
        if (haveGM)
            tytrans = tytrans + ncread(wname,'ty_trans_nrho_gm',[1 1 Ti ti],[xL yL 1 1])*1e9/rho0;
        end
        mhflux(:,ti) = mhflux(:,ti) + rho0*Cp*T(Ti)*nansum(tytrans,1)';
    end
end
save([outD model sprintf('_output%03d',output) '_SurfaceVars.mat'],'mhflux','-append');

% $$$ 
% $$$ %% Latitude-depth overturning:
% $$$ PSI = zeros(yL,TL,tL);
% $$$ 
% $$$ % Ignoring tri-polar (wrong >60N).
% $$$ for ti=1:tL
% $$$     ti
% $$$     PSI(:,:,ti) = squeeze(nansum(ncread(wname,'ty_trans_nrho',[1 1 1 ti],[xL yL TL 1]),1));
% $$$ end
% $$$ save([outD model sprintf('_output%03d',output) '_Tpsi.mat'],'PSI');

% $$$     
% $$$     for TI = 1:TL
% $$$         tytrans = ncread(wname,'ty_trans_nrho',[1 1 TI ti],[xL yL 1 1]);
% $$$         for yi=1:yL
% $$$             inds = lat >= late(yi) & lat <= 
% $$$             PSI(yi,TI,ti) = sum(tytrans(


% $$$ %% Heat Function From tendencies -------------------------------------------------
% $$$ HFSWH    = zeros(yL,TL+1,tL); % W due to SW redistribution
% $$$ HFVDS    = zeros(yL,TL+1,tL); % W due to vdiffuse_sbc.
% $$$ HFRMX    = zeros(yL,TL+1,tL); % W due to rivermix.
% $$$ HFPME    = zeros(yL,TL+1,tL); % W due to P-E.
% $$$ HFFRZ    = zeros(yL,TL+1,tL); % W due to frazil.
% $$$ HFETS    = zeros(yL,TL+1,tL); % W due to eta_smoothing.
% $$$ HFSUB    = zeros(yL,TL+1,tL); % W due to submesoscale.
% $$$ HFVDF    = zeros(yL,TL+1,tL); % W due to vdiffusion
% $$$ HFKNL    = zeros(yL,TL+1,tL); % W due to KPP non-local
% $$$ HFADV    = zeros(yL,TL+1,tL); % W due to advection
% $$$ HFTEN    = zeros(yL,TL+1,tL); % W due to tendency
% $$$ 
% $$$ tph = zeros(yL,1);
% $$$ for yi=1:yL
% $$$     tph(yi) = all(lat(:,yi) == latv(yi));
% $$$ end
% $$$ tph = ~tph;
% $$$ tplast = find(tph==1,1,'first');
% $$$ 
% $$$ for ti=1:tL
% $$$     for ii=1:TL
% $$$         TEN = area.*ncread(wname,'temp_tendency_on_nrho',[1 1 ii ti],[xL yL 1 1]);TEN(isnan(TEN)) = 0;
% $$$         ADV = area.*ncread(wname,'temp_advection_on_nrho',[1 1 ii ti],[xL yL 1 1]);ADV(isnan(ADV)) = 0;
% $$$         SUB = area.*ncread(wname,'temp_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]);SUB(isnan(SUB)) = 0;
% $$$         PME = area.*ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 ii ti],[xL yL 1 1]);PME(isnan(PME)) = 0;
% $$$         RMX = area.*ncread(wname,'temp_rivermix_on_nrho',[1 1 ii ti],[xL yL 1 1]);RMX(isnan(RMX)) = 0;
% $$$         VDS = area.*ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL 1 1]);VDS(isnan(VDS)) = 0;
% $$$         SWH = area.*ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1]);SWH(isnan(SWH)) = 0;
% $$$         VDF = area.*ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1]);VDF(isnan(VDF)) = 0;
% $$$         KNL = area.*ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]);KNL(isnan(KNL)) = 0;
% $$$         FRZ = area.*ncread(wname,'frazil_on_nrho',[1 1 ii ti],[xL yL 1 1]);FRZ(isnan(FRZ)) = 0;
% $$$         ETS = area.*ncread(wname,'temp_eta_smooth_on_nrho',[1 1 ii ti],[xL yL 1 1]);ETS(isnan(ETS)) = 0;
% $$$         
% $$$         % Tripolar region:
% $$$         for yi=yL:-1:(tplast-1)
% $$$             sprintf('Calculating heat function lat %03d of %03d, time %03d of %03d, temp %03d of %03d',yL-yi+1,yL,ti,tL,ii,TL)
% $$$             inds = lat >= latv(yi);
% $$$             HFTEN(yi,ii+1,ti) = HFTEN(yi,ii,ti) + sum(TEN(inds));
% $$$             HFADV(yi,ii+1,ti) = HFADV(yi,ii,ti) + sum(ADV(inds));
% $$$             HFSUB(yi,ii+1,ti) = HFSUB(yi,ii,ti) + sum(SUB(inds));
% $$$             HFPME(yi,ii+1,ti) = HFPME(yi,ii,ti) + sum(PME(inds));
% $$$             HFRMX(yi,ii+1,ti) = HFRMX(yi,ii,ti) + sum(RMX(inds));
% $$$             HFVDS(yi,ii+1,ti) = HFVDS(yi,ii,ti) + sum(VDS(inds));
% $$$             HFSWH(yi,ii+1,ti) = HFSWH(yi,ii,ti) + sum(SWH(inds));
% $$$             HFVDF(yi,ii+1,ti) = HFVDF(yi,ii,ti) + sum(VDF(inds));
% $$$             HFKNL(yi,ii+1,ti) = HFKNL(yi,ii,ti) + sum(KNL(inds));
% $$$             HFFRZ(yi,ii+1,ti) = HFFRZ(yi,ii,ti) + sum(FRZ(inds));
% $$$             HFETS(yi,ii+1,ti) = HFETS(yi,ii,ti) + sum(ETS(inds));
% $$$         end
% $$$         %Non-tripolar region:
% $$$         HFTEN(1:(tplast-2),ii+1,ti) = HFTEN(1:(tplast-2),ii,ti) + cumsum(sum(TEN(:,1:(tplast-2)),1),2,'reverse')' + sum(sum(TEN(:,(tplast-1):end),1),2);
% $$$         HFADV(1:(tplast-2),ii+1,ti) = HFADV(1:(tplast-2),ii,ti) + cumsum(sum(ADV(:,1:(tplast-2)),1),2,'reverse')' + sum(sum(ADV(:,(tplast-1):end),1),2);
% $$$         HFSUB(1:(tplast-2),ii+1,ti) = HFSUB(1:(tplast-2),ii,ti) + cumsum(sum(SUB(:,1:(tplast-2)),1),2,'reverse')' + sum(sum(SUB(:,(tplast-1):end),1),2);
% $$$         HFPME(1:(tplast-2),ii+1,ti) = HFPME(1:(tplast-2),ii,ti) + cumsum(sum(PME(:,1:(tplast-2)),1),2,'reverse')' + sum(sum(PME(:,(tplast-1):end),1),2);
% $$$         HFRMX(1:(tplast-2),ii+1,ti) = HFRMX(1:(tplast-2),ii,ti) + cumsum(sum(RMX(:,1:(tplast-2)),1),2,'reverse')' + sum(sum(RMX(:,(tplast-1):end),1),2);
% $$$         HFVDS(1:(tplast-2),ii+1,ti) = HFVDS(1:(tplast-2),ii,ti) + cumsum(sum(VDS(:,1:(tplast-2)),1),2,'reverse')' + sum(sum(VDS(:,(tplast-1):end),1),2);
% $$$         HFSWH(1:(tplast-2),ii+1,ti) = HFSWH(1:(tplast-2),ii,ti) + cumsum(sum(SWH(:,1:(tplast-2)),1),2,'reverse')' + sum(sum(SWH(:,(tplast-1):end),1),2);
% $$$         HFVDF(1:(tplast-2),ii+1,ti) = HFVDF(1:(tplast-2),ii,ti) + cumsum(sum(VDF(:,1:(tplast-2)),1),2,'reverse')' + sum(sum(VDF(:,(tplast-1):end),1),2);
% $$$         HFKNL(1:(tplast-2),ii+1,ti) = HFKNL(1:(tplast-2),ii,ti) + cumsum(sum(KNL(:,1:(tplast-2)),1),2,'reverse')' + sum(sum(KNL(:,(tplast-1):end),1),2);
% $$$         HFFRZ(1:(tplast-2),ii+1,ti) = HFFRZ(1:(tplast-2),ii,ti) + cumsum(sum(FRZ(:,1:(tplast-2)),1),2,'reverse')' + sum(sum(FRZ(:,(tplast-1):end),1),2);
% $$$         HFETS(1:(tplast-2),ii+1,ti) = HFETS(1:(tplast-2),ii,ti) + cumsum(sum(ETS(:,1:(tplast-2)),1),2,'reverse')' + sum(sum(ETS(:,(tplast-1):end),1),2);
% $$$ % $$$ 
% $$$ % $$$         for yi=(tplast-2):-1:1
% $$$ % $$$             sprintf('Calculating heat function lat %03d of %03d, time %03d of %03d, temp %03d of %03d',yL-yi+1,yL,ti,tL,ii,TL)
% $$$ % $$$             HFTEN(yi,ii+1,ti) = HFTEN(yi,ii,ti) + nansum(nansum(TEN(:,yi:end),1),2);
% $$$ % $$$             HFADV(yi,ii+1,ti) = HFADV(yi,ii,ti) + nansum(nansum(ADV(:,yi:end),1),2);
% $$$ % $$$             HFSUB(yi,ii+1,ti) = HFSUB(yi,ii,ti) + nansum(nansum(SUB(:,yi:end),1),2);
% $$$ % $$$             HFPME(yi,ii+1,ti) = HFPME(yi,ii,ti) + nansum(nansum(PME(:,yi:end),1),2);
% $$$ % $$$             HFRMX(yi,ii+1,ti) = HFRMX(yi,ii,ti) + nansum(nansum(RMX(:,yi:end),1),2);
% $$$ % $$$             HFVDS(yi,ii+1,ti) = HFVDS(yi,ii,ti) + nansum(nansum(VDS(:,yi:end),1),2);
% $$$ % $$$             HFSWH(yi,ii+1,ti) = HFSWH(yi,ii,ti) + nansum(nansum(SWH(:,yi:end),1),2);
% $$$ % $$$             HFVDF(yi,ii+1,ti) = HFVDF(yi,ii,ti) + nansum(nansum(VDF(:,yi:end),1),2);
% $$$ % $$$             HFKNL(yi,ii+1,ti) = HFKNL(yi,ii,ti) + nansum(nansum(KNL(:,yi:end),1),2);
% $$$ % $$$             HFFRZ(yi,ii+1,ti) = HFFRZ(yi,ii,ti) + nansum(nansum(FRZ(:,yi:end),1),2);
% $$$ % $$$             HFETS(yi,ii+1,ti) = HFETS(yi,ii,ti) + nansum(nansum(ETS(:,yi:end),1),2);
% $$$ % $$$         end
% $$$     end
% $$$ end
% $$$ save([outD model sprintf('_output%03d',output) '_HFunc.mat'],'HFSWH','HFVDS','HFRMX','HFPME','HFFRZ', ...
% $$$      'HFETS','HFSUB','HFVDF','HFKNL','HFADV','HFTEN');

end

% $$$ %% Swap in non-NaN'd lon/lat:
% $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/mat_data/';
% $$$ model = 'MOM025';
% $$$ outputs = [8:12];
% $$$ 
% $$$ load([base 'old/' model sprintf('_output%03d_BaseVars.mat',2)]);
% $$$ region = 'Global';
% $$$ for output = 8:12
% $$$     save([base model sprintf('_output%03d_BaseVars.mat',output)], ...
% $$$          'lon','lat','lonu','latu','area','-append');
% $$$ end
