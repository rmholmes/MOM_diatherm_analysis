% This script processes the heat budget and associated variables in
% MOM025 or MOM01 simulations and save's into .mat files. This version
% saves only a specified region (as determined by Heat_Budget_Mask)
% and transports/fluxes into and out of this region.

% $$$ baseL = '/short/e14/rmh561/mom/archive/';
baseL = '/srv/ccrc/data03/z3500785/';

% $$$ model = 'MOM01';
% $$$ baseD = [baseL 'MOM01_HeatDiag/'];

% $$$ model = 'MOM025_nipoall';
% $$$ baseD = [baseL 'MOM_HeatDiag_nipoall/'];

model = 'MOM025';
baseD = [baseL 'MOM_HeatDiag/'];

rstbaseD = baseD;
outD = [baseD 'mat_data/'];

region = 'Pacific'; % See Heat_Budget_Mask

% $$$ post = 'ocean/'; % For ACCESS-OM2 output coulpled;
post = ''; % For MOM-SIS.

haveRedi = 0;
haveGM = 0;
% $$$ for output=[19:-1:1]
    output = 8;
% $$$     output = 19;
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

% Temperature grid  -----------------------------------------
Cp = 3992.1; % J kg-1 degC-1
rho0 = 1035; % kgm-3

T = ncread(wname,'neutral');
Te = ncread(wname,'neutralrho_edges');
TL = length(T);dT = T(2)-T(1);

[mask_t,mask_Ny,mask_Nx,mask_Sx,mask_Sy,mask_Wx,mask_Wy] = ...
    Heat_Budget_Mask('Pacific',gname,fname,wname,outD,model);

save([outD model sprintf('_output%03d',output) '_' region '_BaseVars.mat'], ...
     'T','Te','TL','dT','Cp','rho0','time','time_snap','tL', ...
     'z','zL','lon','lat','area','xL','yL','latv_u','lonv_t','latv_t','lonv_u', ...
     'lonu','latu','mask_t','mask_Ny','mask_Nx','mask_Sx','mask_Sy','mask_Wx','mask_Wy','-v7.3');

%% Calculate volume integrated budget from online T-binned values -----------------------------------------------------------------------------------------------------------
V      = zeros(TL+1,tL); % Volume of water (m3) above temperature T
H      = zeros(TL+1,tL); % Heat content (J) above temperature T
Vsnap  = zeros(TL+1,tL+1); % Volume from snapshots (m3)
Hsnap  = zeros(TL+1,tL+1); % Heat content from snapshots (J)
Temp   = zeros(zL,tL); % Temperature as a function of depth (degC)
SWH    = zeros(TL+1,tL); % W due to SW redistribution
VDS    = zeros(TL+1,tL); % W due to vdiffuse_sbc.
RMX    = zeros(TL+1,tL); % W due to rivermix.
PME    = zeros(TL+1,tL); % W due to P-E.
FRZ    = zeros(TL+1,tL); % W due to frazil.
ETS    = zeros(TL+1,tL); % W due to eta_smoothing.
SUB    = zeros(TL+1,tL); % W due to submesoscale.
VDF    = zeros(TL+1,tL); % W due to vdiffusion
KNL    = zeros(TL+1,tL); % W due to KPP non-local
if (haveRedi)
    K33    = zeros(TL+1,tL); % W due to K33
    RED    = zeros(TL+1,tL); % W due to Redi diffusion
end
if (haveGM)
    NGM    = zeros(TL+1,tL); % W due to GM
end
ADV    = zeros(TL+1,tL); % W due to advection
TEN    = zeros(TL+1,tL); % W due to tendency
SFW    = zeros(TL+1,tL); % surface volume flux into ocean (m3s-1)
TENMON = zeros(TL+1,tL); % W due to tendency from Offline Monthly
JBS    = zeros(TL+1,tL);  % m3s-1 out of Pacific North
JSP    = zeros(TL+1,tL);  % m3s-1 out of Pacific South
JITF   = zeros(TL+1,tL);  % m3s-1 out of Pacific West
QBS    = zeros(TL+1,tL);  % W out of Pacific North
QSP    = zeros(TL+1,tL);  % W out of Pacific South
QITF   = zeros(TL+1,tL);  % W out of Pacific West

%Do IC for Vsnap and Hsnap:
for zi = 1:zL
    sprintf('Doing snapshot IC, depth %02d of %02d',zi,zL)
    %Temperature snapshot:
    tempsnap = ncread(rnameT,'temp',[1 1 zi rstti],[xL yL 1 1]);
    tempsnap(tempsnap==0) = NaN; %This is included because the
                                 %restarts don't have any NaNs in
                                 %temp, just lots of 0s. 
    if (found_rst)
        Volsnap = ncread(rnameZ,'rho_dzt',[1 1 zi rstti],[xL yL 1 1]).*area/rho0;
    else
        Volsnap = ncread(rnameT,'dzt',[1 1 zi rstti],[xL yL 1 1]).*area;
    end
    
    % Mask in Pacific-only:
    tempsnap(~mask_t) = NaN;
    Volsnap(~mask_t) = NaN;
    
    %Accumulate sums:
    for Ti=1:TL
        inds = find(tempsnap>=Te(Ti) & tempsnap<Te(Ti+1));
        Vsnap(Ti,1) = Vsnap(Ti,1)+nansum(Volsnap(inds));
        Hsnap(Ti,1) = Hsnap(Ti,1)+nansum(Volsnap(inds).*tempsnap(inds)*rho0*Cp);
    end
    inds = find(tempsnap>=Te(TL+1));
    Vsnap(TL+1,1) = Vsnap(TL+1,1)+nansum(Volsnap(inds));
    Hsnap(TL+1,1) = Hsnap(TL+1,1)+nansum(Volsnap(inds).*tempsnap(inds)*rho0*Cp);
end
%Integrate to get to T'>T:
Vsnap(:,1) = flipud(cumsum(flipud(Vsnap(:,1))));
Hsnap(:,1) = flipud(cumsum(flipud(Hsnap(:,1))));

%Do Eulerian budget calculations:
for ti=1:tL
    for zi=1:zL
        sprintf('Doing Eul Bud. time %03d of %03d, depth %02d of %02d',ti,tL,zi,zL)

        temp = ncread(fname,'temp',[1 1 zi ti],[xL yL 1 1]);
        tempsnap = ncread(sname,'temp',[1 1 zi ti],[xL yL 1 1]);
        Vol = ncread(fname,'dzt',[1 1 zi ti],[xL yL 1 1]).*area;
        Volsnap = ncread(sname,'dzt',[1 1 zi ti],[xL yL 1 1]).*area;

        % Mask in Pacific-only:
        tempsnap(~mask_t) = NaN;
        Volsnap(~mask_t) = NaN;
        temp(~mask_t) = NaN;
        Vol(~mask_t) = NaN;

        %Calculate T(z):
        areaNaN = area;
        areaNaN(isnan(temp)) = NaN;
        Temp(zi,ti) = squeeze(nansum(nansum(temp.*area,1),2)./nansum(nansum(areaNaN,1),2));

        %Tendency from Monthly snapshots:
        TENf = area.*ncread(hname,'temp_tendency',[1 1 zi ti],[xL yL 1 1]);
        TENf(~mask_t) = NaN;
        
        %Accumulate sums:
        for Ti=1:TL
            inds = find(tempsnap>=Te(Ti) & tempsnap<Te(Ti+1));
            Vsnap(Ti,ti+1) = Vsnap(Ti,ti+1)+nansum(Volsnap(inds));
            Hsnap(Ti,ti+1) = Hsnap(Ti,ti+1)+nansum(Volsnap(inds).*tempsnap(inds)*rho0*Cp);
            inds = find(temp>=Te(Ti) & temp<Te(Ti+1));
            V(Ti,ti) = V(Ti,ti)+nansum(Vol(inds));
            H(Ti,ti) = H(Ti,ti)+nansum(Vol(inds).*temp(inds)*rho0*Cp);
            TENMON(Ti,ti) = TENMON(Ti,ti)+nansum(TENf(inds));            
        end
        inds = find(tempsnap>=Te(TL+1));
        Vsnap(TL+1,ti+1) = Vsnap(TL+1,ti+1)+nansum(Volsnap(inds));
        Hsnap(TL+1,ti+1) = Hsnap(TL+1,ti+1)+nansum(Volsnap(inds).*tempsnap(inds)*rho0*Cp);
        inds = find(temp>=Te(TL+1));
        V(TL+1,ti) = V(TL+1,ti)+nansum(Vol(inds));
        H(TL+1,ti) = H(TL+1,ti)+nansum(Vol(inds).*temp(inds)*rho0*Cp);
        TENMON(TL+1,ti) = TENMON(TL+1,ti)+nansum(TENf(inds));            
    end
    %Integrate to get to T'>T:
    Vsnap(:,ti+1) = flipud(cumsum(flipud(Vsnap(:,ti+1))));
    V(:,ti) = flipud(cumsum(flipud(V(:,ti))));
    Hsnap(:,ti+1) = flipud(cumsum(flipud(Hsnap(:,ti+1))));
    H(:,ti) = flipud(cumsum(flipud(H(:,ti))));
    TENMON(:,ti) = flipud(cumsum(flipud(TENMON(:,ti))));
end

save([outD model sprintf('_output%03d',output) '_' region '_HBud.mat'], ...
     'Vsnap','Hsnap','V','H','Temp','TENMON','-v7.3');

for ti=1:tL
    ii = TL;
    sprintf('Calculating water-mass heat budget time %03d of %03d, temp %03d of %03d',ti,tL,ii,TL)
    TEN(ii,ti) = nansum(nansum(mask_t.*area.*ncread(wname,'temp_tendency_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    ADV(ii,ti) = nansum(nansum(mask_t.*area.*ncread(wname,'temp_advection_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    SUB(ii,ti) = nansum(nansum(mask_t.*area.*ncread(wname,'temp_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    PME(ii,ti) = nansum(nansum(mask_t.*area.*ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    RMX(ii,ti) = nansum(nansum(mask_t.*area.*ncread(wname,'temp_rivermix_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    VDS(ii,ti) = nansum(nansum(mask_t.*area.*ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    SWH(ii,ti) = nansum(nansum(mask_t.*area.*ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    VDF(ii,ti) = nansum(nansum(mask_t.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    KNL(ii,ti) = nansum(nansum(mask_t.*area.*ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    if (haveRedi)
    K33(ii,ti) = nansum(nansum(mask_t.*area.*ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    RED(ii,ti) = nansum(nansum(mask_t.*area.*ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveGM)
    NGM(ii,ti) = nansum(nansum(mask_t.*area.*ncread(wname,'neutral_gm_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    FRZ(ii,ti) = nansum(nansum(mask_t.*area.*ncread(wname,'frazil_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    ETS(ii,ti) = nansum(nansum(mask_t.*area.*ncread(wname,'temp_eta_smooth_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    SFW(ii,ti) = nansum(nansum(mask_t.*ncread(wname,'mass_pmepr_on_nrho',[1 1 ii ti],[xL yL 1 1])/rho0,1),2);
    
    txtrans = ncread(wname,'tx_trans_nrho',[1 1 ii ti],[xL yL 1 1])*1e9/rho0;
    tytrans = ncread(wname,'ty_trans_nrho',[1 1 ii ti],[xL yL 1 1])*1e9/rho0;
    txtrans = txtrans+ncread(wname,'tx_trans_nrho_submeso',[1 1 ii ti],[xL yL 1 1])*1e9/rho0;
    tytrans = tytrans+ncread(wname,'ty_trans_nrho_submeso',[1 1 ii ti],[xL yL 1 1])*1e9/rho0;
    if (haveGM)
    txtrans = txtrans+ncread(wname,'tx_trans_nrho_gm',[1 1 ii ti],[xL yL 1 1])*1e9/rho0;
    tytrans = tytrans+ncread(wname,'ty_trans_nrho_gm',[1 1 ii ti],[xL yL 1 1])*1e9/rho0;
    end
    
    JBS(ii,ti) = nansum(tytrans(mask_Ny==1));
    JSP(ii,ti) = -nansum(tytrans(mask_Sy==1)) - nansum(txtrans(mask_Sx==1));
    JITF(ii,ti) = -nansum(txtrans(mask_Wx==1));
    QBS(ii,ti) = JBS(ii,ti)*rho0*Cp*T(ii);
    QSP(ii,ti) = JSP(ii,ti)*rho0*Cp*T(ii);
    QITF(ii,ti) = JITF(ii,ti)*rho0*Cp*T(ii);
    
for ii=TL-1:-1:1
    sprintf('Calculating water-mass heat budget time %03d of %03d, temp %03d of %03d',ti,tL,ii,TL)
    TEN(ii,ti) = TEN(ii+1,ti) + nansum(nansum(mask_t.*area.*ncread(wname,'temp_tendency_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    ADV(ii,ti) = ADV(ii+1,ti) + nansum(nansum(mask_t.*area.*ncread(wname,'temp_advection_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    SUB(ii,ti) = SUB(ii+1,ti) + nansum(nansum(mask_t.*area.*ncread(wname,'temp_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    PME(ii,ti) = PME(ii+1,ti) + nansum(nansum(mask_t.*area.*ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    RMX(ii,ti) = RMX(ii+1,ti) + nansum(nansum(mask_t.*area.*ncread(wname,'temp_rivermix_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    VDS(ii,ti) = VDS(ii+1,ti) + nansum(nansum(mask_t.*area.*ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    SWH(ii,ti) = SWH(ii+1,ti) + nansum(nansum(mask_t.*area.*ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    VDF(ii,ti) = VDF(ii+1,ti) + nansum(nansum(mask_t.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    KNL(ii,ti) = KNL(ii+1,ti) + nansum(nansum(mask_t.*area.*ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    if (haveRedi)
    K33(ii,ti) = K33(ii+1,ti) + nansum(nansum(mask_t.*area.*ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    RED(ii,ti) = RED(ii+1,ti) + nansum(nansum(mask_t.*area.*ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveGM)
    NGM(ii,ti) = NGM(ii+1,ti) + nansum(nansum(mask_t.*area.*ncread(wname,'neutral_gm_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    FRZ(ii,ti) = FRZ(ii+1,ti) + nansum(nansum(mask_t.*area.*ncread(wname,'frazil_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    ETS(ii,ti) = ETS(ii+1,ti) + nansum(nansum(mask_t.*area.*ncread(wname,'temp_eta_smooth_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    SFW(ii,ti) = SFW(ii+1,ti) + nansum(nansum(mask_t.*ncread(wname,'mass_pmepr_on_nrho',[1 1 ii ti],[xL yL 1 1])/rho0,1),2);

    txtrans = ncread(wname,'tx_trans_nrho',[1 1 ii ti],[xL yL 1 1])*1e9/rho0;
    tytrans = ncread(wname,'ty_trans_nrho',[1 1 ii ti],[xL yL 1 1])*1e9/rho0;
    txtrans = txtrans+ncread(wname,'tx_trans_nrho_submeso',[1 1 ii ti],[xL yL 1 1])*1e9/rho0;
    tytrans = tytrans+ncread(wname,'ty_trans_nrho_submeso',[1 1 ii ti],[xL yL 1 1])*1e9/rho0;
    if (haveGM)
    txtrans = txtrans+ncread(wname,'tx_trans_nrho_gm',[1 1 ii ti],[xL yL 1 1])*1e9/rho0;
    tytrans = tytrans+ncread(wname,'ty_trans_nrho_gm',[1 1 ii ti],[xL yL 1 1])*1e9/rho0;
    end
    
    JBS(ii,ti) = JBS(ii+1,ti) + nansum(tytrans(mask_Ny==1));
    JSP(ii,ti) = JSP(ii+1,ti) - nansum(tytrans(mask_Sy==1)) - nansum(txtrans(mask_Sx==1));
    JITF(ii,ti) = JITF(ii+1,ti) - nansum(txtrans(mask_Wx==1));
    QBS(ii,ti) = QBS(ii+1,ti) + nansum(tytrans(mask_Ny==1))*rho0*Cp*T(ii);
    QSP(ii,ti) = QSP(ii+1,ti) + (- nansum(tytrans(mask_Sy==1)) - nansum(txtrans(mask_Sx==1)))*rho0*Cp*T(ii);
    QITF(ii,ti) = QITF(ii+1,ti) + (- nansum(txtrans(mask_Wx==1)))*rho0*Cp*T(ii);
end
end

save([outD model sprintf('_output%03d',output) '_' region '_HBud.mat'],'SWH','VDS','RMX','PME','FRZ', ...
     'ETS','SUB','VDF','KNL','ADV','TEN','SFW','JBS','JSP','JITF','QBS','QSP','QITF','-append');
if (haveRedi)
    save([outD model sprintf('_output%03d',output) '_' region '_HBud.mat'],'K33','RED','-append');
end
if (haveGM)
    save([outD model sprintf('_output%03d',output) '_' region '_HBud.mat'],'NGM','-append');
end

